# pgr 重复标记（Repeat Masking）方案

> 本文档是 pgr 重复标记的总体方案：现状命令（e-kmer/s-kmer/trf 及其实现管道）、
> 命令命名规划、遮蔽版计划（Dfam 全库 + pgi/lastz 比对）、以及与 RepeatMasker
> 的差异。用户用法见 [docs/rept.md](../../docs/rept.md)；FastK 底层分析见
> [[../references/fastk.md]]；RepeatMasker 源码梳理见附录 A（依据仓库内
> `RepeatMasker/` 目录，open-4.2.4）。

## 1. 场景与现状

### 1.1 背景与动机

完整跑一遍 RepeatMasker（基于 Smith-Waterman 与 RepBase/Dfam 库）非常昂贵且耗时。
pgr 的思路是用 FastK 系工具做**快速近似**：不做逐条 repeat 的分类注释，只回答
"基因组上哪些区间是重复的"，然后把区间喂给 `pgr fa mask` 做 soft/hard masking。
这适合大规模基因组的快速屏蔽；需要注释级结果（family/class 标签）时仍应使用 RepeatMasker。

外部工具依赖（均须在 `$PATH`）：

*   FastK + Profex（FastK 套件）
*   spanr
*   trf（仅 `rept trf`）

### 1.2 命令分工与检测闭环

| 命令 | 重复类型 | 原理 | 输入 |
| :--- | :--- | :--- | :--- |
| `pgr rept e-kmer` | 散在重复（interspersed） | 与重复库做 k-mer 富集比对 | 基因组 + repeat 库（Dfam/RepBase/TnCentral） |
| `pgr rept s-kmer` | 基因组内重复（无库） | 自身 k-mer 深度比较 | 仅基因组 |
| `pgr rept trf` | 串联重复（tandem） | trf 的周期搜索 | 仅基因组 |

三者输出**同一种格式**：runlist JSON（`{"chr": "start-end,start-end,..."}`），
可直接作为 `pgr fa mask --runlist` 的输入，因此检测结果与屏蔽步骤天然闭环：

```bash
# 检测 + 屏蔽闭环
pgr rept s-kmer genome.fa -o repeats.json
pgr fa mask genome.fa --runlist repeats.json -o masked.fa        # soft-mask（小写）
pgr fa mask genome.fa --runlist repeats.json --hard -o masked.fa # hard-mask（N）
```

> **命名说明**：`pgr rept` 组已于 2026-08-03 迁移落地，原 `pgr pl ir` /
> `pgr pl rept` / `pgr pl trf` 分别更名为 `pgr rept e-kmer` / `pgr rept s-kmer` /
> `pgr rept trf`（命名规则见 §1.3）。

### 1.3 命令命名规划

#### 设计过程

1. **问题**：`ir` / `rept` 缩写脱离上下文、难以理解；且未来还有两个命令
   （库 + align 遮蔽版、自身 + align 自比对），加上现有的共 4 个组合，需要
   成体系的命名。
2. **头部命令**：仿照 `pgr align` / `pgr sd`，建 `pgr rept` 组，子命令承载
   4 个组合。
3. **维度取舍**：先讨论前缀表达哪个维度。机制（kmer/align）作前缀的直觉是
   "技术直白"，但机制是可扩展维度（未来可能加 HMM 等），且 pgr 没有机制进
   命令名的先例（`pgr sd search --engine pgi|lastz` 的引擎是参数）。最终选
   **对象（库/自身）作前缀**。
4. **形式**：单字母前缀 + 机制单词后缀——`e-kmer` / `e-align` / `s-kmer` /
   `s-align`，避免单字母命令（`e`/`s` 单独作命令名）的不可读问题，也比
   `repeat-lib` 这类全词组合短。
5. **库文件作位置参数**：`e-*` 前缀已声明"要用库"，库文件直接作为第一个
   位置参数（沿用原 `pl ir <repeat> <infile>` 的接口），不引入 `--lib`。

> **多库场景（2026-08-03 决策）**：`e-*` 命令只接受单个库位置参数；同时
> 使用多个库（如 RepBase + Dfam）时**不引入 `--lib` 多值**，而是对每个库
> 分别执行一次命令、再用 `spanr merge` 合并各 runlist（与 ir + trf 合并
> 同一模式，示例见 [docs/rept.md](../../docs/rept.md) e-kmer 一节）。
> 理由：把多库合并进单 FastK 表会破坏 `--keep-index` 的缓存失效逻辑
> （缓存 key 需按多文件 mtime 组合判断）、丧失每库独立调参与诊断能力；
> 代价只是基因组多扫描几次，FastK 很快，对遮蔽场景可接受。若未来出现
> "很多库 × 超大基因组"的真实需求，再考虑 `--lib` 可重复作为升级。

#### 命名规则

重复检测有 4 个组合：

* 库 + kmer → `pgr rept e-kmer <repeats> <genome>`（原 `pgr pl ir`）
* 自身 + kmer → `pgr rept s-kmer <genome>`（原 `pgr pl rept`）
* 库 + align → `pgr rept e-align <repeats> <genome>`（未来遮蔽版）
* 自身 + align → `pgr rept s-align <genome>`（未来，`scripts/pgr-repeat.sh`）

加上 `pgr rept trf`（串联重复）共 5 个命令。
命名规则：前缀 `e` = 外部库 / `s` = 自身（对象），后缀 `kmer` / `align` = 机制。
库文件是 `e-*` 的位置参数（沿用原 `pl ir <repeat> <infile>` 的接口），无需 `--lib`。
`trf` 不参与 e/s 前缀——它封装 TRF 工具、输入只有基因组，保留工具名即可。

#### 为什么选 e / s 前缀

*   **对象决定输入**：用户选命令时第一个要判断的就是"我有重复库吗？"——
    `e`（要库，需要准备库文件）/ `s`（纯基因组，不需要库）。这是最先、
    最容易确认的事；而"kmer 还是 align"是第二步的、偏实现的问题。
*   **对象是稳定维度，机制是可扩展维度**：库/自身是二元的、不会变；
    kmer/align 只是当前两种机制，未来可能加 HMM、profile 比对等——扩展
    维度放在后缀或参数里，命令名才稳定。如果前缀表达机制，每加一种机制
    就要改命令族。
*   **与 pgr 风格一致**：`pgr sd` 命令名（search / align / cluster）表达
    "做什么"，机制用 `--engine pgi|lastz`；`pgr dist` 也按对象（hv / pgi /
    seq）分命令。pgr 没有"机制进命令名"的先例——机制差异（kmer 快但粗、
    align 慢但准）靠帮助文本和文档传达就够了。
*   **字母自然**：external / self 首字母，配合单词后缀后整体可读。

> 状态：2026-08-03 已迁移。`pgr rept e-kmer / s-kmer / trf` 落地，`pl` 下
> 已移除 ir/rept/trf；`e-align` / `s-align` 随遮蔽版计划实现。

### 1.4 术语澄清：SD 序列不是真正的 repeats

`pgr sd`（BISER 移植，见 [[sd.md]]）检测的**分段重复（segmental duplications, SD）**
在重复标记语境里容易混淆，需要澄清：

*   SD 是祖先复制事件产生的**旁系同源（paralogous）共享片段**（如 T2T-CHM13 标准：
    > 1 kb 且 identity > 90%），它们虽然"序列重复出现"，但**不是转座子等真正的
    重复元件（repeats）**。
*   `pgr rept e-kmer/s-kmer/trf` 检测的是重复序列本身（转座子、rRNA 基因簇、串联重复等）；
    `pgr sd` 检测的是旁系同源片段。两者目的不同，**不要把 `pgr sd` 的输出当成
    repeat masking 的结果**。
*   实践中 SD 在比对/组装中会造成假比对（旁系同源片段会被多处匹配），因此检测出 SD 后，
    下游流程通常会把它们排除或特殊处理；但**先被屏蔽的从来不是 SD，而是真正的重复元件**——
    BISER 的输入就要求预先 soft-mask 重复序列（RepeatMasker/TRF 等），SD 恰恰是
    "屏蔽重复元件后仍剩余的高相似旁系同源片段"，BISER 找的就是它们。顺序是：
    屏蔽 repeats 在前 → 检测 SD 在后 → 排除 SD 在下游比对中，三者并不矛盾。
*   推论：若屏蔽后还要做 SD 搜索（对应 BISER 输入假设），屏蔽应**只用 `pgr rept e-kmer` +
    `pgr rept trf`**（≈ T2T-CHM13 的 TRF + RepeatMasker 预处理），**不要用 `pgr rept s-kmer`
    （自比较）**——它会把 SD 本身也当作"重复"屏蔽掉，屏蔽完 SD 搜索就找不到目标了。
    注意 `e-kmer` 需要重复库（Dfam/RepBase），无库时该组合退化为只用 `trf`：

    ```bash
    # SD 搜索前的正确屏蔽：IR + TRF
    pgr rept e-kmer genome.fa repeats.fa -o ir.json # 散在重复（需重复库）
    pgr rept trf genome.fa -o trf.json              # 串联重复
    spanr merge ir.json trf.json -o mask.json # 合并区间
    pgr fa mask genome.fa --runlist mask.json -o masked.fa
    ```

### 1.5 检测管道实现

#### 1.5.1 e-kmer / s-kmer：FastK → Profex → spanr

共享管道在 `src/libs/pl/repeat.rs`：

1.  **FastK**：
    *   `e-kmer`：跑两次——先用 `-t` 对 repeat 库建表（`-Nrepeat`），再对基因组用 `-p:repeat`
        生成相对该表的 count profile（`-Ngenome`）。
    *   `s-kmer`：只跑一次，`-p` 自比较生成基因组自身的 profile。
2.  **Profex per chr**：`pgr fa size` 得到染色体列表后，对每条染色体跑
    `Profex -z genome <sn>`，解析输出中 `start-end`（rept 还会按 `depth` 过滤，`min_depth=2`），
    写成 `<chr>:start-end` 的 `.rg` 文件（`run_profex_per_chr`）。
3.  **spanr 区间处理**（`run_repeat_spanr_pipeline`）：

    ```
    spanr cover <rg files>
        | spanr span --op fill   -n <fk>   # 填 k-mer 之间的孔
        | spanr span --op excise -n <min>  # 切掉过短的碎片
        | spanr span --op fill   -n <ff>   # 合并邻近片段
        -o <outfile>
    ```

默认参数：`kmer=17`、`fill-kmer=2`、`fill-fragment=10`；`e-kmer` 的 `min-len=300`，`s-kmer` 的 `min-len=100`。

**库表缓存**：`e-kmer` 默认每次在临时目录对重复库建 FastK 表（`repeat.ktab` +
隐藏分块 `.repeat.ktab.N`），用完即删。`--keep-index`（与
`pgr align pgi --keep-index` 同款）把整组表原子写到库文件旁
（`<库>.repeat.k<k>.ktab` + `.complete` 标记），后续运行直接
`FastK -p:<前缀>` 读缓存（验证过 `-p:` 接受路径、零复制）；库文件
mtime 变新时缓存自动失效重建。

#### 1.5.2 trf：trf → 解析 → spanr

`src/cmd_pgr/pl/trf.rs`：按染色体拆分 FASTA，逐条跑
`trf <chr>.fa <match> <mismatch> <delta> <pm> <pi> <minscore> <max_period> -d -h -ngs`，
用 `parse_trf_output`（`src/libs/pl/repeat.rs`）解析 `.dat`（少于 15 列的短行跳过），
再 `spanr cover` 合并输出。默认参数对应 TRF 常用设置（match=2、mismatch=7、delta=7、
pm=80、pi=10、minscore=50、max_period=2000）。

### 1.6 临时文件与 FastK 库文件清理

用户曾担心 FastK 会在工作目录生成一批库文件（`*.ktab.*`）。该问题已由
`src/libs/pl/ctx.rs` 的 `PipelineCtx` 内建解决：

*   管道启动时创建 `tempfile::TempDir`（前缀 `pgr_rept_e_` / `pgr_rept_s_` / `pgr_rept_trf_`）；
*   `enter()` 把 CWD 切进 tempdir，此后 FastK 的 `genome.ktab.*` / `repeat.ktab.*`、
    Profex 的 `prof.*.txt/.rg`、trf 的 `.dat` 全部落在 tempdir 内；
*   ctx drop 时 TempDir 自动删除，`CwdGuard` 保证出错时 CWD 也能恢复。

实测（2026-08-03，MG1655）：跑完 `pgr rept s-kmer` 后 `/tmp` 无新增残留，tempdir 也不存在；
FastK `-P` 默认丢到 /tmp 的排序块由 FastK 自身清理。因此**无需**在代码里额外做删除动作。

### 1.7 实测记录

| 命令 | 基因组 | 耗时 | 区间数 | 覆盖 (bp / %) | 备注 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `s-kmer` | MG1655 | ~0.35s | ~150 | — | 含 rRNA 3941442-3946950 等 |
| `trf` | MG1655 | ~1.6s | 84 | 18,768 / 0.40% | 串联重复，与 e-kmer 零重叠 |
| `e-kmer` (tncentral) | MG1655 | ~0.8s | 48 | 56,969 / 1.23% | 原核 IS 专库，最敏感 |
| `e-kmer` (repbase) | MG1655 | ~4s | 38 | 42,763 / 0.92% | |
| `e-kmer` (dfam) | MG1655 | ~2s | 39 | 42,386 / 0.91% | |

soft-mask 衔接验证：`aaaaaatgcgcggtcagaa` 等区间正确转为小写。

三库 + trf 对 MG1655 的完整执行脚本、结果表与注意事项（FastK 并行
SIGSEGV 等）见 [docs/rept.md](../../docs/rept.md) 的 "Example run:
E. coli MG1655" 一节（用户文档）。

## 2. 遮蔽版计划（Dfam 全库）

### 2.1 决策

**只实现遮蔽版重复标记**：把 Dfam consensus 全库与目标基因组做比对，输出
runlist 区间交给 `pgr fa mask` 做遮蔽。不做 family/class 注释、不做
`.out/.tbl` 报告、不构建自有重复库。

理由：

*   RepeatMasker 的价值一半在库（Dfam 几十年积累），一半在注释后处理（二十年
    打磨）——两者都不追：**短期无法达到 Dfam 的库积累，注释后处理是泥潭**。
*   直接用 Dfam consensus FASTA 当候选集，用 pgr **原生比对命令**
    `pgr align pgi`（或 `pgr align lastz`）找区间——`pgr sd search --engine pgi`
    已验证 pgi 用于重复区间检测的模式；遮蔽管道 `fa mask` 已存在。
*   遮蔽版只回答一个问题："基因组上哪些区间是重复的"，质量取决于比对敏感度
    （验证见 §2.5）。

### 2.2 方案

#### 核心思路：Dfam 全库 + 一套通用比对

RepeatMasker 完整流程 = 搜索（库 vs 基因组）+ 注释后处理（碎片整合、边界精修、
family/class、K2P、报告，详见附录 A.5）。我们只做**搜索侧的遮蔽**。

**简化点：不按物种分类。** RepeatMasker 的物种逻辑只影响两件事——库筛选
（FamDB 按 `-species` 从 Dfam 取子集）和搜索配方选择（按 primates/rodentia
跑专项 stage），两者都是为注释精度服务的。对纯遮蔽目标，直接**把 Dfam
全库 consensus 当候选集**，一套通用搜索参数跑完：

*   不做 Taxonomy / FamDB / 物种库解析；
*   不按谱系选配方，RepeatMasker 的 17 个配方（附录 A.3）收敛为 1–2 套通用
    参数（可参考 `general_search_parameters` 的量级，放宽 min-len/identity）；
*   这正好对应 RepeatMasker 自己的 `-lib` 模式——官方支持"用户给全库 +
    `general_search_parameters` 一套流程"的简化路径（附录 A.2）。

**代价**：全库比对会让近缘家族的 consensus 互相 hit（cross-family），造成
一定过度遮蔽；物种特异配方（如灵长类年轻 Alu 的"先切除再补搜"）的灵敏度
也会丢失。对"宁可多遮不漏"的遮蔽目标可接受，但要在验证里量化（§2.5）。

遮蔽版明确不做 RepeatMasker 的注释后处理，边界见 §3.2。

#### 现状保留为对照与兜底

`ir + trf + fa mask` 继续可用：方案落地后与之对比覆盖区间；低复杂度缺口
（polyA 等，§2.4）由 `rept trf` 兜底。

### 2.3 实现步骤

新命令 `pgr rept e-align <repeats> <genome>`，输出 runlist JSON（与
e-kmer/s-kmer/trf 同构，可直接喂 `pgr fa mask`）。管道整体在
`PipelineCtx` 临时目录内执行；`--keep-index` 沿用 e-kmer 的缓存约定，
把 `.pgi` 索引缓存到输入文件旁（见 §1.6）。

#### 2.3.1 处理管道（算法细节）

```
genome.fa + repeats.fa
  │ 1. 建索引（.pgi，--keep-index 可缓存）
  ▼
pgr align pgi <genome> <repeats> -o hits.psl
  │ 2. 比对：ref = genome（target），query = repeats（候选）
  ▼
PSL 过滤（Rust 内）：identity ≥ min-identity 且 block ≥ min-len
  │ 3. 按 matches / (matches+mismatches+ins+del) 计算每条 alignment 的 identity
  ▼
target 侧 .rg（基因组坐标，1-based inclusive）
  │ 4. pgr psl to-range --target-coords 语义
  ▼
spanr cover → spanr span --op excise -n <min> → spanr span --op fill -n <ff>
  │ 5. 合并重叠、切短碎片、填邻近孔（同 e-kmer 的 spanr 管道）
  ▼
runlist JSON → pgr fa mask
```

要点：

1.  **比对方向**：ref = 基因组、query = 库。理由：遮蔽要的是基因组坐标，
    PSL 的 target 侧正是 ref；`pgr psl to-range --target-coords` 直接给出
    基因组区间。query 索引按 pgi 约定 memory-map（库约 50 MB，无压力）。
    注意 `pgr align pgi` 的 PSL 输出是"每链一个 block"，不是逐 hit 列表。
2.  **身份过滤**：转座子拷贝与 consensus 的 identity 通常 70–90%（`sd search`
    的 >90% 是给分段重复调的，不适用）。identity 在 Rust 内从 PSL 记录的
    matches/mismatches/ins/del 计算，不依赖外部工具（复用 `pgr::libs::psl`
    的记录解析）。`--min-identity` 初始 0.70、`--min-len` 初始 50，均为
    可调参数。
3.  **区间合并**：`spanr cover` 合并重叠块 → `spanr span --op excise -n
    <min-len>` 切掉过短碎片 → `spanr span --op fill -n <fill-fragment>`
    合并邻近片段（fill-fragment 默认 10，与 e-kmer 一致）。不做
    RepeatMasker 的碎片整合/边界精修（§3.2）。
4.  **库**：Dfam consensus FASTA 全库直接作 query（不做物种筛选，§2.2）；
    下载与准备见 [docs/rept.md](../../docs/rept.md)。多库沿用 §1.3 决策
    （每库跑一次 + `spanr merge`）。

#### 2.3.2 参数初始值与调参空间

`pgr align pgi` 参数与 RepeatMasker 配方（附录 A.3）不是一一对应，但方向
可比。初始值如下，全部进 CLI 透传，实施后按 §2.5 验证调整：

| pgi 参数 | 含义 | 初始值 | 方向 | 对应 RM 配方考量 |
| :--- | :--- | :--- | :--- | :--- |
| `-k/--smer/--window` | syncmer 种子参数 | 40 / 8 / 5（pgi 默认） | 待验证 | 种子敏感度（RM 靠 minmatch 档） |
| `-f/--freq` | 高频 k-mer 跳过阈值 | 10（默认）→ 放宽 100 | **放宽** | 高拷贝家族（Alu 百万拷贝）种子会超频被丢 |
| `-c/--min-span` | 链最小种子跨度 | 85（默认）→ 50 | 放宽 | RM 短重复配方 minmatch 低 |
| `-s/--max-gap` | 链内种子最大 gap | 1000（默认） | 保持 | |
| `--band` | 对角带半宽 | 128（默认） | 保持 | RM bandwidth 14（RMBlast 单位不同） |
| `--merge-gap` | 共线链合并 gap | 5000（默认） | 保持 | |
| `--min-shared` | 共享种子最小长度 | k（greedy 默认）→ 16 | 放宽 | 高分歧拷贝共享种子少 |
| `--workflow` | greedy / tube | greedy | 对比 | FastGA tube 的 plen floor=12 更敏感 |
| `-p` | 线程数 | 8 | — | |

过滤参数（CLI 默认，实施后标定）：

| 参数 | 含义 | 初始值 |
| :--- | :--- | :--- |
| `--min-identity` | PSL identity 下限 | 0.70 |
| `--min-len` | 最小 block/碎片长度（bp） | 50 |
| `--fill-fragment` | 邻近片段合并孔宽（bp） | 10 |

#### 2.3.3 代码结构

遵循分层原则：`src/cmd_pgr/rept/e_align.rs` 只做 clap 解析与参数转换；
管道逻辑放 `src/libs/pl/repeat.rs`（与 e-kmer 同文件，新增
`run_align_repeat_pipeline`）：

```rust
pub struct AlignRepeatOpts {
    pgr: PathBuf,            // 自身二进制（子进程调用 align pgi）
    abs_repeat: String,      // 库 FASTA 绝对路径
    abs_infile: String,      // 基因组 FASTA 绝对路径
    abs_outfile: String,
    keep_index: bool,
    // pgi 透传参数
    kmer: usize, smer: usize, window: usize,
    freq: usize, min_span: usize, max_gap: usize, band: usize,
    merge_gap: usize, min_shared: usize, workflow: String,
    // 过滤参数
    min_identity: f64, min_len: usize, fill_fragment: usize,
    parallel: usize,
}
```

子进程调用：`pgr align pgi <genome> <repeats> -o hits.psl`，之后 Rust 读
hits.psl 过滤 → 写 target .rg → `spanr` 合并 → runlist。PSL 记录解析复用
`pgr::libs::psl`。

#### 2.3.4 实施顺序

1.  骨架：`e_align.rs` + `run_align_repeat_pipeline`，参数按 §2.3.2 初始值，
    跑通 MG1655 + 现有三个库（冒烟，确认管道与坐标方向）。
2.  真核验证与调参（§2.5）：拟南芥/玉米 + Dfam 全库，扫描关键参数
    （`-f`、`--min-shared`、`--workflow`、`--min-identity`），确定默认值。
3.  收尾：`docs/rept.md` 的 e-align 一节、集成测试（`tests/cli_rept.rs`）。

工作量比完整 RepeatMasker 小一个数量级。遮蔽版需要的能力映射见附录 A.7。

### 2.4 关键风险

*   **比对敏感度**：k-mer（k=17）对高分歧拷贝会漏；pgi 的 syncmer 种子对
    70% identity 的拷贝同样不轻松。这是决定遮蔽质量的核心，必须实测。
*   **全库 cross-family 假阳性**：不做物种筛选后，保守的转座子区域可能让
    基因/其他序列被误遮蔽（over-masking）。遮蔽场景可接受，但需在验证中
    对比"全库 vs 物种库"的遮蔽量差异。
*   **低复杂度缺口**：RepeatMasker 默认屏蔽 low complexity（polyA、卫星、
    homopolymer）。现有 `e-kmer` 只管库内散在重复，`trf` 覆盖串联重复，polyA 这类
    不一定被覆盖。这是遮蔽质量上更实际的差距，与用 k-mer 还是比对无关。
*   **验证基准**：E. coli 几乎无转座子，无参考价值。需用拟南芥/玉米等
    转座子丰富基因组，与 RepeatMasker 的 masked 输出对比 recall。

### 2.5 验证实验（实施前调参）

方向已定（§2.1），验证的目的是评估比对敏感度、确定过滤参数：

1.  取转座子丰富基因组（拟南芥或玉米），先跑 RepeatMasker 得到 masked
    参考 runlist（沿用 §RepeatMasker (reference) 的 gff→runlist 流程）；
2.  Dfam 全库作 query，跑 `pgr rept e-align` 骨架；**参数扫描**：
    `-f ∈ {10, 50, 100}`、`--min-shared ∈ {12, 16, 40}`、
    `--workflow ∈ {greedy, tube}`、`--min-identity ∈ {0.60, 0.70, 0.80}`，
    其余保持默认；每次记录耗时、hits 数、覆盖区间；
3.  对比数据：
    *   覆盖区间 vs 现有 `e-kmer` 的差异；
    *   与 RepeatMasker masked 输出的 recall（`spanr statop`，同 §1.8 流程）；
    *   时间与内存；
    *   （可选）全库 vs 按物种取库的遮蔽量差异，评估 over-masking 代价。
4.  依据 recall / over-masking / 耗时确定最终默认值，写回 §2.3.2 与
    `pgr rept e-align` 帮助文本。

> 命令命名见本文 §1.3（`pgr rept` 组的 2×2 组合形式）。

## 3. 与 RepeatMasker 的差异

### 3.1 现状命令的局限

*   **无分类注释**：不输出 family/class 标签，只给区间。
*   **k-mer 敏感度**：依赖与库共享的精确 k-mer，分化较远的拷贝会漏检或碎成小片段；
    fill 步骤只能桥接短孔，无法恢复长距离分化的拷贝。
*   **输出是区间而非序列**：mask 后的序列需另行用 `pgr fa mask` 生成。
*   **依赖外部工具**：FastK / Profex / spanr（trf 还需 trf）。

### 3.2 遮蔽版的边界

遮蔽版明确不做 RepeatMasker 的注释后处理：

*   碎片整合（`cycleReJoin`）：遮蔽不在乎拼回完整元件，只在乎别漏区域；
*   边界精修：区间级覆盖即可；
*   family/class 注释、K2P %div、`.out/.tbl` 报告：全部不做。

## 4. 待办

*   `e-kmer` 需要用户自备重复库（Dfam/RepBase/TnCentral，下载与准备见
    [docs/rept.md](../../docs/rept.md)），本机缺库，端到端测试待补。
*   若未来要接近 RepeatMasker 能力，可考虑对检测出的区间补一步 family 注释（如对区间
    重跑库比对），但目前无此需求，不做推测性设计。

## 附录 A：RepeatMasker 源码梳理

> 2026-08-03 依据仓库内 `RepeatMasker/` 目录源码（open-4.2.4，Arian Smit &
> Robert Hubley）梳理，作为 §2 方案的源码证据。

### A.1 概览

RepeatMasker 筛查 DNA 序列中的**散在重复**（interspersed repeats）与**低复杂度
序列**，输出详细注释（family/class、%div/%del/%ins）以及 masked 序列。

*   **实现语言**：Perl（主驱动 + 模块化库），比对本体交给外部搜索引擎
    （RMBlast / crossmatch / NHMMER / ABBLAST），自身只做流程编排与后处理。
*   **主要文件**：
    *   `RepeatMasker`（267 KB）：主驱动——分片、多阶段搜索、汇总 `.cat`。
    *   `ProcessRepeats`（364 KB）：后处理——碎片整合、注释、距离估计、输出。
    *   `DupMasker`：片段重复（segmental duplication）检测。
    *   `RepeatProteinMask`：转座子编码蛋白的蛋白级比对。
    *   `DateRepeats`：按 %div 估计重复插入时间。
    *   `TRF.pm` / `TRFResult.pm`：串联重复（tandem repeat）封装。
    *   `Libraries/RepeatAnnotationData.pm`：family/class 注释数据。
*   **搜索引擎抽象**：`SearchEngineI.pm` 定义接口，
    `CrossmatchSearchEngine` / `WUBlastSearchEngine` / `HMMERSearchEngine` /
    `NCBIBlastSearchEngine`（RMBlast = 特制 blastn）四个实现。
*   **库来源**：Dfam（`DFAM.pm`）/ RepBase（`RepbaseEMBL.pm`）经 FamDB 按物种
    生成 consensus/HMM 库；也可用 `-lib` 直接给 FASTA。库 header 解析
    family/class（如 `>XXX#LTR/ERV`）。

### A.2 主驱动 RepeatMasker

**顶层流程**：

1. **参数解析**：Getopt::Long；配置优先级 CLI > 环境变量 > 配置文件
   （`RepeatMaskerConfig.pm`）。
2. **分片**：`fragmentSize = 60000`、`overlapLen = 2000`（`-frag` 可覆盖，
   且必须 ≥ 2×overlap）。分片边界会产生 edge effect，后处理专门处理。
3. **库解析**：FamDB 按物种取库，或 `-lib` 自定义；`refineableHash.dat`
   记录可精修的元素 id。
4. **并行**：`-parallel N` fork 子进程，每个 batch 独立跑搜索。
5. **每 batch 搜索**：`runSearchStages`（blast 类引擎）或
   `runHMMERSearchStages`（NHMMER）。
6. **汇总**：各 batch 的 `.cat` 拼接成总 `.cat`（含版本/引擎/库 header +
   `## RAW Annotations:` 段）。
7. **后处理**：调用 `ProcessRepeats`（`-nopost` 跳过，供手动分步跑）。

**多阶段搜索（runSearchStages）**：搜索不是一次 blastn，而是**一串配方化搜索 +
迭代切除**：

1. **TRF 阶段**：`runTRFStage` 找 perfect simple repeats（`-nolow` 跳过）。
2. **用户库/物种库**：用 `general_search_parameters` 配方（`-cutoff` 改
   minscore），输出 `*.tmp.custom`。**注意：`-lib` 自定义库时根本不走物种
   分类，直接用这一套通用配方**——这就是正文 §2.2"全库 + 统一参数"简化
   所对应的官方路径。
3. **物种特异阶段**：按分类学（primates/rodentia…）跑对应配方——如灵长类
   先 `cut_young_sines_in_primates` 切 Alu，再 `mask_*` 遮蔽其余 SINE；
   每阶段后 `postProcessSearch` 把已找到的 hit **从当前序列中切除**
   （`excise`），下一阶段只搜剩余序列。

`runHMMERSearchStages` 的注释给出了完整阶段顺序（哺乳类）：sinecutlib →
shortcutlib → cutlib →（人类再扫一次 sinecutlib）→ shortlib → longlib →
mirs.lib → mir.lib → retro.lib → l1.lib → simple.lib → at.lib。每个阶段都
是"搜索 → masklevel 过滤 → 结果筛选 → 切除 → 写 cat"。

**"两轮搜索"的真实含义**：早期版本是两遍 blastn（高严格度找 anchor + 低严格度
补漏），现版本是"多阶段配方 + 切除迭代"——每次切除后剩余序列变少，后续阶段
只处理未覆盖区域。

### A.3 搜索配方（getSearchRecipes）

共 17 个配方，每个定义一套完整搜索参数：

| 配方 | 用途要点 | minscore | minmatch | matrix | bandwidth | masklevel | excise |
| :--- | :--- | ---: | :--- | :--- | ---: | ---: | :--- |
| `perfect_simple_repeats` | 完美简单重复（TRF 替代） | 180 | [8,9,10,11] | simple1 | 1 | 1 | 1 |
| `general_search_parameters` | 通用库搜索（默认路径） | 225 | [8,9,11,13] | 20p##g | 14 | 90 | 0 |
| `cut_young_sines_in_primates` | 灵长类年轻 SINE 切除 | 1200 | [7,8,10,12] | 14p##g | 20 | 1 | 1 |
| `mask_young_sines_in_primates` | 灵长类年轻 SINE 遮蔽 | 1500 | [7,8,10,12] | 14p##g | 20 | 80 | 0 |
| `mask_sines_in_primates` | 灵长类其余 SINE | 225 | [7,8,8,9] | 20p##g | 14 | 80 | 0 |
| `mask_sines_in_non_primate_mammals` | 非灵长哺乳类 SINE | 225 | [6,7,8,10] | 18p##g | 14 | 1 | 1 |
| `general_full_length_repeats` | 全长元件（大带宽跨大 gap） | 300 | [9,10,11,13] | 18p##g | 40 | 1 | 1 |
| `complete_3end_of_young_line1s` | 年轻 LINE1 的 3' 端 | 300 | [9,10,11,13] | 18p##g | 40 | 90 | 1 |
| `older_ALUs_in_primates` | 古老 Alu | 800 | [7,8,10,12] | 14p##g | 20 | 1 | 0 |
| `more_ALUs_in_primates` | 其余 Alu | 400 | [7,8,9,11] | 18p##g | 14 | 10 | 0 |
| `short_repeats_and_satellites_rodents` | 啮齿类短重复/卫星 | 210 | [7,8,9,10] | 25p##g | 14 | 90 | 0 |
| `short_repeats_and_satellites` | 通用短重复/卫星 | 225 | [7,8,10,12] | 20p##g | 14 | 90 | 0 |
| `long_interspersed_repeats` / `ancient_repeats` / `tough_ancient_repeats` / `retroviruses` / `tough_line1s_in_eutheria` | 长散在/古老/反转录病毒等分级搜索 | 各异 | 各异 | 18p##g | 各异 | 各异 | 各异 |
| `simple_repeats_again` / `simple_repeats_flanking` | 低复杂度补搜与侧翼 | 各异 | 各异 | simple 系矩阵 | 小 | 高 | 0 |

`minmatch` 是 4 个值的数组，按严格度档位选（`selectParameter`，受 `-s`/
`-q`/`-qq` 速度档影响）；`matrix` 形如 `20p##g.matrix`（`#` 数量按 GC 含量
调整，见 `runStage` 的矩阵选择逻辑）；`masklevel` 配合 RMBlast 的
complexity-adjusted scoring；`raw=1` 时用 basic scoring。

### A.4 搜索引擎与 RMBlast 参数

`search()` 把配方翻译成引擎参数：

*   `NCBIBlastSearchEngine` → `rmblastn`：
    *   `-word_size 14`（默认；`minmatch` 经 `-word_size` 传入）；
    *   `-dust no`（低复杂度过滤交给 masklevel，而不是 blast 的 dust）；
    *   `-min_raw_gapped_score` = minscore；
    *   `-xdrop_ungap / -xdrop_gap_final / -xdrop_gap` 由 minscore 推导
        （refine 阶段带宽传 `-1` 放宽）；
    *   `-matrix`（配方矩阵）、`-gapopen / -gapextend`（由 gap init/ext 换算）；
    *   默认 complexity-adjusted score mode（`raw=0`），`raw=1` 时 basic mode。
*   **失败重试**：bandwidth > 14 时缩到 14，==4 时缩到 1，minmatch < 10 时
    增大，仍失败则放弃该 batch。
*   crossmatch/WUBlast 用各自的矩阵目录（`Matrices/crossmatch/`、
    `Matrices/wublast/aa/`、`Matrices/ncbi/nt/`）。

### A.5 后处理 ProcessRepeats

输入 `.cat`（raw annotations），输出 `.out` / `.tbl` / `.masked` /
`.align` / GFF / HTML。核心是 `processSequence` 的多 cycle 流水线：

*   **cycle 1**：`cycleReJoin` —— 把被打断的转座子碎片按打分连成链
    （left/right linked hit 结构）。
*   **cycle 2**：移除 edge effect 注释（分片边界产物）、移除 masklevel 违规、
    卫星重复重命名、构建 DNA 转座子等价结构（`%chainBeg/%chainEnd`）。
*   **cycle 3–5**：进一步去碎片化（CYCLE5 维护最近 21 个 join 打分做决策）、
    LTR/LINE 配对（`scoreLINEPair`：按 query/consensus 重叠与名字兼容性打分）、
    等价类间传播。
*   **边界精修**：对 `refineableHash` 中的元素，用更宽松带宽的搜索
    （`-xdrop_gap_final` 放宽）重对齐边界，`replaceRMFragmentChainWithRefinement`
    替换原链。
*   **统计**：`calcKimuraDivergence`（K2P，CpG 位点校正）给出
    %div/%del/%ins；simple/low_complexity 不记 divergence。
*   **输出**：`generateOutput`（.out 标准表）、`generateTableOutput`（.tbl
    汇总）、masked 序列（默认 N，`-xsmall` 小写）、`-a` 时 source alignments。

### A.6 辅助程序

*   `DupMasker`：片段重复检测（独立于库搜索的另一条路径）。
*   `RepeatProteinMask`：用 RepeatPeps 蛋白库跑蛋白比对，找转座子编码蛋白。
*   `DateRepeats`：根据 %div 估计插入时间（转座子年代学）。
*   `TRF.pm`：封装 TRF 4.09+，`runTRFStage` 用它替代 simple consensus 搜索。

### A.7 与 pgr 的对应关系

**遮蔽版需要的能力**（正文 §2.3）：

| 步骤 | pgr 现状 | 判断 |
| :--- | :--- | :--- |
| 库-基因组比对 | `pgr align pgi`（原生）或 `pgr align lastz` | 高可行，预计比 RMBlast 快一个数量级 |
| 区间合并/覆盖 | `spanr cover / merge / fill` | ✅ 已有 |
| 输出遮蔽 | `pgr fa mask --runlist` | ✅ 已有 |
| 低复杂度兜底 | `pgr rept trf` | 已有（缺口见正文 §2.4） |

**遮蔽版明确不做**（对应 RepeatMasker 的注释后处理，A.5）：

| 步骤 | 说明 |
| :--- | :--- |
| 碎片整合（cycleReJoin） | 遮蔽不在乎拼回完整元件，只在乎别漏区域 |
| 边界精修 | 区间级覆盖即可 |
| family/class 注释 | 不做 |
| K2P %div | 不做 |
| `.out/.tbl` 报告 | 不做 |

**其余启示**：

*   搜索侧 = 配方化 blastn + 迭代切除，用 pgi/lastz 替换理论上可行且更快；
*   真正的工程量和价值在 `ProcessRepeats` 的 6+ cycle 后处理——这正是遮蔽版
    不做的部分；
*   低复杂度（Simple_repeat / Low_complexity）在 RepeatMasker 里由
    TRF + simple.lib 覆盖，对应 pgr 的 `rept trf`，是遮蔽质量的实际差距。
