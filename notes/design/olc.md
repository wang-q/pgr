# OLC 组装设计：多 k unitig 层的 overlap → layout → consensus

> 状态：**已实现（2026-08-12，M1–M4 全部落地）**。命令：
> `pgr asm ovlp` / `layout` / `cns` / `olc`；库：`libs/olc/`
> （overlap / layout / consensus）。合成基因组端到端验证通过
> （30× 无错 reads → contigs 全部为基因组精确子串，最长覆盖 97.5%；
> 低覆盖 6× 出现重复区经典错装，符合预期）。
> **真实数据验证（Lambda，2026-08-12）**：见 §12——抓出并修复左向延伸
> 坐标下溢 bug；多 k 冗余 2.4×（v1 待消减）；参考菌株差异会误报
> "嵌合"，需 reads 侧验证。
> 需求来源：用户裁定"不对 reads 做 OLC，把不同 k 各自生成的 unitigs 当
> 伪 reads，在 unitig 层做 OLC 拼接"（`todo.md` §3、`references/canu.md` §8）。
> 参考源码：`canu-2.3/`（bogart + utgcns）、`wgs-8.3rc2/`（AS_BAT + AS_CNS），
> 分析笔记见 `references/canu.md` / `references/celera.md`。

## 1. 目标与范围

为 pgr 增加一条完整的 OLC 组装流程：

```
reads ──(多 k 各生成 unitigs)──> 伪 reads ──(overlap)──> PAF
  ──(layout)──> 路径 ──(consensus)──> contigs.fasta
```

**核心裁定（用户，2026-08-12）**：

* 不做 reads 级 OLC——reads 先用 DBG（`pgr asm unitig`）压缩成 unitigs；
* 多 k（默认 21/51/81）各生成一套 unitigs，合并后当伪 reads；
* unitig 语义 = 最大无分支路径（bcalm graph3 移植，无气泡），OLC 只处理
  unitig 间 overlap，不引入"平行路径选哪条"的启发式；
* 气泡/孤儿合并不处理（既有裁定）。

**成功标准**：合成基因组 reads → `pgr asm olc` 输出 contigs 覆盖完整基因组
（identity 100%，因 overlap 精确）；Lambda 真实 reads 冒烟测试输出合理
（contig 数 / N50 与 tadpole contig 同量级）；1701+ 测试全绿、fmt/clippy 干净。

## 2. 为什么 unitig 层 OLC 可行

* **数据量**：unitigs 数远小于 reads（宏基因组下少 1~2 个数量级），
  all-pairs overlap 成本从不可行变可行；
* **规避气泡**：unitig 无分支，重叠只发生在 unitig 间，天然无
  "平行路径投票"问题；
* **多 k 互补**：小 k 连通性好（低覆盖区、重复边界），大 k 特异性强
  （区分重复/菌株）；不同 k 的 unitigs 共享精确子串，overlap 证据天然存在；
* **精确性**：unitig 序列来自 DBG 固实 k-mer，无测序错误，overlap 可做
  全精确验证（不做 Myers/edlib 扩展，参考项目的高噪声机制不需要）。

## 3. 参考实现对照（详见 references/canu.md、celera.md）

| 环节 | Celera 8.3rc2 | Canu 2.3 | pgr 借鉴 |
|---|---|---|---|
| overlap | OlapFromSeedsOVL（k=9 seed + banded DP） | MHAP / overlapInCore（k-mer 哈希 + Myers） | seed→verify 骨架同源；pgr 是精确版（无错误模型） |
| layout | bogart（AS_BAT）：BestOverlapGraph + greedy 双向延伸 + 覆盖度证据 repeat split | bogart（同源）：互惠 best edge 种子 + 单向延伸 + markRepeatReads | greedy best-edge + 互惠种子；repeat 用"双定位"思想（v1） |
| consensus | AS_CNS：MA 列投票（BaseCallMajority） | utgcns：template stitch + edlib 重比对 + POA-DAG bestPath | 精确 overlap 下缝合即可；列投票 = 将来的鲁棒化方向 |

## 4. 管线设计（四阶段）

### S0 伪 reads 生成（复用，不新写）

每个 k 跑 `pgr asm unitig`（`libs/asm/assemble.rs::assemble_unitigs`，
bcalm graph3 压缩语义，默认 k=31 / solid ≥3）。产出物：unitig FASTA。

**命名**：`asm unitig` 的输出名恒为 `unitig_<id>`，多 k 合并必然撞名。
OLC 阶段统一重命名：`<tag>:<name>`，tag 默认取输入文件 stem
（仅保留 `[A-Za-z0-9_.-]`，空则用文件序号）——确定性且可回溯到 k。
（`pgr asm olc` 驱动器内部直接用 `k<k值>` 作 tag。）

### S1 overlap 检测（新 `libs/olc/overlap.rs`）

**算法**：seed → verify（与 `libs/map.rs` MapIndex 同构，精确版）：

1. 建索引：对全部 unitigs 做 canonical k-mer 索引
   （`MapIndex` 形态：`keys: Vec<u8>` 打包 FastK 字节 + `payloads:
   Vec<u64>` 存 `(cid<<32)|pos`，`radix_sort_bytes` 排序），seed k 默认 17
   （`--overlap-k`，≤ 目标最小 unitig 长度，越界则自动降到 min(k, len)）；
2. 候选：对每条 unitig q 查**边界 k-mer**——5' 端窗口 (0..k) 与
   3' 端窗口 (n-k..n)，各查正链与反互补（canonical 索引丢失方向，
   验证时双侧都试）；命中 → (cid, tpos)；
3. 验证：从 seed 对齐处向两端逐碱基扩展，得**最大精确 overlap**
   （含 seed 的 q∩t 最长精确段，长度 L ≥ k）；
4. 分类（按 q/t 覆盖关系）：
   * `dovetail`：q 5'/3' 端与 t 3'/5' 端重叠（两端各留出 >0 的非重叠段），
     或 q 完全包含于 t（contain，长度 ≥ L）；
   * `contain`：q ⊂ t 或 t ⊂ q——不参与延伸，留作共识覆盖证据；
5. 输出 PAF（复用 `libs/paf/record.rs` 12 列 + `ov:A:D|C` tag），
   去重（同一对多 seed 命中取最长）、排除自身（q==t 及回文 rc）。

**并行**：rayon 按 unitig 并行查询；索引构建与 `asm map` 同路径。

### S2 layout（新 `libs/olc/layout.rs`）

**算法**：bogart 风格 greedy best-edge 路径延伸（简化为无 mate、无气泡）：

1. 只取 dovetail overlap 建**有向图**：node = unitig，edge = 一端 → 另一端
   （q 3' 端 → t 5' 端 / q 5' 端 → t 3' 端，方向由 overlap 坐标推出）；
2. 每 node 两端各选 best edge（最长 overlap L，平局按 (target, L) 字典序
   保证确定性）；**互惠**要求：seed 的 3' best 必须同时是对方相应端的
   best（Canu 互惠种子思想，防错装）；
3. greedy：按 unitig 长度降序取未放置 seed，沿 3' best 单向延伸，
   目标已放置 / 无 edge / 目标被标记 repeat 即停；反向延伸同理
   （通过 rc 复用同一逻辑）；
4. **repeat 标记（v0 简版）**：某端的 top2 edge 长度 ≥ 0.9×best 且指向
   不同 node → 该端标记 repeat，禁止从它延伸（Canu `markRepeatReads`
   的"双定位"思想的单元化近似；覆盖度证据版本留 v1）；
5. 输出 layout TSV（每 contig 一行一步）：
   `contig_id step unitig_name strand q_start q_end overlap_len`，
   `strand` 为 unitig 在 contig 中的方向，`q_start/q_end` 是其在
   contig 坐标系中的 0-based 区间。

### S3 consensus（新 `libs/olc/consensus.rs`）

overlap 全精确 ⇒ consensus = 沿 layout **精确缝合**：

1. 按路径顺序取每步 unitig 的方向片段，与上一步 overlap 部分对齐
   （坐标由 layout 记录），追加非重叠后缀；
2. 输出 FASTA：`>contig_<id>,len=...,cov=...`（cov = 路径上 unitigs
   平均覆盖深度，近似 = 参与步数），70 列换行，与 `asm contig` 输出风格一致；
3. `--min-contig-len` 过滤短 contig。

**列投票留 v1**：若未来引入错配 overlap 或真实数据暴露 junction 不一致，
再加 AS_CNS `BaseCallMajority` 式逐列投票 + min-coverage 修剪
（Canu `consensusNoSplit` 语义），复用 `asm map` + `sam to-rg` +
`rg coverage` 的回放设施（`references/canu.md` §8.3 已论证）。

## 5. 命令设计

新增 `pgr asm` 三个叶子命令 + 一个驱动器（四层：`libs/olc/*` 管逻辑，
`cmd_pgr/asm/*` 薄壳）：

| 命令 | 输入 → 输出 | 逻辑 |
|---|---|---|
| `pgr asm ovlp` | unitig FASTA(s) → PAF | `libs/olc/overlap.rs` |
| `pgr asm layout` | PAF + unitig FASTA → layout TSV | `libs/olc/layout.rs` |
| `pgr asm cns` | layout TSV + unitig FASTA → contigs FASTA | `libs/olc/consensus.rs` |
| `pgr asm olc` | reads → contigs FASTA（驱动器） | 内部组合 S0–S3，阶段间走内存 |

`pgr asm olc` 参数：

```text
pgr asm olc <infiles>... -o contigs.fa \
    --kmer 21,51,81          # 逗号分隔，默认 21,51,81
    --min-count-seed 3       # 透传 asm unitig
    --overlap-k 17           # S1 seed k
    --min-overlap 34         # 最短接受的 overlap 长度
    --min-contig-len 500     # 输出过滤
    --keep-dir DIR           # 调试：落地中间文件（unitigs/ovlp/layout）
```

阶段命令天然可独立测试与组合（也支持用户自己跑
`asm unitig` → `asm ovlp` → `asm layout` → `asm cns` 的管道形态）。

## 6. 数据结构与格式

### overlap（PAF）

* 12 列标准 PAF（q 名 = `<tag>:<name>`，q 长 = unitig 长），
  `matches = block_length = L`，`mapq = 255`；
* tag：`ov:A:D`（dovetail）/ `ov:A:C`（contain）——不做 `cg:Z`（无错配无 CIGAR）。

### layout（TSV，无表头）

```text
contig_1<TAB>0<TAB>k21:unitig_5<TAB>+<TAB>0<TAB>2410<TAB>0
contig_1<TAB>1<TAB>k51:unitig_12<TAB>+<TAB>2410<TAB>4370<TAB>180
```

`q_start/q_end` 为该 unitig 在 contig 中的区间；`overlap_len` = 与上一步的
overlap（第 0 步恒 0）。同 contig 内区间连续（`q_end[i] == q_start[i+1]`）。

## 7. 现有基础设施复用映射

| 环节 | 复用 | 用途 |
|---|---|---|
| S0 | `libs/asm/assemble.rs::assemble_unitigs` | unitig 生成（命令层 `asm unitig` 已包装） |
| S1 | `libs/map.rs`（MapIndex 形态 + `canonical_keys` + radix） | canonical k-mer 种子索引 |
| S1 | `libs/kmer/key.rs::Kmer` | 边界 k-mer 编解码 / rc / canonical |
| S1 | `libs/nt::rev_comp` | 方向验证 |
| S1 | `libs/paf/record.rs` | PAF 写出 |
| S1 | `libs/ds/radix_sort.rs::radix_sort_bytes` | 索引排序 |
| S2 | `libs/ds/dsu.rs`（仅若需要连通分量） | 布局分组（v0 可不用） |
| S3 | `libs/fmt/seq.rs::SeqReader` | unitig FASTA 读取 |
| 全部 | `libs/io.rs` reader/writer、`cmd_pgr/args.rs` 标准参数 | I/O 与 CLI 一致性 |
| 驱动 | `libs/asm/assemble.rs` + 上述各 libs | 内存组合，无中间文件 |

**不引入新依赖**（AGENTS.md 硬性要求）；k-mer 表示统一用 FastK 字节键
（`design/kmer.md` §12 的唯一表示），与 `pgr kmer`/`pgi`/`asm map` 同套。

## 8. 验证计划

### 单元测试（libs/olc/）

* overlap：构造两个已知精确 suffix/prefix overlap 的 unitigs → PAF 记录
  的坐标/L/方向正确；contain 与 dovetail 分类正确；rc 方向正确；
  重复 k-mer（polyA 区）不产生错误重叠；
* layout：线性链 / 分支（bubble）→ 只走一条；repeat unitig → 路径断开；
  确定性（相同输入两次运行逐字节一致）；
* consensus：缝合正确性（含跨多 k 的 contain 不引入重复碱基）。

### 集成测试（tests/cli_asm_olc.rs）

* 合成基因组（随机 ~2 kb 序列）→ 生成多份 reads（子串 + rc，覆盖 ~20×）
  → `pgr asm olc` → contigs 与基因组逐段精确一致（identity 100%）；
* 阶段管道形态（`asm unitig` ×3 k → `asm ovlp` → `asm layout` →
  `asm cns`）与驱动器输出一致；
* Lambda 真实 reads 冒烟（`tests/bbtools/Lambda/R1.2k.fq.gz` 等）：
  contig 数 / N50 合理、无 panic。

### 验收门

每里程碑 `cargo fmt` + `cargo clippy --all-targets -- -D warnings` +
`cargo test` 全绿；命令注册与 docs 由 `cli_consistency.rs` 约束。

## 9. 里程碑

| 里程碑 | 内容 | 验证 |
|---|---|---|
| M1 | `libs/olc/overlap.rs` + `pgr asm ovlp` | ✅ 5 单测 + 3 集成测试 |
| M2 | `libs/olc/layout.rs` + `pgr asm layout` | ✅ 5 单测（线性/反向/分支/互惠/contain） |
| M3 | `libs/olc/consensus.rs` + `pgr asm cns` | ✅ 4 单测（含不一致 overlap 友好报错） |
| M4 | `pgr asm olc` 驱动器 + 集成测试 + `docs/asm.md` + todo/笔记更新 | ✅ 合成基因组重建 + 确定性 + 阶段管道等价 |

### 实现说明（与设计定稿的偏差）

* layout 的互惠检查实现在**连接端**（target 的 junction end 的 best edge
  指回当前 unitig），而非自由端——线性链因此可连续延伸；
* 阶段命令 `layout` 需要 unitig FASTA（不止 PAF），用于长度与命名校验，
  `cns` 同理；命名逻辑抽到 `cmd_pgr/asm/common.rs` 三命令共用；
* 驱动器的 unitig 命名 `k<k>:unitig_<id>`（不是 `<stem>:`），阶段管道用
  文件 stem；两者互不冲突（驱动器 `--keep-dir` 产物可直接喂阶段命令）；
* `asm unitig` 新增内存版 `assemble_unitigs_buf`（`assemble.rs` 最小重构：
  核心逻辑抽出 `assemble_unitigs_core`，写盘版行为不变）。

## 10. 不做 / 待决

* **不做**：气泡/孤儿合并（用户裁定）；允许错配的 overlap（unitig 精确
  假设，DBG 错误在固实阈值下已排除；真实数据暴露再议）；scaffolding
  （无配对需求）；列投票 consensus（v1 待真实数据）。
* **待决（数据驱动）**：repeat breaking 的覆盖度证据阈值（Canu
  `SPURIOUS_COVERAGE_THRESHOLD=6` / `ISECT_NEEDED_TO_BREAK=15` 的
  单元化版本）；contain unitig 是否参与 consensus 投票；短 unitig
  （< seed k）的 overlap 缺失处理。
* **v1 素材来源（2026-08-12）**：repeat breaking 的覆盖度阈值与实现路径
  有两个成熟参考——`references/skesa.md` §7.1（`FilterLowAbundanceNeighbors`
  fraction=0.1 多层过滤 + 可逆性检查）与 `references/metaMDBG.md` §9
  （渐进丰度过滤 t=1.1/10% 步长 + RepeatRemover 的桥接 reads 证据，
  pgr 用 `asm map` + `sam to-rg` + `rg coverage` 回放即等价设施）；
  多 k 反馈（SKESA clean_reads / metaMDBG unitig 反馈）为 v2 候选。
* **参考**：`canu-2.3/src/bogart/`、`wgs-8.3rc2/src/AS_BAT/` 源码随取随用；
  不引入其代码/依赖（Canu EOL）。

## 11. 相关文档

* `references/canu.md`（Canu OLC 源码分析 + §8 设计意图 + §8.5 实现后理解回写）
* `references/celera.md`（Celera 8.3rc2 源码分析 + 对照，§9 已按实现更新）
* `references/skesa.md` §7.1 + `references/metaMDBG.md` §9（v1 素材来源：
  fork 过滤/丰度过滤/桥接 reads repeat breaking）
* `design/kmer.md` §11/§12（k 范围、FastK 字节键唯一表示）
* `design/fq-assemble.md` §8（`asm unitig` 语义与 L: 边）
* `todo.md` §3（多 k unitig OLC 挂账项，本项目承接）

## 12. 真实数据验证：Lambda（2026-08-12）

数据：`tests/bbtools/Lambda/R1.fq.gz` + `R2.fq.gz`（SRR5042715，108 bp ×
20k 对 = 40×，Illumina PE）；参考 `BBTools-40.01/resources/lambda.fa.gz`
（NC_001416.1，48,502 bp）；基线 golden `tadpole_contigs31.fasta.gz`
（BBTools 同 reads 组装，48,214 bp / N50 1199 / 最长 4258）。

### 12.1 抓出的 bug（已修复）

**左向延伸坐标下溢**：layout 坐标回填对 prepend 的首个 step 用占位
`q_end=0`，`prev_end − overlap_len` 下溢 panic（真实数据触发，合成数据
只测了右向延伸）。修复：首步坐标从自身长度算起；overlap > 前步末端改为
友好报错（零 panic 策略）。回归测试 `seed_extends_both_directions` /
`inconsistent_overlap_is_error`。

### 12.2 结果

| 实验 | 输入 | k | contigs | N50 | 最长 | 完美贴回参考 | 参考覆盖（正链） |
|---|---|---|---|---|---|---|---|
| A | 原始 reads 40× | 21,51,81 | 52 | 3409 | 19035 | 40/52 | 65.5% |
| B | 纠错 reads 9×（merge.ecco） | 21,51,81 | 202 | 459 | 2129 | 187/202 (92.6%) | 86.0% |
| C | 纠错 reads 9× | 21,31,41 | 248 | 469 | 2129 | 228/248 (91.9%) | 82.6% |
| D | 原始 reads 40× | 21,31,41 | 67 | 2233 | 8282 | 51/67 (76%) | 62.8% |

### 12.3 结论与教训

1. **长 contig ≠ 嵌合**：A 的最长 contig（19,035 bp）是**单个 k81 unitig**
   （cov=1.0），前 1708 bp 匹配参考（ref 29307–31015），随后跳出的序列
   **在 reads 中实锤存在**（100 bp 探针 fwd/rc 均命中）且两侧都是参考
   匹配区——是相对 NC_001416 的**菌株插入变异**（~1.3 kb at ref 31015），
   不是错装。教训：参考菌株与 reads 不同源时，"完美贴回"是错误判据，
   验证需 reads 侧证据（unitig 由 solid k-mer 建成本身即内部一致性证据）。
2. **多 k 冗余 2.4×**（A 总长 116 kb vs 基因组 48.5 kb）：不同 k 的 unitigs
   覆盖同一区域、contain 重叠被排除在延伸外 → 输出重复。v1 需消减
   （contain 去重，见 §13）。
3. **覆盖度 vs 纯度权衡**：40× 原始 reads 出长 contig 但变异区/重复区
   干扰贴回；9× 纠错 reads 纯度更高（92.6% 贴回）但碎片化（N50 459）。
   宏基因组真实数据的推荐路径待定（等数据）。
4. **k 选择**：108 bp reads 下 21/51/81 优于 21/31/41（原始 40×）——
   大 k 特异性在重复区更稳；k 应随读长自适应（设计默认 21/51/81 面向
   更长 reads）。
5. **合成数据验证的盲区**：合成测试全部是右向延伸链，漏了左向——真实
   数据验证的价值再次体现。

### 12.4 reads 回贴验证（2026-08-12，确认全长 contig 正确）

预过滤后 OLC 把 Lambda 拼成**单条 48,387 bp contig**（≈ 48,502 参考）。
reads 回贴验证（`asm map`，完美匹配）：

* 40,000 reads 中 **34,697（86.7%）完美贴回** contig；对 NC_001416 参考
  只有 34,069（85.2%）——OLC contig 多捕获 628 条，正是参考缺失的
  变异区 reads；
* 覆盖剖面：平均深度 77.4、中位 78，**无 ≥50 bp 零深度缺口**、
  **无 ≥100 bp 低覆盖（<5）长段**——整条 contig 连续 reads 支持；
* 变异区（contig 1708–3000，相对 NC_001416 的插入）平均深度 64.5
  （min 21）——变异序列是 reads 实锤，不是错装。

结论：unitig 级预过滤 + 单条全长 contig 的正确性由 reads 侧证据确认；
后续 OLC 验证（宏基因组）沿用此口径（参考不匹配时看 reads 回贴而非
贴回率）。

## 13. v1 待办（真实数据驱动）

* **多 k 冗余消减：已完成（2026-08-12）**——两级：
  * 输出级：consensus 丢弃完全包含于更长 contig 的 contig（含 rc）；
  * unitig 级（`filter_contained`，布局前）：剔除被更长 unitig 完全包含
    的 unitig。Lambda 实测：unitigs 90→22（-76%），overlaps 386→50，
    layout 从 16 条碎片**合并为 1 条全长基因组 contig**（48,387 bp ≈
    48,502），16 条旧 contig 全部是它的子串（内容零丢失）。注意：
    过滤会改变 greedy 路径选择（这正是目的——多 k 冗余曾打断互惠链），
    "内容保留"而非"布局不变"。
* **repeat breaking 覆盖度证据**：桥接 reads 回放（`asm map` + `sam
  to-rg` + `rg coverage`），阈值参考 SKESA fraction / metaMDBG 语义；
  需 reads 侧验证口径（参考菌株不匹配时不能只用贴回率）。
* 真实宏基因组数据验证 + 调参。
