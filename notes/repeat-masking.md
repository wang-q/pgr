# pgr 重复标记（Repeat Masking）场景

> 本文件记录 pgr 的重复标记方案：`pgr pl ir` / `pgr pl rept` / `pgr pl trf` 三个检测命令、
> 与 `pgr fa mask` 的衔接、底层 FastK → Profex → spanr 管道机制、以及临时文件清理策略。
> 写作时间：2026-08-03。

## 1. 背景与动机

完整跑一遍 RepeatMasker（基于 Smith-Waterman 与 RepBase/Dfam 库）非常昂贵且耗时。
pgr 的思路是用 FastK 系工具做**快速近似**：不做逐条 repeat 的分类注释，只回答
"基因组上哪些区间是重复的"，然后把区间喂给 `pgr fa mask` 做 soft/hard masking。
这适合大规模基因组的快速屏蔽；需要注释级结果（family/class 标签）时仍应使用 RepeatMasker。

外部工具依赖（均须在 `$PATH`）：

*   FastK + Profex（FastK 套件）
*   spanr
*   trf（仅 `pl trf`）

FastK 本身的源码分析见 [[fastk.md]]，本文不重复其内部机制。

## 2. 命令分工

| 命令 | 重复类型 | 原理 | 输入 |
| :--- | :--- | :--- | :--- |
| `pgr pl ir` | 散在重复（interspersed） | 与重复库做 k-mer 富集比对 | 基因组 + repeat 库（Dfam/RepBase/TnCentral） |
| `pgr pl rept` | 基因组内重复（无库） | 自身 k-mer 深度比较 | 仅基因组 |
| `pgr pl trf` | 串联重复（tandem） | trf 的周期搜索 | 仅基因组 |

三者输出**同一种格式**：runlist JSON（`{"chr": "start-end,start-end,..."}`），
可直接作为 `pgr fa mask --runlist` 的输入，因此检测结果与屏蔽步骤天然闭环：

```bash
# 检测 + 屏蔽闭环
pgr pl rept genome.fa -o repeats.json
pgr fa mask genome.fa --runlist repeats.json -o masked.fa        # soft-mask（小写）
pgr fa mask genome.fa --runlist repeats.json --hard -o masked.fa # hard-mask（N）
```

## 3. 术语澄清：SD 序列不是真正的 repeats

`pgr sd`（BISER 移植，见 [[design/sd.md]]）检测的**分段重复（segmental duplications, SD）**
在重复标记语境里容易混淆，需要澄清：

*   SD 是祖先复制事件产生的**旁系同源（paralogous）共享片段**（如 T2T-CHM13 标准：
    > 1 kb 且 identity > 90%），它们虽然"序列重复出现"，但**不是转座子等真正的
    重复元件（repeats）**。
*   `pgr pl ir/rept/trf` 检测的是重复序列本身（转座子、rRNA 基因簇、串联重复等）；
    `pgr sd` 检测的是旁系同源片段。两者目的不同，**不要把 `pgr sd` 的输出当成
    repeat masking 的结果**。
*   实践中 SD 在比对/组装中会造成假比对（旁系同源片段会被多处匹配），因此检测出 SD 后，
    下游流程通常会把它们排除或特殊处理；但**先被屏蔽的从来不是 SD，而是真正的重复元件**——
    BISER 的输入就要求预先 soft-mask 重复序列（RepeatMasker/TRF 等），SD 恰恰是
    "屏蔽重复元件后仍剩余的高相似旁系同源片段"，BISER 找的就是它们。顺序是：
    屏蔽 repeats 在前 → 检测 SD 在后 → 排除 SD 在下游比对中，三者并不矛盾。
*   推论：若屏蔽后还要做 SD 搜索（对应 BISER 输入假设），屏蔽应**只用 `pgr pl ir` +
    `pgr pl trf`**（≈ T2T-CHM13 的 TRF + RepeatMasker 预处理），**不要用 `pgr pl rept`
    （自比较）**——它会把 SD 本身也当作"重复"屏蔽掉，屏蔽完 SD 搜索就找不到目标了。
    注意 `ir` 需要重复库（Dfam/RepBase），无库时该组合退化为只用 `trf`：

    ```bash
    # SD 搜索前的正确屏蔽：IR + TRF
    pgr pl ir genome.fa repeats.fa -o ir.json # 散在重复（需重复库）
    pgr pl trf genome.fa -o trf.json          # 串联重复
    spanr merge ir.json trf.json -o mask.json # 合并区间
    pgr fa mask genome.fa --runlist mask.json -o masked.fa
    ```

## 4. 检测管道实现

### 4.1 ir / rept：FastK → Profex → spanr

共享管道在 `src/libs/pl/repeat.rs`：

1.  **FastK**：
    *   `ir`：跑两次——先用 `-t` 对 repeat 库建表（`-Nrepeat`），再对基因组用 `-p:repeat`
        生成相对该表的 count profile（`-Ngenome`）。
    *   `rept`：只跑一次，`-p` 自比较生成基因组自身的 profile。
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

默认参数：`kmer=17`、`fill-kmer=2`、`fill-fragment=10`；`ir` 的 `min-len=300`，`rept` 的 `min-len=100`。

### 4.2 trf：trf → 解析 → spanr

`src/cmd_pgr/pl/trf.rs`：按染色体拆分 FASTA，逐条跑
`trf <chr>.fa <match> <mismatch> <delta> <pm> <pi> <minscore> <max_period> -d -h -ngs`，
用 `parse_trf_output`（`src/libs/pl/repeat.rs`）解析 `.dat`（少于 15 列的短行跳过），
再 `spanr cover` 合并输出。默认参数对应 TRF 常用设置（match=2、mismatch=7、delta=7、
pm=80、pi=10、minscore=50、max_period=2000）。

## 5. 临时文件与 FastK 库文件清理

用户曾担心 FastK 会在工作目录生成一批库文件（`*.ktab.*`）。该问题已由
`src/libs/pl/ctx.rs` 的 `PipelineCtx` 内建解决：

*   管道启动时创建 `tempfile::TempDir`（前缀 `pgr_rm_` / `pgr_rept_` / `pgr_trf_`）；
*   `enter()` 把 CWD 切进 tempdir，此后 FastK 的 `genome.ktab.*` / `repeat.ktab.*`、
    Profex 的 `prof.*.txt/.rg`、trf 的 `.dat` 全部落在 tempdir 内；
*   ctx drop 时 TempDir 自动删除，`CwdGuard` 保证出错时 CWD 也能恢复。

实测（2026-08-03，MG1655）：跑完 `pgr pl rept` 后 `/tmp` 无新增残留，tempdir 也不存在；
FastK `-P` 默认丢到 /tmp 的排序块由 FastK 自身清理。因此**无需**在代码里额外做删除动作。

## 6. 与 RepeatMasker 的差异（局限）

*   **无分类注释**：不输出 family/class 标签，只给区间。
*   **k-mer 敏感度**：依赖与库共享的精确 k-mer，分化较远的拷贝会漏检或碎成小片段；
    fill 步骤只能桥接短孔，无法恢复长距离分化的拷贝。
*   **输出是区间而非序列**：mask 后的序列需另行用 `pgr fa mask` 生成。
*   **依赖外部工具**：FastK / Profex / spanr（trf 还需 trf）。

## 7. 实测记录

| 命令 | 基因组 | 耗时 | 区间数 | 备注 |
| :--- | :--- | :--- | :--- | :--- |
| `pl rept` | MG1655 | ~0.35s | ~150 | 含 rRNA 3941442-3946950 等 |
| `pl trf` | MG1655 | ~1.6s | 84 | 串联重复 |
| `pl ir` | — | — | — | 本机无 Dfam/RepBase 库，未实测 |

soft-mask 衔接验证：`aaaaaatgcgcggtcagaa` 等区间正确转为小写。

## 8. 待办 / 注意

*   `ir` 需要用户自备重复库（Dfam/RepBase/TnCentral，下载与准备见
    [docs/repeats.md](../docs/repeats.md)），本机缺库，
    端到端测试待补。
*   若未来要接近 RepeatMasker 能力，可考虑对检测出的区间补一步 family 注释（如对区间
    重跑库比对），但目前无此需求，不做推测性设计。
