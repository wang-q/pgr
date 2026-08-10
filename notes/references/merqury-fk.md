# MERQURY.FK-1.2：基于 FastK 的组装质量评估（源码/README 分析）

> 2026-08 整理，源自对仓库内 `MERQURY.FK-1.2/`（Gene Myers + Arang Rhie，
> First: Feb 2021, Current: Aug 2021）的分析。MerquryFK 是原版
> [Merqury](https://github.com/marbl/merqury)（meryl 后端）的 **FastK 后端
> 重写**：把 R/Java/shell 脚本集重构为一批 C 命令行工具，显著提速。
> 关联：[fastk.md](fastk.md)（FastK 计数）与 [kmer.md](../design/kmer.md)
> （pgr kmer 基础设施与缺口）。

## 1. 概况与依赖

- **定位**：用 k-mer 谱评估**组装质量**（QV、完整性、copy-number 谱、倍性、
  三人家系分型），输入是 reads 与组装序列的 FastK k-mer 表。
- **依赖**：所有工具吃 FastK 的 `.ktab` + `.hist`（必须 `-t`/`-t1` 全量表）；
  绘图依赖 R + 包（argparse/ggplot2/scales/viridis/cowplot）。
- **通用约定**：已知类型参数的后缀可省略（`foo` → 自动找 `foo.fasta/fa/fastq/fq`）；
  选项参数位置任意。

## 2. 工具清单

| 工具 | 作用 | 输入 |
|---|---|---|
| `HAPmaker` | 三人家系 hap-mer 表：父母**特有**、**子代继承**、可靠（非错误）的 k-mer | `mat/pat/child.ktab` → `*.hap.ktab` |
| `CNplot` | reads vs assembly 的 **copy-number 谱图**（line/fill/stack，png/pdf） | `reads.ktab` + asm |
| `ASMplot` | assembly k-mer 谱：在/不在 asm1/asm2 的四种集合 | `reads.ktab` + asm1[/asm2] |
| `HAPplot` | 单倍型 blob 图：contig 大小 ∝ 长度，(x,y) = 母/父 hap-mer 数 | 两个 `.hap.ktab` + asm |
| `MerquryFK` | 汇总全部评估（见 §3） | `reads.ktab` + asm1[/asm2]（+ hap 表则 trio） |
| `KatComp` | 两个数据集的 k-mer 谱比较（KAT 家族） | 两个 `.ktab` |
| `KatGC` | k-mer 覆盖度 vs GC 分析（KAT 家族） | `.ktab` |
| `PloidyPlot` | 倍性图（SmudgePlot 改进版） | 两个 `.ktab` |

## 3. MerquryFK 的核心输出（QV / completeness）

```text
OUT.<asm>.spectra-cn.(ln+fl+st).(png|pdf)   # asm vs reads 的 copy-number 谱
OUT.<asm>.qv                                # 每条 scaffold 的错误率与 QV
OUT.<asm>_only.bed                          # asm 特有（reads 不支持）k-mer 位置
OUT.spectra-asm.(...).(png|pdf)             # assembly 谱（含/不含集合）
OUT.qv                                      # 整体 QV
OUT.completeness.stats                      # solid read k-mer 被 asm/并集覆盖的比例
# trio 模式：phased_block.bed/.stats/.blob、block.N、continuity.N 等
```

- QV 表格列：`Assembly Only / Total / Error % / QV`。精确公式
  （`MerquryFK.c`，逐 scaffold）：`miss` = asm 中不被 reads 覆盖的 k-mer 数，
  `tots` = asm 全部 k-mer 数；先按 k-mer 缺失比例推到**单碱基错误率**
  `err = 1 - pow(1 - miss/tots, 1/KMER)`，再 `QV = -10·log10(err)`。
  输出行 `miss \t tots \t err(%.4f) \t qv(%.1f)`——**QV 用的是单碱基错误率，
  不是裸的 miss/tots 占比**（后者会被 k 次方放大为明显更高的"错误率"）。
- completeness = reads 中 solid（高计数）k-mer 被组装覆盖的比例，由
  `Logex` 计算。
- **SOLID 阈值自动推断**（MerquryFK.c/HAPmaker.c 的 `reliable_cutoff` 同款）：
  对 reads 的 `.hist` 从 `low+1` 起扫描 `for (k=low+1; hist[k] < hist[k-1]; k++)`，
  停在直方图**不再严格下降**处 = 错误峰下滑结束/固体峰上升起点，取该 `k` 为
  `SOLID_THRESH`；`SOLID_COUNT` 再累加 `[SOLID_THRESH, high]` 的计数，作为
  completeness 的分母。这是一个可复用的"从 k-mer 谱自动分错误峰/固体峰"模式。
- 内部组合 CNplot/ASMplot/HAPplot；`-k` 保留 `.cni/.asmi/.hpi` 中间数据
  供重复绘图。
- **MerquryFK 本质是驱动脚本**：它自己不实现计数/集合并，而是通过 `system()`
  串联 **FastK**（给 asm 建表 `-t1 -p`、产相对 profile `-p:<reads>`）、**Logex**
  （`A&.B[%d-]`、`A|+B` 等集合表达式）、**Fastrm**（清理临时文件）。assembly
  输入是 FASTA，由内部 FastK 现建表；`reads` 必须已是 `.ktab`+`.hist`。

## 4. 与 pgr 的关联

1. **输入 = FastK k-mer 表 + 直方图**：MerquryFK 全家吃 `.ktab` + `.hist`
   （`-t1` 全量）——再次印证 [kmer.md](../design/kmer.md) §9 的缺口判断：
   pgr 已有 `KmerTable`（表），**缺直方图（.hist）**；补上直方图后，
   pgr 就具备了 MerquryFK 家族的输入基础。
2. **KatGC 是 anchr 2_fastk 直接使用的工具**（`KatGC -x1.9 -s Table-<k> ...`，
   见 [kmer.md](../design/kmer.md) §9.1）——若替换 anchr 的 2_fastk，
   KatGC 需要 pgr 原生实现（k-mer 表 + 直方图 + GC 计算），或保留外部调用。
3. **组装质量评估是 pgr 潜在新功能**：QV / completeness 是纯计算
   （可 Rust 实现）；绘图部分 pgr 有 `plot` 模块（非 R 依赖），但
   MerquryFK 的 R 图（ggplot2 风格）对齐成本高，需按需求取舍。
4. **HAPmaker（三人家系 hap-mer）**：本质是 k-mer 集合运算
   （父母特有 ∩ 子代继承），与 FASTK `Logex` 同族——pgr 目前无此能力，
   属 §9 缺口 3（Logex）的应用场景。

## 5. 局限

- 依赖 R 做图（与 pgr 无 R 依赖的风格冲突，需权衡）；
- 与 FastK 生态强绑定（`.ktab`/`.hist` 格式），pgr 若用 KmerTable 需自建
  等价输入；
- **README 标注 First Feb 2021 / Current Aug 2021，这只是文档编写时间，
  不代表项目停更**（参考 [kmer.md](../design/kmer.md) §2.3 的 FASTK 核对
  教训：README 2021-04 的快照 = FASTK-1.2，但上游 2022-12 / 2023-06 /
  2024-10 仍有提交）。本地快照无 ChangeLog/git 历史，上游活跃度无法判断。
- **source quirk**：`HAPmaker.c` 的头部注释是 `"Command line utility to
  produce CN-spectra plots"`——明显是从 `CNplot.c` 复制粘贴漏改（实际实现是
  hap-mer 建表，Usage 为 `mat pat child [.ktab]`），读源码时别被注释误导。

---

*参考来源: 本项目源码 `MERQURY.FK-1.2/`（README + MerquryFK.c/KatGC.c 等）*
