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
- **依赖**：除 KatComp/KatGC 外的工具吃 FastK 的 `.ktab` + `.hist`（必须
  `-t`/`-t1` 全量表）；KatComp/KatGC **实际只流式读 `.ktab`，不读 `.hist`**
  （见 §3b，README 的通用表述更保守）；绘图依赖 R + 包
  （argparse/ggplot2/scales/viridis/cowplot）。
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
| `KatComp` | 两个数据集的 k-mer 谱**叉积** 3D 热图/等高线（KAT 家族） | 两个 `.ktab` |
| `KatGC` | k-mer 覆盖度 vs GC 的 3D 热图/等高线（KAT 家族） | `.ktab` |
| `PloidyPlot` | README 仅存 stub：**已被最新版 SmudgePlot 取代，无独立源码文件** | — |

> 注：KatComp/KatGC 与 CNplot/ASMplot 的绘图选项 `-l/-f/-s` **语义不同**：
> CNplot/ASMplot 的 `-l/-f/-s` = line/fill/stack 谱线；KatComp/KatGC 的
> `-l/-f/-s` = **contour（等高线）/heat（热图）/combo（热图+等高线叠加）**，
> 且两者**无 `-k` 保留中间数据**（README 称重算代价可接受）。二者也不走
> `.hist`+Logex，而是**流式扫 `.ktab` 表**直接累加一个二维矩阵再喂 R。

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
  （`MerquryFK.c` 的 `scan_asm`，逐 scaffold）：`miss` = asm 中不被 reads
  覆盖的 k-mer 数，`tots` = asm 全部 k-mer 数（来自 asm 自身 profile
  `aprof[x]!=0` 的位置，reads 相对 profile `rprof[x]==0` 计为 miss）；
  先按 k-mer 缺失比例推到**单碱基错误率**
  `err = 1 - pow(1 - miss/tots, 1/KMER)`，再 `QV = -10·log10(err)`。
  输出行 `miss \t tots \t err(%.4f) \t qv(%.1f)`——**QV 用的是单碱基错误率，
  不是裸的 miss/tots 占比**（后者会被 k 次方放大为明显更高的"错误率"）。
  - **source quirk**：逐 scaffold 的 `.qv`（`scan_asm`）写的是 `err` **分数**
    （0~1）；而整体 `OUT.qv`（main 里）同一列写的是 `100.*err` **百分比**。
    两处表头都叫 "Error %" 但量纲不一致，读输出时别被误导。
- completeness = reads 中 solid（高计数）k-mer 被组装覆盖的比例，由
  `Logex` 计算：单倍体 `A&.B[<thresh>-]`（asm ∩ reads 的高计数）；双倍体
  用 `A&.D[thresh-]`/`B&.D[thresh-]`/`C&.D[thresh-]` 分别求 asm1/asm2/并集
  覆盖 solid reads 的比例，分母统一为 `SOLID_COUNT`。
- **SOLID 阈值自动推断**（MerquryFK.c 的 `SOLID_THRESH` 与 HAPmaker.c 的
  `reliable_cutoff`，注意**并非同款**）：
  - MerquryFK：对 reads 的 `.hist` 从 `low+1` 起
    `for (k=low+1; hist[k] < hist[k-1]; k++)`（**严格下降**），停在直方图
    **不再严格下降**处 = 错误峰下滑结束/固体峰上升起点，取该 `k` 为
    `SOLID_THRESH`；随后 `SOLID_COUNT` 累加 `[SOLID_THRESH, high]` 的计数
    （completeness 的分母）。
  - HAPmaker：用 **`<=`**（非严格，`hist[k] <= hist[k-1]`），即停在第一个
    不再下降或持平处，判据更"提前"；且它**不算 SOLID_COUNT**，只取阈值。
  - 两者扫描前都先 `Modify_Histogram(hist,low,high,1)` 把直方图归一化到
    **unique 计数模式**（`hist[i]` = 出现 i 次的**不同** k-mer 数，而非实例
    数；切换时 `*i` / `/i` 互转，边界桶藏在 `hist[high+1/+2]`）。
  - 这是一个可复用的"从 k-mer 谱自动分错误峰/固体峰"模式：**先归一化到
    unique 模式 → 沿低覆盖度方向找"不再单调下降"的拐点**。
- **CN 谱的构造**（`cn_plotter.c`）：调一次 `Logex -H1000` 对 asm 表 A 与
  reads 表 R 求 **7 个集合**直方图：`B-A`（reads-only）、`B&.A[1..4]`
  （reads ∩ asm 拷贝数 =1..4）、`B&.A[5-]`（≥5）、`A-B`（asm-only）；
  标签 `read-only / 1 / 2 / 3 / 4 / >4`。核心思想 = 把 reads 的 k-mer 按其
  在 assembly 中的**拷贝数**分桶——这正是 copy-number 谱的精髓。
- **ASM 谱的构造**（`asm_plotter.c`）：单 asm 用 `B-A / B&.A / A-B`；
  双 asm 用 `C-#(A|B)`（reads 不在任一）、`C&.(A-B)`/`C&.(B-A)`（只在
  asm1/asm2）、`C&.#(A&B)`（共有）等 7 组集合，标签 `read-only /
  asm1-only / asm2-only / shared`。
- **相位块（phased block）算法**（`MerquryFK.c` `phase_blocks`+`merge_blocks`）：
  扫 asm profile + 母/父相对 hap profile，把**同极性重叠的 hap-mer** 合并成
  mark 站点；再按 `ANCHOR_LENGTH=20000` 与 `ANCHOR_MARK=5`（×KMER）挑出
  可靠块，把不可靠块在相反极性可靠块之间**以纯度最小化**切分。输出
  `phased_block.bed`（列 Scaffold/Start/End/Phase/Purity/Switches/Markers）与
  `.phased_block.stats`（#Blocks/Sum/Min/Avg/N50/Max）。这是**相感知组装 QC**
  的参考实现。
- 内部组合 CNplot/ASMplot/HAPplot；`-k` 保留 `.cni/.asmi/.hpi` 中间数据
  供重复绘图。
- **MerquryFK 本质是驱动脚本**：它自己不实现计数/集合并，而是通过 `system()`
  串联 **FastK**（给 asm 建表 `-t1 -p`、产相对 profile `-p:<reads>`）、**Logex**
  （`A&.B[%d-]`、`A|+B` 等集合表达式）、**Fastrm**（清理临时文件）。assembly
  输入是 FASTA，由内部 FastK 现建表；`reads` 必须已是 `.ktab`+`.hist`。

## 3b. KatComp / KatGC 的实现（pgr 可借鉴的流式扫表）

- 二者都**直接流式扫 `.ktab` 表**（`Open_Kmer_Stream`），不做直方图/Logex：
  - **KatComp**：两表 `T`/`U` **并行前缀归并**（`GoTo_Kmer_Index` 按前缀分片，
    逐元素 `mycmp` 求交集/差集），累加 `plot[KF1][KF2]` **叉积矩阵**
    （`JMAX×HMAX`，默认 HMAX=JMAX=1000，`-X/-Y` 可给上限）。
  - **KatGC**：扫单表，对每个 k-mer 算**GC 含量** `gcontent`（用预计算的
    `GC[256]`/`GCR[256]` **字节查找表**逐字节查 C/G 数，比逐碱基判断快），
    累加 `plot[GC][覆盖度]`（GC∈[0,KMER]，覆盖度∈[0,HMAX]）。
  - 峰值语义：都从 x=2 起找"首个不单调下降点"再沿上升段扫全局最大
    （与 CNplot 一致），`-x` 倍率默认 2.1。
  - 输出矩阵（坐标带 `.5` 半格中心）到 `.kx`/`.kgc`，再 `Rscript` 画
    contour/heat/combo。

## 3c. HAPmaker 三阶段管线（三人家系 hap-mer 建表）

`reliable_cutoff`（见 §3）在此被**调用两次**，各作用在不同表上：
1. `Logex -h '%s = A-B' '%s = B-A'` 先求母/父**特有** k-mer 表
   （`mat-pat` / `pat-mat`），各自算 `MSOLID`/`PSOLID`；
2. `Logex -h '%s = C&.A[MSOLID-]' '%s = C&.B[PSOLID-]'` 取**子代继承**的
   特有 k-mer（用上一步阈值过滤），再各自重算可靠性阈值；
3. `Logex '%s.hap = A[MSOLID-]' '%s.hap = B[PSOLID-]'` 用第二步阈值过滤出
   **最终 hap-mer 表** `<mat>.hap.ktab` / `<pat>.hap.ktab`。
三步全是 Logex 集合表达式 + 阈值过滤，无独立计数逻辑。

## 4. 与 pgr 的关联

> 对照 [kmer.md](../design/kmer.md) §9/§10（2026-08 已大改：直方图/.kgc/.pkp
> 均已落地，勿再引用旧的"缺直方图"结论）。

1. **输入基础已部分就绪**：MerquryFK 全家吃 `.ktab` + `.hist`（`-t1` 全量）。
   对照 [kmer.md](../design/kmer.md) §9/§10.2：pgr **现已实现 FASTK 字节兼容
   直方图**（`hist.rs` → `pgr kmer hist`，`.hist` 布局同 FastK）、**KatGC 兼容
   `.kgc`**（`gc.rs`）与 **profile `.pkp`**。旧的"缺直方图"缺口已闭合；但
   MerquryFK 还吃 **`.ktab`**（pgr 用自有的内存版 `.pkt`，不做 FASTK `.ktab`
   磁盘分桶，见 kmer.md §9 表）——要驱动外部 MerquryFK，reads 需由 FastK 产
   `.ktab`（或 pgr 只产出 `.hist` 供 `Histex`/GenomeScope 用）。
2. **KatGC 已原生实现**：anchr 2_fastk 直接用的 `KatGC -x1.9 -s Table-<k>`
   （[kmer.md](../design/kmer.md) §9.1）——pgr `kmer gc`（`gc.rs`，`.kgc`
   KatGC 兼容、实测逐行一致）已覆盖其**矩阵计算**；差异仅在渲染（pgr 用
   `--tex` heat 图，无 R）与输入（pgr 从自家表/序列算，不读 `.ktab`）。
3. **组装质量评估仍是 pgr 潜在新功能**：QV 公式（§3）是纯标量计算可 Rust
   实现；CN/ASM 谱本质是 **Logex 集合运算 + 直方图**（pgr 有直方图，仍缺
   Logex——kmer.md §9 缺口 3）；per-scaffold QV 与相位块依赖 **profile**
   （pgr 有 `.pkp`，但自建格式、不做 FASTK `.prof` 兼容）。可借鉴的算法：
   QV 单碱基错误率公式、CN 谱"按 assembly 拷贝数分桶"、KatGC 的 GC **字节
   查找表**、SOLID 拐点检测（unique 模式 + 不单调下降）。
4. **HAPmaker（三人家系 hap-mer）**：本质是 Logex 集合运算（父母特有 ∩
   子代继承，见 §3c）——pgr 目前仍无 Logex（kmer.md §9 缺口 3），此能力未
   落地。
5. **二进制格式对照**（对 pgr 有直接参考价值，均已在 §3/§3b 核实）：
   - `.ktab`：stub + `.<root>.ktab.1..N` 分片；canonical k-mer，表项 =
     `kbyte(=ceil(k/4))` 字节序列 + 2 字节 count；前缀压缩索引 `4^(4*ibyte)`。
   - `.hist`：`int32 k | low | high | int64 ilowcnt | int64 max_inst |
     int64[high-low+1]`，unique/instance 双模式。
   - `.prof`：`.pidx.N`（每序列偏移）+ `.prof.N`（RLE 压缩 u16 逐碱基 count），
     供 per-scaffold QV / 相位块用。
   这正是 kmer.md §10.2 三种格式（`.pkt`/`.pkp`/`.hist`）的**对照系**：
   pgr 的 `.pkt`/`.pkp` 自建单文件不兼容 FASTK，`.hist` 刻意兼容。

## 5. 局限

- 依赖 R 做图（与 pgr 无 R 依赖的风格冲突，需权衡）；
- **输入强绑定 FastK `.ktab`**：pgr 的 `.pkt` 自建格式不兼容，无法直接喂
  外部 MerquryFK；`.hist` 已兼容（可被 `Histex`/GenomeScope 读），但
  MerquryFK 的 `reads` 仍需 `.ktab`+`.hist` 成对存在。
- **PloidyPlot 无源码**：README 仅作为"已被 SmudgePlot 取代"的占位说明，
  不要当作可实现工具。
- **README 标注 First Feb 2021 / Current Aug 11, 2021，这只是文档编写时间，
  不代表项目停更**（参考 [kmer.md](../design/kmer.md) §2.3 的 FASTK 核对
  教训：README 2021-04 的快照 = FASTK-1.2，但上游 2022-12 / 2023-06 /
  2024-10 仍有提交）。本地快照无 ChangeLog/git 历史，上游活跃度无法判断。
- **source quirk**：`HAPmaker.c` 的头部注释是 `"Command line utility to
  produce CN-spectra plots"`——明显是从 `CNplot.c` 复制粘贴漏改（实际实现是
  hap-mer 建表，Usage 为 `mat pat child [.ktab]`），读源码时别被注释误导。
- **source quirk**：逐 scaffold 的 `.qv` 与整体 `OUT.qv` 的 "Error %" 列
  量纲不一致（分数 vs 百分比，见 §3），解析时需区分。

---

*参考来源: 本项目源码 `MERQURY.FK-1.2/`（README + MerquryFK.c/KatGC.c 等）*
