# Cactus 分析笔记

本文档旨在总结 Cactus 项目的核心算法与架构，涵盖渐进式比对（Progressive Alignment）、泛基因组比对（Pangenome Alignment）及下游分析流程，为 `pgr` 项目提供参考。

## 1. Progressive Cactus (渐进式比对)

基于 `doc/progressive.md` 整理。

**Progressive Cactus** 是 Cactus 软件包的核心组件，用于对数百个脊椎动物级别的基因组进行多重序列比对（MSA）。

*   **适用场景**：不同物种（Cross-species）的全基因组比对。
*   **不适用场景**：同一物种内的样本比对（应使用 Minigraph-Cactus 泛基因组流程）。
*   **核心输出**：HAL 格式（Hierarchical Alignment Format），包含所有输入序列及重建的祖先序列。

### 1.1 核心原理：渐进式比对 (Progressive Alignment)

Cactus 采用**自底向上（Bottom-up）**的策略，依据输入的系统发生树（Phylogenetic Tree）进行比对分解。

1.  **输入分解**：用户必须提供一个 Newick 格式的系统发生树。
2.  **迭代比对**：
    *   从树的叶节点开始，找到亲缘关系最近的两个基因组（Sibling）。
    *   将这两个基因组进行两两比对。
    *   根据比对结果，推断并重建它们的**祖先基因组（Ancestor）**。
    *   这个祖先基因组将作为新的叶节点，参与上一层的比对。
    *   重复此过程，直到到达树的根节点（Root）。
3.  **参数自适应**：利用树的分支长度（Branch Lengths）来动态调整 `LastZ`（两两比对工具）的参数。分支越短（亲缘关系越近），参数越灵敏，比对速度越快且不失准确性。

### 1.2 接口与输入格式

运行命令的基本格式：
```bash
cactus <jobStorePath> <seqFile> <outputHal>
```

#### SeqFile (序列文件)

这是核心配置文件，包含两部分信息：
1.  **Newick 树**：定义物种间的进化关系。
2.  **名称-路径映射**：定义每个物种对应的 FASTA 文件路径。

示例：
```text
((Human:0.1,Chimp:0.1)Anc1:0.2,Gorilla:0.3)Anc0;
Human /path/to/human.fa.gz
Chimp /path/to/chimp.fa
Gorilla /path/to/gorilla.fa
```
*   `*` 前缀：可标记某个基因组为参考质量（Reference Quality），即该基因组可作为外群（Outgroup）。

#### 掩盖 (Masking)

*   **Soft-masking**：输入序列必须进行软掩盖（Soft-masking，重复序列用小写字母表示）。推荐使用 `RepeatMasker`。
*   **Hard-masking**：强烈不推荐（用 N 替换重复序列），会导致大量比对丢失。
*   **预处理**：Cactus 默认使用 `red` 或 `lastz` 进行预处理掩盖，以加速比对。

### 1.3 输出格式与工具

#### HAL (Hierarchical Alignment Format)

*   **特点**：以图结构存储多重比对，包含祖先序列，支持高效的随机访问。
*   **工具**：`halStats`（统计信息）, `hal2maf`（转换格式）。

#### MAF (Multiple Alignment Format)

虽然 HAL 是内部存储格式，但下游分析通常需要 MAF。Cactus 提供了 `cactus-hal2maf` 工具来高效生成 MAF。

*   **`cactus-hal2maf` 的改进**：相比旧的 `hal2maf`，它解决了碎片化问题，支持分布式计算，并利用 TAFFY 进行块归一化（Normalization）。
*   **关键选项**：
    *   `--refGenome`: 指定参考基因组（Reference），输出的 MAF 将以此为基准。
    *   `--outType single`: 生成单拷贝（Single-copy）MAF，过滤掉复杂的 paralogy，适合某些特定分析。
    *   `--outType consensus`: 生成共识序列。

### 1.4 Reference Module (构建参考序列)
基于 `cactus-master/reference` 源码分析。
Cactus 的一个关键步骤是从图结构中重建线性的祖先序列或参考序列。这一过程由 C 语言实现的 `reference` 模块处理。

*   **核心问题 (The Reference Problem)**: 给定一组无序的、相互比对的序列片段（Blocks），如何确定它们的最佳线性顺序和方向，以形成一条连续的染色体序列。
*   **算法策略**:
    *   **Matching Algorithms**: 使用图匹配算法（如最大权重完美匹配）来连接片段的端点。
    *   **Adjacency Scoring**: 计算“邻接得分”（Z-score），基于序列间的支持度来判断两个片段是否应该相邻。
    *   **Top-down Construction**: 自顶向下地在每个“花朵”（Flower，Cactus 的递归分解单元）中构建参考路径。

### 1.5 Preprocessor Module (预处理模块)
基于 `cactus-master/preprocessor` 源码分析。
在进行昂贵的比对计算前，Cactus 提供了一系列 C/Python 工具来清洗和标准化输入数据，确保比对的稳定性和效率。

*   **过滤 (Filtering)**:
    *   `cactus_redPrefilter`: 过滤掉过短（默认 <1kb）或低复杂度（如单碱基重复）的 Contig。这对减少 Red (Repeats Detector) 的误报和计算负担至关重要。
*   **掩盖 (Masking)**:
    *   `cactus_softmask2hardmask`: 将软掩盖（小写字母）转换为硬掩盖（N）。
    *   `lastzRepeatMasking`: 包含脚本用于基于 LastZ 或 BED 区间进行特定的掩盖处理。
*   **标准化 (Sanitization)**:
    *   `cactus_sanitizeFastaHeaders`: 规范化 FASTA 标题，截断空格，甚至处理特殊的 Minigraph/GFA 命名格式（去除 `#` 前缀），防止下游工具因 ID 解析错误而崩溃。

## 2. Minigraph-Cactus Pangenome (泛基因组比对)

基于 `doc/pangenome.md` 整理。

**Minigraph-Cactus** 是一种专为**同物种内（Within-species）**或**近缘物种**设计的泛基因组构建流程。它解决了 Progressive Cactus 在处理群体水平数据时过于敏感、难以捕捉结构变异（SV）的问题。

> **相关格式说明**: 关于 Minigraph 和 Cactus 使用的图格式（GFA 1.0, GFA 1.1, rGFA）的详细规范，请参考 [GFA Format](gfa.md)。

*   **适用场景**：构建人类、酵母等物种的泛基因组图谱（Pangenome Graph）。
*   **核心理念**：结合 `minigraph` 的结构变异构建能力和 `Cactus` 的碱基级比对能力。

### 2.1 核心流程 (Pipeline)

该流程包含五个主要阶段，通常通过 `cactus-pangenome` 命令一键运行，但也支持分步执行：

1.  **Minigraph Construction (`cactus-minigraph`)**
    *   使用 `minigraph` 构建初始的 SV 图（GFA 格式）。
    *   从参考基因组开始，逐个添加样本，仅保留结构变异（>50bp），忽略细微差异。
    *   **特点**：快速构建骨架，但这时的图不包含碱基级的比对细节。

2.  **Graph Mapping (`cactus-graphmap`)**
    *   将所有输入样本的序列重新映射（Map）回上述构建的 SV 图。
    *   这一步确立了每个样本在图中的大致路径。

3.  **Graph Splitting (`cactus-graphmap-split`)**
    *   **分治策略**：将全基因组图按染色体拆分成多个子图（Sub-problems）。
    *   **目的**：降低内存消耗，实现并行计算。

4.  **Base-level Alignment (`cactus-align`)**
    *   **核心步骤**：在每个子图中，使用 Cactus 的算法进行精细的碱基级比对。
    *   **填补细节**：填补 `minigraph` 忽略的小变异（SNP、Indel），生成完整的比对图。

5.  **Graph Joining & Indexing (`cactus-graphmap-join`)**
    *   将各染色体的比对结果合并。
    *   生成多种下游分析所需的索引和格式（GFA, VCF, GBZ, etc.）。

### 2.2 关键输出格式

Minigraph-Cactus 的强大之处在于其丰富的输出格式，支持各种下游工具：

*   **GFA (Graphical Fragment Assembly)**: 标准的图格式（GFA 1.1），包含 Walk 线（W-lines）表示单倍体路径。
*   **GBZ (GBWTGraph)**: 高度压缩的只读图格式，专为 `vg giraffe` 设计，支持大规模路径存储。
*   **VCF (Variant Call Format)**: 传统的变异位点格式。支持嵌套变异（Nested Variants），默认使用 `vcfbub` 进行平铺（Flattening）处理。
*   **ODGI (.og)**: 优化的动态图接口格式，非常适合可视化（`odgi viz`, `odgi draw`）。

### 2.3 关键特性

*   **Reference Handling**: 必须指定一个参考样本（Reference Sample），该样本在图中是无环的（Acyclic），作为坐标系统的基准。
*   **Contig Filtering**: 默认过滤掉无法归类到特定染色体的小 Contig。可通过 `--permissiveContigFilter` 放宽限制。
*   **Haplotype Sampling**: 新版支持通过单倍体采样（Haplotype Sampling）替代传统的频率过滤（Filter graphs），显著提升了稀有变异的保留率和 `vg giraffe` 的比对性能。

## 3. 下游分析实战 (Hackathon 2023 笔记)

基于 `doc/sa_refgraph_hackathon_2023.md` 整理。
本节聚焦于利用 Minigraph-Cactus 构建的泛基因组图谱进行下游分析，特别是 Read Mapping 和变异检测。

### 3.1 序列比对 (Mapping Reads to the Graph)

#### Short Read Mapping (`vg giraffe`)
*   **适用性**：短读长（Short Reads）比对的首选工具。
*   **输出格式**：GAM 或 GAF（Graph Alignment Format）。
*   **比对策略**：
    1.  **传统方法**：比对到等位基因频率过滤后的图（Allele-Frequency Filtered Graph, `.d2.gbz`）。
    2.  **个人泛基因组方法（推荐）**：
        *   使用 `kmc` 从 Reads 中提取 k-mers。
        *   利用 k-mers 在未过滤的图（`.gbz`）上动态采样出“个人泛基因组”（Personal Pangenome）。
        *   **优势**：保留样本特有的稀有变异，同时排除无关的常见变异，提高下游分析准确性。

#### Long Read Mapping (`GraphAligner`)
*   **适用性**：目前长读长（Long Reads, 如 HiFi）比对的推荐工具（`vg giraffe` 长读长支持尚在开发中）。
*   **注意事项**：
    *   需要先将 `.gbz` 转换为 `.gfa`，并使用 `--vg-algorithm` 选项以保持坐标一致性。
    *   输出 GAM 格式时需指定 `-x vg`。
    *   旧版本 `GraphAligner` 可能不输出 Mapping Quality，导致后续 `vg pack` 需要关闭质量过滤 (`-Q0`)。

#### Surjecting (图比对投影)
*   **功能**：将图上的比对结果（GAF/GAM）投影（Surject）回线性的参考基因组路径（如 GRCh38）。
*   **输出**：标准的 BAM 文件。
*   **用途**：使泛基因组比对结果兼容现有的线性分析工具（如 DeepVariant, GATK, samtools）。
*   **命令**：`vg surject`。也可以在 `vg giraffe` 中直接使用 `-o bam` 一步到位。

### 3.2 变异检测 (Genotyping and Variant Calling)

#### Small Variants (DeepVariant)
*   **流程**：
    1.  使用 `vg surject` 生成线性 BAM 文件。
    2.  从图中提取参考序列 FASTA (`vg paths`)。
    3.  运行 DeepVariant（WGS 模型）。
*   **注意**：DeepVariant 是针对线性参考基因组设计的，但经过训练可以很好地处理来自 `vg giraffe` 的 surjected BAM。

#### Structural Variants (SV) with `vg call`
*   **原理**：基于图的覆盖度（Coverage）来检测变异。
*   **流程**：
    1.  **Pack**：使用 `vg pack` 从比对文件（GAF/GAM）生成覆盖度索引（`.pack`）。
    2.  **Call**：使用 `vg call` 结合 `.pack` 和 Snarls（图的拓扑结构文件）进行变异检测。
*   **区分**：
    *   *Genotyping*：确定图中已有的变异在样本中是否存在。
    *   *Calling*：确定 Reads 中的变异（可能不在图中）。

#### Structural Variants (SV) with `PanGenie`
*   **特点**：**不依赖 Read Alignment**。
*   **原理**：利用 HMM 模型，结合图的路径和 Reads 的 k-mers 分布，推断最可能的单倍体组合。
*   **输入**：
    *   经过预处理的 VCF（Phased, Sequence-resolved, Non-overlapping, Diploid）。
    *   原始 Reads (FASTQ)。
*   **优势**：速度快，能有效利用图中的单倍体信息。

### 3.3 可视化 (Visualization)
*   **Panacus**: 用于统计和可视化泛基因组图谱覆盖度的工具（Histgrowth curves），可展示样本多样性。
*   **Bandage-NG**: 图结构可视化（需要 GFA 格式）。
*   **ODGI**: 强大的图操作和可视化工具（`odgi viz` 1D 可视化）。

## 4. 动态更新比对 (Updating Alignments)

基于 `doc/updating-alignments.md` 和 `doc/cactus-update-prepare.md` 整理。
Cactus 支持在不重新计算整个比对的情况下，对 HAL 格式的比对结果进行增删改。这对于维护大型比对项目非常有价值。

### 4.1 核心工具: `cactus-update-prepare`

这是官方推荐的高级封装工具，用于生成更新比对所需的一系列命令（Preprocessing -> Alignment -> HAL Update）。它不直接执行，而是输出脚本供用户分步运行。

*   **Warning**: 在执行任何更新操作前，务必备份 HAL 文件。

#### 1. Adding to a Node (添加为子节点)
*   **命令**: `cactus-update-prepare add node ...`
*   **场景**: 将新基因组直接挂载到现有的祖先节点下。
*   **原理**: 仅需重新计算该祖先节点的比对块（Block）。
*   **底层调用**: `halReplaceGenome`。

#### 2. Adding to a Branch (拆分分支)
*   **命令**: `cactus-update-prepare add branch ...`
*   **场景**: 在父节点和子节点之间插入一个新的祖先节点，将新基因组挂在该新祖先下。
*   **原理**: 适合系统发生树拓扑结构发生变化的情况。需要推断新的祖先序列。
*   **底层调用**: `halAddToBranch`。

#### 3. Replacing a Genome (替换基因组)
*   **命令**: `cactus-update-prepare replace ...`
*   **场景**: 基因组组装版本更新。
*   **原理**: 本质是“删除旧版本” + “添加新版本（Add to Node）”。
*   **底层调用**: `halRemoveGenome` + `halReplaceGenome`。

### 4.2 底层 HAL 命令详解

如果需要手动控制或理解 `cactus-update-prepare` 的输出，可以参考以下底层命令：

#### 4.2.1 删除基因组
*   **命令**：`halRemoveGenome <hal file> <genome to delete>`
*   **限制**：只能删除叶节点（Leaf Genome）。
*   **注意**：HAL 文件大小不会自动减小。若需压缩体积，需使用 `h5repack` 或 `halExtract`。

#### 4.2.2 验证比对
*   **命令**：`halValidate --genome <genome> <hal file>`
*   **建议**：每次修改 HAL 文件后都应运行验证，确保文件结构完整。

## 5. 对 `pgr` 项目的启示

`pgr` 作为一个现代化的基因组分析工具箱，在设计 `join` 和 `refine` 模块时，可以借鉴 Cactus 的以下思想：

### 5.1 关于 `pgr fas join` 的改进方向
目前 `pgr fas join` 仅仅是基于坐标的线性堆叠。参考 Cactus：
*   **引入树的概念**：虽然不需要完整的重建祖先，但在合并多个 Pairwise Alignment 时，应优先合并亲缘关系近的物种。
*   **参考序列引导**：`join` 操作本质上就是 Cactus 中的“以参考基因组为锚点”的投影过程。我们需要确保 Gap 的插入不会破坏已有的对齐结构。

### 5.2 关于 `pgr fas refine` 的改进方向
`refine` 对应于 Cactus 中的局部比对优化。
*   **引入 POA (Partial Order Alignment)**：Cactus 在生成祖先序列时，本质上是在做 MSA。`pgr` 引入 `spoa` (SIMD POA) 正是符合这一趋势的正确路线。
*   **图的视角**：线性的 MSA 容易丢失结构变异信息。通过 POA 图，我们可以更准确地处理插入和缺失。

### 5.3 替代 `multiz` 的完整路径
Cactus 明确指出传统的 `multiz` 流程（基于 LastZ -> Chain -> Net）在处理大规模基因组时存在局限。
`pgr` 的演进路线应该是：
1.  **Pairwise Alignment**: 继续完善 `lastz` / `chain` / `net` 模块（基础）。
2.  **Multiple Alignment**: 用 `fas join` + `fas refine (POA)` 来替代 `multiz` 的 `tba` 流程。
3.  **Graph Alignment**: 未来可考虑引入类似 Minigraph 的图构建能力（Pangenome 方向）。

### 5.4 关于图线性化与 Scaffold 排序
`cactus-master/reference` 模块解决的 "Reference Problem" 对于 `pgr` 未来处理更复杂的组装或共识生成非常有启发：
*   **不仅仅是 MSA**: `pgr fas consensus` 目前更多是基于列的 MSA（如 POA）。但在处理大尺度结构变异或碎片化组装时，我们需要决定 Block 的顺序。
*   **图匹配的应用**: 如果 `pgr` 需要实现 Contig 的排序（Scaffolding）或从泛基因组图中提取新的线性参考序列，Cactus 使用的“最大权重匹配”+“邻接评分”策略是一个标准的算法范式。

---
*文档生成时间：2026-02-07*
