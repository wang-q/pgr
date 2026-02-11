# Cactus 分析笔记

本文档旨在总结 Cactus 项目的核心算法与架构，涵盖渐进式比对（Progressive Alignment）、泛基因组比对（Pangenome Alignment）及下游分析流程，为 `pgr` 项目提供参考。

## 1. Progressive Cactus (渐进式比对)

基于 `doc/progressive.md` 整理。

**Progressive Cactus** 是 Cactus 软件包的核心组件，用于对数百个脊椎动物级别的基因组进行多重序列比对（MSA）。

*   适用场景：不同物种（Cross-species）的全基因组比对。
*   不适用场景：同一物种内的样本比对（应使用 Minigraph-Cactus 泛基因组流程）。
*   核心输出：HAL 格式（Hierarchical Alignment Format），包含所有输入序列及重建的祖先序列。

### 1.1 核心原理：渐进式比对 (Progressive Alignment)

Cactus 采用自底向上（Bottom-up）的策略，依据输入的系统发生树（Phylogenetic Tree）进行比对分解。

1.  输入分解：用户必须提供一个 Newick 格式的系统发生树。
2.  迭代比对：
    *   从树的叶节点开始，找到亲缘关系最近的两个基因组（Sibling）。
    *   将这两个基因组进行两两比对。
    *   根据比对结果，推断并重建它们的祖先基因组（Ancestor）。
    *   这个祖先基因组将作为新的叶节点，参与上一层的比对。
    *   重复此过程，直到到达树的根节点（Root）。
3.  参数自适应：利用树的分支长度（Branch Lengths）来动态调整 `LastZ`（两两比对工具）的参数。分支越短（亲缘关系越近），参数越灵敏，比对速度越快且不失准确性。

### 1.2 接口与输入格式

运行命令的基本格式：
```bash
cactus <jobStorePath> <seqFile> <outputHal>
```

#### SeqFile (序列文件)

这是核心配置文件，包含两部分信息：
1.  Newick 树：定义物种间的进化关系。
2.  名称-路径映射：定义每个物种对应的 FASTA 文件路径。

示例：
```text
((Human:0.1,Chimp:0.1)Anc1:0.2,Gorilla:0.3)Anc0;
Human /path/to/human.fa.gz
Chimp /path/to/chimp.fa
Gorilla /path/to/gorilla.fa
```
*   `*` 前缀：可标记某个基因组为参考质量（Reference Quality），即该基因组可作为外群（Outgroup）。

#### 掩盖 (Masking)

*   Soft-masking：输入序列必须进行软掩盖（Soft-masking，重复序列用小写字母表示）。推荐使用 `RepeatMasker`。
*   Hard-masking：强烈不推荐（用 N 替换重复序列），会导致大量比对丢失。
*   预处理：Cactus 默认使用 `red` 或 `lastz` 进行预处理掩盖，以加速比对。

### 1.3 输出格式与工具

#### HAL (Hierarchical Alignment Format)

*   特点：以图结构存储多重比对，包含祖先序列，支持高效的随机访问。
*   工具：`halStats`（统计信息）, `hal2maf`（转换格式）。

#### MAF (Multiple Alignment Format)

虽然 HAL 是内部存储格式，但下游分析通常需要 MAF。Cactus 提供了 `cactus-hal2maf` 工具来高效生成 MAF。

*   `cactus-hal2maf` 的改进：相比旧的 `hal2maf`，它解决了碎片化问题，支持分布式计算，并利用 TAFFY 进行块归一化（Normalization）。
*   关键选项：
    *   `--refGenome`: 指定参考基因组（Reference），输出的 MAF 将以此为基准。
    *   `--outType single`: 生成单拷贝（Single-copy）MAF，过滤掉复杂的 paralogy，适合某些特定分析。
    *   `--outType consensus`: 生成共识序列。

### 1.4 Reference Module (构建参考序列)
基于 `cactus-master/reference` 源码分析。
Cactus 的一个关键步骤是从图结构中重建线性的祖先序列或参考序列。这一过程由 C 语言实现的 `reference` 模块处理。

*   核心问题 (The Reference Problem): 给定一组无序的、相互比对的序列片段（Blocks），如何确定它们的最佳线性顺序和方向，以形成一条连续的染色体序列。
*   算法策略:
    *   Matching Algorithms: 使用图匹配算法（如最大权重完美匹配）来连接片段的端点。
    *   Adjacency Scoring: 计算“邻接得分”（Z-score），基于序列间的支持度来判断两个片段是否应该相邻。
    *   Top-down Construction: 自顶向下地在每个“花朵”（Flower，Cactus 的递归分解单元）中构建参考路径。

### 1.5 Preprocessor Module (`preprocessor/`)
*   功能: 输入数据的预处理，清洗和屏蔽。
*   核心组件:
    *   `cactus_redPrefilter.c`: 过滤掉极其重复的序列（"Red" filter），减少后续比对计算量。
    *   `cactus_softmask2hardmask.c`: 将软屏蔽（小写字母）转换为硬屏蔽（N），防止比对器在重复区域产生不可靠比对。
    *   `cactus_sanitizeFastaHeaders.c`: 清洗 FASTA 标题，确保 ID 唯一且符合格式要求（无空格等）。
    *   `cactus_lastzRepeatMask.py`: 调用 Lastz 进行重复屏蔽的封装脚本。
*   启示: `pgr` 在处理大规模基因组时，也需要类似的“卫士”模块，确保输入数据的质量和一致性，特别是 ID 的标准化和重复序列的处理。

### 1.6 Bar Module (`bar/`)
*   功能: 核心的多序列比对 (MSA) 引擎。负责在 Cactus 图的 "Flower"（局部子图）内部计算高质量的碱基级比对。
*   核心组件:
    *   `bar.c`: 主入口。负责调度 "Flower" 的比对任务，支持 OpenMP 多线程并行处理。
    *   `poaBarAligner.c`: 集成了 **abPOA** (Adaptive Banded Partial Order Alignment) 算法。这是 Cactus 目前默认的高效 MSA 引擎。
    *   `rescue.c`: 实现 "Rescue" (拯救) 机制。用于找回那些因覆盖度过滤而被丢弃、但对维持图连通性或特定物种覆盖度至关重要的序列片段。
    *   `endAligner.c`: 处理 "End" (末端/邻接) 的比对，确保 Block 边界的拓扑一致性。
*   启示:
    *   MSA 算法: 验证了 abPOA 是泛基因组 MSA 的首选方案。`pgr` 目前集成的 SPOA/POA 方向正确。
    *   Rescue 策略: 在构建一致性序列时，不仅仅是简单的多数投票，还需要考虑“拯救”低频但重要的变异或单倍型，这对 `pgr fas refine` 很有参考价值。
    *   并行粒度: Cactus 在 Flower 级别（相当于 `pgr` 的 Block 或 Cluster 级别）进行并行，这比在序列级别并行更高效。

### 1.7 Caf Module (`caf/`)
*   功能: **Cactus Annealing/Alignment Filter** (Cactus 退火/比对过滤器)。它是连接两两比对 (Pairwise Alignment) 和最终泛基因组图结构的桥梁。
*   核心流程:
    1.  **Annealing (退火/捏合)**: 也就是 `stCaf_anneal`。将来自 Lastz/Blast 的两两比对结果（PAF/CIGAR）作为“约束”，将不同的序列（Threads）在比对位置“捏合”（Pinch）在一起，形成 Pinch Graph 中的 Block。
    2.  **Melting (熔化/松弛)**: 也就是 `stCaf_melt`。这是一种过滤机制。根据覆盖度（Degree）、物种支持数或进化树覆盖度（Tree Coverage），将那些支持度不够的 Block "熔化"（拆开）。这用于去除错误的或偶然的比对噪音。
    3.  **Topology Construction**: 最终确立泛基因组的拓扑结构，为后续的 `bar` 阶段（精细 MSA）做准备。
*   核心组件:
    *   `caf.c`: 模块入口，协调整个流程。
    *   `annealing.c`: 实现将比对转化为图结构的逻辑。
    *   `melting.c`: 实现基于图属性的过滤逻辑。
    *   `stPinchGraphs` (依赖库): 定义了 "Pinch Graph" 数据结构，这是 Cactus 处理比对的核心抽象：序列是 "Thread"，比对是 "Pinch/Block"。
*   与 Bar 模块的区别: `caf` 负责粗粒度的拓扑构建和过滤（决定哪些序列应该在一起），而 `bar` 负责细粒度的序列比对（决定在一起的序列具体如何对齐）。
*   启示:
    *   两阶段比对策略: `pgr` 可以借鉴这种“先拓扑 (Caf) 后序列 (Bar)”的分层策略。先用快速方法确定大概的图结构，清洗噪音后，再用昂贵的 POA/MSA 算法优化细节。
    *   Melting 机制: 在构建一致性序列或图时，不应只接受所有比对，而应有一个基于统计或生物学意义（如系统发生树）的“反悔/过滤”机制。

### 1.8 API Module (`api/`)
这是 Cactus 的核心数据结构定义层，所有其他模块（Bar, Caf, Reference）都基于此 API 操作图数据。

*   核心对象模型:
    *   **Flower (花)**: Cactus 图的递归核心单位。一个 Flower 包含一组 Sequence, Chains, Blocks，以及嵌套的子 Flower。它代表了泛基因组的一个局部或层级。
    *   **Disk (磁盘)**: 负责对象的序列化和持久化（KV Store）。Cactus 的设计允许处理超大基因组，因为对象可以按需从 Disk 加载和卸载。
    *   **EventTree (事件树)**: 关联的系统发生树（物种树）。
    *   **Sequence (序列)**: DNA 序列片段。
    *   **Block (块)**: 无空隙的多序列比对块 (Colinear Alignment Block)。
    *   **Chain (链)**: 一系列连接的 Blocks，类似于 Net 中的链。
    *   **Group (组)**: 连接 Blocks 的“节点”结构。
    *   **Cap (帽)**: 序列在 Block 边界的“端点”或“方向”。

### 1.9 Setup Module (`setup/`)
负责泛基因组构建的初始化阶段。
*   `setup.c`: 读取输入的 FASTA 文件和 Newick 树，构建初始的 "Root Flower"。
*   `.complete` 标记: 识别文件名中的 `.complete` 后缀，将其标记为完整的染色体（Telomere-to-Telomere），否则视为碎片（Fragment）。这影响后续的端点处理逻辑。

### 1.10 关键依赖 Submodules
*   `pinchesAndCacti`: 定义了 **Pinch Graph**。这是构建过程中的动态图结构（可变），而 API 中的 Flower/Block 是构建完成后的静态结构（相对稳定）。Caf 模块主要在操作 Pinch Graph。
*   `sonLib`: Benedict Paten 团队的基础 C 库，提供了类似 GLib 的容器（List, Hash, Tree）、异常处理和 I/O 工具。Cactus 深度依赖此库。
*   `hal`: 层次化比对格式 (Hierarchical Alignment Format) 的 C++ 实现库。Cactus 最终通常输出 HAL 格式。

### 1.11 对比分析: Cactus vs UCSC (PGR)
我们的 `pgr` 项目已经实现了 UCSC 标准的 Chain/Net 库。Cactus 的设计深受 UCSC 体系影响（Benedict Paten 曾在 UCSC 工作），但将其从“线性”扩展到了“图”。

*   **Block (块)**
    *   UCSC/pgr (Linear): **Alignment Block**。查询序列和目标序列的无空隙匹配片段 (Pairwise)。
    *   Cactus (Graph): **Block**。无空隙的多序列比对块 (MSA)。可能有 >2 条序列穿过同一个 Block。
    *   关系: Cactus Block 是 UCSC Block 的多序列泛化。在两两比对模式下，它们是等价的。

*   **Chain (链)**
    *   UCSC/pgr (Linear): **Chain**。一系列共线性的 Blocks，允许中间有 Gap。
    *   Cactus (Graph): **Chain**。连接 Blocks 的路径。
    *   关系: 概念基本一致。都是描述“共线性”结构。

*   **Hierarchy (层级)**
    *   UCSC/pgr (Linear): **Net (Fill/Gap)**。基于参考序列的层级结构。Gap 中可以填充下一层级的 Fill。
    *   Cactus (Graph): **Flower (Nested)**。基于图拓扑的递归结构。一个 Flower 可以包含子 Flower。
    *   关系: 核心对应点。UCSC Net 的 Gap 相当于 Cactus 的 Nested Flower。UCSC Net 的 Fill 相当于 Cactus 的 Chain/Block。

*   **Connectivity (连接性)**
    *   UCSC/pgr (Linear): **Coordinates**。通过坐标隐式连接。
    *   Cactus (Graph): **Group / Link / Cap**。通过指针和拓扑结构显式连接。
    *   关系: Cactus 显式建模了连接性（Adjacency），不仅靠坐标，这对于处理倒位、重复和复杂重排至关重要。

#### 深度解析：从 Net 到 Flower
*   Net 的局限: UCSC Net 是以参考序列为中心 (Reference-centric) 的。所有的层级关系（Level 1, Level 2...）都是为了解决“如何在参考基因组上展示其他序列”的问题。它本质上是把复杂的图结构“投影”到了一条直线上。
*   Flower 的进化: Cactus Flower 是无中心 (Reference-free) 的。它不依赖单一参考序列。
    *   当我们需要“投影”到参考序列时（例如生成 Net 文件），Cactus 实际上是遍历图，选择一条路径作为参考，将其他路径作为 Gap/Fill 挂载上去。
    *   递归: Cactus 通过将复杂的局部图（如由重复序列引起的纠缠）折叠成一个 Group/Nested Flower，从而在上一层保持图的简洁。这与 Net 中将复杂区域留作 Gap，在下一层 Net 中再详细展示 Fill 的思想异曲同工。

#### 对 `pgr` 的启示
1.  兼容性: `pgr` 现有的 Chain/Net 模块非常重要。它们是连接 Graph world (Cactus/GFA) 和 Linear world (Browser/IGV) 的桥梁。
2.  数据结构演进: 如果 `pgr` 未来要处理真正的泛基因组构建（不仅仅是操作现有格式），可能需要引入类似 `Flower` 的递归图结构，或者至少支持 GFA 的 Segment/Path 模型。
3.  算法迁移: 既然 Net Gap 和 Nested Flower 拓扑上同构，我们可以尝试将 Cactus 的一些递归算法（如求 Pinch Graph）映射到 Net 的递归处理上。

## 2. Python Wrapper Layer (`src/cactus/`)
`src/cactus` 是 Cactus 的 Python 胶水层，核心作用是利用 **Toil** 引擎编排分布式工作流。

### 2.1 核心流程 (Pipelines)
*   `progressive/`: **Progressive Cactus** 流程。
    *   解析进化树，自底向上调度。
    *   核心文件: `cactus_progressive.py`。
*   `refmap/`: **Minigraph-Cactus Pangenome** 流程。
    *   基于图的构建流程：`minigraph` (构图) -> `graphmap` (映射) -> `split` (拆分) -> `align` (比对) -> `join` (合并)。
    *   核心文件: `cactus_pangenome.py`, `cactus_graphmap.py`。
*   `pipeline/`: 通用工作流逻辑。
    *   `cactus_workflow.py`: 实现了 "Consolidated" 任务，将 Setup, Caf, Bar, Reference 打包在一个 Job 中执行，减少 I/O 开销。

### 2.2 模块封装 (Wrappers)
对应底层的 C 模块，负责参数准备和二进制调用：
*   `bar/`, `caf/`, `reference/`: 封装核心算法模块。
*   `preprocessor/`: 封装 `cactus_redPrefilter`, `lastzRepeatMask` 等。
*   `blast/`: 封装 `lastz` 或 `blast` 进行局部比对。

### 2.3 Refmap Pipeline (`src/cactus/refmap/`)
这是 Minigraph-Cactus 泛基因组构建的核心流程，采用了 "Graph-Map-Split-Align" (图-映射-拆分-比对) 的分治策略，极大地提高了扩展性。

*   流程总览 (`cactus_pangenome.py`):
    1.  Minigraph Construction (`cactus_minigraph.py`):
        *   使用 `minigraph -xggs` 迭代地构建 SV (结构变异) 图骨架。
        *   只捕获大片段变异 (>50bp)，速度快，作为整个泛基因组的“骨架”。
    2.  Graph Mapping (`cactus_graphmap.py`):
        *   将所有输入基因组序列映射回 SV 图。
        *   生成 PAF 文件，描述每条序列大致位于图的哪个位置。这一步不进行碱基级对齐，只做定位。
    3.  Splitting (`cactus_graphmap_split.py`):
        *   分治核心: 调用 `rgfa-split` 工具，根据 PAF 映射信息，识别出图中的连通分量或染色体级别的独立区域。
        *   将巨大的全基因组图拆解为多个独立的、较小的**子图 (Sub-problems)**。
        *   每个组件包含：
            1.  一部分图的拓扑结构（来自 Minigraph）。
            2.  所有映射到该区域的序列片段（FASTA序列）。
    4.  Batch Alignment (`cactus-align --batch`):
        *   对每个拆分后的组件，运行标准的 Cactus 流程 (`setup` -> `caf` -> `bar` -> `reference`)。
        *   关键点: 这里利用了 Minigraph 的拓扑作为约束，比纯粹的从头比对更准确且高效。
    5.  Joining (`cactus_graphmap_join.py`):
        *   将各个组件的比对结果（HAL 或 VG 格式）合并。
        *   生成最终格式：GFA (文本图), GBZ (压缩图索引), VCF (变异位点)。

*   对 `pgr` 的启示:
    *   分治策略: 处理人类全基因组级别的任务时，必须拆分。Minigraph 提供了一个绝佳的拆分基准。
    *   骨架优先: 先构建粗糙的骨架（Minigraph），再填充细节（Cactus/POA），这比一步到位更可行。
    *   PAF 的重要性: PAF 是连接不同工具（Minigraph -> Cactus）的通用语言，`pgr` 必须完善对 PAF 的支持。

### 2.4 关键脚本详解 (Implementation Details)

#### 2. `cactus_graphmap_split.py` (Splitting Logic)
*   **功能**: 将全基因组 GFA/PAF 切分为以染色体为单位的片段，便于并行处理。
*   **核心逻辑**:
    *   **Heuristic Contig Selection**:
        *   Regex: 支持通过正则表达式提取特定染色体（如 `chr*`）。
        *   Size: 按长度排序，优先选择大片段作为 Reference Contigs。
        *   Dropoff: 当 Contig 长度骤降时截断（`dropoff` 参数）。
    *   **`split_gfa` Function**:
        *   调用外部工具 `rgfa-split` 执行实际拆分。
        *   构建命令参数：`-Q` (Uniqueness), `-n`/`-T` (Specificity filters/thresholds)。
        *   **Hack**: 对 PAF 输出进行重命名处理（`get_faidx_subpath_rename_cmd`），将空格替换为制表符，以适配后续工具。
*   **关键代码**:
    ```python
    # src/cactus/refmap/cactus_graphmap_split.py
    cmd = ['rgfa-split', '-p', paf_path, '-b', out_prefix, ...]
    # Specificity filters
    for cv in query_cov_vals: cmd += ['-n', str(cv)]
    for tv in query_thresh_vals: cmd += ['-T', str(tv)]
    ```

#### 3. `cactus_hal2chains.py` (Format Conversion)
*   **功能**: 将 HAL 格式转换为 UCSC Chain/Net 格式，用于跨物种比对可视化。
*   **依赖工具**: `axtChain`, `faToTwoBit`, `pslPosTarget` (UCSC Utils)。
*   **工作流**:
    *   `check_tools`: 验证环境变量中是否存在 UCSC 工具。
    *   `get_genomes`: 从 HAL 文件解析物种树和基因组列表。
    *   `hal2chains_all`: 并行执行转换，支持 `halSynteny` 或 `halLiftover` 模式。
*   **输出**: `.chain.gz`, `.bigChain.bb`, `.bigChain.link.bb`。

#### 4. `dnabrnnMasking.py` (Repeat Masking)
*   **功能**: 使用 `dna-brnn` 识别并屏蔽 Alpha 卫星序列（Alpha Satellites）。
*   **Masking Actions**:
    *   `softmask`: 小写屏蔽（调用 `cactus_fasta_softmask_intervals.py`）。
    *   `hardmask`: 替换为 `N`。
    *   `clip`: 剪切序列（利用 `bedtools subtract` 和 `samtools faidx`），并重命名子序列（如 `chr1_sub_9_15`）以避免特殊字符问题。
*   **BedTools 集成**: 使用 `bedtools sort` 和 `bedtools merge` 合并屏蔽区间。

### 2.5 启示
1.  **资源管理**: 精细计算每个步骤的内存/CPU 需求（参考 `cactus_workflow.py`），这对大规模计算至关重要。
2.  **任务整合**: 将频繁交互的步骤合并（Consolidated），减少中间文件落盘。
3.  **统一配置**: 使用 XML 集中管理所有算法参数，便于调优。
4.  **外部工具集成**: 高效利用现有的 C/C++ 二进制工具（如 `rgfa-split`, `dna-brnn`），通过 Python 胶水代码进行编排，而非重写所有逻辑。
5.  **数据清洗 (Sanitization)**: 在工具链传递过程中，主动处理特殊字符（如将空格替换为制表符，重命名复杂 ID），确保下游工具（如 Assembly Hubs）的兼容性。

## 3. Minigraph-Cactus 泛基因组构建架构

基于 `doc/pangenome.md` 和相关源码结构整理。

Minigraph-Cactus 是一种混合型（Hybrid）泛基因组构建流程，专为**同物种（Within-species）**或**近缘物种**群体设计。它结合了两种核心技术的优势：

*   **Minigraph**: 擅长快速构建复杂的结构变异（SV）骨架，但忽略细微序列差异。
*   **Cactus (Progressive)**: 擅长进行高精度的碱基级多重序列比对（MSA），但难以直接处理全基因组规模的复杂重排。

**核心设计哲学**: "先骨架，后细节" (Skeleton First, Details Later)。
**核心模式**: "分治与映射" (Divide, Map, and Conquer)。

### 3.1 核心工作流 (Core Pipeline)

整个架构采用经典的 **Map-Reduce** 风格，分为五个主要阶段：

#### Phase 1: 骨架构建 (Skeleton Construction)
*   **工具**: `minigraph` (由 `cactus_minigraph.py` 封装)
*   **过程**:
    1.  从参考基因组开始，构建初始图。
    2.  迭代地将其他样本映射到图上。
    3.  仅当发现大片段结构变异（如 >50bp）时，修改图拓扑增加新节点。
    4.  忽略小的变异（SNVs, Indels），保持骨架的简洁。
*   **输出**: **rGFA** (Reference GFA) 格式的 SV 图。
*   **意义**: 确立泛基因组的整体拓扑，解决最困难的大片段重排问题。

#### Phase 2: 映射定位 (Graph Mapping)
*   **工具**: `minigraph` (由 `cactus_graphmap.py` 封装)
*   **过程**:
    1.  将所有输入样本序列（包括参考基因组）重新映射（Map）回 Phase 1 生成的骨架图。
    2.  此步只做定位，不进行详细比对。
*   **关键参数**: `--mapCores`, `--delFilter`, `--minQueryUniqueness` (用于处理重复序列)。
*   **输出**: **GAF/PAF** 文件。描述每条序列大致位于图的哪个位置。

#### Phase 3: 图拆分 (Splitting)
*   **工具**: `rgfa-split` (由 `cactus_graphmap_split.py` 封装)
*   **过程**:
    1.  利用 PAF 映射信息，将巨大的全基因组图拆解为多个独立的、较小的**染色体级组件 (Components)**。
    2.  每个组件包含一部分图拓扑和对应的序列片段。
*   **策略**: 见前文 `2.4` 节关于 `cactus_graphmap_split.py` 的详细源码分析（正则匹配、大小排序、Dropoff）。
*   **意义**: 将内存需求从 TB 级降低到 GB 级，实现大规模并行化。

#### Phase 4: 局部比对 (Batch Alignment)
*   **工具**: `cactus-align` (调用 `cactus_consolidated`)
*   **过程**:
    1.  对每个拆分后的组件，独立运行标准的 Cactus 流程 (`setup` -> `caf` -> `bar` -> `reference`)。
    2.  在此阶段，Minigraph 的粗糙骨架被细化，所有序列细节（包括 SNVs）都被比对和解析。
*   **核心差异**: 相比传统的 Progressive Cactus，这里的比对受到 Minigraph 骨架的**约束**，因此更稳健，不易在重复区域迷失。

#### Phase 5: 结果合并 (Joining)
*   **工具**: `cactus-graphmap-join`
*   **过程**:
    1.  收集所有组件的局部比对结果（HAL 或 VG 格式）。
    2.  将它们缝合回一个完整的泛基因组图。
    3.  生成索引（GBZ, Snarls）以便下游工具（如 `vg`）使用。
*   **输出**:
    *   `.gfa.gz`: 完整的文本图。
    *   `.gbz`: 压缩的图索引（Giraffe 格式）。
    *   `.vcf.gz`: 导出的变异位点。
    *   `.hal`: 包含所有序列的 HAL 文件。

### 3.2 关键技术点 (Key Technical Points)

1.  **分级比对 (Hierarchical Alignment)**:
    *   Level 1 (Minigraph): 处理 >50bp 的结构变异。
    *   Level 2 (Cactus/POA): 处理 <50bp 的序列变异。
    *   这种分级策略完美解决了“既要看森林（SV），又要看树木（SNV）”的矛盾。

2.  **参考坐标系 (Reference Coordinate)**:
    *   虽然是泛基因组，但 Cactus 依然保留了基于参考基因组的坐标系（rGFA），这使得结果更容易与现有的基因组浏览器（UCSC, IGV）兼容。

3.  **波前算法 (Wavefront Algorithm)**:
    *   在最新的 Cactus 版本中，`WFA` (Wavefront Alignment) 被引入以替代部分的 `LastZ` 或 `POA`，特别是在处理长序列比对时，大幅降低了内存消耗。

## 4. 下游分析实战

基于生成的 HAL 或 GFA 文件，可以进行多种分析：

### 4.1 格式转换与导出
*   **HAL to MAF**: `cactus-hal2maf`。生成标准的 MAF 格式，用于 PhyloP/PhastCons 进化保守性分析。
*   **HAL to VCF**: `hal2vcf`。将变异导出为 VCF，用于群体遗传学分析。
*   **HAL to Chain/Net**: `cactus_hal2chains.py`。生成 UCSC Chain/Net，用于 LiftOver 坐标转换。
*   **GFA to VCF**: `vg deconstruct`。基于图路径导出 VCF，通常比 `hal2vcf` 更准确，尤其是在复杂 SV 区域。

### 4.2 覆盖度与深度分析
*   **HAL Coverage**: `halStats --coverage`。统计每个基因组在图中的覆盖比例。
*   **Depth**: 使用 `vg depth` 或 `odgi depth` 在 GFA 图上计算每个节点的深度，识别拷贝数变异 (CNV)。

### 4.3 比较基因组学
*   **PhyloP**: 利用 MAF 文件运行 `phyloFit` 和 `phyloP`，计算每个碱基的进化加速或保守分数。
*   **CAT (Comparative Annotation Toolkit)**: 利用 HAL 文件将参考基因组的注释（Genes, Transcripts）投影到其他物种上。

## 5. 动态更新比对 (Dynamic Pangenome)

Minigraph-Cactus 的一个重大突破是支持**动态更新**，无需重新运行整个流程。

*   **场景**: 已经构建了 100 个人的泛基因组，现在要新加 10 个人。
*   **流程**:
    1.  **Map**: 将新样本映射到现有的 GFA 图上 (`minigraph -M`).
    2.  **Align**: 仅对新样本和相关区域进行局部比对。
    3.  **Augment**: 将新比对结果并入图中。
*   **优势**: 极大地降低了维护成本，使得“生长型”泛基因组成为可能。

## 6. 对 `pgr` 项目的启示

### 6.1 架构层面
1.  **拥抱 GFA/rGFA**: `pgr` 目前主要基于 Chain/Net/AXT。虽然这些是经典格式，但 GFA 是未来的标准。建议 `pgr` 增加对 GFA (GFA 1.0/1.1) 的读写支持，作为内部的高级数据模型。
2.  **引入分层思想**: 不要试图用一种算法解决所有问题。
    *   对于大片段重排，参考 Minigraph 的图映射策略。
    *   对于细节比对，继续优化现有的 POA/SPOA 模块。
3.  **重视中间格式**: Cactus 成功的关键之一是 HAL。它既是存储格式，也是分析格式。`pgr` 是否需要定义自己的二进制索引格式（类似 `.bs` 的扩展），以支持图结构的快速随机访问？

### 6.2 算法层面
1.  **启发式拆分**: 学习 `cactus_graphmap_split.py` 的逻辑，在处理大基因组时，先用 PAF 做快速定位和拆分，再并行处理。
2.  **图对其 (Graph Alignment)**: 探索将序列映射到图（Sequence-to-Graph mapping）的算法，而不仅仅是序列对序列（Sequence-to-Sequence）。这是解决复杂变异（如 MHC 区域）的终极方案。

### 6.3 工程层面
1.  **Python + C/Rust**: Cactus 使用 Python 做胶水（Toil Workflow），C 做核心（Bar/Caf）。`pgr` 使用 Rust 实现核心，这是一个巨大的优势（性能+安全）。但可以考虑提供 Python binding，方便集成到类似 Toil 的工作流中。
2.  **工具链复用**: 不要重造轮子。对于 `LastZ`, `Minigraph`, `WFA` 等成熟工具，可以直接调用或集成库，专注于将它们串联起来解决特定问题。

## 7. 参考文献
*   **Paper**: "Progressive Cactus is a multiple-genome aligner for the thousand-genome era" (Nature, 2020).
*   **Paper**: "Minigraph-Cactus: Constructing the pangenome graph" (Nature Biotechnology, 2023).
*   **Repo**: https://github.com/ComparativeGenomicsToolkit/cactus
