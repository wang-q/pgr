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
        *   分治核心: 根据 PAF 映射结果，将巨大的泛基因组拆解为独立的染色体或组件（Components）。
        *   每个组件包含一部分图结构和对应的序列片段。这使得后续的昂贵比对可以并行化。
    4.  Batch Alignment (`cactus-align --batch`):
        *   对每个拆分后的组件，运行标准的 Cactus 流程 (`setup` -> `caf` -> `bar` -> `reference`)。
        *   关键点: 这里利用了 Minigraph 的拓扑作为约束，比纯粹的从头比对更准确且高效。
    5.  Joining (`cactus_graphmap_join.py`):
        *   将各个组件的比对结果（HAL/VG）合并。
        *   生成最终格式：GFA (文本图), GBZ (压缩图索引), VCF (变异位点)。

*   对 `pgr` 的启示:
    *   分治策略: 处理人类全基因组级别的任务时，必须拆分。Minigraph 提供了一个绝佳的拆分基准。
    *   骨架优先: 先构建粗糙的骨架（Minigraph），再填充细节（Cactus/POA），这比一步到位更可行。
    *   PAF 的重要性: PAF 是连接不同工具（Minigraph -> Cactus）的通用语言，`pgr` 必须完善对 PAF 的支持。

### 2.4 启示
1.  资源管理: 精细计算每个步骤的内存/CPU 需求（参考 `cactus_workflow.py`），这对大规模计算至关重要。
2.  任务整合: 将频繁交互的步骤合并（Consolidated），减少中间文件落盘。
3.  统一配置: 使用 XML 集中管理所有算法参数，便于调优。

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
*   **工具**: `minigraph`
*   **输入**: 多个样本的基因组组装 (FASTA)。
*   **过程**: 
    1.  从参考基因组（Reference）开始。
    2.  迭代地将其他样本映射到当前的图上。
    3.  仅当发现大于一定长度（如 >50bp）的结构变异时，才修改图拓扑增加新节点。
*   **输出**: **rGFA** (Reference GFA) 格式的 SV 图。
*   **意义**: 确立了泛基因组的整体拓扑结构，解决了最困难的大片段重排和拷贝数变异问题。

#### Phase 2: 映射定位 (Graph Mapping)
*   **工具**: `minigraph` (with `-dist` option) 或 `cactus-graphmap` 封装
*   **过程**: 
    *   将所有输入样本的序列（包括参考序列）重新映射（Map）回 Phase 1 构建的 SV 骨架图。
    *   这一步**不进行**详细的比对（Alignment），只进行**定位**（Mapping）。
*   **输出**: **PAF** (Pairwise Alignment Format)。
    *   描述了每条输入序列在图上的路径（大致坐标）。
*   **关键点**: PAF 在这里充当了连接“线性序列世界”和“图世界”的通用胶水。

#### Phase 3: 图拆分 (Graph Splitting)
*   **工具**: `cactus-graphmap-split`
*   **过程**:
    *   根据 PAF 映射信息，识别出图中的连通分量或染色体级别的独立区域。
    *   将巨大的全基因组图拆解为多个独立的、较小的**子图 (Sub-problems)**。
    *   每个子图包含：
        1.  一部分图的拓扑结构（来自 Minigraph）。
        2.  所有映射到该区域的序列片段（FASTA序列）。
*   **意义**: 实现了任务的**并行化**。解决了内存瓶颈，使得处理人类全基因组成为可能。

#### Phase 4: 局部精细比对 (Base-level Alignment)
*   **工具**: `cactus-align` (调用 `cactus` 核心库: Setup -> Caf -> Bar -> Reference)
*   **过程**:
    *   对每个拆分后的子图，独立运行标准的 Cactus MSA 流程。
    *   **约束**: 利用 Minigraph 的拓扑作为“软约束”或“引导”，指导 Cactus 进行比对。
    *   **填充**: Cactus 会在骨架节点之间，通过 **abPOA** (Adaptive Banded POA) 算法计算细微的 SNP 和 Indel，完善比对细节。
*   **输出**: 局部的 **HAL** 或 **VG** 格式文件。

#### Phase 5: 合并与索引 (Joining & Indexing)
*   **工具**: `cactus-graphmap-join`
*   **过程**:
    *   将所有子图的比对结果“缝合”回一个完整的泛基因组图。
    *   处理边界连接问题。
*   **输出**: 最终的交付格式。
    *   **GFA 1.1**: 包含 Walk (W-line) 的完整图。
    *   **GBZ**: 压缩的图索引（供 `vg giraffe` 使用）。
    *   **VCF**: 导出的变异位点文件。

### 3.2 数据流与关键格式 (Data Flow)

架构中各个组件通过标准文本格式进行解耦：

```mermaid
graph TD
    FASTA[Input FASTA] --> Minigraph
    Minigraph -->|rGFA| Skeleton[SV Skeleton Graph]
    FASTA --> Mapper
    Skeleton --> Mapper
    Mapper -->|PAF| MapInfo[Mapping Info]
    Skeleton --> Splitter
    MapInfo --> Splitter
    Splitter -->|Sub-FASTA + Sub-GFA| BatchAlign[Batch Alignment (Cactus)]
    BatchAlign -->|Local HAL/VG| Joiner
    Joiner -->|GFA 1.1| FinalGraph
    Joiner -->|VCF| Variants
    Joiner -->|GBZ| Index
```

*   **rGFA**: 用于表示骨架，强调基于参考坐标（Stable Coordinate）。
*   **PAF**: 用于传递映射信息，极其灵活，易于解析和分割。
*   **GFA 1.1**: 最终图格式，使用 Walk (W) 表达单倍型，比 Path (P) 更紧凑且适合泛基因组。

## 4. 下游分析实战 (Hackathon 2023 笔记)

基于 `doc/sa_refgraph_hackathon_2023.md` 整理。
本节聚焦于利用 Minigraph-Cactus 构建的泛基因组图谱进行下游分析，特别是 Read Mapping 和变异检测。

### 4.1 序列比对 (Mapping Reads to the Graph)

#### Short Read Mapping (`vg giraffe`)
*   适用性：短读长（Short Reads）比对的首选工具。
*   输出格式：GAM 或 GAF（Graph Alignment Format）。
*   比对策略：
    1.  传统方法：比对到等位基因频率过滤后的图（Allele-Frequency Filtered Graph, `.d2.gbz`）。
    2.  个人泛基因组方法（推荐）：
        *   使用 `kmc` 从 Reads 中提取 k-mers。
        *   利用 k-mers 在未过滤的图（`.gbz`）上动态采样出“个人泛基因组”（Personal Pangenome）。
        *   优势：保留样本特有的稀有变异，同时排除无关的常见变异，提高下游分析准确性。

#### Long Read Mapping (`GraphAligner`)
*   适用性：目前长读长（Long Reads, 如 HiFi）比对的推荐工具（`vg giraffe` 长读长支持尚在开发中）。
*   注意事项：
    *   需要先将 `.gbz` 转换为 `.gfa`，并使用 `--vg-algorithm` 选项以保持坐标一致性。
    *   输出 GAM 格式时需指定 `-x vg`。
    *   旧版本 `GraphAligner` 可能不输出 Mapping Quality，导致后续 `vg pack` 需要关闭质量过滤 (`-Q0`)。

#### Surjecting (图比对投影)
*   功能：将图上的比对结果（GAF/GAM）投影（Surject）回线性的参考基因组路径（如 GRCh38）。
*   输出：标准的 BAM 文件。
*   用途：使泛基因组比对结果兼容现有的线性分析工具（如 DeepVariant, GATK, samtools）。
*   命令：`vg surject`。也可以在 `vg giraffe` 中直接使用 `-o bam` 一步到位。

### 4.2 变异检测 (Genotyping and Variant Calling)

#### Small Variants (DeepVariant)
*   流程：
    1.  使用 `vg surject` 生成线性 BAM 文件。
    2.  从图中提取参考序列 FASTA (`vg paths`)。
    3.  运行 DeepVariant（WGS 模型）。
*   注意：DeepVariant 是针对线性参考基因组设计的，但经过训练可以很好地处理来自 `vg giraffe` 的 surjected BAM。

#### Structural Variants (SV) with `vg call`
*   原理：基于图的覆盖度（Coverage）来检测变异。
*   流程：
    1.  Pack：使用 `vg pack` 从比对文件（GAF/GAM）生成覆盖度索引（`.pack`）。
    2.  Call：使用 `vg call` 结合 `.pack` 和 Snarls（图的拓扑结构文件）进行变异检测。
*   区分：
    *   *Genotyping*：确定图中已有的变异在样本中是否存在。
    *   *Calling*：确定 Reads 中的变异（可能不在图中）。

#### Structural Variants (SV) with `PanGenie`
*   特点：不依赖 Read Alignment。
*   原理：利用 HMM 模型，结合图的路径和 Reads 的 k-mers 分布，推断最可能的单倍体组合。
*   输入：
    *   经过预处理的 VCF（Phased, Sequence-resolved, Non-overlapping, Diploid）。
    *   原始 Reads (FASTQ)。
*   优势：速度快，能有效利用图中的单倍体信息。

### 4.3 可视化 (Visualization)
*   **Panacus**: 用于统计和可视化泛基因组图谱覆盖度的工具（Histgrowth curves），可展示样本多样性。
*   **Bandage-NG**: 图结构可视化（需要 GFA 格式）。
*   **ODGI**: 强大的图操作和可视化工具（`odgi viz` 1D 可视化）。

## 5. 动态更新比对 (Updating Alignments)

基于 `doc/updating-alignments.md` 和 `doc/cactus-update-prepare.md` 整理。
Cactus 支持在不重新计算整个比对的情况下，对 HAL 格式的比对结果进行增删改。这对于维护大型比对项目非常有价值。

### 5.1 核心工具: `cactus-update-prepare`

这是官方推荐的高级封装工具，用于生成更新比对所需的一系列命令（Preprocessing -> Alignment -> HAL Update）。它不直接执行，而是输出脚本供用户分步运行。

*   Warning: 在执行任何更新操作前，务必备份 HAL 文件。

#### 1. Adding to a Node (添加为子节点)
*   命令: `cactus-update-prepare add node ...`
*   场景: 将新基因组直接挂载到现有的祖先节点下。
*   原理: 仅需重新计算该祖先节点的比对块（Block）。
*   底层调用: `halReplaceGenome`。

#### 2. Adding to a Branch (拆分分支)
*   命令: `cactus-update-prepare add branch ...`
*   场景: 在父节点和子节点之间插入一个新的祖先节点，将新基因组挂在该新祖先下。
*   原理: 适合系统发生树拓扑结构发生变化的情况。需要推断新的祖先序列。
*   底层调用: `halAddToBranch`。

#### 3. Replacing a Genome (替换基因组)
*   命令: `cactus-update-prepare replace ...`
*   场景: 基因组组装版本更新。
*   原理: 本质是“删除旧版本” + “添加新版本（Add to Node）”。
*   底层调用: `halRemoveGenome` + `halReplaceGenome`。

### 5.2 底层 HAL 命令详解

如果需要手动控制或理解 `cactus-update-prepare` 的输出，可以参考以下底层命令：

#### 5.2.1 删除基因组
*   命令：`halRemoveGenome <hal file> <genome to delete>`
*   限制：只能删除叶节点（Leaf Genome）。
*   注意：HAL 文件大小不会自动减小。若需压缩体积，需使用 `h5repack` 或 `halExtract`。

#### 5.2.2 验证比对
*   命令：`halValidate --genome <genome> <hal file>`
*   建议：每次修改 HAL 文件后都应运行验证，确保文件结构完整。

## 6. 对 `pgr` 项目的启示

`pgr` 作为一个现代化的基因组分析工具箱，在设计泛基因组模块时，应深入借鉴 Cactus 的架构经验。

### 6.1 泛基因组架构启示
1.  **分治是必须的 (Divide & Conquer)**: 处理大规模基因组时，必须设计拆分（Splitting）机制。`pgr` 现有的 Chain/Net 结构天然适合作为拆分依据。
2.  **分层构建 (Layered Construction)**: 不要试图一步到位。先解决 SV（Structure），再解决 SNP（Sequence）。
    *   `pgr` 可以考虑引入类似的“骨架图”概念。
3.  **PAF 的核心地位**: PAF 不仅仅是比对结果，更是流程控制的中间语言。`pgr` 应加强对 PAF 的操作支持（目前已有基础，需强化 split/filter/join 功能）。
4.  **Minigraph 集成**: 考虑直接调用 `minigraph` 二进制或移植其核心逻辑来构建初始图。
5.  **图模型**: 逐步从“线性参考中心”模型（Chain/Net/MAF）向“图模型”（GFA/Walks）过渡。

### 6.2 现有模块改进方向
1.  **`pgr fas join`**:
    *   引入树的概念：虽然不需要完整的重建祖先，但在合并多个 Pairwise Alignment 时，应优先合并亲缘关系近的物种。
    *   参考序列引导：`join` 操作本质上就是 Cactus 中的“以参考基因组为锚点”的投影过程。确保 Gap 的插入不会破坏已有的对齐结构。
2.  **`pgr fas refine`**:
    *   对应于 Cactus 中的局部比对优化。
    *   Cactus 证明了 abPOA 在局部 MSA 中的效率和准确性。`pgr` 引入 `spoa` (SIMD POA) 是正确的路线。通过 POA 图，可以更准确地处理插入和缺失。
3.  **替代 `multiz` 的完整路径**:
    *   `pgr` 的演进路线应该是：
        1.  Pairwise Alignment: 继续完善 `lastz` / `chain` / `net` 模块（基础）。
        2.  Multiple Alignment: 用 `fas join` + `fas refine (POA)` 来替代 `multiz` 的 `tba` 流程。
        3.  Graph Alignment: 未来可考虑引入类似 Minigraph 的图构建能力（Pangenome 方向）。
4.  **图线性化与 Scaffold 排序**:
    *   `cactus-master/reference` 模块解决的 "Reference Problem" 对于 `pgr` 未来处理更复杂的组装或共识生成非常有启发。
    *   如果 `pgr` 需要实现 Contig 的排序（Scaffolding）或从泛基因组图中提取新的线性参考序列，Cactus 使用的“最大权重匹配”+“邻接评分”策略是一个标准的算法范式。

### 6.3 目录结构建议
参考 Cactus，`pgr` 的泛基因组模块可以组织为：
*   `pgr-skeleton`: 负责骨架构建（调用 minigraph 或类似算法）。
*   `pgr-map`: 负责序列到骨架的映射。
*   `pgr-split`: 负责任务拆分。
*   `pgr-align`: 负责局部精细比对（利用现有的 `fas consensus` / POA 模块）。
*   `pgr-join`: 负责结果合并。

## 7. 参考文献
*   Cactus Documentation: `doc/pangenome.md`
*   Minigraph Paper: Li, H. et al.
*   GFA Specification: `doc/gfa.md`

---
*文档更新时间：2026-02-11*
