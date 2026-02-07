# GFA (Graphical Fragment Assembly) Format

GFA 是用于描述序列图（Sequence Graphs）的通用文本格式，广泛应用于基因组组装和泛基因组分析。它通过节点（Segments）和边（Links）来表示序列及其连接关系。

## 1. 版本与变体 (Versions & Dialects)

目前主要流通的版本包括：

*   **GFA 1.0**: 最基础且广泛支持的版本。定义了 Segment (S), Link (L), Path (P), Containment (C)。
*   **GFA 1.1**: 在 1.0 基础上增加了 **Walk (W)** 行，专门用于表示泛基因组中的单倍型路径（Haplotypes）。这是 Minigraph-Cactus 输出的标准格式。
*   **GFA 2.0**: 对 1.0 的泛化和重设计，支持更复杂的层级结构和间隙（Gaps），但与 1.0 **不兼容**。目前主流泛基因组工具（vg, odgi）主要支持 GFA 1.x。
*   **rGFA (Reference GFA)**: GFA 1.0 的一个特定子集/约定。它要求图中的节点必须基于一个线性参考序列进行坐标定义（通过特定的 Tag）。**Minigraph** 原生输出此格式。

## 2. GFA 1.0/1.1 核心语法

GFA 是 Tab 分隔的文本文件。

### Header (H)
```gfa
H	VN:Z:1.0
```

### Segment (S)
定义图中的节点（序列片段）。
```gfa
S	s1	ACGT	LN:i:4
S	s2	*	LN:i:1000  # * 表示序列未在文件中显式列出
```
*   字段: `Type`, `Name`, `Sequence`, `Optional Tags`

### Link (L)
定义节点间的连接（边）。
```gfa
L	s1	+	s2	+	0M
```
*   字段: `Type`, `From_Segment`, `From_Orient`, `To_Segment`, `To_Orient`, `Overlap`
*   `0M` 表示无重叠连接（CIGAR 格式）。

### Path (P) - GFA 1.0
定义一条路径（有序的节点列表）。
```gfa
P	path1	s1+,s2+	*
```
*   字段: `Type`, `Path_Name`, `Segment_Names`, `CIGARs`

### Walk (W) - GFA 1.1
相比 Path，Walk 更适合描述泛基因组中的样本和单倍型信息。
```gfa
W	sample1	1	chr1	0	1000	>s1<s2>s3
```
*   字段: `Type`, `SampleID`, `HaplotypeIdx`, `SeqID`, `Start`, `End`, `Path`
*   Path 格式为 `>s1` (forward) 或 `<s1` (reverse) 的连续字符串。

## 3. rGFA (Reference GFA) 特性

Minigraph 生成的 rGFA 使用特定的 Tags 来标记节点在参考基因组上的位置：

*   `SN:Z:name`: 参考序列名称 (Stable Name)
*   `SO:i:offset`: 偏移量 (Stable Offset)
*   `SR:i:rank`: 等级 (Stable Rank, 0表示主骨架)

示例：
```gfa
S	s1	ACGT	SN:Z:chr1	SO:i:0	SR:i:0
```
这意味着 `s1` 对应于参考序列 `chr1` 的起始位置。

## 4. 变异图数据模型 (Variation Graph Data Model)

基于 Erik Garrison 的 "Untangling graphical pangenomics" 总结。

### 核心三要素
1.  **Nodes (节点)**: DNA 序列片段。通常有数字 ID。
2.  **Edges (边)**: 允许的连接。因为 DNA 是双链的，边是**双向 (Bidirectional)** 的。这意味着有四种连接类型：
    *   `+/+` (Forward to Forward)
    *   `+/-` (Forward to Reverse complement)
    *   `-/+` (Reverse to Forward)
    *   `-/-` (Reverse to Reverse)
3.  **Paths (路径)**: 图中的游走 (Walks)。Path 具有双重身份：
    *   **Genomes**: 表示参考基因组或单倍型。
    *   **Alignments**: 表示外部序列比对到图上的路径。

### 核心操作：Edit
"Edit" 是变异图的基础操作。它指将一个新的序列比对（Alignment）“烘焙”到图中，从而扩展图的结构（例如引入新的 SNP 或 Indel 节点）。

## 5. 格式对应关系 (Formats Mapping)

vg 生态系统倾向于同时维护 Schema (Protobuf, 二进制) 和 Text (文本) 格式。

| 概念 | 文本格式 (Text) | 二进制/Schema 格式 (Protobuf) | 说明 |
| :--- | :--- | :--- | :--- |
| **Graph** | **GFA** (1.0/1.1) | **.vg** (subset of GFA) | 图的拓扑结构和参考路径 |
| **Alignment** | **GAF** (Graph Alignment Format) | **GAM** (Graph Alignment/Map) | 序列比对到图上的结果 |

*   **GAF**: 类似于 PAF/SAM，是 GAM 的文本表示。
*   **GAM**: 类似于 BAM，是 GAF 的二进制表示。

## 6. 工具生态 (Ecosystem)

### Minigraph
*   生成 **rGFA**。
*   构建基于结构的泛基因组图，保留结构变异。

### Minigraph-Cactus
*   生成 **GFA 1.1** (包含 W lines)。
*   结合了 Minigraph 的结构图和 Cactus 的碱基级比对。

### vg (Variation Graph)
*   "瑞士军刀"。
*   支持 GFA 1.0/1.1 导入导出。
*   常用命令 `vg convert`：
    *   `vg convert -g in.gfa -f > out.gfa`: 转换/标准化 GFA。
    *   `vg convert -W -g in.gfa > out.gfa`: 将 GFA 1.1 (W lines) 转换为 GFA 1.0 (P lines)，这对不支持 W line 的工具（如旧版 odgi）很有用。

### odgi (Optimized Dynamic Genome/Graph Implementation)
*   **论文**: "ODGI: understanding pangenome graphs" (Guarracino et al., 2022).
*   **核心理念**: 提供一套高效的命令行工具，用于操作大规模泛基因组图。
*   **通用坐标系统**: 利用图中的 Paths 作为通用坐标系统，支持将图操作转换为标准的 BED/PAF 格式，从而与现有的生物信息学工具链（如 SAMtools）互操作。
*   **核心格式**: `.og` (Optimized Graph)，一种高效的二进制格式，支持快速加载和随机访问。
*   **常用功能**:
    *   `odgi build`: GFA -> OG
    *   `odgi view`: OG -> GFA
    *   `odgi viz`: 可视化
    *   `odgi paths`: 提取路径信息
    *   `odgi stats`: 统计图的指标

### 评估与比较 (Evaluation & Comparison)

### Gretl
*   Rust 编写的 GFA 统计评估工具。
*   需要数值型 Node ID 和拓扑排序。

### pancat
*   **论文**: "Pairwise graph edit distance characterizes the impact of the construction method on pangenome graphs" (Dubois et al., 2025, Bioinformatics).
*   **核心功能**: 计算两个泛基因组图之间的**成对图编辑距离 (Pairwise Graph Edit Distance)**。
*   **原理**: 通过将相同的基因组“穿线” (Threading) 过两个图，来量化由于构建方法（如 Minigraph-Cactus 的参考顺序或 PGGB 的参数）导致的拓扑结构差异。
*   **应用**: 评估构建偏差，比较不同构建工具的一致性。
*   **语言**: Rust (`rs-pancat-compare`)。

### PGGB (PanGenome Graph Builder)
*   **论文**: "Building pangenome graphs" (Garrison et al., 2024, Nature Methods).
*   **核心理念**: 无偏差（Bias-free）的泛基因组构建。不依赖单一参考基因组，而是使用全对全（All-to-All）比对来构建变异图。
*   **流程**:
    1.  **wfmash**: 基于 Wavefront Alignment 算法的高效全对全比对器。
    2.  **seqwish**: 从比对结果中诱导（Induce）出变异图（Variation Graph）。
    3.  **smoothxg**: 对图进行平滑（Smoothing）和规范化，通过局部主成分分析（PCA）优化图的拓扑结构，使其更适合可视化和分析。
*   **特点**: 能够捕捉复杂的结构变异和重组事件，适合构建包含多个高质量基因组的图。
*   **格式**: 输出 GFA。

### TwoPaCo
*   **算法**: 双遍法构建压缩 de Bruijn 图 (cDBG)。
*   **特点**: 专为多个完整基因组设计，内存效率高（基于 Bloom Filter）。
*   **格式**: 支持直接输出 **GFA 1.0** 和 **GFA 2.0**。
*   **应用**: 适合无参考序列的泛基因组构建。

## 7. 对 pgr 的启示

1.  **格式支持**: 如果 `pgr` 未来涉及泛基因组操作，应优先支持 **GFA 1.0** 和 **GFA 1.1**。对于比对数据，应考虑支持 **GAF** 格式，它是图比对的文本标准。
2.  **转换能力**: 实现类似 `vg convert` 的 W-line <-> P-line 转换功能可能很有价值。
3.  **图线性化**: 参考 rGFA 的设计，如何在图中保留线性坐标系统对于下游分析（如变异调用）至关重要。
4.  **数据模型**: 理解 Path 的双重性（基因组 vs 比对）对于设计内部数据结构很重要。如果涉及图构建，需考虑 "Edit" 操作的实现。

---

# 附录：Untangling graphical pangenomics (译文)

> **原文**: [Untangling graphical pangenomics](https://ekg.github.io/2019/07/09/Untangling-graphical-pangenomics)
> **作者**: Erik Garrison
> **日期**: 2019年7月9日

五年来，我一直致力于为泛基因组参考系统构建一个简单、通用的图形模型：**变异图 (variation graph)**。
我并不孤单。
这项工作始于 [GA4GH](https://www.ga4gh.org) 数据工作组 (GA4GH-DWG)，该小组花费了数年时间辩论并提出了各种图形泛基因组模型。
在我基于变异图模型构建了一个工作系统 [vg](https://github.com/vgteam/vg) 之后，该小组的许多成员随后加入了我，继续其开发，组成了 [vgteam](https://github.com/vgteam)。
虽然我们努力产生的数据格式从未被 GA4GH 正式批准为标准，但它们现在已被广泛使用。

这个小组继续构建和完善基于变异图模型的生物信息学方法的参考实现。
我们要一起证明，[这些方法可以显著改善读取映射到变异等位基因的效果](https://europepmc.org/articles/pmc6126949)，尽管当参考图是从针对线性参考的比对检测到的小变异构建时，这种优势是微妙的，但[当图是从代表所有类型和规模变异的全基因组组装构建时，这种优势变得巨大](https://www.biorxiv.org/content/10.1101/654566v1.abstract)。
这些方法以及开发它们所吸取的教训，对于构建统一的数据模型至关重要，这些模型用于[利用第三代测序技术快速生成的大型基因组的从头组装](https://vertebrategenomesproject.org)。

### 图基因组的模型 (A model for graph genomes)

变异图在泛基因组数据结构中结合了三种类型的元素。
我们有 DNA 序列（节点），允许它们之间的连接（边），以及作为穿过图的路径的基因组（路径）。
节点具有标识符（通常是数字），路径具有名称（文本字符串）。
为了反映这些图的基因组用途，做出了一个让步。
它们是双向的，代表 DNA 的两条链，因此位置是指节点的正向或反向互补方向。
这意味着有四种边 (+/+, +/-, -/-, -/+)，每种边都暗示了其自身的反向互补。

有很多方法可以可视化这些图，但也许最具启发性的是我基于 graphviz 的 dot 开发的第一个方法。
这个渲染图（来自我的论文 [Graphical pangenomics](https://doi.org/10.17863/CAM.41621)）展示了一个变异图的片段，该图是基于 HLA 基因 H-3136 的 GRCh38 ALT 单倍型的渐进比对构建的。
我们没有任何倒置边或循环，但我们可以看到一个在 [VCF](https://vcftools.github.io/specs.html) 中无法表达的特征：嵌套变异。
在上面，我们看到了图形结构的核心。
在下面，路径由它们穿过的节点标识符序列表示。

![H-3136](http://ekg.github.io/assets/H-3136.dot.png)

### 基于 Schema 的变异图数据格式 (Schema based data formats for variation graphs)

与 Adam Novak（他帮助使其成为双向的）一起，我首先在[一个简短的 protobuf schema](https://github.com/vgteam/libvgio/blob/f6f93ddeeb97a3977ad8109f16bf031da718e40a/deps/vg.proto) 中实现了这个模型。
通过将此 schema 编译成一组类库，我们直接从 schema 构建了基于图的基因组学的 API 和一组数据类型。
这种方法非常符合 GA4GH-DWG 社区讨论的活动（例如：[1](https://github.com/ga4gh/ga4gh-schemas/issues/275) [2](https://github.com/ga4gh/ga4gh-schemas/issues/444)），旨在构建基于 schema 的交换格式，而不是基本的文本格式。

该 schema 包括基本的图实现（.vg 格式，相当于 [GFAv1](https://github.com/GFA-spec/GFA-spec) 的一个子集）和 *Alignment* 对象定义，允许我们在重测序中使用这些图。
其他人已经扩展了这个 schema 以支持特定的应用程序，如基因分型和图结构分解。
还有更多人在他们自己的项目中使用了它。
特别值得注意的是两个基于图的读取比对器，[GraphAligner](https://github.com/maickrau/GraphAligner) 和 [PaSGAL](https://github.com/ParBLiSS/PaSGAL)，它们实现了针对序列图的高效启发式和精确长读比对。
这些方法以及 [vg](https://github.com/vgteam/vg) 输出 GAM（图形比对/映射）格式。
GAM 将成对序列比对记录的概念扩展到序列比对到图的情况。
与 vg 中的所有数据格式一样，GAM 既有文本（JSON）也有二进制（protobuf）版本。
反映 .vg 与 GFA 的等价性，GAM 在语义上等价于最近提出的新生的基于行的 [GAF](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf) 格式。
vg schema 和一个用于读写其中文件的独立 API 现在位于 [libvgio](https://github.com/vgteam/libvgio) 中。

在 vg 中，*Path* 数据类型起着双重作用。
它既代表穿过图的行走，也代表序列到变异图的碱基级比对（它是 Alignment 对象的组件之一）。
我们在 GFA/GAF 中看到了相同的模式。
GFA 中的路径步骤有 cigar 字符串，GAF 中的比对记录也可以有。
然而，重要的是要注意这允许我们做什么。
我们可以支持变异图的“基本”操作，即*编辑 (edit)*。
在这里，我们将左侧的图 (*G*) 与包含 SNP 等位基因的比对 *z* 进行增强，产生右侧的图：

![Edit](http://ekg.github.io/assets/example_vg_construction_edit_only.png)

通过此操作和比对器，我们可以逐步构建变异图（如在 `vg msga` 中）。
如果我们构建一个发出作为 *Paths* 集的基因型的变异调用器，那么我们可以直接使用变异调用器的输出来使用相同的功能扩展参考系统（如在 `vg call` 中）。

### 泛基因组模型选择 (Pangenomic model choice)

变异图模型的设计旨在尽可能简单。
它对图结构或坐标没有任何断言。
这种设计是有意的，因为简单性允许通用性。
几乎任何类型的序列图都可以用作变异图的基础。
在 [vg](https://github.com/vgteam/vg) 中，我们实现了从 VCF 文件和参考、de Bruijn 图、字符串图、多序列比对、RNA 剪接图和全基因组比对开始的变异图构造函数。
如果你能构建一个序列图，vg 就可以将其用作参考系统。

尽管从线性参考和 VCF 文件构建的图是我们的主要关注点，但在我们[关于变异图读取比对的论文](https://europepmc.org/articles/pmc6126949)中，我们证明了我们可以使用 vg 将来自保留菌株的长 PacBio 读取比对到从其他六个从头 *S. cerevisiae* 组装构建的全基因组比对图（由 [Cactus](https://www.nature.com/articles/nbt.4227) 制作）。
我的论文涵盖了许多其他例子，包括淡水病毒宏基因组、细菌泛基因组、结构变异图、酵母剪接图和人类肠道微生物组。

这种通用性的代价是增加了开发时间，因为我们需要学习和发现更多的东西。
然而，随着推动 vg 的合作的进展，我们学会了如何处理这种复杂性。
我们发现，虽然变异图是一个理想的数据集成系统，但我们可能需要构建基因组图的技术转换，以实现针对它们的高效读取比对。
一般来说，这些图是更大图的严格子集，复杂度较低。
但也可能涉及复杂区域的展开或复制，以减少当许多变异在邻近发生时可能发生的路径爆炸。
我们保留了技术图和基础图之间的映射，在各自的最佳应用中使用每一个。
如果说这总是很容易，那就是贬低了我们的工作。
我们的成功表明，这是可能的，而且比那些选择 *先验地* 限制其泛基因组模型的小组所预期的更权宜。

### 路径为基因组图提供坐标和记忆 (Paths provide coordinates and memory to genome graphs)

许多基因组学研究人员在面对基因组图时，立即担心不再有稳定的坐标将注释放置在图上或将其与已知基因组联系起来。
图只是多序列比对的表示。
[比对是反复无常的](https://lh3.github.io/2014/07/25/on-the-graphical-representation-of-sequences)，并且可能会根据评分参数、顺序或输入序列的结构及其重叠而改变。
我们不应该相信相同的输入序列会导致相同的图。

在 GA4GH-DWG 中，我们花了大约一年时间（2014-2015）辩论试图解决这个问题的不同参考泛基因组数据模型。
许多参与者提出了在图结构中直接保留预先存在的参考坐标系的模型。
然而，很少有人提供工作实现来支持他们的论点，因此这些对话大多针对假设性问题。

作为回应，我提出了 vg 本身的第一个原型，并证明了我能够使用它将读取比对到从整个 [千人基因组计划](http://www.internationalgenome.org) 第3阶段发布构建的图形基因组模型。
由于其通用性，变异图模型可以消费任何其他提出的数据模型，这使得它具有直接的实际用途，以推动[不同图模型之间的比较评估，这是 GA4GH-DWG 的主要研究成果](https://www.biorxiv.org/content/10.1101/101378v1.abstract)。

变异图以尽可能简单的方式解决泛基因组坐标问题。
只要它们嵌入覆盖其大部分空间的参考路径，使得没有图位置距离参考位置“远”，我们就可以推导出图上位置的参考相对坐标。
如果我们以不同的方式从相同的序列构建图，不同参考路径中坐标之间的关系会发生变化，但我们不会获得或丢失任何坐标。
这个解决方案是不优雅的（我们没有得到一个单一的、不可变的坐标层级），但是实用的且极其灵活。

基因组图只是线性序列的比对，即使它们嵌入在图中，这些序列仍然是线性的。
我们可以利用此属性将相对于基因组图获得的结果投影到线性序列，例如通过将图中的比对“满射 (surjecting)”到覆盖大部分图的参考路径子集中（从 GAM 到 BAM）。
我的论文中的许多分析都使用了此功能，这是我们正在进行的古代 DNA 工作的关键方面。
在这种情况下，图只是一个黑盒子，用于在比对过程中减少参考偏差。
为了说明这个概念，下图显示了针对图及其对应于其嵌入路径的子图的节点空间的满射函数。
这种投影对于具有相同嵌入路径的任何一对图都有效。

![Surjection](http://ekg.github.io/assets/surject-example.png)

此功能还解决了基因组图的另一个关键缺点。
从序列及其之间的连接构建的图是无记忆的，并且允许现实中极不可能存在的等位基因重组。
变异图中的路径恰好解决了这个问题，并为对象提供了长距离结构，如果它仅仅是一个图形模型，这种结构将会丢失。

虽然路径存储起来很昂贵，但在给定物种中，给定基因座的大多数单倍型将是相似的，我们可以利用这一点将[基因组压缩成单倍型索引，如 GBWT](https://arxiv.org/abs/1805.03834)，它每存储一个碱基对的序列仅使用一小部分比特，同时提供线性时间的单倍型匹配和提取功能。
这个结果表明，存储我们用于构建泛基因组的所有序列可能是站得住脚的，甚至是可取的，即使这些序列的数量达到数百万。

### 图形泛基因组学的下一阶段 (The next phase of graphical pangenomics)

现在，进入我在这个子领域的经历的第五年，我们正处于一个转折点。
vg 和其他图形基因组方法远非仅仅是减少针对小变异的参考偏差的一种奇特方式，它们正准备成为管理和理解即将到来的大量完整、大型脊椎动物（包括人类）基因组的关键工具集。
为了解决这一需求，我在过去一年的大部分时间里致力于一系列工具，使我们要能够构建和操作代表大量大型基因组的全基因组比对的变异图。

第一个工具，[seqwish](https://github.com/ekg/seqwish)，消费由 [minimap2](https://github.com/lh3/minimap2) 在一组序列上生成的比对，并产生一个变异图（GFA 格式），该图无损地编码所有序列（作为路径）及其碱基对精确比对（在图拓扑本身中）。
这种方法比类似的方法如 [Cactus](https://github.com/glennhickey/progressiveCactus) 快几个数量级，并且比基于 de Bruijn 图的方法如 [SibeliaZ](https://github.com/medvedevgroup/SibeliaZ) 更灵活。
虽然它还只是一个原型，但我现在已经能够可靠地将 seqwish 应用于真核生物基因组集合，如人类、[medaka](https://twitter.com/erikgarrison/status/1127984636452274177) 和丽鱼 (cichlids)。
为了测试，我在笔记本电脑上几分钟内就构建了酵母的泛基因组，这在以前需要在大型计算节点上花费一天时间。
尽管如此，它的性能还有很大的提升空间，并且通过对其数据结构进行完全的基于范围的压缩（使用 Heng Li 很棒的 [implicit interval tree](https://github.com/lh3/cgranges)），它将能够扩展到从数百个人类基因组构建无损泛基因组。
因为它完全尊重其输入的比对结果，seqwish 可以充当泛基因组构建流程中一个独立的内核。
一旦它得到充分优化，难点将在于为给定的泛基因组结构化和 [选择最佳的比对集](https://github.com/natir/fpa)。

vg 的一个主要缺点是其内存中的图形模型，它每输入一个图的碱基可能消耗高达 100 字节的内存。
虽然对于博士项目的起步和原理验证方法来说是可以接受的，但这对于我们日常应用 vg 的问题来说是不可扩展的。
我们不得不通过各种方式来解决这个问题，例如将图划分为单个（人类规模）染色体，使用更大内存的机器，或将图转换为静态的简洁索引，如 [xg](https://github.com/vgteam/xg)。
去年冬天，在 [Matsue 举行的 NBDC/DBCLS BioHackathon](http://2018.biohackathon.org/) 上，Jordan Eizenga 和我决定通过构建一套 [动态简洁变异图模型](https://github.com/vgteam/sglib) 来解决这个问题。
这些模型提供了一致的 [HandleGraph](https://github.com/vgteam/libhandlegraph) C++ 接口，并且可以与任何编写为在该接口上工作的算法互换使用，[现在已经有很多这样的算法了](https://github.com/vgteam/vg/tree/master/src/algorithms)。
我自己在这个主题上的工作成果是 [odgi](https://github.com/vgteam/odgi)，它为使用 seqwish 构建的大型图提供高效的后处理。
目前，该方法支持图拓扑排序、修剪、简化和 kmer 索引等困难步骤。
我计划扩展它以支持读取比对和渐进式图构建，以补充 seqwish 和其他 [新兴方法如 minigraph](https://github.com/lh3/minigraph)。

为了说明这两个工具，考虑来自 [Yeast Population Reference Panel](https://yjx1217.github.io/Yeast_PacBio_2016/data/) 的一组基因组。
在这里，仅取 *cerevisiae* 基因组，我们运行 minimap2，过滤比对以去除 <10kb 的短比对，然后使用 seqwish 诱导变异图，并使用 [Bandage](https://github.com/ekg/yeast-pangenome/blob/master/steps.sh) 渲染输出图像：

![Seqwish yeast Bandage visualization](http://ekg.github.io/assets/seqwish10kbyeast.png)

这显示了一个相对开放的图，有一些塌陷和一些明显未比对的染色体末端（可能是因为我们的长度过滤器）。

使用 Bandage 渲染可能需要极长的时间，并且生成的图很难解释，因为不容易查看图中不同嵌入路径之间的关系。
然而，使用 `odgi viz`，我们可以获得嵌入的染色体如何相互关联以及如何与图拓扑关联的图像。
这种线性时间渲染方法在图像底部显示图拓扑，使用矩形图案显示边（挂在下方）连接图的排序序列空间（分隔图像顶部和底部的黑线）中的两个位置的位置。
路径显示在此拓扑之上，y 轴上的每个位置显示一条路径。
布局是非线性的，仅显示给定路径触及图中的哪些位置。
但是，由于基因组图是从线性序列构建的，它们具有流形线性属性，我们通常可以应用线性直觉来解释它们。

在下图中，基因组从上到下排列：S288c, DBVPG6765, UWOPS034614, Y12, YPS128, SK1, DBVPG6044。
染色体顺序部分源于 seqwish 的初始节点 ID 顺序分配。
我们立即看到 UWOPS034614（一个高度分化的菌株）中的一条染色体已被重新排列到其他染色体中。
同样值得注意的是，频繁的结构变异和 CNV 似乎嵌入在染色体末端。
这证实了 [描述这些组装的论文](https://www.nature.com/articles/ng.3847) 中的另一个发现，其作者观察到亚端粒区域是结构变异的热点。

![Seqwish yeast odgi viz](http://ekg.github.io/assets/seqwish_yeast_l10k.dg.png)

如果你好奇这是如何工作的，[我已经记录了为整个集合执行此操作的步骤](https://github.com/ekg/yeast-pangenome/blob/master/steps.sh)。
在接下来的几周里，我将更多地报告这些方法及其在人类基因组组装集合中的应用。

### 与其他泛基因组方法合作 (Working with other pangenomic methods)

我写这篇文章的部分动机是 [Heng Li 发布了一个类似的工具链，用于构建泛基因组参考图](https://lh3.github.io/2019/07/08/on-a-reference-pan-genome-model)。
这种方法是 GA4GH 图基因组对话期间出现的关于坐标系和参考图的想法的第一个实际实现。
Heng 提出我们需要一种新的数据模型来在泛基因组中编码稳定的坐标。
这个模型（在 rGFA 中实现）类似于当时讨论的“侧图 (side graph)”，也类似于 Seven Bridges Genomics 使用的层级模型。
它用指示其起源的信息注释 GFA 元素，这些起源来自于一个逐步构建的坐标层级。
实际上，rGFA 将能够支持我们在 GFA 中维护的全方位语义，但它通过这些注释简化了将稳定坐标系关联到图中的过程。
rGFA 为交换路径相对坐标提供了正式规范，我们长期以来一直在 vg 生态系统中的工具使用的各种索引中缓存这些坐标，这对于泛基因组模型的用户可能非常重要。

Heng 还提出了一个原型渐进式泛基因组构建算法 (minigraph)，该算法构建图以仅包含相对于已建立的参考或正在扩展的图在结构上新颖的序列。
minigraph 是 Heng 对简单、高效且概念清晰的生物信息学方法的执着追求的完美体现。
其性能特别令人振奋（我测得它比 seqwish 流程快 5-10 倍），并且只要它基于 minimap2 链接 (chaining)，我期望生成的图在捕捉序列之间的大规模关系方面具有高质量。
我很高兴现在有另一种构建泛基因组的方法可以扩展到我正在处理的问题规模，并且我希望在我自己的工作流程中使用 minigraph 的输出。

然而，我想指出 minigraph 产生的内容与我在这里列出的泛基因组模型之间存在明显的区别。
该方法的用户应该理解，这不是在泛基因组中记录基因组集合的通用解决方案，而是一种根据给定的渐进比对模型推导出与新颖序列相关的坐标层级的方法。
我担心，如果作为主要的泛基因组模型被采纳，目前形式的 minigraph 方法将延续参考偏差的长周期，这种偏差自重测序方法建立以来一直主导着基因组学的许多方面。

Minigraph 泛基因组相对于其输入序列是有损的。
除非这些序列之间没有小变异，否则不可能构建一个包含一组样本中所有序列的 minigraph。
包含在图中的序列将是顺序依赖的。
这将使我们面临各种变异的参考偏差，这些变异不够大或不够分化，不足以阻碍渐进比对算法。
触发此情况的确切配置对于最终用户可能是不透明的。
在这些图中，我们将很难知道何时应该期望在气泡的这一侧或那一侧找到变异，或者很难根据对图的读取比对确定给定的变异以前从未见过。

在基于 minigraph 的图中，我们没有任何小变异。
这简化了事情，并允许我们在图中使用遗留算法（序列到序列的链接和映射），但我们的结果表明，这将在针对排除在图之外的变异的比对准确性方面付出代价。
这也意味着在 minigraph 参考的背景下处理小变异将需要 VCF 的泛化，并将保留与该格式相关的所有复杂性。
其效果将是将基因分型的难度从比对（如在图基因组中）转移回变异调用（如在线性基因组中）。

Heng 提出的坐标模型将提供一个稳定的层级结构，但前提是图是以相同的方式、由相同的算法构建的，或者图是先前图的严格扩展。
对基础参考的更改将使图不兼容。
由于这种不灵活性，我不认为我们应该将我们使用的坐标系嵌入到图的结构中。
相反，两者应该是独立的，允许在研究期间根据需要构建和使用各种图结构，即使我们维护一个共同的坐标空间。
提供坐标（在变异图术语中）所需的一切只是一组主要覆盖图的路径。
这表明了这些模型如何协同工作。
使用 minigraph，我们可以构建一个基础 rGFA，为图提供覆盖坐标空间。
然后，将这些坐标空间作为路径嵌入到我们的图中，我们可以随心所欲地装饰、扩展或重建图，使用对我们分析有用的任何变异和基因组。
这一切都很容易，因为我们已经在读写相同的格式 (GFA)，尽管由于算法限制，这种交换可能并不总是双向的。

为了使这一切具体化，我留给你们 HLA 中 DRB1-3123 的 GRCh38 ALT 的两个版本。
该基因在人类群体中高度分化，我在 vg 的开发和探索过程中长期将其用作测试。
上面的图是由 minigraph 生成的，使用了比默认值更短 (2kb) 的比对阈值，而下面的图是由 seqwish 生成的，使用了相同的 2kb 比对过滤器。
它们在大小（minigraph 在图中产生多 20% 的序列）和复杂性（seqwish 在拓扑上更复杂）方面有所不同，但暗示了相同的结构。

![minigraph DRB1-3123](http://ekg.github.io/assets/DRB1-3123.minigraph.png)
![seqwish DRB1-3123](http://ekg.github.io/assets/DRB1-3123.seqwish.png)

我期待着与 Heng 和社区的其他成员合作，以了解构建和使用这些对象的最佳方式。
我相信这将需要竞争性但支持性和以参与者为导向的项目，在这些项目中，我们不仅构建，而且使用泛基因组来驱动基因组推断。
我们在一起，我相信我们会弄清楚如何实现这些想法并构建可扩展和有用的泛基因组方法。
重要的是我们学会读写相同的数据类型。
而且，至关重要的是，我希望这些数据类型足够通用，以支持研究人员可能需要做的所有事情。
激动人心的时刻即将来临！

#### 修正 (corrections)

2019-07-10：我已经更新了这篇文章，以正确区分 rGFA 和 minigraph，并解决变异图编辑图中的一个错误。
