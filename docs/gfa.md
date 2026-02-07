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
*   专注于 GFA 1.0 的高效处理和可视化。
*   核心格式 `.og`。
*   通常需要 GFA 1.0 输入（或者通过 `vg convert` 将 W 转为 P）。

### Gretl
*   Rust 编写的 GFA 统计评估工具。
*   需要数值型 Node ID 和拓扑排序。

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
