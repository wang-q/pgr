# nwk cut

`pgr nwk cut` 的目标是：给定一棵 Newick 树（系统发生树 / 层次聚类树），按照用户指定的规则把叶子节点切分成一组互不重叠的分组（partition），并以稳定、可复用的表格格式输出。

它关注的是“从树导出扁平聚类结果”，而不是“从数据构建树”。因此命名为 `cut` 比 `cluster` 更准确：树本身已经表达了层次结构，我们要做的是在树上选取一条切割规则并导出分组。

本文偏设计稿：描述 `pgr nwk cut` 的背景、输入输出约定、算法模式与选参思路，并对比相关生态工具。

## 适用场景

在实际分析中，经常会遇到这样的需求：

- 已经有一棵树（例如系统发生树、基于距离矩阵的层次聚类树、或某种推断得到的 dendrogram）。
- 希望在某个阈值下把叶子分组（比如得到“簇”用于下游统计、注释、画图、或与其他方法比较）。
- 切割规则可能不止一种：按高度切、按簇数切、按簇内最大两两距离（直径）切、要求每个簇必须是单系群（clade）、或者禁止跨越低支持度的边等。

`pgr nwk cut` 旨在提供一套与现有生态对齐但更“树友好”的切割方式：

- 对齐 R 的 `cutree()`（在 dendrogram 上切一刀得到分组）。
- 对齐 TreeCluster（在系统发生树上按生物学常用约束得到分组）。
- 与 `pgr clust` 区分：`pgr clust` 主要是从相似度/距离矩阵或图结构“构建聚类”；而 `pgr nwk cut` 是从“已有树”导出分组。

## 输入与输出

### 输入

- 输入树：Newick 格格式（建议与 `pgr nwk` 其它命令一致，支持文件或 stdin）。
- 分支长度：用于距离/高度相关方法（例如按 root distance、max pairwise distance 等）。
- 分支支持度（可选）：若树节点/边上携带支持度（例如 bootstrap），可作为“不可跨越”的约束条件。

### 输出

输出建议采用 TreeCluster 风格的 TSV（便于与既有工具互操作）：

```
SequenceName	ClusterNumber
leafA	        1
leafB	        1
leafC	        -1
...
```

- `SequenceName`：叶子标签（leaf label）。
- `ClusterNumber`：簇编号。
  - 非单例簇：从 1 开始递增编号。
  - 单例簇（singleton）：输出 `-1`（与 TreeCluster 一致，便于下游约定）。

也可以考虑提供替代输出格式（例如 `--format list|pair`），但默认 TSV 最通用。

## 核心概念：切树并导出 partition

不管采用何种规则，`cut` 的结果都可以理解为：

1. 在树上选择一组“切断点”（cut edges / cut nodes）。
2. 切断后，树被分成若干个互不相交的连通分量（component）。
3. 每个连通分量包含若干叶子；这些叶子构成一个输出簇。

不同算法的差异主要在于“切断点如何确定”。

## 选择阈值/簇数：扫描与准则

在 `cut` 场景里，用户常见的选择有两类：

- 直接指定簇数 `K`（类比 R `cutree(k=...)`）。
- 指定阈值 `t`（距离/高度/直径等），由阈值决定切割后的簇数。

当你不确定 `K` 或 `t` 应该取多少时，更稳妥的策略通常是“扫描 + 选点”，而不是一次性拍脑袋给出某个值。

### 扫描（scan）

建议提供（或在实现时优先考虑）一种扫描输出的能力：给定一组候选 `t` 或 `K`，对每个候选值计算并输出摘要指标，便于画曲线或人工挑选。

常见摘要指标包括：

- 簇数（总簇数 / 非单例簇数）
- 单例数量（singleton count）
- 最大簇大小、簇大小分布分位数
- 由支持度阈值导致的强制切断数量（若启用 `--support`）

直观用法是：先扫描得到一张表，再结合领域知识（例如期望簇规模、希望减少单例、或更保守地对待低支持区域）选择折中点。

### 选择准则（criterion）

如果需要自动选点，可以把“选择准则”显式做成一个可选项，并明确它与 `mclust` 这类模型选择的区别：

- `mclust` 的 BIC 依赖显式概率模型（高斯混合）与可计算的参数复杂度；`cut` 是在既定树上导出 partition，不天然对应同一个 BIC 语义。
- 因此在 `cut` 中，更合适的是提供若干“规则驱动/统计摘要驱动”的准则，例如：
  - 最大化非单例簇数（TreeCluster 的 `argmax_clusters` 属于这一类）
  - 最小化单例数量（在簇规模有意义时）
  - 约束最大簇大小/最小簇大小后再最大化某个目标（更贴近实践）

实现层面上，可以先从“扫描并输出摘要表”做起；自动选择可以作为扫描之上的薄层逻辑叠加，避免把一个难以解释的单一分数当作唯一答案。

一个可参考的现成思路是 TreeCluster 的无阈值模式 `-tf argmax_clusters`：对 `t'∈[0,t]` 的一组候选阈值运行同一种切割方法，选择“非单例簇（size>1）数量最多”的阈值作为输出。它本质上是把“扫描”内置化，再用一个简单可解释的准则做自动选点。

### 手肘规则（elbow）

手肘规则是一种常用的启发式：当你扫描一系列 `K`（或阈值 `t`）并计算某个指标时，曲线往往呈现“先明显改善，后收益递减”的形态；手肘点就是从“改善很快”过渡到“改善变慢”的拐点。

在 `pgr nwk cut` 的语境中，手肘规则更适合作为“扫描之后的人工选点方法”，而不是一个强约束的自动决策。

- **对 `K` 的手肘**：当命令支持 `--k <K>` 时，可以让用户扫描不同 `K`，再观察诸如单例比例、最大簇大小、或某个簇内距离摘要随 `K` 的变化趋势。
- **对阈值 `t` 的手肘**：在 TreeCluster 风格方法中更常见。随着 `t` 变大，切割会变“松”，簇数（尤其是非单例簇数）通常会快速下降并逐渐进入平台期。平台开始处常是一个实用的手肘点。

实践建议：

- 先扫描得到一张表（包含 `t/K`、簇数、非单例簇数、单例数、簇大小分布等）。
- 画出 `非单例簇数` 或 `单例数` 的曲线，优先找“平台开始处”，并结合业务期望（例如希望减少单例但不希望出现过大的超级簇）选择最终参数。

## 计划支持的模式与算法（设计稿）

这里列出与 TreeCluster / `cutree()` 对齐的常见模式，作为 `pgr nwk cut` 的设计目标。具体参数命名可在实现时再统一到 `pgr` 的命令风格中。

### 1) 按簇数切：`--k <K>`

等价于 R 的 `cutree(hc, k=K)`：

- 从根开始，逐步把某些内部节点“展开”，直到得到 K 个子树。
- 输出这 K 个子树的叶子集合。

直观解释：你不关心阈值是多少，只想要固定数量的分组。

### 2) 按高度/距离切：`--height <H>` 或 `--root-dist <D>`

等价于 R 的 `cutree(hc, h=H)`（对于 `hclust`，高度是合并高度；对于系统发生树，更自然的是 root distance 或 branch length）。

- 对系统发生树：可以定义为“从根向下累计分支长度，到达阈值就切断”。
- 对 dendrogram：可直接对应合并高度。

### 3) TreeCluster 风格：按簇内约束切（距离/直径/单系）

TreeCluster 的核心价值在于：除了“横向切一刀”，还支持“簇必须满足某种约束”。

典型方法包括：

- **max / max_clade**：每个簇内叶子两两距离的最大值（直径）不超过阈值。
  - `max_clade` 额外要求簇必须是 clade（即某个内部节点的整棵子树）。
- **avg_clade / med_clade**：簇内叶子两两距离的平均值/中位数不超过阈值，并要求 clade。
- **single_linkage_cut / single_linkage_union**：单链接思想在树上的变体。
- **length / length_clade**：切断所有长度大于阈值的边（或约束每个簇内最大边长）。
- **root_dist / leaf_dist_{min,max,avg}**：以根到叶的距离统计量为基准切树。

这些方法的共同点：仍然输出叶子的分组，但“切断策略”不再是单纯的水平线，而是根据树拓扑和距离统计量自适应决定。

### 4) 支持度约束：`--support <S>`

TreeCluster 的做法是：当某条边（或节点）支持度低于阈值时，视为“不可跨越”，相当于强制切断（或者把该边长度视为无限大）。

在 `pgr nwk cut` 中也建议提供类似选项：

- 使聚类结果对低支持区域更保守。
- 对下游解释更友好（不会把低支持连接当作可靠证据）。

## 与相关工具的关系与区别

### 与 R `hclust + cutree()`

- **相同点**：都是“树 → 叶子分组”。
- **不同点**：
  - `cutree()` 面向的是 `hclust` 产生的 dendrogram；`pgr nwk cut` 面向 Newick 树（系统发生树或一般 Newick）。
  - `pgr nwk cut` 计划支持 TreeCluster 风格的生物树约束（clade、支持度阈值、树上单链接等），这超出了 `cutree()` 的常规用法。

### 与 TreeCluster

- **相同点**：目标与输出格式高度一致（叶子 → 簇，单例为 -1 的约定也可沿用）。
- **不同点**：
  - TreeCluster 是独立工具；`pgr nwk cut` 将融入 `pgr` 的 Newick 工具链，便于与 `pgr nwk prune/reroot/subtree/...` 串联。
  - `pgr` 复用自己的树数据结构与遍历/查询能力，便于后续扩展（例如与 `pgr nwk stat/distance/topo` 联动）。

### 与 `pgr clust`

- **`pgr clust`**（例如 [clust/mod.rs](file:///c:/Users/wangq/Scripts/pgr/src/cmd_pgr/clust/mod.rs)、[mcl.rs](file:///c:/Users/wangq/Scripts/pgr/src/cmd_pgr/clust/mcl.rs)）更偏“从数据/图出发做聚类”：
  - 输入往往是相似度/距离表或图。
  - 输出是簇或簇对关系。
- **`pgr nwk cut`** 是“从树出发导出分组”：
  - 输入是 Newick 树。
  - 输出是叶子分组。

两者可以互补：例如先用某种方法构树，再用 `nwk cut` 导出不同阈值下的分组；或用 `clust` 在图上聚类后与树切割结果对比。

## 典型用法（建议）

下面示例为“拟定接口”，最终以实现时的参数为准：

```bash
# 1) 固定簇数（类比 R cutree(k=...)）
pgr nwk cut tree.nwk --k 20 > clusters.tsv

# 2) 按 root distance 切割
pgr nwk cut tree.nwk --root-dist 0.03 > clusters.tsv

# 3) TreeCluster 风格（max_clade），并启用支持度阈值
pgr nwk cut tree.nwk --method max-clade -t 0.045 --support 70 > clusters.tsv
```

## 与 `pgr nwk` 工具链的协作

`pgr nwk` 已包含大量对树的操作与分析（见 [nwk/mod.rs](file:///c:/Users/wangq/Scripts/pgr/src/cmd_pgr/nwk/mod.rs)）。`cut` 的定位是 “analysis/ops 的桥梁”：

- 你可以先用 `reroot/prune/subtree/replace/...` 规范化树，再 `cut` 导出分组。
- 也可以先 `cut` 得到分组，再回到 `nwk` 其它命令进行统计、对比或可视化。

## 工作流与工具链协作

为了保持命令的专注与正交性，我们推荐以下“生成-评估”分离的工作流：

### 1. 生成 (Generation)

使用 `pgr nwk cut` 专注于从树生成分组（Partition）。

- 它只负责“切”，不负责“评”。
- 支持多种策略（k, height, max_clade 等）和参数扫描。
- 输出标准 TSV 格式。

### 2. 评估 (Evaluation)

评估聚类质量通常需要参考标准（Ground Truth）或与其他聚类结果对比。这部分逻辑放入独立的 `pgr clust` 或 `pgr nwk` 命令中：

- **通用指标 (`pgr clust eval` / `compare`)**：
  - 输入：两个聚类结果 TSV（或一个结果 + 一个参考）。
  - 输出：ARI (Adjusted Rand Index), AMI (Adjusted Mutual Information), V-Measure, Fowlkes-Mallows 等。
  - 适用场景：当你已知样本的真实分类，或者想比较两种切割参数的差异度时。

- **树相关指标 (`pgr nwk metrics`)**：
  - 输入：树文件 + 聚类结果。
  - 输出：Parsimony score, Likelihood (需配合序列), Silhouette score (基于树上距离矩阵) 等。
  - 适用场景：没有真实分类，需要评估聚类在树结构上的紧密性或分离度。

### 推荐工作流示例

#### 1. 经典系统发育分析
```bash
# 1. 扫描不同参数，生成多个聚类结果
pgr nwk cut input.nwk --method max-clade --scan 0.01,0.05,0.10 > partitions.tsv

# 2. (可选) 如果有真实分类 metadata.tsv，评估哪个阈值最好
pgr clust eval partitions.tsv metadata.tsv --metric ari

# 3. 选定最佳阈值，生成最终聚类
pgr nwk cut input.nwk --method max-clade -t 0.05 > final_cluster.tsv

# 4. 可视化或提取子树
pgr nwk subset input.nwk --list final_cluster.tsv --cluster-id 1 > cluster1.nwk
```

#### 2. 层次聚类（hclust）接入
从距离矩阵出发，经由 hclust 生成树，再进行切分与评估（参见 [hclust.md](file:///c:/Users/wangq/Scripts/pgr/docs/hclust.md)）。

```bash
# 1. 准备 PHYLIP 距离矩阵 (若有 pair TSV，先转为 phylip)
pgr mat to-phylip pairs.tsv -o matrix.phy

# 2. 生成层次聚类树 (Ward 方法，启用叶序优化)
pgr mat hclust matrix.phy --method ward --optimal-ordering > tree.nwk

# 3. 切分 (按高度阈值切，或按 K 切)
pgr nwk cut tree.nwk --height 0.05 > clusters.tsv

# 4. 评估 (计算 Cophenetic 相关系数与 Silhouette)
pgr nwk metrics tree.nwk --metrics cophenet --dist matrix.phy > fit.tsv
pgr nwk metrics tree.nwk --part clusters.tsv --metrics silhouette > sil.tsv
```
