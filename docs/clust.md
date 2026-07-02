# pgr clust - Clustering Algorithms

## 概述 (Overview)

`pgr clust` 模块提供了一系列用于序列、基因组特征和一般数据的聚类算法。这些工具旨在处理生物信息学中常见的距离矩阵、相似性网络和特征向量。

命令按输入数据类型分为两类：
1.  **Tree-based**: 基于距离矩阵构建系统发育树或层级结构 (`hier`, `nj`, `upgma`)。
2.  **Flat**: 基于图或向量直接生成分组 (`cc`, `mcl`, `k-medoids`, `dbscan`)。

## 算法列表 (Algorithms)

### MCL (Markov Cluster Algorithm)

- **原理**：通过在图上模拟随机游走（Random Walk），通过交替执行“扩展（Expansion）”和“膨胀（Inflation）”操作，使强连接区域内的流更加集中，弱连接区域的流逐渐消失，最终自然分割出模块。
- **命令**：`pgr clust mcl`
- **特点**：基于流模拟的图聚类。
- **适用场景**：**生物网络**（如 SSN）、蛋白家族检测、模块发现。
- **优势**：对噪声鲁棒，能处理复杂的网络结构。

### Connected Components (CC)

- **原理**：图论中的基础概念，寻找图中所有互相可达的节点集合。相当于在给定阈值下，将所有直接或间接相连的样本归为一个簇。
- **命令**：`pgr clust cc`
- **特点**：最基础的连通分量。
- **适用场景**：极高相似度阈值下的快速去重。
- **优势**：极快（线性复杂度）。

### K-Medoids

- **原理**：类似 K-Means 的迭代优化，但中心点（Medoid）必须是数据集中的实际样本。通过最小化所有点到其最近中心点的距离之和（Dissimilarity）来更新中心。
- **命令**：`pgr clust k-medoids`
- **特点**：类似 K-Means，但中心点必须是实际样本（Medoid）。
- **适用场景**：抗噪声需求，或仅有**距离矩阵**（非欧氏空间）的情况。
- **优势**：对异常值鲁棒，结果可解释（中心是真实样本）。

### DBSCAN

- **原理**：基于密度的聚类。从任意点出发，若其 $\epsilon$ 邻域内的点数超过 `min_samples`，则形成核心点并扩展簇；密度不足的区域被视为噪声。
- **命令**：`pgr clust dbscan`
- **特点**：基于密度的聚类，需指定邻域半径 `eps` 和最小点数 `min_samples`。
- **适用场景**：**非凸形状**簇，密度不均匀分布，需**离群点检测**。
- **优势**：不需要指定簇数 K，能识别噪声。
- **文档**：[clust-dbscan.md](clust-dbscan.md)

### UPGMA

- **原理**：非加权组平均法。一种自底向上的层次聚类，每次合并距离最近的两个簇，新簇与其他簇的距离计算为所有成员间距离的算术平均。假设进化速率恒定（分子钟）。
- **命令**：`pgr clust upgma`
- **特点**：层次聚类（平均链接），输出**有根树**。
- **适用场景**：假设**分子钟**（超度量）的系统发育分析。
- **优势**：生成层级结构，树高有明确的距离意义。

### NJ (Neighbor-Joining)

- **原理**：邻接法。通过最小化树的总枝长（基于 Q 矩阵校正距离），迭代地合并“净发散度”最小的节点对。不假设分子钟，允许不同分支有不同的进化速率。
- **命令**：`pgr clust nj`
- **特点**：基于距离矩阵的构树算法，输出**无根树**。
- **适用场景**：一般加性距离（无需分子钟假设），构建进化树。
- **优势**：速度快，对不同演化速率鲁棒。

### Hierarchical Clustering

- **原理**：通用的自底向上（Agglomerative）聚类框架。通过不同的链接准则（如 Ward 最小方差、Complete 最远距离）合并簇，构建完整的树状层级。
- **命令**：`pgr clust hier` (别名 `hclust`)
- **特点**：通用层次聚类，支持 `single`, `complete`, `average`, `ward` 等方法。
- **实现状态**：已实现（$O(N^2)$ NN-Chain 优化）。
- **价值**：提供通用层级结构视图（不限于生物演化），配合 `clust cut` 可灵活获取不同粒度的分组。
- **文档**：[clust-hier.md](clust-hier.md)

### Tree Cut

- **原理**：对已有的 Newick 树（系统发生树或层次聚类树），按指定规则切分为扁平的聚类分组（Partition）。支持按簇数 (`--k`)、高度 (`--height`)、簇内直径 (`--max-clade`)、动态切割 (`--dynamic-tree`/`--dynamic-hybrid`) 等多种切割策略。
- **命令**：`pgr clust cut`
- **特点**：从已有树导出分组，不重建聚类；支持参数扫描 (`--scan`) 与代表点选择 (`--rep`)。
- **适用场景**：已有树结构（来自 `clust hier`、`clust upgma`、`clust nj` 或外部工具），需要在不同阈值下切分并评估。
- **文档**：[clust-cut.md](clust-cut.md)

## 评估与分析 (Evaluation)

这些命令不产生聚类，而是评估聚类或树的质量。

- **Tree-based Evaluation**
  - **命令**：`pgr nwk eval` (计划中)
  - **定位**：树结构的多维度评估。
  - **功能**：几何紧密性（Silhouette）、分类纯度（Purity）、演化一致性（Discordance）。

- **Partition-based Evaluation**
  - **命令**：`pgr clust eval`
  - **定位**：通用聚类质量评估（支持有/无 Ground Truth）。
  - **功能**：ARI, AMI, V-Measure (外部); Silhouette, Davies-Bouldin (内部)。
  - **文档**：[clust-eval.md](clust-eval.md)

## 计划中 (Planned)

GMM、HDBSCAN、Louvain/Leiden 等算法已列入路线图。

## 不建议实现 / 暂无计划 (Not Planned)

这些算法虽然经典，但在生物大数据场景下存在局限性，暂不作为核心功能引入。

- **K-Means**
  - **原因**：虽然速度快，但假设簇是球形且方差相等，且质心（Centroid）通常不是真实的样本点，缺乏生物学解释性（如无法直接作为代表序列）。
  - **替代**：`K-Medoids`（已实现），其中心点（Medoid）必须是真实样本，且支持任意距离矩阵，更适合生物序列分析。

- **Bisecting K-Means**
  - **原理**：自顶向下的分裂式聚类。初始将所有点视为一簇，每次选择 SSE 最大的簇进行二分 K-Means 分裂，直到达到 K 个簇。
  - **原因**：虽然能生成树状结构（二叉树），但继承了 K-Means 的局限性（需欧氏距离、质心非真实样本）。生物学构树通常偏好自底向上的 Agglomerative 方法（如 UPGMA/NJ）。

- **Affinity Propagation (AP)**
  - **原理**：基于消息传递机制，所有点相互竞争成为代表点（Exemplar）。不需要指定簇数，但计算复杂度高。
  - **原因**：时间与空间复杂度较高 ($O(N^2)$)，难以处理大规模生物序列数据（如 >1万条序列）。
  - **替代**：对于小规模数据寻找代表点，推荐使用 `K-Medoids`；对于自动确定簇数，推荐 `DBSCAN` 或 `MCL`。

- **Spectral Clustering (谱聚类)**
  - **原理**：利用拉普拉斯矩阵的特征向量进行降维，然后在低维空间进行 K-Means 聚类。本质上是寻找图的最小归一化割（Normalized Cut）。
  - **原因**：这就涉及构建拉普拉斯矩阵并进行特征分解，计算开销大 ($O(N^3)$)。
  - **替代**：`MCL` 在生物网络聚类中通常能提供类似甚至更好的效果，且扩展性更好。

- **Mean Shift**
  - **原理**：基于密度的爬山算法。通过不断将点移动到其邻域的密度中心（均值漂移），最终收敛到局部密度峰值（模态）。
  - **原因**：计算复杂度高，且带宽参数（bandwidth）难以自适应选择。
  - **替代**：`DBSCAN` 或 `GMM` 通常能覆盖其密度估计的需求。

- **OPTICS**
  - **原理**：通过生成一个可达距离图（Reachability Plot），对数据点进行排序，从而在一次运行中捕获所有可能的密度层级。解决了 DBSCAN 对全局 `eps` 敏感的问题。
  - **原因**：其核心思想（层级密度聚类）已被 **HDBSCAN** 更好地继承和自动化；OPTICS 的结果（可达距离图）需要复杂的后处理才能得到明确的簇。
  - **替代**：推荐使用更现代、参数更少且自动化程度更高的 `HDBSCAN`。

- **Biclustering (双聚类)**
  - **原因**：同时对行和列进行聚类（如 Spectral Co-Clustering），主要用于基因表达谱分析等特定矩阵子块挖掘场景，与 `pgr` 当前专注的“样本分组”目标差异较大。
  - **替代**：若需对特征（列）进行聚类，可转置矩阵后使用标准聚类；若需寻找共表达模块，建议使用专门的表达谱分析工具（如 WGCNA）。

- **BIRCH**
  - **原理**：基于聚类特征树（CF Tree）的增量聚类。通过单次扫描构建一棵高度压缩的树，树节点存储簇的统计摘要（Sum, SquareSum），极适合超大规模数据集。
  - **原因**：强依赖于欧氏空间的统计特性（计算质心和半径），不适用于生物序列的复杂距离度量（如 Edit Distance）；且对簇形状有限制。
  - **替代**：对于大规模向量，`K-Means (MiniBatch)` 是更通用的选择；对于大规模序列，推荐 `MCL`（图聚类）或 `CD-HIT/MMseqs2`（贪心聚类）。

## 算法详细说明 (Detailed Descriptions)

### GMM (Gaussian Mixture Models) [计划中]

引入 GMM 的动机：
- **软聚类 (Soft Clustering)**：不同于 K-means 的硬划分，GMM 给出样本属于某簇的概率，适合处理边界模糊的生物学分类（如亚种、基因家族过渡态）。
- **非球形簇**：通过协方差矩阵建模，适应不同形状和大小的簇（K-means 假设簇是等方差球形）。
- **生成式模型**：可用于密度估计（Density Estimation）和异常检测（Outlier Detection）。

**规划接口**：
```bash
# 从 CSV/TSV 向量输入进行 GMM 聚类
pgr clust gmm input.tsv --k 5 --cov full > clusters.tsv

# 输出包含：ID, Cluster, PosteriorProb (后验概率)
```

### 模型选择 (Model Selection)

如何确定聚类的簇数（K）或最佳模型复杂度？

- **BIC (Bayesian Information Criterion)** [计划中]：
  - 在 GMM 中，BIC 权衡了对数似然（拟合度）与参数数量（复杂度）。
  - `pgr` 可提供 `clust gmm --scan-k 2..20`，自动计算并输出 BIC 曲线，辅助用户选择最佳 K（通常是 BIC 最低点或手肘点）。
- **Silhouette / Calinski-Harabasz** [部分支持]：基于几何距离的评估指标，适用于 K-means 或一般距离聚类（`clust eval` 已支持距离矩阵版 Silhouette；树上 Silhouette 计划在 `pgr nwk eval` [计划中] 中实现）。

## 大规模数据策略 (Two-stage / Representative Strategy)

对于 $N > 20,000$ 的大规模数据，全连接层次聚类的内存 ($O(N^2)$) 和计算 ($O(N^2)$) 开销急剧增加。

**内存估算 (f32 Condensed Matrix)**:
- **1 GiB**: ~23,000 点
- **10 GiB**: ~73,000 点
- **32 GiB**: ~130,000 点
- **64 GiB**: ~185,000 点

**结论**: 即使在 64G 内存的高配服务器上，处理 $N=200k$ 也已接近极限。

**推荐策略**：采用“两步法”，结合快速聚类与精细构树。
1.  **预聚类/压缩**: 使用线性或近线性算法（如 `pgr clust k-medoids`、`pgr clust mcl` 或外部工具 `mmseqs2`）将数据压缩为 $K$ 个代表点（$K \approx 5000 \sim 10000$）。
2.  **层次聚类**: 提取代表点之间的距离矩阵，运行 `pgr clust hier` 构建骨架树。

**工作流示例**:
```bash
# 1. 快速聚类选出代表点 (k=5000)
pgr clust k-medoids all_data.tsv --k 5000 --format pair > clusters.tsv

# 2. 提取代表点列表
cut -f1 clusters.tsv | sort -u > representatives.list

# 3. 提取代表点的子矩阵
pgr mat subset all_data.tsv --list representatives.list -o sub_matrix.phy

# 4. 对代表点构树
pgr clust hier sub_matrix.phy --method ward > backbone.nwk
```

## 推荐工作流

### 场景 A：蛋白家族挖掘 (Graph-based)

```bash
# 1. 序列比对构建网络 (如 mmseqs/blast -> pair.tsv)
# 2. MCL 聚类
pgr clust mcl pairs.tsv --inflation 2.0 > families.tsv
```

### 场景 B: 层次聚类参数扫描与评估 (Workflow)

结合 `clust cut` 的扫描能力与 `clust eval` 的批量评估，寻找最佳切分阈值。

```bash
# 1. 生成层次聚类树
pgr clust hier matrix.phy --method ward > tree.nwk

# 2. 扫描阈值并评估内部指标 (Silhouette)
# clust cut 在 scan 模式下输出长表，直接传给 clust eval
pgr clust cut tree.nwk --height 1.0 --scan 0,1.0,0.05 | \
    pgr clust eval - --format long --matrix matrix.phy > evaluation.tsv

# 3. 分析 evaluation.tsv 选择最佳阈值 (如 Silhouette 最大处)
# 假设最佳阈值为 0.45
pgr clust cut tree.nwk --height 0.45 > final_clusters.tsv
```

## 输入输出格式约定 (File Formats)

`pgr clust` 系列命令涉及多种数据格式，为了便于与其他工具交互，约定如下标准格式。

### 1. 分区文件 (Partition File)

用于表示聚类结果（即样本到簇的映射）。支持三种格式，通过 `--format` 参数指定。

#### Pair Format (默认, `--format pair`)
最通用的长表格式，每行表示一个样本所属的簇。
- **结构**：`ClusterID <tab> Item`
- **特点**：易于解析，支持流式处理。
- **示例**：
  ```text
  # Numeric ID
  1	GeneA
  1	GeneB
  2	GeneC

  # Representative as ID
  GeneA	GeneA
  GeneA	GeneB
  GeneC	GeneC
  ```

#### Cluster Format (`--format cluster`)
宽表格式，每行代表一个簇，包含该簇的所有成员。
- **结构**：`Item1 <space/tab> Item2 ...`
- **特点**：人类可读性好，适合查看聚类结果。行号（1-based）即为 ClusterID。
- **示例**：
  ```text
  GeneA GeneB
  GeneC
  ```

#### Long Format (Batch, `--format long`)
用于批量评估的专用格式，通常由 `pgr clust cut --scan` 生成。
- **结构**：`Group <tab> ClusterID <tab> Item`
- **Group 列**：用于标识不同的参数组合或切割方法。格式通常为 `Method=Value`（如 `height=0.5`）。
  - `pgr clust eval` 会保留此列作为评估结果的标识符。
- **示例**：
  ```text
  height=0.1	1	GeneA
  height=0.1	2	GeneC
  height=0.2	1	GeneA
  height=0.2	1	GeneC
  ```

### 2. 距离矩阵 (Distance Matrix)

用于 `clust hier`, `nj`, `upgma` 以及 `eval --matrix`。

#### PHYLIP Format (Relaxed)
- **结构**：
  - 第一行：样本数量 $N$。
  - 后续 $N$ 行：`Name <space> Dist1 <space> Dist2 ...`
- **特点**：标准生物信息学格式。`pgr` 支持“宽松”格式（名字与数据间可有任意空白）。
- **示例**：
  ```text
  3
  A  0.0 0.1 0.5
  B  0.1 0.0 0.5
  C  0.5 0.5 0.0
  ```

### 3. 坐标/特征向量 (Coordinates / Feature Vector)

用于 `clust eval --coords` (Davies-Bouldin Index) 或未来可能的 `kmeans/gmm`。

#### FeatureVector Format
- **结构**：`Name <tab> Val1,Val2,Val3...`
- **分隔符**：名字与向量间用 **制表符**，向量数值间用 **逗号**。
- **示例**：
  ```text
  GeneA	1.2,0.5,3.3
  GeneB	1.1,0.6,3.1
  ```
- **兼容性**：此格式与 `pgr dist vector` 的输出一致。
