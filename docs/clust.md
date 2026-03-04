# clust

`pgr clust` 模块专注于“从数据/图出发构建聚类”，与 `pgr nwk cut`（从树导出聚类）互补。

本文档规划 `pgr clust` 的核心能力、适用场景、以及引入 GMM（高斯混合模型）等高级统计聚类方法的思路，并探讨模型选择（Model Selection）在生物序列分析中的应用。

## 核心定位

- **输入**：相似度矩阵、距离矩阵、图结构（边列表）、或高维向量（embeddings）。
- **输出**：
  - **扁平分组 (Flat Partition)**：如 `clusters.tsv`，每行一个样本及其所属簇 ID。
  - **树状结构 (Hierarchy/Tree)**：如 `tree.nwk`，通过 Newick 格式表达样本间的层级关系（仅 `hier`/`hdbscan` 支持）。
- **目标**：发现数据中的自然结构（Structure Discovery），不依赖预先存在的树。

## 算法概览与实现状态 (Overview & Status)

本节汇总 `pgr` 支持或规划中的聚类算法，按实现状态分类说明。

### ✅ 已实现 (Implemented)

这些算法已在 `pgr clust` 中可用，适合直接投入生产流程。

- **MCL (Markov Cluster Algorithm)**
  - **原理**：通过在图上模拟随机游走（Random Walk），通过交替执行“扩展（Expansion）”和“膨胀（Inflation）”操作，使强连接区域内的流更加集中，弱连接区域的流逐渐消失，最终自然分割出模块。
  - **命令**：`pgr clust mcl`
  - **特点**：基于流模拟的图聚类。
  - **适用场景**：**生物网络**（如 SSN）、蛋白家族检测、模块发现。
  - **优势**：对噪声鲁棒，能处理复杂的网络结构。

- **Connected Components (CC)**
  - **原理**：图论中的基础概念，寻找图中所有互相可达的节点集合。相当于在给定阈值下，将所有直接或间接相连的样本归为一个簇。
  - **命令**：`pgr clust cc`
  - **特点**：最基础的连通分量。
  - **适用场景**：极高相似度阈值下的快速去重。
  - **优势**：极快（线性复杂度）。

- **K-Medoids**
  - **原理**：类似 K-Means 的迭代优化，但中心点（Medoid）必须是数据集中的实际样本。通过最小化所有点到其最近中心点的距离之和（Dissimilarity）来更新中心。
  - **命令**：`pgr clust k-medoids`
  - **特点**：类似 K-Means，但中心点必须是实际样本（Medoid）。
  - **适用场景**：抗噪声需求，或仅有**距离矩阵**（非欧氏空间）的情况。
  - **优势**：对异常值鲁棒，结果可解释（中心是真实样本）。

- **DBSCAN**
  - **原理**：基于密度的聚类。从任意点出发，若其 $\epsilon$ 邻域内的点数超过 `min_samples`，则形成核心点并扩展簇；密度不足的区域被视为噪声。
  - **命令**：`pgr clust dbscan`
  - **特点**：基于密度的聚类，需指定邻域半径 `eps` 和最小点数 `min_samples`。
  - **适用场景**：**非凸形状**簇，密度不均匀分布，需**离群点检测**。
  - **优势**：不需要指定簇数 K，能识别噪声。

- **UPGMA**
  - **原理**：非加权组平均法。一种自底向上的层次聚类，每次合并距离最近的两个簇，新簇与其他簇的距离计算为所有成员间距离的算术平均。假设进化速率恒定（分子钟）。
  - **命令**：`pgr clust upgma`
  - **特点**：层次聚类（平均链接），输出**有根树**。
  - **适用场景**：假设**分子钟**（超度量）的系统发育分析。
  - **优势**：生成层级结构，树高有明确的距离意义。

- **NJ (Neighbor-Joining)**
  - **原理**：邻接法。通过最小化树的总枝长（基于 Q 矩阵校正距离），迭代地合并“净发散度”最小的节点对。不假设分子钟，允许不同分支有不同的进化速率。
  - **命令**：`pgr clust nj`
  - **特点**：基于距离矩阵的构树算法，输出**无根树**。
  - **适用场景**：一般加性距离（无需分子钟假设），构建进化树。
  - **优势**：速度快，对不同演化速率鲁棒。

### 📅 计划中 (Planned)

这些算法已列入路线图，旨在补全统计聚类与大规模向量分析能力。

- **Hierarchical Clustering**
  - **原理**：通用的自底向上（Agglomerative）聚类框架。通过不同的链接准则（如 Ward 最小方差、Complete 最远距离）合并簇，构建完整的树状层级。
  - **命令**：`pgr clust hier` (别名 `hclust`)
  - **计划内容**：支持 `ward`, `complete` 等 linkage 方法，输出 Newick 树。
  - **价值**：提供通用层级结构视图（不限于生物演化），配合 `nwk cut` 可灵活获取不同粒度的分组。

- **GMM (Gaussian Mixture Models)**
  - **原理**：假设数据由 $K$ 个高斯分布混合而成。使用 EM（期望最大化）算法迭代估计每个高斯分量的参数（均值、协方差）及每个样本属于各分量的后验概率。
  - **命令**：`pgr clust gmm`
  - **计划内容**：支持**软聚类**（概率输出）与 **BIC** 模型选择。
  - **价值**：适合**椭球形簇**与密度估计，解决 K-Means 仅适应球形簇的限制；BIC 可辅助确定最佳 K。

- **HDBSCAN**
  - **原理**：结合层次聚类与 DBSCAN。通过构建基于密度的层次树（Condensed Tree），并根据簇的稳定性（Stability）在不同层级自动提取最佳簇，无需全局 $\epsilon$。
  - **命令**：`pgr clust hdbscan`
  - **scikit-learn 对应**：`HDBSCAN`
  - **计划内容**：层次化 DBSCAN，无需手动指定 `eps`。
  - **价值**：DBSCAN 的现代高级版，**自动适应不同密度的簇**，参数更少且更稳健。

- **Louvain / Leiden**
  - **原理**：基于模块度（Modularity）优化的社区发现算法。Louvain 贪心地最大化模块度；Leiden 改进了 Louvain 的局部合并策略，保证连通性并加速收敛。
  - **命令**：(待定)
  - **计划内容**：社区发现算法。
  - **价值**：比 MCL 更适合**超大规模网络**的层次化结构探索。

### 🚫 不建议实现 / 暂无计划 (Not Planned)

这些算法虽然经典，但在生物大数据场景下存在局限性，暂不作为核心功能引入。

- **K-Means**
  - **原因**：虽然速度快，但假设簇是球形且方差相等，且质心（Centroid）通常不是真实的样本点，缺乏生物学解释性（如无法直接作为代表序列）。
  - **替代**：`K-Medoids`（已实现），其中心点（Medoid）必须是真实样本，且支持任意距离矩阵，更适合生物序列分析。

- **Bisecting K-Means**
  - **原理**：自顶向下的分裂式聚类。初始将所有点视为一簇，每次选择 SSE 最大的簇进行二分 K-Means 分裂，直到达到 K 个簇。
  - **原因**：虽然能生成树状结构（二叉树），但继承了 K-Means 的局限性（需欧氏距离、质心非真实样本）。生物学构树通常偏好自底向上的 Agglomerative 方法（如 UPGMA/NJ）。
  - **替代**：`clust hier` (Hierarchical Clustering) 或 `HDBSCAN`。

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

### GMM (Gaussian Mixture Models) [规划中]

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

- **BIC (Bayesian Information Criterion)** [规划中]：
  - 在 GMM 中，BIC 权衡了对数似然（拟合度）与参数数量（复杂度）。
  - `pgr` 可提供 `clust gmm --scan-k 2..20`，自动计算并输出 BIC 曲线，辅助用户选择最佳 K（通常是 BIC 最低点或手肘点）。
- **Silhouette / Calinski-Harabasz** [部分支持]：基于几何距离的评估指标，适用于 K-means 或一般距离聚类（`nwk metrics` 已支持树上 Silhouette）。

## 实现分析与对比 (Implementation Analysis)

通过与 scikit-learn 的源码对比，总结 `pgr` 当前实现的特点与未来优化方向。

### DBSCAN

- **scikit-learn 实现**：
  - **核心**：使用 `NearestNeighbors` 模块（基于 BallTree/KDTree）加速邻域查询。
  - **优化**：支持稀疏矩阵；通过 `n_jobs` 并行化；核心逻辑部分使用 Cython 加速。
  - **适用性**：能处理数百万量级的数据（如果维度不高）。
- **pgr 实现**：
  - **核心**：基于 `ScoringMatrix` 的朴素 $O(N^2)$ 距离遍历；`region_query` 为线性扫描。
  - **特点**：代码简洁，无需额外空间索引库；输出包含“代表点对”等生物学便利功能。
  - **局限**：缺乏空间索引，在大规模（>1万点）或高维数据上性能不如 sklearn。
- **未来方向**：对于大规模向量输入，需引入空间索引（如 R-tree/KD-tree）或并行化邻域搜索。

### 层次聚类 (UPGMA / NJ vs Agglomerative)

- **scikit-learn (AgglomerativeClustering)**：
  - **定位**：通用统计聚类，输出 Linkage Matrix（$N-1$ 次合并记录）。
  - **优化**：
    - 使用 **MST (最小生成树)** 加速 Single Linkage ($O(N^2)$)。
    - 使用 **Heap (堆)** 结构加速 Ward/Average/Complete Linkage 的最近邻查找。
  - **输出**：不直接生成 Newick 树，需转换。
- **pgr (UPGMA / NJ)**：
  - **定位**：生物系统发育构树，直接输出 **Tree** 对象和 **Newick** 格式。
  - **实现**：
    - **UPGMA**：动态维护距离矩阵（HashMap），每次迭代寻找最小值，复杂度 $O(N^3)$。
    - **NJ**：经典的 Neighbor-Joining 实现，计算净发散度与 Q 矩阵，输出无根树。
  - **优势**：逻辑直观，原生支持生物学所需的枝长计算与树操作。
  - **局限**：未采用 Heap 优化，在大规模数据（>5000 序列）上速度慢于优化过的 Linkage 算法。
- **未来方向**：实现通用的 `clust hier` 时，应参考 sklearn 的 Heap 优化思路，而非简单复用 UPGMA 的 $O(N^3)$ 逻辑。

## 推荐工作流

### 场景 A：蛋白家族挖掘 (Graph-based)

```bash
# 1. 序列比对构建网络 (如 mmseqs/blast -> pair.tsv)
# 2. MCL 聚类
pgr clust mcl pairs.tsv --inflation 2.0 > families.tsv
```

### 场景 B：基于 Embedding 的亚型分类 (Vector-based)

```bash
# 1. 计算序列 embedding (如 k-mer profile 或 ESM) -> vectors.tsv
# 2. GMM 聚类并扫描最佳 K (基于 BIC)
pgr clust gmm vectors.tsv --scan-k 2..15 --cov diag > bic_scores.tsv

# 3. 选定 K=8 进行聚类
pgr clust gmm vectors.tsv --k 8 --out-prob > soft_clusters.tsv
```

## 实现路线图

1. **基础图聚类**：已完成 MCL、CC、DBSCAN、K-Medoids。
2. **系统发育构树**：已完成 UPGMA、NJ。
3. **向量支持**：建立读取稠密向量/矩阵的基础设施（进行中）。
4. **统计聚类**：引入 GMM 实现，支持 BIC 模型选择（规划中）。
5. **层次聚类扩展**：实现通用的 Hierarchical Clustering (hclust) 和 HDBSCAN（规划中）。
