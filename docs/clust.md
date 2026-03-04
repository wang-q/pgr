# clust

`pgr clust` 模块专注于“从数据/图出发构建聚类”，与 `pgr nwk cut`（从树导出聚类）互补。

本文档规划 `pgr clust` 的核心能力、适用场景、以及引入 GMM（高斯混合模型）等高级统计聚类方法的思路，并探讨模型选择（Model Selection）在生物序列分析中的应用。

## 核心定位

- **输入**：相似度矩阵、距离矩阵、图结构（边列表）、或高维向量（embeddings）。
- **输出**：样本的分组（Partition）或软聚类概率（Soft Clustering）。
- **目标**：发现数据中的自然结构（Structure Discovery），不依赖预先存在的树。

## 功能规划

## 功能概览

### 1. 基于图的聚类 (Graph Clustering)

适用于稀疏相似度网络（如 sequence similarity network, SSN）。

- **MCL (Markov Cluster Algorithm)** [已实现]：基于流模拟的图聚类，适合发现紧密连接的模块（蛋白家族检测标准方法）。
  - 命令：`pgr clust mcl`
- **Connected Components (CC)** [已实现]：最基础的连通分量，适合极高阈值下的快速去重。
  - 命令：`pgr clust cc`
- **Louvain / Leiden** [规划中]：社区发现算法，适合大规模网络的层次化结构探索。

### 2. 基于距离/向量的统计聚类 (Statistical Clustering)

当数据以连续坐标（如 k-mer profile, embedding）或稠密距离矩阵存在时，统计聚类提供基于分布的建模能力。

- **Hierarchical Clustering** [规划中]：层次聚类，支持 Ward/Average/Complete 等方法，输出 dendrogram 树。
  - 命令：`pgr clust hier`（别名 `hclust`）
  - 详情：参见 [clust-hier.md](file:///c:/Users/wangq/Scripts/pgr/docs/clust-hier.md)
- **K-Medoids** [已实现]：类似于 K-Means，但中心点必须是实际样本（Medoid），对噪声更鲁棒，支持任意距离矩阵。
  - 命令：`pgr clust k-medoids`
- **DBSCAN** [已实现]：基于密度的聚类，能发现任意形状的簇并识别噪声点。
  - 命令：`pgr clust dbscan`
- **K-Means** [规划中]：最基础的欧氏空间聚类，适合大规模向量数据。

#### GMM (Gaussian Mixture Models) [规划中]

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

### 3. 模型选择 (Model Selection)

如何确定聚类的簇数（K）或最佳模型复杂度？

- **BIC (Bayesian Information Criterion)** [规划中]：
  - 在 GMM 中，BIC 权衡了对数似然（拟合度）与参数数量（复杂度），是选择 GMM 组件数的标准方法。
  - 对于 K-Means，也可通过假设簇内分布为高斯来近似计算 BIC（X-Means 算法思路）。
  - **注意**：BIC 不直接适用于 K-Medoids（因其缺乏明确的概率似然函数），除非强行假设簇内分布。
- **Silhouette / Calinski-Harabasz** [部分支持]：基于几何距离的评估指标。
  - **推荐用于 K-Medoids**：对于 K-Medoids，建议优先使用 Silhouette 系数或 Elbow Method（手肘法）来评估 K，而非 BIC。
  - `nwk metrics` 已支持树上 Silhouette；后续可扩展至向量/距离矩阵的 Silhouette 计算。

## scikit-learn 映射与参考

| pgr 命令 | 状态 | scikit-learn 对应 | 适用场景 | 备注 |
| :--- | :--- | :--- | :--- | :--- |
| `clust cc` | ✅ 已实现 | `connected_components` (scipy) | 简单连通分量 | 快速去重 |
| `clust mcl` | ✅ 已实现 | (无直接对应) | 图/网络聚类 | 马尔可夫流模拟 |
| `clust k-medoids` | ✅ 已实现 | `KMedoids` (sklearn-extra) | 任意距离矩阵 | 抗噪 |
| `clust dbscan` | ✅ 已实现 | `DBSCAN` | 密度聚类, 噪声识别 | 需调参 eps/min_samples |
| `clust gmm` | 📅 规划中 | `GaussianMixture` | 连续向量, 软聚类 | 需处理协方差矩阵奇异性 |
| `clust kmeans` | 📅 规划中 | `KMeans` | 欧氏空间, 硬聚类 | 基础基准 |

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

1. **基础图聚类**：优先实现 MCL 与 CC（已在 [clust/cc.rs](file:///c:/Users/wangq/Scripts/pgr/src/cmd_pgr/clust/cc.rs) 与 [clust/mcl.rs](file:///c:/Users/wangq/Scripts/pgr/src/cmd_pgr/clust/mcl.rs) 原型中）。
2. **向量支持**：建立读取稠密向量/矩阵的基础设施。
3. **统计聚类**：引入 GMM 实现（可基于 `linfa` 或自研简化版 EM 算法），重点支持 BIC 模型选择。
