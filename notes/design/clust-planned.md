# clust 计划中算法

> **实现状态注记**：本文档列出尚未实现的聚类算法规划。当前 `pgr clust` 已实现 hier/dbscan/mcl/k-medoids/cc/nj/upgma + cut/eval。

## 1. GMM (Gaussian Mixture Models)

- **原理**：假设数据由 $K$ 个高斯分布混合而成。使用 EM（期望最大化）算法迭代估计每个高斯分量的参数（均值、协方差）及每个样本属于各分量的后验概率。
- **命令**：`pgr clust gmm`
- **计划内容**：支持**软聚类**（概率输出）与 **BIC** 模型选择。
- **价值**：适合**椭球形簇**与密度估计，解决 K-Means 仅适应球形簇的限制；BIC 可辅助确定最佳 K。

## 2. HDBSCAN

- **原理**：结合层次聚类与 DBSCAN。通过构建基于密度的层次树（Condensed Tree），并根据簇的稳定性（Stability）在不同层级自动提取最佳簇，无需全局 $\epsilon$。
- **命令**：`pgr clust hdbscan`
- **scikit-learn 对应**：`HDBSCAN`
- **计划内容**：层次化 DBSCAN，无需手动指定 `eps`。
- **价值**：DBSCAN 的现代高级版，**自动适应不同密度的簇**，参数更少且更稳健。

## 3. Louvain / Leiden

- **原理**：基于模块度（Modularity）优化的社区发现算法。Louvain 贪心地最大化模块度；Leiden 改进了 Louvain 的局部合并策略，保证连通性并加速收敛。
- **命令**：(待定)
- **计划内容**：社区发现算法。
- **价值**：比 MCL 更适合**超大规模网络**的层次化结构探索。
