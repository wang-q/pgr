# clust eval [设计中]

`pgr clust eval` 提供通用的聚类质量评估与比较功能。它主要关注**外部有效性（External Validity）**，即通过与 Ground Truth 或其他聚类结果的对比来量化一致性。

## 定位与场景

- **定位**：通用聚类评估工具，不依赖树结构。
- **互补**：
  - `pgr nwk eval`：关注树结构与分组的一致性（几何/演化）。
  - `pgr clust eval`：关注两个分组（Partition）之间的一致性（统计/信息论）。
- **场景**：
  - **算法对比**：比较 MCL 与 K-Medoids 在同一数据集上的结果差异。
  - **基准测试**：将聚类结果与已知的标准分类（Ground Truth）对比，计算准确性。
  - **参数调优**：比较不同参数（如 `eps` 或 `inflation`）下聚类结果的稳定性。

## 核心指标

### 1. 基于配对 (Pair-counting)
*关注样本对在两个分区中是否保持同组/异组关系。*

- **ARI (Adjusted Rand Index)**:
  - 最常用的聚类一致性指标。
  - 范围 `[-1, 1]`，0 表示随机，1 表示完全一致。
  - 对簇大小不平衡鲁棒。
- **RI (Rand Index)**: 基础的一致性比例（未校正随机性）。

### 2. 基于信息论 (Information Theoretic)
*关注两个分区所共享的信息量。*

- **AMI (Adjusted Mutual Information)**:
  - 校正了随机性的互信息。
  - 范围 `[0, 1]`。
  - 相比 ARI，更适合簇数量较多或簇大小极度不平衡的场景。
- **V-Measure**:
  - **Homogeneity (同质性)**: 每个簇是否只包含某一个类的成员？（类似 Precision）
  - **Completeness (完整性)**: 某一个类的所有成员是否都被分到了同一个簇？（类似 Recall）
  - V-Measure 是两者的调和平均。

### 3. 基于集合匹配 (Set Matching)
*关注簇与类之间的最佳匹配。*

- **Jaccard Index**: 集合重叠度。
- **F1 Score**: 基于 Precision 和 Recall 的综合指标。

## 输入与输出约定

### 输入
- **Partition 1 (`-1` / `--p1`)**: 第一个分组文件（TSV）。
- **Partition 2 (`-2` / `--p2`)**: 第二个分组文件（TSV）。
  - 格式：`Item <tab> ClusterID`。
  - 支持 `pgr` 标准的 `cluster` 格式（每行一个簇）或 `pair` 格式（每行一个成员）。
- **Merge Strategy**: 自动取两个文件样本的交集进行评估。

### 输出
- **TSV 格式**，包含所有计算的指标。

## 典型用法

```bash
# 比较聚类结果与 Ground Truth
pgr clust eval --p1 clustering_result.tsv --p2 ground_truth.tsv > eval.tsv
# 输出: ARI, AMI, Homogeneity, Completeness, V-Measure

# 比较两个算法的结果
pgr clust eval --p1 mcl_clusters.tsv --p2 dbscan_clusters.tsv
```

## 实现备注（技术细节）

- **列联表 (Contingency Table)**:
  - 所有指标的基础是构建 $R \times C$ 的列联表（混淆矩阵），其中 $n_{ij}$ 表示同时属于 Partition 1 的第 $i$ 簇和 Partition 2 的第 $j$ 簇的样本数。
  - 稀疏优化：对于大数量簇，列联表应使用稀疏矩阵或 Hash Map 存储。
- **性能策略**:
  - 核心算法复杂度通常为 $O(N)$（构建列联表）。
  - AMI 计算涉及大量对数运算，需注意数值稳定性。

## 实施计划

- [ ] **CLI 搭建**: 支持读取两个 Partition 文件并对齐样本。
- [ ] **核心算法**:
    - 实现列联表构建。
    - 实现 ARI, AMI, Homogeneity, Completeness, V-Measure。
- [ ] **验证**:
    - 对比 `scikit-learn.metrics` 的结果。
