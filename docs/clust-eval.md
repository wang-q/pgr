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

## 现有工具参考 (Prior Art)

### 1. Scikit-learn (`sklearn.metrics`)
- **地位**：Python 数据科学事实标准。
- **支持**：ARI, AMI, V-Measure, Fowlkes-Mallows, Silhouette。
- **特点**：API 设计简洁（`metric(labels_true, labels_pred)`），所有 Adjusted 指标均为默认配置。`pgr` 的核心算法逻辑将与 sklearn 对齐以确保结果可比性。

### 2. R (`mclust`, `fpc`)
- **地位**：统计学与生物信息学常用。
- **特点**：`mclust` 专注于 ARI；`fpc` 混合了内部有效性（距离统计）和外部有效性。

### 3. ClustEval (Bioinformatics)
- **特点**：专为大规模生物聚类设计，强调处理数百万序列的能力。
- **启示**：对于生物数据（如 Gene Families），簇的数量可能极大（>10k）。`pgr` 在实现列联表时需采用稀疏策略（HashMap），避免 $O(K^2)$ 的内存消耗。

## 实现备注（技术细节）

- **Scikit-learn 借鉴**:
  - **参考路径**:
    - 核心指标逻辑：[sklearn/metrics/cluster/_supervised.py](file:///c:/Users/wangq/Scripts/pgr/scikit-learn-main/sklearn/metrics/cluster/_supervised.py)
    - EMI 算法：[sklearn/metrics/cluster/_expected_mutual_info_fast.pyx](file:///c:/Users/wangq/Scripts/pgr/scikit-learn-main/sklearn/metrics/cluster/_expected_mutual_info_fast.pyx)
  - **列联表 (Contingency Table)**: 
    - 核心实现：利用稀疏矩阵（COO/CSR）构建二维直方图。
    - 逻辑：统计 `(true_label, pred_label)` 对的频次。
    - Rust 方案：`HashMap<(u32, u32), u32>`（最通用）或 CSR（Compressed Sparse Row，若 ID 已映射为紧凑整数，可大幅节省内存）。
  - **ARI 高效计算**:
    - 避免 $O(N^2)$ 遍历所有样本对。
    - 利用列联表：$ \sum_{ij} \binom{n_{ij}}{2} $ 计算同簇对数量。
  - **数值稳定性**: 计算 Entropy 和 MI 时，需处理 `x * log(x)` 当 `x=0` 的情况（应为 0），并使用 `f64` 避免精度溢出。
  - **EMI (Expected Mutual Information)**:
    - AMI 的核心难点。涉及超几何分布的期望值。需仔细移植 `_expected_mutual_info_fast` 的逻辑。

- **性能策略**:
  - **输入对齐**: 两个输入文件可能包含不完全重叠的样本。第一步必须是**取交集**并**按样本名排序**，生成对齐的 Label 数组。
  - **算法复杂度**: 构建列联表为 $O(N)$。基于列联表的指标计算通常为 $O(K_1 \times K_2)$（稀疏情况下为 $O(\text{NonZero})$）。

## 实施计划

- [ ] **CLI 搭建**: 支持读取两个 Partition 文件并对齐样本。
- [ ] **核心算法**:
    - 实现列联表构建。
    - 实现 ARI, AMI, Homogeneity, Completeness, V-Measure。
- [ ] **验证**:
    - [ ] **Perfect Matches**: ID 重命名不影响结果（ARI=1.0）。
    - [ ] **Non-consecutive Labels**: 非连续 ID（如 `0, 4`）不影响结果。
    - [ ] **Homogeneity/Completeness**: 验证单侧完美情况（H=1 vs C=1）。
    - [ ] **Integer Overflow**: 确保在大样本（N > 65536）下计数器不溢出（使用 `u64/usize`）。
    - [ ] **Random Baseline**: 随机分区的 Adjusted 指标应接近 0。
