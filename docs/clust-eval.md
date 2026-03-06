# clust eval

`pgr clust eval` 提供通用的聚类质量评估与比较功能。它支持**外部有效性（External Validity）**（与 Ground Truth 对比）和**内部有效性（Internal Validity）**（基于数据本身的几何/统计特性）。

## 设计哲学

`pgr` 采用**组件化**的设计哲学，将聚类的**生成（Generation）**与**评估（Evaluation）**分离。这与 Python 包 `clusteval` 的一体化设计不同：

*   **Python `clusteval`**: `fit()` 方法内部自动执行 Grid Search（尝试不同的 $k$ 或 $\epsilon$），计算内部指标（如 Silhouette），并返回最优结果。
*   **`pgr` Workflow**:
    1.  **生成**: 使用 `pgr nwk cut --scan` 或 `pgr clust dbscan --scan` 生成一系列候选聚类方案（Partitions）。
    2.  **评估**: 使用 `pgr clust eval` 批量计算这些方案的评估指标。
    3.  **决策**: 用户根据指标（如 Silhouette 峰值、Elbow 点）选择最优参数。

- **互补**：
  - `pgr nwk eval` [计划中]：关注树结构与分组的一致性（几何/演化）。
  - `pgr clust eval`：关注分组的统计有效性，支持外部（两分组对比）与内部（单分组+矩阵/坐标/树）评估。
- **场景**：
  - **算法对比**：比较 MCL 与 K-Medoids 在同一数据集上的结果差异。
  - **基准测试**：将聚类结果与已知的标准分类（Ground Truth）对比，计算准确性。
  - **参数调优**：比较不同参数（如 `eps` 或 `inflation`）下聚类结果的稳定性。

这种设计使得评估工具可以独立于聚类算法存在，支持任意来源的聚类结果。

## 核心指标 (Core Metrics)

聚类评估指标通常分为两类：**外部有效性**（依赖 Ground Truth 或参考分区）和**内部有效性**（仅依赖数据本身的几何/统计特性）。

### 1. 外部有效性 (External Validity)
*用于比较两个聚类结果的一致性，或评估聚类结果与真实分类的吻合度。*

#### 1.1 基于配对 (Pair-counting)
*关注样本对在两个分区中是否保持同组/异组关系。*

- **ARI (Adjusted Rand Index)**
  - **定义**：校正了随机性的 Rand Index。
  - **原理**：统计同时在两个分区中属于同组或异组的样本对数量，并减去随机分配下的期望值。
  - **范围**：`[-1, 1]`。1 表示完全一致，0 表示随机水平，负值表示比随机更差。
  - **优点**：
    - **可解释性强**：0 为随机基线，直观。
    - **对称性**：`ARI(A, B) == ARI(B, A)`。
  - **缺点**：对簇的内部结构（如形状）不敏感。
  - **适用**：簇大小不平衡、簇数量较多的通用场景。

- **RI (Rand Index)**
  - **定义**：正确分类的样本对比例。
  - **范围**：`[0, 1]`。
  - **缺点**：未校正随机性。随着簇数量增加，随机分区的 RI 也会趋近于 1，导致区分度降低。通常**不推荐**单独使用。

- **Homogeneity (同质性)**
  - **定义**：每个簇是否只包含某一个类的成员？（类似 Precision，要求簇够“纯”）。
  - **范围**：`[0, 1]`。

- **Completeness (完整性)**
  - **定义**：某一个类的所有成员是否都被分到了同一个簇？（类似 Recall，要求簇够“全”）。
  - **范围**：`[0, 1]`。

- **V-Measure**
  - **定义**：同质性 (Homogeneity) 和完整性 (Completeness) 的调和平均。
  - **范围**：`[0, 1]`。
  - **缺点**：未校正随机性。在样本量小或簇数量多时，得分偏高。
  - **适用**：需要分析聚类误差来源（是分得太碎还是混得太杂）时。

- **FMI (Fowlkes-Mallows Index)**
  - **定义**：Precision 和 Recall 的几何平均数。
  - **原理**：$FMI = \sqrt{\frac{TP}{TP+FP} \times \frac{TP}{TP+FN}}$。
  - **范围**：`[0, 1]`。
  - **适用**：当对 Precision 和 Recall 均有要求时。

#### 1.2 基于信息论 (Information Theoretic)
*关注两个分区所共享的信息量（熵）。*

- **AMI (Adjusted Mutual Information)**
  - **定义**：校正了随机性的互信息 (Mutual Information)。
  - **原理**：基于熵（Entropy）计算两个分区的共享信息，并减去随机期望。
  - **范围**：`[0, 1]`。1 表示完全一致，0 表示随机。
  - **优点**：
    - 对簇数量极多（甚至接近样本数）的情况更鲁棒。
    - 能捕捉非线性的复杂关系。
  - **适用**：小样本、多簇（Large K）场景。

- **NMI (Normalized Mutual Information)**
  - **定义**：标准化的互信息。
  - **原理**：$NMI = \frac{MI(U, V)}{\sqrt{H(U) \cdot H(V)}}$（几何平均）或 $\frac{2 \cdot MI}{H(U) + H(V)}$（算术平均）。`pgr` 采用算术平均。
  - **范围**：`[0, 1]`。
  - **缺点**：未校正随机性。
  - **适用**：簇大小分布均衡的场景。

- **MI (Mutual Information)**
  - **定义**：两个分区之间的互信息量。
  - **原理**：$MI(U, V) = \sum \sum P(u,v) \log \frac{P(u,v)}{P(u)P(v)}$。
  - **范围**：`[0, +∞)`。
  - **缺点**：难以直接解释，受分区熵影响大。通常作为 AMI/NMI 的中间计算步骤。

#### 1.3 基于集合匹配 (Set Matching)
*关注簇与类之间的最佳匹配关系。*

- **Jaccard Index**: 两个集合交集与并集的比率。用于衡量特定簇的重叠度。
- **Precision/Recall**: 基于对（Pair-wise）的精确率与召回率。

---

### 2. 内部有效性 (Internal Validity)
*用于在没有 Ground Truth 的情况下，评估聚类结果本身的质量（紧密度与分离度）。*

#### 2.1 基于距离 (Distance-based)
*需要距离矩阵 (`--matrix`) 或系统发育树 (`--tree`)。*

- **Silhouette Coefficient (轮廓系数)**
  - **原理**：对每个样本 $i$，计算其与同簇样本的平均距离 $a(i)$ 和与最近异簇样本的平均距离 $b(i)$。$s(i) = (b - a) / \max(a, b)$。
  - **范围**：`[-1, 1]`。
    - 接近 1：样本聚类良好（离同簇近，离异簇远）。
    - 0：样本位于簇边界。
    - 负值：样本可能分错簇了。
  - **优点**：直观，兼顾凝聚度和分离度。
  - **缺点**：计算复杂度高 ($O(N^2)$)，大规模数据需优化。

- **Dunn Index**
  - **原理**：最小簇间距离与最大簇内直径之比。
  - **范围**：`[0, +∞)`。**越大越好**。
  - **优点**：简单直观。
  - **缺点**：对噪声极其敏感（因为基于 min/max）。

- **C-Index**
  - **原理**：比较簇内距离之和与整个数据集中最小的 $N_W$ 个距离之和（$N_W$ 为簇内对数）。
  - **范围**：`[0, 1]`。**越小越好**。
  - **缺点**：计算复杂度高 ($O(N^2 \log N)$)，需要对所有成对距离排序。

- **Hubert's Gamma**
  - **原理**：距离矩阵与二值聚类矩阵（0=同簇，1=异簇）之间的相关性。
  - **范围**：`[-1, 1]`。**越大越好**（注意定义中 Y=1 为异簇，通常需根据具体实现确认符号方向，`pgr` 实现中越大表示区分度越好）。

- **Kendall's Tau**
  - **原理**：距离矩阵与聚类矩阵的秩相关系数。
  - **范围**：`[-1, 1]`。**越大越好**。

#### 2.2 基于坐标 (Coordinate-based)
*需要坐标矩阵 (`--coords`)。适用于欧几里得空间数据。*

- **Davies-Bouldin Index (DBI)**
  - **原理**：计算每对簇的“相似度”（簇内散度之和 / 簇心距离），取每个簇最差（最大）相似度的均值。
  - **范围**：`[0, +∞)`。**越小越好**。
  - **优点**：计算比 Silhouette 快。
  - **适用**：评估基于质心的聚类算法。

- **Calinski-Harabasz Index (CH)**
  - **原理**：簇间离散度 (BGSS) 与簇内离散度 (WGSS) 之比。
  - **范围**：`[0, +∞)`。**越大越好**。
  - **优点**：计算快。

- **PBM Index**
  - **原理**：基于总离散度、簇内离散度和簇心最大距离的组合指标。
  - **范围**：`[0, +∞)`。**越大越好**。

- **Ball-Hall Index**
  - **原理**：各簇平均离散度的均值。
  - **范围**：`[0, +∞)`。**越小越好**（越紧凑）。

- **Xie-Beni Index**
  - **原理**：簇内紧凑度与簇间分离度（最小簇心距离）的比值。
  - **范围**：`[0, +∞)`。**越小越好**。

- **Wemmert-Gancarski Index**
  - **原理**：基于相对距离（点到所属簇心距离 / 点到最近其他簇心距离）的紧凑度指标。
  - **范围**：`[0, 1]`。**越大越好**。

### 3. SciPy/Scikit-learn 兼容性与借鉴

`pgr` 的评估指标设计深度借鉴了 `scikit-learn.metrics` 模块，力求在概念和数值上与其保持一致，同时针对生物信息学场景进行优化。

- **混淆矩阵优化**：借鉴 `sklearn.metrics.cluster.pair_confusion_matrix`。
- **距离计算**：支持 PHYLIP 距离矩阵和直接基于 Newick 树的距离计算（Patristic Distance）。

### 4. 指标选择指南

| 场景 | 推荐指标 | 理由 |
| :--- | :--- | :--- |
| **有 Ground Truth** | ARI, AMI | 校正了随机性，结果可信。 |
| **关注聚类纯度** | V-Measure | 可以分别查看 Homogeneity（纯度）和 Completeness（完整性）。 |
| **无 Ground Truth (距离)** | Silhouette | 直观反映几何质量。 |
| **无 Ground Truth (坐标)** | Davies-Bouldin, CH | 计算效率较高。 |
| **簇数量极大** | AMI | 比 ARI 更稳定。 |

### 5. 高级选项 (Advanced Options)

- **`--no-singletons`**
  - **功能**：在评估时，排除 Ground Truth (P2) 中的单例（Singleton，即簇大小为 1 的样本）。
  - **适用场景**：Ground Truth 包含大量未分类或独特的样本（Singletons），而聚类算法主要关注识别群组结构。如果不排除这些单例，可能会因为算法将它们合并或拆分而导致评分不合理地降低。
  - **对齐**：此选项与 `TreeCluster` 的 `score_clusters.py` 中的 `-ns` / `--no_singletons` 功能一致。

## 典型工作流 (Workflows)

### 场景 A: 有 Ground Truth（外部评估）

比较算法生成的聚类结果与已知分类：

```bash
# 比较 result.tsv 和 truth.tsv
pgr clust eval result.tsv truth.tsv
# 输出: ARI, AMI, NMI, FMI, V-Measure...
```

### 场景 B: 无 Ground Truth（内部评估）

#### 1. 使用距离矩阵 (Silhouette, Dunn, etc.)
```bash
# 1. 准备距离矩阵
pgr nwk distance tree.nwk --mode phylip > dist.mat

# 2. 评估
pgr clust eval result.tsv --matrix dist.mat
```

#### 2. 直接使用树文件 (无需生成矩阵)
```bash
# 直接基于树计算距离 (Patristic Distance)
pgr clust eval result.tsv --tree tree.nwk
```

#### 3. 使用坐标/向量 (Davies-Bouldin, CH, etc.)
```bash
# 输入特征向量
pgr clust eval result.tsv --coords vectors.tsv
```

### 场景 C: 批量扫描与评估

结合 `nwk cut --scan` 生成多组阈值结果并批量评估。

```bash
# 1. 扫描树切割阈值，输出长表 (Group, Cluster, ID)
pgr nwk cut tree.nwk --scan 0.01,1.0,0.01 --mode leaf-dist-min > partitions.tsv

# 2. 批量评估内部指标 (直接传入树文件)
pgr clust eval partitions.tsv --format long --tree tree.nwk > scores.tsv

# 3. 查看结果 (找出 Silhouette 最高的阈值)
cat scores.tsv | sort -k2 -nr | head
```

## 输入与输出约定

### 输入
- **单次比较模式**：
  - **Partition 1 (`<p1>`)**: 第一个分组文件（TSV）。
  - **Partition 2 (`<p2>`)**: 第二个分组文件（TSV，可选）。
  - 若提供 `<p2>`，计算外部指标（ARI/AMI）。
  - 若不提供 `<p2>` 且提供了 `--matrix/--tree/--coords`，计算内部指标。
  - 支持 `cluster` / `pair` 格式（通过 `--format` 指定）。

- **批量评估模式 (Batch Mode)**：
  - **Partition (`<p1>`)**: 包含多个分组方案的长表文件（TSV）。
  - 必须指定 `--format long`。
  - **列定义**：
    1. `Group`: 分组标识（如阈值、参数）。
    2. `ClusterID`: 簇 ID。
    3. `SampleID`: 样本 ID。
  - 数据必须按 `Group` 列排序或聚集（程序会按 Group 逐块处理）。
  - 通常与 `pgr nwk cut --scan` 的输出直接对接。
  - 支持 `Group` 列包含 `Method=Value` 格式的元数据（如 `height=0.01`）。

### 输出
- **TSV 格式**，包含所有计算的指标。
- **单次模式 (External)**：
  - 只有一行表头 + 一行数据。
  - 列顺序：`ari`, `ami`, `homogeneity`, `completeness`, `v_measure`, `fmi`, `nmi`, `mi`, `ri`, `jaccard`, `precision`, `recall`。
- **单次模式 (Internal)**：
  - 两行：Metric Name + Value。
  - 若提供 `--matrix`/`--tree`，列顺序：`silhouette`, `dunn`, `c_index`, `gamma`, `tau`。
  - 若提供 `--coords`，列顺序：`davies_bouldin`, `calinski_harabasz`, `pbm`, `ball_hall`, `xie_beni`, `wemmert_gancarski`。
- **批量模式 (Batch)**：
  - 一行表头 + 多行数据（每组一行）。
  - 第一列为 `Group`。
  - 后续列取决于提供的参数（如果同时提供多种参数，列将按以下顺序拼接）：
    1. 有 `<p2>`：包含所有 External 指标。
    2. 有 `--matrix`/`--tree`：包含所有基于距离的 Internal 指标。
    3. 有 `--coords`：包含所有基于坐标的 Internal 指标。

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

#### ClustEval 组件概览与行为
- `silhouette.fit`: 针对聚类数 `k` 做网格搜索，计算每个 `k` 的 Silhouette 得分并取最大值；输出 `labx`（最佳 k 的标签）和 `score` 表（`cluster_threshold`, `clusters`, `score`）。
- `dbindex.fit`: 计算 Davies–Bouldin Index（越小越好），对不同 `k` 取最小值；输出最佳 `labx` 和 `score` 表（`clusters`, `score`）。
- `derivative.fit`: 基于层次聚类的合并高度二阶差分（Elbow），选择加速度最大的 `k`；输出 `labx`（无显式评分）。
- `dbscan.fit`: 对 `eps` 扫描（默认 0.1..5，指定分辨率），用 Silhouette 选最优 `eps`；输出 `labx` 和扫描曲线（`eps`, `silscores`, `sillclust`）。
- `hdbscan.fit`: 调用 `hdbscan`，输出 `labels_`、`probabilities_`、`cluster_persistence_`、`outlier_scores_` 等信息；提供树/凝缩树绘图。
- `coord2density`: 用 `KernelDensity` 计算坐标密度（可视化辅助手段）。
- `plot_dendrogram`/`bubblegrid`: 绘图辅助（树切线、矩阵气泡图）。

#### 与 `pgr` 的差异与可采纳点
- 输入约定差异：ClustEval 偏向“坐标/观测矩阵”作为输入；`pgr` 当前 `clust eval` 偏向“Partition vs Partition”外部有效性评估。两类评估应该并存：
  - 外部有效性（当前主线）：`ARI/AMI/V-Measure`，输入为两个分区。
  - 内部有效性（可作为补充）：`Silhouette/DBIndex`，输入为坐标或距离矩阵 + 单个分区（或直接做 k/eps 网格搜索）。
- 扫描与评分：ClustEval 把“扫描 + 评分 + 选最优”打包在一起；`pgr` 更推荐“扫描（生成表）”与“决策/评估（独立命令）”解耦，便于组合与审计。
- HDBSCAN/DBSCAN：算法本身属于“聚类”范畴，评分是“评估”。`pgr` 侧更适合在 `clust` 命令里实现算法，在 `clust eval`/`nwk metrics` 里补内外部指标。
- 可视化：文档保留参考，CLI 仅输出 TSV

## 实现与采纳要点

- 技术基线（sklearn）：
  - 指标与算法参考路径：[sklearn/metrics/cluster/_supervised.py](file:///c:/Users/wangq/Scripts/pgr/scikit-learn-main/sklearn/metrics/cluster/_supervised.py)，EMI：[expected_mutual_info_fast.pyx](file:///c:/Users/wangq/Scripts/pgr/scikit-learn-main/sklearn/metrics/cluster/_expected_mutual_info_fast.pyx)
  - 列联表：COO/CSR 构建二维直方图；Rust 侧用 `HashMap<(u32,u32),u32>` 或 CSR
  - ARI 高效计算：用列联表求同簇对（∑ C(n_ij,2)），避免 $O(N^2)$
  - 数值稳定性：处理 `x*log(x)` 的零值；使用 `f64`
  - EMI 实施：移植 `_expected_mutual_info_fast` 的核心逻辑
  - **Silhouette**：采用分块计算（Chunking）策略，避免生成全量 $N \times N$ 距离矩阵，降低内存消耗；支持预计算距离矩阵（Precomputed）。
  - **Davies-Bouldin**：仅在提供坐标矩阵时启用，或者实现基于 Medoid 的变体以支持距离矩阵。
- 性能策略：
  - 对齐输入样本（交集 + 排序）；列联表构建 $O(N)$；指标计算稀疏下按非零计
  - 内部有效性指标（Silhouette）利用 Rayon 并行计算各样本的 $a(i)$ 和 $b(i)$。
- 采纳方案（clusteval 融合）：
  - 外部有效性：默认输出 `ARI/AMI/V-Measure`，与 `nwk cut --scan` 联动选阈值
  - 内部有效性：
    - 不再试图构建一个万能的 `clust eval` 命令。
    - **基础评估**：`pgr clust eval` 增加 `--matrix` 或 `--coords` 参数，计算单分区指标。
    - **扫描集成**：在 `pgr nwk cut --scan` 和 `pgr clust dbscan --scan` 中直接集成评估逻辑，输出含指标的扫描表。
  - 算法整合：`DBSCAN/HDBSCAN` 保留在 `pgr clust`，评估统一在 `clust eval`/`nwk metrics`
  - 可视化：文档保留参考，CLI 仅输出 TSV
- 代码参考：
  - Silhouette：[silhouette.py](file:///c:/Users/wangq/Scripts/pgr/clusteval-2.2.7/clusteval/silhouette.py)
  - DBIndex：[dbindex.py](file:///c:/Users/wangq/Scripts/pgr/clusteval-2.2.7/clusteval/dbindex.py)
  - DBSCAN 扫描与评分：[dbscan.py](file:///c:/Users/wangq/Scripts/pgr/clusteval-2.2.7/clusteval/dbscan.py)
  - HDBSCAN：[hdbscan.py](file:///c:/Users/wangq/Scripts/pgr/clusteval-2.2.7/clusteval/hdbscan.py)

## 测试策略 (Testing Strategy)

参考 `scikit-learn` 和 `clusteval` 的测试方法，我们将采用以下策略确保 `pgr clust eval` 的正确性与鲁棒性。

### 1. 对照测试 (Exactness against Reference)
*目标：确保核心算法实现与业界标准（scikit-learn）完全一致。*

- **测试数据生成**：
  - 使用 Python 脚本生成多组典型聚类结果（包括高一致性、随机、完全不一致）。
  - 调用 `sklearn.metrics` 计算预期指标（ARI, AMI, V-Measure, Silhouette, DBIndex）。
  - 将输入分区与预期分数保存为测试用例（JSON/TSV）。
- **Rust 集成测试**：
  - 读取测试用例，运行 `pgr clust eval`。
  - 断言计算结果与 `sklearn` 的误差在 `1e-10` 范围内。
  - *注意*：AMI 的计算依赖于 `log` 底数（通常为 `e` 或 `2`）和列联表构建方式，需确保参数对齐（sklearn 默认 `log_e`）。

### 2. 不变性测试 (Invariance)
*目标：确保指标仅依赖于分区的数学结构，而非表达形式。*

- **标签置换 (Label Permutation)**：
  - 将分区中的 Cluster ID 随机重命名（如 `1->A, 2->B` 变为 `1->B, 2->A`）。
  - 断言 ARI/AMI/V-Measure 等指标结果**完全不变**。
- **样本顺序 (Sample Order)**：
  - 打乱输入文件的行顺序（保持 ID 对应关系）。
  - 结果应完全不变。
- **标签缩放 (Label Scaling)**：
  - 使用非连续整数或大整数作为 Cluster ID（如 `1, 100, 10000`）。
  - 结果应完全不变。

### 3. 边界条件 (Boundary Conditions)
*目标：处理极端情况，避免 Panic 或 NaN。*

- **完全一致 (Perfect Match)**：
  - 输入两个完全相同的分区。
  - 预期：ARI=1.0, AMI=1.0, V-Measure=1.0。
- **单簇 (Single Cluster)**：
  - 所有样本都属于同一个簇。
  - 预期：ARI=0.0, AMI=0.0。
- **全单例 (All Singletons)**：
  - 每个样本自成一簇（簇数 = 样本数）。
  - 预期：ARI=0.0, AMI 视归一化方法而定（通常接近 0 或 1，需查阅 sklearn 定义）。
- **空输入 (Empty Input)**：
  - 0 个样本。
  - 预期：返回错误提示或特定的空值，不应 Panic。

### 4. 随机基线 (Random Baseline)
*目标：验证 Adjusted 指标的归一化特性。*

- **随机分区**：
  - 生成两个完全独立的随机分区（样本量 N > 1000）。
  - 预期：ARI 和 AMI 应接近 0.0（允许微小正负波动）。
  - *注意*：非 Adjusted 指标（如 RI, V-Measure）在随机情况下通常 > 0。

### 5. 内部有效性特有测试
- **Silhouette**:
  - 验证单样本簇（Cluster size = 1）的处理（通常定义为 0）。
  - 验证距离矩阵对角线为 0。
- **DBSCAN 噪声**:
  - 验证噪声点（通常标记为 -1 或 Unclassified）在评估时的处理方式（作为独立簇还是忽略）。`sklearn` 通常将噪声视为独立簇或忽略，需明确 `pgr` 策略（建议：视为独立单例或统一为一个特殊簇，需在文档中明确）。

### 6. 验证清单 (Verification Checklist)
- Perfect Matches：ID 变更不影响结果（ARI=1.0）
- Non-consecutive Labels：非连续 ID 不影响结果
- Homogeneity/Completeness：单侧完美情况验证（H=1 vs C=1）
- Random Baseline：随机分区的 Adjusted 指标接近 0
- 与 sklearn 对齐：在小规模数据上交叉验证 ARI/AMI/V-Measure 的一致性

## 实施计划

### 阶段 1：外部有效性 MVP [已完成]
- [x] CLI：`pgr clust eval <p1> <p2> -o eval.tsv`
- [x] 算法：ARI, AMI, V-Measure, FMI, NMI, MI, RI, Jaccard, Precision, Recall。

### 阶段 2：内部有效性（指标库） [已完成]
- **基于距离 (`--matrix`, `--tree`)**：
  - [x] Silhouette Coefficient
  - [x] Dunn Index
  - [x] C-Index
  - [x] Hubert's Gamma
  - [x] Kendall's Tau
- **基于坐标 (`--coords`)**：
  - [x] Davies-Bouldin Index
  - [x] Calinski-Harabasz Index
  - [x] PBM Index
  - [x] Ball-Hall Index
  - [x] Xie-Beni Index
  - [x] Wemmert-Gancarski Index

### 阶段 3：扫描与集成 [已完成]
- [x] 批量评估模式：`--format long` 支持读取多组分区。
- [x] 输入兼容：与 `pgr nwk cut --scan` 输出直接对接。

### 阶段 4：树结构整合 [已完成]
- [x] 直接树支持：`--tree <FILE>` 参数，无需预先生成距离矩阵。
- [x] 统一接口：`DistanceMatrix` trait 适配 PHYLIP 矩阵和 Newick 树。
