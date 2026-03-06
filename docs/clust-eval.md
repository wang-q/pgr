# clust eval

`pgr clust eval` 提供通用的聚类质量评估与比较功能。它支持**外部有效性（External Validity）**（与 Ground Truth 对比）和**内部有效性（Internal Validity）**（基于数据本身的几何/统计特性）。

## 设计哲学

`pgr` 采用**组件化**的设计哲学，将聚类的**生成（Generation）**与**评估（Evaluation）**分离。这与 Python 包 `clusteval` 的一体化设计不同：

*   **Python `clusteval`**: `fit()` 方法内部自动执行 Grid Search（尝试不同的 $k$ 或 $\epsilon$），计算内部指标（如 Silhouette），并返回最优结果。
*   **`pgr` Workflow**:
    1.  **生成**: 使用 `pgr nwk cut --scan` 或 `pgr clust dbscan --scan` 生成一系列候选聚类方案（Partitions）。
    2.  **评估**: 使用 `pgr clust eval` 批量计算这些方案的评估指标。
    3.  **决策**: 用户根据指标（如 Silhouette 峰值、Elbow 点）选择最优参数。

* **互补**：
  * `pgr nwk eval` [计划中]：关注树结构与分组的一致性（几何/演化）。
  * `pgr clust eval`：关注分组的统计有效性，支持外部（两分组对比）与内部（单分组+矩阵/坐标/树）评估。
* **场景**：
  * **算法对比**：比较 MCL 与 K-Medoids 在同一数据集上的结果差异。
  * **基准测试**：将聚类结果与已知的标准分类（Ground Truth）对比，计算准确性。
  * **参数调优**：比较不同参数（如 `eps` 或 `inflation`）下聚类结果的稳定性。

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

- **Jaccard Index (Pair-wise)**
  - **定义**：在所有被任一分区视为同组的样本对（TP + FP + FN）中，两个分区均视为同组（TP）的比例。
  - **公式**：$J = \frac{TP}{TP + FP + FN}$。
  - **意义**：即集合 $S_1$（P1 中同组对）与 $S_2$（P2 中同组对）的 Jaccard 相似系数。

- **Precision (Pair-wise)**
  - **定义**：在所有被预测分区（P1）视为同组的样本对中，真实分区（P2）也视为同组的比例。
  - **公式**：$P = \frac{TP}{TP + FP}$。

- **Recall (Pair-wise)**
  - **定义**：在所有真实分区（P2）视为同组的样本对中，预测分区（P1）也视为同组的比例。
  - **公式**：$R = \frac{TP}{TP + FN}$。

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

- **Homogeneity (同质性)**
  - **定义**：每个簇是否只包含某一个类的成员？（类似 Precision，要求簇够“纯”）。
  - **原理**：基于条件熵 $H(C|K)$ 计算。若 $H(C|K)=0$，则同质性为 1。
  - **范围**：`[0, 1]`。

- **Completeness (完整性)**
  - **定义**：某一个类的所有成员是否都被分到了同一个簇？（类似 Recall，要求簇够“全”）。
  - **原理**：基于条件熵 $H(K|C)$ 计算。若 $H(K|C)=0$，则完整性为 1。
  - **范围**：`[0, 1]`。

- **V-Measure**
  - **定义**：同质性 (Homogeneity) 和完整性 (Completeness) 的调和平均。
  - **范围**：`[0, 1]`。
  - **缺点**：未校正随机性。在样本量小或簇数量多时，得分偏高。
  - **适用**：需要分析聚类误差来源（是分得太碎还是混得太杂）时。

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

### 4. 指标选择指南

| 场景 | 推荐指标 | 理由 |
| :--- | :--- | :--- |
| **有 Ground Truth (通用)** | ARI, AMI | 校正了随机性，结果可信。 |
| **有 Ground Truth (纯度)** | V-Measure | 可以分别查看 Homogeneity（纯度）和 Completeness（完整性）。 |
| **有 Ground Truth (精确匹配)** | Jaccard, F1/FMI | 关注具体簇或配对的重叠程度（而非整体分布）。 |
| **无 Ground Truth (距离)** | Silhouette | 直观反映几何质量，兼顾凝聚度与分离度。 |
| **无 Ground Truth (距离相关性)** | Gamma, Tau | 评估聚类结构与原始距离矩阵的相关程度。 |
| **无 Ground Truth (坐标)** | Davies-Bouldin, CH | 计算效率较高，适合大规模数据。 |
| **无 Ground Truth (坐标-紧凑性)** | PBM, Xie-Beni | 对簇的紧凑性有更严格的惩罚。 |
| **簇数量极大** | AMI | 比 ARI 更稳定。 |

## 典型工作流 (Workflows)

### 场景 A: 有 Ground Truth（外部评估）

比较算法生成的聚类结果与已知分类：

```bash
# 比较 result.tsv 和 truth.tsv
pgr clust eval result.tsv --truth truth.tsv
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
  - **Other (`--other`)**: 第二个分组文件（TSV，可选，用于外部评估；`--truth` 为同义别名）。
  - 若提供 `--other`，计算外部指标（ARI/AMI）。
  - 可选：`--no-singletons` 在评估前从 `--other` 中排除簇大小为 1 的样本（仅对外部评估生效）。
  - 若不提供 `--other` 且提供了 `--matrix/--tree/--coords`，计算内部指标。
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
    1. 有 `--truth`：包含所有 External 指标。
    2. 有 `--matrix`/`--tree`：包含所有基于距离的 Internal 指标。
    3. 有 `--coords`：包含所有基于坐标的 Internal 指标。
