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

#### 与 ROC/PR 的关系（为什么这里不画“曲线”）

不少用户会把 “同簇/不同簇” 的配对指标联想到 ROC/PR 曲线：把每个样本对 `(i, j)` 当成二分类（同簇=正类），然后随阈值变化画出曲线。

这里需要区分两种输入类型：

- **Partition（分区结果）**：输入是一个固定的聚类划分（`cluster` 或 `pair`）。它对每个样本对只给出 0/1（同簇或不同簇），没有可调阈值，因此只能对应 ROC/PR 空间中的一个点，而不是一条曲线。
- **Scored pairs / 可扫阈值的过程**：如果你对每个样本对还有一个连续分数（相似度/距离/共聚类概率），或者聚类过程本身可随阈值连续切分（例如层次聚类的 cut height），才会自然产生 ROC/PR 曲线与 AUC。

`pgr clust eval` 的定位是 “Partition vs Partition” 的一致性评估，因此核心输出是 ARI/AMI/V-Measure 这类对两个分区的整体比较；若你确实需要 ROC/PR 的曲线视角，通常意味着要从“带阈值”的来源（距离阈值、树高阈值、或相似度阈值）生成一系列分区，再逐点计算对应的 TP/FP 等统计量。

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
- **Partition 1 (`<p1>`)**: 第一个分组文件（TSV，位置参数 1）。
- **Partition 2 (`<p2>`)**: 第二个分组文件（TSV，位置参数 2）。
  - 支持与其他命令一致的两种分组表示（`cluster` / `pair`），需通过 `--format` 参数指定（默认为 `pair`）：
    ```text
    * cluster: Each line contains points of one cluster. items separated by tabs.
    * pair: Each line is "Representative/Label <tab> Member".
    ```
  - 评估内部会将输入统一规整为 “每个样本一个标签” 的映射：`Item -> Label`。
    - 对于 `pair` 输入，第一列（Representative/Label）作为标签。
    - 对于 `cluster` 输入，行号（自增整数）作为标签。
- 样本对齐：取两个文件样本的交集进行评估。

### 输出
- **TSV 格式**，包含所有计算的指标。
  - 列名约定：`ari`, `ami`, `homogeneity`, `completeness`, `v_measure`
  - 若后续增加指标（如 `ri`, `jaccard`, `f1`），将以新增列的方式扩充，保持现有列不变

## 典型用法

```bash
# 比较聚类结果与 Ground Truth
pgr clust eval clustering_result.tsv ground_truth.tsv -o eval.tsv
# 输出: ARI, AMI, Homogeneity, Completeness, V-Measure

# 比较两个算法的结果
pgr clust eval mcl_clusters.tsv dbscan_clusters.tsv
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
- CLI：`pgr clust eval <p1> <p2> -o eval.tsv`（位置参数 + 统一 `-o`）
- 算法：构建列联表（交集对齐），实现 `ARI/AMI/V-Measure/Homogeneity/Completeness`
- 输入兼容：`cluster/pair` 两种格式，需通过 `--format` 显式指定（默认 `pair`）
- 输出：TSV 列包含上述指标，列名与顺序固定

### 阶段 2：内部有效性（指标库） [进行中]
- **核心算法 (`libs/clust/eval.rs`)**：
  - [x] 实现 `silhouette_score(partition, distance_matrix)`：支持 NamedMatrix。
  - [ ] 实现 `davies_bouldin_score(partition, coordinates)`：支持坐标输入。
- **CLI 增强 (`pgr clust eval`)**：
  - [x] 新增参数：
    - `--matrix <file>`: 输入距离矩阵（PHYLIP）。
  - [ ] 新增参数：
    - `--coords <file>`: 输入坐标矩阵（TSV，用于 DBIndex）。
    - `--methods <list>`: 指定计算指标（默认 `ari,ami`，可选 `silhouette,dbindex`）。
  - 逻辑：
    - 若提供 `<p1>` 和 `<p2>`：计算外部指标（ARI/AMI）。
    - 若提供 `<p1>` 和 `--matrix`：计算内部指标（Silhouette）。

### 阶段 3：扫描与集成（各命令独立支持）
- **层次聚类 (`pgr nwk cut`)**：
  - 增强 `--scan` 模式，支持 `--eval-matrix <file>` 和 `--eval-methods silhouette`。
  - 在遍历阈值切分树时，直接计算指标并追加到输出 TSV 的列中。
- **扁平聚类 (`pgr clust dbscan`)**：
  - 实现 `--scan <eps_range>`。
  - 集成内部指标计算，输出“参数-指标”扫描表。

### 阶段 4：数值与性能
- **内存优化**：对于大规模矩阵，避免全量加载，支持流式读取或分块计算。
