# pgr clust cut

`pgr clust cut` 用于将 Newick 树（系统发生树或层次聚类树）切分为扁平的聚类分组（Partition）。

与 `pgr clust`（从数据构建聚类）不同，`cut` 关注于“从已有树结构导出分组”。它支持多种生物学与统计学切割规则，并提供稳定、可复用的表格输出。

本文描述该命令的算法模式、选参思路及输入输出约定。

## 适用场景与设计理念

在实际分析中，我们经常已经有一棵树（系统发生树或层次聚类树），并希望在某个阈值下把叶子切分为不同的分组（Partition）。切割规则可能不止一种：按高度、按簇数、按簇内最大距离（直径）、要求单系群（clade）等。

`pgr clust cut` 旨在提供一套高效、规范且功能全面的切割工具：

- **算法核心**：
  实现了基于簇数 (`--k`) 和高度 (`--height`) 的基础切割，逻辑与 `SciPy.cluster.hierarchy` 及 R `cutree` 保持一致；同时完整移植了 TreeCluster 的生物学约束算法（如 `--max-clade`, `--med-clade` 等），专为系统发生树优化。

- **性能与体验升级**：
  - **高性能**：基于 Rust 实现，无外部依赖，处理大规模树更高效。
  - **标准化**：统一了不同来源算法的术语差异（例如统一使用 `--height`, `--single-linkage`），降低认知负担。
  - **可组合**：作为 `pgr` 工具链的一部分，可直接与 `clust eval` 等命令配合进行聚类评估。

## 支持的模式与算法

`pgr clust cut` 提供了丰富的切割算法，涵盖了从简单的阈值切割到复杂的生物学约束聚类。以下是详细的算法定义与复杂度分析。

### 1. 按簇数切 (`--k <K>`)

- **定义**：将树分割成 $K$ 个簇，使得这 $K$ 个簇是由 $K-1$ 次切割产生的。切割顺序基于节点的高度（距最远叶子的距离），优先切割高度最大的节点。
- **复杂度**：$O(N \log N)$，其中 $N$ 是叶子数量。需要对所有内部节点按高度排序。
- **适用场景**：探索性分析，当你只想要固定数量的分组，而不关心具体的距离阈值时。

### 2. 按高度切 (`--height <H>`)

- **定义**：切断树中所有高度（距最远叶子的距离）大于 $H$ 的边。
  - 对于任意结果簇 $C$，其中的所有节点 $u$ 满足 $height(u) \le H$。
  - 等价于 `SciPy` 的 `fcluster(criterion='distance')` 或 R 的 `cutree(h=H)`。
- **复杂度**：$O(N)$。只需一次后序遍历（Post-order Traversal）。
- **适用场景**：适用于超度量树（Ultrametric Tree），其中高度严格代表时间或遗传距离。

### 3. 按根距离切 (`--root-dist <D>`)

- **定义**：切断树中所有距离根节点路径长度大于 $D$ 的边。
  - 对于任意结果簇 $C$ 的根节点 $r_C$，满足 $dist(root, r_C) \le D$。
  - 一旦某条路径累积长度超过 $D$，该路径即被切断，不再向下延伸。
- **复杂度**：$O(N)$。只需一次前序遍历（Pre-order Traversal）。
- **适用场景**：系统发生树分析，定义从共同祖先（根）演化特定时间后的分化群。

### 4. 按最大簇内直径切 (`--max-clade <T>`)

- **定义**：将叶子划分为互不重叠的簇 $\{C_1, C_2, ..., C_m\}$，使得对于每个簇 $C_i$：
  1. **单系性 (Monophyly)**：$C_i$ 中的叶子必须在原树 $T$ 中构成一个单系群（Clade）。
  2. **直径约束**：$\max_{u, v \in C_i} dist(u, v) \le T$。
  3. **支持度约束**（若指定 `--support`）：$C_i$ 内部的任意两点路径上，不能包含支持度低于阈值的边。
- **算法**：TreeCluster "Max Clade" 算法。采用高效的自底向上（Bottom-Up）直径计算与自顶向下（Top-Down）贪心选择。
- **复杂度**：$O(N)$。避免了 $O(N^2)$ 的全距离矩阵计算。
- **适用场景**：病毒分型、OTU 划分等需要严格控制簇内差异度的场景。

### 5. 按平均簇内距离切 (`--avg-clade <T>`)

- **定义**：类似于 `--max-clade`，但约束条件改为：
  1. **单系性**。
  2. **平均距离约束**：$\frac{1}{|C_i|(|C_i|-1)} \sum_{u, v \in C_i, u \neq v} dist(u, v) \le T$。
  3. **支持度约束**。
- **算法**：TreeCluster "Avg Clade" 算法。自底向上维护子树内的距离和与节点数。
- **复杂度**：$O(N)$。
- **适用场景**：相比最大距离，平均距离对个别离群点（Outlier）更鲁棒。

### 6. 按中位数簇内距离切 (`--med-clade <T>`)

- **定义**：类似于 `--max-clade`，但约束条件改为：
  1. **单系性**。
  2. **中位数距离约束**：$median(\{dist(u, v) \mid u, v \in C_i, u \neq v\}) \le T$。
  3. **支持度约束**。
- **算法**：TreeCluster "Med Clade" 算法。采用自底向上合并排序列表的方式计算中位数。
- **复杂度**：$O(N^2 \log N)$（最坏情况）。相比前两种方法，计算开销显著更大。
- **注意**：不建议用于超大规模树（如 >10k 叶子），除非确实需要中位数鲁棒性。

### 7. 按簇内总枝长切 (`--sum-branch <T>`)

- **定义**：类似于 `--max-clade`，但约束条件改为：
  1. **单系性**。
  2. **总枝长约束**：簇 $C_i$ 对应的最小生成子树的总枝长（Phylogenetic Diversity, PD） $\le T$。
  3. **支持度约束**。
- **生物学意义**：总枝长对应 **Phylogenetic Diversity (PD)**，代表该簇所包含的进化历史总量。
- **算法**：TreeCluster "Sum Branch Clade" 算法。
- **复杂度**：$O(N)$。
- **注意**：
    - PD 是一个**广延量**（随样本数增加而单调增加），不像直径或平均距离那样是“强度量”。
    - 因此，作为切割阈值时，它倾向于把紧密的大簇切碎（因为累积枝长很容易超标），而保留松散的小簇。
    - 除非有特定的生物学理由（如“限制每个 OTU 的最大进化潜能”），否则通常不建议作为首选切割标准。

### 8. 按叶子距离切 (`--leaf-dist-max/min/avg <T>`)

- **定义**：基于簇根节点到叶子的距离进行切割。
  - **Max Leaf Dist** (`--leaf-dist-max <T>`): 切割树，使得簇根到任意叶子的最大距离 $\le T$。
    - 等价于 `root_dist(max_depth - T)`。
    - 类似于 `--height`，但适用于非超度量树（Non-ultrametric Tree），以最远叶子为基准对齐。
  - **Min Leaf Dist** (`--leaf-dist-min <T>`): 切割树，使得簇根到任意叶子的最小距离 $\le T$。
    - 等价于 `root_dist(min_depth - T)`。
  - **Avg Leaf Dist** (`--leaf-dist-avg <T>`): 切割树，使得簇根到所有叶子的平均距离 $\le T$。
    - 等价于 `root_dist(avg_depth - T)`。
- **复杂度**：$O(N)$。需要预先遍历树计算深度统计量。
- **适用场景**：非超度量树（如病毒树），其中“时间”不是统一的，需要以采样时间（叶子）为基准回溯。

### 9. 按最大边长切 (`--max-edge <T>`)

- **定义**：切断树中所有长度大于 $T$ 的边。
  - 对于任意结果簇 $C$，其中的所有边 $e$ 满足 $length(e) \le T$。
  - 这种方法在图论中也被称为 **Single Linkage Clustering**（单链接聚类）：只要两点间存在一条由“短边”（$\le T$）构成的路径，它们就在同一簇。
  - 别名：`--single-linkage <T>`。
- **复杂度**：$O(N)$。只需一次遍历。
- **适用场景**：
  - 去除长枝吸引（Long Branch Attraction）造成的影响。
  - 快速识别紧密相连的群组，忽略稀疏连接。

### 10. 按不一致系数切 (`--inconsistent <T>`)

- **定义**：基于节点与其子树的“不一致性”进行切割。
  - 对于树上每个非叶子节点 $i$，计算其不一致系数 $I_i$。
  - $I_i = \frac{height(i) - \text{mean}(H)}{\text{std}(H)}$，其中 $H$ 是节点 $i$ 下方 $d$ 层（`--deep`）内的所有合并高度集合。
  - **SciPy 参考**: `scipy.cluster.hierarchy.inconsistent`。
  - 若 $I_i > T$，则切断该节点。这意味着该节点的合并高度显著高于其子树的合并高度，暗示这里是一个自然的聚类边界。
- **复杂度**：$O(N)$。
- **适用场景**：当树的整体演化速率不均匀，寻找“自然”聚类边界。

### 11. 动态树切割 (`--dynamic-tree`) [已实现]

- **定义**：参考 R 语言 `dynamicTreeCut` 包的 `cutreeDynamicTree` 算法 ([dynamicTreeCut/R/cutreeDynamic.R](file:///c:/Users/wangq/Scripts/pgr/dynamicTreeCut/R/cutreeDynamic.R))。
- **原理**：自顶向下的递归算法。
  1.  首先基于全局高度进行初步切割。
  2.  对每个初步簇，分析其内部结构（高度分布）。
  3.  如果一个簇内部包含显著的子结构（即存在“高度差”和“子簇大小”满足条件的分割点），则将其进一步递归拆分。
- **参数**：
  - `--min-cluster-size <N>` (默认 20): 最小簇大小。
  - `--deep-split`: 启用更激进的拆分（对应 R 中的 `deepSplit=TRUE`）。
  - `--max-tree-height <H>` (可选): 最大合并高度（默认 99% 树高）。
- **输入**：仅需树结构（Dendrogram），无需距离矩阵。
- **适用场景**：
  - 只需要树结构，追求快速、自动化的切割。
  - 适合处理那些“大簇里套小簇”的嵌套结构。

### 12. 混合动态切割 (`--dynamic-hybrid`) [已实现]

- **定义**：参考 R 语言 `dynamicTreeCut` 包的 `cutreeHybrid` 算法。
- **原理**：自底向上的两阶段算法。
  1.  **Core Detection (核心检测)**:
      - 实现 R `cutreeHybrid` 的自底向上（Bottom-Up）算法，识别满足紧密度（Core Scatter）和分离度（Gap）要求的“基本簇”。
  2.  **PAM-like Reassignment (二次分配)**: 利用原始距离矩阵，将第一阶段未分配的对象（Outliers/Singletons）尝试吸附到最近的核心簇中（Medoid-based assignment）。
      - 默认行为与 R 保持一致 (`pamStage=TRUE`, `pamRespectsDendro=TRUE`, `respectSmallClusters=TRUE`)。
- **输入**：必须提供树结构 + 原始距离矩阵 (`--matrix`)。
- **参数**：
  - `--dynamic-hybrid <N>`: 启用混合切割，N 为最小簇大小。
  - `--matrix <FILE>`: 距离矩阵文件（PHYLIP 格式）。
  - `--max-pam-dist <D>`: PAM 分配的最大距离阈值（默认等于切割高度 `cutHeight`）。如果未分配点到最近 Medoid 的距离超过此值，则保持未分配。
- **适用场景**：
  - 对聚类边界的准确性要求极高。
  - 需要利用距离矩阵信息来修正树结构中可能存在的微小误差或不确定性。
  - 能够有效识别并处理离群点（Outliers）。

### 13. 支持度过滤 (`--support <S>`)

- **定义**：作为上述所有方法的**预处理步骤**。
  - 遍历树中所有边，若某条边的支持度值 $< S$，则将其长度视为 $+\infty$（无限长）。
  - **默认行为**：对于没有明确支持度值的节点（如多叉树解析产生的内部节点），`pgr` 默认其支持度为 100（完全可信）。
- **效果**：任何跨越低支持度边的聚类尝试都会因距离/高度超标而被阻止，从而强制在低支持度处切断。

## 输入与输出

### 输入

- **输入树**：Newick 格式（单棵树）。
- **分支长度**：用于距离/高度相关方法（例如按 root distance、max pairwise distance 等）。
- **分支支持度（可选）**：若树节点/边上携带支持度（例如 bootstrap），可作为“不可跨越”的约束条件。

### 输出

输出建议采用与 `pgr clust dbscan` 相同的格式（便于与既有工具互操作）：

```text
* cluster: Each line contains points of one cluster. The first point is the representative.
* pair: Each line contains a (representative point, cluster member) pair.
```

- `pair` 格式：每行包含 (Representative, Member)。
  - Representative：簇的代表点。
  - Member：簇成员。
  - 单例簇（singleton）：自己是自己的代表点。

**代表点选择 (`--rep`)**：
适用于 `cluster` 和 `pair` 两种格式：
- `root` (默认)：距离根节点最近的成员（字母序作为 Tie-break）。
- `medoid`：Medoid，即到簇内其他成员距离之和最小的成员。
- `first`：字母序第一个成员。

## 工作流与工具链协作

为了保持命令的专注与正交性，我们推荐以下“生成-评估”分离的工作流：

### 1. 生成 (Generation)

使用 `pgr clust cut`：
- 它只负责“切”，不负责“评”。
- 支持多种策略（k, height, max_clade 等）和参数扫描。
- 输出标准 TSV 格式。

### 2. 评估 (Evaluation)

评估聚类质量通常需要参考标准（Ground Truth）或与其他聚类结果对比。这部分逻辑放入独立的 `pgr clust` 或 `pgr nwk` 命令中：

- **通用指标 (`pgr clust eval` / `compare`)**：
  - 输入：两个聚类结果 TSV（或一个结果 + 一个参考）。
  - 输出：ARI (Adjusted Rand Index), AMI (Adjusted Mutual Information), V-Measure 等。
  - 适用场景：当你已知样本的真实分类，或者想比较两种切割参数的差异度时。

- **树相关指标 (`pgr nwk metrics`)**：
  - 输入：树文件 + 聚类结果。
  - 输出：Parsimony score, Silhouette score (基于树上距离矩阵) 等。
  - 适用场景：没有真实分类，需要评估聚类在树结构上的紧密性或分离度。

### 推荐工作流示例

#### 1. 经典系统发育分析
```bash
# 1. 扫描不同参数，生成多个聚类结果
# pgr clust cut input.nwk --method max-clade --scan 0.01,0.05,0.10 > partitions.tsv

# 2. 选定最佳阈值，生成最终聚类
pgr clust cut input.nwk --method max-clade -t 0.05 > final_cluster.tsv

# 3. 可视化或提取子树
pgr nwk subset input.nwk --list final_cluster.tsv --cluster-id 1 > cluster1.nwk
```

#### 2. 层次聚类（hclust）接入
从距离矩阵出发，经由 hclust 生成树，再进行切分与评估。

```bash
# 1. 生成层次聚类树
pgr clust hier matrix.phy --method ward > tree.nwk

# 2. 切分 (按高度阈值切)
pgr clust cut tree.nwk --height 0.05 > clusters.tsv

# 3. 评估 (计算 Cophenetic 相关系数与 Silhouette)
# pgr nwk metrics tree.nwk --part clusters.tsv --metrics silhouette > sil.tsv
```

#### 3. SciPy 风格分析 (不一致系数)
对于演化速率不均匀的树，使用不一致系数可以找到更自然的聚类边界。

```bash
# 使用不一致系数切割 (默认 depth=2)
pgr clust cut tree.nwk --inconsistent 1.5 > clusters.tsv
```

## 选择阈值/簇数：扫描与准则

在 `cut` 场景里，用户常见的选择有两类：

- 直接指定簇数 `K`（类比 R `cutree(k=...)`）。
- 指定阈值 `t`（距离/高度/直径等），由阈值决定切割后的簇数。

当你不确定 `K` 或 `t` 应该取多少时，推荐使用 `--scan` 进行参数扫描。

### 扫描（Scan）

`pgr` 提供显式的扫描能力：适用于所有基于数值参数的方法（如 `--k`, `--height`, `--max-clade`, `--inconsistent` 等）。

**用法**：
`pgr clust cut ... --scan <start>,<end>,<steps>`
（注：扫描仅针对方法的**主阈值参数**。例如对于 `--inconsistent`，扫描的是系数阈值 `T`，而深度 `--deep` 保持固定为用户指定值或默认值）

**输出指标表**：
| Group | Clusters | Singletons | Non-Singletons | MaxSize |
| :--- | :--- | :--- | :--- | :--- |
| height=0.01 | 500 | 480 | 20 | 5 |
| height=0.02 | 300 | 200 | 100 | 15 |
| ... | ... | ... | ... | ... |

- **Non-Singletons**: 即 TreeCluster `argmax_clusters` 试图最大化的指标。
- **MaxSize**: 辅助判断是否存在“超级大簇”（under-clustering）。

### 扫描模式的输出格式

当启用 `--scan` 时，`--format` 参数将被忽略。输出行为如下：

1.  **标准输出 (`stdout` 或 `-o`)**：始终输出详细的分区表（Long format / Tidy Data）。
    - 列定义：`Group`, `ClusterID`, `SampleID`。
    - `Group` 列格式为 `Method=Value`（例如 `height=0.5`, `max-clade=0.02`），便于区分不同的切割参数。
    - 这种格式可以直接作为 `pgr clust eval --format long` 的输入进行批量评估。
2.  **统计输出 (`--stats-out`)**：若指定，将摘要统计表（阈值, 簇数, 单例数, 非单例数, 最大簇大小）写入该文件。

示例：
```bash
# 1. 仅输出详细分区表（用于后续分析或评估）
pgr clust cut tree.nwk --max-clade 0.5 --scan 0,0.5,0.01 > partitions.tsv

# 2. 同时保存统计信息（用于快速检视）
pgr clust cut tree.nwk --max-clade 0.5 --scan 0,0.5,0.01 -o partitions.tsv --stats-out stats.tsv
```

### 与 `pgr clust eval` 的联动

`pgr clust cut` 与 `pgr clust eval` 通过 Long Format 完美配合，支持两种评估模式：

#### 1. 批量内部评估 (Batch Internal Evaluation)
不需要 Ground Truth，使用距离矩阵或坐标评估所有扫描生成的阈值。

```bash
# 生成所有阈值的分区，并直接通过管道传给 eval 进行 Silhouette 评估
pgr clust cut tree.nwk --max-clade 0.5 --scan 0,0.5,0.01 | \
    pgr clust eval - --format long --matrix dist.phy > evaluation.tsv
```

#### 2. 针对性外部评估 (Targeted External Evaluation)
如果你手头有 Ground Truth，通常不需要评估所有阈值（计算量大且无必要）。推荐流程：

1. 先用 `--scan` 快速定位几个有意义的候选阈值区间（例如手肘点附近）。
2. 对少数候选阈值，分别运行一次 `pgr clust cut` 生成分区，再用 `pgr clust eval` 计算 ARI/AMI/V-Measure 等外部一致性指标。

示例：

```bash
# 1) 扫描阈值，先看摘要趋势
pgr clust cut tree.nwk --max-clade 0.5 --scan 0,0.5,0.01 > scan.tsv

# 2) 选定阈值后生成分区
pgr clust cut tree.nwk --max-clade 0.12 > pred.tsv

# 3) 与 ground truth 对比（Partition vs Partition）
pgr clust eval pred.tsv truth.tsv -o eval.tsv
```

### 选点策略参考

当你不确定最佳阈值时，可以使用 `--scan` 生成数据，并参考以下两种常用策略进行决策：

#### 策略 1：最大化非单例簇 (Max Non-Singletons)

- **原理**：寻找一个阈值，使得生成的簇中“非单例簇（Non-Singleton Clusters）”的数量最多。
- **适用性**：当你期望得到尽可能多有意义的（包含 >1 个成员）聚类结果，同时避免过度切碎（导致大量单例）或欠切分（导致巨大簇）时。
- **操作**：观察扫描结果表中的 `Non-Singletons` 列，选择其最大值对应的阈值。

#### 策略 2：手肘规则 (Elbow Rule)
这是数据分析中的通用策略。

- **原理**：观察阈值与簇数量（或单例数）的变化曲线，寻找“拐点”。
  - **陡峭下降期**：随着阈值放松，簇数量迅速减少（大量微小簇合并）。
  - **平缓平台期**：簇数量变化趋于稳定。
  - **拐点（手肘）**：即从“陡峭”转变为“平缓”的点，通常对应着数据内在的自然结构。
- **操作**：
  1. 运行扫描：`pgr clust cut ... --scan ... > scan.tsv`
  2. 观察变化率：若阈值从 $T_1$ 增至 $T_2$ 时簇数剧烈变化，而从 $T_2$ 增至 $T_3$ 时变化平缓，则 $T_2$ 可能是最佳切点。
  3. 可视化：将 `scan.tsv` 导入绘图工具辅助判断。

#### 策略 3：基于评估指标 (Evaluation Metrics)
这是最严谨的策略，通过 `pgr clust eval` 计算聚类质量指标。

- **原理**：直接计算分区的内部有效性（如 Silhouette）或外部一致性（如 ARI，如果有 Ground Truth）。
- **操作**：结合 `pgr clust eval` 使用。
  ```bash
  # 生成所有候选分区的详细列表
  pgr clust cut ... --scan ... > partitions.tsv
  # 批量评估
  pgr clust eval partitions.tsv --format long --matrix dist.phy
  ```

## 现有工具参考 (Prior Art)

`pgr clust cut` 的设计吸收了多个领域的最佳实践：

- **SciPy (`scipy.cluster.hierarchy`)**:
  - 提供了 `fcluster` 函数，支持按高度 (`distance`)、簇数 (`maxclust`) 和不一致系数 (`inconsistent`) 切割。
  - `pgr` 复用了其 `height` 和 `inconsistent` 的定义。

- **R (`dynamicTreeCut`)**:
  - 提供了 `cutreeDynamic` (Tree) 和 `cutreeHybrid` (Hybrid) 方法。
  - 引入了“自适应递归切割”和“核心检测+二次分配”的思路。
  - `pgr` 完整实现了其 `Dynamic Tree` 和 `Hybrid` 算法，提供了高性能的 Rust 版本。

- **TreeCluster**:
  - 专为系统发生树设计，引入了 `Max Clade` (直径)、`Avg Clade` 等约束。
  - 解决了非超度量树的切割问题。
  - `pgr` 完整实现了其核心算法集。

- **R (`cutree`)**:
  - 提供了最基础的 `h` (高度) 和 `k` (簇数) 切割。
