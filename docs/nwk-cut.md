# pgr nwk cut

`pgr nwk cut` 的目标是：给定一棵 Newick 树（系统发生树 / 层次聚类树），按照用户指定的规则把叶子节点切分成一组互不重叠的分组（partition），并以稳定、可复用的表格格式输出。

它关注的是“从树导出扁平聚类结果”，而不是“从数据构建树”。因此命名为 `cut` 比 `cluster` 更准确：树本身已经表达了层次结构，我们要做的是在树上选取一条切割规则并导出分组。

本文偏设计稿：描述 `pgr nwk cut` 的背景、输入输出约定、算法模式与选参思路，并对比相关生态工具。

## 适用场景

在实际分析中，经常会遇到这样的需求：

- 已经有一棵树（例如系统发生树、基于距离矩阵的层次聚类树、或某种推断得到的 dendrogram）。
- 希望在某个阈值下把叶子分组（比如得到“簇”用于下游统计、注释、画图、或与其他方法比较）。
- 切割规则可能不止一种：按高度切、按簇数切、按簇内最大两两距离（直径）切、要求每个簇必须是单系群（clade）、或者禁止跨越低支持度的边等。

`pgr nwk cut` 旨在提供一套与现有生态对齐但更“树友好”的切割方式：

- **对齐 R `cutree()`**：在 dendrogram 上切一刀得到分组。
- **对齐 SciPy `fcluster()`**：支持按距离 `distance` 或簇数 `maxclust` 导出扁平聚类。
- **对齐 TreeCluster**：在系统发生树上按生物学常用约束得到分组。
- **与 `pgr clust` 区分**：`pgr clust` 主要是从相似度/距离矩阵或图结构“构建聚类”；而 `pgr nwk cut` 是从“已有树”导出分组。

## 功能对照表

为了方便从其他工具迁移，以下是 `pgr` 与主流工具的功能对照：

### vs SciPy (`cluster.hierarchy`)

| SciPy Criterion (`fcluster`) | `pgr nwk cut` 参数 | 说明 | 状态 |
| :--- | :--- | :--- | :--- |
| `maxclust` | `--k <N>` | 指定生成的簇数量 | ✅ 已实现 |
| `distance` | `--height <N>` | 指定切割的高度/距离阈值 | ✅ 已实现 |
| `inconsistent` | `--inconsistent <T>` | 基于不一致系数切割（需配合深度参数） | ✅ 已实现 |
| `monocrit` | - | 基于自定义单调统计量切割 | ❌ 未计划 |

> 注：SciPy 的 `cut_tree` 函数主要对应 `maxclust` (n_clusters) 和 `distance` (height)。

### vs TreeCluster

| TreeCluster Method | `pgr nwk cut` 参数 | 说明 | 状态 |
| :--- | :--- | :--- | :--- |
| `max_clade` | `--max-clade <N>` | 簇内最大两两距离（直径） | ✅ 已实现 |
| `root_dist` | `--root-dist <N>` | 根到叶子的最大距离 | ✅ 已实现 |
| `single_linkage` | `--height <N>` | 等同于按高度切割 | ✅ 已实现 |
| `avg_clade` | `--avg-clade <N>` | 簇内平均两两距离 | ✅ 已实现 |

## 输入与输出

### 输入

- **输入树**：Newick 格式（单棵树）。
- **分支长度**：用于距离/高度相关方法（例如按 root distance、max pairwise distance 等）。
- **分支支持度（可选/规划中）**：若树节点/边上携带支持度（例如 bootstrap），可作为“不可跨越”的约束条件。

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

## 核心概念：切树并导出 partition

不管采用何种规则，`cut` 的结果都可以理解为：

1. 在树上选择一组“切断点”（cut edges / cut nodes）。
2. 切断后，树被分成若干个互不相交的连通分量（component）。
3. 每个连通分量包含若干叶子；这些叶子构成一个输出簇。

不同算法的差异主要在于“切断点如何确定”。

## 选择阈值/簇数：扫描与准则

在 `cut` 场景里，用户常见的选择有两类：

- 直接指定簇数 `K`（类比 R `cutree(k=...)`）。
- 指定阈值 `t`（距离/高度/直径等），由阈值决定切割后的簇数。

当你不确定 `K` 或 `t` 应该取多少时，更稳妥的策略通常是“扫描 + 选点”，而不是一次性拍脑袋给出某个值。

### 扫描（Scan）

`pgr` 计划提供显式的扫描能力：适用于所有基于数值参数的方法（如 `--k`, `--height`, `--max-clade`, `--inconsistent` 等）。它将遍历一组候选值，对每个值计算并输出摘要指标。

**用法规划**：
`pgr nwk cut ... --scan <start>,<end>,<steps>`
（注：扫描仅针对方法的**主阈值参数**。例如对于 `--inconsistent`，扫描的是系数阈值 `T`，而深度 `--deep` 保持固定为用户指定值或默认值）

**与 TreeCluster 的对应关系**：
如果要复刻 TreeCluster 的默认扫描行为，可设置：
- `start = 0`
- `end = <T>` (用户指定的阈值)
- `step = <T> / 1000` (或自定义步长，1000 是 TreeCluster 的硬编码默认值)

**输出指标表（示例）**：
| Threshold | Clusters | Singletons | Non-Singletons | MaxSize |
| :--- | :--- | :--- | :--- | :--- |
| 0.01 | 500 | 480 | 20 | 5 |
| 0.02 | 300 | 200 | 100 | 15 |
| ... | ... | ... | ... | ... |

- **Non-Singletons**: 即 TreeCluster `argmax_clusters` 试图最大化的指标。
- **MaxSize**: 辅助判断是否存在“超级大簇”（under-clustering）。

### 自动选点与 TreeCluster 策略

TreeCluster 提供了一种“无阈值”（Threshold-Free）模式（通过 `-tf argmax_clusters` 启用），其本质是基于扫描的自动选点。

**TreeCluster `argmax_clusters` 的实现细节**：
1.  **输入**：用户指定一个上限阈值 $T_{max}$（即 `-t` 参数）。
2.  **扫描**：在 $[0, T_{max}]$ 区间内均匀取 `NUM_THRESH` (默认为 1000) 个点作为候选阈值。
3.  **评估**：对每个候选阈值执行聚类，计算 **非单例簇（Non-Singleton Clusters）的数量**。
4.  **决策**：选择该指标最大的阈值作为最佳阈值；若有多个，取第一个。

**`pgr` 的设计思路**：
`pgr` 计划将这一过程解耦为“扫描”与“决策”两步：
1.  **扫描 (`--scan`)**：显式输出区间内各阈值的指标表（TSV）。
2.  **决策**：用户既可以复刻 TreeCluster 的策略（选非单例簇最多的点），也可以结合“手肘法”或业务约束（如最大簇大小限制）选择更优的点。

### 手肘规则（Elbow Rule）与选点建议

当你不确定最佳阈值时，可以使用 `--scan` 生成数据，然后利用“手肘规则”进行辅助决策。

**什么是手肘点？**
想象你绘制一条曲线，横轴是阈值（或 K），纵轴是某个关键指标（如“簇数量”）。曲线通常呈现“L”型：
1.  **陡峭下降期**：随着阈值放松，指标迅速下降（例如大量微小的簇合并）。
2.  **平缓平台期**：指标变化很小，甚至不再变化。
3.  **拐点（手肘）**：即从“陡峭”转变为“平缓”的那个点。这个点通常对应着数据内在的自然结构。

**实战操作**：
1.  **运行扫描**：
    `pgr nwk cut tree.nwk --max-clade 0.5 --scan 0,0.5,0.01 > scan.tsv`
2.  **观察变化率**：
    - 假设 `Threshold` 从 0.05 增至 0.06 时，`Clusters` 减少了 100 个（变化剧烈）。
    - 而从 0.06 增至 0.07 时，只减少了 5 个（变化平缓）。
    - 那么 `0.06` 很可能就是最佳切点。
3.  **可视化**：将 `scan.tsv` 导入 Excel 或使用 Python/R 绘图，能更直观地识别拐点。

## 支持的模式与算法

### 1. 按簇数切 (`--k <K>`) [已实现]

等价于 R 的 `cutree(hc, k=K)` 或 SciPy 的 `fcluster(..., criterion='maxclust')`。

- **逻辑**：从根开始，优先分割高度（距最远叶子的距离）最大的节点，直到树被分割成 `K` 个子树。
- **适用场景**：你不关心阈值是多少，只想要固定数量的分组。

### 2. 按高度切 (`--height <H>`) [已实现]

等价于 R 的 `cutree(hc, h=H)` 或 SciPy 的 `fcluster(..., criterion='distance')`。

- **逻辑**：任何高度（距最远叶子的距离）大于 `H` 的节点都会被切断；高度小于等于 `H` 的节点形成簇。
- **适用场景**：适用于超度量树（Ultrametric Tree），其中高度代表时间或遗传距离。

### 3. 按根距离切 (`--root-dist <D>`) [已实现]

- **逻辑**：模拟在时间轴上的切割。从根节点出发，累积路径长度，一旦分支距离根节点的距离超过 `D` 则切断。
- **适用场景**：系统发生树分析，定义从共同祖先（根）演化特定时间后的分化群。

### 4. TreeCluster 风格：按最大簇内直径切 (`--max-clade <T>`) [已实现]

这是 **TreeCluster** 的核心算法（`Method: max_clade`）。

- **逻辑**：确保每个簇内的**最大成对距离（直径）**不超过阈值 `T`。同时隐含了单系群（Clade）约束（即簇必须是树上的完整子树）。
- **算法**：采用高效的自底向上（Bottom-Up）直径计算与自顶向下（Top-Down）贪心选择，避免了 $O(N^2)$ 的全距离矩阵计算。
- **适用场景**：病毒分型、OTU 划分等需要严格控制簇内差异度的场景。

### 5. SciPy 风格：按不一致系数切 (`--inconsistent <T>`) [已实现]

这是 SciPy `fcluster(..., criterion='inconsistent')` 的默认方法。

- **逻辑**：不一致系数（Inconsistent Coefficient）用于检测某个合并事件（节点）是否比其子树内的合并事件显著更“突兀”。
- **计算公式**：
  对于树上每个非叶子节点 $i$，考虑它以及它下方 $d$ 层（`--deep`，默认 2）内的所有子节点的合并高度集合 $H$。
  $$ I_i = \frac{height(i) - \text{mean}(H)}{\text{std}(H)} $$
  如果 $I_i > T$，则认为该节点是聚类边界，予以切断。
- **参数**：
  - `--inconsistent <T>`: 阈值，通常在 0.8 ~ 3.0 之间。
  - `--deep <D>`: 计算系数时的回溯深度，SciPy 默认为 2。
- **适用场景**：当树的整体演化速率不均匀，或者你想寻找“自然”聚类边界而不是强制切断时。

### 6. 更多 TreeCluster 变体 [规划中]

- **`avg_clade`**：簇内平均成对距离不超过阈值。
- **`med_clade`**：簇内中位数成对距离不超过阈值。
- **`single_linkage`**：树上的单链接聚类。

### 7. 支持度过滤 (`--support <S>`) [规划中]

- **逻辑**：当某条边（或节点）支持度低于阈值时，视为“不可跨越”，相当于强制切断。
- **用途**：防止聚类跨越不可靠的进化分支。

## 工作流与工具链协作

为了保持命令的专注与正交性，我们推荐以下“生成-评估”分离的工作流：

### 1. 生成 (Generation)

使用 `pgr nwk cut`：
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
# 1. 扫描不同参数，生成多个聚类结果 (规划中支持 --scan)
# pgr nwk cut input.nwk --method max-clade --scan 0.01,0.05,0.10 > partitions.tsv

# 2. 选定最佳阈值，生成最终聚类
pgr nwk cut input.nwk --method max-clade -t 0.05 > final_cluster.tsv

# 3. 可视化或提取子树
pgr nwk subset input.nwk --list final_cluster.tsv --cluster-id 1 > cluster1.nwk
```

#### 2. 层次聚类（hclust）接入
从距离矩阵出发，经由 hclust 生成树，再进行切分与评估。

```bash
# 1. 生成层次聚类树
pgr clust hier matrix.phy --method ward > tree.nwk

# 2. 切分 (按高度阈值切)
pgr nwk cut tree.nwk --height 0.05 > clusters.tsv

# 3. 评估 (计算 Cophenetic 相关系数与 Silhouette)
# pgr nwk metrics tree.nwk --part clusters.tsv --metrics silhouette > sil.tsv
```

#### 3. SciPy 风格分析 (不一致系数)
对于演化速率不均匀的树，使用不一致系数可以找到更自然的聚类边界。

```bash
# 使用不一致系数切割 (默认 depth=2)
pgr nwk cut tree.nwk --inconsistent 1.5 > clusters.tsv
```

## 与相关工具的关系与区别

### 与 R `hclust + cutree()`

- **相同点**：都是“树 → 叶子分组”。
- **不同点**：
  - `cutree()` 面向 `hclust` 产生的 dendrogram；`pgr nwk cut` 面向 Newick 树。
  - `pgr` 支持 TreeCluster 风格的生物树约束（直径、单系），性能更高。

### 与 TreeCluster

- **相同点**：目标与输出格式高度一致（叶子 → 簇）。
- **不同点**：
  - TreeCluster 是 Python 工具；`pgr` 是 Rust 实现，速度更快，且无外部依赖。
  - `pgr nwk cut` 融入了 `pgr` 工具链，可直接与 `prune`, `reroot` 等命令配合。

## 开发计划 (Roadmap)

### 第一阶段：核心功能完善 [✅ 已完成]
- [x] 实现基础切割：`--k`, `--height`, `--root-dist`.
- [x] 实现 TreeCluster 核心：`--max-clade` (diameter).
- [x] 输出格式对齐：支持 `cluster` (一行一簇) 和 `pair` (代表点-成员) 格式。
- [x] 代表点选择：支持 `root` (距离根最近), `medoid` (中心点), `first` (字母序).

### 第二阶段：高级准则与 SciPy 对齐 [✅ 部分完成]
- [x] **Inconsistent Coefficient**:
    - 实现 `calculate_inconsistency(node, depth)` 算法。
    - 添加 `--inconsistent <T>` 和 `--deep <D>` 参数。
    - 验证与 SciPy `fcluster(..., criterion='inconsistent')` 的结果一致性（因 Tie-breaking 略有差异，已添加回归测试）。
- [x] **扫描模式 (Scan Mode)**:
    - 实现 `--scan <start,end,step>` 参数。
    - 输出包含 (Threshold, ClusterCount, SingletonCount) 的摘要表。

### 第三阶段：评估与整合 [📅 待定]
- [x] 更多 TreeCluster 变体：`avg_clade` [已实现], `med_clade` [跳过].
- [ ] 整合到 `pgr clust` 统一评估流程。

### 已取消特性
- **支持度过滤 (`--support`)**：聚类树通常不带支持度信息，系统发生树场景暂不作为 MVP 核心。
