# clust hier: 层次聚类

`pgr clust hier`（别名 `hclust`）提供通用的层次聚类（dendrogram）生成能力，支持 `single/complete/average/ward.D2` 等方法，输出 Newick 形式，便于下游 `nwk cut`。

## 背景与定位
- **归属**：`clust` 模块，与 `k-medoids`、`mcl` 等并列。
- **目标**：统计意义的 dendrogram（合并高度表达链接准则的代价），不强制“演化/分子钟”语义。
- **与 pgr 现有能力协同**：
  - 构树：`clust upgma`（有根、超度量）与 `clust nj`（加性、无根）已存在
  - 切分：`docs/nwk-cut.md` 的切树分组
  - 评估：`docs/nwk-metrics.md` 的树上指标（silhouette/直径/最近簇间距）

## 与 UPGMA/NJ 的关系
- 共同点：都以距离矩阵为输入，输出树状结构；均可配合 `nwk cut` 得到扁平分组。
- 与 UPGMA 的关系：
  - R `hclust(method="average")` 等价“平均链接”；UPGMA 是在“超度量（分子钟）”假设下的专用版本，输出有根且严格超度量的树，分支长度具有“时间/演化”意义。
  - 结论：两者链接更新一致，但语义不同；UPGMA 更偏系统发育场景，`clust hier` 更偏统计聚类。
  - 参考实现：CLI [upgma](file:///c:/Users/wangq/Scripts/pgr/src/cmd_pgr/clust/upgma.rs)，库 [clust::upgma](file:///c:/Users/wangq/Scripts/pgr/src/libs/clust/upgma.rs)
- 与 NJ 的关系：
  - NJ（Neighbor-Joining）通过 Q 矩阵最小化总树长，生成“加性最短树”，不属于链接更新范式，输出通常为无根树。
  - 在一般加性距离下，NJ比UPGMA更鲁棒；若距离是超度量，UPGMA/hclust-average与NJ在拓扑上通常一致（无根视角）。
  - 参考实现：CLI [nj](file:///c:/Users/wangq/Scripts/pgr/src/cmd_pgr/clust/nj.rs)，库 [clust::nj](file:///c:/Users/wangq/Scripts/pgr/src/libs/clust/nj.rs)

## 方法与算法要点
- `single/complete/average`：标准链接更新（Lance–Williams 框架），合并高度为链接准则对应的距离/代价。
- `ward.D2`：
  - 概念：最小化簇内平方误差（总类内方差，SSE）的增加量；常用且效果稳健。
  - 更新（平方距离版本，n为簇大小）：
    - 设合并簇 `u∪v` 与第三簇 `w` 的平方距离：
    - `d(u∪v,w)^2 = [ (n_u+n_w) d(u,w)^2 + (n_v+n_w) d(v,w)^2 − n_w d(u,v)^2 ] / (n_u+n_v+n_w)`
  - 若输入是非平方距离：可先平方进行更新，合并高度需要时取平方根或按 SSE 增量定义输出。
  - 距离前提：理论上要求欧氏或近欧氏距离；在一般生物学距离上可用，但统计意义的“方差最小化”解释会削弱。

## 输出与约定
- 输出 Newick dendrogram：
  - 分支长度表示合并高度（链接代价或 SSE 增量的相应量纲处理）。
  - 不保证严格 ultrametric（除非数据满足相应条件），但满足 `nwk cut --height` 的使用需求。
- 数值格式：统一六位小数，去除尾随零；与 `nwk distance` 的约定一致（见 [distance.rs](file:///c:/Users/wangq/Scripts/pgr/src/cmd_pgr/nwk/distance.rs)）。

## 推荐工作流
- 生成树：
  - 近分子钟/超度量场景：`clust upgma` 输出有根超度量树
  - 一般加性距离场景：`clust nj`
  - 通用层次聚类分析或需要 `ward.D2`：`clust hier --method ward.D2`
- 切分与评估：
  - 切分：`pgr nwk cut --height H` 或按 TreeCluster 风格阈值/约束
  - 无 Ground Truth：`pgr nwk metrics`（silhouette/直径/最近簇间距）
  - 有 Ground Truth：`pgr clust eval/compare`（ARI/AMI/V-Measure）
- 参考文档：
  - 切分：[nwk-cut.md](file:///c:/Users/wangq/Scripts/pgr/docs/nwk-cut.md)
  - 评估：[nwk-metrics.md](file:///c:/Users/wangq/Scripts/pgr/docs/nwk-metrics.md)

## 实现规划与优化分析 (Implementation & Optimization)

参考 `scikit-learn` 实现，`pgr clust hier` 的实现应关注以下优化点，以支撑大规模生物数据分析。

### 核心数据结构优化
- **Heap (堆)**：
  - 适用：Ward, Average, Complete, Weighted Linkage。
  - 原理：与其在每次迭代中 $O(N^2)$ 扫描距离矩阵寻找最小值，不如维护一个优先队列（Min-Heap）。虽然更新距离仍需 $O(N)$，但查找最小值降为 $O(1)$，总体性能在稀疏或受限场景下更优。
  - `scikit-learn` 参考：[`scikit-learn-main/sklearn/cluster/_hierarchical_fast.pyx`](file:///c:/Users/wangq/Scripts/pgr/scikit-learn-main/sklearn/cluster/_hierarchical_fast.pyx) 中的堆实现。
- **MST (最小生成树)**：
  - 适用：**Single Linkage** (最近邻)。
  - 原理：Single Linkage 聚类等价于求最小生成树（MST）。使用 Prim 或 Kruskal 算法可在 $O(N^2)$ (稠密) 或 $O(E \log E)$ (稀疏) 内完成，显著快于通用 Linkage 的 $O(N^3)$。
  - `scikit-learn` 参考：[`scikit-learn-main/sklearn/cluster/_agglomerative.py`](file:///c:/Users/wangq/Scripts/pgr/scikit-learn-main/sklearn/cluster/_agglomerative.py) 中的 `_single_linkage_tree` 函数。
- **Union-Find (并查集)**：
  - 配合 MST 使用，用于快速合并簇和标记标签。

### 空间与时间复杂度权衡
- **稠密矩阵 (Dense Matrix)**：
  - 现状：`pgr` 目前主要处理 PHYLIP 距离矩阵，属于稠密矩阵。
  - 策略：对于 $N < 10,000$，朴素的 $O(N^2)$ 存储和 $O(N^3)$ 计算是可接受的（且利于 SIMD 优化）。
  - 优化：对于更大规模，必须避免全矩阵存储。
- **稀疏/受限连接 (Connectivity Constraints)**：
  - 场景：图像像素聚类或基于 KNN 图的聚类。
  - `scikit-learn` 特性：支持 `connectivity` 参数（稀疏矩阵），限制只有相邻节点可以合并。这能将复杂度从 $O(N^3)$ 降至 $O(N \log N)$ 甚至 $O(N)$。
  - `pgr` 规划：未来可支持从 `pair.tsv`（稀疏边列表）直接构建 Linkage，而不强制转为全距离矩阵，从而支持超大规模序列聚类。

## 现有 Rust 生态参考
- **kodama** ([`kodama-master/`](file:///c:/Users/wangq/Scripts/pgr/kodama-master/))：
  - 实现了现代层次聚类算法（NN-chain），性能对标 `fastcluster`。
  - 核心接口 `linkage` 接受 Condensed Matrix（上三角压缩），输出 Stepwise Dendrogram。
  - 提供了完整的 `Method` 枚举（Single, Complete, Average, Ward 等）。
  - **决策**：`pgr` 将参考 `kodama` 的 NN-chain 算法实现自己的逻辑，保持对核心数据结构的完全控制（如适配稀疏输入）。
  - **价值**：利用 `kodama` 的测试用例（[`kodama-master/tests/`](file:///C:\Users\wangq\Scripts\pgr\kodama-master\src\test.rs)）和基准测试（[`kodama-master/benches/`](file:///c:/Users/wangq/Scripts/pgr/kodama-master/benches/)）来验证 `pgr` 实现的正确性与性能。
- **linfa-hierarchical** ([`linfa-master/algorithms/linfa-hierarchical/`](file:///c:/Users/wangq/Scripts/pgr/linfa-master/algorithms/linfa-hierarchical/))：
  - 提供了符合 `linfa` 生态的 `Transformer` 接口。
  - 内部直接调用 `kodama`，并增加了对 Similarity Kernel 的支持（自动转为 Distance）。
  - **借鉴**：参考其清晰的参数校验（`ParamGuard`）和从 Stepwise Dendrogram 到 Flat Clusters 的后处理逻辑（[`linfa-hierarchical/src/lib.rs`](file:///c:/Users/wangq/Scripts/pgr/linfa-master/algorithms/linfa-hierarchical/src/lib.rs) 中的 `clusters` HashMap 维护）。

### 阶段性实现路线

#### Phase 1: MVP (Primitive Implementation) - **已完成 (Completed)**
- **状态**：已在 `src/libs/clust/hier.rs` 中实现。
- **特性**：
  - 实现了基于 `CondensedMatrix`（上三角压缩矩阵）的存储，节省 50% 内存。
  - 实现了通用的 Lance-Williams 更新公式，支持 7 种 Linkage 方法：
    - `Single`, `Complete`, `Average` (UPGMA), `Weighted` (WPGMA), `Centroid` (UPGMC), `Median` (WPGMC), `Ward` (Ward's D2)。
  - 复杂度：$O(N^3)$ 时间，$O(N^2)$ 空间。
  - 验证：单元测试覆盖了各种 Linkage 方法的正确性（对比预期的新簇距离）。
- **代码位置**：
  - 核心逻辑：[`src/libs/clust/hier.rs`](file:///c:/Users/wangq/Scripts/pgr/src/libs/clust/hier.rs)
  - 矩阵结构：[`src/libs/pairmat.rs`](file:///c:/Users/wangq/Scripts/pgr/src/libs/pairmat.rs) (CondensedMatrix)

#### Phase 2: 性能优化 (NN-chain) - **待办 (Planned)**
- **目标**：将时间复杂度降至 $O(N^2)$，支撑 $N \approx 5000$ 到 $10000$ 的数据。
- **算法**：NN-chain (Nearest-neighbor chain) 算法。
  - **适用性**：Ward, Average, Complete, Weighted (要求空间具有可还原性/Reducibility，Single Linkage 不适用此法但可用 MST)。
  - **核心逻辑**：
    1.  维护一条“最近邻链”：$x_1 \to x_2 \to ... \to x_k$，其中 $x_{i+1} = \text{NN}(x_i)$。
    2.  当链出现“互为最近邻” ($x_{k-1} = \text{NN}(x_k)$ 且 $x_k = \text{NN}(x_{k-1})$) 时，合并 $(x_{k-1}, x_k)$。
    3.  合并后，从链中移除这两个点，继续寻找新的最近邻。
  - **优势**：避免了全局搜索最小值，利用局部性原理加速。
- **参考**：直接参考 `kodama` 的实现逻辑，特别是其 `chain.rs` 和 `linkage.rs` 中的状态管理。

#### Phase 3: 大规模扩展 (Sparse/Linear) - **待办 (Planned)**
- **目标**：支持超大规模数据 ($N > 10000$)。
- **策略**：
  - **稀疏输入**：不再构建全/压缩矩阵，直接基于邻接表 (`HashMap<usize, Vec<(usize, f32)>>`) 进行计算。
  - **算法**：
    - Single Linkage: 使用 MST (Prim/Kruskal) + Union-Find。
    - 其他 Linkage: 引入 Connectivity Constraints，仅计算稀疏图上的连通分量。
  - **内存优化**：实现 $O(N)$ 内存的 SLINK/CLINK 算法。

## 校验与提示：
  - 若方法为 `average` 且用户预期“演化意义”→提示使用 `upgma` 更合适。
  - 若方法为 `ward.D2` 且输入非欧氏距离→提示统计解释的偏差风险。

## 保留 UPGMA 的原因
- 语义清晰：有根、超度量、分支长度可解释为时间/演化距离。
- 生物流程稳定：与系统发育工具链更自然协作（`upgma/nj → nwk cut → nwk metrics`）。
- 用户认知与可用性：独立入口降低心智负担，避免 `method` 选择歧义。

## CLI 设计（规划）

### 命令概览
- 名称：`pgr clust hier`（可见别名 `hclust`、`hc`、`linkage`）
- 作用：从距离矩阵生成层次聚类树（dendrogram），输出为 Newick，便于后续 `nwk cut`。
- 归属：`clust` 模块，与 `k-medoids` 等并列。

### 输入
- 矩阵文件：PHYLIP 距离矩阵（标准或宽松格式）
- 格式转换：若手头是 pair TSV（三列 `name1  name2  distance`），请先使用 `pgr mat to-phylip` 转换为 PHYLIP；统一入口减少歧义，便于与 `clust upgma/nj` 一致
- 名称来源：自动从输入解析；无需额外标签文件
- 复杂度：Phase 1 朴素实现 O(n³)；Phase 2 优化至 O(n²)

### 主要参数
- `--method {single|complete|average|weighted|centroid|median|ward|ward.D2}`：链接/准则选择（默认 `ward.D2`）。命名与 SciPy linkage 对齐：
  - `single`：最近点（Nearest）
  - `complete`：最远点（Farthest/Voor Hees）
  - `average`：UPGMA（算术平均）
  - `weighted`：WPGMA（加权平均）
  - `centroid`：UPGMC（质心距离，欧氏）
  - `median`：WPGMC（质心平均）
  - `ward` / `ward.D2`：Ward 方差最小化（欧氏）
- `--outfile/-o <path>`：输出文件路径（默认 `stdout`，即打印到屏幕）。如需写入文件，可用 `-o tree.nwk` 或使用 `> tree.nwk` 重定向。
- `--optimal-ordering`：启用叶序优化，使相邻叶的距离之和最小，提升树的直观性（参考 SciPy linkage 的 `optimal_ordering`）
- `--tie {alpha|size}`：并列时的确定性规则（默认 `alpha`）
  - `alpha`：按名称字典序打破并列
  - `size`：先比较簇大小，再按名称

### 输出
- 默认输出：Newick dendrogram，分支长度表示合并高度
- 数值格式：统一六位小数、移除尾随零；与 `nwk distance` 的约定一致（见 [`src/cmd_pgr/nwk/distance.rs`](file:///c:/Users/wangq/Scripts/pgr/src/cmd_pgr/nwk/distance.rs)）

### 示例
```bash
# 先将 pair TSV 转为 PHYLIP
pgr mat to-phylip pairs.tsv -o matrix.phy

# Ward.D2（PHYLIP 输入，默认 Newick 输出）
pgr clust hier matrix.phy --method ward.D2 > tree.nwk

# Average/complete/single（PHYLIP 输入）
pgr clust hier matrix.phy --method average > tree.nwk

```

### 注意事项
- 距离前提：Ward.D2 理论依赖欧氏或近欧氏距离；在一般生物学距离上可用，但“最小总类内方差”的统计解释会减弱
- 语义差异：
  - hier 的合并高度是链接/准则的代价；不保证 ultrametric（除非数据满足相应条件）
  - 若需要“有根、超度量、演化意义”的分支长度，请使用 `clust upgma`；一般加性距离建议使用 `clust nj`
- 稳定性：并列合并通过 `--tie` 选项保证确定性；名称字典序作为默认 Tie-break
- 实现约定：`ward.D2` 内部自动按“平方距离”完成更新并返回“距离量纲”的分支长度；用户无需提供或区分 `D` 与 `D^2`
- 方法特性：
  - `centroid/median` 可能产生非单调的合并高度（inversion），属于算法特性；输出仍为合法 Newick，但高度的直觉性较 `average/ward` 略弱
  - `optimal-ordering` 会改变叶子的输出顺序以提升可读性，不改变树的拓扑与分支长度

## 与 SciPy 的映射与差异
- 方法映射：与 SciPy `linkage` 的 `method` 集合对齐，`ward` 等价 `ward.D2`（内部按平方距离更新）；`average` 等价 UPGMA，`weighted` 等价 WPGMA，`centroid/median` 等价 UPGMC/WPGMC。
- 输入差异：SciPy 接受“condensed 距离向量”或“观测矩阵”，pgr 统一使用 PHYLIP 距离矩阵；如需从 pair TSV 转换，请使用 `pgr mat to-phylip`。
- 输出差异：SciPy 返回 `(n-1)×4` 的 linkage 矩阵 Z；pgr 输出 Newick 树，直接用于 `nwk cut / to-dot / to-forest`。普通用户无需关心 Z；若需与 SciPy 互操作，请在 Python 端继续使用 Z 与 `fcluster/cophenet`。
- 叶序优化：`--optimal-ordering` 对齐 SciPy 的 `optimal_ordering` 行为，仅影响叶子顺序，保持拓扑与分支长度不变。
- 平切（flat clustering）：SciPy 的 `fcluster` 提供 `criterion='distance'|'maxclust'|...`；在 pgr 中分别对应 `nwk cut --height H` 与 `nwk cut --k K`，其它 `monocrit/inconsistent` 等准则暂不引入。
- 评估指标：SciPy 有 `cophenet`（共生相关系数）；pgr 建议在 `nwk metrics` 中加入 cophenetic 相关系数作为树质量评估的补充（与 silhouette/直径/最近簇间距并列）。

### 用户提示
- 新手路径（推荐）：`mat to-phylip → clust hier --method ward → nwk cut --height → nwk metrics → nwk 可视化`
- 互操作与审计：若需要逐步核对合并过程或在 Python 端进一步平切/统计，请使用 SciPy 的 linkage 矩阵与工具；pgr 侧保持 Newick 为主，减少心智负担。

### 示例映射
- SciPy linkage（Ward）:
  - Python: `Z = linkage(y, method='ward', optimal_ordering=True)`
  - pgr: `pgr mat to-phylip pairs.tsv -o matrix.phy` → `pgr clust hier matrix.phy --method ward --optimal-ordering > tree.nwk`
- SciPy fcluster（按距离平切）:
  - Python: `labels = fcluster(Z, t=0.05, criterion='distance')`
  - pgr: `pgr nwk cut tree.nwk --height 0.05 > clusters.tsv`
- SciPy fcluster（按簇数平切）:
  - Python: `labels = fcluster(Z, t=20, criterion='maxclust')`
  - pgr: `pgr nwk cut tree.nwk --k 20 > clusters.tsv`
- SciPy cophenet:
  - Python: `c, dists = cophenet(Z, Y)`
  - pgr（规划）：`pgr nwk metrics tree.nwk --metrics cophenet --dist matrix.phy > metrics.tsv`

### scikit-learn 映射
- AgglomerativeClustering (Ward):
  - Python: `model = AgglomerativeClustering(linkage='ward').fit(X)`
  - pgr: `pgr clust hier matrix.phy --method ward > tree.nwk`（需先计算距离矩阵）
- AgglomerativeClustering (Average/Complete/Single):
  - Python: `model = AgglomerativeClustering(linkage='average').fit(X)`
  - pgr: `pgr clust hier matrix.phy --method average > tree.nwk`
- 差异说明:
  - scikit-learn 侧重于直接输出聚类标签（`labels_`），`pgr` 侧重于生成树结构（Newick）。
  - 若需在 `pgr` 中获得标签，请配合 `nwk cut` 使用。

### 与工具链协作
- 构树：`pgr clust hier` → 生成 dendrogram
- 切分：`pgr nwk cut --height H` → 导出分组（也可用 TreeCluster 风格约束方法）
- 评估：
  - 无 Ground Truth：`pgr nwk metrics`（silhouette、簇内直径、最近簇间距）
  - 有 Ground Truth：`pgr clust eval/compare`（ARI/AMI/V-Measure）
- 可视化：`pgr nwk to-dot/to-forest` → 图形/LaTeX 展示
