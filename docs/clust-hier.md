# clust hier: 层次聚类

`pgr clust hier`（别名 `hclust`）提供通用的层次聚类（dendrogram）生成能力，支持 `single/complete/average/ward.D2` 等方法，输出 Newick 形式，便于下游 `nwk cut`。

## 背景与定位
- **归属**：`clust` 模块，与 `k-medoids`、`mcl` 等并列。
- **目标**：统计意义的 dendrogram（合并高度表达链接准则的代价），不强制“演化/分子钟”语义。
- **与 pgr 现有能力协同**：
  - 构树：`clust upgma`（有根、超度量）与 `clust nj`（加性、无根）已存在
  - 切分：`docs/nwk-cut.md` 的切树分组
  - 评估：`docs/nwk-eval.md` 的树上指标（几何/分类/演化/地理多维度评估）

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
  - 内部评估（无 Ground Truth）：`pgr clust eval --matrix ...` (Silhouette) 或 `pgr nwk eval` (树结构评估)
  - 外部评估（有 Ground Truth）：`pgr clust eval` (ARI/AMI/V-Measure)
- 参考文档：
  - 切分：[nwk-cut.md](file:///c:/Users/wangq/Scripts/pgr/docs/nwk-cut.md)
  - 评估：[nwk-eval.md](file:///c:/Users/wangq/Scripts/pgr/docs/nwk-eval.md)

## SciPy 实现借鉴与对比 (Insights from SciPy)
通过深入分析 `scipy.cluster.hierarchy` 源码（基于 Cython 的高性能实现），`pgr` 吸收了以下关键设计思想：

1.  **Generic Clustering Algorithm (Heap 优化)**:
    - **背景**: NN-Chain 算法仅适用于可归约方法（Ward, Average, Complete, Single, Weighted），无法处理 Centroid 和 Median。
    - **SciPy 方案**: 在 `fast_linkage` (`_hierarchy.pyx`) 中实现了 Daniel Müllner (2011) 的算法。该算法结合 `neighbor` 数组和 Binary Heap，将所有方法的复杂度统一优化至 $O(N^2 \log N)$ 甚至 $O(N^2)$。
    - **pgr 借鉴**: 目前 `pgr` 对 Centroid/Median 使用 $O(N^3)$ 朴素实现。未来计划移植该 Heap 算法，消除性能短板。

2.  **Ward 方法的数值稳定性与效率**:
    - **SciPy 实现**: `ward` 更新公式在内部计算时涉及平方和开方（`sqrt`），这在大量迭代中可能积累浮点误差，且计算开销较大。
    - **pgr 优化**: `pgr` 采用全程平方距离运算（Internal Squared Euclidean），仅在最终输出时开方。这避免了中间步骤的精度损失和 `sqrt` 开销，使得 Ward 方法的性能与 Average 方法完全持平（而在许多其他库中 Ward 通常更慢）。

3.  **生态一致性**:
    - **Flat Clustering**: `pgr nwk cut` 的设计与 SciPy `fcluster` 的 `criterion='distance'|'maxclust'` 保持概念一致。
    - **Cophenetic Correlation**: 确认将 `cophenet` 引入 `pgr nwk eval`，作为衡量树对原始距离矩阵拟合优度的核心指标。

4.  **Optimal Leaf Ordering (OLO)**:
    - **背景**: 标准层次聚类算法生成的树，左右子树的顺序是任意的。这导致在绘制热图（Heatmap）时，相似的行/列可能不相邻，视觉效果杂乱。
    - **SciPy 方案**: `scipy.cluster.hierarchy.optimal_leaf_ordering`。
    - **算法**: Bar-Joseph et al. (2001) 的动态规划算法。在不改变树拓扑结构的前提下，通过旋转内部节点，最小化相邻叶子之间的距离之和。
    - **pgr 借鉴**: 计划在 `pgr nwk order` 中实现此功能（`--olo` 或 `--optimal`），作为聚类后的标准优化步骤，显著提升下游可视化（`pgr mat plot` 或外部工具）的效果。

## 实现规划与优化分析 (Implementation & Optimization)

### 核心数据结构优化
- **Heap (堆) - Generic Clustering Algorithm**:
  - 适用：所有方法，特别是 **Centroid** 和 **Median**（不可归约，无法使用 NN-chain）。
  - 原理：维护一个距离最近邻的优先队列。这是 Daniel Müllner (2011) 提出的 "Generic Clustering Algorithm"。
  - SciPy 参考：`fast_linkage` in `_hierarchy.pyx`。
  - `pgr` 规划：作为 Phase 4 的一部分，替换目前的 Primitive $O(N^3)$ 实现，统一所有方法的性能基线。
- **MST (最小生成树)**:
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
- **状态**：已在 `src/libs/clust/hier.rs` 中实现，并集成到 CLI `src/cmd_pgr/clust/hier.rs`。
- **特性**：
  - 实现了基于 `CondensedMatrix`（上三角压缩矩阵）的存储，节省 50% 内存。
  - 实现了通用的 Lance-Williams 更新公式，支持 7 种 Linkage 方法：
    - `Single`, `Complete`, `Average` (UPGMA), `Weighted` (WPGMA), `Centroid` (UPGMC), `Median` (WPGMC), `Ward` (Ward's D2)。
  - 复杂度：$O(N^3)$ 时间，$O(N^2)$ 空间。
  - 验证：单元测试覆盖了核心算法，集成测试覆盖了 CLI 功能。

#### Phase 2: 性能优化 (NN-chain) - **已完成 (Completed)**
- **状态**：已在 `src/libs/clust/hier.rs` 中实现 NN-chain 算法。
- **特性**：
  - **算法**：NN-chain (Nearest-neighbor chain) 算法。
  - **适用性**：Ward, Average, Complete, Weighted (空间具有可还原性/Reducibility)。
  - **复杂度**：时间复杂度优化至 $O(N^2)$。
  - **自动调度**：`linkage` 函数自动根据 Method 选择最佳算法（Reducible 方法用 NN-chain，其它用 Primitive）。
  - **验证**：
    - 单元测试验证了 NN-chain 与 Primitive 算法输出的一致性（包括 ID 映射和拓扑）。
    - 基准测试证明了显著的性能提升。

**Benchmark Results (Average & Ward):**

| N | Method | Primitive $O(N^3)$ | NN-Chain $O(N^2)$ | Speedup |
|---|---|---|---|---|
| 100 | Average | ~300 µs | ~63 µs | ~4.7x |
| 200 | Average | ~2.1 ms | ~248 µs | ~8.5x |
| 400 | Average | ~15.6 ms | ~975 µs | ~16x |
| | | | | |
| 100 | Ward | ~315 µs | ~70 µs | ~4.5x |
| 200 | Ward | ~2.3 ms | ~266 µs | ~8.6x |
| 400 | Ward | ~15.8 ms | ~1.0 ms | ~15.8x |

*注：Ward Linkage 在优化后（平方距离更新）性能与 Average Linkage 几乎持平。*

#### Phase 3: 大规模数据策略 (Two-stage / Representative) - **推荐 (Recommended)**
参见 `docs/clust.md` 中的“大规模数据策略”章节。

#### Phase 4: 性能与正确性优化 (Pending)
通过分析 `kodama`、`scikit-learn` 和 `scipy` 实现，确定以下优化方向：
1.  **Generic Clustering Algorithm (Heap)**:
    - 目标：优化 **Centroid** 和 **Median** 方法。
    - 方案：参考 SciPy 的 `fast_linkage` 实现（基于 Müllner 2011），引入 Binary Heap 维护最近邻距离。这将把这两个方法的复杂度从 $O(N^3)$ 降至 $O(N^2 \log N)$。
    - 优先级：中（除非用户有大量 Centroid/Median 聚类需求）。
2.  **Ward/Centroid 平方距离优化 (已完成)**:
    - 改进：在算法开始时一次性将距离矩阵平方，使用简化版 Lance-Williams 更新，仅在输出时开方。
    - 效果：消除了每次迭代中的 `sqrt` 调用，使得 Ward Linkage 的性能与 Average Linkage 持平（基准测试证实）。
2.  **In-place 接口 (已完成)**:
    - 引入 `linkage_inplace`，允许消耗输入的 `CondensedMatrix`（避免克隆），节省 $O(N^2)$ 内存复制开销。
3.  **Chain 循环优化 (已分析)**:
    - 分析：`kodama` 使用了高效的 `ActiveList` (双向链表) 来跳过非活跃节点。
    - 结论：虽然这能将寻找最近邻的复杂度从 $O(N)$ 降至 $O(K)$，但鉴于 Condensed Matrix 的顺序访问对 CPU 缓存非常友好，且当前实现在 $N=400$ 时仅需 ~1ms，引入链表的跳跃访问可能收益有限甚至负优化（对于中小规模数据）。
    - 决策：暂不实施，直至 profiling 显示 NN 搜索成为显著瓶颈。
4.  **MST 算法 (已分析)**:
    - 分析：`scikit-learn` 和 `kodama` 对 Single Linkage 使用 MST 算法。
    - 结论：对于稠密矩阵，NN-Chain 和 Prim MST 都是 $O(N^2)$。当前的 NN-Chain 实现通用且足够高效。MST 主要优势在于处理稀疏图输入（Phase 3 范畴）。
    - 决策：维持现状。


#### Phase 5: 测试覆盖率增强 (已完成)
参考 `kodama` 和 `scikit-learn` 的测试策略，已增加以下测试以提升稳健性：
1.  **Fuzzing / Randomized Testing (Kodama)**:
    - 目标：验证 NN-Chain 算法与 Primitive 算法在大量随机输入下的一致性。
    - 状态：已实现 `test_nn_chain_fuzzing`，循环测试 20 个不同大小（$N=10 \sim 105$）的随机矩阵，验证输出 Step 数量和合并距离的一致性（包括 Ward 方法）。
2.  **Monotonicity Check (Sklearn)**:
    - 目标：验证生成的 Dendrogram 是否满足单调性（除了 Centroid/Median 方法）。
    - 状态：已实现 `test_monotonicity`，断言所有单调方法的 `steps[i].distance <= steps[i+1].distance`。
3.  **Edge Cases (Kodama)**:
    - 目标：验证极小输入的处理 ($N=0, 1, 2$)。
    - 状态：已实现 `test_edge_cases`，确保空输入或单点输入正确返回空结果。

#### Phase 6: 基准测试增强 (已完成)
参考 `kodama` 和 `scikit-learn` 的基准测试策略，已实施了以下测试：
1.  **多尺度性能曲线 (Scalability)**:
    - 验证了 NN-Chain 算法在 $N=1000 \sim 4000$ 范围内的 $O(N^2)$ 扩展性。
    - **Ward** 与 **Average** 的性能曲线几乎重合，证明了平方距离优化的有效性。
    - $N=4000$ 时耗时约 0.18s，推算 $N=20000$ 时约需 5s，完全可接受。
2.  **方法间对比 (Method Comparison)**:
    - 在 $N=1000$ 时，所有方法（Single, Complete, Average, Weighted, Ward）的耗时高度一致（~6.0ms）。
    - 表明核心算法框架的效率主导了计算，具体距离公式的差异对性能影响微乎其微。

**最新 Benchmark 数据 (Average & Ward):**

| N | Primitive $O(N^3)$ | NN-Chain $O(N^2)$ |
|---|---|---|
| 100 | ~0.3 ms | ~0.06 ms |
| 400 | ~16 ms | ~0.9 ms |
| 1000 | (未测) | ~5.3 ms |
| 2000 | (未测) | ~29.0 ms |
| 4000 | (未测) | ~174 ms |

#### Phase 7: 真实分布与效果验证 (Planned)
参考 `linfa-hierarchical` 的 `test_blobs` 测试，计划增加以下内容以验证算法的统计有效性：
1.  **真实分布测试 (Blobs Test)**:
    - 目标：验证算法能否正确聚类具有明显几何结构的合成数据（Statistical Correctness）。
    - 计划：在 `tests/` 中添加集成测试，生成两个高斯分布的簇（Blob A 和 Blob B），计算距离矩阵，运行 `clust hier`，验证生成的 Newick 树是否将两个簇的点分在不同的主分支上。
    - **注意**：由于 `pgr nwk cut` 命令已实现，该测试可以被激活。
2.  **输入预处理文档 (已完成)**:
    - 目标：澄清输入要求。
    - 状态：已在 `mat-transform.md` 和 `clust-hier.md` 中更新，并增强了 `pgr mat transform` 功能（支持对角线归一化），确保用户能正确地将相似度转换为距离。

## CLI 设计

### 命令概览
- 名称：`pgr clust hier`（可见别名 `hclust`、`hc`、`linkage`）
- 作用：从距离矩阵生成层次聚类树（dendrogram），输出为 Newick，便于后续 `nwk cut`。
- 归属：`clust` 模块，与 `k-medoids` 等并列。

### 输入
- 矩阵文件：PHYLIP 距离矩阵（标准或宽松格式）
- 格式转换：若手头是 pair TSV（三列 `name1  name2  distance`），请先使用 `pgr mat to-phylip` 转换为 PHYLIP；统一入口减少歧义，便于与 `clust upgma/nj` 一致。
- 距离/相似度转换：`clust hier` 仅接受**距离矩阵**（越小越相似）。如果输入是相似度矩阵（如 BLAST Identity, Alignment Score），请先使用 `pgr mat transform` 进行转换（如 `--op inv-linear --max 100` 或 `--op log`）。
- 名称来源：自动从输入解析；无需额外标签文件

### 主要参数
- `--method {single|complete|average|weighted|centroid|median|ward}`：链接/准则选择（默认 `ward`）。命名与 SciPy linkage 对齐。
- `--outfile/-o <path>`：输出文件路径（默认 `stdout`，即打印到屏幕）。

### 输出
- 默认输出：Newick dendrogram，分支长度表示合并高度
- 数值格式：统一六位小数、移除尾随零；与 `nwk distance` 的约定一致（见 [`src/cmd_pgr/nwk/distance.rs`](file:///c:/Users/wangq/Scripts/pgr/src/cmd_pgr/nwk/distance.rs)）

### 示例
```bash
# 先将 pair TSV 转为 PHYLIP
pgr mat to-phylip pairs.tsv -o matrix.phy

# Ward (PHYLIP 输入，默认 Newick 输出)
pgr clust hier matrix.phy --method ward > tree.nwk

# Average/complete/single (PHYLIP 输入)
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
- 叶序优化：`pgr` 推荐 `pgr nwk order --nd` (Ladderize) 以换取极高的性能，且可视化效果通常足够好。
- 平切（flat clustering）：SciPy 的 `fcluster` 提供 `criterion='distance'|'maxclust'|...`；在 pgr 中分别对应 `nwk cut --height H` 与 `nwk cut --k K`，其它 `monocrit/inconsistent` 等准则暂不引入。
- 评估指标：SciPy 有 `cophenet`（共生相关系数）；pgr 建议在 `nwk eval` 中加入 cophenetic 相关系数作为树质量评估的补充。
- 内部实现细节：
  - SciPy 的 `ward` 更新公式在内部进行平方和开方（`sqrt`）；`pgr` 采用全程平方距离运算（仅输出时开方），避免了中间步骤的精度损失和开方开销，理论上更高效。
  - SciPy 的 `fast_linkage` 使用了 Heap 优化；`pgr` 目前对非 NN-chain 方法使用朴素实现，未来可借鉴此优化。

### 用户提示
- 新手路径（推荐）：`mat to-phylip → clust hier --method ward → nwk cut --height → nwk metrics → nwk 可视化`
- 互操作与审计：若需要逐步核对合并过程或在 Python 端进一步平切/统计，请使用 SciPy 的 linkage 矩阵与工具；pgr 侧保持 Newick 为主，减少心智负担。

### 示例映射
- SciPy linkage（Ward）:
  - Python: `Z = linkage(y, method='ward', optimal_ordering=True)`
  - pgr: `pgr mat to-phylip pairs.tsv -o matrix.phy` → `pgr clust hier matrix.phy --method ward > tree.nwk` → `pgr nwk order tree.nwk --nd > ordered.nwk`
- SciPy fcluster（按距离平切）:
  - Python: `labels = fcluster(Z, t=0.05, criterion='distance')`
  - pgr: `pgr nwk cut tree.nwk --height 0.05 > clusters.tsv`
- SciPy fcluster（按簇数平切）:
  - Python: `labels = fcluster(Z, t=20, criterion='maxclust')`
  - pgr: `pgr nwk cut tree.nwk --k 20 > clusters.tsv`
- SciPy cophenet:
  - Python: `c, dists = cophenet(Z, Y)`
  - pgr: `pgr nwk eval tree.nwk --dist matrix.phy > metrics.tsv`

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
- 切分：`pgr nwk cut --height H` → 导出分组
- 评估：
  - 无 Ground Truth：`pgr nwk eval`（几何/分类/演化/地理多维度评估）
  - 有 Ground Truth：`pgr clust eval`（ARI/AMI/V-Measure）
- 可视化：`pgr nwk to-dot/to-forest` → 图形/LaTeX 展示
