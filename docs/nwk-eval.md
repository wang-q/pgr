# nwk eval [设计中]

`pgr nwk eval` 的目标是建立一个**多维度的树评估框架**。除了基于树拓扑的几何指标外，重点引入**生物学语境**，评估基因树（Gene Tree）与物种树（Species Tree）或分类单元（Taxonomy）的一致性。

它旨在回答四个核心问题：
1.  **几何紧密性**：从图论角度，分组是否紧密且分离？（Silhouette, Diameter）
2.  **分类一致性**：分组是否对应自然的生物分类单元（如科、属、种）？（Taxonomic Purity）
3.  **演化一致性**：基因树的分组与公认的物种演化历史是否冲突？（Discordance）
4.  **地理/性状一致性**：分组是否对应特定的地理区域或表型特征？（Trait Purity/Entropy）

## 设计目标与范围

- **多层次评估**：支持从纯数学指标到引入外部生物学知识（分类表、参考树、地理分布）的深度评估。
- **基因树 vs 物种树**：承认并量化由 ILS、HGT 或基因重复/丢失引起的不一致性。
- **性状映射**：支持将叶子节点映射到分类层级或地理区域，计算分组的“纯度”和“熵”。
- **距离定义**：优先采用分支长度（Patristic distance）；当长度缺失时，以边数作为距离替代。

## 输入与输出约定

### 输入

- **Target Tree (`-t`)**: 待评估的树（通常是基因树）。Newick 格式。
- **Partition (`-p`)**: 待评估的分组（可选；若不提供，则评估整棵树或根据 `--cut` 动态生成）。支持 `cluster` (列表) 和 `pair` (对) 格式。
- **External Context (可选)**:
  - **Trait Map (`--traits`)**: TSV 文件，用于分类或地理信息。格式 `LeafName <tab> Trait1,Trait2...`。
  - **Reference Tree (`--ref`)**: 参考树（通常是物种树），用于计算拓扑差异。
  - **Original Matrix (`--dist`)**: 原始 PHYLIP 距离矩阵，用于计算 Cophenetic 相关性。

### 输出

- **TSV 格式**，包含多组列：
  - `Basic`: Size, Diameter, AvgDist.
  - `Geom`: Silhouette, Separation.
  - `Trait`: Purity, Entropy, DominantTrait. (复用分类学指标逻辑)
  - `Phylo`: RF-Distance (to Ref), ConflictScore.
  - `Fit`: CopheneticCorrelation.

## 指标详细定义

### 1. 几何/拓扑指标 (Geometric/Topological)
*无需外部信息，仅基于输入树的边长和拓扑。*

设树上任意两个叶子的距离为 `d(x, y)`。

#### 簇内紧密性（Cohesion）
- **簇内平均两两距离**: `mean(d_intra) = 平均{ d(x, y) | x, y ∈ C, x≠y }`
- **簇直径**: `diameter(C) = max{ d(x, y) | x, y ∈ C }`

#### 簇间分离度（Separation）
- **最近簇间距离**: `d_min_inter(Ci) = min_{j≠i} min_{x∈Ci, y∈Cj} d(x, y)`

#### Silhouette（基于树距离）
对样本 `x`：
- `a(x)`：与同簇其他成员的平均距离
- `b(x)`：对所有其他簇，取“与该簇所有成员的平均距离”的最小值
- `s(x) = (b(x) - a(x)) / max(a(x), b(x))`
- **聚合**: 计算全局平均值和每簇的均值/中位数。单例的 `s(x)` 仅由 `b(x)` 决定（`a(x)=0`）。

#### Cophenetic 相关系数（树拟合度）
衡量树结构对原始距离矩阵的保真度。需提供 `--dist`。
- **定义**: 树上任意两叶子节点的距离（即 LCA 高度的 2 倍，或路径长度）称为 Cophenetic 距离。
- **计算**: 计算原始距离矩阵 $D$ 与 Cophenetic 距离矩阵 $C$ 之间的 Pearson 相关系数 $r$。
- **SciPy 参考**: `scipy.cluster.hierarchy.cophenet`。
- **应用**: 评估不同聚类方法（如 UPGMA vs NJ vs Ward）对数据的拟合优度。$r$ 越接近 1，表示树结构越能真实反映原始数据。

### 2. 性状/分类指标 (Trait/Taxonomic Consistency)
*需提供 `--traits` (或兼容 `--tax`)。评估分组与外部标签（分类、地理、表型）的一致性。*

- **Purity**: 簇内最优势标签的占比。
  - 分类示例：9 个 *E. coli*，1 个 *S. enterica* -> Purity = 0.9。
  - 地理示例：9 个 *Asia*，1 个 *Europe* -> Purity = 0.9。
- **Entropy**: 标签分布的香农熵。`H(C) = - sum(p_i * log2(p_i))`。衡量簇内标签的混乱程度。
- **LCA Rank Consistency** (仅限分类): 如果提供层级信息，评估 LCA 是否对应特定层级。

### 3. 系统发育指标 (Phylogenetic Discordance)
*需提供 `--ref`。评估基因树局部结构与物种树的冲突。*

- **Local RF Distance**: 簇内子树与参考树对应子集的 Robinson-Foulds 距离。
- **Monophyly Check**: 基因树上的簇成员，在物种树上是否也聚集成单系群？
  - 若基因树聚类但物种树分散 -> 可能暗示 HGT 或 LBA（长枝吸引）。

## 典型用法 (Use Cases)

```bash
# 场景 A: 纯几何评估 (无外部信息)
pgr nwk eval tree.nwk --part clusters.tsv > geom_eval.tsv
# 输出: ClusterID, Size, Silhouette, Diameter

# 场景 B: 性状/地理一致性验证
# traits.tsv: LeafName <tab> Region
pgr nwk eval tree.nwk --part clusters.tsv --traits location.tsv > geo_eval.tsv
# 输出: ..., Purity, DominantTrait, Entropy

# 场景 C: 基因树质量控制 (Reference Tree)
pgr nwk eval gene_tree.nwk --ref species_tree.nwk > phylo_eval.tsv
# 输出: Global_RF_Dist, Cluster_Conflict_Score

# 场景 D: 原始距离拟合度 (Cophenetic)
pgr nwk eval tree.nwk --dist matrix.phy --metrics cophenet > fit.tsv
```

## 实现备注（技术细节）

- **代码复用与协同**：
  - **距离计算**：
    - 基础：复用 `Tree::get_distance` (来自 [distance.rs](file:///c:/Users/wangq/Scripts/pgr/src/cmd_pgr/nwk/distance.rs))。
    - 优化：对于 `avg_clade` (AvgDist) 和 `max_clade` (Diameter)，直接复用 `libs/phylo/tree/stat.rs` 中的 $O(N)$ 自底向上聚合算法（已在 `pgr nwk cut` 中验证）。
  - **拓扑比较**：
    - RF 距离：核心逻辑应复用或提取自 [cmp.rs](file:///c:/Users/wangq/Scripts/pgr/src/cmd_pgr/nwk/cmp.rs)。建议将 `compute_metrics` 封装为 `libs::phylo::cmp` 模块。
  - **单系性检查**：
    - 复用 `Tree::is_monophyletic` (已在 [label.rs](file:///c:/Users/wangq/Scripts/pgr/src/cmd_pgr/nwk/label.rs) 和 [subtree.rs](file:///c:/Users/wangq/Scripts/pgr/src/cmd_pgr/nwk/subtree.rs) 中使用)。

- **性能策略**：
  - 优先使用基于遍历的聚合算法，避免构建 $O(N^2)$ 全距离矩阵。
  - 对于 Silhouette，针对大树考虑采样。

- **数值格式**：统一到六位小数，移除尾随零。

## 实施计划 (Roadmap)

### Phase 1: 几何核心 (Geometric Core)
- [ ] **CLI 搭建**: 支持 `-t`, `-p`, `--dist`。
- [ ] **核心指标**: Size, Diameter, AvgDist, MinInterDist。
- [ ] **Silhouette**: 实现树上距离计算逻辑。
- [ ] **Cophenetic**: 实现 Pearson 相关系数。

### Phase 2: 分类学扩展 (Taxonomic Extension)
- [ ] 解析 `--tax` 文件（KV 映射）。
- [ ] 实现 `Purity` 和 `Entropy` 指标。

### Phase 3: 参考树对比 (Reference Comparison)
- [ ] 引入 `phylo::cmp` 模块。
- [ ] 实现 RF 距离和单系性检查。

### Phase 4: 文档与高级功能
- [ ] 完善文档，添加 Benchmark。
- [ ] 支持 NCBI Taxonomy Dump。
