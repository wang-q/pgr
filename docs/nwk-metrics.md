# nwk metrics

`pgr nwk metrics` 的目标是：在已有 Newick 树与分组（partition）的基础上，计算“树结构友好”的聚类质量指标，用于评估分组的紧密性与分离度，补全 `pgr nwk cut` 的下游分析链。

它与 `pgr clust eval/compare` 的通用对比指标（ARI/AMI 等）互补：当没有 Ground Truth 或更关注树拓扑和分支长度时，使用 `nwk metrics` 更合适。

## 设计目标与范围

- 聚焦“树相关评估”，不引入概率模型（如高斯混合的 BIC），不做序列似然构建。
- 复用树上的距离定义：优先采用分支长度；当长度缺失时，以边数作为距离替代。
- 与 `pgr` 现有工具链风格对齐：简单的 TSV 输出、可管道化、参数直观。

参考：
- 距离计算复用 [distance.rs](file:///c:/Users/wangq/Scripts/pgr/src/cmd_pgr/nwk/distance.rs) 的逻辑与约定。
- 生成分组参考 [nwk-cut.md](file:///c:/Users/wangq/Scripts/pgr/docs/nwk-cut.md) 的输出格式与工作流。

## 输入与输出约定

### 输入

- 树文件：Newick 格式（文件或 stdin）。
- 分组文件：TSV，两列：
  - `SequenceName`：叶子标签。
  - `ClusterNumber`：簇编号；单例可为 `-1`（与 TreeCluster 约定一致）。
- 原始距离矩阵（可选，用于 cophenetic）：PHYLIP 距离矩阵（文件或 stdin），参数 `--dist <matrix.phy>`

可选参数（规划）：
- 仅评估指定簇：`--clusters <id1,id2,...>`
- 排除单例：`--exclude-singleton`

### 输出

- 全局指标（单行或少量行）：整体 silhouette 平均值、簇内直径的均值/最大值、最近簇间距离的均值/最小值、非单例簇数、单例数等。
- 全局拟合度（可选）：`CopheneticCorrelation`（当提供 `--dist` 时输出）
- 每簇指标（多行）：簇大小、簇内平均两两距离、簇直径、该簇到其他簇的最近距离、该簇的 silhouette 均值/中位数。
- 格式：TSV，字段固定，便于与其它命令串联。

## 指标定义

设树上任意两个叶子的距离为 `d(x, y)`。当分支长度存在时为长度和；否则为边数。

### 簇内紧密性（Cohesion）

- 簇内平均两两距离：`mean(d_intra) = 平均{ d(x, y) | x, y ∈ C, x≠y }`
- 簇直径：`diameter(C) = max{ d(x, y) | x, y ∈ C }`

### 簇间分离度（Separation）

- 最近簇间距离：`d_min_inter(Ci) = min_{j≠i} min_{x∈Ci, y∈Cj} d(x, y)`

### Silhouette（基于树距离）

对样本 `x`：
- `a(x)`：与同簇其他成员的平均距离
- `b(x)`：对所有其他簇，取“与该簇所有成员的平均距离”的最小值
- `s(x) = (b(x) - a(x)) / max(a(x), b(x))`

聚合：
- 全局 silhouette：`mean_x s(x)`
- 每簇 silhouette：均值与中位数
- 单例处理：单例的 `a(x)=0`；其 `s(x)` 仅由 `b(x)` 决定。可通过 `--exclude-singleton` 在全局统计中剔除。

### Cophenetic 相关系数（树拟合度）

衡量树对原始距离矩阵的保真度。给定原始距离矩阵 `D(i,j)` 与从树派生的 cophenetic 距离 `C(i,j)`，计算两者的皮尔逊相关系数：

- `C(i,j)`：树的 cophenetic 距离（在层次树上等价于两样本首次合并的高度；在一般 Newick/非超度量场景下，采用树上路径距离作为近似）
- 输出：`CopheneticCorrelation`（`r ∈ [-1, 1]`，实践中期望接近 1）
- 依赖：需要提供原始 PHYLIP 距离矩阵 `--dist <matrix.phy>`

## 计算流程（概述）

1. 解析树与分组；建立叶名到节点的映射。
2. 为每个簇：
   - 计算簇内平均两两距离与直径（必要时按簇内叶集做距离查询与缓存）。
   - 计算与其他簇的最近簇间距离（最小跨簇叶-叶距离）。
3. 为每个样本：
   - 计算 `a(x)` 与 `b(x)` 并得到 `s(x)`；聚合为全局与每簇的 silhouette。
4. 若提供 `--dist`：计算 `CopheneticCorrelation` 并输出。
5. 输出 TSV：全局指标 + 每簇指标；数值统一格式化（六位小数，去除尾随零）。

## 典型用法（规划）

```bash
# 评估分组（排除单例），输出全局与每簇指标
pgr nwk metrics tree.nwk --part clusters.tsv --exclude-singleton > metrics.tsv

# 仅评估指定簇，并输出 silhouette 相关列
pgr nwk metrics tree.nwk --part clusters.tsv --clusters 1,3,5 --metrics silhouette > sil.tsv

# 计算 cophenetic 相关系数（需原始 PHYLIP 距离矩阵）
pgr nwk metrics tree.nwk --dist matrix.phy --metrics cophenet > fit.tsv
```

## 与工具链协作（推荐工作流）

- 生成：使用 `pgr nwk cut` 扫描参数并导出分组（见 [nwk-cut.md](file:///c:/Users/wangq/Scripts/pgr/docs/nwk-cut.md)）。
- 评估：
  - 无 Ground Truth：用 `pgr nwk metrics` 查看树上的紧密性与分离度，选择合理参数（如手肘、平台点）。
  - 有 Ground Truth：用 `pgr clust eval/compare` 计算 ARI/AMI/V-Measure 等通用指标。
- 可视化与下游：用 `pgr nwk subset / to-dot / to-forest` 挑选簇并查看结构。

## 实现备注（规划）

- 距离来源：复用 `Tree::get_distance`，当长度为 0 时以边数替代（参见 [distance.rs](file:///c:/Users/wangq/Scripts/pgr/src/cmd_pgr/nwk/distance.rs) 中的约定）。
- 性能策略：小规模直接查询；中大规模按簇内叶集做批量/缓存；必要时引入并行。
- 数值格式：统一到六位小数，移除尾随零，便于比较与可视化。
- cophenetic 计算：在层次树/超度量树上取首次合并高度；在一般 Newick/非超度量场景下采用树上路径距离近似，输出相关系数以衡量保真度。
