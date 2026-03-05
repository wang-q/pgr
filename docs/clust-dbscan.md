# clust dbscan: 密度聚类与规划

`pgr clust dbscan` 提供基于密度的聚类（DBSCAN）。当前实现以“距离矩阵”为输入，输出标准的 `cluster`/`pair` 两种分区格式，便于与其它命令互操作。

## 现状概览

- 输入与输出
  - 输入：成对距离 `.tsv`（lower is better），通过 `pairmat::ScoringMatrix::from_pair_scores` 读入
  - 输出：`cluster`（一行一簇，首元素为代表点）或 `pair`（代表点-成员对）
- CLI 参数
  - `infile`、`--format {cluster|pair}`、`--same`、`--missing`、`--eps`、`--min_points`、`-o/--outfile`
- 代码参考
  - CLI：[clust/dbscan.rs](file:///c:/Users/wangq/Scripts/pgr/src/cmd_pgr/clust/dbscan.rs)
  - 算法库：[libs/clust/dbscan.rs](file:///c:/Users/wangq/Scripts/pgr/src/libs/clust/dbscan.rs)

## 规划目标

### 扫描与评分
- 新增 `--scan <start>,<end>,<steps>`：参考 `nwk cut --scan` 风格，对主参数 `eps` 进行扫描，输出 TSV 摘要表（不绘图）
- 评分列（可选）：
  - `Silhouette`（距离矩阵版）
  - `DBIndex`（Davies–Bouldin）
- 建议输出列：
  - `Epsilon`, `Clusters`, `Noise`, `Silhouette`, `DBIndex`
- 自动选优（可选）：`--opt-eps {silhouette|max-clusters|min-noise}`，在扫描后直接选点并输出该分区

**用法规划**：
`pgr clust dbscan ... --scan <start>,<end>,<steps>`
（注：扫描仅针对方法的**主阈值参数**，此处为 `eps`。例如 `min_points` 保持固定为用户指定值或默认值）

**输出指标表（示例）**：
| Epsilon | Clusters | Noise | Silhouette | DBIndex |
| :--- | :--- | :--- | :--- | :--- |
| 0.10 | 25 | 12 | 0.42 | 1.85 |
| 0.12 | 18 | 20 | 0.47 | 1.72 |
| ... | ... | ... | ... | ... |

### 代表点与比例参数
- 代表点（`cluster`/`pair` 通用）：维持现状（medoid），并在 `pair` 中输出 `(medoid, member)`
- 新增 `--min-pct <0..1>`：按样本比例折算为 `min_points`，与 `--min_points` 二选一

### 评分实现（距离矩阵版本）
- Silhouette
  - 对每个点 i：簇内平均距离 `a(i)`；到其它簇的最小平均距离 `b(i)`；`s(i) = (b-a)/max(a,b)`；总体取平均
- DBIndex
  - 每簇散度（簇内到中心的平均距离）；簇间中心距；计算最大比值并平均
- 位置建议：`libs::clust::metrics` 或 `libs::metrics`，供扫描与 `clust eval` 复用

### 互操作与职责分离
- 算法侧（本命令）：负责 DBSCAN 聚类与扫描 TSV 输出
- 评估侧（`clust eval`）：外部有效性（`ARI/AMI/V-Measure`）为主；内部有效性（`Silhouette/DBIndex`）作为补充
- 与树工具协作：不直接涉及 `nwk`，但输出的 `cluster/pair` 可用于后续评估或可视化

## 性能与边界

- 复杂度
  - DBSCAN：从距离矩阵出发，整体约 `O(N^2)`；扫描为 `steps × O(N^2)`
  - 评分计算也需 `O(N^2)`（平均距离/中心距）
- 建议
  - 缩小 `start,end`：可据距离分布的分位数设范围（如 `p10..p90`）
  - 合理的 `steps`（默认 100），规模较大时降低分辨率
  - 清晰文档提示中大型数据的计算成本

## 测试计划

- 单元测试
  - 距离矩阵版 Silhouette/DBIndex 的正确性（小矩阵、可手算）
  - `--min-pct` 与 `--min_points` 的互斥与折算
  - 噪声计数与聚类数的统计正确性
- 集成测试
  - `--scan-eps` 输出 TSV 字段与排序一致性
  - `--opt-eps silhouette` 在简单数据上的选点合理性
- Fuzz
  - 随机小矩阵，验证扫描输出与聚类稳定性

## 使用示例（规划）

```bash
# 1) 基本聚类（pairwise 距离输入）
pgr clust dbscan pairs.tsv --eps 0.15 --min_points 3 -o clusters.tsv

# 2) 扫描 eps 并输出评分曲线（TSV）
pgr clust dbscan pairs.tsv --scan 0.05,0.5,100 -o scan.tsv

# 3) 自动基于 Silhouette 选优 eps 并直接输出分区
pgr clust dbscan pairs.tsv --scan 0.05,0.5,100 --opt-eps silhouette -o best.tsv

# 4) 使用比例表达 min_points
pgr clust dbscan pairs.tsv --eps 0.15 --min-pct 0.02 -o clusters.tsv

# 5) 输出 pair 格式，便于后续评估
pgr clust dbscan pairs.tsv --eps 0.15 --min_points 3 --format pair -o pairs.out.tsv
```

## 相关参考

- 本仓库
  - CLI：[clust/dbscan.rs](file:///c:/Users/wangq/Scripts/pgr/src/cmd_pgr/clust/dbscan.rs)
  - 库：[libs/clust/dbscan.rs](file:///c:/Users/wangq/Scripts/pgr/src/libs/clust/dbscan.rs)
  - 评估文档：[clust-eval.md](file:///c:/Users/wangq/Scripts/pgr/docs/clust-eval.md)
- 生态
  - ClustEval DBSCAN 扫描与评分（Silhouette）：`clusteval/dbscan.py`

