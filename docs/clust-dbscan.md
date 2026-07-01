# clust dbscan: 密度聚类

`pgr clust dbscan` 提供基于密度的聚类（DBSCAN）。当前实现以“距离矩阵”为输入，输出标准的 `cluster`/`pair` 两种分区格式，便于与其它命令互操作。

## 现状概览

- 输入与输出
  - 输入：成对距离 `.tsv`（lower is better），通过 `pairmat::ScoringMatrix::from_pair_scores` 读入
  - 输出：`cluster`（一行一簇，首元素为代表点）或 `pair`（代表点-成员对）
- CLI 参数
  - `infile`、`--format {cluster|pair}`、`--same`、`--missing`、`--eps`、`--min_points`、`-o/--outfile`
- 代码参考
  - CLI：[clust/dbscan.rs](../src/cmd_pgr/clust/dbscan.rs)
  - 算法库：[libs/clust/dbscan.rs](../src/libs/clust/dbscan.rs)

## 规划

`--scan`/`--opt-eps`/`--min-pct` 等尚未实现的参数扫描与评分功能规划详情移至 [notes/design/dbscan-planned.md](../notes/design/dbscan-planned.md)。

## 使用示例

```bash
# 基本聚类（pairwise 距离输入）
pgr clust dbscan pairs.tsv --eps 0.15 --min_points 3 -o clusters.tsv

# 输出 pair 格式，便于后续评估
pgr clust dbscan pairs.tsv --eps 0.15 --min_points 3 --format pair -o pairs.out.tsv
```

## 相关参考

- 本仓库
  - CLI：[clust/dbscan.rs](../src/cmd_pgr/clust/dbscan.rs)
  - 库：[libs/clust/dbscan.rs](../src/libs/clust/dbscan.rs)
  - 评估文档：[clust-eval.md](clust-eval.md)
- 生态
  - ClustEval DBSCAN 扫描与评分（Silhouette）：`clusteval/dbscan.py`

