# pgr clust dbscan

`pgr clust dbscan` 提供基于密度的聚类（DBSCAN）。当前实现以“距离矩阵”为输入，输出标准的 `cluster`/`pair` 两种分区格式，便于与其它命令互操作。

## 现状概览

- 输入与输出
  - 输入：成对距离 `.tsv`（lower is better）
  - 输出：`cluster`（一行一簇，首元素为代表点）或 `pair`（代表点-成员对）
- CLI 参数
  - `infile`、`--format {cluster|pair}`、`--same`、`--missing`、`--eps`、`--min-points`、`-o/--outfile`

## 规划

`--scan`/`--opt-eps`/`--min-pct` 等参数扫描与评分功能尚未实现。

## 使用示例

```bash
# 基本聚类（pairwise 距离输入）
pgr clust dbscan pairs.tsv --eps 0.15 --min-points 3 -o clusters.tsv

# 输出 pair 格式，便于后续评估
pgr clust dbscan pairs.tsv --eps 0.15 --min-points 3 --format pair -o pairs.out.tsv
```

## 相关参考

- 评估文档：[clust-eval.md](clust-eval.md)

