# `rg runlist` 命令行基准：pgr（SpanIndex 二分）vs rgr（IntSpan 交集）

> 目的：对比 `pgr rg runlist`（按 runlist 过滤 `.rg` 行）与外部 `rgr
> runlist` 的耗时与内存。2026-08-04 实测。

## 数据

* runlist `rl.json`：8 条染色体 × 100 Mb，154k spans（约覆盖 30%）
* `target.20k.rg`：20,000 条随机查询区间

## 正确性验证

`overlap` 与 `superset` 两种 op 在 20k 目标上，pgr 与 rgr 输出 `sort` 后
`diff` 为空（逐行一致）。

## 结果（20k 目标，5 次取均值）

| op | pgr `rg runlist` | rgr `runlist` | 加速 | RSS（pgr / rgr） |
| :--- | ---: | ---: | ---: | ---: |
| overlap | 14.6 ± 0.4 ms | 1.210 ± 0.004 s | **~83×** | 15.7 / 9.8 MB |
| superset | 15.0 ± 0.6 ms | 8.595 ± 0.025 s | **~588×** | — |

## 优化说明

初版 `rg runlist` 用 `IntSpan::intersect` 判定（对 19k-span 集合是线性
O(n) 扫描），仅比 rgr 快 ~1.3×。改为复用 `SpanIndex`（有序 span 数组 +
两次二分，O(log n + k)）后：

* overlap：覆盖数 > 0；
* non-overlap：覆盖数 == 0；
* superset：覆盖数 == 区间长度（完全包含）。

三个 op 都从单次 `overlap()` 查询派生，每行 O(log n + k)。rgr 的
superset 走 `diff`（O(n·m) 旧实现），在 19k-span 集合上每行数毫秒，所以
差距最大（~588×）。

## 结论

与 prop 同源：性能差来自"线性交集 vs 二分重叠段"。`rg runlist` 已与
`rg prop` 一样用 SpanIndex，三种过滤 op 均达毫秒级；相对 rgr 快两个
数量级。
