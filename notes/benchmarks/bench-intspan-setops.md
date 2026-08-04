# IntSpan 集合运算与批量构建基准（criterion）

> 目的：把 `intersect`/`union`/`diff`/`xor` 线性化与 `from_pairs` 批量构建
> 的收益固化为库级 criterion 基准，直接对比新旧实现（旧实现用公开 API
> 在 bench 内重建）。2026-08-04 实测。

## 基准文件

`benches/intspan_setops_benchmark.rs`：

* `setops`：5k / 20k span 的两条 runlist（100 Mb 染色体，长度 100–2000），
  四个集合运算 × 新旧实现。
* `construction`：10k / 100k 个随机 1 bp span（互不重叠的 adversarial
  输入），`from_pairs`（排序构建）vs `add_pair` 循环（逐个插入）。

```bash
cargo bench --bench intspan_setops_benchmark
```

## 结果（median，release，measurement-time 2s）

### setops

| op | n=5k 旧 | n=5k 新 | 加速 | n=20k 旧 | n=20k 新 | 加速 |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| intersect | 672 µs | 29.7 µs | **23×** | 11.45 ms | 102.8 µs | **111×** |
| union | 1.13 ms | 26.1 µs | **43×** | 15.2 ms | 87.9 µs | **173×** |
| diff | 1.51 ms | 20.9 µs | **72×** | 16.9 ms | 76.2 µs | **221×** |
| xor | 2.00 ms | 84.0 µs | **24×** | 32.8 ms | 445.8 µs | **74×** |

加速随 n 放大（O(n·m) → O(n+m)）；20k 时 74–221×，远高于 CLI 基准的
5–8×（CLI 时间被 runlist JSON 加载占据）。

### construction（随机 1 bp 稀疏区间，adversarial）

| n | `from_pairs` | `add_pair` 循环 | 加速 |
| :--- | ---: | ---: | ---: |
| 10k | 124 µs | 2.24 ms | **18×** |
| 100k | 1.66 ms | 173 ms | **105×** |

`add_pair` 逐个插入无序区间是 O(n²)（VecDeque 中间搬移），`from_pairs`
排序 + 单遍合并为 O(n log n)，差距随 n 放大。

## 结论

库级基准证实：集合运算线性化（20k span 时 ~100–220×）与批量构建
（100k 时 ~105×）都是数量级收益，且随规模增长；这与 CLI 基准
（bench-rg-prop.md、interval-overlap.md）互补——CLI 看端到端，criterion
看纯库层。旧的 O(n·m)/O(n²) 实现保留在 bench 内作基线，便于后续回归
监控。
