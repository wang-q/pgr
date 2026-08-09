# fa 逐字节统计 SIMD 基准（nt_simd）

> 2026-08-09。`libs/nt_simd.rs`（count_valid / count_n / masked_bitmap）矢量化
> 的分派：单一 AVX2 手写路径 + `is_x86_feature_detected!` + portable `wide`
> 回退（复杂判定如 N 家族 12 值在非 AVX2 平台回退标量，见代码注释）。
> 机器：AMD Ryzen 9 7945HX（x86_64，AVX2），stable 1.97，release。
> 输入：随机 DNA（大写 ACGT + 5% 小写/N/IUPAC/无效字节），长度 1 MB / 10 MB。

## count_valid（`fa size --no-ns`：A/C/G/T/U 大小写计数）

| 规模 | scalar（NT_VAL 查表） | wide | avx2 | avx2/scalar |
|---|---:|---:|---:|---:|
| 1 MB | 513 µs | 182 µs | 35.7 µs | ~14.4× |
| 10 MB | 5.11 ms | 1.81 ms | 358 µs | ~14.3× |

wide 相对标量 ~2.8×（128-bit 拆分 + 判等掩码开销）。

## count_n（N 家族：IUPAC 歧义码 + N + X，12 个小写化值）

| 规模 | scalar（NT_VAL 查表） | wide | avx2 | avx2/scalar |
|---|---:|---:|---:|---:|
| 10 MB | 2.81 ms | 3.99 ms（慢 42%，回退标量） | 430 µs | ~6.5× |

`wide` 12 次判等开销超过标量 filter（2026-08-09 实测），非 AVX2 平台回退
标量——**不可类推**：同是 12 值判等，`count_bases` 的 wide 却有 ~7.5×，
因为它的标量基准（12 值链式 match）慢得多。

## masked_bitmap（`fa masked` 掩码位图，默认模式 = 小写 ∪ N 家族）

| 规模 | scalar | avx2 | avx2/scalar |
|---|---:|---:|---:|
| 10 MB | 6.41 ms | 425 µs | ~15.1× |

## count_bases（`fa count`：A/C/G/T/N 五类计数，U→T、IUPAC/X→N）

| 规模 | scalar | wide | avx2 | wide/scalar | avx2/scalar |
|---|---:|---:|---:|---:|---:|
| 1 MB | 4.66 ms | 624 µs | 98.7 µs | ~7.5× | ~47× |
| 10 MB | 46.9 ms | 6.26 ms | 987 µs | ~7.5× | ~47.5× |

标量版为 12 值 N 家族 match 分支（每字节最多十几次比较），AVX2 固定
17 次向量判等 + popcount，故加速比高于 count_n；`wide` 同样有 ~7.4×
（首版按 count_n 推断"wide 无收益"走标量，实测后纠正为 wide 路径）。
实现 2026-08-09，`nt_simd::count_bases` + `fasta::stat::count_bases` 转发。

## 结论

- AVX2 对纯字节扫描/比较型统计 ~6–15×，与 POA 的 DP 向量化量级一致；
  `wide` 回退对 ≤5 值判等仍 ~2.8×，复杂判定（12 值）放弃向量化。
- 单基因组（5 MB）CLI 上 `fa size --no-ns` / `fa masked` 均为 ~10 ms 级
  （I/O 主导），SIMD 收益在 4 万 cohort 批量处理时累积（每基因组 ~4 ms →
  ~0.4 ms 量级）。
