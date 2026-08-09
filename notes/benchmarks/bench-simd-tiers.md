# SIMD 三级回退速度对比（AVX2 / wide128 / 标量）

> 2026-08-09 汇总。全部为 criterion 实测（AMD Ryzen 9 7945HX，release），
> 数据来自 `target/criterion/` 存档（`byte_stat_*`、`poa_*`、
> `norm_*`、`twobit_from_dna`）。wide 一律为 **128-bit 类型**
> （SSE2/NEON 原生）；标量为同一函数的逐元素参考。

## 对比表

| 函数 | 规模 | 标量 | wide(128) | AVX2 | wide/标量 | AVX2/标量 | AVX2/wide |
|---|---|---:|---:|---:|---:|---:|---:|
| `count_valid` | 1 MB | 513 µs | 182 µs | 35.7 µs | 2.8× | 14.4× | 5.1× |
| `count_bases` | 1 MB | 4.66 ms | 624 µs | 98.7 µs | 7.5× | 47× | 6.3× |
| `count_n` | 10 MB | 2.81 ms | 回退标量* | 430 µs | — | 6.5× | — |
| `masked_bitmap` | 10 MB | 6.40 ms | 回退标量* | 425 µs | — | 15.1× | — |
| POA 对齐 120 bp | — | 904 µs | 150 µs | 111 µs | 6.0× | 8.1× | 1.35× |
| POA 对齐 600 bp | — | 23.4 ms | 3.44 ms | 1.94 ms | 6.8× | 12.1× | 1.8× |
| `norm_l2` | 10005 | 5.86 µs | 747 ns | 无手写** | 7.8× | — | — |
| `twobit` 分类 | 10 MB | 53.5 ms | 6.75 ms | 5.43 ms | 7.9× | 9.9× | 1.24× |
| `cigar` 分类 | 40 M 列 | 138 ms（旧双趟） | 无 wide*** | 9.6 ms | — | ~14× | — |

* `count_n`/`masked_bitmap`：12 值判等在 wide 上实测比标量慢 42%，回退标量。
** `norm_l2`：无 AVX2 手写路径，wide 双累加器为主实现。
*** `cigar` 分类：只有 AVX2 + 标量，未实现 wide。

## 解读

- **纯逐字节统计类**（count_valid/count_bases）：AVX2 比 wide 快 5–6×——
  256-bit 一次 movemask 处理 32 列，128-bit 只有一半吞吐。
- **DP/打包类**（POA、twobit）：AVX2 只比 wide 快 1.2–1.8×——128-bit 已
  接近依赖链/打包型任务的硬件瓶颈，加宽收益递减。
- **wide 的定位**：主流非 AVX2 平台（SSE2/NEON）拿到 6–8× 级收益；AVX2
  是上限；标量保证任何 CPU 可跑（三级原则见
  `design/simd-optimization.md` §5 第 7 条）。
- 完整 `from_dna`（含块合并+打包）端到端 39.6 ms，说明串行部分仍是
  剩余成本（分类本身 AVX2 仅 5.4 ms）。
