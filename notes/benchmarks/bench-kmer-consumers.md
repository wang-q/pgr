# Kmer 采样消费方分析（2026-08-09）

> 盘点所有使用 k-mer 采样的命令，量化"窗口最小型"与"哈希草图型"各自
> 的占比，确定 simd-minimizers 思路（分块窗口最小 / 8 路并行哈希）的
> 推广价值。perf 实测（Ryzen 9 7945HX，release，mg1655 单基因组场景）。

## 消费方汇总

| 消费方 | 采样类型 | 采样占比 | 端到端 | 借鉴价值 |
|---|---|---:|---:|---|
| `pgi build` | syncmer（窗口最小） | ~7.5% | 850 ms | ✅ 分块法已落地（64→53 ms，-16.6%）；SIMD 上限 ~5% |
| `align rest` prefilter | syncmer/minimizer（窗口最小） | <5% | 22 s | ❌ LASTZ 主导 |
| `dist mini` | minimizer（窗口最小） | 8% | 90 ms | ❌ 小场景 + I/O 主导 |
| `dist frac` | FracMinHash（**哈希草图**） | ~43% | 133 ms | ⚠️ 占比最高；需 8 路并行哈希，非窗口最小 |
| `dist mash` | Mash（哈希草图） | 分散 | 107 ms | ❌ memcmp 主导 |
| `dist hv` / `pgi to-hv` | syncmer/minimizer + HV 编码 | 小 | — | ❌ RNG 主导 |

## 结论

- **窗口最小型**（分块法可推广）：pgi build 是占比最高处，已优化；
  align rest / dist mini 采样占比 <8%，被 LASTZ / I/O 主导，不值得
  继续 SIMD。
- **哈希草图型**（frac/mash）：`dist frac` 采样 43% 是最高单点，但性质
  是逐 k-mer 哈希 + 阈值筛选——要加速需借鉴 simd-minimizers 的
  **8 路并行滚动哈希**（ntHash 思路的 rapidhash 变体），批量 cohort
  两两 dist 场景有潜力。
- **hv**：RNG 主导，维持原结论。

## 基线基准

- `benches/syncmer_benchmark.rs`：syncmer_dna 基线（64 ms，分块法后 53 ms）
  + rolling_hashes_only（~6 ms，占 9%——热点在窗口最小，非哈希）。
- `benches/fracminhash_benchmark.rs`：`seq_fracminhash` 基线（mg1655，
  k21）：scale=1000 89.5 ms、scale=100 92.2 ms（scale 影响小，成本在
  滚动 key + rapidhash）——8 路并行哈希试点的 A/B 基准。
- `benches/minimizer_benchmark.rs`：`seq_mins` 基线（mg1655，k21/w5）：
  rapid 46.6 ms、fx 48.1 ms（窗口最小型，若分块法推广到 minimizer 时
  A/B 用）。

## frac 8 路并行试点（2026-08-09，证伪）

细分：fracminhash 87 ms 中 canonical 打包（rolling + rc_key）占 84 ms
（96.6%），rapidhash 仅 ~3 ms。据此实现 8 路并行滚动（8 个相邻 k-mer
同步推进，u128 独立依赖链 → ILP），**实测 87 → 282 ms（慢 3.2×）**：
每 key 全程滚动 = 8 倍滚动工作量，ILP 抵消不了冗余；且 u128 移位 +
mask 未充分并行。已回退。教训：simd-minimizers 的 8 路收益来自
**SIMD 打包 + 并行哈希**（packed_seq），不是"8 个 u128 标量滚动"；
pgr 的 u128 key + rapidhash 结构不匹配，此方向关闭。
