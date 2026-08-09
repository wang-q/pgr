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

## simd-minimizers 参考分析（含 packed-seq）

SEA 2025 论文 "SimdMinimizers: Computing random minimizers, fast" 的 Rust
实现（源码在项目根 `simd-minimizers-master/`）。

### 它是什么

SIMD 加速的随机 minimizer / syncmer 计算库：

- 序列切成 8 个 chunk 并行处理（`packed_seq::iter_bp`，一次 SIMD 处理
  8 个碱基）；
- 32-bit **ntHash** 滚动哈希（查表 + xor + rol，SIMD 宽度内无依赖）；
- **滑动窗口最小**用分块 prefix/suffix 法（Bender-Farach 式），不用
  单调队列；
- canonical 版本：左/右最小 + TG 偏好；也支持 closed/open syncmers。

依赖：packed-seq 5.0.0（**底层用 `wide`——pgr 已在用**）+ seq-hash
0.2.0 + ensure_simd；主路径要求 AVX2/NEON（`target-cpu=native`），有
`scalar` feature 但慢。

### "人基因组 4 秒"的来源（不是端到端）

`bench/src/bin/paper.rs` 的 `bench_human_genome`：

- `read_human_genome`（读 FASTA + `PackedSeqVec::from_ascii` 2bit 打包）
  在计时**外**；
- `time()` 只包住从**已打包 2bit 输入**开始的 minimizer 计算：人基因组
  3.1 Gbp → 775 MB（数据量 1/4），约 4 s（单线程）。

快的三个 SIMD 部件：`iter_bp` 8 路碱基迭代 + ntHash 8 路并行滚动 +
SIMD 分块滑动窗口最小。

### 与 pgr 的关键差异（为什么借鉴处处碰壁）

| | simd-minimizers | pgr |
|---|---|---|
| 输入 | 2bit 打包内存（不计打包时间） | ASCII FASTA / gz（需解压+打包） |
| 哈希 | 32-bit ntHash（查表+xor，SIMD 友好） | u128 2bit key + rc_key + rapidhash（字节哈希，依赖链） |
| 窗口最小 | SIMD 分块 prefix/suffix | 单调队列（已改分块法） |

pgr 的采样语义锚定在 u128 canonical key + rapidhash：8 路并行滚动已证伪
（见下节）；换哈希则破坏采样集合/兼容性。结论：思路借鉴，依赖不引入。

### packed-seq 评估（参考价值 > 使用价值）

packed-seq 5.0.0：2bit 打包序列内存格式 + `iter_bp`（SIMD 8 路碱基
迭代/转置）+ ASCII↔打包转换。

- **参考**：`iter_bp` 的 8 路 SIMD 打包/迭代技巧可用于 pgr 的 2bit 场景
  （pgi 索引构建、twobit 打包——后者已 SIMD 化），是批量打包实现样板。
- **不直接引入**：① pgr 的 2bit 是 twobit 标准文件格式，packed_seq 是
  内存格式，不兼容；② 新依赖（seq-hash/ensure_simd，违反约束）；
  ③ 覆盖不了 pgr 核心瓶颈（u128 滚动 + rapidhash）。
