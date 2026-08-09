# 热点 profiling 实测（2026-08-09）

> 目的：验证"下一步优化空间"的假设（gzip I/O 主导？纯计算热点在哪？），
> 工具 perf（`cpu/cycles/P`，997 Hz）+ flamegraph，release 构建；
> 机器 AMD Ryzen 9 7945HX（Zen 4，AVX2/AVX-512）。火焰图产物在
> `/tmp/{size,kmer,pgibuild}-flame.svg`（临时，未入库）。

## 场景 1：`fa size --no-ns`（50 MB gz，/tmp/big.fa.gz）

| 符号 | self |
|---|---:|
| `zlib_rs::inflate::inflate_fast_help_avx2` | 42.8% |
| `__memset_avx512_unaligned_erms` | 35.0% |

- **zlib-rs 的 inflate 内部已是 AVX2**，逐字节统计（nt_simd）未进前列；
  另一次计时中纯解压占 gz 总耗时的 ~81%。
- 结论：gz 输入下 I/O（解压 + 输出缓冲清零）主导，继续扣统计 SIMD
  边际收益很小；inflate 的"SIMD 化"已被依赖库吃掉。

## 场景 2：`rept s-kmer`（mg1655.fa.gz，**单 contig**）

| 符号 | self |
|---|---:|
| `table_profiles`（嵌套 Vec collect） | 78.5% |
| `radix_sort::partition_at` | 10.4% |

- `table_profiles` = 每个 k-mer 窗口一次 `partition_point` 二分查表
  （`src/libs/kmer/profile.rs`）：~460 万窗口 × ~23 次 u128 比较，全是对
  ~73 MB 排序表（远超 L3）的随机访问。perf stat：**cache-misses 1.31 亿 /
  41% miss rate**，即每窗口 ~28 次 cache miss，与二分比较次数吻合。
- 单 contig → "按序列 par_iter" 并行无效；热点是**内存访问模式（数据
  结构）而非 SIMD**。优化候选：哈希表（O(1) 期望，miss 降一个数量级）
  或按 key 高位分桶 + 桶内小二分（首级 O(1) 索引、桶 ~70 条可入 L1）。
- 源码对照（2026-08-09）：pgr 的二分查表是迁移简化，FastK 原版
  （`FASTK-master/`）用前缀索引 + 分桶（`libfastk.c` `_Kmer_Stream` 的
  `index[]`、`split.c` `Split_Table` 的 1-byte 前缀块表），查询 O(1) 定位
  + 桶内小范围；详见 [[../design/kmer.md]] §3.5。
- radix_sort 已是并行实现，占比低。

### 优化落地（2026-08-09，排序合并）

按 FastK 结构先试前缀索引，**隔离基准证伪**（`kmer_lookup`，MG1655
4.6 M 窗口）：全局 `partition_point` 1.125 s vs 前缀桶 `binary_search`
1.195 s——73 MB 表随机访问延迟主导，比较次数 23→7 无收益（前缀桶的
首次访问仍是随机 DRAM/L3 miss）。结论：**逐窗口查表模式本身是瓶颈，
换查找结构没用**。

最终改为 FastK 的排序合并路线（`src/libs/kmer/profile.rs`）：收集全部
窗口 key（rayon 并行）→ `radix_sort_u128_par` → 与排序去重的 `table.keys`
线性归并一次写回。criterion 对比（同机同参数）：

| 基准 | 旧（逐窗口二分） | 新（排序合并） | 提升 |
|---|---:|---:|---:|
| self_profiles_mg1655 | 1.43 s | 250 ms | ~5.2× |
| relative_profiles_mg1655 | 1.41 s | 270 ms | ~5.4× |
| `rept s-kmer` 整命令 | 1.67 s | 0.50 s | ~3.4× |

语义由 `sort_merge_matches_binary_search` 对照测试 + 全量集成测试保证。

## 场景 3：`pgi build`（mg1655.fa.gz）

| 符号 | self |
|---|---:|
| `pgi::execute`（memmove 写索引 13.9%） | 24.1% |
| `build_from_seqs`（内部 15.8% + hash_one 2.6%） | 22.4% |

- 无单一主导热点：syncmer 采样 + HashMap 去重 + 排序 + 写文件均摊；
  总耗时 ~0.85 s。无明显向量化切入点。

### 优化落地 2：MAF→PAF 的 CIGAR 生成（2026-08-09）

`maf_block_to_paf`（`src/libs/paf/maf_import.rs`）对每个 block 调用
`cigar_from_alignment` + `cs_from_alignment`，40 M 列 `maf to-paf` 实测
占 53.4%（两次逐列扫描 + cg:Z 格式化 11.8%）。优化：`classify_alignment`
（AVX2 + 标量回退）一次生成 I/D/=/X 四掩码，两个生成函数共享，并用
`trailing_zeros/ones` 位运算跳扫 match run。0.55 s → 0.347 s（~37%），
输出逐字节一致；剩余瓶颈 cg:Z 格式化 ~13.6%。

函数级拆分（`benches/cigar_benchmark.rs`，40 M 列，95% match）：

| 测量 | 耗时 | 说明 |
|---|---:|---|
| old_two_pass | 138 ms | 旧实现（两次逐列扫描 + 输出） |
| new_classify_scan | 88 ms | 新全流程，~1.57× |
| new_classify_only | 9.6 ms | 纯 SIMD 分类，~7× 级 |

口径教训：纯计算函数（count_bases 等）微基准直接给 10–50×；cigar 的
另一半是字符串输出（~79 ms，不可 SIMD），函数级仅 1.57×、端到端 37%。
百分比 vs 倍数 = 优化覆盖面差异，不是 SIMD 失效。

## 结论

1. gzip 主导假设成立，且 inflate 已被依赖库（zlib-rs AVX2）优化。
2. **下一个真实热点是 `table_profiles` 的 cache miss（数据结构/访问
   模式），不是 SIMD 能解决的**——与"向量化三步模式"的边界一致。
3. `rept s-kmer` 场景从 ~1.7 s 里可挖的 ~1.3 s 在查表；哈希/分桶化后
   预计 5–20×，比再扣一个逐字节 SIMD 更值。
