# SIMD 优化经验笔记（向量化方法论）

> 2026-08-09 建立。总结 HV / POA / fa 逐字节统计三轮向量化优化的共同
> 模式、适用边界与教训，作为后续热点优化的检查清单。基准数据见
> [[../benchmarks/bench-simd-hv-jaccard.md]]、
> [[../benchmarks/bench-nt-simd.md]] 与 `benches/poa_benchmark.rs`。

## 1. 总体方式（三步模式）

HV → POA → nt_simd 三轮完全同构，固定为三步：

1. **热点选择**：数据量大、逐元素纯计算、无跨元素依赖（profile /
   场景驱动，不拍脑袋）。
2. **双路径实现**：AVX2 手写 intrinsic 为主路径（x86-64），`wide`
   portable 回退，`is_x86_feature_detected!` 运行时分派；无 SSE4.1
   中间档、无 SIMDe（"HV 式"，见 [[hv.md]]）。
3. **验证**：与标量逐位一致 + 随机数据对照测试 + criterion 基准落盘。

## 2. 适用边界（经验教训）

| 热点类型 | 效果 | 实例 |
|---|---|---|
| 纯字节扫描 / 无依赖链 DP / f32 向量运算 | 6–15× | nt_simd 6.5–15×、POA 8.7–12.3×、norm_l2 ~7.8× |
| RNG 生成主导 | 仅 ~2×；AVX-512 跳步 RNG 才 2.9–4.5× | HV 编码（见 [[../benchmarks/bench-simd-hv-jaccard.md]] §5） |
| 哈希函数本身 | 1.3×，不值得 | `dist mash` murmur3（已裁定不做） |
| 依赖链（滚动 2-bit 打包、单调队列） | 不适用，需变体 | kmer `canonical_keys`、syncmer |
| I/O 主导（gzip 解压等） | 向量化被淹没 | 见 §3 实测 |

结论：**先确认热点不是 I/O / RNG / 依赖链主导，再套三步模式。** 合成
随机数据基准能证明"峰值"提速，但真实场景的收益取决于输入形态与瓶颈
分布——基准结论与场景结论要分开记录。

## 3. 实测：I/O 可能主导（2026-08-09）

50.6 MB 随机 DNA（大写 ACGT + 小写 + N），`gzip -9` 后 23.3 MB，
AMD Ryzen 9 7945HX、stable 1.97、release：

| 场景 | 耗时 |
|---|---:|
| `fa size --no-ns` 纯文本 | 87 ms |
| `fa size --no-ns` gzip 输入 | 207 ms |
| 纯 `gzip -dc` 解压 | 168 ms（占 81%） |

项目大量 `.gz` 输入下，单线程 inflate 是主导瓶颈；继续在逐字节统计上
扣 SIMD 的边际收益很小。**下一步先 profile 真实场景（4 万 cohort /
mg1655 全流程）确认 CPU 分布，再选点。**

2026-08-09 已用 perf 实测三个场景（`fa size` gz / `rept s-kmer` /
`pgi build`），数据与火焰图见
[[../benchmarks/bench-profile-hotspots.md]]：gzip 假设成立且 zlib-rs
inflate 内部已是 AVX2；`rept s-kmer` 的 79% 在 `table_profiles` 的
`partition_point` 查表（cache miss 41%），**真实下一个热点是数据结构
（访问模式），不是 SIMD**。

## 4. 后续候选热点

### 第一梯队：可完全套用三步模式，风险低

- `nt::complement` / `rev_comp`（`src/libs/nt.rs`）：逐字节 NT_COMP
  查表，与 nt_simd 同构（pshufb 或 eq_any）；消费方遍布 pbit 压缩/解压、
  pgi build、paf graph builder、`fa rc`、kmer profile。注意 pbit 内为短段
  调用，收益需按场景评估。
- `translate::translate`（`src/libs/translate.rs`）：3 碱基 → 1 aa 查表，
  数据并行但输出紧凑，内存带宽主导，预期 2–4×。
- `twobit::from_dna` 位打包：16 字节 → 4 字节并行；2bit 命令 I/O 主导，
  收益存疑。

### 第二梯队：依赖链需变体，收益大但复杂度高

- kmer `canonical_keys` / `rolling_kmer_keys`（`src/libs/kmer/mod.rs`）：
  2-bit 滚动打包为纯依赖链，可"分块并行打包 + 边界修正"；`pgi build` /
  `rept s-kmer` 核心，4 万 cohort 场景下值得。
- `rc_key` 已有 4-base 查表，可继续展开到 8-base。

### 第三梯队：非向量化，但可能才是真瓶颈

- gzip 并行解压：多成员 gzip 流检测边界后并行 inflate，或换
  zlib-ng / libdeflate。**受"不引入新依赖"约束，需用户决策。**
- 单文件多 contig 的 rayon 并行（`dist mash` 序列级并行等待同类场景）。

### 已评估不做

- `dist mash` murmur3 SIMD（1.3×，风险大于收益）。
- HV AVX-512 跳步 RNG（只保留在 `benches/hv_benchmark.rs` 作对照）。

## 5. 建议流程

1. 端到端 profile 真实场景（perf / 火焰图），确认 CPU 时间分布。
2. 若 `.gz` 主导 → 压缩层解压优化优先；若解压后 CPU 主导 →
   `rev_comp` / `complement` 是套用现有模式最自然、风险最低的下一块。
3. 每个候选按三步模式走：热点确认 → 双路径实现 → 位一致 + 基准落盘。
