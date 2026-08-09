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

- gzip 并行解压：**已裁定不做**（2026-08-09——程序常被 shell 包裹并行
  执行，pgr 侧 `fa` 保持单线程；zlib-rs inflate 内部已是 AVX2）。
- 单文件多 contig 的 rayon 并行（`dist mash` 序列级并行等待同类场景）。

### 已评估不做

- `dist mash` murmur3 SIMD（1.3×，风险大于收益）。
- HV AVX-512 跳步 RNG（只保留在 `benches/hv_benchmark.rs` 作对照）。

## 5. SIMD 优化设计原则（步骤清单，2026-08-09 定稿）

### 阶段 0 · 立项前

1. **约束检查**：fa 路径保持单线程（不引入 rayon）、不引入新依赖、
   只动该动的地方。
2. **消费方核实**：全库搜索所有生产调用方与量级（搜索勿用 `head`
   截断——`cigar_from_alignment` 曾因此误判"无消费方"）；无消费方
   的函数不做。
3. **场景 profile**：对真实场景 perf 确认该函数在 CPU 分布中的占比；
   先排除 I/O 主导（gz 解压、memset）、RNG/哈希主导、依赖链（滚动
   打包、单调队列）——SIMD 无效或覆盖有限。

### 阶段 1 · 热点与语义

4. **隔离验证**：微基准单独测目标函数确认是热点（`canonical_keys`
   实测 2.95 ms 后直接证伪，省下实现成本）。
5. **先基准后修改（2026-08-09 twobit 教训）**：改动前先用 criterion 为
   现状实现建立基线（可复现、存档），确认后再开始实现；实现后在同一
   基准上 A/B。**禁止先改后补基准**——twobit `from_dna` 直接实现、事后
   才补 `twobit_benchmark.rs`，函数级加速比只能事后补测，缺干净基线。
6. **语义锚定**：写下与标量逐位一致的契约（U→T、X→N、
   `eq_ignore_ascii_case` 精确行为），作为测试基准。

### 阶段 2 · 实现

7. **数据结构先行**：先检查访问模式再谈向量化——FastK 对照显示
   "排序合并替代逐窗口查表"（~5×）远大于 SIMD 收益；查表型热点先
   考虑数据结构。
8. **三级回退链（2026-08-09 用户定稿）**：
   1. **AVX2 手写**（`is_x86_feature_detected!` 运行时检测）第一优先；
   2. 无 AVX2 → **`wide` 固定 128-bit 类型**（`u8x16`/`i32x4`/`f32x4`）：
      128-bit 在 x86_64 是 SSE2、aarch64 是 NEON，均为平台原生宽度，
      且 **wide 的 128-bit 分支先匹配 `sse2`，不受编译期 avx2 影响**
      （反汇编实证：`u8x16` 在 `+avx2` 编译下仍 0 条 ymm；`u8x32`
      则有 3 条）——老 CPU 因此可兜底。**nt_simd 已改造**
      （2026-08-09，`u8x16` + `SimdPath::Scalar` 显式兜底，性能与 256-bit
      时代一致）；**linalg/poa/hv 已改完（2026-08-09 晚）**。教训：有跨
      chunk 累加依赖链的函数（linalg norm/dot 等），128-bit 化必须配
      **双累加器**（8 元素块拆两个独立向量），否则依赖链变长慢 ~2×
      （实测 norm 0.72→1.45µs，双累加恢复 747ns）；逐 chunk 独立无累加链
      的函数（nt_simd 统计）无此问题。
   3. 无法检测/无 SIMD 平台 → **纯标量**最终兜底。
   无 SSE4.1 中间档、无 SIMDe；wide 256-bit 类型（`u8x32` 等）**禁用**。
9. **wide 必须实测且受编译前提约束**：
   * 宽/标量收益不能按函数类推：`count_n` wide 慢 42%（回退标量）、
     `count_bases` wide 7.5×（保留）——同模式结果相反，必须实测。
   * **历史教训（已由 128-bit 原则取代）**：wide 的 256-bit 类型
     （`u8x32`）由编译时 `cfg(target_feature="avx2")` 决定，全局开 avx2
     编译时变真 AVX2，无 AVX2 CPU 上 SIGILL（反汇编实证）。改用 128-bit
     类型后此问题消失；`nt_simd` 曾用 `#[cfg]` 双轨门控 + 纯标量兜底，
     128-bit 化后可简化（仍保留纯标量作为第三级）。
   * **NEON 维度（2026-08-09 验证）**：wide 在 `aarch64 + neon` 下
     `u8x16 = uint8x16_t`（源码确证），`cargo check --target
     aarch64-apple-darwin` 通过（x86 特有代码均有 cfg 门控）——**编译与
     NEON 指令可用可保证**；但① 无手写 NEON 路径，受益上限 128-bit
     级（~2–4× 预期，非 AVX2 的 10–47×）；② 无 aarch64 硬件/CI，wide
     在 NEON 上是否正收益未实测；③ 32 位 ARM 不在 wide 支持范围
     （编译失败）。若 NEON 是发布目标，需加 aarch64 CI + 实测基准。
   * **必须保留纯标量兜底（2026-08-09 加固）**：每个 SIMD 函数都要有
     不依赖任何 SIMD 的最终回退。nt_simd 的 `count_valid`/`count_bases`
     曾删掉标量、wide 是唯一非 AVX2 路径——若编译开 avx2（wide 变真
     AVX2）则无 AVX2 CPU 上 SIGILL。已加固：标量函数
     `#[cfg(target_feature = "avx2")]`、wide 函数及 `u8x32` import
     `#[cfg(not(target_feature = "avx2"))]`，两种编译 + aarch64 均验证。
     **linalg / hv / poa::simd 仍只有 wide 回退，待加固。**
10. **量化输出占比**：输出型函数（字符串构建、格式化）先量化——SIMD
   只覆盖分类/计算部分；输出主导时预期给"百分比"而非"倍数"
   （cigar：SIMD 分类 7×、函数级 1.57×、端到端 37%）。

### 阶段 3 · 验证

11. **位一致对照**：随机（含边界字符）+ 真实数据，与标量/旧实现逐位
    对照；`cargo test` 的 debug 模式必须过（能抓到 release 不炸的
    shift overflow）。
12. **双口径基准**：函数级 micro-bench（倍数）+ 端到端命令（百分比），
    两者都记录并说明口径。
13. **工程门禁**：fmt / clippy `-D warnings` / 全量测试 clean。

### 阶段 4 · 落盘

14. **结果记录**：基准数据、口径、实现位置写入 `notes/benchmarks/`；
    决策与方法论更新进本文件。
15. **不做也记录**：证伪/暂缓候选写清原因与触发条件（`rev_comp`、
    cigar 无消费方、murmur3、gzip），避免重复立项。

## 6. 全库 SIMD 候选清单（2026-08-09 扫描）

按三步模式适用条件（纯数据并行、无依赖链、非 I/O/RNG 主导）全库扫描，
`fa` 路径保持单线程（SIMD 为指令级并行，不引入 rayon，与用户 2026-08-09
裁定一致）。

### 第一梯队：与 nt_simd 完全同构，低风险，~10× 预期

1. `fasta::stat::count_bases`（`pgr fa count`，`src/libs/fasta/stat.rs`）：
   逐字节 5 类统计（A/C/G/T/N，IUPAC→N，U→T，X→N），= count_valid +
   count_n 的组合，复用 `eq_any_lower` 骨架。**已实现（2026-08-09）**：
   `nt_simd::count_bases`（AVX2 + `wide` 回退，实测 wide ~7.5×、AVX2
   ~47×——首版按 count_n 推断"wide 无收益"走标量，实测后纠正），基准见
   [[../benchmarks/bench-nt-simd.md]]），`fa count` 接入。
2. `nt::complement` / `rev_comp`（`src/libs/nt.rs`）：NT_COMP 查表，
   pshufb 或 eq_any；消费方 pbit/pgi/paf graph/`fa rc`/chain `to_axt`/loc，
   **已评估暂缓（2026-08-09）**：perf 实测 `fa rc` 50 MB 时 rev_comp
   （内联进 `rc::execute`）仅 ~14%，memset 40% + 读 5% 主导；pbit 为短段、
   `to_axt` 低频中等规模——当前消费方都不值得 SIMD 化。

### 第二梯队：模式需变体（位图 + 串行 run 扫描），已确认值得

3. `paf::cigar::cigar_from_alignment` / `cs_from_alignment`
   （`src/libs/paf/cigar.rs`）：逐列 `=`/`X`/`I`/`D` 分类 + run 合并，
   同 `fa masked` 的"SIMD 位图 + 串行扫描"结构。**已确认是热点
   （2026-08-09）**：`maf_block_to_paf`（`src/libs/paf/maf_import.rs`）对
   每个 block 调用两者，`pgr maf to-paf` 40 M 列实测占 53.4%
   （两次逐列扫描 + 格式化 11.8%）；MAF→PAF 是 4 万 cohort 管线核心。
   **已实现（2026-08-09）**：`classify_alignment`（AVX2 + 标量回退）一次
   生成 I/D/=/X 四掩码，`maf_block_to_paf` 共享后分别 `scan_cigar_ops` /
   `scan_cs`；两扫描用 `trailing_zeros/ones` 位运算直接跳 match run。
   40 M 列 `maf to-paf` 0.55 s → 0.347 s（~37%），输出逐字节一致；
   剩余瓶颈为 cg:Z 字符串格式化（~13.6%，非 SIMD 范畴）。

### 第三梯队：曾判 I/O 主导，profiling 翻案后已完成（2026-08-09 晚）

4. `twobit::from_dna` 位打包（`src/libs/fmt/twobit.rs`）：4 碱基→1 字节
   SIMD + 位图追踪 N/小写块。**新证据：`pbit create`（mg1655 自压缩）
   perf 实测 `Blocks::from_dna` 占 40.6%**（其次 LZ 索引 18.5%、inflate
   15.3%）——pbit 编码路径的核心热点，不是 I/O 主导。**已实现
   （2026-08-09）**：`classify_dna` 三级（AVX2 / 128-bit `wide` / 标量测试
   参考）生成 N/小写位图 + 2bit 码向量，`blocks_from_mask` 位运算合并块，
   `pack_codes` 标量打包；`pbit create` 83 ms → 58 ms（~30%），from_dna
   占比 40.6% → 不在前列（LZ 索引成为新主导 40%）。三级 + avx2/aarch64
   配置验证通过；pbit 输出含 HashMap 随机种子导致非确定（既有行为，未改）。

### 明确排除（有据可依）

- `pgi/wave.rs`（Myers wave front）：diagonal 推进依赖链；
- `pbit/lz_diff.rs`：哈希匹配主导；
- `pgi/build.rs` syncmer：滚动哈希 + 单调队列依赖链；
- kmer 滚动打包：已证伪（2.95 ms，占比 ~1%，见
  [[../benchmarks/bench-profile-hotspots.md]]）；
- `alignment/variation.rs`：规模小；
- gzip 解压：用户已裁定不做。

### 实施顺序建议

先 `count_bases`（骨架现成）→ perf 验证 `fa rc` / PAF 生成里 rev_comp 与
cigar 的真实占比，决定第二梯队 → 其余按数据说话。
