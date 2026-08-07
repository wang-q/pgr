# SIMD / HV / Jaccard 基准（自 hnsm 迁移）

> 2026-08-06 自 `hnsm/benches/{simd,jaccard,hd}.rs` 迁移至 `benches/`，
> 对齐 pgr 现役实现（`libs::linalg` / `libs::hv`，wide 1.6.0）。
> 机器：AMD Ryzen 9 7945HX（32 核，x86_64），stable 1.97.0，release 构建。

## 1. SIMD：L2 范数（`benches/simd_benchmark.rs`）

向量长度 10005（8-lane SIMD 的余数分支随之验证）。四个实现：
`map` / `fold` / pgr `linalg::norm_l2`（wide SIMD）/ `nalgebra::DVector::norm`。

| 实现 | 耗时 |
|---|---:|
| norm_map | 5.617 µs |
| norm_fold | 5.610 µs |
| **norm_simd_pgr** | **723.9 ns** |
| norm_nalgebra | 1.354 µs |

SIMD 比 map/fold 快 **~7.8×**，比 nalgebra 快 ~1.9×。
历史（Ryzen 7 8745HS）：map 6.31 µs / simd 810 ns / nalgebra 1.56 µs——趋势一致。

## 2. HV 编码（`benches/hv_benchmark.rs`，hv_d = 4096）

对比 pgr 现役 `hash_hv_bit`（位操作 SIMD）、`hash_hv_i8`（i8 累加 SIMD）
与三个标量 + RNG 对照（RapidRng / StdRng / SmallRng，i16 累加）。
seed 集合固定种子（StdRng seed 42，RapidHashSet）。

### small（1000 seeds）

| 实现 | 耗时 |
|---|---:|
| **hash_hv_i8** | **447.8 µs** |
| hash_hv_bit | 677.1 µs |
| 标量 + RapidRng | 906.9 µs |
| 标量 + SmallRng | 936.0 µs |
| 标量 + StdRng | 1.039 ms |

### medium（10000 seeds）

| 实现 | 耗时 |
|---|---:|
| **hash_hv_i8** | **4.452 ms** |
| hash_hv_bit | 6.818 ms |
| 标量 + RapidRng | 9.150 ms |
| 标量 + SmallRng | 9.376 ms |
| 标量 + StdRng | 10.445 ms |

i8 实现比 bit 实现快 ~1.53×、比最快标量对照快 ~2×。
历史（2026-01-30）：i8 421 µs / 4.20 ms，lib(bit) 670 µs / 6.73 ms——一致。

### AVX2 快速路径（2026-08-07，Ryzen 9 7945HX = Zen 4；主实现）

`hash_hv_bit` / `hash_hv_i8` 以 **AVX2（256-bit）为主实现**：x86-64 运行时
检测 `avx2`，无则降级到上表 wide 可移植路径；输出与串行逐位一致，由
`test_hash_hv_bit_serial_vs_simd` / `test_hash_hv_i8_serial_vs_simd` 及
avx2 变体测试保证。AVX-512 实现**只保留在 `benches/hv_benchmark.rs` 作
参考对照，不参与运行时分派**（作者决策：AVX-512 相对 256 无优势，实测
Zen 4 两者持平，见下表参考列）。

**平台策略**：AVX2 只是 x86-64 的可选加速；非 x86-64 目标（如 aarch64，
wide 自动映射到 NEON）与不支持 AVX2 的 x86 CPU 走 wide 可移植路径，
语义逐位一致。已用 `cargo check --target aarch64-apple-darwin --lib`
交叉编译验证（2026-08-07）。

| 场景 | bit（AVX2 主路径） | i8（AVX2） | AVX-512 参考 bit | AVX-512 参考 i8 | i8 旧 wide | bit 旧 wide |
|---|---:|---:|---:|---:|---:|---:|
| n=10k, D=4096 | **1.14 ms** | 2.09 ms | 1.48 ms | 2.08 ms | 4.33 ms | 6.69 ms |
| n=100k, D=4096 | **11.4 ms** | 20.9 ms | 14.8 ms | 20.8 ms | 43.4 ms | 66.8 ms |
| n=10k, D=16384 | —（≈参考值） | — | 5.91 ms | 8.30 ms | 17.4 ms | 26.9 ms |

AVX2 主路径实测：bit ±1 编码相对旧 wide bit **~5.9×**、相对旧 wide i8
**~3.8×**；i8 保语义 ~2.1×。AVX-512 参考列（同机测量）与 AVX2 持平，
佐证"256 作为主实现"的决策。速度来自两点：

1. **RapidRng 跳步 + 块主序**：RapidRng 状态是常数步长计数器，输出
   j = mix(seed + j·SECRET0, …)，因此可以把 HV 分块常驻寄存器、遍历全部
   seed，每 seed 每 32 维只做 1 次 RNG（旧 bit 路径的 next_u32 语义，逐位
   一致），且不同 seed 的 mix 相互独立、天然 ILP；
2. **256-bit 位展开 + 延迟 −N 偏移**：`vpsrlvd + vpand + vpslld + vpaddd`
   四个 8-lane 寄存器一次覆盖 32 维（AVX2 intrinsics，无 wide 的冗余
   展开/装载）；因 ±1 数值围绕 0 平衡，每种子只累加 2b，−N 每块末尾一次
   减去（每组少 1 个向量 op）。

**限速步骤分解**（n=10k、D=4096，AVX2 主路径逐组件测量，2026-08-07）：

| 路径 | 完整 | 纯拆分/展开 | 纯 RNG（独立） |
|---|---:|---:|---:|
| bit（主路径，deferred 前） | 1.41 ms | **1.37 ms** | 0.71 ms |
| i8 | 2.10 ms | 1.49 ms | 2.85 ms |

bit 路径中**位拆分/展开是主导**（1.37/1.41 ms，~97%），RNG 只占独立
0.71 ms 且与 SIMD 展开在不同执行端口上**重叠**（标量端口 vs 向量端口）；
采用延迟 −N 后完整耗时降到 **1.14 ms**。i8 路径 RNG 调用是 bit 的 4 倍
（每 8 维一次 vs 每 32 维一次），RNG 独立成本更高，但同样与展开重叠，
总耗时介于两者之间。

**累加器宽度实验（2026-08-07）**：作者提示 bit 路径数值围绕 0 平衡、
i16 应够用。实测两种方案：

* **i16 lane 分段累加**（每段 8192 种子，段内值域确定地落在 [−8192, 8192]，
  无溢出）反而慢 **~2.5×**（3.50 ms vs 1.40 ms）：16-bit 变量移位
  （`vpsrlvw`）延迟更高，且 32 维块从 4 条独立累加链降到 2 条，链式
  依赖主导；
* **延迟 −N 偏移（采用）**：每种子只累加 2b、每块末尾一次性减 N——每组
  少 1 个向量 op，4 条链不变，bit 主路径 **1.40 → 1.14 ms**（n=10k，
  ~1.23×）、n=100k **14.0 → 11.4 ms**。

i8 路径确认需要 i32：字节均值 −0.5（直流偏置，不上下平衡），260 万种子
实测各维 ±1.3e6，i16 必溢出（见 [[dist-cohort-validation.md]] §2）。

**幅度与区分度无关（数值模拟，2026-08-07）**：作者提出 i8 每维幅度大、
同样 D=4096 下基因组间差距更大。模拟（N=10k、D=4096，共享 2000/500/50）
显示：bit ±1 与零均值 i8（±127）的 Jaccard 估计误差**完全相同**
（幅度是比值的公共因子，分子分母同时缩放、抵消）；真正的 i8（均值
−0.5）误差反而最大（0.16–0.19，直流偏置引入 ~0.25·N_A·N_B 二次噪声底）。
区分度由共享种子占比与 D 决定，与编码幅度无关。

### 迁移中发现并修复的性能退化

`hash_hv_i8` 初版迁移用 `bytes.map(|b| b as i8 as i32)`（8 元素标量转换）
构造 SIMD 向量，导致 small 447 → 1789 µs（~4.3× 退化，medium 同样）。
修复：wide 无 u8→i8→i32 数值 lane 转换，改用
`u8x16 → u16x8::from_u8x16_low → i32x8::from_u16x8` 零扩展链 +
`(x << 24) >> 24` 算术移位还原有符号语义（等价 `b as i8 as i32`）。
`hash_hv_bit` 的 u32→i32 转换改 `bytemuck::cast`（0/1 bit pattern 重解释，零开销）。
正确性由 `test_hash_hv_i8_serial_vs_simd` 等对照测试保证。

## 3. 集合 Jaccard（`benches/jaccard_benchmark.rs`）

每个集合 5005 个 u64（0..u16::MAX），105 个集合，随机取两集合比较。

| 实现 | Jaccard 耗时 |
|---|---:|
| **rapidinlinehash** | **20.55 µs** |
| rapidhash | 20.66 µs |
| BTreeSet | 38.24 µs |
| HashSet（SipHash） | 54.46 µs |
| tinyset | 54.41 µs |
| nohash | 59.98 µs |

rapidhash 系最快，比 std HashSet 快 ~2.6×；intersection/union/access
（btree 38.4/47.6 µs、hashset 55.4/40.8 µs、btree_access 4.7 ns）。
历史（msvc）：rapidhash 22.7 µs 最快、hashset 75.6 µs——结论一致。

## 4. 结论

- SIMD 加速真实有效：norm ~7.8×、HV i8 编码 ~2× over 标量。
- pgr 现役选择（`norm_l2` wide SIMD、`hash_hv_i8`）均为基准最优/次优，
  迁移到 wide 后性能与 std::simd 时代持平。
- 集合场景若日后需要，rapidhash 的 HashSet 比 std 默认快 ~2.6×。

## 5. 遗留疑虑（2026-08-07 已解决）

HV 编码矢量化提速不足的根因是 **RNG 生成主导耗时**（n=10k、D=4096：
rng-only 2.68 ms vs 总 4.36 ms），8-lane AVX2 只加速了位展开，无法减少
RNG 调用次数（i8 每 8 维一次 128-bit mix），也消除不了串行 next_u64 的
依赖链。解法见 §2 AVX-512 小节：跳步 RNG + 块主序 + 512-bit 展开，
实测 bit ±1 编码 2.9–4.5×、i8 保语义 2.1×（详见 §2 表格）。
