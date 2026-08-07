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

> 下表为 **wide 可移植路径**（x86 无 AVX2 时的降级路径；2026-08-06
> 迁移后、AVX2 落地前实测）。2026-08-08 同机复测：AVX2 主路径自动分派下
> small 为 bit 113.1 µs / i8 209.5 µs（medium 见下方 AVX2 表）。

### small（1000 seeds，wide 可移植路径）

| 实现 | 耗时 |
|---|---:|
| **hash_hv_i8** | **447.8 µs** |
| hash_hv_bit | 677.1 µs |
| 标量 + RapidRng | 906.9 µs |
| 标量 + SmallRng | 936.0 µs |
| 标量 + StdRng | 1.039 ms |

### medium（10000 seeds，wide 可移植路径）

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
| n=10k, D=4096 | **1.11 ms** | 2.11 ms | 1.44 ms | 2.14 ms | 4.33 ms | 6.69 ms |
| n=100k, D=4096 | **11.2 ms** | 21.0 ms | 14.4 ms | 21.2 ms | 43.4 ms | 66.8 ms |
| n=10k, D=16384 | 4.49 ms | 8.45 ms | 5.70 ms | 8.62 ms | 17.4 ms | 26.9 ms |

> "旧 wide"列 2026-08-08 已由显式基准（`hash_hv_bit_wide` /
> `hash_hv_i8_wide`）直接复测：bit 6.64 ms、i8 4.38 ms（与表列一致）。

> 2026-08-08 升级：`hash_hv_bit` 主实现改为**每 64 维消耗完整 64 位**
> （低 32 位 → 前 32 维、高 32 位 → 后 32 维），RNG 调用减半，bit 1.15 →
> **1.11 ms**（n=100k 11.5 → 11.2 ms）；AVX-512 参考同步。三处实现
> （AVX2 / wide 回退 / serial 测试参考）逐位一致，测试通过；高低 32 位
> 随机性统计检查通过（各 ~50% 位密度、lo&hi 交叉 ~25%）。

AVX2 主路径实测：bit ±1 编码相对旧 wide bit **~6.0×**、相对旧 wide i8
**~3.9×**；i8 保语义 ~2.1×。AVX-512 参考列（同机测量）与 AVX2 持平，
佐证"256 作为主实现"的决策。速度来自两点：

> 2026-08-08 复测（同机，Ryzen 9 7945HX / stable 1.97）：64 位用满升级后
> 上表为现役数字；D=16384 的 AVX2 实测 bit 4.49 ms / i8 8.45 ms，仍快于
> 同维 AVX-512 参考（5.70 / 8.62 ms）。

1. **RapidRng 跳步 + 块主序**：RapidRng 状态是常数步长计数器，输出
   j = mix(seed + j·SECRET0, …)，因此可以把 HV 分块常驻寄存器、遍历全部
   seed，每 seed 每 64 维只做 1 次 RNG（一次 u64 输出用满：低 32 位 → 前
   32 维、高 32 位 → 后 32 维；2026-08-08 升级，旧版为每 32 维 1 次取
   低 32 位），且不同 seed 的 mix 相互独立、天然 ILP；
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
采用延迟 −N 后完整耗时降到 **1.14 ms**（2026-08-08 64 位用满后再降至
**1.11 ms**）。i8 路径 RNG 调用是 bit 的 8 倍（每 8 维一次 vs 每 64 维
一次），RNG 独立成本更高，但同样与展开重叠，总耗时介于两者之间。

**展开替代路线实测否决（2026-08-08，n=10k）**："减少每 bit 指令数"的
两个替代方案均更慢：i16 lane 3.64 ms（见下方累加器实验，慢 ~3.3×）；
**pshufb 4-bit 查表 5.21 ms（慢 ~4.7×）**——字节级 LUT 展开需要标量
nibble 提取 + `vpmovzxbd` 字节→i32 扩展，每 bit op 数反而更多。srlv
家族（variable shift → and → shift-left → add）已接近该方向最优；剩余
杠杆在广播/依赖链（const 实验显示其占 ~98% 成本）。

**累加器宽度实验（2026-08-07）**：作者提示 bit 路径数值围绕 0 平衡、
i16 应够用。实测两种方案：

* **i16 lane 分段累加**（每段 8192 种子，段内值域确定地落在 [−8192, 8192]，
  无溢出）反而慢 **~2.5×**（3.50 ms vs 1.40 ms）：16-bit 变量移位
  （`vpsrlvw`）延迟更高，且 32 维块从 4 条独立累加链降到 2 条，链式
  依赖主导；
* **延迟 −N 偏移（采用）**：每种子只累加 2b、每块末尾一次性减 N——每组
  少 1 个向量 op，4 条链不变，bit 主路径 **1.40 → 1.14 ms**（n=10k，
  ~1.23×）、n=100k **14.0 → 11.4 ms**。

**2026-08-08 64 位框架复测**（`hash_hv_bit_i16`，n=10k）：每 64 维 4 个
i16 16-lane 寄存器（bits 0-15/16-31/32-47/48-63）+ 延迟 −N + 末尾
i16→i32 扩展转存，实测 **3.64 ms vs i32 主实现 1.11 ms，慢 ~3.3×**——
链 8→4 条 + `vpsrlvw` 延迟 + 转存开销（i32 路径无此步）。08-07 的 2.5×
（链 4→2）在新框架下扩大为 3.3×（链 8→4），"bit 数值平衡所以 i16 够用"
在数学上成立但工程上更慢，i32 + 延迟 −N 维持。

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

### 64 位用满升级（2026-08-08）

原生 64 位 RNG
（RapidRng / wyrand / splitmix64 / xoshiro256++ / PCG64）一次 mix 输出
64 位，旧版只取低 32 位是"信息没用上"而非"计算浪费"（64 位 mix 比 32 位
原生 RNG 的一次输出还便宜：RapidRng u64 0.51 ns vs Pcg32 u32 1.02 ns；
MT19937 / PCG32 / ChaCha12 的 u64 反而是两个 u32 拼出来的，见 rapidrand
表：Pcg32 u64 ≈ 2×u32、ChaCha12 u64 ≈ 2×u32）。把高 32 位也用上（每
64 维 1 次 RNG）可减半 RNG 调用——实测快 ~3%（1.1509 → 1.1145 ms，
n=100k 11.5 → 11.2 ms），与"RNG 免费"（const）持平；高低 32 位随机性
统计检查通过（各 ~50% 位密度、lo&hi 交叉 ~25%，独立平衡）。`hash_hv_bit`
无生产消费者（FASTA 路径走 i8、`.hv` v2 走稀疏），无兼容性负担，已直接
采纳：AVX2 / wide 回退 / serial 测试参考三处同步，逐位一致测试通过；
AVX-512 参考同步。i8 路径每 8 维用满 8 字节，无此问题。

### RNG 候选对比（2026-08-08）

为验证"换更快的 RNG 是否值得"，在 `hv_benchmark.rs` 中以与
`hash_hv_bit_avx2` 逐指令一致的宏生成变体（仅 RNG 行不同，RNG 输出统一
`black_box` 包裹以防编译器对常量 RNG 做 LICM 广播折叠），同机测量
（AVX2 64 位框架——每 64 维 1 次 RNG，与 2026-08-08 升级后的主实现同
口径；n=10k/100k，D=4096）。经典 RNG（MT19937 / LCG / PCG）无法塞进
O(1) 跳步框架的另行说明：

**O(1) 跳步候选（counter+mix 结构，块主序 AVX2 bit）**：

| RNG | n=10k | n=100k | vs 基线 |
|---|---:|---:|---:|
| **hash_hv_bit（现役主实现，无 black_box）** | **1.1145 ms** | **11.222 ms** | 1.00× |
| 常量（RNG 免费，广播保留） | 1.0935 ms | 10.934 ms | ~0.98× |
| seed 原值不混合（仅测量） | 1.1025 ms | 11.099 ms | ~0.99× |
| wyrand（HyperGen 编码 RNG） | 1.1227 ms | 11.254 ms | ~1.01× |
| rapid（宏版，同公式 + black_box） | 1.1266 ms | 11.283 ms | ~1.01× |
| splitmix64 | 1.1556 ms | 11.666 ms | ~1.04× |
| **PCG（O(log j) 跳步）** | **3.013 ms** | **30.14 ms** | **~2.7×** |
| **LCG/MINSTD（O(log j) 跳步）** | **10.69 ms** | **106.6 ms** | **~9.6×** |

**标量经典 RNG 对照（每 seed 串行流，i16 累加，n=10k）**：

| RNG | n=10k | vs RapidRng |
|---|---:|---:|
| **RapidRng（现役）** | **9.01 ms** | 1.00× |
| wyrand | 8.90 ms | ~0.99× |
| xoshiro256++（SmallRng） | 9.23 ms | ~1.02× |
| PCG32 | 9.50 ms | ~1.05× |
| xorshift64* | 9.57 ms | ~1.06× |
| splitmix64 | 9.92 ms | ~1.10× |
| ChaCha12（StdRng） | 10.38 ms | ~1.15× |
| **MT19937** | **25.82 ms** | **~2.9×** |

> **三组实现（AVX2 / wide / 标量，n=10k）**：bit 1.11 / 6.64 / 9.01 ms
> （AVX2/wide ~6.0×、AVX2/标量 ~8.1×）；i8 2.11 / 4.38 / 11.69 ms
> （AVX2/wide ~2.1×、AVX2/标量 ~5.5×）。wide vs 标量仅 1.4–2.7×
> （SIMD 只加速累加、RNG 串行主导）——"8 lane 只有约两倍"的观察即来自
> 此；AVX2 的跳步 RNG + 块主序额外拉开差距（见 hv.md §2.1）。

解读：

* **广播 + 依赖链是绝对主体**：常量 RNG（计算完全免费）与主实现只差
  ~2%，mix 计算本身约占 2%——换任何 RNG 的收益上限都 <3%。
* **快速候选全持平**：wyrand / rapid / 常量差异 <3%；splitmix64 慢 ~4%
  （64 位框架下 RNG 调用减半，其更重的 mix 影响被稀释）。换 RNG 无收益，
  维持原决策。
* **经典 RNG 全部落选**：MT19937 标量慢 ~2.9×（每 seed 需初始化 624 项
  状态数组 + temper，短流场景吃亏）；LCG/PCG 的跳步是 O(log j)，在块主序
  内层每 chunk 重算，慢 2.7×（PCG）/ 9.6×（LCG，MINSTD 模除更贵）。
  "O(1) 跳步（counter+mix）"是块主序框架的硬前提，经典 LCG/MT 结构不满足。
* **0.063 ms 是 LICM 假象**：不加 black_box 时，编译器把常量 RNG 的广播
  （`set1`）整体折叠/提升，测得 63 µs——不代表任何真实 RNG 可达。修正后
  常量与基线仅差 ~2%，进一步优化的真正方向是削减每 seed 广播/依赖链
  （非换 RNG）。

基准变体（`bit_avx2_rng_*`）保留在 `benches/hv_benchmark.rs` 供复现。

### 采样哈希吞吐（2026-08-08）

HyperGen 的 t1ha2 与 pgr 现有 minimizer 哈希对比（10k 个 21-mer，
一次性哈希）——

| 哈希 | 10k × 21-mer | vs rapidhash |
|---|---:|---:|
| **rapidhash（pgr 现役）** | **16.85 µs** | 1.00× |
| fxhash | 22.91 µs | ~1.36× |
| wyhash | 26.91 µs | ~1.60× |
| t1ha2（HyperGen FracMinHash 采样） | 29.71 µs | ~1.76× |
| murmurhash3 | 65.06 µs | ~3.86× |

pgr 现役 rapidhash 是最快采样哈希，t1ha2 慢 ~1.76×——不引入。

### 粒度换算（2026-08-08）

rapidrand 等外部基准报告的是**单次 draw**
（`next_u64()` 一次，纳秒级：RapidRng 0.51 ns、SmallRng 1.13 ns、
ChaCha20 13.24 ns）；pgr 基准报告的是**完整 HV 编码**（n=10k × D=4096
= 64 万次 draw（64 位用满后）+ 广播 + 向量展开，毫秒级）。两者交叉印证：
纯 RNG 独立成本 0.71 ms ÷ 128 万 ≈ **0.55 ns/draw**（旧 32 位框架分解），
与 rapidrand 的 0.51 ns 一致。标量对照中 RNG 间的差异被初始化/累加/内存
稀释（rapidrand 里 SmallRng 是 RapidRng 的 2.2×，在完整编码中只剩 ~2%）。

### 稀疏投影（2026-08-08，`.hv` v2 路径）

> 定位：稀疏编码是历史"尽全力提高速度"压力下的产物，非有意设计，当前
> 是 `.hv` v2 实际路径但**应视为待重新审视的候选**（见 hv.md §2.7 定位
> 说明与 §5.4 立项）。理论归属稀疏随机投影（Achlioptas 2003 / Li 2006 /
> feature hashing 2009），无直接产品先例。

`hash_hv_sparse`（splitmix64，每 seed 更新 s 个随机维度）性能
（n=10k、D=4096）：**s=3 0.055 ms**（s=1 0.022 / s=8 0.155 / s=16
0.365 / s=64 1.622 ms，近似随 s 线性）——比稠密 bit（1.11 ms）快
**~20×**、比 i8（2.11 ms）快 ~38×。n=100k、s=3：0.554 ms。**s 精度
扫描**（50 组独立集合对，N=3000/shared=500/D=4096，s=1..4096）：MAE
完全平坦（0.010–0.013）——**s 不是精度杠杆**（误差由 D 决定），只影响
速度；s=3 是历史默认无依据，速度最优为 s=1。固化于
`test_hash_hv_sparse_s_error_scan`。

**D 解耦实验（2026-08-08，s=1）**：编码时间与 D 几乎无关
（D=4096/16384/65536 → 0.022/0.023/0.024 ms，O(n·s)），而精度
`MAE ∝ 1/√D`（0.0155 / 0.0062 / 0.0030）——稀疏可用大 D 免费换精度，
稠密（O(n·D)）做不到。**大 s 是错误用法**：s=64 时 1.62 ms 比稠密
1.11 ms 慢（随机内存访问 vs 连续 SIMD），且无精度收益。详见 hv.md
§2.7 决策。

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
