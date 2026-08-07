# HV（Hypervector）设计笔记

> 2026-08-07 新增。来源：
> * HyperGen 论文：Xu et al., *Bioinformatics* 2024, 40(7), btae452，
>   DOI 10.1093/bioinformatics/btae452（PDF 见
>   `~/Downloads/Bioinformatics - 2024 - HyperGen compact and efficient genome s.pdf`）；
> * HyperGen 参考代码：仓库根目录 `Hyper-Gen-main/`（Rust，MIT，v0.2.2）。
>
> 本文汇总论文算法、代码实现，并对照 pgr 现有 HV 状态
> （`src/libs/hv.rs`、`dist hv` / `pgi to-hv`、`benches/hv_benchmark.rs`、
> [[../benchmarks/bench-simd-hv-jaccard.md]] 与
> [[../benchmarks/dist-cohort-validation.md]]），给出差异、启示与待验证项。

## 1. HyperGen 论文总览

### 1.1 定位

HyperGen 是一个**基于超维计算（HDC）的基因组草图 + ANI 估计**工具，目标是
大规模基因组集合（数据库搜索、聚类、物种分类）的快速粗筛。与 Mash /
Dashing 2 / Sourmash 等“离散 hash 集合”草图不同，HyperGen 把采样后的
k-mer hash 集合编码成高维整数向量（sketch HV），从而：

* 草图体积 O(D)（与 k-mer 集合大小 N 无关），比 Mash / Dashing 2 小 1.8–2.7×；
* 距离计算变成向量点积，可 SIMD / GEMM 化；
* 论文宣称 sketch 速度比 Mash 快 ~1.7×、搜索比 Dashing 2 快最多 4.3×，
  峰值内存 ~1 GB（FastANI / Skani 在 GTDB MAGs 规模 OOM）。

### 1.2 算法三步

#### 1.2.1 FracMinHash 采样

与 Mash 的 MinHash（保留最小的 N 个 hash）不同，HyperGen 用 FracMinHash：

```
S_k(A) = { h(x) | x ∈ A, h(x) ≤ M/S }
```

其中 M 为 hash 值域上限，S 为 scaled factor（默认 1500）。FracMinHash 对
大小差异悬殊的集合给出**无偏** Jaccard 估计（MinHash 有偏），代价是采样集
通常更大——这正是后续 HV 编码要压缩的部分。默认 canonical k-mer、k=21。

#### 1.2.2 HDC 编码（关键步骤）

每个被采样的 k-mer hash 值作为种子喂给伪随机数生成器（代码里用
`WyRng::seed_from_u64(hash)`），生成 D 维二进制向量 hv ∈ {0,1}^D，
转成 ±1 后逐位累加：

```
H = Σ_i (hv_i × 2 − 1)        （等价于：hv 初值 −N，每位置位 +2）
```

默认 D=4096。要点：

* 每个 k-mer 影响**所有** D 维（稠密、全维更新）；
* 随机位保证不同 k-mer 的向量准正交；
* 聚合结果 H ∈ [−N, N]^D（代码用 i16 存，值域限制见 §2.4）。

#### 1.2.3 Jaccard / ANI 估计

利用准正交性：

```
|S_k(A)| = ‖H_A‖² / D
|S_k(A) ∩ S_k(B)| = H_A · H_B / D
J = dot / (‖H_A‖² + ‖H_B‖² − dot)          （D 约掉，代码直接算原始点积/范数）
ANI = 1 + (1/k)·ln(2J / (1+J))
```

后一条即 Mash 的距离公式（ANI = 1 − Mash distance）。L2 范数在 sketch
阶段预计算并存储，比较时只算点积。

### 1.3 默认参数与论文结论

* 默认 `k=21, S=1500, D=4096, seed=123`；D>4096 后误差不再显著下降；
* **S 越小（采样越密）误差越大**：聚集进 HV 的向量越多，正交性被破坏得越
  厉害——这是论文 §3.2.1 的核心观测，也是理解我们饱和问题的钥匙（§5.1）；
* 细菌数据集（Bacillus cereus / E. coli）D=4096 时 MAE 0.37/0.37、Pearson
  0.986/0.952，为所有 sketch 工具中最低档；
* GTDB MAGs 搜索：sketch 阶段 130.4 s / 1.0 GB，单查询 0.3 s / 0.9 GB
  （对照 Dashing 2 sketch 5632 s / 8.9 GB）；
* GPU fast 模式（CUDA）比 CPU 再快 1.8–2.7×。

## 2. 代码实现梳理（Hyper-Gen-main）

### 2.1 模块结构

| 文件 | 职责 |
|---|---|
| `src/main.rs` | CLI 入口，sketch / dist / search（search 未实现，TODO） |
| `src/utils.rs` | clap 参数、glob 收集 FASTA、进度条、bincode 读写 sketch、ANI 输出 |
| `src/sketch.rs` | 并行 sketch：`extract_kmer_hash`（FracMinHash 采样）→ 编码 → 压缩 |
| `src/hd.rs` | `encode_hash_hd`（标量 i16）、`encode_hash_hd_avx2`（SIMD）、无损量化 + bitpack 压缩/解压 |
| `src/dist.rs` | L2 范数、点积、Jaccard→ANI；对称/非对称全对距离 |
| `src/types.rs` | `FileSketch` / `SketchParams` / `SketchDist`，mm_hash64 等 |
| `src/sketch_cuda.rs` | CUDA k-mer 哈希内核（feature 门控） |

### 2.2 关键实现细节

* **采样**：`needletail` 的 `canonical_kmers`（k=21）+ `t1ha2_atonce(kmer, seed)`
  哈希，保留 `h < u64::MAX/scaled` 的 canonical hash 进 `HashSet<u64>` 去重；
  另有 `mm_hash64`（MurmurHash 变体）供 CUDA 路径用。
* **编码**：
  * 标量：`hv` 初值 `vec![-(N as i16); D]`，每 seed 按 64 维一块取
    `WyRng::next_u64()`，逐位把 `((bit << 1) as i16)` 累加进 i16；
  * AVX2：每批 4 个 seeds × 每 64 维一块，`_mm256_shuffle_epi8` 把 4 个 u64
    重组为 4×16-bit lanes，逐位 `srl + and + hadd` 得 4 维计数，结果与标量
    逐位一致（单元测试 `test_simd_hd_enc` 断言相等）。
* **压缩（无损）**：`compress_hd_sketch` 先找覆盖 `[min,max]` 的最小位宽
  （6→16 逐级试），再加偏移用 `bitpacking::BitPacker8x` 按 256 元素块压成
  位串；`hv_quant_bits` 存进文件头，读取时 `decompress_hd_sketch` 还原。
  论文称再省 ~1.3×。
* **距离**：`compute_pairwise_ani`（标量；另有 AVX2 `_mm256_madd_epi16`
  版本）算原始 i32 点积，`jaccard = dot/(norm2_r+norm2_q−dot)`，
  `ani = 1 + ln(2/(1/J+1))/k`，clamp [0,1]×100、NaN→0；对称模式只算上三角。
* **输出**：按 ANI 降序、过滤 `ani_threshold`（默认 85.0）写 TSV。
* **文件格式**：整批 genome 的 `Vec<FileSketch>` 直接 `bincode::serialize`
  成一个文件（含 ksize/scaled/seed/hv_d/hv_quant_bits/hv_norm_2/file_str/hv）。

### 2.3 工程观察

* CLI `sketch` 输入是**目录**（glob `*.fna/*.fa/*.fasta`），`dist` 输入是两个
  sketch 文件；`search` 子命令已在 clap 里占位但无实现；
* `if_compressed` 在 CLI 里硬编码 `true`（TODO）；
* `SketchDist::default` 的 `hv_d=1024` 与 sketch 默认 4096 不一致，但 dist
  实际从文件头读 hv_d，默认值不参与计算（仅 ksize/hv_d 一致性断言）；
* 依赖：`needletail`（k-mer）、`wyhash`/`t1ha`（哈希/RNG）、`bitpacking`、
  `bincode`、`rayon`、`rand`、`xxhash-rust`、`bloom`、`ndarray`、
  `cudarc`（可选，CUDA）。

### 2.4 代码侧风险（对 pgr 有参考意义）

* **i16 值域**：每维值域 [−N, N]，N > 32767 即溢出（`-N as i16` 直接截断）。
  论文评测为细菌规模（5 Mb / S=1500 ≈ 3.3k），无虞；人类规模或小 S 会爆。
* **sketch 文件无版本/magic**：bincode 直接反序列化，跨版本不兼容且无自描述
  （pgr 的 `.hv` 有 `PGV1` magic + version，见 §3.2）。
* 搜索未实现：论文的“GEMM 加速搜索”只存在于实验代码，仓库里是 TODO。

## 3. pgr 现状

### 3.1 `src/libs/hv.rs`（486 行）

| 函数 | 编码 | 说明 |
|---|---|---|
| `hash_hv_bit` | 稠密位编码，i32 | hv 初值 −N，每 32 维一个 `RapidRng::next_u32`，位 0/1 加 2；AVX2 8-lane |
| `hash_hv_i8` | 稠密 i8 累加，i32 | hv 初值 0，每 8 维一个字节（i8）；RNG 调用少（1 u64/8 维） |
| `hash_hv_sparse` | 稀疏 ±1，i32 | splitmix64 派生，每 k-mer 只更新 `s`（默认 3）个随机维度 |
| `hv_norm_l2_sq` / `hv_cardinality` / `hv_dot` | — | wide SIMD 范数；cardinality=‖H‖²/D；dot 按 √D 归一 |
| `calc_distances` | — | jaccard / containment / mash 多口径输出 |
| `load_hv_from_fasta(_syncmer)` | i8 | FASTA → minimizer / closed syncmer 集合 → `hash_hv_i8` |

注意：pgr 的稠密编码用 **RapidRng**（对照 HyperGen 的 WyRng），且累加用 i32
（无 HyperGen 的 i16 上限问题，但存储 2×）。

### 3.2 消费链与 `.hv` 格式

* `pgr dist hv`：FASTA 路径（minimizer/syncmer → `hash_hv_i8` →
  `calc_distances`）；`.hv` 路径（`pgi to-hv` 产物直接比较，稀疏余弦）。
* `pgr pgi to-hv`：把 `.pgi` 的 unique k-mer keys 投影成**稀疏** HV；
  `.hv` v2 格式：`PGV1` magic + version 2 + `k/dim/sparse/n_kmer/name` +
  i32 数组。稀疏投影 + 存储 n_kmer 是 v2 的关键修复（见下）。
* `pgr dist hv a.hv b.hv`：余弦相似度 → `inter = cos·√(n1·n2)`，
  集合大小用文件头存储的真实 n_kmer。

### 3.3 基准与验证（已有记录）

* [[../benchmarks/bench-simd-hv-jaccard.md]]（2026-08-06，Ryzen 9 7945HX）：
  `hash_hv_i8` 1k/10k seeds 447.8 µs / 4.452 ms，`hash_hv_bit` 677 / 6.818 ms，
  标量 RNG 对照 ~2× 更慢；`norm_l2` SIMD 7.8×；集合 Jaccard rapidhash 最快。
  遗留疑虑：HV 编码 SIMD 相对标量仅 ~1.5–2×（[[../todo.md]] §4）。
* [[../benchmarks/dist-cohort-validation.md]]（10 株 E. coli × 45 对）：
  * 稠密 i8 在 260 万 seeds 下**饱和**（各维 ±1.3e6，dim 无关），containment
    饱和为 1.0、与身份率 Spearman ≈ 0（初测 −0.05）；
  * **稀疏 v2 修复**：与 `dist pgi` mash 排序 Spearman 0.969、45 对 0.12 s
    （~71× 提速）、共享 k-mer 计数平均误差 2.39%；
  * `dist seq`（k=8 syncmer 草图）仍是与身份率最贴近的草图层（ρ=0.616）。

### 3.4 推导发现（2026-08-07，待实测确认）

FASTA 路径 `hash_hv_i8` + `calc_distances` 存在**量纲不匹配**：

* i8 字节均值 ≈ −0.5（0..255 → −128..127），所以每维 H 的期望 ≈ −N/2；
* E[H_i²] ≈ N·E[b²] + N(N−1)·E[b]² = N·5390 + N(N−1)/4，于是
  `hv_cardinality = ‖H‖²/D` ≈ N·6140（不是 N）；
* 点积被 ~N²/4 的“直流项”主导，`inter = hv_dot = dot/D` 与 cardinality 的
  缩放不一致，Jaccard 随 N 增大趋向只依赖集合大小的常数（等大时 ~0.5），
  mash 失去区分度。

数值模拟（D=4096，N=3000，shared=500）：reported jaccard 0.147 vs 真值
0.091；N 更大时偏差更严重。稀疏 `.hv` 路径（余弦 + 真实 n_kmer）不受影响。
`hash_hv_bit`（±1 编码）无此问题——该路径建议改用 bit 编码或稀疏投影，
并用两株 E. coli 对照 `dist seq` / `dist pgi` 实测确认后再定。

## 4. HyperGen vs pgr 对照

| 维度 | HyperGen | pgr 现状 |
|---|---|---|
| 采样 | FracMinHash（S=1500，canonical k=21） | minimizer / closed syncmer（FASTA 侧）；`.pgi` unique k-mers（索引侧） |
| 编码 | 稠密 ±1 累加 i16（WyRng，4 seeds×64 bits SIMD） | 稠密 bit/i8 累加 i32（RapidRng）；稀疏 ±1（v2 修复） |
| 维度 | 4096 | 4096 默认（`dist hv` / `pgi to-hv`） |
| 距离 | Jaccard→ANI（预存 L2 范数 + 点积） | jaccard / containment / mash 多口径；稀疏时余弦→inter 估计 |
| 压缩 | 无损最小位宽量化 + BitPacker8x | 无（`.hv` 原始 i32，4096 维 16 KB/文件） |
| 批量 | 目录→一个 sketch 文件，ref×query 全对 | `.hv` 单文件，两两比较；FASTA 侧同 `dist seq` 参数 |
| 搜索 | 论文宣称 GEMM；仓库未实现 | `dist hv` / `dist pgi` / `dist seq` 均已落地 |
| GPU | CUDA fast 模式（feature 门控） | 无 |
| 值域 | i16，N ≤ 32767 | i32，无此限制 |
| 文件格式 | bincode，无 magic/版本 | `PGV1` magic + version（v2） |

## 5. 启示与后续方向

### 5.1 饱和问题的理论解释

稠密 ±1（或 i8）编码下，每个 k-mer 更新所有 D 维，N 个随机向量聚合后各维
标准差 ~√N，共享信号与无关噪声的比值随 N 增大而恶化。HyperGen 用
FracMinHash 把 N 控制在 ~基因组长度/S（E. coli ~3.3k）来保住正交性；
论文“S 越小误差越大”的观测与此同源。pgr 的稀疏投影（每 k-mer 只碰 s=3 维）
是**另一条正交性控制路线**，已在 `.hv` v2 实证有效（ρ=0.969）。

### 5.2 SIMD 提速疑虑（2026-08-07 已解决）

根因是 **RNG 生成主导耗时**（n=10k、D=4096：RNG-only 2.68 ms vs 总
4.36 ms），8-lane AVX2 只加速位展开、不减少 RNG 调用（i8 每 8 维一次
128-bit mix），也消除不了串行 next_u64 的依赖链。解法：

1. **RapidRng 跳步 + 块主序**：RapidRng 状态是常数步长计数器，输出
   j = mix(seed + j·SECRET0, …)，因此 HV 分块可常驻寄存器、遍历全部 seed，
   每 seed 每 32 维只做 1 次 RNG，且不同 seed 的 mix 相互独立 → ILP；
2. **AVX2（256-bit）位展开**：`vpsrlvd + vpand + vpslld + vpsubd`
   四个 8-lane 寄存器一次覆盖 32 维（intrinsics，避免 wide 的冗余展开）。

落地：`hash_hv_bit` / `hash_hv_i8` 以 **AVX2 为主实现**（x86-64 运行时
检测 `avx2`，无则降级 wide 可移植路径；输出与串行逐位一致），实测
bit ±1 编码 1.14 ms（n=10k、D=4096，vs 旧 wide bit 6.69 ms ~5.9×、
vs 旧 wide i8 4.33 ms ~3.8×）、i8 保语义 2.09 ms（~2.1×）。
累加用 i32，但利用"±1 数值围绕 0 平衡"做**延迟 −N 偏移**：每种子只累加
2b、每块末尾一次减 N（每组少 1 个向量 op）。曾试 i16 lane 分段累加
（每段 8192 种子、段内值域确定不溢出），实测反而慢 ~2.5×（16-bit 变量
移位延迟高、独立累加链从 4 条降到 2 条），未采用；i8 路径确认必须 i32
（字节均值 −0.5 直流偏置，260 万种子实测 ±1.3e6）
（详见 [[../benchmarks/bench-simd-hv-jaccard.md]] §2）。
**AVX-512 仅保留在 `benches/hv_benchmark.rs` 作参考对照，不参与运行时
分派**（作者决策：AVX-512 不优于 256；Zen 4 实测两者持平，参考列数字
见 [[../benchmarks/bench-simd-hv-jaccard.md]] §2）。

**平台策略**：AVX2 仅是 x86-64 的运行时可选加速，不引入平台依赖——非
x86-64 目标（wide 在 aarch64 自动用 NEON）和不支持 AVX2 的 x86 CPU 都走
可移植 wide 路径，所有路径输出逐位一致；
`cargo check --target aarch64-apple-darwin --lib` 交叉编译通过（2026-08-07）。

### 5.3 pbit 决策 B 的参考方案

[[pbit.md]] 决策 B（HV sketch 内嵌）触发条件是“无源 FASTA、仅归档、需距离
粗筛”——正是 HyperGen 的完整场景。触发时可直接参考：FracMinHash 采样 +
稠密 ±1 编码（i16 或 i32）+ 无损量化 bitpack + 预存范数 + 点积 ANI，并
注意 §2.4 的 i16 值域限制与 batch 全对输出。

### 5.4 可借鉴工程点

* 无损量化 + bitpack：pgr `.hv` 现在存原始 i32，可无损压到 ~6–16 bit/维；
* 预计算 L2 范数 / 存储 n_kmer（pgr 稀疏版已做后者）；
* 对称模式只算上三角、输出按 ANI 排序 + 阈值过滤（`dist` 输出习惯可对齐）；
* 批量搜索的 GEMM 化（论文 future 方向，pgr `linalg` 已是 wide SIMD）。

### 5.5 待验证 / 待办

1. **FASTA `dist hv` 路径量纲问题**（§3.4）：先实测确认，再决定改
   `hash_hv_bit` / 稀疏投影；
2. HV 编码 SIMD 深挖（todo.md §4）：**已解决（2026-08-07，见 §5.2）**；
3. 若决策 B 立项，按 §5.3 做设计。

## 参考

* [[pbit.md]]（决策 B 与 `.hv` 消费链）
* [[../benchmarks/bench-simd-hv-jaccard.md]]（HV 编码基准）
* [[../benchmarks/dist-cohort-validation.md]]（饱和问题与稀疏 v2 验证）
* [[../todo.md]]（§4 HV SIMD 疑虑）
* `Hyper-Gen-main/`（参考代码）与 `~/Downloads/...btae452.pdf`（论文）
