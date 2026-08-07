# HV 外部参考（HyperGen + hdlib）

> 整理于 2026-08-07（自 `design/hv.md` §4 迁出），后续将随更多
> HDC / 基因组草图文献持续扩充。目的：了解超维计算（HDC）基因组草图的
> 主流做法与参数选择，为 pgr 的 HV 实现（`src/libs/hv.rs`、
> [[../design/hv.md]]）提供外部参考。本文是背景材料；pgr 的实现决策以
> [[../design/hv.md]] §1/§2 为准。
>
> 来源：
> * HyperGen 论文：Xu et al., *Bioinformatics* 2024, 40(7), btae452
>   （PDF 见 `~/Downloads/Bioinformatics - 2024 - HyperGen compact and efficient genome s.pdf`）；
> * HyperGen 参考代码：仓库根目录 `Hyper-Gen-main/`（Rust，MIT，v0.2.2）；
> * hdlib 参考代码：仓库根目录 `hdlib-2.0.0/`（VSA 通用库，Python，JOSS 2023）。

## 1. HyperGen 论文算法总览

**定位**：基于超维计算（HDC）的基因组草图 + ANI 估计，面向大规模基因组
集合的快速粗筛。草图体积 O(D)（与 k-mer 集合大小 N 无关），距离计算变成
向量点积（可 SIMD / GEMM 化）；论文宣称 sketch 比 Mash 快 ~1.7×、搜索比
Dashing 2 快最多 4.3×、峰值内存 ~1 GB。

**算法三步**：

1. **FracMinHash 采样**：保留 `h(x) ≤ M/S` 的 canonical k-mer（默认
   k=21、S=1500）。对大小差异悬殊的集合仍能给出**可校正**的
   Jaccard / containment 估计（MinHash 对大小悬殊集合有偏；偏差性质与
   校正见 [[../design/hv.md]] §2.6），代价是采样集更大。
2. **HDC 编码**：每个被采样 k-mer 的 hash 作为种子，用 `WyRng` 生成
   D 维二进制向量，转 ±1 后逐位累加（`H = Σ(hv×2−1)`，默认 D=4096）。
   每个 k-mer 影响所有 D 维（稠密、全维更新）。
3. **Jaccard / ANI**：`|S|=‖H‖²/D`、`|A∩B|=H_A·H_B/D`、
   `J = dot/(‖H_A‖²+‖H_B‖²−dot)`、`ANI = 1 + ln(2J/(1+J))/k`
   （即 Mash 公式的 ANI = 1 − Mash distance）。L2 范数预计算存储。

**默认参数与结论**：`k=21, S=1500, D=4096, seed=123`；D>4096 后误差不再
显著下降；**S 越小（采样越密）误差越大**（聚集向量越多，正交性被破坏，
与 pgr 饱和问题同源）；细菌数据集 D=4096 时 MAE 0.37、Pearson 0.95+；
GTDB MAGs 搜索 sketch 130.4 s / 1.0 GB、单查询 0.3 s / 0.9 GB；GPU
fast 模式再快 1.8–2.7×。

## 2. HyperGen 代码实现梳理（Hyper-Gen-main）

| 文件 | 职责 |
|---|---|
| `src/main.rs` | CLI 入口，sketch / dist / search（search 未实现，TODO） |
| `src/utils.rs` | clap 参数、glob 收集 FASTA、bincode 读写 sketch、ANI 输出 |
| `src/sketch.rs` | 并行 sketch：FracMinHash 采样 → 编码 → 压缩 |
| `src/hd.rs` | `encode_hash_hd`（标量 i16）、AVX2 版（4 seeds×64 bits）、无损量化 + bitpack |
| `src/dist.rs` | L2 范数、点积、Jaccard→ANI；对称/非对称全对距离 |
| `src/types.rs` | `FileSketch` / `SketchParams` / `SketchDist`，mm_hash64 等 |
| `src/sketch_cuda.rs` | CUDA k-mer 哈希内核（feature 门控） |

关键细节：

* 采样：`needletail` canonical_kmers + `t1ha2_atonce`，保留
  `h < u64::MAX/scaled` 进 `HashSet` 去重；
* 编码：标量 i16（hv 初值 −N，每 64 维一个 u64 随机位）；AVX2 版每批
  4 seeds × 64 位，`shuffle + srl + and + hadd`，与标量逐位一致；
* 压缩：找覆盖 `[min,max]` 的最小位宽（6→16）加偏移，
  `BitPacker8x` 按 256 元素块压位串，读取时解压（论文称再省 ~1.3×）；
* 距离：i32 点积 → jaccard → ANI，clamp [0,1]×100、NaN→0；对称模式
  只算上三角；输出按 ANI 降序、过滤阈值（默认 85.0）；
* 文件格式：整批 `Vec<FileSketch>` 直接 bincode（无 magic/版本）。

**代码侧风险（对 pgr 有参考意义）**：

* **i16 值域**：每维值域 [−N, N]，N > 32767 溢出；论文评测为细菌规模
  （5 Mb / S=1500 ≈ 3.3k）无虞，人类规模或小 S 会爆（pgr 用 i32 规避）；
* **sketch 文件无版本/magic**：跨版本不兼容且无自描述（pgr `.hv` 有
  `PGV1` magic + version）；
* 搜索未实现：论文宣称的 GEMM 加速搜索只存在于实验代码。

## 3. 对照表（HyperGen vs pgr）

| 维度 | HyperGen | pgr 现状 |
|---|---|---|
| 采样 | FracMinHash（S=1500，canonical k=21） | minimizer / closed syncmer（FASTA 侧）；`.pgi` unique k-mers（索引侧） |
| 编码 | 稠密 ±1 累加 i16（WyRng，4 seeds×64 bits SIMD） | 稠密 bit/i8 累加 i32（RapidRng，跳步块主序）；稀疏 ±1（v2） |
| 维度 | 4096 | 4096 默认（`dist hv` / `pgi to-hv`） |
| 距离 | Jaccard→ANI（预存 L2 范数 + 点积） | jaccard / containment / mash 多口径；稀疏时余弦→inter |
| 压缩 | 无损最小位宽量化 + BitPacker8x | 无（`.hv` 原始 i32） |
| 批量 | 目录→一个 sketch 文件，ref×query 全对 | `.hv` 单文件两两比较；FASTA 侧同 `dist seq` 参数 |
| 搜索 | 论文宣称 GEMM；仓库未实现 | `dist hv` / `dist pgi` / `dist seq` 已落地 |
| GPU | CUDA fast 模式（feature 门控） | 无 |
| 值域 | i16，N ≤ 32767 | i32，无此限制 |
| 文件格式 | bincode，无 magic/版本 | `PGV1` magic + version（v2） |

## 4. hdlib 参考（VSA 通用库，Python）

hdlib（Cumbo et al., JOSS 2023，仓库根目录 `hdlib-2.0.0/`）是通用
VSA/HDC 库：`Space`（随机向量容器）、`Vector`（binary/bipolar 随机向量 +
cosine/hamming/euclidean 距离）、`arithmetic`（bind/bundle/permute）、
`model.graph`（图编码，Poduval 2022）、`model.classification/clustering`。
`examples/pangenome/minimizers.cpp` 只有朴素的 minimizer 提取
（O(n·w) 字符串比较），无借鉴价值。值得留意/评估的点：

* **bind + permute 编码结构（未立项）**：`bind(node_u, permute(node_v))`
  编码有向边、bundle 全部边编码整图——把"顺序/邻接/结构"压进 HV 的标准
  VSA 手法。pgr 当前 HV 是**无序集合语义**（Jaccard 的前提），不需要；
  若未来出现"顺序敏感草图"（k-mer 沿基因组顺序、contig 邻接、pangenome
  图分类）需求，这套 bind/permute 是现成方案。
* **带权重/多重度的 bundle（可考虑）**：同一向量在 bundle 中出现多次 →
  结果偏向主导向量。pgr 现在聚合是去重集合（无权重），重复内容已被证明
  干扰排序（见 [[../benchmarks/dist-cohort-validation.md]] §2 的 k=40 集合
  受重复影响）；若未来要保留 k-mer 多重度（类似 Dashing 2 SetSketch），
  加权聚合是一条路，但需先定义"权重=什么"（拷贝数？覆盖度？）。
* **binary/bipolar 双类型与三种距离**：对 ±1 向量，cosine 与 hamming
  等价（cos = 1 − 2·hamming/D）；我们主路径已用 cosine（.hv v2）与
  Jaccard/Mash，无新增。
* 随机向量生成与 seed 复现：与我们一致，无新增。
