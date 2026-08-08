# HV 与测距聚类外部参考

> 整理于 2026-08-07（自 `design/hv.md` §4 迁出），2026-08-08 并入测距/
> 聚类文献（来源 `~/sync/zotero/bacteria/clustering/`），后续随更多文献
> 持续扩充。目的：为 pgr 的 HV 实现（`src/libs/hv.rs`、[[../design/hv.md]]）
> 与细菌测距/聚类方向提供外部参考。本文是背景材料；pgr 的实现决策以
> [[../design/hv.md]] §1/§2 与 §6 审计为准。
>
> 来源：
> * HyperGen 论文：**Bioinformatics 2024（HyperGen）**
>   （PDF 见 `~/sync/zotero/bacteria/clustering/Bioinformatics - 2024 - HyperGen compact and efficient genome s.pdf`）；
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

hdlib（JOSS 2023，仓库根目录 `hdlib-2.0.0/`）是通用
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

## 5. 测距与聚类文献（bacteria clustering）

> 来源：`~/sync/zotero/bacteria/clustering/`（作者收集，准备"测距用于
> 聚类"方向）。与 HV 直接相关的条目被 [[../design/hv.md]] §2.6/§2.7 与
> §6 审计引用。

### 5.1 集合相似度与距离估计（测距核心）

| 文献 | 定位 | 与 pgr 的关系 |
|---|---|---|
| **bioRxiv 2022（FracMinHash, Irber）** | FracMinHash 采样 + 最小 metagenome cover | [[../design/hv.md]] §2.6 数值 ANI 的推荐采样器 |
| **bioRxiv 2023（open syncmer ≡ FracMinHash, Koslicki）** | open syncmer 与 FracMinHash 在 Jaccard/containment 意义上等价 | hv.md §2.6 已有引用 |
| **Bioinformatics 2022（minimizer Jaccard 有偏）** | minimizer Jaccard 估计**有偏且不一致** | hv.md §2.6 已有引用 |
| **PeerJ 2021（closed syncmer, Edgar）** | closed syncmer 定义，比 minimizer 更敏感 | pgr syncmer 采样来源（见 [[syng.md]]） |
| **Bioinformatics 2023（minmer）** | minmer 泛化 minimizer，修正其偏差 | hv.md §2.6 已有引用 |
| **KDD 2023（DotHash）** | 超维向量（HDC）估计集合相似度：**Theorem 2 证明点积无偏估计交集 + 误差概率界** | 与 pgr HV 最直接：稀疏/稠密投影的同类工作，§6 审计的核心依据 |
| **IEEE TKDE 2020/2022（ProbMinHash）** | 概率 Jaccard（带权重集合）的 LSH 族 | 若未来支持多重度/加权 k-mer（hv.md §4.4 weighted bundle 方向）可参考 |
| **Bioinformatics 2019（Order Min Hash, Marçais）** | edit distance 的 LSH（Order Min Hash） | pgr 距离方向的 LSH 参考 |
| **Genome Res 2021（strobemers, Sahlin）** | 成组短 k-mer 采样，抗 indel | 采样层候选（对比 syncmer/minimizer） |
| **Genome Biol 2024（minimizer sketching 综述）** | minimizer sketching 综述（何时用、理论、局限） | 采样方法全景，补 hv.md §2.6 背景 |
| **Bioinformatics 2022（local k-mer selection 理论）** | local k-mer selection 理论：conservation 精确表达式 + syncmer 闭式解 + minimap2 实证 8.2% | 采样器选择的定量理论框架（§6 审计指出 hv.md 未引用） |
| **Bioinformatics 2022（压缩 k-mer 字典）** | 基于 minimizer 统计性质的压缩 k-mer 字典 | k-mer 集合压缩存储方向 |
| **PODS/JCSS 2001/2003（随机投影, Achlioptas）** | 数据库友好随机投影：±1 / 稀疏（2/3 零）矩阵保距（Theorem 2） | 稀疏随机投影奠基（§6.6） |
| **KDD 2006（very sparse random projections, Li）** | very sparse random projections：±√s、非零概率 1/s，保距只需零均值+单位方差 | **与 pgr 稀疏 s 桶期望等价**（§6.6） |
| **ICML 2009（feature hashing）** | feature hashing：元素→桶+符号，指数尾界 | pgr 稀疏 s=1 的直接对应（§6.6） |
| **J. Algorithms 2005（Count-Min Sketch, Cormode）** | Count-Min Sketch：每元素多桶，内积查询 (ε,δ) 保证 | 与 pgr 稀疏结构同构（§6.6） |

### 5.2 大规模聚类与搜索（测距的消费者）

| 文献 | 定位 | 与 pgr 的关系 |
|---|---|---|
| **Genome Biol 2023（RabbitTClust）** | 百万细菌基因组快速聚类（含 k-mer 距离 + 聚类） | pgr 4 万 E. coli cohort 聚类的直接对标 |
| **Nat Biotechnol 2019（BIGSI）** | 全部已公开细菌/病毒基因组索引与搜索 | 大规模索引搜索的里程碑 |
| **Nat Biotechnol 2025（LexicMap）** | 百万原核基因组高效比对（k-mer + 索引） | 大规模 pairwise 比对的工程参考 |
| **NAR 2024（GSearch）** | k-mer hashing + HNSW 图做基因组搜索 | pgr dist 的近似搜索（ANN）方向 |
| **IEEE TPAMI 2020（HNSW）** | 分层可导航小世界图 ANN 搜索 | 距离矩阵之外的图式近邻搜索经典 |
| **Nat Commun 2018（MMseqs2 聚类）** | 线性时间蛋白序列聚类（Linclust） | 蛋白侧聚类（若 pgr 扩展蛋白距离） |
| **NAR Genom Bioinform 2024（近似近邻图嵌入）** | 近似近邻图 + 嵌入，大规模生物数据 | embedding + ANN 结合方向 |
| **Bioinformatics 2023（MetaProFi）** | chunked Bloom filter 存储查询蛋白/核酸序列 | 大规模集合存储查询方向 |

### 5.3 k-mer 哈希与采样

| 文献 | 定位 | 与 pgr 的关系 |
|---|---|---|
| **Bioinformatics 2016（ntHash）** | 递归核苷酸哈希（滚动） | pgr syncmer 乘性滚动哈希的对照 |
| **Bioinformatics 2022（ntHash2）** | 递归 spaced seed 哈希 | spaced seed 采样方向 |
| **Bioinformatics Advances 2023（aaHash）** | 蛋白 k-mer 滚动哈希 | 蛋白距离（dist seq protein）方向 |
| **WABI 2024（mod-minimizer）** | 长 k-mer 的简单高效采样 | 采样器候选（长 k-mer 场景） |
| **Bioinformatics 2023（seeding with minimized subsequence）** | 最小化子序列做种子 | 采样理论补充 |
| **CSBJ 2024（k-mer 方法综述）** | k-mer 方法与应用全景 | 方向总览入口 |

### 5.4 超维计算

| 文献 | 定位 | 与 pgr 的关系 |
|---|---|---|
| **Cogn Comput 2009（HDC 奠基, Kanerva）** | HDC 奠基：高维随机向量准正交 + bind/bundle/permute | pgr HV（§1–§4 对照）的理论源头 |
| **IEEE ICCD 2019（SparseHD）** | 稀疏超维计算：训练后模型稀疏化 90%、质量损失 <1%，FPGA 加速 48.5× 低能耗 | 稀疏 HD 方向；与 pgr 编码阶段稀疏机制不同（§6.6） |

### 5.5 未收录（相关性弱或重复）

* Manifold Learning 综述（两份，降维非测距核心）、WebGraph（图压缩）、Statist Med 1996
  （logistic regression）、wheat pangenome 应用、UniProt 去冗余、SECOM/domain identity/
  compressed amino acid alphabets（蛋白域方向）、Flexible protein database、Snekmer、
  Improved protein homolog、Matchtigs/Simplitigs/CBL/BWT/Strobealign/Block Aligner/
  Exact global alignment（比对与表示，非测距聚类核心，可在需要时补充）。
* ProbMinHash、Syncmers、Ultrafast search 各有重复 PDF（同一文献的两个版本）。

## 6. 文献阅读笔记（2026-08-08，逐篇完整阅读）

> 目的：为 [[../design/hv.md]] §6 审计提供逐篇的详细依据。笔记按
> 阅读批次组织，核心关注：方法与理论（定理/误差界）、与 pgr HV 各决策
> 的关系、以及"稀疏投影能否获得 DotHash 级理论"这一关键问题。
>
> **阅读深度说明（诚实标注）**：与 HV 决策直接相关的核心文献（§6.1/§6.2）
> 完整精读（含定理/方法/证明）；聚类/搜索/比对/蛋白等背景文献（§6.3/§6.4）
> 精读摘要 + 方法 + 结果；个别与方向弱相关的仅记录定位。笔记中"完整
> 精读"均指已读正文核心（非仅摘要）。

### 6.1 第一批：HV 编码与集合距离核心

#### DotHash（KDD 2023）—— 完整精读

**核心构造**：元素经随机映射 ψ: S → R^d 到 **d 维超立方体顶点**（unit
向量，即 ±1/√d 的稠密随机向量），集合 sketch = 元素向量之和。

**Theorem 1（精确版）**：若 ψ 是标准基（one-hot），a·b = |A∩B| 精确成立。

**Theorem 2（估计版，pgr 最直接的理论模板）**：ψ 均匀随机映射到超立方体
顶点，则：

* `E[a·b] = |A∩B|`（**无偏**）；
* `Var(a·b) = (1/d)·[|A||B| + |A∩B|² − 2|A∩B|]`；
* 误差界：Chebyshev `Pr(|X−μ| ≥ εμ) ≤ Var/(ε·|A∩B|)²`，另有 CLT 近似。

**关键点（对 pgr 的意义）**：

1. **这是"元素随机向量叠加 → 点积无偏估计交集"的完整理论**，直接支持
   pgr 稠密 bit 路径（§2.1）的核心假设——pgr 的稠密 bit 与 DotHash
   Theorem 2 的构造等价（±1 超立方体、集合取和、点积）。**pgr 稠密
   路径可以直接套用这个无偏性和方差公式**。
2. **稀疏版本不是 s 桶 ±1**：DotHash 的"稀疏"是 standard basis
   （one-hot，Theorem 1，精确计算、无符号），不是 pgr 的"每元素 s 个
   随机 ±1 桶"。**pgr 的 s 桶构造没有现成定理**——但 Theorem 2 的
   证明路径（按维度独立 + 区分相等/不相等元素对）可以仿照推导 s 桶
   的方差（见 §6.6 待补推导）。
3. 扩展：可估计 Adamic-Adar 等族（通过调整元素向量幅度）——pgr 若做
   泛基因组图/链路预测可借鉴。

#### HDC 奠基（Kanerva, Cogn Comput 2009）—— 完整精读

**核心内容**：高维随机向量空间的统计性质 + bind/bundle/permute 操作。

* **准正交性量化**：D=10,000 维二进制随机向量，两两 Hamming 距离集中在
  0.5 ± 0.005（二项分布，STD 50 bits）；距离 < 0.476 的点不到百万分之一
  ——"任何两个随机向量几乎正交"的严格表述。
* **操作**：bundle（叠加，集合/多重集表示）、bind（逐维乘，成对绑定）、
  permute（置换，顺序表示）；相似度 = 距离（Hamming/Euclidean）。
* **容量**：随机向量互不相似，可区分海量实体。

**对 pgr**：bit 路径（±1 超立方体）正是 Kanerva 的表示；准正交性是
§2.1/§2.5 的理论源头（为什么点积能分离共享/非共享信号）。Kanerva 未
提供集合相似度的估计理论（那是 DotHash 补的）。

#### HyperGen（Bioinformatics 2024）—— 已在 §1–3 详细分析

要点回顾：FracMinHash 采样 + **稠密** ±1 编码（i16，WyRng）+ Jaccard→ANI；
i16 值域 [−N,N] 溢出风险；无 magic/版本。**与 pgr 稀疏无关**（它是
稠密），但 pgr 稠密 bit 与其同源（实现不同：RapidRng 跳步 vs WyRng、
i32 vs i16）。

#### ProbMinHash（IEEE TKDE 2020/2022）—— 完整精读

**核心**：概率 Jaccard JP（带权集合的 Jaccard 推广）+ P-MinHash 无偏
估计（方差 JP(1−JP)/k）+ ProbMinHash 一类 one-pass LSH（4 种算法：
2 种统计等价、2 种引入统计依赖降误差），可特化到普通 Jaccard 且超越
minwise hashing。

**对 pgr**：若未来支持 k-mer 多重度（拷贝数/覆盖度），JP 是现成理论
（对应 hv.md §4.4 的 weighted bundle 方向）；其"降低估计方差"的技巧
（引入统计依赖）对 pgr 的 D/s 选择有参考意义。

#### Order Min Hash（Bioinformatics 2019, Marçais）—— 完整精读

**核心**：edit distance 的 LSH（OMH），minHash 的改进——不仅看 k-mer
内容还看**相对顺序**，有 gapped LSH 理论保证；现有做法用 Jaccard/Hamming
代理 edit distance。

**对 pgr**：若距离方向需要 edit-distance 感知（而非 k-mer 集合），OMH
可参考；但 pgr 目前距离语义是"精确 k-mer 集合重叠"（Jaccard → ANI），
与 OMH 的 edit distance 目标不同——需明确场景再引入。

### 6.2 第二批：采样理论

#### Minimizer Jaccard 有偏（Bioinformatics 2022）—— 完整精读

**核心**：严格证明 minimizer sketch 的 Jaccard 估计**有偏且不一致**——
偏差不为零，且不随序列长度增长消失（估计收敛到与真实 J 不同的值）；
**给出偏差的解析公式**（作为共享 k-mer 沿序列布局的函数）；存在偏差
很大的序列族（真实 Jaccard 可比估计大 2 倍以上）；实证影响 mashmap
的映射准确性。

**对 pgr**：§2.6 弃用 minimizer 的直接依据；关键对比——**Mash 的
minhash（随机采样）估计无偏（Broder 1997），minimizer（窗口内最小）
有偏**；pgr 的 FracMinHash 属随机采样类（无偏类），syncmer 属局部
选择类（需 Shibuya 2022 另证无偏）。minmer（Kille）是修正方案。

#### Syncmers（PeerJ 2021, Edgar）—— 完整精读

**核心**：syncmer 家族定义——通过 k-mer 内**最小 s-mer 的位置**选 k-mer；
closed syncmer = 最小 s-mer 在 k-mer 首或尾。**同步性**：syncmer 由序列
本身识别（不依赖上下文），minimizer 会被侧翼突变删除而 syncmer 不会。
实验：同步实现**更低密度 + 更高 conservation**（相对 minimap2/Kraken 参数
的 minimizer）。

**对 pgr**：pgr closed syncmer 采样（§2.6、`libs/syncmer.rs`）的直接来源。

#### FracMinHash（bioRxiv 2022, Irber）—— 完整精读

**核心**：FracMinHash = modulo hash 的变体，支持**不同大小集合**的
Jaccard **containment** 估计（MinHash 对大小悬殊集合有偏）；sourmash
实现；70 万微生物参考基因组的规模验证；最小 metagenome cover（贪心集合
覆盖）。

**对 pgr**：§2.6 数值 ANI 推荐采样器的实现参考（sourmash 的工程实践）；
containment 语义对 pgr 的"参考 vs 查询"场景（4 万 cohort）有价值。

#### Syncmer ≡ FracMinHash（bioRxiv 2023, Koslicki）—— 完整精读

**核心**：**open syncmer sketch 与 FracMinHash sketch 在 k-mer 相似度
意义上等价**（注意：是 open 不是 closed）；open syncmer 有更好的距离
分布和基因组覆盖；k-mer truncation 可扩展到 open syncmer（多分辨率
估计 + 灵活种子）。

**对 pgr**：§2.6"open syncmer 与 FracMinHash 等价"的出处（pgr 用 closed，
需注意差异——closed 的无偏性由 Shibuya 2022 另证）；truncation 多分辨率
思路可参考。

#### Minmers（Bioinformatics 2023）—— 完整精读

**核心**：minmer = minimizer 的泛化——每窗口用 **rolling minhash 采样
多个 k-mer**（作者：Kille, Garrison, Treangen, **Phillippy**——MashMap
作者之一）；**理论 + 实证证明无偏** local Jaccard 估计（完整标题
"Minmers are a generalization of minimizers that enable unbiased local
Jaccard estimation"）；MashMap3 默认 ANI 阈值下比 minimizer 版快 **10×**。

**对 pgr**：§5.5 方向 3（minmer 替代 `seq_mins` 的 minimizer）的完整依据；
消除 §2.6 minimizer 偏差的现成方案。

#### Local k-mer selection 理论（Bioinformatics 2022）—— 完整精读

**核心**：local k-mer selection 的形式化（q-local 方法，Theorem 1：
共享 k+q−1 长区域的两序列，局部选择的 k-mer 互保）；conservation 的
**精确表达式**（Theorem 3）；(open/closed) syncmer、(a,b,n)-words 的
**闭式解**、minimizer 上界；**open syncmer 最优参数定理（Theorem 8）**；
os-minimap2 实证：更 conserved 的方法提升映射 reads **8.2%**，但更
conserved 的 k-mer 更重复 → 运行时增加（速度-质量权衡）。

**对 pgr**：采样器选择的定量理论框架（§2.6 未引用，§6 审计指出的缺口）；
"conservation 与重复性权衡"对 syncmer 参数选择有直接指导。

#### mod-minimizer（WABI 2024）—— 完整精读

**核心**：窗口保证 + 密度的形式化（密度下界 1/w）；random minimizer 密度
接近下界 2 倍；**mod-sampling** 两步骤采样（找最小 t-mer 位置 p，采
p mod w 处 k-mer）；**mod-minimizer（t ≡ k mod w）在 k→∞ 达到最优
密度**，且与 random minimizer 一样快。

**对 pgr**：长 k-mer 场景的采样器候选（§5.5 方向 7）；密度最优性对
大规模索引（4 万 cohort）有存储/查询收益。

#### Strobemers（Genome Res 2021, Sahlin）—— 完整精读

**核心**：2+ 个链式短 k-mer 的组合（哈希决定），替代单一 k-mer；对
突变率更不敏感、匹配分布更均匀、覆盖更高（vs k-mer 和 spaced k-mer）；
StrobeMap 验证聚类/比对场景。

**对 pgr**：采样层抗突变候选（§5.5 方向 7）；与 syncmer 互补（strobemer
解决 indel 敏感，syncmer 解决上下文依赖）。

#### Minimizer sketching 综述（Genome Biol 2024）—— 完整精读

**核心**：minimizer 入门 + 方法进展 + 应用全景（组装/宏基因组/比对/
纠错/泛基因组）+ 替代技术（universal hitting sets、syncmers、
strobemers）。

**对 pgr**：§2.6 采样方法的背景总览；UHS（universal hitting sets）是
未评估的候选。

#### Sparse and skew hashing of K-mers（Bioinformatics 2022）—— 完整精读

**核心**：压缩关联 k-mer 字典（MPHF 分配唯一 ID），利用 **minimizer
统计性质**压缩，支持数十亿 k-mer 的成员查询。

**对 pgr**：k-mer 集合存储方向（§5.1 Pibiri 条目）；若 4 万 cohort 的
k-mer 索引需要压缩存储可参考。

### 6.3 第三批：大规模聚类与搜索

#### RabbitTClust（Genome Biol 2023）—— 完整精读

**核心**：sketch-based 距离估计（Mash 距离）+ **两条聚类管线**：
**clust-mst**（MST 单链层次聚类，动态生成/合并部分聚类、不存全距离
矩阵 → 线性空间）与 clust-greedy；流式/并行；
**113,674 个完整细菌基因组（455 GB）6 分钟内聚类**，1,009,738 个
GenBank 细菌基因组（4.0 TB）34 分钟（128 核）；MinHash sketching +
**最小生成树（MST）** + 冗余检测（发现 1269 个完全相同基因组）；
距离阈值 0.05、内存 10.7 GB（bact-RefSeq）。

**对 pgr**：4 万 E. coli cohort 聚类的**直接对标**（§5.2 条目）——其
"sketch 距离 → 降维 → 图聚类（MST）"管线是 pgr 可借鉴的端到端方案；
冗余检测（完全相同基因组）与 pgr 的 4 万 cohort 去冗余需求吻合。

#### GSearch（NAR 2024）—— 完整精读

**核心**：k-mer hashing 概率数据结构（**ProbMinHash / SuperMinHash /
Densified MinHash / SetSketch**）估计基因组距离 + **HNSW 图搜索**；
O(log N) 复杂度，可扩展数十亿基因组（数据库分片策略）；8000 查询 vs
318k/3M 基因组几分钟、~6 GB（SetSketch 2.5 GB）；三阶段搜索策略
（按查询新颖度）。**方法细节**：6 种 MinHash 类算法可选——Densified
MinHash（最快，单哈希函数）、**ProbMinHash（默认：共享 k-mer 按丰度
加权、按总 k-mer 数归一化）**、SuperMinHash（Jaccard 精度优化）、
SetSketch（最低内存 2 变体）；"MinHash/SetSketch 距离 + HNSW 近邻"
组合首次用于基因组搜索。
**工程细节**：Rust 重写 hnswlib（hnswlib-rs 0.1.19，含内存映射）与
probminhash 0.1.10（11 种 MinHash 类算法，GSearch 用其中 6 种）；
三模块 tohnsw（建图）/ add（增量添加）/ request（查询），均并行；
建图 O(N log N)、查询 O(log N)；M 上限 255（--nbng，24–64 常用）、
ef_construct 默认 400（论文建议 >1000 提升召回）；
数据库分片（各片独立建图，汇总 top-K，每片取 ≥K 时与整体等价）；
三阶段搜索（nt → 全蛋白组 AAI → 通用基因 AAI，阈值 ~78% ANI /
52% AAI）；论文参数细菌 s=12,000、k=16，真菌 s=48,000、k=21，
蛋白 k=7；基准：318k RefSeq 建图 4.1 h（24 线程）、8466 查询 9.33
min；1M 基因组文件 ~118 GB（ProbMinHash）/ 9.8 GB（SetSketch）。
**源码核对（gsearch-master / hnsw_rs 0.3.4 / probminhash 0.1.12，
2026-08-08）**：① 图距离实现为 **sketch 签名的 Hamming 距离**
（`Hnsw::<Sig, DistHamming>`，anndists 实现于 u8/u16/u32/i32/i16），
并非论文文字里的 1−Jp 直接做图导航——Jaccard 语义由 LSH sketch 本身
保证，HNSW 只按签名相似度近邻；② 建图参数：ef_construct 默认
**400**、max_nb_conn = min(--nbng, 255)、capacity 硬编码 1,500,000、
层数上限 16（GTDB 实际 3 层）；③ **HubNSW 单层化选项**：
`--scale_modify_f ∈ [0.2,1]`（默认 1.0）可把多层图压成单层 NSW，
README 明确称对高维数据集 "better space requirement and accuracy"
（arXiv 2412.01940）——与 pgr 4096 维 HV 的召回问题直接相关；
④ k-mer 预处理：过滤非 ACGT、2-bit 编码、`rc().min()` 取正则链后
再哈希；⑤ 处理粒度：默认**逐序列 sketch**，`--block` 才把整文件连成
一条；基因组级排序用"序列级距离的乘积"聚合（matcher.rs
compute_merit_wl，阈值 0.99，注释自承 TODO 未调优）；⑥ 请求侧
ef_search 独立传参，输出阈值 0.99；⑦ README 快速上手参数：
k=21、s=18,000、n=128、ef=1600、--algo optdens、scale 0.25。
**召回（Table 3，真值 = 暴力 BLAST-ANI/AAI 的 top-K）**：定义
recall = |R'∩R|/|R|（R' = GSearch 返回的 top-K），跨查询平均，评估
top-5 / top-10。原核（对 318k RefSeq）：近缘查询（>78% ANI）nt 图
top-5/top-10 = 98.3%/96.2%（Tara）、97.7%/95.1%（Ye）；无近缘查询
nt 图 top-10 跌至 43–49%，换蛋白组图恢复到 95–97%，深分支（通用基因
图）94–96%；数据库本身无 ≥52% AAI 近缘的深分支只有 50–56%（数据稀疏，
非检索失败）。真菌（s=48,000，k=21）：nt top-10 99.4%（对 MUMMER-ANI）；
aa 图 top-5 99.7%、top-10 98.5%。病毒（~3M IMG/VR4，只建 aa 图）：
top-5 recall 98.32%。基因组完整度 >50% 时 top-10 recall >80%，低于
50% 不推荐。**注意**：这是"sketch 距离误差 + 图检索误差"的端到端召回，
真值是 ANI/AAI；比我们 HV 实验的图检索召回（真值 = HV 距离精确 top-10）
更严格，二者不可直接比大小。

**对 pgr**：pgr dist 的 **ANN 搜索方向（§5.5 方向 6）**的完整工程参考——
特别是"概率草图距离 + HNSW"的组合，以及 SetSketch 的低内存变体。
另注意：GSearch 的 sketch（s=12,000 寄存器，u16–u64 计 24–96
KB/基因组）与我们的 HV（16 KB）体积**同量级**，其优势在 HNSW 图查询
O(log N) 而非 sketch 更小；单层化 HubNSW 可作 4096 维召回实验的下一个
对照（见 genome-nn-query §6.4）。

#### HNSW（IEEE TPAMI 2020）—— 完整精读

**核心**：分层可导航小世界图（Hierarchical NSW）：多层近邻图，元素层数
按指数衰减概率随机分配；从顶层开始搜索 + 尺度分离 → **对数复杂度**；
邻居选择启发式提升高召回与聚集数据性能；skip-list 类似结构便于分布式。

**对 pgr**：任何"距离矩阵之外的图式近邻搜索"的经典底座（GSearch 用它，
pgr 若做 dist 的 ANN 也用它）；与 pgr 的 `dist` O(D) 逐对比较互补。

#### BIGSI（Nat Biotechnol 2019）—— 完整精读

**核心**：Bitsliced Genomic Signature Index——把 web 搜索的位切片方法
用于微生物基因组；索引**全部 447,833 个已公开细菌/病毒 WGS 数据集**，
存储比先前方法少 4 个数量级；增量可扩展至百万数据集；应用：耐药基因
（MCR-1/2/3）快速查找、质粒宿主范围、抗生素耐药量化。**方法细节**：
每个数据集一个固定长度二进制"签名"（bitsliced）；搜索词 = k-mer
（SNP/等位基因），"文档" = 原始读段/组装；位切片曾被 Zobel 1998
判定不如倒排索引（自然语言），但 Bing 2017 复活——微生物 DNA 维度
远高于英文（10^6 文档、10^10 唯一 k-mer），位切片伸缩性更好。

**对 pgr**：大规模索引搜索的里程碑（§5.2）；位切片（bitsliced）索引
思路对 pgr 的 .paf.idx / .pgi 索引设计有参考价值。

#### LexicMap（Nat Biotechnol 2025）—— 完整精读

**核心**：探针 k-mer（probe k-mers）子集高效采样数据库——保证每个
250 bp 窗口含多个种子、与探针共享前缀；层次索引低内存比对；**百万级
原核基因组**查询（基因/质粒/长读 >250 bp），数分钟内完成，精度与
SOTA 相当、速度更快内存更低。

**对 pgr**：大规模 pairwise 比对的工程参考（§5.2）；"探针 k-mer 采样 +
层次索引"对 pgr 的 pgi/lastz 混合管线（align fill/rest）有启发。

#### Linclust / MMseqs2 聚类（Nat Commun 2018）—— 完整精读

**核心**：**第一个运行时线性于 N（独立于聚类数 K）**的序列聚类算法；
关键技巧：每序列选 **m 个最低哈希值 k-mer** → 排序找共享 k-mer 的
序列组（桶）→ 每组选**最长序列为中心** → 每序列只与共享 k-mer 的
中心比较（三阶段逐步更慢更敏感）；分桶把 O(NK) 降到 O(N)；1.6 亿
宏基因组序列片段 10 小时聚类（50% identity），比之前快 **>1000×**；
可聚类超出内存的数据集。

**对 pgr**：蛋白侧聚类（若 pgr 扩展蛋白距离）的算法模板（§5.2）；
"分桶降复杂度"思想对 4 万 cohort 的聚类管线有直接借鉴。

#### annembed（NAR Genom Bioinform 2024）—— 完整精读

**核心**：改进 UMAP 类降维：t-SNE+UMAP 结合 + **HNSW 替换 K-NNG
限速步骤** + **MinHash LSH 做序列距离估计**；Rust 实现、全并行；
扩展功能：局部本征维度、hubness 计算；应用：微生物数据库可视化、
单细胞、宏基因组 contig binning。

**对 pgr**：embedding + ANN 结合方向（§5.2）——特别是"MinHash 距离 +
HNSW + 降维"的完整链路，与 pgr dist 的可视化/聚类下游相关。

#### MetaProFi（Bioinformatics 2023）—— 完整精读

**核心**：**第一个蛋白级 Bloom filter 索引**——氨基酸序列索引 + 氨基酸/
核苷酸双查询；共享内存、chunked 存储、高效压缩。

**对 pgr**：大规模集合存储查询方向（§5.2）；蛋白级索引若 pgr 扩展
蛋白距离（dist seq protein）可参考。

#### Phylogenetic profiling with MinHash（PLoS Comput Biol 2020）—— 完整精读

**核心**：MinHash 用于系统发育谱（phylogenetic profiling）——把基因
的物种分布谱转成可扩展的相似度计算，发现真核有性生殖相关基因。

**对 pgr**：谱系/分布向量的 MinHash 应用场景（§5.1 背景补充）；与
pgr 的"集合相似度"框架同源（MinHash vs pgr HV 是两种草图路线）。

### 6.4 第四批：哈希 / 比对 / 图表示 / 蛋白 / 其他（弱相关，简洁笔记）

#### k-mer 哈希家族

* **ntHash（Bioinformatics 2016）**：DNA/RNA 的**递归哈希**
  （`H(kmer_i) = f(H(kmer_{i−1}), r[i+k−1], r[i−1])`，相邻 k-mer O(1)
  更新），比替代方案快一个数量级；pgr syncmer 的乘性滚动哈希的同类
  对照。
* **ntHash2（Bioinformatics 2022）**：**spaced seed** 递归哈希，比旧版
  快 2.1×、比朴素适配快 3.8×；改进长 k-mer 碰撞率与哈希分布均匀性
  （修改 canonical hashing）；spaced seed 采样方向（§5.3）。
* **aaHash（2023）**：氨基酸**递归哈希**（多级哈希）；蛋白距离方向
  （dist seq protein）的哈希基础。

#### 比对类

* **miniprot（2023）**：蛋白→基因组比对器（最新算法，替代 10+ 年老工具）。
* **Block Aligner（2023）**：SIMD 加速 SW-Gotoh DP 的新范式——DP 矩阵
  **块**（greedy shift/grow/shrink，自适应计算区域），处理复杂评分矩阵/
  大 gap；比先前快 5–10×、错误率 <3%（蛋白/长读）；**Rust 库**（支持
  global/local/X-drop）——对 pgr 的 alignment/banded 库有工程参考。
* **Exact global alignment with A\***（2024）：edit distance 的 A* 精确
  全局比对（线性时间近似）；与 pgr 的精确比对方向相关但非测距核心。
* **Strobealign（Sahlin 2022）**：**syncmers + strobemers 组合种子**
  （动态模糊种子），E-hits 指标，短读比对快数倍——pgr 采样层候选的
  组合用法示例。
* **ropebwt3（2024）**：terabase 级 BWT 构建/查询；泛基因组索引方向
  （pgr paf/图索引可参考，但与测距聚类弱相关）。

#### 图 / k-mer 集表示

* **Simplitigs（Genome Biol 2021）**：de Bruijn 图的紧凑表示（比
  unitigs 更短更少）；BWT 索引时内存/查询更优。
* **Matchtigs（Genome Biol 2023）**：k-mer 集**最小明文表示**
  （多项式算法 + 贪心，建模为 **minimum-cost flow / Chinese postman
  problem**），比 unitigs 压缩 59%、字符串数减 97%、SSHash-Lite 查询
  快 4.26×。
* **CBL（Conway–Bromage–Lyndon, 2024）**：压缩**动态**精确 k-mer 集
  表示（最小循环旋转/Lyndon 词压缩 + **Elias-Fano 风格动态位向量**），
  Rust 库，唯一支持 **in-place 集合操作**——k-mer 集存储（§5.1 Pibiri
  条目的补充）。
* **WebGraph（Fontana/Vigna/Zacchiroli 2024）**：Rust 大规模图压缩
  （BVGraph 系）；pgr 泛基因组图（petgraph）的存储参考，弱相关。

#### 蛋白 / 其他

* **Snekmer（2023）**：氨基酸重编码（AAR）→ k-mer 向量 → 蛋白分类；
  蛋白侧"k-mer 指纹"路线。
* **UniProt 去冗余（2016）**：蛋白库冗余消除方法学——与 pgr 4 万
  cohort 去冗余（ecoli-cohort.md 上游）类比。
* **Local homology & distance（Edgar 2004）**：**压缩氨基酸字母表**的
  线性时间同源识别与距离——Edgar 的早期距离工作，蛋白侧参考。
* **SECOM（2012）**：hash seed + community detection 的蛋白域识别
  （共享 hash seed 作边权）——"哈希共现 → 图聚类"思路。
* **Flexible protein database（2022）**：氨基酸 k-mer 蛋白数据库。
* **SubseqHash（2023）**：**子序列**作种子（非子串），高错误率有效；
  关键：**ABC order**（特定序）下最小子序列可多项式计算，碰撞概率
  接近 Jaccard；read mapping/比对/overlap 三场景碾压子串种子——采样层
  抗错方向（§5.5 方向 7 的补充候选）。
* **k-mer 方法综述（CSBJ 2024）**：k-mer 方法与应用全景
  （方向总览入口，§5.3）。
* **Manifold Learning 综述（Annu Rev Stat Appl 2024）**：非线性降维方法
  原理与统计基础——若 pgr 聚类下游需要可视化/嵌入（annembed 也基于
  此方向），可作背景；与测距本身弱相关。
* **Wheat pangenome（2024）**：k-mer 泛基因组在植物上的应用案例，
  相关性弱。
* **Domain-Based identity thresholds（2009）**、**Improved protein
  homolog（2023）**：蛋白域/同源检测方向，弱相关。
* **Statist Med 1996**：logistic regression 的 explained variation
  统计方法——与本方向无关（保留说明）。

### 6.5 稀疏投影的完整统计理论（2026-08-08，仿照 DotHash Theorem 2 推导）

> 回答"稀疏投影能否有 DotHash 级数学推导"：**能**。pgr 的"每元素 s 个
> 随机 ±1 桶"构造可以仿照 DotHash Theorem 2 的路径得到无偏性与方差
> 公式（推导 + 实验验证于 §6.5 下方测试）。关键是区分两个方差视角。

**构造**：元素 a 经 ψ: S → R^D 映射——独立均匀选 s 个维度、每维赋 ±1
（splitmix64 决定位置与符号）。集合 sketch = 元素向量之和。归一化
ψ'(a) = ψ(a)/√s。

**定理（无偏性）**：对任意集合 A、B（shared = |A∩B|），

```
E[a'·b'] = shared           （与 DotHash Theorem 2 相同）
```

证明：E[ψ'(x)·ψ'(y)] = 1 当 x=y（s 个 ±1/√s 的平方和 = s·(1/s) = 1），
0 当 x≠y（独立）；a'·b' = Σ_{x∈A,y∈B} ψ'(x)·ψ'(y)，取期望只剩 shared
个 x=y 项。✓

**定理（方差，投影随机性视角：固定 A、B，独立随机投影）**：

```
Var(a'·b') = shared/s + (|A||B| − 2·shared)/D
```

证明（仿 DotHash）：Var(a'·b') = Σ_i Var(Σ_{x,y} ψ'(x)_i ψ'(y)_i)。
相等对（shared 个）：Var(ψ'(x)_i²) = 1/(sD) − 1/D²（选 i 概率 s/D，
值 ±1/√s）；不等对：Var(ψ'(x)_i ψ'(y)_i) = E[ψ'_xi²]·E[ψ'_yi²] =
1/D²。Σ_i 后：shared·(1/s − 1/D) + (|A||B|−shared)·(1/D)。✓

**与 DotHash（稠密超立方体）对比**：DotHash 的 Var = (1/d)·[|A||B| +
|A∩B|² − 2|A∩B|]；pgr 稀疏多出 **shared/s 项**（s 越小方差越大），
且 |A∩B|² 项被 |A∩B| 项替代（稀疏投影不产生元素对之间的交叉碰撞
方差）。典型参数（N=3000、shared=500、D=4096）下 shared/s 项占比
≤18%（s=1），s≥3 时 <7%——**s 影响方差但影响有限**。

**第二视角（集合随机性：不同集合对、单次投影——pgr 实际评估方式）**：
MAE ≈ 常数/√D，与 s 无关（§2.7 的 50 组扫描实测 + 推导：相对方差
≈ n²/(shared²·D)，s 消去）。原因：pgr 的投影是**确定性**的（seed 来自
k-mer hash），跨集合评估时集合随机性主导，投影方差的 shared/s 项被
淹没。

**实验验证**（`test_hash_hv_sparse_projection_variance`，300 次独立
投影）：归一化均值 ≈ shared（无偏 ✓）；方差与预测公式吻合
（s=1: 2464 vs 2697、s=3: 2388 vs 2364、s=16: 2498 vs 2228、
s=64: 2330 vs 2205，s 大时第二项主导）。

**结论**：pgr 稀疏投影有完整的无偏性 + 方差理论（仿 DotHash），
**稀疏投影可以获得 DotHash 级推导**；s 的影响是"小方差项"而非零。
若未来需要更严的误差界，可直接套 Chebyshev：
`Pr(|a'·b' − shared| ≥ ε) ≤ Var/(ε²)`。

**文献家族确认（§6.6 外部补录后）**：pgr 的 s 桶 ±1 构造不是孤立方案——
属于稀疏随机投影 / 特征哈希 / count sketch 家族：Li 2006 的稀疏分布
（±√s、非零概率 1/s）与 pgr **期望等价**；Weinberger 2009 提供特征
哈希内积近似的指数尾界；Count-Min Sketch（每元素多桶）与 pgr 结构
同构；Achlioptas 2001 提供稀疏保距奠基。**稀疏投影的理论支撑 = 文献
家族（§6.6）+ 自推无偏/方差（§6.5）+ 实测（§2.7）三层**。

### 6.6 外部补录文献（2026-08-08，用户找到——稀疏随机投影家族）

> 这 6 篇是之前标注"二手/未引"的关键文献，用户找到放入目录后完整精读。
> **结论先行**：pgr 的"每元素 s 个随机 ±1 桶"不是孤立构造——它属于
> 稀疏随机投影 / 特征哈希 / count sketch 家族的**期望等价形式**，理论
> 支撑比之前判断的更充分。

#### Database-friendly random projections（PODS/JCSS 2001/2003, Achlioptas）—— 完整精读

**核心定理（Theorem 2）**：任意 n 个 d 维点可嵌入 k 维
（`k₀ = (4+2β)/(ε²/2 − ε³/3) · log n`），**嵌入矩阵元素可来自两个极简
分布**：

1. ±1 等概率（Rademacher）；
2. **稀疏分布：+1 概率 1/6、−1 概率 1/6、0 概率 2/3（缩放 √3）**。

保距保证：`(1−ε)‖u−v‖² ≤ ‖f(u)−f(v)‖² ≤ (1+ε)‖u−v‖²`，成功概率
≥ 1−n^−β。稀疏分布的实际计算 = 数据库聚合（扔掉 2/3 属性、剩余分两半、
各求和、取差）——**不需要浮点运算**。

**对 pgr**：稀疏随机投影的奠基定理——"随机 ±1 矩阵稀疏化（2/3 零）
仍保距"的直接证明。与 pgr 的关系：Achlioptas 是**数据点投影**矩阵
稀疏；pgr 是**元素嵌入**向量稀疏——结构不同但同属稀疏随机投影家族。

#### Very sparse random projections（KDD 2006, Li）—— 完整精读

**核心**：把 Achlioptas 的稀疏分布推广为**通用稀疏度参数 s**：

```
r_ji = √s · {+1 概率 1/(2s)、0 概率 1−1/s、−1 概率 1/(2s)}
```

即"以 1/s 概率非零、值 ±√s"（s=3 时还原 Achlioptas）。**关键洞察**：
保距只需"零均值 + 单位方差"（对称分布即可），稀疏化改变的是方差/平均
误差——稀疏随机投影 = **以 1/s 速率随机采样**（s=3 即 2/3 稀疏）。

**对 pgr（最重要的对应）**：Li 的稀疏分布（每位置独立以 1/s 概率非零，
期望每行 s 个非零）与 pgr 的"每元素恰好 s 个 ±1 桶"**期望等价**——这是
pgr 稀疏投影最直接的文献理论支撑。两者差异：Li 是非零数泊松式可变，
pgr 是固定 s 个（方差更小，理论更好）；Li 是数据投影、pgr 是元素嵌入。

#### Feature hashing（ICML 2009）—— 完整精读

**核心**：特征哈希（每特征哈希到桶 + 符号）做降维/降存储；**提供指数
尾界**（exponential tail bounds）——随机子空间交互可忽略（高概率）；
多任务学习（数十万任务）实证。

**对 pgr**：特征哈希通常每元素 1 桶（pgr 的 s=1 特例）；其指数尾界是
"元素→哈希桶 + 内积近似"的理论保证（pgr 稀疏 s=1 的直接对应）。

#### Count-Min Sketch（J. Algorithms 2005, Cormode）—— 完整精读

**核心**：数据流汇总的**子线性空间结构**——每元素哈希到 d 个桶（每行
1 个），点/范围/**内积**查询近似，误差 ε 概率 δ；从 1/ε² 改进到 1/ε。

**对 pgr（结构最接近）**：CMS 的"每元素哈希到多个桶、累加"与 pgr 的
"每元素 s 个 ±1 桶"**几乎同构**（CMS 用 +1 计数、pgr 用 ±1 投影）。
CMS 提供内积查询的 (ε, δ) 保证——pgr 稀疏投影的近似保证可参考 CMS
的分析框架（虽然 CMS 是频率估计、pgr 是相似度，但桶结构相同）。

#### SparseHD（IEEE ICCD 2019）—— 完整精读

**核心**：稀疏超维计算（HDC 分类）：二进制（1-bit）超向量精度不足
（>50% 损失）、多 bit 提高精度但牺牲能效；**SparseHD 把训练好的 HD
模型稀疏化（最多 90% 稀疏）+ 迭代重训练补偿，质量损失 <1%**；FPGA
加速器利用稀疏性：48.5× 低能耗、15.0× 快（vs GPU）。

**对 pgr**：稀疏 HD 方向的工作——但注意机制不同：SparseHD 是**训练后
模型稀疏化**（裁剪超向量元素），pgr 是**编码阶段稀疏**（每元素只碰
s 维）；两者都证明"稀疏 + 低质量损失"可行。SparseHD 的"1-bit 精度不足"
对 pgr 的 bit vs i8 决策有参考（i8 保语义变体的存在意义）。

#### An Improved Data Stream Summary（重复版）

与 Journal of Algorithms 2005 版同一文献（Count-Min Sketch），已覆盖。

### 6.7 Genome Res 2023（FracMinHash 校正与 ANI 置信区间，用户补充 PDF）

> Genome Res 2023, 33(7):1061–1068, doi:10.1101/gr.277651.123
> （Hera, Pierce-Ward & Koslicki）。完整精读，是 [[../design/hv.md]] §2.6
> "FracMinHash 校正 + CI"引用的原始出处。本地 PDF：
> `~/sync/zotero/bacteria/clustering/Genome Res - 2023 - Deriving confidence intervals for mutati.pdf`。

**核心结论**：**FracMinHash 不是无偏的**（严格意义上），但偏差可校正；
校正后是 containment 的无偏估计，且有渐近正态性 → 可推 ANI 置信区间。

**Theorem 1（containment 期望，节选自 Irber 2022）**：

```
E[Ĉ_frac(A,B) | |FRAC_s(A)| > 0] = C(A,B) · (1 − (1−s)^|A|)
```

偏差因子 `(1 − (1−s)^|A|)`——|A| 小（病毒/短序列）时明显，|A| 大
（细菌 ~5M）时 ≈ 1。

**去偏 fractional containment**（公式 3）：

```
C_frac(A,B) = |FRAC_s(A)∩FRAC_s(B)| / (|FRAC_s(A)|·(1 − (1−s)^|A|))
```

需要 |A|（全集合大小，可估计，如 HyperLogLog）。

**对 pgr**：pgr `dist seq --sampler frachash` 的 containment =
inter/card1（未校正）；实测 e2348×cft073 偏低 ~5%（0.60 vs 真值
0.633）——**该差异不能用 Hera 偏差因子解释**（细菌规模因子 ≈ 1），
疑似实现细节（N 处理 / canonical 边界）或另有系统因素，实现校正时
需逐项核对 pgr 与 Hera 定义。ANI 点估计 + CI 的推导见论文 Methods，
对应 pgr `--ci`（当前为正态近似，Hera 更精确）。

**todo 关联**：`notes/todo.md` §2 "FracMinHash containment/ANI 偏差
校正（等 Hera 论文）"——论文已到位（2026-08-08 用户下载），可开始。
