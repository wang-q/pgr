# 百万级基因组最近邻查询：方法与工程调研

> 目的：梳理文献中已实现"百万级基因组最近邻查询/聚类/搜索"的工作，
> 归纳技术路线与用户交互形态，为 pgr 后续设计同类功能的用户界面做
> 前置调研。来源：`references/hv.md` §5.2 / §6.3 逐篇精读笔记，
> 细节与完整阅读记录见该文档。

## 1. 工作全景

| 工作 | 出处 | 功能 | 规模 | 方法核心 |
|---|---|---|---|---|
| **GSearch** | NAR 2024 | 最近邻查询 | 318k/3M 基因组，查询几分钟 | ProbMinHash/SetSketch 距离 + HNSW 图 |
| **RabbitTClust** | Genome Biol 2023 | 聚类 / 去冗余 | 113k 细菌 6 分钟；1M GenBank 34 分钟（128 核） | Mash 距离 + MST / greedy |
| **BIGSI** | Nat Biotechnol 2019 | 索引搜索（存在性/定位） | 447k 细菌/病毒 WGS 数据集 | bitsliced 二进制签名索引 |
| **LexicMap** | Nat Biotechnol 2025 | 序列比对查询 | 百万级原核基因组，数分钟 | 探针 k-mer 采样 + 层次索引 |
| **MMseqs2 (Linclust)** | Nat Commun 2018 | 蛋白序列聚类 | 1.6 亿序列 10 小时 | 低哈希 k-mer 分桶降复杂度 |
| **HNSW** | IEEE TPAMI 2020 | ANN 搜索底座 | —（算法） | 分层可导航小世界图，O(log N) |
| **annembed** | NAR Genom Bioinform 2024 | 可视化 / 嵌入 | 微生物数据库 | MinHash 距离 + HNSW + 降维 |
| **MetaProFi** | Bioinformatics 2023 | 蛋白索引查询 | — | 蛋白级 Bloom filter 索引 |
| Mash（背景） | Genome Res 2016 | 全 RefSeq pair 距离 | <30 CPU h | bottom-k MinHash |
| sourmash（背景） | bioRxiv 2022 | 大规模宏基因组分析 | — | FracMinHash |

## 2. 逐工作总结

### GSearch（NAR 2024）—— 最近邻查询

- **方法**：k-mer 概率哈希（6 种可选：Densified MinHash / **ProbMinHash
  （默认，共享 k-mer 按丰度加权）** / SuperMinHash / SetSketch 低内存变体）
  估计基因组距离，**HNSW 图做近邻搜索**——"MinHash 距离 + HNSW 近邻"
  组合首次用于基因组搜索；O(log N)，数据库可分片扩展至数十亿。
- **规模**：8000 查询 vs 318k / 3M 基因组，数分钟，内存 ~6 GB（SetSketch
  2.5 GB）。
- **交互形态**：给定查询基因组 → 返回按距离排序的近邻列表（含可调
  搜索策略，按查询新颖度分三阶段）。
- **技术路线（论文 Methods）**：① sketch：ProbMinHash3a（默认）/
  Densified MinHash / SuperMinHash / SetSketch，距离 = 1−Jp（加权）
  或 1−J；**距离只按需计算，不做 all-vs-all**，HNSW 图节点即各基因组
  sketch；② 建图 O(N log N)：逐点贪心插入 + 插入后反向更新邻居列表；
  M（节点最大度数）20–200、ef_construct >1000（GTDB 65,703 基因组
  实际 3 层：65,180/519/4）；③ 查询 O(log N)：从顶层贪心下探到 0 层；
  ④ 实现：Rust 重写 hnswlib（hnswlib-rs 0.1.19）+ probminhash 0.1.10，
  三模块 tohnsw / add / request，支持增量添加新基因组与内存映射；
  ⑤ 扩展：数据库随机分片、各片独立建图，汇总各片 top-K 排序
  （每片取 ≥K 时与整体 top-K 等价）；⑥ 三阶段搜索：nt 图 → 全蛋白组
  AAI → 通用基因 AAI（切换阈值 ~78% ANI / 52% AAI），弥补 k-mer
  Jaccard 在远缘下失准；论文参数：细菌 s=12,000、k=16，真菌
  s=48,000、k=21，蛋白 k=7。
- **基准**：GTDB v207（65,703 基因组）建图 1.3 h（24 线程）/
  0.27 h（128 线程），数据库文件 3.0 GB；318k RefSeq（2 TB 数据）
  建图 4.1 h（24 线程）、文件 15 GB，8466 查询 9.33 min / ~16 GB
  （对比 Dashing v1/v2 与 BinDash 21–42 min）；查询复杂度 O(log N)
  实测随 N 增长符合 log 拟合；1M 基因组数据库文件 ~118 GB
  （ProbMinHash）/ 9.8 GB（SetSketch）。
- **局限**：病毒基因组 nt 图过稀易陷局部最优 → 官方推荐 aa 级；
  4 种 sketch 均要求 LSH/度量性质（非度量的 FastANI/FracMinHash
  距离被明确排除）；未来方向 = spaced k-mer / tensor sketch、
  GPU、O(N·c) 建图。

### RabbitTClust（Genome Biol 2023）—— 聚类 / 去冗余

- **方法**：sketch 距离（Mash）+ 两条聚类管线：**clust-mst**（最小生成树
  单链层次聚类，动态合并、不存全距离矩阵 → 线性空间）与 clust-greedy；
  流式/并行；附带冗余检测（识别完全相同基因组）。
- **规模**：113,674 个完整细菌基因组（455 GB）6 分钟；1,009,738 个
  GenBank 基因组（4.0 TB）34 分钟（128 核），内存 10.7 GB（阈值 0.05）。
- **交互形态**：输入基因组集合 → 输出聚类划分 + 冗余组。

### BIGSI（Nat Biotechnol 2019）—— 索引搜索

- **方法**：把 web 搜索的**位切片（bitsliced）**方法用于微生物基因组；
  每个数据集一个固定长度二进制签名，搜索词 = k-mer（SNP/等位基因）。
  微生物 DNA 维度（10⁶ 文档、10¹⁰ 唯一 k-mer）下位切片伸缩性优于
  倒排索引。
- **规模**：索引全部 447,833 个已公开细菌/病毒 WGS 数据集，存储比先前
  方法少 4 个数量级；增量可扩展至百万。
- **交互形态**：给定 k-mer / 序列特征 → 返回包含它的数据集列表（存在性
  查询，如耐药基因 MCR-1/2/3 的宿主范围）。

### LexicMap（Nat Biotechnol 2025）—— 序列比对查询

- **方法**：**探针 k-mer（probe k-mers）**子集高效采样——保证每个
  250 bp 窗口含多个种子、与探针共享前缀；层次索引低内存比对。
- **规模**：百万级原核基因组查询（基因/质粒/长读 >250 bp），数分钟，
  精度与 SOTA 相当、速度更快内存更低。
- **交互形态**：给定查询序列 → 返回数据库中的比对命中（近似 pairwise
  比对，而非仅距离）。

### MMseqs2 / Linclust（Nat Commun 2018）—— 蛋白序列聚类

- **方法**：每序列选 m 个最低哈希 k-mer → 排序分桶（共享 k-mer 的序列
  成组）→ 组内选最长序列为中心 → 每序列只与共享 k-mer 的中心比较
  （三阶段逐步更慢更敏感）。**分桶把 O(NK) 降到 O(N)**——第一个运行时
  线性于 N 的序列聚类算法。
- **规模**：1.6 亿宏基因组序列片段 10 小时聚类（50% identity），比之前
  快 >1000×；可处理超内存数据集。
- **交互形态**：输入序列集合 → 输出聚类（Linclust 为 MMseqs2 的
  easy-cluster 管线，命令行式）。

### HNSW（IEEE TPAMI 2020）—— ANN 搜索底座

- **方法**：分层可导航小世界图——多层近邻图，元素层数按指数衰减概率
  分配；顶层搜索 + 尺度分离 → **O(log N)**；邻居选择启发式提升召回。
- **意义**：GSearch / annembed 的底层（也是 pgr 若做 ANN 的候选底座）；
  与 pgr `dist` 的 O(D) 逐对比较互补。

### annembed（NAR Genom Bioinform 2024）—— 可视化 / 嵌入

- **方法**：t-SNE+UMAP 结合 + **HNSW 替换 K-NNG 限速步骤** + **MinHash
  LSH 做序列距离**；Rust 实现、全并行；附带局部本征维度、hubness 计算。
- **意义**："MinHash 距离 + HNSW + 降维"完整链路，与 pgr dist 的可视化/
  聚类下游相关。

### MetaProFi（Bioinformatics 2023）—— 蛋白索引

- **方法**：蛋白级 Bloom filter 索引（氨基酸序列 + 氨基酸/核苷酸双查询），
  共享内存、chunked 存储。
- **意义**：蛋白距离扩展时的存储查询参考。

## 3. 技术路线归纳

百万级场景的共同前提：**先 sketch（降维）再索引/聚类**，避免全量
pairwise 比较。归纳为三条路线：

1. **sketch 距离 + ANN 图搜索**（GSearch / HNSW / annembed）：适合
   "给定一个查询基因组，返回 K 个最近邻"的**查询型**交互；图索引
   增量可更新，O(log N)。
2. **sketch 距离 + 图聚类**（RabbitTClust / MMseqs2）：适合"给一组
   基因组，划分聚类 / 去冗余"的**批处理型**交互；MST 单链或分桶
   降复杂度。
3. **索引搜索**（BIGSI / LexicMap / MetaProFi）：适合"某个 k-mer /
   序列特征出现在哪些基因组"的**定位型**交互（存在性 / 比对命中）；
   bitsliced / 探针 k-mer / Bloom filter 各有取舍。

## 4. 用户交互形态（UI 设计参考）

| 工作 | 输入 | 输出 | 交互类型 |
|---|---|---|---|
| GSearch | 查询基因组 | 排序近邻列表 | 查询（在线） |
| RabbitTClust | 基因组集合 | 聚类 / 冗余组 | 批处理（离线） |
| BIGSI | k-mer / 序列特征 | 含该特征的数据集 | 定位（在线） |
| LexicMap | 查询序列 | 比对命中 | 查询（在线） |
| MMseqs2 | 序列集合 | 聚类 | 批处理（离线） |
| annembed | 序列集合 | 嵌入 / 可视化坐标 | 批处理（离线） |

## 5. 场景：物种内聚类选参考 + PBit 归档（核心用例）

> 设计目标是和 pgr 其他组件配合，端到端跑通一个物种的压缩归档流程。
> 本场景串起 `dist`（聚类/选参考）→ `align pgi`（双序列比对）→
> `pbit`（归档），其中**聚类、构树、剪枝由独立项目 Necom 承担**
> （`~/Scripts/necom`），pgr 不重复实现。

### 5.1 工作流

1. **输入**：某物种已有的一系列基因组（FASTA 集合）。
2. **距离矩阵**：对全部基因组两两 sketch 距离（`dist` 家族：
   mini/mash/frac），输出 pair TSV。
3. **聚类 / 构树 / 剪枝（Necom）**：pair TSV →
   `necom mat to-phylip`（转 PHYLIP）→ `necom clust`（hier / dbscan /
   k-medoids / mcl / nj / upgma）聚类或构树 → `necom cut`（按
   height / K / clade size / dynamic 剪枝成扁平分区）→ 输出分组。
4. **挑选 reference**：每组选代表基因组（组内中心 / 最长 / 最全），
   成为 pbit 的多参考集合。
5. **路由 + 比对**：其余 sample 按组归属路由到对应 reference，对
   reference–sample 对做双序列比对（`align pgi` → PAF）。
6. **归档**：reference + sample + PAF 交给 `pgr pbit create`
   （多参考 `-r` 可重复、PAF 经 TSV 第 3 列/`--paf` 传入，CIGAR 编码
   增强压缩）→ 产出 `.pbit` 归档。

### 5.2 现有组件对接点

| 环节 | 组件 | 现状 |
|---|---|---|
| 距离矩阵 | `dist mini/mash/frac` | 已有：草图距离，输出 pair TSV（name1 name2 dist） |
| 聚类 / 构树 / 剪枝 | **Necom**（`~/Scripts/necom`） | 已有：`mat to-phylip` / `clust hier/dbscan/k-medoids/mcl/nj/upgma` / `cut simple/dynamic/hybrid` / `nwk` / `eval`；**pgr 不重复实现聚类** |
| 比对 | `align pgi` | 已有：双序列比对输出 PAF |
| 归档 | `pgr pbit create` | 已有：多参考 `-r` + 样本 `-i` + PAF 列，CIGAR/LZ-diff 编码 |
| 选参考 | Necom 分组 + 手工/脚本 | 组内选代表（中心/最长/最全），可由管线命令生成 |
| 自动路由 | `pbit` | **缺口**：路由手动（TSV 第 4 列，默认参考 0）；`design/pbit.md`
  决策点 1 明确"自动路由留待多样性 cohort 数据证明收益后再做"——本场景
  即触发条件 |

### 5.3 UI 候选形态（待后续讨论）

- **批处理管线**：一个命令串起全部步骤，如
  `pgr pl pbit-cluster <genomes...> -o out.pbit`（内部：`dist` 输出
  pair TSV → **Necom** 聚类/剪枝 → 选参考 → `align pgi` →
  `pbit create`）；
- **分步命令**：`pgr dist` 保证 pair TSV 可直接喂
  `necom mat to-phylip`（格式对齐），Necom 分组输出作为 `pbit create`
  的 TSV 路由来源（自动生成）；
- 渐进式：先对齐 `dist` 输出 ↔ Necom 输入格式，pbit 自动路由作为
  消费端增强。

### 5.4 待定问题

- sketch 选型：frac（无偏数值）vs mini（快）vs mash（兼容）做聚类；
- 聚类阈值 / 方法由 Necom 提供（`clust` 多算法 + `cut` 多准则），
  场景中选型待定；
- reference 挑选准则（中心 vs 最长 vs 最全）；
- 路由失败处理（离群 sample 不属任何簇时：新建参考 / 警告跳过）；
- 与 `.hv` 路径的分工（超大规模初筛用 hv，候选对再精确聚类）。
- `dist` pair TSV 与 Necom 输入格式的精确对齐（列序、缺失值、自比较
  行），必要时加 `--necom` 兼容输出。

## 6. HV 的存储与最近邻检索（SQLite 偏好）

> 2026-08-08 调研。问题：有没有现成数据库/工具可以直接做 HV 矢量
> 最近邻（用户偏好 SQLite 内置向量比较）。先纠正一个前提，再列方案。

### 6.1 HV 数据形态（前提纠正）

- HV = 固定 D 维稠密 `Vec<i32>`（D=4096）：`.hv` 文件 = 4096 × 4 B
  + 头部 ≈ **16 KB/基因组**（与 `pbit-decisions.md` 的 16 KB/样本一致）。
- "bit（±1）"是每维 ±1 的**编码方式**，不是打包位向量；稀疏路径
  s=1 只保证值稀疏，文件仍按稠密 i32 存储。
- 规模账：1 万 ≈ 160 MB；4 万 ≈ 640 MB；百万 ≈ 16 GB。
- 距离语义：inter = H_A·H_B/D（点积恢复共享 k-mer 数）；余弦 =
  inter/√(card₁·card₂)，与 Jaccard 单调 → 近邻**排序**可直接用
  cosine/dot；精确 Jaccard 可由 cosine + 范数（card = ‖H‖²/D）恢复。

### 6.2 SQLite 方案对比

| 方案 | 许可 | 向量类型 | 检索方式 | 备注 |
|---|---|---|---|---|
| **sqlite-vec**（asg017） | MIT / Apache-2.0 | float[N] / int8[N] / bit[N] | `vec0` 虚拟表 KNN（`match … order by distance limit k`），**精确扫描、无 ANN** | Rust crate `sqlite-vec`；metric 支持 l2/cosine/l1/hamming；k ≤ 4096；HV 存 float[4096] 与 i32 等体积（16 KB），HV 值域（√N 量级）在 f32 精确整数范围内 |
| **sqlite-vector**（rqlite） | **Elastic 2.0**（生产需商业授权） | Float32 / Float16 / BFloat16 / Int8 / UInt8 | 普通表 BLOB + 精确扫描，可量化 | 许可与 MIT 的 pgr 不兼容，**排除** |
| 自研：SQLite 存 BLOB + pgr SIMD 扫描 | MIT | i32 原样 | 线性扫描（`hv.rs` 已有 SIMD 点积） | 零新依赖；≤10 万规模与 sqlite-vec 等价 |

### 6.3 超大规模（>10 万）候选

- **usearch 深入评估（2026-08-08）**：
  * **本质**：单文件 C++11 HNSW（核心 ~3K SLOC）+ Rust `cxx` FFI
    绑定。`build.rs` 用 `cxx_build` 编译 C++17（需 g++/clang++/MSVC），
    默认 `-ffast-math -O3`；可选 OpenMP/numkong（SIMD）feature。
    Apache-2.0，支持 f32/f16/i8/二进制向量与 cosine/inner 等度量，
    O(log N)、索引可持久化。
  * **用户负担**：pgr 目前全 Rust 依赖（无 C++ 构建链）；加 usearch =
    源码构建必须带 C++17 工具链 + 编译时间显著变长；索引文件格式与
    行为绑定外部库，2.x 活跃迭代、pre-1.0 语义，升级有破坏性风险。
  * **掌控性**：跨 FFI 黑盒，图结构/内部状态不可见、难调试；cxx
    bridge 维护成本。
  * **技术面**：HNSW 在 4096 维的召回/延迟高度依赖实现——rust-cv
    `hnsw` 0.11 实测 N=30k 时 recall@10 上限仅 ~0.92（ef=400）、
    ef=10 仅 0.73（`benchmarks/bench-hv-ann-recall.md`）；同数据下
    GSearch 同源 `hnsw_rs` 0.3.4 可达 0.974–0.990，但查询慢 7–10×
    （`benchmarks/bench-hv-ann-hubnsw.md`），HubNSW 单层化（scale 0.2）
    只带来 0.4–1.8 pp 的微弱改善。GSearch 的 sketch（s=12,000 寄存器，
    u16–u64 计 24–96 KB）与我们的 16 KB HV **同量级**，其优势在 HNSW
    图查询 O(log N) 而非 sketch 更小（SetSketch 低内存变体才明显更省，
    ~2.5 GB/318k 基因组）。
  * **官方技术路线（README；无同行评议论文，引用为软件引用
    `@software{Vardanian_USearch}`，Zenodo DOI 10.5281/zenodo.7949416）**：
    算法即 HNSW（与 FAISS 相同，论文 = Malkov & Yashunin, IEEE TPAMI
    2020）；工程主打单文件 C++11（~3K SLOC）+ 零强制依赖 + 用户自定义
    metric（JIT）+ 多语言绑定；支持 f32→bf16/f16/float8/i8/单比特量化、
    uint40_t（4B+ 条目）、mmap 磁盘视图、多索引并行、谓词过滤、索引内
    聚类与语义 join。官方明确划分精确/近似：HNSW 用于**百万级**，小集合
    用 `exact=True` 的 SIMD 暴力扫描（README 示例 10k×1024d）。**立场
    差异**：README 反对"量化模型 + 降维"路线（"only sometimes
    reliable，会改变数据统计特性、分布漂移需重调"），主张高精度计算 +
    低精度存储；这与我们 >30k 先 PCA 降维的建议动机不同（我们是为绕开
    4096 维 HNSW 召回天花板），记录以免文档自相矛盾。
- 文献参照：GSearch = SetSketch/ProbMinHash + HNSW（查询型）；
  RabbitTClust = Mash + MST（批处理聚类型，§2）。

### 6.4 推荐与待定

- **ANN 召回实测结论（2026-08-08，`benchmarks/bench-hv-ann-recall.md`）**：
  * N=1k：HNSW 召回 0.98（ef≥10），比精确快 ~3×；
  * N=10k：ef=20 召回 0.90、0.46 ms/查询（精确 10 ms）；
  * N=30k：召回上限 0.92（ef=400 + ef_c=200），ef=10 仅 0.73；
    查询 0.65–1.5 ms vs 精确 30 ms（快 20–46×）；
  * 召回上限由查询 ef 决定，构建 ef_c 64→200 仅提升 0.6–1.2 pp——
    4096 维高维诅咒是主因，不是图质量（rust-cv 实现结论）。
- **hnsw_rs / HubNSW 对照（2026-08-08，`benchmarks/bench-hv-ann-hubnsw.md`）**：
  * N=30k：`hnsw_rs`（GSearch 同源）ef=200–400 召回 0.974–0.990，
    但查询 10.4–18 ms（精确 24.3 ms），只快 2.3–6.5×；HubNSW 单层化
    （scale=0.2）与多层召回几乎持平（中高 ef 略好 0.4–1.8 pp），
    不是关键变量；
  * N=10k：ef=50 召回 0.96–0.97、2.9 ms（精确 8.4 ms，快 ~3×）；
  * 结论修正：此前"4096 维 HNSW 召回天花板 <0.92"是 rust-cv 特定实现
    的结论；换成召回优先的 hnsw_rs 后 ANN 收益大幅缩水（接近精确扫描
    量级），降维仍是大规模下的首选。
- **距离层标定（2026-08-08，`benchmarks/bench-hv-ani-calibration.md`）**：
  HV 距离 vs skani ANI 在 135 个真实基因组上：仅 ANI 90–98% 区间中等
  可靠（Spearman 0.5–0.6），**≥99% 近缘与 <85% 远缘失效**（ρ≈0.38 /
  0.05）；D=16384 只改善中远缘，不救近缘；Mash 同种内 ρ=−0.97、
  ANI-truth recall@10 = 0.76 vs HV 0.62。**含义：物种内（≥98% ANI）
  聚类 / 选参考用 `dist mash` / `dist frac`；HV 定位为嵌入 / 粗筛 /
  查询路由（85–98% 带），不做 ANI 级精排。**
- **推荐路线**（据此更新）：
  * **≤10k**：精确扫描即够——SQLite + sqlite-vec（float[4096] +
    cosine，命令形态候选 `pgr hv db` / `pgr hv nn`），或零依赖的
    自研 SIMD 扫描（`hv.rs` 已有 `hv_dot`/`calc_distances`，SQLite 只
    存 BLOB/元数据）；实测精确扫描 10k 规模 ~10 ms/查询，批处理聚类
    完全够用。
  * **10k–30k 且接受 ~0.9 召回**：HNSW 可考虑（0.5–1.3 ms/查询），
    候选 = 自研纯 Rust HNSW（算法 ~几百行；`hnsw` crate 召回受限但快、
    `hnsw_rs` 召回高但慢，两版实现都是参考）或 usearch（接受 C++ 工具
    链）；若想要 0.96+ 召回，用 hnsw_rs 的 ef=50（10k 时 2.9 ms）；
    若已有 Necom 树/分组，先用代表路由 + clade 内图（§6.5 实测 30k
    同延迟下召回 +12–16 pp）。
  * **>30k 或百万级**：4096 维 HNSW 召回/延迟权衡不佳——rust-cv 召回
    受限（<0.92），hnsw_rs 召回虽高但查询已接近精确量级（30k 时
    2.3–6.5×），构建分钟–小时级（N=100k ≈ 4–5 min，N=1M ≈ 45 min
    单线程；hnsw_rs 构建再慢 ~2×）；**应先降维**
    （如 PCA 到 256–512 维）再评估 ANN，或回 SQLite 精确扫描 +
    分桶/倒排路线。
- 待验证：sqlite-vec 4096 维 float 查询延迟（目标机实测）；HV→f32 无损
  （按最大 N 校验值域）；cosine 排序 vs 精确 Jaccard 的一致性。
- Necom `mat from-vector`（f32，euclid/cosine/jaccard，O(N²)）是**批处理
  聚类消费者**，与查询型 NN 互补：小规模聚类可直接喂，查询型用本节方案。

### 6.5 外部知识辅助检索（系统发育 / 性状路由）—— 设想与评估

> 设想：HNSW 的入口点 / 分层代表目前由算法随机或贪心决定；能否用额外
> 知识（生物学性状、系统发育树距离）**人工指定更好的中间查询中介**，
> 改善查询过程？2026-08-08 讨论并实验（结论见本节末尾）。

**拆解"代表性点"的三层含义与可行性**：

1. **换顶层入口点**：基本无用。顶层只有极少数节点（GTDB 规模实测
   3 层、顶层 4 个节点），贪心下探到第 0 层后入口影响很小；HubNSW
   实验（`benchmarks/bench-hv-ann-hubnsw.md`）把多层压成单层后召回
   几乎不变，也侧面说明 4096 维下瓶颈不在入口选择。
2. **外部知识查询路由（最有希望）**：不用一棵全局 HNSW，而是按
   系统发育 / 性状分片（如按 Necom 聚类出的 clade），每个 clade 有
   代表节点（树距离或性状选的 medoid）和局部索引；查询先路由到最近
   的 R 个 clade，再在局部做 HNSW 或精确扫描。这是 GSearch 数据库
   分片策略（随机分片 → 汇总 top-K）的"知识化"升级，论文已证明分片
   汇总与整体等价；GSearch 三阶段搜索（nt → 全蛋白组 AAI → 通用基因
   AAI，按新颖度路由）也是外部知识参与路由的文献先例。
3. **外部知识重排（成本最低、最可控）**：ANN 照旧（如 ef=50 拿
   ~0.9 召回），对 top-k 候选用系统发育距离 / 性状做二次排序或软加权。
   先验错误时不会漏掉真近邻（候选已被 ANN 框住），只影响排序。

**风险与前提**：

- 先验必须与 HV 距离一致：HV 是 k-mer 内容相似度，系统发育 / 性状
  大致相关，但存在反例（水平基因转移、质粒获得、基因组大小差异 →
  "树近但 k-mer 远"或反之）。先验若作**硬路由规则**（只搜某 clade），
  判断错误时召回直接塌掉——GSearch 论文中病毒 nt 图陷局部最优就是
  k-mer 距离本身在稀疏数据上不可靠的先例。先验应作软信号或两级级联
  （路由 + 局部图内仍按 HV 精排），不宜替代距离。
- 增量维护：树 / 分组索引需要与数据库同步；HNSW 图本身已编码距离
  结构，知识路由只有在"图容易迷路"（远缘、稀疏、高维）时才体现价值。
- 工程现实：核心场景中 Necom 聚类时反正要建树，把树的分组与代表固化为
  查询路由表几乎是免费的；真正的问题是新基因组进库时树要不要重建。

**实验结论（2026-08-08，`benchmarks/bench-hv-ann-clade.md`）**：合成
cohort 注入 16 个 clade 结构（clade 内 2,048 共享核心 + 256 全局核心），
对比全局 HNSW vs 16 棵 clade 内 HNSW + 代表路由（R=1/2/4），recall@10
对全局精确 top-10：
- 30k 时路由在**同延迟**下召回高 12–16 pp（ef=10：0.940 vs 0.822；
  ef=20：0.982 vs 0.840；ef=50：0.984 vs 0.884）；全局卡在 0.908
  （ef≥100），路由卡在 0.988。收益机制 = 有效 N 缩小 16×，总构建成本
  相近（30k：41.8 s vs 16×2.3 s）。
- 10k 时全局 ef≥20 已 0.994，路由无必要；R>1 在完美先验下召回不变、
  只增延迟，其价值在兜底。
- 先验质量是决定因素：误路由查询 recall ≈ 0（只共享 256 全局核心），
  期望召回 ≈ (1−m)·R₁（m = 误路由率，近似线性）——先验应作软路由
  （R>1 / 候选重排），不宜硬性只搜一个 clade。
- 落地形态：Necom 聚类/构树的分组与代表天然是路由表，
  "代表精确路由 → clade 内 HNSW/精确"是 10k–30k 的更优选择。

## 7. 证据清单与验证计划（对标 GSearch 的验证体系）

> GSearch（NAR 2024）的可靠性验证分四层，我们逐层对照，把缺的证据
> 补成可执行实验。目标：链路上每一步都有"距离/检索/聚类/归档"的量化
> 依据，而不是沿用假设。

### 7.1 GSearch 的验证体系（可对标的方法库）

**A. 距离估计层**：① RMSE vs 真 Jaccard（6 种 sketch 算法互比，
m=12,000、k=16；蛋白 k=7）——SuperMinHash/ProbMinHash 理论方差最小、
Dashing 类 HLL 在小 Jaccard（75–78% ANI）时误差最大；② 相关性标定
（Spearman ρ）：ProbMASH-ANI vs BLAST-ANI/FastANI ≥78% ANI 时
ρ≈0.964、ProbMASH-AAI vs BLAST-AAI（95%>AAI>52%）ρ=0.90、通用基因
AAI ρ=0.939——**阈值由此定出**；③ mergeability（分段 sketch 可合并性）；
④ 度量性质（4 种 sketch 均 metric，HNSW 建图前提，非度量距离被排除）；
⑤ 基因组完整度 >50% 时 top-10 recall >80%。

**B. 检索层**：① 端到端 recall@5/10 vs 暴力 BLAST-ANI/AAI 真值，按
新颖度分层（近缘/无近缘/深分支，Table 3）；② **阈值过滤**——超过
0.9850（nt）/0.9720（aa）/0.9800（病毒）的匹配剔除后再算 recall，
把"距离估计失效区"与"图检索召回"分离；③ 查询留出（Tara 8466 / Ye
1000 均为库外）；④ 多真值交叉验证（BLAST-ANI、OrthoANI、FastANI、
MUMMER-ANI、BLAST-AAI）；⑤ M / ef_construct 参数敏感性。

**C. 端到端实用性**：① 与 GTDB-Tk 分类一致率（301 基因组 87.1%，
差异归因于污染/分类不一致）；② 单查询逐案核对（S13–S17，对比
Sourmash/Mash/Dashing）；③ 完整度与新颖度分层的召回矩阵。

**D. 性能/可扩展性**：① 建图/查询 wall time vs 线程数（24/128）与
分片数，O(log N) 拟合；② 内存/文件体积（ProbMinHash vs SetSketch）；
③ 同参数公平对比 Mash/Dashing/BinDash/Sourmash/skani/GTDB-Tk/
PhageCloud；④ hnswlib-rs vs 原版 hnswlib 三维对照；⑤ 分片 top-K
等价性（证明 + 实测）。

### 7.2 我们链路的证据现状与缺口

| 链路步骤 | 已有证据 | 缺口（对照 7.1） | 补证据实验（按优先级） |
|---|---|---|---|
| ① sketch/HV 距离 | dist mash 与 Mash 逐位一致；SIMD/HV 性能；合成数据 Jaccard 基准；RNG 检验；**P1 完成（2026-08-08）**：135 真实基因组上 HV 与 skani ANI 标定（`benchmarks/bench-hv-ani-calibration.md`） | frac/mini 未标定；无完整度/长度鲁棒性；sampler/k/D 扫描未做 | **P1 剩余**：完整度鲁棒性（模拟删 contig）；frac/mini 同 cohort 标定；sampler/k 扫描（对标 A①②⑤） |
| ② HV 最近邻检索 | 合成数据 recall@10（rust-cv/hnsw_rs/HubNSW/路由）；**P1 完成（2026-08-08）**：真实数据上以 ANI 为真值的**排序** recall@10（HV 0.62 vs Mash 0.76，`benchmarks/bench-hv-ani-calibration.md`） | 图检索（HNSW/路由）尚未在真实 HV 向量上以 ANI 真值重测；未用阈值过滤分离"距离误差 vs 检索误差" | **P1 剩余**：把 cohort 的 HV 向量喂 `hv_ann_clade`/`hv_ann_hubnsw` 类 bench，以 ANI top-10 为真值测图检索召回（对标 B①–④） |
| ③ 聚类/构树/剪枝（Necom） | 算法本身（MST 等） | 聚类结果 vs 已知分类/GTDB 的一致率；距离误差→聚类误差的传播 | **P2**：真实 cohort 聚类 vs GTDB 标签的 ARI/一致率；扰动距离矩阵看聚类稳定性（对标 C①） |
| ④ 选参考 + 比对（pgi） | pgi 对齐正确性/性能基准 | 参考选择策略对下游收益的量化；比对精度 vs ANI/AAI 真值 | **P2**：不同选参考策略（中心/最长/随机）→ pgi → pbit 压缩率对比；pgi 距离 vs ANI 标定（对标 A②C③） |
| ⑤ PBit 归档 | pbit 压缩基准（孤立） | 端到端收益：聚类→参考→比对→pbit vs 朴素方案 | **P2**：端到端压缩率 vs 单参考/无参考/直接 gzip（对标 C） |
| ⑥ 参数层 | D=4096 沿用 hypergen；采样/哈希理论笔记 | 参数（D、k、采样密度）对距离方差的实证曲线 | **P3**：真实数据上 D/k 扫描 → Jaccard/ANI 误差曲线（对标 A①） |

### 7.3 依赖与前提

- **P1 两项需要真实基因组 cohort**：本地无 FASTA 数据（`~/Scripts/genomes`
  只有 assembly 元数据表）。候选：RefSeq 按 `Escherichia.assembly.tsv`
  抽 50–200 个基因组（同种内 + 跨种，含已知 GTDB/RefSeq 分类），或
  用户已有的测序/下载数据。下载需用户侧网络。
- FastANI/skani 可用 `pgr` 之外的系统工具（用户已装 Mash；FastANI 可
  由用户安装）或 `pgr dist mash` 的 ANI 变换（已有 ProbMASH-ANI 公式）。
- 每项实验完成后更新本节状态与 `benchmarks/` 文档，结论回流 §6.4。

**P1 进度小结（2026-08-08）**：① 的距离标定已完成主体（HV + Mash vs
skani ANI，9,045 对）；关键发现——HV 只在 ANI 90–98% 中等可靠
（Spearman 0.5–0.6），≥99% 近缘与 <85% 远缘失效，D=16384 只救中远缘
不救近缘；Mash 同种内 ρ=−0.97，recall@10 比 HV 高 14 pp。② 的排序
recall@10 已测，图检索部分待补。详见
`benchmarks/bench-hv-ani-calibration.md`。

### 7.4 实验执行计划（真实数据，2026-08-08 起）

> 数据资产（全部本地，无需下载）：132,572 个 QC 组装；15,574 个 NR
> 代表（每物种 NR.lst）；每基因组 .msh sketch 与每物种 mash.nr.tsv；
> 全局 mash.dist.tsv（680 rep）+ tree.nwk + groups.tsv（height 0.4
> 聚类）；genome.taxon.tsv（name→species/genus/family/order/class）；
> assembly 元数据（N50/contig/大小）；minhash 与 bac120 两套树；
> 蛋白簇与 bac120 标记基因。工具：pgr（release）、mash、skani、necom、
> hmmsearch/mafft/trimal/FastTree、python3+numpy/scipy/pandas。
>
> 执行规则：每项实验先写方法（本表 + 对应 benchmarks 文档），跑完更新
> 状态与结论；结果回流 §6.4。状态列：⬜ 待做 / 🔄 进行 / ✅ 完成。

| # | 实验 | 方法要点 | 产出证据 | 状态 |
|---|---|---|---|---|
| 1 | frac/mini vs ANI | 同 cohort 跑 `dist frac`/`dist mini`（--list-files），与 skani ANI 求 Spearman/RMSE/recall@10，和 HV/Mash 并列 | dist 家族谁最贴近 ANI；"frac 用于 ANI"建议是否成立 | ✅ frac≈Mash（ρ0.97–0.99）；mini≈HV 近缘失效（详见 `benchmarks/bench-hv-ani-calibration.md` 补充节） |
| 2 | 完整度鲁棒性 | 子集基因组删 contig 至 90/70/50%，重算 HV/Mash/skani ANI，测距离漂移 | 完整度→距离误差曲线（对标 GSearch A⑤） | ✅ 近缘对 50% 完整度 HV +43%/Mash +84%，ANI 稳定；中等对稳定（详见 `benchmarks/bench-completeness-robustness.md`） |
| 3 | sampler/k/D 扫描 | 子集上扫 minimizer k/w、frac scale、syncmer、稀疏 s、D | 参数→误差曲线；默认参数有据 | ✅ frac 默认 s=1000 已近上限（RMSE≈1.15 ANI 点）；mini k/w/hasher 影响小；syncmer HV 不稳定（详见 `benchmarks/bench-parameter-scan.md`） |
| 4 | 长度/大小偏差 | 距离残差 vs N50/contig 数/总长（元数据现成） | HV 归一化对大小是否敏感 | ✅ 种内误差由碎片化驱动（N50 低/contig 多 → 误差大），大小差异也有贡献（详见同上） |
| 5 | 距离 CI | `dist frac --ci` 与自助法，看 CI 与误差关系 | 单对距离可靠性区间 | ✅ frac CI 对 skani ANI 覆盖率仅 8.4%（CI=采样误差，非金标准区间；详见 `benchmarks/bench-taxonomy-ci-tree.md`） |
| 6 | 真实 HV 图检索 | cohort 真实 HV 向量（pgi to-hv 或复用 hv.tsv）喂 HNSW，以 ANI top-10 为真值 | 图检索层真实数据召回（P1 ② 收尾） | ✅ 全局 HNSW recall_HV≥0.993、recall_ANI 0.664=精确；差距全在距离层（详见 `benchmarks/bench-hv-ann-real.md`） |
| 7 | 真实 clade 路由 | 用 cohort 自身 mash 距离聚类（necom hier/cut）做 clade，代表路由 + clade 内 HV 检索 | §6.5 真实验证（收益/误路由代价） | ✅ 物种硬路由 R=1 反降（0.70）；小 clade 是失败前提，需 clade≥K 成员（详见同上） |
| 8 | E. coli NR 全量 HV | 2,115 NR 基因组 HV 建库 + 精确 top-k 延迟 | 万级规模账 | ✅ 494 规模实测：精确 1.17 ms、HNSW 0.31 ms（recall_HV 0.985）；外推 2,115≈5 ms / 15k≈37 ms（详见 `benchmarks/bench-scale-and-pbit.md`） |
| 9 | 全 NR HV 可行性 | 15,574 NR 建库时间/内存估算 | 万级上限 | ✅ 估算：建库 ~3 CPU·时（0.7 s/基因组，可并行）；向量 249 MB；精确扫描 ~37 ms/查询（外推自 #8） |
| 10 | SQLite vs SIMD | BLOB+SIMD vs sqlite-vec 真实 HV 延迟（需装 sqlite-vec） | §6.2 实证 | ⬜（等安装） |
| 11 | Necom 聚类 vs 物种 | cohort 距离矩阵 → necom clust → ARI/NMI vs species 标签 | 聚类一致性（对标 C①） | ✅ mash K16 ARI 0.68/HV 0.65，K10 最优（mash 0.74/HV 0.57）（详见 `benchmarks/bench-clustering-validation.md`） |
| 12 | 聚类稳定性 | 距离加噪/自助重聚类，测再现度 | 聚类对距离误差敏感性 | ✅ ≤20% 噪声 ARI≥0.73，40% 崩至 0.36（详见同上） |
| 13 | groups.tsv 一致性 | groups.tsv 成员 vs species 标签/ANI 分布 | 现成分组能否当 clade/路由键 | ✅ 仅 13 个科/目级大组，物种纯度 0.03，不能当物种级路由键（详见同上） |
| 14 | 参考→pgi→pbit 端到端 | 子集（5–10 基因组）多参考 vs 单参考 vs 无参考压缩率 | 选参考收益量化（§7.2④） | 🔄 发现并修复 pbit `--paf` 跨组装命名 bug（含回归测试）；naive create 会静默丢数据；CIGAR 路径需 cg:Z PAF（pgr 转换缺口）+ 长共线比对（3 条约束，详见 `benchmarks/bench-scale-and-pbit.md` #14b）；真实压缩率待转换器/对齐器补齐 |
| 15 | ANI 物种阈值 | 同种/异种 ANI 分布找 ~95% 边界 | CLI/文档阈值实证 | ✅ 95% 阈值实证成立（同种误伤 0.8%/异种漏判 11.4%）（详见同上） |
| 16 | pgi 距离 vs ANI | 子集 dist pgi vs skani ANI 标定 | pgi 有偏结论量化 | ✅ 总体 ρ=−0.92，近缘段（≥95%）ρ=−0.71 弱（详见 `benchmarks/bench-pgi-calibration.md`） |
| 17 | pgi to-hv 一致性 | pgi to-hv 与 FASTA 直算 HV 对比 | 两条 HV 路径等价 | ✅ pgi→HV 保距（ρ=0.97 vs pgi 距离）；与 FASTA 直算参数不同不可直接比（详见同上） |
| 18 | 树一致性 | minhash 树 vs bac120 树 cophenetic 相关 | 距离树 vs 标记基因树吻合度 | ✅ 137 物种对 ρ=0.57；近缘段 ρ≈0.3–0.4 弱（详见同上） |
| 19 | 标记基因路由 | bac120 蛋白做快速先验，测路由准确率 | §6.5 生物学路由键选择 | ⬜ |
| 20 | dist mash 全量计时 | NR 子集两两计时 | 距离矩阵成本 | ✅ 0.13 ms/对（8 线程）：2,115 NR ≈ 5 min，15,574 ≈ 4.3 CPU·时 |
| 21 | HV 确定性 | 重复运行 dist hv 输出比对 | 可复现性 | ✅ 值完全确定，仅并行写盘行序不同（排序后逐行一致） |
| 22 | 结果回流 | 全部结论更新 §6.4/§7 + 用户文档 | 决策闭环 | 🔄 §6.4/§7 已随实验持续更新；用户文档（docs/*.md）阈值/默认参数建议待与语言问题一起处理 |

**执行顺序**：第一批 1、2、5、13、15、18、21（纯本地、快）；第二批 3、
4、6、7、11、12、16、17；第三批 8、9、10（等 sqlite-vec）、14、19、
20。每批跑完更新状态列并在 `benchmarks/` 落结果文档。

**执行日志（2026-08-08 续）**：
- 已完成 19/22（#1–#8、#11–#18、#20、#21），见上表状态。
- 本轮：① #19 bac120 标记基因路由准确率（hmmsearch 产物现成）；
  ② #14 扩展——高分歧样本（E. albertii）与双参考的 pbit 边际成本；
  ③ #8 全量——E. coli NR 2,115 基因组真实 HV 建库 + 检索延迟；
  ④ #10 sqlite-vec 侧视 crate 缓存是否可用，不可用则保持待装。
- 完成标准：每项有 benchmarks/ 结果文档 + 状态列更新。

**执行日志（2026-08-08 第三轮）**：
- #19 ✅：8×bac120 标记蛋白 aa 最近邻路由准确率 0.756（ANI 金标准上限
  0.800、HV 路由 0.822），详见 `benchmarks/bench-marker-routing.md`。
- #14 🔄：naive pbit create 静默丢数据（contig 名不匹配）→ 深挖发现
  `append_sample_with_paf` 命名门 bug，已修复 + 回归测试（617 测试通过）；
  CIGAR 路径 3 条约束（需 cg:Z PAF、段全覆盖、段内目标）记录于
  `benchmarks/bench-scale-and-pbit.md` #14b；真实压缩率待转换器补齐。
- #8 全量 2,115：未跑（494 规模实测 + 外推已覆盖规模账，见 #8 状态）。
- #10：sqlite-vec 未缓存，保持待装。
