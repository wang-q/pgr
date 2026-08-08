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
  + 头部 ≈ **16 KB/基因组**（与 `pbit.md` §附录的 16 KB/样本一致）。
- "bit（±1）"是每维 ±1 的**编码方式**，不是打包位向量；稀疏路径
  s=1 只保证值稀疏，文件仍按稠密 i32 存储。
- 规模账：1 万 ≈ 160 MB；4 万 ≈ 640 MB；百万 ≈ 16 GB。
- 距离语义：inter = H_A·H_B/D（点积恢复共享 k-mer 数）；余弦 =
  inter/√(card₁·card₂)，与 Jaccard 单调 → 近邻**排序**可直接用
  cosine/dot；精确 Jaccard 可由 cosine + 范数（card = ‖H‖²/D）恢复。

### 6.2 SQLite 方案对比

| 方案 | 许可 | 向量类型 | 检索方式 | 备注 |
|---|---|---|---|---|
| **sqlite-vec**（asg017） | MIT / Apache-2.0 | float[N] / int8[N] / bit[N] | `vec0` 虚拟表 KNN（`match … order by distance limit k`），**精确扫描、无 ANN** | **核心是 C 扩展**（用户明确不采用，未评测）；metric 支持 l2/cosine/l1/hamming；k ≤ 4096；HV 存 float[4096] 与 i32 等体积（16 KB），HV 值域（√N 量级）在 f32 精确整数范围内 |
| **sqlite-vector-rs**（quinnjr） | MIT / Apache-2.0 | Float4 / Float8 / Int1 / Int2 / Int4 / Float2 | PGVector-like `vector` 虚拟表 + **usearch HNSW（ANN）** | **纯 Rust 封装**（usearch 核心经 cxx FFI，C++ 源码构建）；`knn_match(distance, ?) LIMIT k`；**2,088 真实 HV 实测**（`benchmarks/bench-scale-and-pbit.md` #10b）：构建 1.6 s、查询 1.58 ms、recall@10 1.000、DB 69.9 MB（HNSW 图 +35 MB）。0.1.0 短板：library `register()` 是 `todo!()`（需手动注册 vtab）、shadow id 忽略用户 id、`sqlite3-ext-vtab` 需 bundled/static |
| **sqlite-vector**（rqlite） | **Elastic 2.0**（生产需商业授权） | Float32 / Float16 / BFloat16 / Int8 / UInt8 | 普通表 BLOB + 精确扫描，可量化 | 许可与 MIT 的 pgr 不兼容，**排除**（与上面两行是不同项目） |
| 自研：SQLite 存 BLOB + pgr SIMD 扫描 | MIT | i32 原样 | 线性扫描（`hv.rs` 已有 SIMD 点积） | 零新依赖；2,088 实测 DB 35 MB、预取后 2.46 ms/查询（#10a），**SQLite 非瓶颈** |

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
    ef=10 仅 0.73（`design/hv.md`）；同数据下
    GSearch 同源 `hnsw_rs` 0.3.4 可达 0.974–0.990，但查询慢 7–10×
    （`design/hv.md`），HubNSW 单层化（scale 0.2）
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

- **ANN 召回实测结论（2026-08-08，`design/hv.md`）**：
  * N=1k：HNSW 召回 0.98（ef≥10），比精确快 ~3×；
  * N=10k：ef=20 召回 0.90、0.46 ms/查询（精确 10 ms）；
  * N=30k：召回上限 0.92（ef=400 + ef_c=200），ef=10 仅 0.73；
    查询 0.65–1.5 ms vs 精确 30 ms（快 20–46×）；
  * 召回上限由查询 ef 决定，构建 ef_c 64→200 仅提升 0.6–1.2 pp——
    4096 维高维诅咒是主因，不是图质量（rust-cv 实现结论）。
- **hnsw_rs / HubNSW 对照（2026-08-08，`design/hv.md`）**：
  * N=30k：`hnsw_rs`（GSearch 同源）ef=200–400 召回 0.974–0.990，
    但查询 10.4–18 ms（精确 24.3 ms），只快 2.3–6.5×；HubNSW 单层化
    （scale=0.2）与多层召回几乎持平（中高 ef 略好 0.4–1.8 pp），
    不是关键变量；
  * N=10k：ef=50 召回 0.96–0.97、2.9 ms（精确 8.4 ms，快 ~3×）；
  * 结论修正：此前"4096 维 HNSW 召回天花板 <0.92"是 rust-cv 特定实现
    的结论；换成召回优先的 hnsw_rs 后 ANN 收益大幅缩水（接近精确扫描
    量级），降维仍是大规模下的首选。
- **距离层标定（2026-08-08，`design/hv.md`）**：
  HV 距离 vs skani ANI 在 135 个真实基因组上：仅 ANI 90–98% 区间中等
  可靠（Spearman 0.5–0.6），**≥99% 近缘与 <85% 远缘失效**（ρ≈0.38 /
  0.05）；D=16384 只改善中远缘，不救近缘；Mash 同种内 ρ=−0.97、
  ANI-truth recall@10 = 0.76 vs HV 0.62。**含义：物种内（≥98% ANI）
  聚类 / 选参考用 `dist mash` / `dist frac`；HV 定位为嵌入 / 粗筛 /
  查询路由（85–98% 带），不做 ANI 级精排。**
- **推荐路线**（据此更新）：
  * **≤10k**：精确扫描即够——SQLite + pgr SIMD 扫描（零依赖，
    2,088 实测 2.46 ms/查询）或 sqlite-vector-rs vtab（实测 1.58 ms、
    recall 1.000，但多 ~35 MB HNSW 图 + 0.1.0 成熟度风险）
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
- 已验证：sqlite-vector-rs vtab 2,088 真实 HV 延迟/召回（#10b）；HV→f32
  无损（值域 √N ≪ 2²⁴）
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
   实验（`design/hv.md`）把多层压成单层后召回
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

**实验结论（2026-08-08，`design/hv.md`）**：合成
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
| ① sketch/HV 距离 | dist mash 与 Mash 逐位一致；SIMD/HV 性能；合成数据 Jaccard 基准；RNG 检验；**P1 完成（2026-08-08）**：135 真实基因组上 HV 与 skani ANI 标定（`design/hv.md`） | frac/mini 未标定；无完整度/长度鲁棒性；sampler/k/D 扫描未做 | **P1 剩余**：完整度鲁棒性（模拟删 contig）；frac/mini 同 cohort 标定；sampler/k 扫描（对标 A①②⑤） |
| ② HV 最近邻检索 | 合成数据 recall@10（rust-cv/hnsw_rs/HubNSW/路由）；**P1 完成（2026-08-08）**：真实数据上以 ANI 为真值的**排序** recall@10（HV 0.62 vs Mash 0.76，`design/hv.md`） | 图检索（HNSW/路由）尚未在真实 HV 向量上以 ANI 真值重测；未用阈值过滤分离"距离误差 vs 检索误差" | **P1 剩余**：把 cohort 的 HV 向量喂 `hv_ann_clade`/`hv_ann_hubnsw` 类 bench，以 ANI top-10 为真值测图检索召回（对标 B①–④） |
| ③ 聚类/构树/剪枝（Necom） | 算法本身（MST 等） | 聚类结果 vs 已知分类/GTDB 的一致率；距离误差→聚类误差的传播 | **P2 部分完成（2026-08-08）**：真实 E. coli 物种内聚类（mash→UPGMA→cut）形成 7 簇结构（`design/genome-nn-query.md`）；物种级一致率见 #11（mash K10 ARI 0.74）；扰动稳定性见 #12 |
| ④ 选参考 + 比对（pgi） | pgi 对齐正确性/性能基准 | 参考选择策略对下游收益的量化；比对精度 vs ANI/AAI 真值 | **P2 完成（2026-08-08）**：中心/最长/随机三策略实测——**最长参考最优**（delta/gzip 0.520 vs 中心 0.554、随机 0.521），中心参考（距离和最小）反而最差；pgi 距离 vs ANI 标定见 #16（`design/genome-nn-query.md`） |
| ⑤ PBit 归档 | pbit 压缩基准（孤立） | 端到端收益：聚类→参考→比对→pbit vs 朴素方案 | **P2 完成（2026-08-08）**：端到端全链跑通（dist→necom→pgi→pbit），100 样本 × 3 策略 279 对 ≈ 12 分钟；LZ 路径 delta = gzip-9 的 52–55%，to-fa 覆盖率 ≥0.998（`design/genome-nn-query.md`） |
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
`design/hv.md`。

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
| 1 | frac/mini vs ANI | 同 cohort 跑 `dist frac`/`dist mini`（--list-files），与 skani ANI 求 Spearman/RMSE/recall@10，和 HV/Mash 并列 | dist 家族谁最贴近 ANI；"frac 用于 ANI"建议是否成立 | ✅ frac≈Mash（ρ0.97–0.99）；mini≈HV 近缘失效（详见 `design/hv.md` 补充节） |
| 2 | 完整度鲁棒性 | 子集基因组删 contig 至 90/70/50%，重算 HV/Mash/skani ANI，测距离漂移 | 完整度→距离误差曲线（对标 GSearch A⑤） | ✅ 近缘对 50% 完整度 HV +43%/Mash +84%，ANI 稳定；中等对稳定（详见 `design/genome-nn-query.md.md`） |
| 3 | sampler/k/D 扫描 | 子集上扫 minimizer k/w、frac scale、syncmer、稀疏 s、D | 参数→误差曲线；默认参数有据 | ✅ frac 默认 s=1000 已近上限（RMSE≈1.15 ANI 点）；mini k/w/hasher 影响小；syncmer HV 不稳定（详见 `design/hv.md`） |
| 4 | 长度/大小偏差 | 距离残差 vs N50/contig 数/总长（元数据现成） | HV 归一化对大小是否敏感 | ✅ 种内误差由碎片化驱动（N50 低/contig 多 → 误差大），大小差异也有贡献（详见同上） |
| 5 | 距离 CI | `dist frac --ci` 与自助法，看 CI 与误差关系 | 单对距离可靠性区间 | ✅ frac CI 对 skani ANI 覆盖率仅 8.4%（CI=采样误差，非金标准区间；详见 `design/genome-nn-query.md.md`） |
| 6 | 真实 HV 图检索 | cohort 真实 HV 向量（pgi to-hv 或复用 hv.tsv）喂 HNSW，以 ANI top-10 为真值 | 图检索层真实数据召回（P1 ② 收尾） | ✅ 全局 HNSW recall_HV≥0.993、recall_ANI 0.664=精确；差距全在距离层（详见 `design/hv.md`） |
| 7 | 真实 clade 路由 | 用 cohort 自身 mash 距离聚类（necom hier/cut）做 clade，代表路由 + clade 内 HV 检索 | §6.5 真实验证（收益/误路由代价） | ✅ 135 小 clade 路由反降（R=1 0.70）；**2,088 正向案例**：C=8/R=2 HV 路由 0.942、保 94% 全量 Mash recall；recall≈路由准确率（线性公式定量确认）（详见 `design/hv.md` 补充节） |
| 8 | E. coli NR 全量 HV | 2,115 NR 基因组 HV 建库 + 精确 top-k 延迟 | 万级规模账 | ✅ **2,088 全量实测（2026-08-08）**：精确 5.55 ms、HNSW ef10 0.45 ms（12×）、recall_HV 0.958–0.984；HV vs Mash 真值 recall@10=0.09（种内排序脱钩，详见 `benchmarks/bench-scale-and-pbit.md` #8b） |
| 9 | 全 NR HV 可行性 | 15,574 NR 建库时间/内存估算 | 万级上限 | ✅ 估算：建库 ~3 CPU·时（0.7 s/基因组，可并行）；向量 249 MB；精确扫描 ~37 ms/查询（外推自 #8） |
| 10 | SQLite vs SIMD | BLOB+SIMD vs sqlite-vector-rs 真实 HV 延迟/召回 | §6.2 实证 | ✅ 双路径都实测（`benchmarks/bench-scale-and-pbit.md` #10a/#10b）：BLOB+扫描 2.46 ms/查询、DB 35 MB；sqlite-vector-rs vtab 1.58 ms/查询、recall@10 1.000、DB 69.9 MB；裸 usearch ef≥50 recall≥0.999（sqlite-vec 因 C 核心未评测，用户决定） |
| 11 | Necom 聚类 vs 物种 | cohort 距离矩阵 → necom clust → ARI/NMI vs species 标签 | 聚类一致性（对标 C①） | ✅ mash K16 ARI 0.68/HV 0.65，K10 最优（mash 0.74/HV 0.57）（详见 `design/genome-nn-query.md.md`） |
| 12 | 聚类稳定性 | 距离加噪/自助重聚类，测再现度 | 聚类对距离误差敏感性 | ✅ ≤20% 噪声 ARI≥0.73，40% 崩至 0.36（详见同上） |
| 13 | groups.tsv 一致性 | groups.tsv 成员 vs species 标签/ANI 分布 | 现成分组能否当 clade/路由键 | ✅ 仅 13 个科/目级大组，物种纯度 0.03，不能当物种级路由键（详见同上） |
| 14 | 参考→pgi→pbit 端到端 | 子集（5–10 基因组）多参考 vs 单参考 vs 无参考压缩率 | 选参考收益量化（§7.2④） | ✅ **LZ 内容回退落地**（§8.5 路线 1，119 测试过）：6 样本边际 delta = gzip-9 的 51–81%（近缘 51–56%、E. albertii 81%），~100% 无损；完整参考对 draft 覆盖更好（`benchmarks/bench-scale-and-pbit.md` #14f/#14g） |
| 15 | ANI 物种阈值 | 同种/异种 ANI 分布找 ~95% 边界 | CLI/文档阈值实证 | ✅ 95% 阈值实证成立（同种误伤 0.8%/异种漏判 11.4%）（详见同上） |
| 16 | pgi 距离 vs ANI | 子集 dist pgi vs skani ANI 标定 | pgi 有偏结论量化 | ✅ 总体 ρ=−0.92，近缘段（≥95%）ρ=−0.71 弱（详见 `design/pgi-align.md.md`） |
| 17 | pgi to-hv 一致性 | pgi to-hv 与 FASTA 直算 HV 对比 | 两条 HV 路径等价 | ✅ pgi→HV 保距（ρ=0.97 vs pgi 距离）；与 FASTA 直算参数不同不可直接比（详见同上） |
| 18 | 树一致性 | minhash 树 vs bac120 树 cophenetic 相关 | 距离树 vs 标记基因树吻合度 | ✅ 137 物种对 ρ=0.57；近缘段 ρ≈0.3–0.4 弱（详见同上） |
| 19 | 标记基因路由 | bac120 蛋白做快速先验，测路由准确率 | §6.5 生物学路由键选择 | ✅ 8×bac120 标记蛋白 aa 最近邻路由准确率 0.756（ANI 上限 0.800、HV 0.822）（详见 `design/hv.md`） |
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
  0.800、HV 路由 0.822），详见 `design/hv.md`。
- #14 🔄：naive pbit create 静默丢数据（contig 名不匹配）→ 深挖发现
  `append_sample_with_paf` 命名门 bug，已修复 + 回归测试（617 测试通过）；
  CIGAR 路径 3 条约束（需 cg:Z PAF、段全覆盖、段内目标）记录于
  `benchmarks/bench-scale-and-pbit.md` #14b；真实压缩率待转换器补齐。
- #8 全量 2,115：未跑（494 规模实测 + 外推已覆盖规模账，见 #8 状态）。
- #10：sqlite-vec 未缓存，保持待装。

**执行日志（2026-08-08 第四轮）**：
- 目标：在共线性的完整 E. coli 对上跑通 pgi→PSL→cg:Z PAF→pbit 全链路，
  拿到真实 delta 压缩率（绕开重排对撞上的 CIGAR 约束 3）。
- 方法：候选完整 E. coli 样本逐个 align pgi → PSL 块贪心单调链 → 分块
  difflib 生成 =/X/I/D CIGAR → 修好的 debug pgr 建 pbit → to-fa 覆盖率；
  选覆盖率最高者记录归档大小 vs gzip-9。
- 成功标准：#14 状态转 ✅（带真实压缩率 + 约束说明）；失败则如实记录
  覆盖率下限并保留约束文档。

**执行日志（2026-08-08 第四轮结果）**：
- 候选完整 E. coli 全部重建 0%——CIGAR 逻辑经两个新单元测试验证正确
  （多段无 gap 往返 ✓；跨段删除只丢受影响段 ✓），真实失败根因 =
  pgi PSL 碎片化（重叠块 2–3× 覆盖）无法朴素链化；`--min-shared`/
  `--freq` 参数不改善。结论：#14 收尾为"bug 已修 + 约束精确化 +
  链级 cg:Z 生产者列为后续功能"（`benchmarks/bench-scale-and-pbit.md`
  #14c）。本轮新增 2 个测试（共 117 个 pbit 测试通过）。

**执行日志（2026-08-08 第五轮）**：
- 新尝试：**每条 PSL 链单独一条 PAF 记录**（+ 链按链记录；− 链按
  contig 合并成整 RC 查询 CIGAR），不再跨链合并——此前 0% 重建是
  跨链合并大 CIGAR 撞约束 3 所致，修复二进制后此路径未试过。
- 方法：converter v3 → pbit create（debug pgr）→ to-fa 覆盖率 → 若
  ≥80% 则记录归档大小 vs gzip-9，#14 转 ✅；否则如实记录并结束 #14。

**执行日志（2026-08-08 第五轮结果）**：
- 修复转换器 `re.sub` 误删数字 bug 后 0 畸形 CIGAR，重建仍 0%。
  实测 pgi 链跨度中位数 1,142 bp、p90 3.2 kb——**链粒度 << 段粒度**，
  pbit"单记录全覆盖 4096 段"在真实基因组上无段可编码。两条后续路线
  （pgi 长链化 / pbit 跨记录组装）挂 todo；#14 收尾（bug ✅ + 约束与
  量化根因 ✅ + 真实压缩率待功能）。

**执行日志（2026-08-08 第六轮）**：
- 目标：① #8 全量实测——E. coli NR 2,115 基因组 HV 建库 + 检索延迟
  （此前 494 规模外推）；② 写"证据汇总与设计决策建议"（§8），把
  20 项实验收敛成核心场景每步的决策依据（#22 技术主体）。
- 方法：NR.lst 全量 2,115 → xargs -P 8 并行 pgi build + to-hv → 
  `bench hv_ann_real`（空 ANI，recall_HV + 延迟）→ 记录到
  `benchmarks/bench-scale-and-pbit.md`；同步撰写 §8 汇总。
- 完成标准：#8 状态更新为实测数字；§8 完成并与各 benchmarks 文档
  交叉引用一致。

**执行日志（2026-08-08 第六轮结果）**：
- #8 ✅ 全量实测：2,088 基因组（8 并行流式建库 ~15 min，避免 tmpfs
  配额打爆——首次尝试 41 MB/pgi × 2,115 ≈ 87 GB 超限，改流式后
  峰值 ~330 MB）；精确 5.55 ms、HNSW ef10 0.45 ms（12×）、recall_HV
  0.958–0.984；**HV vs 全量 Mash 真值 recall@10 = 0.09**（种内近同
  株排序脱钩，随机基线 ≈0.05）——种内精细排序 HV 不可用（任意规模）。
- §8 证据汇总与设计决策建议已完成（13 条决策 + 工作流修订 + GSearch
  对标 + 未决项）。

**执行日志（2026-08-08 第七轮）**：
- 新假设：pbit create 的 `--segment-size`（默认 4096）是参数——把段降
  到 1,024/512 bp，pgi 链中位数 1,142 bp 即可全覆盖段，可能绕开
  "链粒度 << 段粒度"的结论，直接拿到真实压缩率。
- 方法：per-chain PAF（make_paf_chains.py，已修复 CIGAR 生成）+
  修好的 debug pgr + `-s 1024`/`-s 512`，对完整 E. coli 对
  （00_3076 vs 00_3230）测 to-fa 覆盖率与归档大小 vs gzip-9。
- 成功标准：覆盖率 ≥80% 则 #14 转 ✅（真实压缩率 + "段大小需按链
  粒度调参"的结论）；否则如实记录并维持 #14 待功能状态。

**执行日志（2026-08-08 第七轮结果）**：
- 段大小调参无效（512 时仅 1 段编码）——最根本约束浮出：**段相位对齐**
  （目标起点 mod 段长 == 0）。新增单元测试
  `test_append_sample_with_paf_indel_breaks_phase`（1 bp 插入后全部下游
  段被跳过）钉死该约束。真实基因组 indel 永久破坏相位 → CIGAR 编码
  需"跨相位组装"设计改动（挂 todo）。#14 终态 = bug ✅ + 三条约束全部
  精确化（118 个 pbit 测试）。

**执行日志（2026-08-08 第八轮）**：
- 目标：补 #7 路由的**正向规模案例**——135 cohort 因 clade 太小路由
  失效；现在用 2,088 E. coli（HV 向量 + 全量 Mash 矩阵现成）聚成成员
  充足的真实 clade，测"**HV 路由先验 + clade 内 Mash 精确检索**"
  （§6.5 推荐形态）的 recall 与路由准确率。
- 方法：mash2115.tsv → necom ward → cut（C=8/16）→ 每 clade 一个代表；
  查询按 HV 点积路由到 top-R clade（R=1/2）→ 路由 clade 内 Mash 精确
  top-10 → recall vs 全量 Mash top-10；同时统计 HV 路由准确率。
- 成功标准：路由准确率与 recall 记录进文档；若 C=16/R=1 时 recall
  ≥0.9 则"路由生效前提 = clade ≥K 成员"得到规模级正向证据。

**执行日志（2026-08-08 第八轮结果）**：
- 2,088 E. coli 正向路由案例完成：C=8/R=2 时 HV 路由准确率 0.942、
  路由后 Mash recall 0.940（搜 25% 库保 94% recall）；**routed recall ≈
  路由准确率**（4 组数据 ±0.004），§6.5 线性容错公式定量确认；
  C=16（小 clade）准确率跌到 0.46——clade ≥K 成员前提再次验证。

**执行日志（2026-08-08 第九轮）**：
- 目标：实现 §8.5 路线 1——**LZ 兜底内容匹配化**：名字匹配失败时，
  按 canonical k-mer 倒排索引找"内容最相似的参考段"，LZ 编码（不改
  归档格式、无需 PAF）。这应解锁 pbit 真实压缩率。
- 设计：`Compressor` 加惰性 `ref_kmer_index: HashMap<u64, Vec<u32>>`
  （canonical 15-mer → 参考段 id）；`best_ref_group(seg)` 投票取最高
  共享 k-mer 的参考段；append_sample / append_sample_with_paf 的名字
  匹配失败回退到它；方向由 LZ 内部 alt 分支处理。
- 验证：① 单元测试——跨组装命名样本（无 PAF）LZ 内容匹配往返；
  ② 真实完整 E. coli 对 pbit create（无 PAF）→ to-fa 覆盖率与归档
  大小 vs gzip-9。
- 完成标准：测试过 + 覆盖率 ≥80% 则 #14 转 ✅（真实压缩率 + 路线 1
  落地）；否则记录失败模式。

**执行日志（2026-08-08 第九轮结果）**：
- **#14 ✅ 路线 1 落地**：canonical k-mer 倒排索引 + `best_ref_group`
  内容匹配，119 个 pbit 测试过（含跨组装无损往返 + 原相位/跨段测试
  更新为"LZ 回退恢复 CIGAR 丢弃的段"）。
- 真实数据：完整近缘对 100% 重建（delta = gzip-9 的 53%）、draft 近缘
  99.99%（57%）、E. albertii 100%（78%）；归档为 2bit 参考 + delta 的
  结构化格式。真实压缩率终于拿到（`benchmarks/bench-scale-and-pbit.md`
  #14f）。

**执行日志（2026-08-08 第十轮）**：
- 目标：① 多样本 pbit 边际成本——ref + 6–8 个跨亲缘样本（完整/draft
  E. coli + E. albertii），测每样本边际 delta 与覆盖率；② 查 draft
  近缘对缺失的 737 bp 是什么（无内容匹配的 contig/质粒？）。
- 方法：增量 create（1→N 样本）比归档大小差 = 边际成本；to-fa 覆盖率；
  draft 缺失部分对照原始 FASTA contig 列表。
- 完成标准：边际成本表 + 缺失 contig 归因记录到
  `benchmarks/bench-scale-and-pbit.md`。

**执行日志（2026-08-08 第十轮结果）**：
- 6 样本边际 delta = gzip-9 的 51–81%（完整近缘 51–56%、draft 56/53%、
  E. albertii 81%）；6 样本归档 6.94 MB（含参考），重建 ≈100%
  （32.38 Mbp）。
- draft 737 bp 缺失归因：3 个 contig 的边缘段（样本特有序列，参考无
  匹配）；换完整参考后仅丢 53 bp——**完整参考覆盖更好**，`to-fa` 覆盖
  率应作归档质量门。

**执行日志（2026-08-08 第十一轮）**：
- 目标：① #10 可做的一半——SQLite 存 HV BLOB + 扫描的查询延迟
  （Python 标准库 sqlite3 + numpy 近似 pgr SIMD 路径，不引新依赖；
  sqlite-vec 侧仍等安装）；② #22 技术内容——把发现整理成"用户文档
  改动清单"（§8.6），语言处理时直接套用。
- 方法：2,088 个真实 HV（i32 4096）写入 SQLite BLOB → 逐查询读全部
  BLOB + 转 f32 + 点积 top-10 → 测延迟；对比纯内存扫描（#8b 5.5 ms）。
- 完成标准：BLOB 路径延迟与"是否可行"结论记录；§8.6 清单完成。

**执行日志（2026-08-08 第十一轮结果）**：
- #10 BLOB 侧 ✅：2,088 向量 SQLite DB 35 MB（入库 0.03 s），预取后
  扫描 2.46 ms/查询、一次性取+转换 68 ms——**SQLite 存储不是瓶颈**，
  ≤10k 精确扫描方案成立；sqlite-vec 对比仍待安装（#10 转 🔄）。
- §8.6 用户文档改动清单完成（6 个文档 × 8 项改动，含依据证据），

**执行日志（2026-08-08 第十二轮）**：
- 目标：#10 换道——用户明确不用 C 核心的 sqlite-vec，改测纯 Rust 的
  **sqlite-vector-rs 0.1.0**（PGVector-like vtab + usearch HNSW；
  MIT/Apache-2.0，与 rqlite 的 Elastic-2.0 "sqlite-vector" 是不同项目）。
- 踩坑记录（防重走）：① library 模式 `register()` 是 `todo!()`，需按
  loadable-extension 入口等价逻辑手动注册 vtab；② `sqlite3-ext-vtab`
  默认（非 static）走运行时 API 表，library 用法直接 SIGSEGV，需
  `bundled`/`static`；主机缺 libsqlite3-dev，故用 bundled；③ vtab
  shadow 表 `id` 是 AUTOINCREMENT，用户提供的 id 被忽略，返回 id 从 1
  起（比对真值要 `id-1`）。
- 结果 ✅（`benchmarks/bench-scale-and-pbit.md` #10b）：2,088 真实 HV，
  裸 usearch HNSW 构建 1.4 s、查询 159–842 µs（ef10–200）、recall@10
  0.984–1.000（ef≥50 ≥0.999）；vtab ef64 端到端构建 1.56 s、查询
  1.58 ms、recall 1.000，DB 69.9 MB（HNSW 图 +35 MB）；重开加载 58 ms
  后 warm 查询 1.55 ms。**#10 转 ✅（22/22）**。
  待语言处理时落地（#22 技术主体就绪）。

**执行日志（2026-08-08 第十三轮）**：
- 目标：**核心工作流端到端验证（P2，§7.2 ③④⑤）**——真实 E. coli
  物种内跑通"dist mash → necom 聚类 → 选参考 → pgi 比对 → pbit 归档"
  全链，量化参考选择策略对压缩率的影响（用户核心场景）。
- 方法：2,088 个 E. coli NR 中 farthest-point 采样 100 个；mash 全对
  距离（4,950 对）→ UPGMA → cut k=7（44/35/9/8/2/1/1）；每簇按
  center（距离和最小）/ longest（总长最大）/ random 各选 1 参考；
  簇内非参考样本 × 参考做 align pgi → psl to-paf → pbit create
  （纯 LZ + --paf CIGAR 双路径）；delta/gzip-9 对比 + to-fa 覆盖率。
- 结果 ✅（`design/genome-nn-query.md.md`）：279 对，
  **longest 最优**（delta/gzip 0.520 vs center 0.554、random 0.521；
  longest < center 71/98 查询）；center 参考（典型 draft、内容覆盖
  少）反而最差；单参考随机性影响 3–4 pp；to-fa 覆盖率 ≥0.998；
  CIGAR 与 LZ 差平均 32 B（相位约束，基本回退 LZ）。
- 补充（用户质疑"参考必须 Complete"）：clade 0/1 配对实验（75 对，
  控制"最长"只变组装级别）——draft-longest 参考压缩率反而更优
  （0.505 vs Complete 0.525，56/75 对），机制 = contig 颗粒对齐 +
  质粒/未装配内容共享；cohort 中 7 簇有 2 簇无 Complete 成员。
  **"参考必须 Complete"的规则不来自压缩率**，而是比对/可解释性/
  一致性（§8.1 决策 14 修正）。

**执行日志（2026-08-08 第十四轮）**：
- 目标：用户要求 pbit **严格无损**（"什么样的东西进去，什么样的
  东西出来"）；此前只核对了碱基计数覆盖率（99.99%），不严格。
- 核对结果：逐碱基（名称+顺序+位置）抽查 12 个真实样本，7 个有
  缺失（6–614 bp）；极端样本 Es_coli_188（66 万 N）缺 20 万 bp。
  根因 = 无参考匹配的段/contig 被静默跳过（内容匹配路径
  `best_ref_group=None` 丢弃段；PAF 路径丢弃整 contig）。
- 修复（v1006）：`DeltaEncoding::Raw = 2`——未匹配段 flate2 压缩
  原文存储（挂 ref_group 0、解码不读参考段），两处跳过分支全部
  替换；`test_raw_fallback_lossless` 钉死往返一致；pbit create
  帮助文本同步。修复后逐碱基核对 10 样本：9/10 完全一致，188 的
  差异全部为简并碱基 → N（用户明确允许的已知边界）。
- 影响：Raw 段占比小，端到端压缩率结论（决策 14）不变；119 个
  pbit 测试 + 全量 621 测试通过，fmt/clippy 干净。

**执行日志（2026-08-08 第十五轮）**：
- 目标：用户确认 pbit 设计意图 = **省掉 PAF 文件依赖**（放入时基于
  PAF 编码、比对信息内嵌、可从归档导出），推荐 PAF 生成链路 =
  `align pgi → chainnet → maf to-paf`。
- 验证：链路跑通，但发现两个真实障碍并修复：
  ① chainnet 默认在序列名前加 basename 前缀（`GCF_xxx.NZ_...`），
  与样本 FASTA contig 名不匹配 → PAF 驱动全部落空（此前"CIGAR 0
  段"的主因其实是这个，不是相位约束）；修复 = 推荐链路加
  `--t-name '' --q-name ''`（空前缀，纯 contig 名）。
  ② 即使名字匹配，CIGAR 路径仍要求"单条记录全覆盖 4kb 段 + target
  不跨段"（相位约束）→ **v1007 格式升级**：SegmentDesc 加 `q_start`、
  `ref_start/ref_end` 改参考文件全局坐标；CIGAR 段级混合编码（PAF
  覆盖部分 CIGAR + 剩余 Raw），不再要求整段对齐。
- 结果 ✅：真实 98.6% 完整对（00_3076 vs 00_3230）CIGAR 1246 段 /
  LZ 250 / Raw 451；to-fa 严格无损；delta/gzip 0.393（纯 LZ 0.539，
  **CIGAR 混合 -14.6 pp**）。
- 新增 `pgr pbit to-paf`：从归档导出内嵌比对（每个 CIGAR 段一条
  PAF，12 列 + cg:Z）；闭环验证 pbit → to-paf → 重新 create 编码
  分布一致 + 严格无损。
- 注：fas_xlsx 有并行审计改动，未触碰（用户指示）。

## 8. 证据汇总与设计决策建议（2026-08-08）

> 目的：把 §7.4 的 20+ 项实验收敛为对核心场景（物种内聚类选参考 +
> PBit 归档）每一步的**决策依据**。每条建议后附证据与来源；全部数据
> 来自真实 Enterobacterales 数据（135–2,115 基因组，skani ANI 金标准）。

### 8.1 决策 → 建议 → 证据

| # | 决策点 | 建议 | 关键证据 | 来源 |
|---|---|---|---|---|
| 1 | 距离度量 | 物种内（≥98% ANI）聚类/选参考用 **`dist mash`/`dist frac`**；HV 定位为嵌入/粗筛/路由（85–98% 带），不做 ANI 精排 | frac≈Mash（ρ0.97–0.99，recall@10 0.76）；HV 近缘 ρ0.38、<85% ρ0.05；HV recall@10 0.62 | `design/hv.md` |
| 2 | frac 参数 | 默认 scale=1000 合理（RMSE≈1.15 ANI 点）；s=100 仅近缘段微增，s=10000 近缘变差 | 30 基因组扫描 | `design/hv.md` |
| 3 | minimizer | mini 近缘缺陷是结构性的（同种 ρ≈0.61），k/w/hasher 无法解决 | 同上 | `design/hv.md` |
| 4 | HV 维度 | D=16384 只救中远缘，近缘段无改善——近缘场景不必升维 | D 对比 | `design/hv.md` |
| 5 | 完整度 | 种内聚类前过滤低质量组装；完整度 <50% 距离膨胀 HV +43%/Mash +84%（ANI 稳定） | 删 contig 实验 | `design/genome-nn-query.md` |
| 6 | 检索 | ≤10k 精确扫描够用（494 规模 1.17 ms/查询）；HNSW 图检索误差可忽略（recall_HV≥0.993），HV→ANI 差距全在距离层 | 真实 HV ANN | `design/hv.md`、`benchmarks/bench-scale-and-pbit.md` |
| 7 | 路由 | clade 须 ≥K 成员；小 clade 硬路由有害（R=1 掉到 0.70）；**生物学先验可用**——8 个 bac120 标记蛋白路由准确率 0.756（ANI 上限 0.800）；2,088 规模 C=8/R=2 时 HV 路由 0.942、保 94% 全量 Mash recall，recall≈路由准确率 | 合成 + 真实（135/2,088）+ 标记实验 | `design/hv.md`、`design/hv.md`、`design/hv.md` |
| 8 | 聚类 | Necom 物种级恢复良好（mash K10 ARI 0.74 / HV 0.57）；≤20% 距离噪声稳定（ARI≥0.73） | 聚类验证 | `design/genome-nn-query.md` |
| 9 | ANI 阈值 | 95% 物种边界实证成立（同种误伤 0.8%、异种漏判 11.4%）；90% 为保守提示 | 物种标签分布 | `design/genome-nn-query.md` |
| 10 | 树 | minhash 树近缘段与 bac120 树弱一致（ρ0.3–0.4）→ 物种级参考拓扑用 bac120 | 两树 cophenetic | `design/genome-nn-query.md` |
| 11 | pgi | pgi 距离近缘段弱（ρ−0.71）；`pgi to-hv` 保距（ρ0.97），可作 pgi→HV 嵌入 | pgi 标定 | `design/pgi-align.md` |
| 12 | PBit | **v1006 起 ACGTN 严格无损**（无匹配段 Raw 存储；唯一有损 = 简并 → N，用户允许）；**v1007 起 CIGAR 支持任意参考区间**（段级混合编码），PAF 驱动真正生效（98.6% 对 delta/gzip 54%→39%）；**v1009 起 `pbit to-paf` 无损还原输入 PAF**（809/809 行逐字段一致：大链按 `paf_id` 合并且 cg/cs/gi/bi 重算、ms 存表；碎链行原样存储；含 PAF 恢复区 delta/gzip 0.448）；**v1010 起 Identity 零载荷**（纯 `=` 段指向参考区间） | pbit 深挖 v1007–v1010，能力状态统一维护于 `design/pbit.md` | `benchmarks/bench-scale-and-pbit.md` #14a–d/#14h/#14k/#14l |
| 13 | 规模 | 2,088 全量实测：精确 5.55 ms、HNSW ef10 0.45 ms（12×）；建库 15 min（8 并行流式）；15,574 外推 ≈ 37 ms | `benchmarks/bench-scale-and-pbit.md` #8b | `benchmarks/bench-scale-and-pbit.md` |
| 14 | **参考选择** | 簇内选参考**不要选距离中心**（典型 draft、内容覆盖少，压缩率最差：longest vs center 差 3.4 pp，71/98 查询 longest 更优）；**"必须 Complete"的规则不来自压缩率**——配对实验 draft-longest 反而 -2 pp（56/75 对，contig 颗粒对齐 + 质粒/未装配内容共享）；Complete 优先的正当理由是比对质量/坐标可解释性/下游一致性；簇内无 Complete 时用最长 draft 或并入相邻簇；单参考随机性影响 3–4 pp，小簇可多参考摊平 | 100 样本 × 3 策略 + Complete/draft 配对 | `design/genome-nn-query.md` |

### 8.2 核心场景工作流（证据修订版）

1. **输入 + 质控**：QC 名单（NWR pass.lst 式）+ N50/完整度门（决策 5）；
2. **距离**：`dist mash`/`dist frac` 两两（决策 1、2），输出 pair TSV；
3. **聚类/构树**：Necom（ward/MST），以物种标签 ARI 校准 K（决策 8）；
4. **选参考**：组内**最长/高完整度优先**（决策 14；不选距离中心），
   **共线性优先**（pbit 约束 3，决策 12）；
5. **比对**：`align pgi` → PSL（注意链粒度限制，决策 12）；
6. **归档**：`pbit create --paf` + **`to-fa` 覆盖率质量门**（决策 12）；
7. **查询/检索**（可选）：HV 嵌入 + SQLite/精确扫描 ≤10k；HNSW 仅在大
   clade 路由下使用（决策 6、7）。

### 8.3 与 GSearch 验证体系的对标

| GSearch 验证层 | 我们的对应证据 |
|---|---|
| A 距离估计（RMSE/相关性/完整度/mergeability/度量性） | frac/mash/HV/mini vs ANI 分层 Spearman+RMSE；完整度；参数扫描（§7.4 #1–5） |
| B 检索召回（分层、阈值过滤、多真值） | ANI 真值 recall@10 分层 + HV 真值图检索误差分离（#6、#1） |
| C 端到端（分类一致率、单查询核对） | Necom vs 物种 ARI；聚类稳定性；标记路由（#11/12/19） |
| D 性能/可扩展（线程/分片/O(log N)/内存） | 494/2,115 规模延迟实测；mash 全量计时；HV 建库成本（#8/9/20） |

### 8.4 未决与后续

- #10 sqlite-vector-rs vtab 实测完成（#10b）；sqlite-vec（C 核心）按
  用户决定不再评测；
- pbit 真实压缩率：依赖现有 chainnet 链路（**pgi 长链链化已明确不做**，
  2026-08-09 用户裁定——自研 chain 不如 UCSC chainnet，见
  `design/pbit.md` §PAF 驱动编码的演进）；"pbit 跨记录组装"仍可选；
- pgr `psl to-paf` 的 cg:Z 生产者（链级，与上一条相关）；
- 用户文档（docs/*.md）阈值与默认参数建议随语言处理一起落地。

### 8.5 pbit 编码演进（已收拢至 design/pbit.md）

pbit 的编码路线与约束现状已统一维护在
[design/pbit.md](../design/pbit.md)（"PAF 驱动编码的演进"章节；三条路线
现状：LZ 内容匹配化 v1006 ✅、跨相位 CIGAR v1007 ✅（+ v1009 to-paf、
v1010 Identity）、pgi 长链链化明确不做）。

### 8.6 用户文档改动清单（发现 → docs/*.md，语言处理时套用）

| 文档 | 改动 | 依据（证据） |
|---|---|---|
| `dist.md`（frac） | ① 注明 ANI 95% CI 只覆盖 FracMinHash 采样误差、不覆盖金标准 ANI（覆盖率仅 8.4%）；② 注明 protein 场景应降 `--scale`（默认 1000 对 ~300 aa 蛋白 sketch 仅 1 元素）；③ 强调 frac/mash 是 ANI 数值估计的推荐命令 | #5、#3、#1 |
| `dist.md`（hv） | ① HV 距离只适合粗分层（85–98% ANI），**种内（≥98%）排序不可靠**（2,088 规模 recall vs Mash 0.09）；② D=16384 不改善近缘分辨率，默认 4096 即可；③ 输出行序不保证稳定（并行写），按对去重使用 | #1、#8b、#21 |
| `dist.md`（mini） | 注明 minimizer 采样近缘分辨率结构性地弱于 frac/mash（同种 ρ≈0.61） | #1/#3 |
| `pbit.md` | ① naive create 只按 contig 名匹配——跨组装样本会**静默丢数据**，必须用 `--paf` 或依赖新 LZ 内容匹配（版本说明）；② 归档后跑 `to-fa` 覆盖率质量门；③ 近缘样本边际 delta ≈ gzip-9 的 51–56%、分歧 ~81%（实测）；④ 完整参考对 draft 样本覆盖更好 | #14a–g |
| `align-pgi.md` | 注明 pgi PSL 链粒度 ~1 kb（中位），与 pbit CIGAR 的 4 kb 段不匹配；建议 CIGAR 路径用长链对齐器或等格式升级 | #14c–e |
| `pgi.md` | 注明 pgi 距离近缘段（≥95% ANI）相关性弱（ρ≈−0.71），适合粗距离/索引，不做 ANI 精排 | #16 |
| 新小节（dist.md 或独立） | ANI 阈值：≥95% 同种（实测误伤 0.8%/漏判 11.4%）；90% 为保守提示 | #15 |

> 状态：清单已定，待语言处理时落地（#22）。

---

## 证据附录：真实数据聚类验证：Necom 聚类 vs 物种标签、距离噪声稳定性（#11/#12）

> 日期：2026-08-08。cohort 同 `design/hv.md`（135 基因组，
> 16 个物种目录/14 个有效物种标签）。对应 `design/genome-nn-query.md`
> §7.4 #11/#12。

## #11 Necom 聚类 vs 物种标签

方法：mash（NWR .msh）与 HV（k21/w5/D4096）距离矩阵 → `necom mat
to-phylip` → `necom clust hier --method ward` → `necom cut simple -k K`
→ `necom eval partition` vs 物种标签（cluster format）。

| 距离 | K | ARI | 同质度 homogeneity | RI |
|---|---|---|---|---|
| mash | 10 | **0.739** | 0.785 | 0.931 |
| mash | 16 | 0.681 | 0.878 | 0.935 |
| mash | 25 | 0.588 | 0.902 | 0.924 |
| mash | 40 | 0.333 | 0.924 | 0.900 |
| HV | 10 | 0.573 | 0.711 | 0.899 |
| HV | 16 | 0.646 | 0.835 | 0.924 |
| HV | 25 | 0.525 | 0.909 | 0.916 |
| HV | 40 | 0.256 | 0.928 | 0.894 |

结论：物种级划分（K≈10–16）两种距离都能较好恢复（ARI 0.57–0.74），
Mash 各 K 下都优于 HV（+0.07–0.17）；K 过大（40）时 ARI 骤降（过度
细分，同质度上升但一致率崩）。注意 E. coli 与 Escherichia sp. 混聚在
一起（sp. 标签实际多为 E. coli），物种标签本身有噪声。

## #12 距离噪声下的聚类稳定性

方法：HV 距离矩阵加乘性高斯噪声（σ = 5/10/20/40%），每水平 3 个种子
重聚类（ward, K=16），与原聚类算 ARI。

| 噪声 σ | ARI（3 次重复均值 ± SD） |
|---|---|
| 5% | 0.918 ± 0.014 |
| 10% | 0.821 ± 0.141 |
| 20% | 0.726 ± 0.110 |
| 40% | 0.363 ± 0.101 |

结论：物种级聚类对 ≤20% 距离噪声稳健（ARI ≥0.73）；40% 噪声破坏聚类。
结合 `design/hv.md`（HV 近缘区间排序噪声大），物种级
划分不受影响，但**种内亚结构（如 E. coli 内部的分群）对距离噪声更
敏感**，后续应专门用 E. coli 子集验证。

## 复现

```bash
# 距离矩阵 → phylip → ward 聚类 → cut → eval：
# necom mat to-phylip mash.pair -o mash.phylip
# necom clust hier --method ward mash.phylip -o tree.nwk
# necom cut simple -k 16 tree.nwk -o clust.tsv
# necom eval partition clust.tsv --other species.pair \
#     --input-format cluster --other-format pair
# 稳定性脚本：/tmp/hv_calib/run_stability.py
```

---

## 证据附录：基因组完整度对 HV / Mash / ANI 距离的影响（#2）

> 日期：2026-08-08。方法：1 个 E. coli 锚点基因组，3 个不同亲缘靶
> （近缘 E. coli ANI≈98% / 中等 E. albertii ANI≈90% / 远缘 Yersinia
> enterocolitica），对靶基因组随机删 contig 至 90/70/50% 完整度
> （seed=42 固定），重算 HV（k21/w5/D4096）、Mash（k21/s10000，与
> NWR 同参数）、skani ANI。对应 `design/genome-nn-query.md` §7.4 #2。

## 结果

### 近缘对（E. coli vs E. coli，完整 ANI ≈ 98.07）

| 完整度 | HV 距离 | Mash 距离 | skani ANI |
|---|---|---|---|
| 100% | 0.0618 | 0.0243 | 98.07 |
| 90% | 0.0653 | 0.0280 | 98.01 |
| 70% | 0.0714 | 0.0336 | 98.03 |
| 50% | 0.0883 | 0.0448 | 98.00 |

距离膨胀：HV +43%（相对），Mash +84%；ANI 几乎不变（−0.07）。

### 中等对（E. coli vs E. albertii，完整 ANI ≈ 90.02）

| 完整度 | HV 距离 | Mash 距离 | skani ANI |
|---|---|---|---|
| 100% | 0.1081 | 0.0867 | 90.02 |
| 90% | 0.1082 | 0.0867 | 90.02 |
| 70% | 0.1082 | 0.0867 | 90.02 |
| 50% | 0.1086 | 0.0876 | 90.02 |

三种距离全部稳定（<1% 变化）。

### 远缘对（E. coli vs Yersinia enterocolitica）

| 完整度 | HV 距离 | Mash 距离 | skani ANI |
|---|---|---|---|
| 100% | 0.2073 | 0.2607 | 无（skani 无比对） |
| 90% | 0.2023 | 0.2680 | 无 |
| 70% | 0.2688 | 0.2767 | 无 |
| 50% | 0.2317 | 0.2800 | 无 |

Mash 单调膨胀（+7%）；HV 非单调、噪声大；skani 在远缘（ANI<85% 且
截断后）无法给出 ANI——与 `design/hv.md` 的覆盖率一致。

## 结论

1. **完整度主要影响 k-mer 计数型距离（HV/Mash），不影响 ANI 本身**：
   近缘对 50% 完整度时 Mash/HV 距离显著膨胀（+84%/+43%），而 skani
   ANI 稳定在 98.0。机理 = 不完整基因组丢失共享 k-mer，Jaccard 低估 →
   距离高估。
2. **敏感性随亲缘距离非单调**：中等对（90% ANI）几乎不受影响（共享
   k-mer 富余），近缘对受影响最大——与"近缘分辨率本就差"叠加，完整度
   差异会在物种内聚类中制造伪距离（不完整株显得更远）。
3. **对设计的含义**：用 Mash/HV 做物种内聚类前应过滤/标记低完整度
   基因组（GSearch 建议完整度 >50%，我们的数据支持该阈值——50% 时
   近缘距离已 +43%）；ANI 层面 skani 对完整度稳健，标定不受影响。
4. 局限：随机删 contig 模拟，真实 MAG 缺失可能有偏（如偏向低覆盖区）；
   后续可换"删最大 contig 子集"或读段层面模拟复核。

## 复现

```bash
# 脚本: /tmp/hv_calib/run_completeness.py（读 gz FASTA → 随机保留
# contig → pgr dist hv / mash sketch+dist / skani dist）
```

---

## 证据附录：真实数据第二批验证：frac CI 校准、groups 分组、ANI 阈值、树一致性

> 日期：2026-08-08。数据与 cohort 同 `design/hv.md`
> （135 基因组，skani ANI 全矩阵）。对应 `design/genome-nn-query.md`
> §7.4 的 #5/#13/#15/#18。

## #5 frac 的 ANI 置信区间 vs skani ANI（覆盖率）

方法：`pgr dist frac --merge --ci` 全 cohort 输出每对 ANI 95% CI
（0–1 标度），与 skani ANI（0–100）对比，CI 乘 100 后计算覆盖率。

| ANI 区间 | n | CI 覆盖率 | CI 宽度中位数 | \|CI 中心 − skani ANI\| 中位数 |
|---|---|---|---|---|
| ≥99% | 122 | 15.6% | 0.10 | 0.21 |
| 95–99% | 1,339 | 21.4% | 0.17 | 0.23 |
| 90–95% | 2,427 | 1.4% | 0.51 | 0.92 |
| 85–90% | 721 | 0.0% | 0.67 | 1.27 |
| <85% | 1,328 | 11.8% | 2.42 | 1.91 |
| 全部 | 5,937 | **8.4%** | 0.53 | 0.96 |

结论：frac 的 CI 覆盖的是 **FracMinHash 估计自身的采样误差**（Hera et
al. 2023 公式），不能当作"覆盖金标准 ANI"的区间——frac 与 skani 的
系统偏差（近缘 ~0.2、远缘 ~1.9 ANI 单位）远大于 CI 宽度。文档/帮助文本
应注明 CI 语义，避免用户误读。

## #13 NWR groups.tsv（minhash 树 height 0.4 分组）与物种一致性

groups.tsv 覆盖 680 个代表基因组，**只有 13 个组**（组大小中位数 34，
最大 169）——在 Enterobacterales + Pasteurellales 全数据集上 height
0.4 切出的是**科/目级大组**，不是物种级 clade。

- 物种纯度：中位数 0.029（几乎每组都混合几十个物种）；属纯度中位数
  0.472（31% 的组单属）。
- mash.dist 组内距离中位数 0.197 vs 组间 0.284——区分度弱。

结论：groups.tsv 不能直接当物种级路由键；#7 的真实 clade 实验要用
更低的切割高度（或物种标签）在 cohort 自身距离上重新聚。

## #15 ANI 物种阈值（物种标签为真值的实证）

skani ANI 按物种标签分层（排除 "sp."）：

| 关系 | n | ANI q05 / q50 / q95 |
|---|---|---|
| 同种 | 1,034 | 96.1 / 98.2 / 99.1 |
| 同属异种 | 3,827 | 83.4 / 91.0 / 97.8 |
| 跨属 | 1,076 | 80.6 / 81.1 / 82.6 |

阈值判定（同种 vs 异种）：

| ANI 阈值 | 异种误判为同种（<阈值占比） | 同种误伤（<阈值占比） |
|---|---:|---:|
| 90 | 25.4% | 0.0% |
| 92 | 68.1% | 0.2% |
| 95 | 88.6% | 0.8% |
| 97 | 93.2% | 19.9% |

结论：~95% 阈值与经典物种边界一致（95% 阈值下同种误伤仅 0.8%、
异种漏判 11.4%）；90–92% 更保守。同种内 ANI 下限可达 91.7（个别
异常株/标注噪声），建议 CLI/文档以 95% 为默认"同种"判据、90% 为
"近似同种"提示。

## #18 minhash 距离树 vs bac120 标记基因树（物种级 cophenetic）

注意：两棵树 tip 命名体系不同（旧版 `Atl_hermannii_...` vs 新版
`At_herm_...`），需按 GCF accession 映射；bac120 树另有物种级
`Species___N` tip（`[S=...]` 注释使 necom 无法解析 condensed 版，
用 order 版 + accession 映射解决）。

- 共有物种：137（minhash 675 / bac120 146），物种对 9,316。
- cophenetic Spearman = **0.568**，Pearson = 0.919（深度分裂主导）。
- 分距离段：minhash 距离 0–0.05（近缘物种对）ρ=0.31；0.05–0.15
  ρ=0.39；0.15–0.4 ρ=0.71；>0.4 ρ=0.005。

结论：minhash 距离树与标记基因树在中远缘一致尚可，**近缘物种对的
排序一致性弱（ρ≈0.3–0.4）**——树层面再次印证 k-mer 距离在近缘下的
分辨率缺陷（与 `design/hv.md` 的距离层结论一致）。
bac120（保守标记基因）树更适合做物种级参考拓扑。

## 复现

```bash
# #5:  pgr dist frac cohort.fa.lst --list-files --merge --ci -o frac.ci.tsv
# #13: groups.tsv + genome.taxon.tsv + mash.dist.tsv（纯 python 分析）
# #15: ani.full.tsv + genome.taxon.tsv（纯 python 分析）
# #18: necom nwk label/distance + accession 映射（analyze_trees4.py）
```

---

## 证据附录：端到端验证：聚类 → 选参考 → 比对 → PBit 归档（2026-08-08）

对应 `design/genome-nn-query.md` §7.2 ③④⑤（P2）与核心工作流
"物种内聚类选参考 + PBit 归档"。目的：量化**参考选择策略**对
归档压缩率的影响，为实施方案提供决策依据。

## 数据与子集

- 数据源：`~/data/Escherichia/`（Enterobacterales NR，E. coli 全库
  51,318 组装；cohort = 2,088 个 E. coli NR 基因组，
  `/tmp/hv_calib/meta2115.tsv` + `mash2115.tsv` 全对距离现成）。
- Pilot：随机 30 个；全量：**farthest-point 采样 100 个**（基于
  mash2115.tsv 距离，seed 42，保证覆盖物种内多样性）。
- 所有样本 FASTA 路径校验存在（`ASSEMBLY/Escherichia_coli/<name>/`）。

## 方法

1. `pgr dist mash --merge --list-files` 全对距离（30 样本 435 对 /
   100 样本 4,950 对；0.01–0.04 mash 距离 = ANI ~96–99%）。
2. `necom mat to-phylip` + `necom clust upgma` + `necom cut simple`：
   30 样本 k=4（23/3/3/1）；100 样本 k=7（44/35/9/8/2/1/1）。
3. 每簇按三种策略各选 1 个参考：
   - **center**：簇内到其他成员 mash 距离和最小；
   - **longest**：FASTA 解压总长最大；
   - **random**：固定 seed 随机。
4. 簇内每个非参考样本 × 参考：
   `pgr align pgi`（PSL）→ `pgr psl to-paf` → `pgr pbit create`
   （纯 LZ 与 `--paf` CIGAR 两条路径）。
5. delta 压缩率 =（pbit 归档 − 参考 self-archive）/ gzip-9 样本大小；
   `pgr pbit to-fa` 抽查覆盖率。

## 结果

### Pilot（30 样本，26 查询 × 3 策略 = 78 对）

| 策略 | LZ delta/gzip（mean） | 备注 |
|---|---|---|
| center | 0.499 | 中心参考最差 |
| longest | **0.466** | 簇内最长参考最优 |
| random | 0.490 | 介于两者 |

### 全量（100 样本，98 查询 × 3 策略 = 279 对）

| 策略 | LZ delta/gzip mean | median | min | max |
|---|---|---|---|---|
| center | 0.5535 | 0.5505 | 0.4486 | 0.6327 |
| longest | **0.5198** | 0.5172 | 0.4334 | 0.6378 |
| random | 0.5213 | 0.5262 | 0.4302 | 0.6276 |

逐查询配对（n=98）：**longest < center 71/98**、random < center 61/98；
longest vs random 无稳定差异（39/49）。按簇看，longest 在 4/5 个
有样本的簇最优或接近最优；random 在 clade 1 偶然更优（0.496 vs
longest 0.532），说明**单参考的随机性影响可达 3–4 pp**。

### 参考特征（解释 longest 优势）

| 策略 | 参考 contig 数（中位） | 参考总长（平均） |
|---|---|---|
| center | 52 | 4.81 Mb |
| longest | **28** | **5.38 Mb** |
| random | 74 | 4.99 Mb |

center 参考（簇内距离和最小）倾向选"典型 draft"，内容覆盖量反而
最少；压缩率主要由**参考的内容覆盖量（总长 × 完整度）**决定，而非
与样本的平均相似度。

### 覆盖与 CIGAR 路径

- to-fa 覆盖率抽查（8/8）：100 版 ≥0.9984（KTE66 draft 极端 0.16%；
  其余 ≤0.02%），完整参考下多为 100%——LZ 内容匹配基本无损。
- CIGAR（`--paf`）与纯 LZ 归档大小差平均仅 32 B（max 143 B）：
  已知段相位约束使 CIGAR 路径基本回退 LZ（#14e），端到端无差异。

### Complete vs draft 参考（配对对照，2026-08-08 补充）

用户提出生产规则："选作参考的应该必须是 Complete"。本 cohort
（farthest-point 100 个 E. coli）中组装级别分布：Contig 48 /
Scaffold 27 / Complete Genome 22 / Chromosome 3（仅 25% 达标），
7 个簇中有 **2 簇完全没有 Complete/Chromosome 成员、1 簇仅 1 个**。

配对实验（控制"簇内最长"变量，只变组装级别）：clade 0/1 各选
complete-longest（48/GF60，Complete Genome，5–9 contigs 含质粒）与
draft-longest（49832/531，Contig，~200 contigs）两个参考，对同一批
75 个查询做 pbit（LZ）：

| 参考 | delta/gzip（mean） | complete 更优的对数 |
|---|---|---|
| complete-longest | 0.5249 | 19/75 |
| draft-longest | **0.5049** | —（draft 更优 56/75） |

即**压缩率维度 Complete 参考不占优**（draft-longest 反而 -2 pp）。
机制假说：① 查询 75% 是 draft，其 contig 颗粒与 draft 参考的 contig
天然对齐（LZ 内容匹配的 `best_ref_group` 命中率高）；② draft 参考含
质粒/未装配片段，与 draft 查询共享更多内容；③ Complete 参考的 4 kb
参考段是染色体上切出的，与查询 contig 边界错位。覆盖率两种参考都
≈100%（抽查），差异在编码效率而非覆盖。

**方案含义**：Complete 参考的正当理由不在压缩率，而在比对质量、
坐标可解释性与下游一致性；生产规则应定为"参考优先 Complete
（Chromosome 次之），若簇内无 Complete 则用簇内最长 draft（压缩率
不差）或并入相邻簇"，而不是"必须 Complete"。

## 结论与决策建议

1. **参考选择影响压缩率 ~3.4 pp（longest vs center）**；在 E. coli
   物种内（ANI 96–99%）最长/高完整度参考是稳定优选的简单启发式。
2. **不要用"距离中心"当参考**：中心参考是典型 draft，内容覆盖少，
   压缩率反而最差。压缩率维度上"最长 draft" ≥ "Complete"（配对
   实验，见上节）——**完整度不是压缩率的充分条件**；但生产规则仍
   建议 Complete 优先（比对/可解释性/一致性理由，见上节方案含义）。
3. 单参考随机性影响 3–4 pp：小簇（<10 成员）时参考选择比簇划分更
   影响结果；多参考（每簇 2–3 个）可摊平随机性（待验证）。
4. 端到端流程（dist → necom → pgi → pbit）全部就绪、耗时可控：
   100 样本 × 3 参考策略 ≈ 12 分钟（8 线程），pgi align ~1.5 s/对、
   pbit ~3.5 s/样本。

## 复现

```bash
# 子集选择（farthest-point，seed 42）
cd /tmp/e2e   # 脚本与中间产物
python3 - <<'PY'   # 见 bench 文档注释；cohort100.tsv/list100.txt 已生成
PY
pgr dist mash --merge --list-files -p 8 list100.txt -o dist100.pair
necom mat to-phylip dist100.pair.clean -o dist100.phylip
necom clust upgma dist100.phylip -o tree100.nwk
necom cut simple tree100.nwk -k 7 -o clust100_k7.tsv
python3 run_e2e_100.py   # 279 任务：align + pbit（LZ/CIGAR）
```
