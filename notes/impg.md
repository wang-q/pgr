# impg 分析笔记

本文档旨在总结 `impg` (implicit pangenome graph) 项目的核心设计、 命令分发架构与关键数据结构，
重点剖析 `src/main.rs` 的组织方式，为 `pgr`项目在泛基因组比对、区间投影与图构建管道方面提供参考。

`impg` 由 Andrea Guarracino、Bryce Kille、Erik Garrison 开发（v0.4.1，Rust 2021 edition），
其核心思想是：**不显式构建泛基因组图**，而是把 all-vs-all 两两比对当作一张"隐式图"，
通过区间树 (coitrees) 在比对网络上投影目标区间，按需提取同源序列。

`pgr` 自身的 pairwise 比对已成熟，并已实现多基因组共享 core 部分的比对； 本笔记主要关注
impg 在**泛基因组部分**（隐式图模型、按需投影、GFA 构建管道、crush）的设计，作为 `pgr`
向泛基因组方向扩展的参考。`pgr` 不打算参考 impg 的 syng syncmer 免比对后端，本文档中该部分一律从略。

## 1. 项目概览

### 1.1 设计哲学：隐式泛基因组图

传统泛基因组工具（pggb、Minigraph-Cactus）走"先物化 GFA，再在图上分析"的路线：
all-vs-all 比对 →seqwish 诱导图 → smoothxg 平滑 → 全图分析。这种路线的代价是：
**即使用户只关心一个位点，也要先 构建整张图**。对于 HLA、C4 这类只想看一个 locus 的场景，
等于把全基因组的图构建开销强加给单点查询。

`impg` 的核心洞察是：**all-vs-all 两两比对本身就是一张图的隐式描述** —— 比对双方的坐标区间是"边"，
序列是"节点"。既然比对文件已经存在，为什么还要物化一遍？于是走第三条路：

- **比对即图**：把 PAF/1ALN/TPA 两两比对视为图的隐式描述，不预先物化 GFA。
- **按需投影**：给定目标区间 `seq:start-end`，在区间树上查找所有重叠比对，把目标坐标 lift 到查询
  序列上，输出 BED/BEDPE/PAF/FASTA/GFA/VCF/MAF。
- **传递闭包**（§4.2）：`-x` 选项递归地把初次结果当作新的查询，找全所有同源片段。
- **不物化图**：除非用户显式要求 `gfa`/`vcf` 输出，否则永远不落到 GFA 文件。

#### 1.1.1 传递闭包：impg 区别于传统 UCSC 工具的核心

**这是 impg 真正的创新点，而非"隐式图"本身。**单条 pairwise 比对只给出 A↔B 的同源关系；但生物学
同源具有传递性：若 A~B、B~C，则 A~C（在 paralog 横跳、基因转换等场景下尤其重要）。要找全 A 的所有
同源片段，必须沿比对网络做 BFS，把间接同源也找出来。

- **UCSC Chain/Net 体系**：每条 Chain 是单条 A↔B 线性映射。`pgr chain lift` 把目标区间沿一条 Chain
  投影到查询序列，本质是单链线性投影。要找"所有同源片段"需要用户手动挑选一组 Chain，且不处理
  "通过第三序列间接同源"的情况。Net 在一定程度上编码了 syntenic 关系，但仍不是图遍历。
- **impg `-x`**：把所有 pairwise 比对当作图的边集，从目标区间出发做 BFS/DFS，自动发现所有直接和
  间接同源片段。`-m` 控制深度，`-d` 控制 hop 内最大 gap（也即单跳能吸收的最大 SV 长度）。
- **Cactus 的 transitive alignment**：Cactus 也做传递比对，但它在构建 Cactus 图时**一次性物化**全部
  传递关系。impg 把传递闭包**延迟到查询时才计算**，按需展开。

**代价**：impg 的传递闭包**不等于**真正的多序列比对。它是"在 pairwise 比对网络上做可达性查询"，
每个 hop 都依赖原比对的质量。若 A~B、B~C 但 A、C 在 B 的同一区段上不兼容（重叠但不一致），
impg 仍会报告 A~C 同源，但不会给出三者的多重比对——后者需要再走 `gfa`/`maf` 输出路径做局部 MSA。
这是"图遍历 vs 多序列比对"的本质区别，pgr 借鉴时需要清醒认识。

#### 1.1.2 隐式图 vs 物化图的适用边界

| 维度     | 隐式图（impg）赢                  | 物化图（pggb/odgi）赢                       |
|----------|-----------------------------------|---------------------------------------------|
| 查询模式 | 单 locus / 窗口化分区 / 动态子集  | 全图统计（coverage、layout、id stats）      |
| 重复访问 | 一次性查询、稀疏查询              | 反复遍历整图（genotyping、variant calling） |
| 内存     | 比对索引（可分片）+ 区间树        | 全图节点 + 路径 + 边                        |
| 预处理   | all-vs-all 比对（O(n²) 不可避免） | all-vs-all 比对 + seqwish + smooth          |
| 表达能力 | 同源关系（无 global 不变量）      | 路径保真、bubble 结构、coverage 向量        |

**关键认识**：impg **没有摆脱 O(n²) all-vs-all 比对**——它只是把"图构建"延迟到了查询时。隐式图的
真正优势不在"省比对"，而在"按需计算图遍历"：单 locus 查询只需扫一遍区间树，不必构建/平滑整张图。
对于 pgr 的 UCSC Chain/Net 体系，**Chain 本身就是隐式图的边**：一条 Chain = 一条 A↔B 边集。pgr
若要做 cohort 级隐式图，无需引入 PAF/wfmash，可直接在 Chain/Net 上构建区间索引，复用已有的 pairwise
成熟基础设施。

#### 1.1.3 能力栈：从索引到基因分型的四层

impg 的 22 个子命令并非平级，而是按**数据流方向**构成四层能力栈。理解这个层次是阅读后续章节的 地图：

```
索引层（前提）      Index — 把上游比对文件全量装入 coitrees（不过滤）
    ↓
查询层              Query / Project — 区间投影（单点）↔ 传递查询（-x BFS 闭包）
                    过滤参数在此层：-d/--merge-distance、--min-result-identity、
                    -l/--min-output-length、--subset-sequence-list、-m/--max-depth
    ↓
图构建层            Graph — 直接构建 / Partition + Lace — 分区构建（窗口化控内存）
    │               Crush — bubble 压缩（图构建后处理）
    │               GraphReport — 图质量评估 / Render — 渲染
    ↓
应用层              Genotype / Infer — 基因分型与等位基因推断
```

**关键认识**：

- **索引是前提且不过滤** — 所有后续能力都建立在"比对已装入区间树"之上。`Index` 命令
  （[`AlignmentOpts`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L4067)）
  只有文件 路径、`--index-mode`、`--unidirectional`、`--trace-spacing`
  等参数，**没有** min_aln/min_match_len/min_mapq 等过滤开关；
  [`Impg::from_multi_alignment_records`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L11561)
  把输入 PAF/1ALN/TPA **全量**装入。比对质量的挑选责任在生成 PAF 的上游工具（wfmash/sweepga），
  不在 impg——这是"比对即图"哲学的体现：索引保留所有边，挑选推迟到查询。
- **挑选发生在查询层** — 真正的"挑选比对"由
  [`QueryOpts`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L4319)控制：
  `-d/--merge-distance`（必填或 `--no-merge`，合并间隔 ≤ D bp 的区间，D 也是单跳能吸收的 最大
  gap/SV）、`--min-result-identity`（最低 gap-compressed identity）、`-l/--min-output-length`、
  `--subset-sequence-list`（只保留指定序列）；传递闭包专属 `-m/--max-depth`（默认 2）、
  `--min-transitive-len`（默认 101）、`--min-distance-between-ranges`（默认 10）。同一份索引可
  服务不同严格度的查询。
- **查询层的两种模式** — 区间投影（单点查找重叠比对）与传递查询（BFS 找全所有同源片段）不是
  平级命令，而是同一个 `Query` 的两种模式（`-x` 开关）。后者是前者的扩展。
- **图构建层的两种形态** — 直接构建（`Graph` 单图）与分区构建（`Partition` 切窗口 + `Lace` 拼回）
  是同一层的两种内存策略，不是独立步骤。`Crush` 是图构建的后处理（压缩 bubble），不是 独立的图构建路径。
  图构建层另有 [`SeqwishOpts`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L2063)的
  `--min-match-len`（默认 23）/`--repeat-max`/`--sparse-factor` 控制传递闭包的比对挑选。
- **应用层依赖图构建** — 基因分型需要图或隐式图作为 evidence 来源，是能力栈的顶端消费者。

pgr 的当前能力覆盖索引层（Chain/PSL 处理）与查询层（`chain lift` 单链投影），但
**图构建层与 应用层是空白**——这正是 §9 启示聚焦的方向。

#### 1.1.4 关键概念名词解释

能力栈涉及的核心术语在此集中定义，后续章节直接引用：

- **隐式泛基因组图（implicit pangenome graph）** — 不显式构造 GFA，把 all-vs-all pairwise 比对
  当作"图的边集"装入区间树，按需通过查询提取同源片段。impg 的核心哲学（§1.1.1）。
- **coitrees** — cache-oblivious interval trees，impg 索引层的物理基础。按 target sequence 分多棵 树（
  [`TreeMap = FxHashMap<u32, Arc<BasicCOITree<QueryMetadata, u32>>>`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/impg.rs#L226)），
  支持高效区间重叠查询。详见 §3.1。
- **区间投影（interval projection）** — 查询层的基础模式：给定目标区间，在区间树上查所有重叠 比对，
  把目标坐标投影到查询序列坐标。等价于"单跳同源查找"（§4.1）。
- **传递闭包（transitive closure）** — 查询层的高阶模式（`-x`）：把 pairwise 比对当图边，从目标
  区间出发做 BFS/DFS，自动发现所有直接和间接同源片段。是 impg 区别于传统 UCSC 工具的核心创新
  （§4.2）。**注意：传递闭包是图遍历，不是多序列比对**（§4.3）。
- **单跳（one hop）** — 传递闭包遍历的一步：从一条序列经一条比对边跳到另一条序列。
  `-d/--merge-distance`限制单跳能吸收的最大 gap/SV 长度。
- **locus（复数 loci）** — 基因组上的一段区域。在 impg 语境下，partition 产出的 locus 是
  **跨多条 序列的同源区段集合**，不是单序列的一段坐标。
- **partition（分区）** — 把整个 cohort 的基因组切成一组互不重叠的 loci，每个 locus 独立处理。
  算法是"批量传递闭包 + masking 去重"，**不是**简单按窗口切分（§6.3）。
- **bubble（气泡）** — 泛基因组图中的变异结构：从一个节点分叉出多条路径再汇合。crush 算法的 处理对象
  （§6.5）。
- **crush** — 图构建后处理：把 seqwish 产出的"bubble 状"结构压缩成更紧凑的变异图。8 阶段流程， 15
  种 `ResolutionMethod`（§6.5）。
- **stage 化管道** — 用冒号分隔的字符串（如 `gfa:cut-n=100:pggb:crush:sort`）描述图构建流水线， 由
  `graph_pipeline.rs` 解析（§6.2、§2.4）。
- **syng 后端** — 免比对后端：从 FASTA/AGC 构建 syncmer GBWT 索引，用 syncmer anchor 链定义同源，
  不跑 wfmash/FastGA。与 alignment 后端（PAF/1ALN/TPA）并列（§1.2、§5）。
- **GFA engine** — 图构建引擎枚举
  （[`GfaEngine`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/lib.rs#L38)）：`Pggb`/
  `Seqwish`/`Poa`/`SyngNative`/`SyngLocal` 5 种，决定基础管道形态（§6.1）。
- **pair-selection** — 避免 all-vs-all 比对的机制：用 Mash KNN（`--sparsify`）或 syng anchor counts
  选对比对，把 N² 降到 N×K（§6.4）。
- **graph feature evidence** — 基因分型的后端中立抽象：候选 haplotype 编码为特征向量（节点/边/ 路径
  presence），既可来自物化 GFA，也可来自隐式图查询（§7.1）。

### 1.2 双后端架构

`impg` 支持两种对齐数据来源，通过同一个 `-a` 参数透明切换：

| 后端               | 输入                        | 工作方式                     | 典型用途        |
|--------------------|-----------------------------|------------------------------|-----------------|
| **Alignment 后端** | PAF/1ALN/TPA                | coitrees 索引，按 CIGAR 投影 | all-vs-all 比对 |
| **syng 后端**      | syng 索引前缀（多 sidecar） | syncmer GBWT 免比对          | FASTA/AGC 构建  |

`main.rs` 中的 [`detect_syng_prefix`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L4604)
函数会根据 `-a` 参数的扩展名自动判定走哪个后端。

### 1.3 源代码组织

`impg` 把库逻辑放在 `src/lib.rs`，把命令行入口与参数解析放在 `src/main.rs`。全 `src/` 共 44 个
`.rs` 文件、约 8.5 万行；其中 `main.rs` 1.57 万行、`resolution.rs` 1.72 万行、`syng.rs` 0.9 万行
三个超大文件占了近一半（这是 impg 最显著的工程反例，详见 §2.1）。`lib.rs` 导出 31 个模块，按
**职责层次**（而非文件清单）可分为四层：

- **索引层（隐式图核心）** — `Impg` struct（单文件 `.impg` 索引）+ `ImpgIndex` trait（统一单/多
  文件接口）+ `MultiImpg`（协调 per-file 子索引）+ `ForestMap`（target_id→tree_offset 反向索引）。
  这是"比对即图"哲学的物理载体：把 all-vs-all 比对装进 coitrees，查询时在区间树上投影。
- **格式抽象层** — `AlignmentRecord` 把 PAF/1ALN/TPA 三格式统一成 8 字段（strand + 文件偏移打包），
  使上层索引与查询代码不感知具体比对格式。`SequenceIndex` + `PanSn`（`sample#haplotype#contig`）
  提供序列命名抽象。
- **图构建层（泛基因组核心）** — `resolution.rs`（crush 算法实体，§6.5）+ `smooth.rs`（smoothxg
  风格块分解 + SPOA 平滑）+ `graph_report.rs`（15+ 维度图质量报告，为 sweep 脚本提供评分依据）+
  `graph_pipeline.rs`（stage 化管道解析）。这一层是 impg 区别于纯投影工具的关键：不仅查询图，
  还能物化图、压缩图、评估图。
- **命令层（`commands/`）** — 各子命令的业务逻辑。`main.rs` 只做参数装配与分发，业务逻辑下沉到
  `commands::<sub>::run`（详见 §2.3 的"瘦分发 + 胖模块"模式）。

**层次依赖方向**：命令层 → 图构建层 → 索引层 → 格式抽象层。syng 后端（6 文件）作为平行分支，
共享索引层与命令层但不进入图构建层（本文档不参考）。

### 1.4 关键依赖

impg 的依赖反映其能力栈，按**能力维度**可分为五组：

- **区间查询能力** — `coitrees`（cache-oblivious interval trees）是"比对即图"哲学的物理基础， 把
  all-vs-all 比对装进内存区间树。这是 pgr 最值得直接借鉴的依赖。
- **图操作能力** — `gfa` + `handlegraph` 提供 GFA 读写与 handlegraph 抽象，支撑 lace、crush 等
  图操作。impg 把图当一等公民处理，但**只在需要时物化**（隐式图路线）。
- **比对能力（aligner 路由）** — `lib_wfa2`（WFA 仿射比对，BiWFA 边界精修）+ `spoa_rs`/`poasta`
  （POA MSA 引擎，用于 `gfa:poa` 与 crush 阶段）。crush 的 median 三档路由（§6.5.2）本质上是 在这些
  aligner 之间做选择。
- **工具链集成（同团队生态）** — `sweepga`/`seqwish`/`allwave`/`gfasort`/`bluntg`/`povu` 都是 Erik
  Garrison 团队维护的泛基因组工具，**作为库直接嵌入而非子进程调用**。这种"工具链库化"模式使 impg
  能深度集成各阶段，但也把外部工具的复杂性内化到二进制中。pgr 若引入类似生态，应评估"库化 vs
  子进程"的权衡。
- **基础设施** — `rayon`/`crossbeam-channel`/`indicatif`（并行/通道/进度条）+ `noodles`/
  `rust-htslib`（BIO 格式/BAM/CRAM）+ `ragc-core`（AGC 序列归档）+ `onecode`/`tpa`/`tracepoints`
  （1ALN/TPA 格式与 tracepoint 编解码）。

### 1.5 主要子 crate 作用清单

§1.4 按能力维度分组，本节列出每个 git 依赖的具体作用与调用方式，便于评估 pgr 是否需要引入等价物。
impg 的子 crate 全部来自 Erik Garrison/Andrea Guarracino 团队的 pangenome 生态，**除 gfaffix 外
均作为库直接嵌入**（与 §1.4 "工具链库化"一致）。

15 个 git/vendor 依赖按能力分组（前 14 个库嵌入，最后 1 个子进程）：

- **图构建生态（6 个，pangenome 团队）**
    - **sweepga**（pangenome/sweepga）— wfmash/FastGA 集成 + KNN sparsification + PAF filter。
      模块：[`pansn`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/lib.rs#L875)、
      [`knn_graph`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/syng_graph.rs#L600)、
      [`paf_filter`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/syng_graph.rs#L825)。
      调用位置：`graph`/`align`/`syng2gfa`。
    - **seqwish**（pangenome/seqwish）— 从 PAF 诱导 GFA。
      入口：[`generate_gfa_seqwish_from_intervals`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/graph.rs#L1036)。
      调用位置：`graph` 的 Seqwish engine。
    - **allwave**（pangenome/allwave）— crush 的
      [`Allwave` ResolutionMethod](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/resolution.rs#L8685)
      + syng-native BiWFA pair sparsification。调用位置：`crush`/`syng`。
    - **gfasort**（pangenome/gfasort）— GFA 排序（Ygs pipeline：path-guided SGD + grooming + 拓扑排序）。
      调用位置：`graph` 的 sort 阶段（[graph.rs#L423](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/commands/graph.rs#L423)）。
    - **bluntg**（pangenome/bluntg）— GFA bluntify（把 link/path overlap 转 0M）。
      调用位置：`syng2gfa --gfa-mode blunt`（[syng2gfa.rs#L3660](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/commands/syng2gfa.rs#L3660)）。
    - **povu**（povu-rs，pangenome/povu）— GFA 解析（`NativeGfa`）+ flubble 检测。
      调用位置：`graph report`（[graph_report.rs#L2135](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/graph_report.rs#L2135)）。
- **比对能力（3 个）**
    - **lib_wfa2**（ekg/lib_wfa2）— WFA 仿射比对（crush 边界精修 + syng-native BiWFA）。
      调用位置：`crush`/`syng`。
    - **spoa_rs**（AndreaGuarracino/spoa-rs）— SPOA POA MSA 引擎。
      调用位置：`gfa:poa` + crush 阶段。
    - **poasta**（crates.io）— POASTA POA MSA 引擎（新算法）。调用位置：crush 路由选项。
- **图抽象（1 个）**
    - **handlegraph**（chfi/rs-handlegraph）— handlegraph 抽象（handle/edge/path 统一接口）。
      调用位置：`lace`/`crush` 图操作基础。
- **格式支持（4 个）**
    - **ragc-core**（AndreaGuarracino/ragc）— AGC 序列归档格式（syng 后端输入）。
    - **onecode**（pangenome/onecode-rs）— 1ALN 比对格式解析（AlignmentRecord 后端之一）。
    - **tpa**（AndreaGuarracino/tpa）— TPA 比对格式解析（AlignmentRecord 后端之一）。
    - **tracepoints**（AndreaGuarracino/tracepoints）— 1ALN/TPA tracepoint 编解码。
- **子进程调用（1 个）**
    - **gfaffix**（vendor/gfaffix）— **唯一子进程调用**（其余 14 个都是库嵌入）。
      通过 [`run_gfaffix`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/graph.rs#L974) 用
      `current_exe().with_file_name("gfaffix")` 定位 sibling 二进制。
      作用：GFA 归一化（图构建后的 normalize 阶段）。

**关键认识**：

- **库化 vs 子进程的不一致** — 14 个 git 依赖中 13 个库化嵌入，唯独 gfaffix 走子进程
  （[`run_gfaffix`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/graph.rs#L974) 定位 sibling
  二进制）。原因可能是 gfaffix 用了不同的 handlegraph 版本，库化会引入 API 冲突。pgr 若引入类似
  生态，应优先全库化，避免二进制分发复杂度。
- **sweepga 是最重的依赖** — 它同时提供 aligner 集成（wfmash/FastGA）、sparsification（KNN graph）、
  PAF filter、PanSN 命名抽象，是 impg 图构建层的核心。pgr 若借鉴 pair-selection，sweepga 的
  `knn_graph` 模块是参考点（impg.md §6.4）。
- **pgr 可暂不引入整个生态** — pgr 第一步聚焦 Chain 传递闭包（pairwise-selection.md §3.1），
  只需 coitrees 等价物（区间树）。sweepga/seqwish/gfasort 等图构建生态是后续阶段的需求。

## 2. main.rs — 命令分发与参数解析（重点）

`src/main.rs` 是 `impg` 二进制的入口，单文件约**613 KB / 1.5 万行+**， 承担了所有 clap 命令定义、
参数解析、stage 解析、命令分发与大量辅助逻辑。这是该项目最显著的工程特点（也是值得反思之处）。

### 2.1 文件顶层结构

`main.rs` 的内部组织遵循"辅助函数 → 命令定义 → 命令分发"的顺序，但比例严重失衡：

- **前 ~4500 行：辅助函数与 stage 解析器** — 大量 `parse_*` 工具函数（度量后缀、merge distance、
  round count）+ `apply_gfa_output_engine_shorthand` 及其一系列 `parse_*_stage` 子函数 （§2.4
  详述）。这部分是 `main.rs` 膨胀的主因——stage 化字符串 DSL 的解析逻辑全塞在这里。
- **~L4707：顶层 `Args` enum** — 所有 22 个子命令的 clap 定义集中在此（§2.2）。
- **~L6168：`run` 函数的巨型 dispatch** — `match args { ... }` 分发到各 `commands::*` 模块。
- **末尾：剩余辅助函数** — `build_graph_config`/`validate_*`/`initialize_threads_and_log` 等。

**组织逻辑反思**：impg 把"参数解析 + stage DSL 解析 + 命令分发"全塞进单文件，导致 `main.rs`
既是入口又是 DSL 解释器又是分发中心。pgr 的 `cmd_pgr/` 按格式/功能分组、每命令独立模块的结构
更优——若未来需要 stage DSL，应独立到 `libs/pipeline_dsl.rs` 而非塞进 `pgr.rs`。

### 2.2 顶层 `Args` enum 与子命令清单

[main.rs#L4707](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L4707) 处的
`#[derive(Parser)] enum Args` 定义了所有子命令。下面按功能归类：

- **索引**
    - `Index` — 构建 `.impg` 索引（single 或 per-file 模式）；逻辑在 `main.rs` 内联 +   `impg.rs`/
      `multi_impg.rs`。
- **查询/投影**
    - `Query` — 把目标区间通过比对网络投影（核心命令）；路由到 `impg_index.rs` 的 trait 方法。
    - `Project` (别名 `projection`) — 投影命令；`projection.rs` + `projection/converter.rs`。
    - `Map` — 把短读段映射到 syng 索引（GAF/PAF/pack/proj 输出）；syng 后端。
- **分割/拼接**
    - `Partition` — 把 cohort 按互不重叠的 locus 切分（传递闭包 + masking 去重，详见 §6.3）；
      `commands/partition.rs`（含 `rehome_singleton_slivers`）。
    - `Lace` — 把多个 per-window GFA/VCF 拼回一张图；`commands/lace.rs`（支持 gzip/bgzf/zstd）。
    - `Refine` — 收紧 loci 边界以最大化 sample/haplotype 支持度；`commands/refine.rs` （含
      transitive BFS/DFS、blacklist）。
- **分析**
    - `Similarity` — 区间内成对相似度/距离 + PCA/MDS；`commands/similarity.rs`
      （Jaccard/cosine/dice/estimated_identity）。
    - `Stats` — 汇总比对统计；`main.rs` 内联。
- **图构建**
    - `Graph` — 从 FASTA 直接构建泛基因组图（不需要预比对）；`commands/graph.rs` + `graph.rs` +
      `smooth.rs`。
    - `NormalizeSelfLoops` — 折叠 blunt GFA 中路径局部的 self-loop 重复单元； `gfa_self_loops.rs`。
    - `Crush` — 解析 bounded bubble，支持多种 POA/POASTA/sweepga 路由； `resolution.rs`（§6.5，15
      种 `ResolutionMethod`）。
    - `Gfa2vcf` (别名 `gfa-to-vcf`, `povu`) — GFA → VCF；`main.rs` 内联 + povu 库。
    - `GraphReport` (别名 `describe-graph`, `describe-gfa`) — 输出图特征 Markdown/JSON/TSV 报告；
      `graph_report.rs`（15+ 维度）。
    - `Syng2gfa` — syng 索引物化为 GFA（`--gfa-mode blunt|raw`，§6.4）；`commands/syng2gfa.rs`
      （`SyngGfaMode::Blunt` 走 bluntg，`Raw` 输出 syng-native overlap 图）。
    - `Render` — 用 gfalook 渲染 1D 图；`commands/render.rs` + `render_bundle.rs`。
    - `Align` — 调用 wfmash/FastGA 跑比对；`commands/align.rs` +
      `commands/mod.rs::create_aligner_adaptive`。
- **syng 后端**
    - `Syng` — 从 FASTA/AGC 构建 syng 索引；`syng.rs` + `syng_parallel.rs`。
    - `SyngRepair` — 重建 `.syng.pstep`/`.syng.spos` 而不重读序列；`syng.rs`。
    - `ReadIndex` — 构建 read-to-syncmer 倒排索引（输出 `.r2s.*` sidecar）；syng 后端。
- **基因分型**
    - `Genotype` (别名 `gt`) — 基因分型命名空间，含 `Cos` (cosigt) 子命令； `commands/genotype.rs` +
      `genotyping.rs`（`SyngCosigtConfig`/`CandidateMode`）。
    - `Infer` — 跨区间/分区输出等位基因 call，支持 stitching；`commands/infer.rs` （`InferTarget`/
      `PartitionDiscoveryConfig`）。

`Query` 命令的 `after_help` 文档
（[main.rs#L4892](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L4892)）
是整个项目最完整的输出格式说明，列出了 `auto/bed/bedpe/paf/gfa/vcf/maf/fasta/fasta+paf/fasta-aln/gbwt`
共 11 种输出与各自的约束。

### 2.3 命令分发模式 (`run` 函数)

[main.rs#L6168](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L6168) 处的
`fn run() -> io::Result<()>` 是命令分发中心。它的模式如下：

```rust
fn run() -> io::Result<()> {
    let args = Args::parse();
    match args {
        Args::Index { common, alignment, sequence } => {
            initialize_threads_and_log(&common);
            let alignment_files = resolve_alignment_files(&alignment)?;
            // ...
            initialize_index(...)?;
            info!("Index created successfully");
        }
        Args::Lace { common, sequence, files, file_list, format, output, compress,
                     fill_gaps, skip_validation, temp_dir, reference } => {
            initialize_threads_and_log(&common);
            // 参数验证 → 路由到 lace::run_gfa_lace 或 lace::run_vcf_lace
        }
        Args::Partition { common, alignment, syng_padding, /* ... 大量字段 */ } => {
            initialize_threads_and_log(&common);
            // 验证 → 解析 engine spec → 路由到 partition::run
        }
        // ... 其余 17 个分支
    }
}
```

每个分支的模式高度一致：

1. `initialize_threads_and_log(&common)` — 设置 rayon 线程池与日志级别。
2. 解析/验证参数（`validate_output_format`、`require_merge_distance`、 `validate_selection_mode`、
   `validate_region_size` 等）。
3. 通过 `engine_cli.parse_engine()` 解析 GFA 引擎简写（如 `pggb:10000`）。
4. 路由到 `commands::<sub>::run(...)` 中的实际实现。

`main.rs` 本身**不包含**子命令的业务逻辑（除了 `Index`/`Lace` 等较简单的命令），
它只负责参数装配与调用 `commands::*` 模块。这种"瘦分发 +胖模块"的边界是项目刻意维持的。

### 2.4 辅助函数分类

`main.rs` 顶部的约 4500 行是大量辅助函数，按职责可归为几类：

#### 参数解析 (`parse_*`)

```rust
fn parse_size(s: &str) -> Result<u64, String>            // 接受 k/m/g 后缀
fn parse_merge_distance(s: &str) -> Result<i32, String>  // i32 边界检查
fn parse_usize_size(s: &str) -> Result<usize, String>
fn parse_round_count(s: &str) -> Result<usize, String>   // "until-done" → usize::MAX
```

注意 `parse_merge_distance` 内部调用 `sweepga::parse_metric_number`， 作者注释说明这是为了"让 impg
和 sweepga 用同样的方式解释用户输入的距离后缀"——这是与依赖库保持一致性的设计选择。

#### GFA 引擎简写解析 (`apply_gfa_output_engine_shorthand`)

这是 `main.rs` 中最复杂的部分之一。用户可以写：

```
-o gfa:pggb                         # 等价于 --gfa-engine pggb
-o gfa:pggb:10000                   # partitioned 模式，窗口 10kb
-o gfa:cut-n=100:pggb               # 终端 N-run 裁剪 + pggb
-o gfa:syng:blunt,k=63,s=8,seed=7   # syng 引擎 + 参数断言
-o gfa:syng:crush                   # syng + crush 阶段
-o gfa:pggb:crush,method=allwave,k-nearest=5
-o vcf:syng                         # syng 引擎 + VCF 输出
```

[main.rs#L164](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L164) 的
`apply_gfa_output_engine_shorthand` 拆解这些冒号分隔的 stage，并分派给一系列 `parse_*_stage` 函数：

| Stage 函数                    | 处理的内容                                                           |
|-------------------------------|----------------------------------------------------------------------|
| `parse_terminal_n_clip_stage` | `cut-n=<bp>` 终端 N-run 裁剪                                         |
| `parse_syng_mask_stage`       | syng syncmer mask 参数                                               |
| `parse_crush_stage`           | bubble crush (`method`/`k-nearest`/`pair-trees`/`polish-rounds` 等)  |
| `parse_smooth_stage`          | smoothxg 平滑 (`target-poa-length`/`max-node-length`/`block-source`) |
| `parse_graph_sort_stage`      | 最终 gfasort pipeline (默认 `Ygs`)                                   |
| `parse_syng_assertion_params` | syng `k/s/seed` 参数断言                                             |

每个 stage 解析器返回 `Option<String>` 表示剩余的 stage 串，串成一条管道。 这种"stage 化字符串
DSL"是 `impg` 的特色，但也使 `main.rs` 体积膨胀。

#### FASTA/FASTQ 流式读取

`read_fasta_records`、`read_fastq_records`、`stream_fasta_query_chunks`、
`stream_fastq_query_chunks`、`stream_query_chunks_with` 等函数为 `Map`命令提供流式读序支持，
避免一次性把所有读段加载到内存。

#### 进程级 stdout 静默

[main.rs#L1883](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L1883) 的
`silence_stdout_for_process`（unix 分支）通过 `RawFd` 重定向屏蔽子进程的 stdout，
用于 `gfaffix`、`odgi` 等外部工具的调用。这是 unix-only 代码，非 unix 平台有空实现（
[main.rs#L1926](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L1926)）。

### 2.5 测试入口

[main.rs#L6151](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L6151) 的
`fn args_command_for_test() -> clap::Command` 暴露完整的命令定义给集成测试，用于断言 help
文本与参数互斥规则。

## 3. 索引层：隐式图核心数据结构

本层是 impg"比对即图"哲学的物理载体——把 all-vs-all 比对装进 coitrees，查询时在区间树上投影。
下列数据结构支撑了 §1.1.3 能力栈的"索引层"，是后续查询层（§4）与图构建层（§6）的前提。

### 3.1 `Impg` 结构

[impg.rs#L394](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/impg.rs#L394)
定义了单文件索引的核心结构：

```rust
pub struct Impg {
    pub trees: RwLock<TreeMap>,                          // target_id -> COITree
    pub seq_index: SequenceIndex,                        // 序列名 <-> ID 映射
    alignment_files: Vec<String>,                        // PAF/1ALN/TPA 文件路径
    pub forest_map: ForestMap,                           // 反向索引
    index_file_path: String,
    pub sequence_files: Vec<String>,
    trace_spacing_cache: RwLock<Vec<Option<i64>>>,       // .1aln/.tpa trace_spacing 懒加载
}
```

其中 `TreeMap = FxHashMap<u32, Arc<BasicCOITree<QueryMetadata, u32>>>` — 每个 target 序列一棵区间树，
节点 metadata 是 `QueryMetadata`。`RwLock` + `Arc`允许查询时无锁共享。

### 3.2 `QueryMetadata` 与 `CigarOp` 紧凑编码

[impg.rs#L165](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/impg.rs#L165) 的 `QueryMetadata`
用 bit-packing 压缩 metadata：

```rust
pub struct QueryMetadata {
    query_id: u32,
    target_start: i32,
    target_end: i32,
    query_start: i32,
    query_end: i32,
    alignment_file_index: u32,
    strand_and_data_offset: u64,   // bit 63 = strand, bit 62 = reversed, 低位 = data offset
    data_bytes: usize,
}
```

`CigarOp` ([impg.rs#L74](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/impg.rs#L74))
把 CIGAR 操作 + 长度压进单个 `u32`：高 3 位是 op (`=`/`X`/`I`/`D`/`M`)，低 29 位是长度。
这种设计使区间树节点极紧凑，能装下全基因组规模的 all-vs-all 比对。

`Impg` 通过 `get_cigar_ops` 按需还原 CIGAR：

- **PAF**：直接从原文件 `data_offset` 处读取 `cg:Z:` 标签的字节。
- **1ALN/TPA**：从 tracepoint 解码，需要 `trace_spacing` 与（可能的）目标序列做 BiWFA 还原。

### 3.3 `SortedRanges` 与区间合并

[impg.rs#L243](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/impg.rs#L243) 的
`SortedRanges` 维护按起点排序的区间集合，支持基于 `min_distance`的合并。`insert`
方法返回"未被现有区间覆盖的新增部分"，是传递闭包查询中"只把新发现区间加入下一轮 BFS/DFS"的关键。

### 3.4 `ImpgIndex` trait 与 `MultiImpg`

[impg_index.rs](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/impg_index.rs) 定义了 `ImpgIndex`
trait：

```rust
pub trait ImpgIndex: Send + Sync {
    fn seq_index(&self) -> &SequenceIndex;
    fn query(&self, ...) -> io::Result<Vec<AdjustedInterval>>;
    fn query_with_cache(&self, ..., cigar_cache: &FxHashMap<...>) -> io::Result<...>;
    fn populate_cigar_cache(&self, ...);
    fn query_transitive_dfs(&self, ...) -> io::Result<...>;
    fn query_transitive_bfs(&self, ...) -> io::Result<...>;
    // ...
}
```

`Impg`（单文件）与 `MultiImpg`（多文件）都实现这个 trait。`MultiImpg` 内部维护
`TreeLocation { index_idx, local_target_id }` 把全局 `target_id`翻译到子索引的本地 ID。`main.rs`
中的命令代码只与 `&dyn ImpgIndex` 打交道，从而对单/多文件透明。

`MultiImpg` 还实现了 staleness 检测：当 `.impg` 索引比源比对文件旧时， 警告并要求 `--force-reindex`。

<!-- crush 算法原 §3.5 内容已移至 §6 图构建层 -->
## 4. 查询层：区间投影与传递闭包

索引层（§3）把比对装进区间树后，查询层负责"在比对网络上回答同源问题"。impg 的 `query`/`project`
命令提供两种模式，本质是同一个区间树查找的两种展开方式。

### 4.1 区间投影（单点查询）

**区间投影**是查询层的基础模式：给定一个目标区间 `[target, start, end]`，在区间树上查找所有
与之重叠的比对记录，把目标区间投影到查询序列坐标。这等价于"单跳同源查找"——只找直接比对到
目标区间的片段，不沿比对网络遍历。

`pgr chain lift` 本质上就是单条 Chain 的区间投影：通过 Chain 的坐标映射把目标区间 lift 到查询 序列。
区别在于 impg 在 all-vs-all 比对网络的并集上做区间树查找，而 pgr chain lift 在用户手动 指定的单条
Chain 上做线性投影。

### 4.2 传递闭包（`-x` BFS 查询）

传递闭包的哲学意义与适用边界已在 §1.1.1 详述——这是 impg 区别于传统 UCSC 工具的核心创新。本节
聚焦查询层的技术实现。

`impg query -x` 把所有 pairwise 比对当作图的边集，从目标区间出发做 BFS/DFS，自动发现所有直接和
间接同源片段。关键参数（完整清单见 §1.1.3 的 `QueryOpts`）：`-d/--merge-distance`（必填，
控制 hop 内最大 gap，也即单跳能吸收的最大 SV 长度）、`-m/--max-depth`（BFS 深度，默认 2）、
`--min-transitive-len`（默认 101，小于此长度的区间不进入下一轮）、`--min-distance-between-ranges`
（默认 10，同序列上传递区间最小间距）。实现上依赖 §3.3 的 `SortedRanges`——每轮 BFS 只把"未被
现有区间覆盖的新增部分"加入下一轮，避免重复遍历。

**注意**：impg 的"挑选比对"全部发生在查询层，**不**在索引层（§1.1.3）。索引全量装入上游比对 文件，
由查询参数控制严格度。pgr 若做 `chain query`，应在查询参数而非索引参数上暴露过滤。

### 4.3 传递闭包 ≠ 多序列比对

§1.1.1 末尾已指出：impg 的传递闭包是"图遍历"而非"多序列比对"——它能找全同源片段，但不给出 多重比对，
后者需要再走 `gfa`/`maf` 输出路径做局部 MSA。

对 pgr 的意义：pgr 已有的 `fas multiz`（多基因组 core 比对）+ `fas consensus`（POA consensus）
恰好可以填补"从同源关系到多重比对"这一步——传递闭包找全片段，`fas` 系列做真正的 MSA。这是 pgr 相对
impg 的天然优势：impg 的 MSA 路径是 per-bubble POA，pgr 的 MSA 路径是成熟的 banded DP 合并。

### 4.4 深度对比：区间投影 vs Chain lift

| 维度     | impg 区间投影/传递闭包              | pgr chain lift                     |
|----------|-------------------------------------|------------------------------------|
| 数据源   | all-vs-all 比对网络（PAF/1ALN/TPA） | UCSC Chain（已 syntenic 净化）     |
| 查找方式 | 区间树查找 + BFS 传递闭包           | 单条 Chain 线性映射                |
| 同源发现 | 自动找全直接 + 间接同源             | 需用户手动选 Chain，不处理间接同源 |
| 代价     | O(n²) all-vs-all 比对预处理         | 直接用现成 Chain，无预处理         |
| 适用场景 | cohort 级"找全所有同源片段"         | 已知 syntenic 关系的坐标转换       |

**对 pgr 的启示**（详见 §9）：pgr 无需重新跑 all-vs-all，可直接在 Chain/Net 上构建区间索引，
复用 impg 的 BFS 传递闭包思路，实现 cohort 级"找全所有同源片段"——这是 pgr 当前 `chain lift`
（单链线性投影）的天然扩展。

## 5. syng 免比对后端（略）

`impg` 还有一个 syng syncmer GBWT 免比对后端（`SyngIndex` / `SyngMatcher` / `SyncmerParams`，6
个 sidecar 文件 `.1khash`/`.1gbwt`/`.syng.names`/`.syng.pstep`/`.syng.spos`/`.syng.meta`，以及
`impg map` 的 GAF/PAF/pack/proj 输出）。需要时直接阅读 `impg-0.4.1/src/syng.rs` 与 `docs/` 下 syng
相关文档。

## 6. 图构建层：GFA 管道与 crush

本层把"查询得到的同源片段集合"物化为 GFA 图，并可选地压缩 bubble。两种物化形态：直接构建
（`graph` 单图）与分区构建（`partition` 切窗口 + `lace` 拼回），是同一层的两种内存策略。`crush`
是图构建的后处理（压缩 bubble），不是独立的图构建路径。

### 6.1 引擎选择

`graph`、`query -o gfa`、`partition -o gfa` 三个命令共享同一套引擎实现，由 `--gfa-engine`
选择 ([lib.rs#L38](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/lib.rs#L38)的 `GfaEngine`
枚举)：

| Engine                   | Pipeline                                    | 用途                  |
|--------------------------|---------------------------------------------|-----------------------|
| `Pggb` (默认)            | sweepga + seqwish + smoothxg 平滑 + gfaffix | 平滑变异图            |
| `Seqwish`                | sweepga + seqwish + gfaffix                 | 原始（未平滑）图      |
| `Poa`                    | 单遍 SPOA                                   | 小区域、快速 MSA 输出 |
| `SyngNative`/`SyngLocal` | syng 锚点 + BiWFA + allwave + seqwish       | 近缘单倍型快速通道    |

### 6.2 stage 化管道

通过 `-o gfa:<stage1>:<stage2>:...` 简写，用户可在引擎前后插入 stage（见 §2.4）：

```
-o gfa:cut-n=100:pggb:crush,method=allwave:sort,pipeline=Ygs
       ^^^^^^^^^^ ^^^^ ^^^^^^^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^
       终端N裁剪  引擎  bubble 解析           最终排序
```

`build_graph_config` 与 `build_engine_opts` 把这些 stage 装配成 `EngineOpts` 结构，传给
`commands::graph::run` / `commands::partition::run` /`commands::syng2gfa::run`。

### 6.3 partitioned 模式与 `lace`：什么是"分区"

**"分区"分的是基因组区段（loci），不是切窗口、不是分序列、不是分样本。**

`Partition` 命令把整个 cohort 的基因组切成一组**互不重叠的 loci**，每个 locus 是
**跨多条序列 的同源区段集合**。目的是控内存：大 cohort（如 100+ 人类基因组）无法一次性构建 GFA，
切成 per-locus 批次，每个 locus 独立构建 GFA，最后用 `lace` 拼回。

**算法**
（[`partition_alignments`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/commands/partition.rs#L158)，
核心循环 L295-L570）：

1. 初始化 `missing_regions` = 所有序列的全长（尚未被任何 partition 覆盖）； `masked_regions` = 空
   （已被 partition 认领的区段）。
2. `select_and_window_sequences`（
   [L714](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/commands/partition.rs#L714)）按
   `--selection-mode`（longest/total/sample/haplotype）选"缺失最多"的序列，切成 `--window-size` bp
   的**初始窗口**（窗口只是起点，不是最终 partition 边界）。
3. 对每个窗口跑**传递闭包 BFS**（`query_transitive_bfs`，与 `query -x` 同一实现），找全所有同源片段。
4. `merge_overlaps`（合并间隔 ≤ `--merge-distance` 的区间）。
5. `mask_and_update_regions`：把这些区间从 `missing_regions` 减去，加入 `masked_regions` （后续
   partition 不会再认领这些区段，保证 loci 互不重叠）。
6. 输出这个 partition（BED/FASTA/GFA/VCF/MAF）。
7. 回到 2，选下一个"缺失最多"的序列。当 `missing_regions` 全空时停止（L906）。
8. `rehome_singleton_slivers`（
   [L27](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/commands/partition.rs#L27)）把贪心
   masking 产生的碎片归并到相邻 partition。

**关键认识**：

- Partition 内部用的就是 `query_transitive_bfs`——**Partition = 批量化的传递闭包查询 + masking 去重**。
- 同一个区段只属于一个 partition（`masked_regions` 保证），partitions 之间互不重叠。
- `--window-size` 是初始窗口大小，但 partition 的实际边界由传递闭包结果 + masking 决定，
  **不是简单按 WINDOW 切分**。
- `:WINDOW` 后缀（如 `pggb:10000`）触发 partitioned 模式：每个 partition 独立构建 GFA， 最后用
  `impg lace` 拼回一张图，可选 `--fill-gaps` 用参考序列填充窗口间空隙。

`lace` 同时支持 GFA 与 VCF 输入，路径名必须遵循 `NAME:START-END` 约定（最后一个 `:` 是分隔符），
坐标驱动重新拼装。

### 6.4 GFA 管道的阶段化组织

GFA 构建管道分散在 8 个源文件中，但按**管道阶段**而非文件来理解更清晰。一条完整的
`gfa:cut-n=100:pggb:crush:sort` 管道经过以下阶段：

1. **引擎选择与配置**（`lib.rs` 的 `GfaEngine` 枚举 + `EngineOpts`）— 5 种引擎 （`Pggb`/`Seqwish`/
   `Poa`/`SyngNative`/`SyngLocal`）决定基础管道形态。`SmoothConfig`定义 post-crush smoothxg pass
   的默认值（`target_poa_lengths=[700,1100]`、`max_node_length=100`，与 pggb 一致）。

2. **stage DSL 解析**（`graph_pipeline.rs` 的 `GraphPipelineSpec`）— 把 `stage,key=value:stage,...`
   字符串解析为类型化的 `Vec<GraphPipelineStage>`，只做语法校验，不决定可执行性。这是 §2.4 stage
   解析器的下游消费方。

3. **管道执行入口**（`commands/graph.rs`）— `build_graph`/`induce_graph_from_alignment`
   （seqwish 传递闭包）+ `run_graph_build_*` 系列函数（poa/pggb/partitioned 四种模式）。
   `GraphBuildConfig` 含 20+ 字段（threads/frequency/min_aln_length/repeat_max/min_match_len/
   adaptive_min_match_len/sparse_factor/transclose_batch/disk_backed 等），是管道行为的 完整参数化。

**避免 all-vs-all 比对的机制**（align 阶段）— impg 的 Graph/Query/Partition 命令在 align
   阶段默认跑 self-PAF（all-vs-all），但提供 4 层参数控制"挑选哪些对比对"，避免 O(n²) 全对比：

    - **`--sparsify <策略>`**（
      [`GraphBuildConfig.sparsify`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/commands/graph.rs#L98)，
      类型 `sweepga::knn_graph::SparsificationStrategy`）— wfmash-only 的 pair-selection
      heuristic。用 Mash sketch 算全对距离矩阵，建 KNN 图，只对 K 近邻边跑比对，把
      N² 降到 N×K。README 描述为 `--sparsify auto # pair-selection heuristic`。
      [main.rs#L3990](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L3990)
      验证逻辑："`--sparsify controls external-aligner pair selection`"。
      [graph.rs#L1153-L1243](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/commands/graph.rs#L1153)
      分支：`None`/`WfmashDensity` 走 wfmash 密度模式，其他变体走 `sweepga_align`（Mash KNN）。
      **这是"连通性达到一定程度就放弃 all-vs-all"的核心参数**。
    - **Syng 引擎**（`--gfa-engine syng`/`syng-local`）— 从 syng syncmer anchor counts 选对，
      [main.rs#L4009](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L4009)：
      "Sparsification is driven by syng anchor density and sweepga::knn_graph"，完全跳过外部比对器。
    - **`--num-mappings <m:n>`**（默认 many:many）— 限制 query:target 维度映射数（如 `1:1` 只留
      最佳映射），在 sweepga plane-sweep 过滤阶段裁剪。
    - **`--sparse-factor 0.0`**
      （[`SeqwishOpts`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L2076)，默认
      0.0=保留全部）— 在 seqwish 图归纳阶段丢弃一定比例的 input matches。另外**Partition 的 `--selection-mode`**（longest/total/sample/haplotype，
   [partition.rs#L714-L906](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/commands/partition.rs#L714)）
   是基于覆盖度的贪心选择：跟踪 `missing_regions`，选"缺失最多"的序列开窗，当所有序列都被
   覆盖时停止（L906）。但 Partition 消费已有 PAF，不跑比对——它控内存的方式是"按 locus 切分 cohort"
   （详见 §6.3），而非"挑选对比对"。
4. **图操作原语**（`graph.rs`）— `unchop_gfa`/`sort_gfa`（gfasort Ygs 集成）/ `build_spoa_engine`/
   `feed_sequences_to_graph`（SPOA）/`terminal_n_clip_span`/`prepare_poa_graph_and_sequences`。
   这些是各阶段共享的低层原语。
5. **crush 阶段**（`resolution.rs`，§6.5）— 由 `parse_crush_stage` 装配 `ResolutionConfig` 后调用
   `resolve_gfa_bubbles`。
6. **smooth 阶段**（`smooth.rs`）— smoothxg 风格块分解 + SPOA 平滑，`SmoothConfig` 的
   `block_source: SmoothBlockSource`（PathOverlap/Flubble/NeighborMergePoasta 三策略）决定块来源。
7. **后处理**（`gfa_self_loops.rs` + `commands/syng2gfa.rs`）— `NormalizeSelfLoops` 折叠 blunt GFA
   路径局部的 self-loop 重复单元；`syng2gfa` 处理 syng → GFA 物化 （`SyngGfaMode::Blunt` 物化精确
   source-spelling 0M paths，`Raw` 输出 syng-native overlap 图）。

**管道外的评估设施** — `graph_report.rs` 虽不参与构建，但为所有 sweep 脚本提供评分依据：
`GraphReport` 输出 15+ 维度（`status`/`failures`/`warnings` + 各维度测量值），序列化为 JSON 供
`c4-crush-cmaes.py` 等优化器消费（§10.2）。这是"把图质量变成可排序数值向量"的关键设施。

### 6.5 crush 算法：bubble 压缩的后处理

crush 是图构建层的后处理步骤——在 GFA 图物化后，压缩"bubble-like"结构
（同一对 boundary 之间的 多条可替换路径）为更紧凑的变异图。实体在
[`resolution.rs`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/resolution.rs)（17169 行，
全项目最大源文件）。

#### 6.5.1 crush 要解决什么问题

pggb 走的是"all-vs-all 比对 → seqwish 诱导图 → smoothxg 全图平滑"的路线。问题是 smoothxg 在
大图上做 all-vs-all 局部 MSA，**代价随图规模非线性增长**。对于 HLA、C4 这种高变 locus，PGGB 跑
一遍要数小时甚至跑不动。

crush 的核心动机记录在
[crush-architecture-spec.md](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/docs/crush-architecture-spec.md)
开头："three days of agents working on crush have been flailing without a clear spec"——
**agent 开发失去方向后才补的规范**。其设计目标明确：

- **避免 all-to-all alignment**，改用**per-bubble 局部 condensation**（这是 crush vs PGGB
  的全部意义）。
- **分钟级**完成单 locus（C4 GRCh38 规模），靠 rayon 并行 + level hierarchy 限制迭代深度。
- **路径保真不变量**：每条输入路径的 byte 序列在输出图中必须逐字节保留（emitted paths are the
  coordinate system）。不允许 lossy representative collapse。
- **PGGB 不是质量 oracle 的字节级目标**，只做"合理压缩"。

#### 6.5.2 算法的 8 阶段（来自 spec，对应代码实现）

1. **Bubble 检测** — POVU/flubble 检测（Lean-proven correct），每轮识别 top-level 非重叠 flubbles。
2. **按 bubble 规模分层路由 aligner** — 这是 crush 的核心设计决策：
    - median < ~1 kb → **sPOA**（短串联重复线性化为 clean linear representation）
    - ~1 kb – ~10 kb → **POASTA**（中等规模快速 POA，仍近线性输出）
    - ≥ ~10 kb → **sweepga**（基因转换/非等位同源重组预期非线性拓扑，sweepga 保非线性结构）
    - **aligner 由"输出允许长什么样"决定**，不是单纯比速度。
3. **局部图构建** — 每个 bubble 独立抽取内部路径 → 喂给 aligner → 得到局部压缩图 → 分配
   `[n+current_size, ...)` 区间的新 node ID（避免与存活 segment 冲突）。
4. **Strike + link** — 删除旧内部 segment，插入新局部图，连接 boundary 到新路径首尾。
5. **批处理 lacing** — 同一轮的非重叠 batch 整体并行处理（rayon），单批更新父图后排序。 非重叠由
   level hierarchy 保证。
6. **Level descent 迭代** — 解决 level-1 后内部可能还有更小 bubble，下降到 level-2 继续解决，
   不重复处理已解决的 level-1。
7. **POA 规整短重复 motif** — sweepga/allwave 输出的局部图中若有紧密高拷贝节点（短串联重复）， 再用
   sPOA/POASTA 局部重对齐线性化，为下游 genotyping 输出干净图。
8. **(可选) 最终 gfaffix 去重 pass** — 若前 7 阶段正确，此步应为 no-op。

**算法产生的不变量**（by construction）：

1. 不存在两个 segment 共享相同序列（若存在则说明 bubble 重叠/boundary 被误吸收/aligner 冗余）。
2. 每条输入路径逐字节保留。
3. 每轮图只缩不涨（segment 数与 bp 都 ≤ 上轮）；若涨说明算法在自我打架。
4. 新旧 node ID 不复用。

#### 6.5.3 15 种 `ResolutionMethod` 是演进史，不是设计

代码中 [`ResolutionMethod`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/resolution.rs#L274)
有 15 个变体。**这不是"15 种可选策略"的设计，而是 8 个月实验迭代留下的考古层**：

- 早期：`Poa`/`Poasta`/`Abpoa`/`StarBiwfa` — 直接 POA 替换的尝试，`StarBiwfa` 现仅作调试。
- 中期：`Allwave`/`Sweepga`/`Wfmash` — 引入 seqwish + SPOA polish 的全量对比对路线，其中 `Wfmash`
  专为与 PGGB 对比保留。
- 成熟期：`Auto`（默认，median 三档路由）+ `Hierarchical`（按深度路由）—— 这两个是 spec §Phase-2
  的代码化，**实际生产中只用这两个**。
- 后期补救：`ChainGreedy`/`ChainPovu`/`TopFlubbleSweepga`/`IterativeMultiLevel`/
  `CoverageMultiBubble`/`MotifLocal` — 针对 C4 难用例的 residual 提出的"更激进窗口策略"。
  `MultiLevelWindowMode`（7 变体）和 `MultiLevelObjectiveMode`（2 变体）只为这些方法服务。

**教训**：impg 自己的 C4 sweep 脚本（§10.2）显示大部分场景只跑 `Auto`/`Hierarchical`/`Sweepga` 3-4
种。**15 种 method 的存在本身就是工程失控的标志**——每个新方法都是"前一个搞不定 C4 时加的 新尝试"，
但没有删除旧方法的纪律。`pgr` 若做类似 bubble 处理，应坚持单一 `Auto` 路由（median 三档），不积累
method zoo。

#### 6.5.4 关键入口与配置

- **入口函数**
  [`resolve_gfa_bubbles`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/resolution.rs#L997)
  接收 GFA 1.0 字符串（仅 `S`/`L`/`P`，link 必须 blunt `0M`），返回 `ResolvedGfa { gfa, stats }`。
- **`ResolutionConfig`** ([L35](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/resolution.rs#L35))：
  `max_iterations`（默认 1，`until-done` = `usize::MAX`）、`method`、
  `auto_spoa_max_traversal_len`/`auto_poasta_max_traversal_len`（三档阈值，设 0 禁用对应档）。
- **关键常量** ([L505-L569](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/resolution.rs#L505))：
  `DEFAULT_AUTO_SPOA_MAX_TRAVERSAL_LEN=1_000`、`DEFAULT_AUTO_POASTA_MAX_TRAVERSAL_LEN=10_000`、
  `DEFAULT_REPLACEMENT_SEQWISH_MIN_MATCH_LEN=311` 等——这些数字直接对应 spec Phase-2 阈值，是 C4
  调参的产物。
- **诊断**：`IMPG_CRUSH_DEBUG_DIR` 环境变量触发 `DEBUG_REPLACEMENT_ID` 等原子计数器输出，被
  `audit_poasta_replacement_cycles.py`（§10.2）消费。

### 6.6 `examples/` 诊断与实验工具（7 个独立可执行程序）

`examples/` 下 7 个 `.rs` 文件，都是 `cargo run --example <name>` 程序，服务于 crush/smooth 算法的
**离线诊断与实验**，不属于正式 CLI。按用途分三类：

- **path-preserving 校验（2 个）** —
  [`compare_gfa_paths.rs`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/examples/compare_gfa_paths.rs)
  比较两个 GFA 路径拼写是否逐字节一致；
  [`validate_gfa_path_sources.rs`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/examples/validate_gfa_path_sources.rs)
  验证路径名 `NAME:START-END` 的拼写与源序列文件一致（含反向互补）。
  **这是任何"路径保真图操作" 都需要的端到端校验**，pgr 若引入 GFA 操作可直接复用这两个 example
  的思路。
- **POVU/POA 诊断（3 个）** — `povu_decomp_report.rs`（POVU flubble 分解 TSV）、
  `neighbor_merge_existing_gfa.rs`（独立调用 smooth 管道）、`poasta_order_driver.rs`（809 行，
  最大的 example，测试 8+ 种序列顺序策略对 POA 输出图质量的影响，对应 `docs/poasta-order-*.md`）。
  后者的"多策略 × 多指标 TSV 输出"实验框架值得借鉴——把"POA 顺序敏感性"这个模糊问题变成可排序的指标表。
- **syng 探测（2 个，本文档不参考）** — `syng_anchor_probe.rs`/`syng_probe.rs`。

### 6.7 `tests/` 测试体系（13 个集成测试 + 2 个验证脚本 + fixture 矩阵）

`impg` 没有传统的 `src/` 内单元测试，所有测试都在 `tests/` 目录下，共 13 个 `test_*.rs`
集成测试 文件 （合计 10866 行）+ 2 个 shell 验证脚本 + 一个中心化 fixture 矩阵。这与 `pgr` 的
`tests/cli_*.rs`约定形成鲜明对比（见本节末尾）。

**测试代码（`tests/*.rs`，13 文件）**按主题分四类：

- **crush 回归** — `test_crush_integration.rs`（1674 行，最大非 syng 测试，C4A 切片 2942 segments/
  465 haplotypes）+ `test_local_compression_testbed.rs`（645 行，消费 fixture 矩阵，13 classes ×12
  methods 笛卡尔积）+ `test_graph_output_crush.rs`（端到端 `gfa:pggb:crush`）。
- **传递闭包与投影** — `test_transitive_integrity.rs`（767 行，传递闭包完整性：非重叠区域保持分离、
  坐标投影准确、双向查询对称）+ `test_gfa_projection.rs`（GAF→GFA projection）。
- **图构建引擎** — `test_graph_seqwish.rs`/`test_graph_poa.rs`/`test_pipeline_integration.rs`。
- **syng/genotype**（本文档不参考）— `test_syng_integration.rs`（5302 行，最大）、
  `test_genotype_validation_suite.rs`、`test_genotype_gfa.rs` 等。

**测试数据（`tests/test_data/`，三层 fixture）**

- **根目录小 fixture** — `a.fa`/`b.fa`/`c.fa` 等 30-80 字节级，用于基础索引查询。
- **`crush/` 子目录（6.9MB）** — crush 算法真实数据：`c4_slice_1500_3000.gfa`（7.2MB C4A 切片）+
  `c4_fragments/`（4 个 C4 子片段：`easy_shared_flank`/`bounded_multi_bubble`/`short_floor`/
  `duplicated_repeat`）。
- **`local_compression/` 子目录（216KB）—— 最有组织的测试数据集**。由
  `scripts/local_compression_testbed.py write-fixtures` 生成，`manifest.json` 中心化管理。
  13 个 fixture class 覆盖所有局部压缩场景（snp_bubble / short_indel / alu_like_insertion /
  tandem_copy_number_loop_cyclic / inversion_like / nested_bubbles 等）。每个 fixture 子目录含
  `input.fa` + `expected_paths.tsv`（路径保真基准）+ `metadata.json`（含 `expected_topology`断言 +
  `allowed_ranges` 拓扑区间 + `known_failure_mode`）+ `notes.md`。tier 标签区分 `ci`（快速 CI 子集）
  与 `local`（完整测试）。

**验证脚本（`tests/validation/`，2 个 bash）** — `battery_syng_vs_paf.sh`（syng vs PAF 批量对比，
11 列 TSV 诊断）+ `compare_syng_vs_paf.sh`（单次对比）。均为诊断工具，退出码恒 0。

**测试基础设施约定**

- 二进制定位：`CARGO_BIN_EXE_impg` 环境变量，回退到 `target/{debug,release}/impg`，再回退到
  `/home/erikg/impg/target/release/impg`（作者机器路径，硬编码）。
- 外部依赖：gfaffix 必须存在；wfmash/samtools 可选（`#[ignore]` 标记）。
- 串行化：syng 测试用 `LazyLock<Mutex<()>>` 全局锁（C 库非线程安全）；C4 fragment 测试用 `Once` +
  `Mutex` 32 线程池串行。
- **已知失败用 `#[ignore]` 标记保持 CI 绿**，文档注释说明复现命令与对应 docs/ 审计条目——这是 "agent
  开发导致 RED 测试堆积"的妥协产物。

**与 `pgr` 测试约定的对比**：`pgr` 用 `tests/cli_<command>.rs` + `PgrCmd` 辅助结构体
（`tests/common/mod.rs`），测试数据放 `tests/<command>/`，强调 Zero Panic。impg 则把测试数据
中心化到 `manifest.json` + tier 分级 + `allowed_ranges` 区间断言——这套做法对"多方法 × 多场景"
的回归矩阵特别有效，`pgr` 若引入 crush/POA 类算法可借鉴。

## 7. 应用层：基因分型与等位基因推断

应用层是能力栈的顶端消费者，依赖图构建层（§6）或隐式图查询（§4）提供 evidence。
impg 的 `genotype`/`infer` 命令实现基于 graph feature 的基因分型，核心设计在
[`docs/genotype-architecture.md`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/docs/genotype-architecture.md)
与 [`docs/infer-design.md`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/docs/infer-design.md)
（README 显式引用）。

### 7.1 基因分型架构（5 步模型）

`genotype-architecture.md` 定义了**后端中立的 5 步模型**——围绕 graph feature evidence 而非
特定图表示：

1. **选 locus** — 确定要分型的基因组区域。
2. **提取 candidate haplotype** — 从图或隐式图查询得到该 locus 的候选单倍型。
3. **表示为 graph feature 向量** — 把候选 haplotype 编码为特征向量（节点/边/路径 presence）。
4. **sample 表示为 coverage/support 向量** — 把样本 reads 投影到特征向量上，得到支持度。
5. **score ploidy-sized 组合** — 在 ploidy 约束下搜索最优等位基因组合。

**关键设计**：graph feature 是后端中立的——既可来自物化 GFA 图，也可来自隐式图查询。这意味着
基因分型不强制要求图构建层，可直接在查询层（§4）的输出上工作。

### 7.2 infer 设计：跨区间 stitching

`infer-design.md` 定义跨区间/分区输出等位基因 call 的 stitching/mosaic 设计。
`pangenome- genotyping-roadmap.md` 给出长期目标管线：

```
panel sequences/pangenome graph → implicit graph backend → sample evidence projection
→ local candidate subwalks → local genotype scoring → recombination/copying inference
→ inferred phased haplotype mosaics
```

### 7.3 对 pgr 的意义

应用层是 pgr 当前完全空白的能力。但 §9 启示将说明：pgr 现阶段不应直接实现完整基因分型管道，
而应先打通"索引 → 查询 → 图构建"三层基础设施。基因分型的 5 步模型可作为未来设计的目标架构参考，
尤其是"后端中立的 graph feature"思路——与 pgr 的隐式图路线契合。

## 8. 对比分析：impg vs pgr

`pgr` 的真正强项是**UCSC 体系的 pairwise 比对处理**（Chain/Net/MAF/AXT/PSL/LAV 全套， 见
[docs/chain.md](file:///Volumes/ExtHome/Scripts/pgr/docs/chain.md)）与**Block FA 多序列比对**
（`fas` 全套子命令 + `libs/poa/` 的 SPOA 移植 + `libs/fas_multiz.rs` 的 multiz 风格 banded DP 合并，
其 `FasMultizMode::Core` 即"多基因组共享 core 比对"）。pairwise 与 core 比对均已成熟。

**泛基因组图部分是 `pgr` 的空白**：[docs/gfa.md](file:///Volumes/ExtHome/Scripts/pgr/docs/gfa.md)
明确写"如果 `pgr` 未来涉及泛基因组操作"，是规划/知识背景文档而非实现；`src/cmd_pgr/` 下无 `gfa`
子命令，`src/libs/` 下无 GFA 模块。本节聚焦"泛基因组图"这一维度对比，作为 §9 启示的依据。

- **pairwise 比对** — `pgr` 成熟（AXT/MAF/PSL/Chain/Net/LAV）；`impg` 成熟（PAF/1ALN/TPA）
- **core 比对** — `pgr` 已实现（`fas multiz --mode core` + `fas consensus` POA）；`impg` cohort
  all-vs-all + 投影间接得到
- **图模型** — `pgr` 显式 Chain/Net（GFA 尚为规划）；`impg` 隐式图（比对网络），按需物化 GFA
- **核心数据结构** — `pgr` Newick 树 / PSL / Chain / Block FA；`impg` coitrees + 紧凑 CIGAR delta
- **比对输入** — `pgr` AXT/MAF/PSL/Chain（UCSC 风格）；`impg` PAF/1ALN/TPA（wfmash 风格）
- **查询模式** — `pgr` 按 coordinate 直接读取 / `psl lift` 线性单链投影；`impg` 区间投影 + 传递闭包
- **泛基因组图构建** — `pgr` 无（gfa.md 规划中）；`impg` 内嵌 sweepga/seqwish/allwave/crush 管道
- **bounded bubble 处理** — `pgr` 无；`impg` crush 算法（POVU + aligner 路由 + polish）
- **免比对后端** — `pgr` 无；`impg` syng syncmer GBWT
- **基因分型** — `pgr` 无；`impg` `genotype cos` + `infer` (cosigt 模型)
- **CLI 风格** — `pgr` `pgr <format> <subcommand>` 多级；`impg` 单级，参数众多
- **代码组织** — `pgr` `cmd_pgr/` 下按格式分组；`impg` `commands/` 分文件，但 `main.rs` 巨大

### 8.1 impg 的局限性（启示的前提）

impg 的设计有很多亮点，但以下局限性是 `pgr` 借鉴时必须清醒认识的——这些局限性正是 §9 启示的 前提。

1. **all-vs-all O(n²) 不可避免** — impg 把"图构建"延迟到查询时，但**没有摆脱 O(n²) all-vs-all 比对**。
   隐式图的真正优势不在"省比对"，而在"按需计算图遍历"（§1.1.2）。pgr 若走隐式图路线，应复用已有的
   UCSC Chain（已经过 syntenic 净化），而非重新跑 all-vs-all。
2. **传递闭包 ≠ 多序列比对** — 传递闭包只做"图遍历可达性查询"，不产出多重比对；间接同源 （A~B、B~C）
   会被报告，但三者的 MSA 需另走 `gfa`/`maf` 输出路径（机制详见 §4.3）。pgr 借鉴时 需清醒认识"图遍历
   vs 多序列比对"的本质区别。
3. **15 method zoo 是工程失控的标志** — crush 的 15 种 `ResolutionMethod` 不是"15 种可选策略"的设计，
   而是 8 个月实验迭代留下的考古层。每个新方法都是"前一个搞不定 C4 时加的新尝试"，但没有删除旧
   方法的纪律（§6.5.3）。pgr 若做类似 bubble 处理，应坚持单一 `Auto` 路由，不积累 method zoo。
4. **C4 难用例至今未完全解决** — `validate-current-c4.md` 的结论 "better than prior SYNG-derived
   C4 outputs on the specific repeat-artifact failure mode, but still not solved and not yet..."
   （§10.3.5）反映 crush 算法在 C4 这一"超难"区域至今未完全解决。这是 pgr 决定"现阶段不需要 crush"
   （§9.3）的实证依据。
5. **agent 评分系统与实际完成度脱节** — `evaluations/` 下大量 `Overall score: 0.00/1.00` 与
   `Rubric underspecified: true`（§10.3.5）反映 agent 评分系统与任务实际完成度的脱节。pgr
   若引入类似 agent 工作流，应确保评分 rubric 明确。
6. **main.rs 巨型化** — 单文件 1.5 万行 + 60 万字符，承担所有 clap 命令定义、参数解析、stage 解析、
   命令分发（§2.1）。这是 impg 最显著的工程反例，pgr 应坚持 `cmd_pgr/` 按格式/功能分组的结构。
7. **single-developer 工作流痕迹** — artifact 路径硬编码 `/home/erikg/impg/data/<RUN_ID>/`、PNG
   上传 `hypervolu.me/~erik/impg/`、测试二进制定位回退到 `/home/erikg/impg/target/release/impg`
   （§6.7）。这些是 single-developer (Erik Garrison) 工作流的产物，pgr 移植时需参数化。
8. **大量 `unwrap_or_else(|e| panic!(...))`** — impg 源码中存在大量 panic 风格错误处理（如
   `get_cigar_ops`、`get_target_sequence_cached`），违反 pgr 的 Zero Panic 原则。pgr 借鉴其算法时
   应改为 `anyhow::Result` + `bail!`。

## 9. 对 `pgr` 的启示（聚焦泛基因组部分）

`pgr` 的 pairwise 比对与多基因组 core 比对已成熟，无需从 impg 借鉴；且 `libs/poa/` 已有完整的
SPOA 移植（`Poa` struct，`add_sequence`/`consensus`/`msa`），被 `fas consensus` 与 `fas refine`
使用——这恰好是 impg crush 算法的 aligner 基础设施。下列启示聚焦**泛基因组部分**，围绕四个关键
问题展开，再附工程层借鉴要点。

### 9.1 路线选择：隐式图 vs 物化图

**答：pgr 应走"Chain/Net 即隐式图"的路线，而非重新构建 PAF-based 隐式图，也非先物化 GFA。**

impg 的隐式图以 PAF/1ALN/TPA all-vs-all 比对为边集；pgr 的强项是 UCSC Chain/Net 体系，
**Chain 本身就是隐式图的边**——一条 Chain = 一条 A↔B 边集（含 indel gap）。pgr 无需引入 PAF/wfmash，
可直接 在 Chain/Net 上构建区间索引，复用已有的 pairwise 成熟基础设施（见 §9.2）。

物化图（pggb/Minigraph-Cactus 路线）的代价是"即使用户只关心一个位点，也要先构建整张图"。pgr 当前
没有 GFA 构建管道（[docs/gfa.md](file:///Volumes/ExtHome/Scripts/pgr/docs/gfa.md) 是规划文档），
**"先物化再分析"对 pgr 是过载的**——应优先实现按需投影，把物化推迟到用户显式要求 `gfa`/`maf`输出时。

§1.1.2 的适用边界表给出判断依据：pgr 的典型场景（UCSC 风格的 locus 查询、cohort 区间投影）落在
"隐式图赢"的区域；只有未来需要全图统计（coverage/layout/variant calling）时才考虑物化。

### 9.2 Chain/Net 体系如何与隐式图结合

**答：把 Chain 索引到 coitrees，复用 impg 的传递闭包 BFS，实现 cohort 级"找全所有同源片段"。**

具体路径：

1. **Chain → 区间树** — 一条 Chain 的每个 block（匹配段）是一条 A[start-end]↔B[start-end] 边。
   把所有 cohort 的 Chain 文件解析后，按 target seq 建立区间树（可直接复用 `coitrees` crate），
   节点 metadata 用 `QueryMetadata` 风格的 bit-packing（§3.2）：query_id/target_start/target_end/
   query_start/query_end/strand/data_offset。Chain 的 CIGAR 信息（match/insertion/deletion run）
   可懒加载，只存文件偏移量。
2. **`chain lift` 升级为传递闭包** — 当前 `pgr chain lift` 是单条 Chain 的线性投影。可在区间树
   之上实现 `chain lift --transitive`（对应 impg `-x`）：从目标区间出发做 BFS，自动发现所有直接
   和间接同源片段，包括通过第三序列的间接同源。`-m` 控制深度，`-d` 控制 hop 内最大 gap。
3. **Net 作为 syntenic 过滤器** — Net 编码了 syntenic 关系（top-level vs alignmentNet vs classNet），
   可作为传递闭包 BFS 的剪枝依据：只在 syntenic Net 内扩展，避免 paralog 横跳引入 噪声。这是 pgr
   相对 impg 的潜在优势——impg 的 PAF 没有 syntenic 注释，只能靠 `-d` 硬阈值。
4. **输出格式** — 复用 impg 的输出矩阵（bed/bedpe/paf/gfa/vcf/maf/fasta），但 pgr 的 MAF 输出
   可直接复用 `libs/fas_multiz.rs` 的 banded DP 合并（`FasMultizMode::Core`），比 impg 的
   per-bubble POA 更适合 core 区段。

**前提设施**：当前 `pgr` 的 Chain 处理是流式的，缺少随机访问能力。第一步需引入 Chain 索引格式
（类似 `.impg` 但基于 Chain），把全 cohort 的 Chain 装进内存区间树。

### 9.3 crush 算法的必要性

**答：现阶段不需要，未来若需要可只移植 median 三档路由核心。**

理由（详见 §6.5）：

- crush 解决的是"已有 GFA 图后如何压缩 bubble"，是图构建**之后**的优化步骤。pgr 当前还没有 GFA
  构建管道，先解决"能否构建图"再考虑"如何压缩图"。
- pgr 的 `fas consensus`（`libs/poa/` SPOA）已能做单 locus 的 MSA consensus，对于"只想看一个位点"
  的场景，隐式图投影 + 局部 MSA 比 crush 更轻量。
- crush 的 8 阶段流程 + 15 method 复杂度对 pgr 是过载的。impg 自己的 C4 sweep（§10.2/§10.3.2）显示
  大部分场景只跑 `Auto`/`Hierarchical`/`Sweepga` 3-4 种，15 种 method 的存在本身就是工程失控的
  标志。若未来确需，可只移植 §6.5.2 的 median 三档路由核心（< 1kb → SPOA、1-10kb → POASTA、≥ 10kb
  → sweepga），跳过 8 个月的 method zoo 考古层。
- `libs/poa/` 已是 crush 算法 aligner 层的现成基础。`fas consensus` 已验证该 POA 在 MSA 场景可用，
  迁移到 graph bubble 场景的门槛低于从零起步。

### 9.4 第一步最小原型

**答：`pgr chain query` —— Chain 索引 + 区间投影 + 传递闭包，输出 BED/MAF。**

具体目标：

1. **Chain 索引格式**（`.pgr.chain.idx`）— 把多个 Chain 文件解析后按 target seq 建立区间树，
   序列化为单文件索引。参考 impg 的 `Index` 命令与 `.impg` 格式，但 metadata 用 Chain block 而非
   PAF record。**索引层不过滤**（遵循 impg §1.1.3 的设计：全量装入，挑选推迟到查询层）。
2. **`pgr chain query <region>`** — 在区间树上查找所有重叠 Chain block，把目标坐标 lift
   到查询 序列，输出 BED（最简输出）。查询参数暴露过滤：`--merge-distance`/`--min-identity`/
   `--min-output-length`/`--subset-sequence-list`（对应 impg `QueryOpts` 的 `-d`/
   `--min-result-identity`/`-l`/`--subset-sequence-list`）。
3. **`pgr chain query <region> --transitive`** — 实现传递闭包 BFS（对应 impg `-x`），找全所有
   同源片段。`--max-depth`/`--max-gap` 控制遍历（对应 impg `-m`/`-d`）。
4. **`pgr chain query <region> --transitive -o maf`** — 传递闭包结果 + 局部 MSA，复用
   `libs/fas_multiz.rs` 输出 MAF。这是"隐式图 + 按需物化"的最小闭环。

**验证标准**：

- 在已知 Chain 文件上，`chain query --transitive` 的结果应包含所有手动挑选 Chain 的 `chain lift`
  结果，且额外包含间接同源片段。
- MAF 输出的序列与源基因组 FASTA 逐字节一致（path 保真不变量）。
- Zero Panic：畸形/二进制 Chain 输入返回友好错误而非 panic。

**暂不实现**：GFA 物化、crush、partitioned 模式、partition/lace/refine 子命令。这些是"按需物化"
路径上的后续步骤，待最小原型验证后再决定是否需要。

**为何 pgr 不需要 `--sparsify`**：impg 的 `--sparsify auto`（§6.4）是为了在**没有现成比对**时 避免
all-vs-all 比对——用 Mash KNN 选对再跑 wfmash。pgr 的场景是**已有 UCSC Chain**（已经过 syntenic
净化），天然避开了 all-vs-all 问题。pgr 的"挑选"发生在查询层（`chain query` 的 `--merge-distance`
等参数），而非比对层。这是 pgr 相对 impg 的根本优势：复用成熟的 pairwise 基础设施，不重新跑
all-vs-all。

### 9.5 工程层借鉴要点

除上述路线决策外，impg 的以下工程实践值得 pgr 借鉴：

1. **区间树 + 紧凑 CIGAR delta** — `coitrees` + `CigarOp` 风格的 bit-packing（§3.2）把全基因组
   all-vs-all 比对装进内存。当前 `pgr` 的 PSL/Chain 处理是流式的，缺少随机访问能力——这是把
   pairwise/core 比对升级为"全 cohort 可查询隐式图"的前提设施。
2. **trait 抽象单/多文件索引** — `ImpgIndex` trait + `MultiImpg`（§3.4）是处理"单大文件 vs
   多小文件"两种部署模式的干净做法。`pgr` 若引入类似的索引层，可让命令代码与索引物理形态解耦。
3. **PAF `cg:Z:` 懒加载** — impg 不把 CIGAR 存进区间树节点，只存 `data_offset` + `data_bytes`，
   查询时按需读取。这对 `pgr` 处理大型 Chain/PSL 是直接可借鉴的内存优化。
4. **thread-local 缓存模式** — impg 用 `thread_local!` 缓存 WFA aligner、1aln/TPA 句柄、目标序列
   片段，避免重复分配。`pgr` 在并行处理 PSL/Chain 时可借鉴同样的模式。
5. **PanSN 命名约定** — impg 全程用 `sample#haplotype#contig` 命名（`#` 分隔），`pgr` 在 `pgr pl`
   流水线中若需要处理群体数据，可采用同一约定以与 pggb/impg/odgi 生态兼容。
6. **stage 化字符串 DSL 的反面教训** — impg 的 `-o gfa:cut-n=100:pggb:crush:sort` 简写表达力强，
   但代价是 `main.rs` 膨胀到 60 万字符。`pgr` 当前坚持 `pgr <format> <subcommand>` 的多级结构
   更易维护，应保持。若未来需要类似的管道组合，可考虑专门的 pipeline 配置文件而非 CLI 简写。
7. **避免 main.rs 巨型化** — impg 把 20 个子命令的 clap 定义与分发全塞进单文件是明显的反例。
   `pgr`的 `src/cmd_pgr/` 按格式/功能分组、每命令独立模块的结构更优，应继续坚持——`main.rs` 只做
   `ArgMatches` 分发，业务逻辑下沉到模块。
8. **Zero Panic 与 AGENTS.md 的契合** — impg 源码中存在大量 `unwrap_or_else(|e| panic!(...))`
   （如 `get_cigar_ops`、`get_target_sequence_cached`），违反了 `pgr` 的 Zero Panic 原则。`pgr`
   在借鉴其算法时应改为 `anyhow::Result` + `bail!`，把错误返回到调用方而非 panic。
9. **POA 基础设施可直接复用** — `pgr` 的 `libs/poa/`（`Poa` struct + `AlignmentParams` +
   `AlignmentType::Global`）已是 crush 算法 aligner 层的现成基础。`fas consensus` 已验证该 POA 在
   MSA 场景可用，迁移到 graph bubble 场景的门槛低于从零起步。

## 10. impg-0.4.1 的 notes/scripts/docs 目录

本章汇总 impg 仓库中三个辅助目录（notes/scripts/docs）的组织与内容，作为前述技术分析的补充 资料。
这些目录反映 impg 的开发工作流与文档实践，对 pgr 的工程层借鉴有参考价值。

### 10.1 `impg-0.4.1/notes/` 目录

`notes/` 是独立于 `docs/` 的小目录，仅 5 个 Markdown 文件：

| 文件                                      | 主题                             | 相关性                  |
|-------------------------------------------|----------------------------------|-------------------------|
| `FAST_MODE_IMPLEMENTATION.md`             | `.1aln` tracepoint + approx      | **相关** — 与 §3.2 同源 |
| `SYNG_NEXT_STEPS.md` (14 KB)              | syng 集成路线图（PR #162）       | 不参考                  |
| `SYNG_OPTION3_PAIRWISE_REFINE.md` (14 KB) | syng discovery + pairwise refine | 不参考                  |
| `SYNG_TRANSITIVE_DESIGN.md` (8 KB)        | syng-seeded transitive 投影      | 不参考                  |
| `SYNG_VS_PAF_VALIDATION.md` (7 KB)        | syng vs PAF 验证（yeast235）     | 不参考                  |

5 个文件中 4 个是 syng 专属笔记，与本文档范围无关；唯一值得关注的是 `FAST_MODE_IMPLEMENTATION.md`，
它是 alignment 后端（.1aln/TPA tracepoint）的性能优化指南，与 §3.2 的 CIGAR 懒加载思路一脉相承，
可视为该机制的进阶实现说明。

### 10.2 `impg-0.4.1/scripts/` 目录

`scripts/` 是 21 个独立可执行脚本（Python 14 + Bash 6 + R 1），几乎都是
**实验驱动器 (experiment driver)**：不属于 impg 二进制，而是封装"调用 impg 子命令 → 跑 sweep/参数
搜索 → 收集时间/RSS/graph-report → 写 TSV/JSON → 渲染"的工作流。重型 artifact 写到仓库外的
`/home/erikg/impg/data/<RUN_ID>/`，仅驱动器本身入库。

**主题分布**：C4 难用例实验矩阵占 10 个（围绕 GRCh38 chr6:31891045-32123783 这一"超难"区域跑
参数/算法矩阵），crush 诊断与通用图 QC 3 个，demo 等价验证 5 个，杂项 3 个。C4 主题占比近半，反映
crush/C4 难用例是 impg 开发后期的主战场。

**核心模式（对 pgr 的直接借鉴价值）**：

1. **driver 外置 + 评分解耦** — 所有 sweep 脚本把 impg 当黑盒 subprocess 调用，自身只做参数装配、
   artifact 落盘、TSV 汇总。`c4-crush-cmaes.py` 明言 "optimizer is deliberately external to
   impg"——候选生成与评分逻辑不入主二进制。`pgr` 设计泛基因组图构建的 benchmark 体系时可直接采用
   这一思路：主二进制只暴露 `--json` 输出，外部脚本驱动参数搜索。
2. **path 保真检查是事实标准** — `compare_gfa_paths`（§6.6）被几乎所有 sweep 脚本调用，做 path 序列
   逐字节保真度检查。任何"路径保真图操作"都需要端到端校验工具，pgr 若引入 GFA 操作可直接复用该
   example 的思路。
3. **graph-report 作为评分依据** — `GraphReport`（§6.4）的 15+ 维度输出 JSON，被 `c4-crush-cmaes.py`
   等优化器消费。把"图质量"变成可排序的数值向量，是从"模糊感觉图好不好"到 "可优化目标"的关键一步。
4. **`local_compression_testbed.py`（89 KB，全目录最大）** — 消费 fixture 矩阵（§6.7），13 classes
   × 12 methods 笛卡尔积，写 scoreboard。这是"多方法 × 多场景"回归矩阵的可执行实现，pgr 若引入
   crush/POA 类算法可借鉴该 fixture matrix 模式。

**单开发者工作流痕迹**：artifact 路径硬编码 `/home/erikg/impg/data/<RUN_ID>/`、PNG 上传
`hypervolu.me/~erik/impg/`、`hprcv2-syng-smoke.py` 等 syng 主题脚本（2 个，本笔记不参考）。pgr 移植
时需参数化路径。

### 10.3 impg-0.4.1/docs 文档结构

`impg-0.4.1/docs/` 是项目的开发文档目录，规模庞大：
**109 个顶层 Markdown 文件 + 2 个子目录 (`designs/`、`evaluations/`)，合计约 3.5 MB / 34619 行**。
它不是面向最终用户的文档（那是 `README.md` 与 `--help`），而是开发过程中的设计笔记、实验报告、bug
诊断与审计记录。下面按主题脉络梳理，突出对 `pgr` 有借鉴价值的模式与教训，而非逐文件罗列。

#### 10.3.1 crush 算法的"设计→实验→诊断"三循环（~64 个文档）

crush 是 docs/中最大的主题，约 64 个文档（21 设计 + 24 实验 + 19 诊断）。这不是"15 种可选策略"
的设计，而是 8 个月迭代留下的考古层——**三类文档形成循环**：写设计→跑实验→发现 bug 写诊断→
改设计→再跑实验。算法成熟度与文档数量成正比，但也反映工程失控风险。

**设计/规范类（~21 个）**的核心是
[`crush-architecture-spec.md`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/docs/crush-architecture-spec.md)：
开篇明言 "three days of agents working on crush have been flailing without a clear spec... No
code change to crush should be proposed without referencing this document"，定义 3 阶段算法
（POVU flubble 检测 → 按中位遍历长度路由 aligner → polish 直到收敛）。其余设计文档记录算法演进
的各个分支：层次化处理（`crush-hierarchical`/`crush-level-descent`/`crush-true-level-descent`）、
wider context（`crush-wider-context-bubbles`）、flank-aware（`flank-aware-crush-design`）、
neighbor-merge（`crush-neighbor-merge-iterate`）、local-compression-testbed
（`local-graph-compression-testbed-design`）等。**演进规律**：每个新设计文档都是"前一个搞不定 C4
时加的新尝试"，但没有删除旧设计的纪律（与 §6.5.3 的 15 method zoo 同源）。

**实验报告类（~24 个）**每个遵循"假设 → 实验设置（binary commit/输出目录/PNG
链接）→ TL;DR →详细结果"格式，多数以 C4 为 canonical 测试区域。**关键发现**：
[`crush-aligner-speed-study.md`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/docs/crush-aligner-speed-study.md)
记录 sPOA 占 99.2% 时间（831s）、POASTA 在相同输入上 9.93s（**84× 加速**），是 sPOA 被
POASTA 全面 取代的起点；`crush-exp-poasta-everywhere` 确认 POASTA 用于所有 bubble 是新最佳；
`crush-vs-pggb-comparison` 显示 crush 当前最佳 36:53/101GB 产出 19836 seg / 553585 bp，PGGB
13:38/64GB 产出 13288 seg / 234524 bp（**2.36× 更多序列**），crush 仍未追平 PGGB 的压缩率。
实验也记录**假设被反驳**的情况（`crush-exp-sweepga-k31`：降低 seqwish-k 不能 rescue 小 bubble；
`crush-exp-min-run-5`：min-run=5 与 min-run=3 mask 掉完全相同节点）——这种"明确记录失败"的纪律
值得借鉴。

**诊断/审计类（~19 个）**是回溯性文档，记录 bug 调查、性能瓶颈、行为审计。关键工具是
`IMPG_CRUSH_DEBUG_DIR` 环境变量（`crush-aligner-failure-trace` 抓 per-replacement 子图，无源码
改动）。`crush-poasta-pass-through-audit` 揭示一个重要设计决策：replacement 接受仅由 path-validity
gating，**无**图质量/compression-ratio/objective 检查守卫——这是"只保不变量、不保质量"的激进选择。
`crush-verify-report` 给出端到端验证：5 个真实 GFA × 3 scale，path 序列保留 ✓，whitespace 减少
18-66%。`crush-crush-handoff`（即 `c4-crush-handoff.md`）记录 branch/PR/binary commit、known-good
artifact、unrelated dirty state，便于下一轮从有用 artifact 继续——这种"工作交接"文档是 agent-assisted
开发的产物。

#### 10.3.2 C4 (HLA) 难用例作为主战场（~17 个文档 + evaluations/ 中 ~20 个）

C4 是 HLA 复合体区域（GRCh38#0#chr6:31891045-32123783，~232kb，HPRCv2 465 paths），impg 团队用作
"难测试用例"——高重复、高 CNV、SV 富集。docs/中 C4 主题文档约 17 个，evaluations/ 中还有 ~20 个
评估，合计占 docs/总量三分之一强，反映 crush/C4 难用例是开发后期的主战场。

**blocker 系列（5 个，编号 01-05b）** — 阻塞 C4 正确处理的疑难问题排查，从 POASTA scale、残留 路由、
全量重跑记分板到完整遍历聚合。`c4-blocker-05b` 记录关键决策：候选只在 cluster 范围 union 覆盖所有
graph path 时发出，bounded homologous node expansion 仅在 range-only union 不足时使用。

**CMA-ES 黑盒优化器（3 篇）** — 用 Python `cma` 包优化 crush 参数。`c4-cmaes-optimizer`
明言 "optimizer is deliberately external to impg"——每个候选作为普通 `impg crush` 命令运行，
`impg graph-report` 后由 wrapper 单独判分。`c4-cmaes-target-shape-results` 把 PGGB target spectrum
加入 objective（bp-weighted node copy-frequency distribution、EMD + TV、node length distribution
distance、excess total segment bp），是从"模糊质量"到"可优化数值向量"的关键一步。

**诊断与修复循环** — `diagnose-and-fix-c4-compound`/`diagnose-residual-two-c4`/
`diagnose-residual-underaligned-c4`/`fix-c4-syng-20260530`/`fix-top-flubble-sweepga` 记录具体 bug
与修复。`fix-top-flubble-sweepga` 揭示两个管道失败模式：min_match_len 大于所有 exact CIGAR run；零
PAF block 仍被接受为有效 replacement。

**集成与扩展** — `integrate-c4-local`/`iterative-multi-level-c4`/`expand-multi-bubble-c4`/
`settle-local-replacement`/`evaluate-low-min-match-c4` 记录算法扩展尝试。

#### 10.3.3 syng 后端文档（~8 个，含 `designs/`）— 略

包括 `designs/syng-integration.md`（集成架构总览）、`syng-gfa-query.md`（local GFA 查询配方，
README 显式引用）、`syng-parallel-construction.md`（并行构建 + 6 个 sidecar）等。需要时直接阅读
`impg-0.4.1/docs/` 下这些文件。

#### 10.3.4 其他主题文档（基因分型 / 图管道 / 外部工具审计，~15 个）

**基因分型与 infer（~7 个）**的核心是 `genotype-architecture.md`：围绕
graph-feature evidence 而非特定图表示，5 步模型（选 locus → 提取 candidate
haplotype → 表示为 graph feature 向量 →sample 表示为 coverage/support 向量 →
score ploidy-sized 组合）。`pangenome-genotyping-roadmap.md`给出长期目标管线：
`panel sequences/pangenome graph → implicit graph backend → sample evidence projection → local candidate subwalks → local genotype scoring → recombination/copying inference → inferred phased haplotype mosaics`。
`infer-design.md`（README 显式引用）定义跨区间/分区输出等位基因 call 的 stitching/mosaic 设计。

**图管道、DSL 与渲染（~3 个）** — `graph-pipeline-dsl.md` 定义 `-o gfa:<stage>:<stage>`
DSL（对应 §2.4 的 stage 解析）。`render-gbz-translation-design.md` 提出一个重要的抽象：
IMPG 是隐式泛基因组 上的翻译系统，scalable object 不是单张物化图，而是在源序列坐标、graph
feature、evidence projection、inferred haplotype 之间不丢身份地移动的能力；Root namespace 是
`source_sequence_id : [0, source_length)`。

**外部工具封装审计（4 个）** — 同行评审系列，审计 aligner 封装（wfmash/FastGA/SweepGA）与 render
封装（gfalook/odgi）的质量。

#### 10.3.5 `evaluations/` 子目录 — agent 工作流产物（~35 个文件 + 2 个嵌套子目录）

测试运行结果与评估脚本输出。**重要特点**：多数文件以 `Task: <id>` + `Evaluator: agent-<N>` +
`Date:` 开头，是 agent 系统自动评分的产物，常见 `Overall score: 0.00 / 1.00` 与
`Rubric underspecified: true/false` 字段。这反映了 impg 团队的 agent-assisted 开发模式：agent
执行任务 → 自动评分 → 滚动迭代。

**C4 系列评估（~20 个）** — 与 §10.3.2 的 C4 主题文档呼应，但格式更结构化（含 metrics.json/
metrics.tsv 数据文件）。`validate-current-c4.md` 的结论 "better than prior SYNG-derived C4 outputs
on the specific repeat-artifact failure mode, but still not solved and not yet..."反映 C4
难用例至今未完全解决。

**implement-* 系列（agent 评分普遍低）** — 4 个 implement 文件中 3 个评分 0.00/1.00，唯一非零分
是 `implement-occurrence-level.md`（0.20/1.00，confidence 0.78）。这反映 agent 评分系统与任务
实际完成度的脱节——rubric underspecified 是常态。

**`local-compression-autopoietic/` — 自循环迭代测试** — "autopoietic"（自生系统）指自动反馈循环：
每轮运行 → 分析候选 → 综合 → 进入下一轮。含 iter-1/iter-2/iter-2b 三轮的 hypothesis/change/
result/next recommendation 滚动状态表（`summary.md`）。这种"滚动状态表"是 agent-assisted 开发中
保持迭代上下文的重要工具。

**`local-compression-testbed-fast/` — 测试床 fast profile** — 含 `scoreboard.json`（每 fixture
×method 的 `fixture_class`/`tier`/`method_family`/`method_parameters`/`command_line`/
`output_gfa_path` 字段）+ `fixture-validation.json`（12 个 validated fixtures）+ `fixtures/`
子目录（每 fixture 一份 input.fa/expected_paths.tsv/metadata.json + 各 method 的 output.gfa）。
这是"多方法 × 多场景"回归矩阵的完整可执行产物，与 §6.7 的测试体系对应。

#### 10.3.6 文档命名约定与特点

通读 109 个文件后，核心模式可归纳为三点：

1. **前缀主题化 + 内容四分类** — 文件用 `主题-子主题.md` 命名（如 `crush-aligner-deep-diag.md`），
   便于 `ls crush-*` 列出同主题。按内容（而非仅文件名）可分为四类：设计/规范（`*-design.md`/
   `*-spec.md`/`*-architecture.md`，前瞻性）、实验报告（`*-exp-*.md`/`*-sweep.md`/`run-*.md`/
   `compare-*.md`，含 hypothesis + TL;DR + 实测数据）、诊断/审计（`*-diag*.md`/`*-audit.md`/
   `diagnose-*.md`/`*-blocker-*.md`，回溯性）、实施/评估（`implement-*.md`/`evaluate-*.md`/
   `validate-*.md`，agent 系统产物）。**前缀不总等于主题**：`c4-crush-handoff.md` 前缀是 c4，
   内容是 crush 交接。
2. **无顶层索引，靠前缀自组织** — `docs/` 顶层无 `README.md` 或索引文件，新读者需要 `ls | sort`
   后按前缀浏览。`README.md` 仅显式引用 3 个文档（`syng-gfa-query.md`/`genotype-architecture.md`/
   `infer-design.md`），其余 100+ 文档是"暗物质"——只对开发者可见。这是 impg 文档体系最大的可改进点。
3. **agent-assisted 工作流痕迹** — `evaluations/` 下多数文件以 `Task: <id>` + `Evaluator: agent-<N>`
    - `Overall score: X.XX / 1.00` + `Rubric underspecified: true/false` 开头；实验文档记录
      `Branch: wg/agent-<N>/<topic>` 与 binary commit hash 便于复现；artifact 路径硬编码
      `/home/erikg/impg/data/<dir>/` 与 `hypervolu.me/~erik/impg/<png>`。这些是
      single-developer (Erik Garrison) + agent-assisted 工作流的产物。

#### 10.3.7 与 `pgr` 文档实践的对比

`pgr` 的 `docs/` 是**面向用户的设计文档**（每格式/每模块一篇，正文英文， AGENTS.md 第 10
行约束），数量受控（~25 篇），每篇对应一个 `pgr <command>`或核心算法。`impg` 的 `docs/` 是
**面向开发者的工作笔记**（英文为主、大量实验报告、无 README 索引、含 agent 评分产物），数量爆炸
（109 篇 + evaluations / 子树），反映其活跃开发节奏与"先写文档再写代码"+"agent-assisted"的工程文化。

- **可借鉴**：impg 的"前缀主题化 + 实验/设计/诊断三分类"命名约定，在 `pgr` 项目规模扩大时值得参考。
- **应避免**：impg 缺少顶层索引，新读者难以切入。`pgr` 的 `AGENTS.md`
  中的"关键设计文档"列表实际充当了 docs/ 的索引，应继续保持。
- **可借鉴**：impg 把"实验结果"（参数 sweep、性能基准、agent 评分）也纳入 docs / 而非丢弃，
  便于回溯。`pgr` 的 `benches/` 当前只有基准代码，没有结果报告，可考虑增加 `benches/results/`
  子目录归档历次基准结果。
- **应避免**：impg 的 `evaluations/` 中大量 `Overall score: 0.00/1.00` 文件， 反映 agent
  评分系统与任务实际完成度的脱节。`pgr` 若引入类似 agent 工作流，应确保评分 rubric 明确（避免
  `Rubric underspecified: true`）。
- **可借鉴**：impg 的 `crush-architecture-spec.md` 模式——当某个算法复杂到"three days of agents
  working on crush have been flailing without a clear spec"时，强制写权威 spec 并要求 "No code
  change should be proposed without referencing this document"。`pgr` 的 `docs/spoa_port.md`
  等移植笔记已有类似性质，可强化。

