# pgr 泛基因组 — PAF 隐式图

**复用已有的 pairwise 比对基础设施，构建 PAF 隐式图，按需回答"哪些序列的哪些区段同源"， 而非物化一张泛基因组图。**

V1-V6（graph-report）已全部完成。综合参考：impg（POA → GFA）、minigraph（chain → CIGAR → rGFA）、
seqwish（DSU 传递闭包）、Cactus（Caf 退火-熔化）。

参考文档：[[impg.md]]、[[seqwish.md]]、[[minigraph.md]]、[[cactus.md]]、[[cactus_lastz.md]]、
[[project-understanding.md]]、[[ecoli-cohort.md]]。

```
索引层 ✅  |  查询层 ✅  |  图构建层 ✅（V4a 粗全局 + V4b 局部精细 + V5 VCF + V6 graph-report）|  应用层 ← 远期
```

---
## 1. 路线决策

### 1.1 四条原则

1. **隐式图优先，粗 GFA 作为可选投影** — 默认用 PAF 索引 + 区间树 + BFS 传递闭包做"按需图遍历"，
   不物化 GFA。粗全局 GFA 是索引的**显式投影**（V4a），数据源仍是 PAF 索引。粒度差异：seqwish
   传递闭包是全局一次性，pgr 是局部按需（每次查询从一个区间出发 BFS）——见 [[seqwish.md]] §5。
2. **复用 pairwise 资产，大 cohort 用 Mash KNN sparsify** — pgr 已有成熟的 pairwise 比对链
   （`pgr lav lastz` → chain → net → axt → maf）。小 cohort + 已有 MAF 直接复用；大 cohort + 无先验
   （如 4 万 E. coli）用 Mash KNN sparsify 把 N² 降到 N×K，传递闭包推断未比对的对。
3. **查询层全量，图构建层粗框架** — 查询层（V1-V3）全量返回同源区段，用户用 `--merge-distance`
   控制粗细；图构建层（V4a+）物化 GFA 时才引入 `--min-var-len`（默认 100）粗框架过滤，对齐
   minigraph。
4. **pipe 友好，两段式 GFA** — V1 默认 PAF 输出、`-o bed` 可选（最 pipe 友好），
   `bed → fa range → fas consensus` 的 MSA 路径已通。GFA 两段式：粗全局（地图）+ 区域精细（碱基级）。

### 1.2 pgr 与 impg 的起点差异

pgr 走向泛基因组时，面对的问题与 impg **完全不同**。impg 的起点是"只有 FASTA，没有 pairwise 比对"，
需要先选对、再比对、再索引。pgr 的起点是"已有成熟的 pairwise 比对基础设施"，需要的是"复用已有资产，
补上缺失的图遍历层"。

| 维度     | impg                       | pgr                              |
|----------|----------------------------|----------------------------------|
| 比对来源 | 从 FASTA 跑 wfmash/sweepga | 已有两序列 MAF（可转 PAF）       |
| 挑选时机 | align 阶段（无先验）       | 可借已有 MAF 先验                |
| 核心问题 | 选哪些对比对               | 复用已有 pairwise，做 PAF 隐式图 |
| 比对工具 | wfmash/FastGA              | pgr 已有 `pgr lav lastz` 全套    |

### 1.3 与 `--sparsify` 的关系：分场景

impg 的 `--sparsify auto` 用 Mash KNN 从 N 个基因组中选 K 个近邻做比对，把 N² 降到 N×K。
**这是隐式图架构避免 N² 爆炸的核心机制**——稀疏比对 + 查询时 BFS 传递闭包推断未比对的对。

pgr 是否需要 sparsify **取决于 cohort 规模和已有资产**：

| 场景                     | 规模          | 已有资产          | 是否需要 sparsify                       |
|--------------------------|---------------|-------------------|-----------------------------------------|
| **小 cohort + 已有 MAF** | 几十基因组    | 现成 pairwise MAF | **不需要** — MAF 里的对已跑过比对       |
| **大 cohort + 无先验**   | 27000 E. coli | 只有 FASTA        | **必需** — 否则 27000² ≈ 3.6 亿对不可行 |

**小 cohort**：不需要选对、不需要 wfmash；"挑选"发生在查询层（`--min-identity` 等参数过滤 PAF
记录）。

**大 cohort**：Mash KNN 把 27000² 降到 27000×K（K≈50，~270 倍缩减）；用 FastGA 比对产 MAF → 转 PAF；
稀疏比对的缺口由查询时 BFS 推断。

> **spanning tree 优化（远期）**：seqwish 在传递闭包前用最大权生成树剪枝，把 N(N-1)/2 边压缩到 N-1
> 边。pgr 查询层 BFS 若性能瓶颈显现，可在加载 PAF 阶段预计算生成树。当前不做，待性能数据出来再评估。

### 1.4 三种图构建路线

| 路线          | 流程                                          | 输出                 |       pgr 采用       |
|---------------|-----------------------------------------------|----------------------|:--------------------:|
| **impg**      | BFS 传递闭包 → POA → GFA → gfaffix → gfasort  | 标准 GFA 1.0         |       V4b 参考       |
| **minigraph** | k-mer seeds → linear chain DP → gchain → rGFA | rGFA 1.0（SN/SO/SR） | V4a 参考（DSU 风格） |
| **pgr**       | MAF → PAF → PafIndex → BFS → POA → MSA/GFA    | MAF/GFA/VCF          |      ✅ 已实现       |

pgr 的独特起点：已有 PAF index + BFS 查询，不需要重新做比对；`libs/poa/` 纯 Rust POA，零外部依赖。

---
## 2. 核心决策

以下决策是后续所有行动的**不变前提**。

### 2.1 用 PAF 作隐式图边集，不用 Chain

PAF 是图的边，Chain/Net 是查询层的 syntenic 过滤器。理由（详见 [[impg.md]] §9.1）：

- Chain 是 star topology（ref↔query_i），做传递闭包时 ref 成为必经枢纽，ref 缺失区段会断开间接同源路径
- Chain 已被 UCSC 流程过滤（score 阈值、syntenic 净化），不是原始比对，丢失了 paralog/低质量区间
- Chain 的 gap-less tBlock 分段结构不适合做图边

**Chain/Net 的正确角色**：PAF 边集提供"所有可能同源"（全量装入）；Chain/Net 提供 syntenic 验证
（已落地为 `--syntenic-filter`）。这是 pgr 独有的优势——impg 没有 UCSC Chain/Net 体系。

### 2.2 PAF 来源：MAF → PAF 转换，不跑新比对

pgr 已有的两序列 MAF 直接转换为 PAF。两序列 MAF 的每个 block 等价于一条 pairwise alignment—— `s`
行给出坐标和链向，可直接映射到 PAF 的 12 列。

### 2.3 索引全量装入，挑选发生在查询层

PAF 索引时不做过滤，所有记录全量装入区间树。过滤参数只在查询时生效。同一份索引可服务不同严格度的
查询。对应 impg 的 `Index` 命令只有文件路径和 index-mode 参数，而 `QueryOpts` 才有过滤开关。

**索引层工程优化（4 万大肠杆菌规模可考虑）**：seqwish 的 `SeqIndex` 用 FM-index 索引序列名 +
`SparseBitVec` 记录序列边界，比 HashMap 省内存；`PosT` 把 offset+方向打包进单 u64。详见 [[seqwish.md]
] §2.1、§2.2。pgr 当前用 HashMap + coitrees 在 4 万大肠杆菌规模尚可，若扩到 HPRC 规模（数百单倍型、
Gb 级），可借鉴此方案。

### 2.4 传递闭包是图遍历，不是多序列比对

传递闭包做"图遍历可达性"，不产出多重比对。找到所有同源片段后，如需 MSA，再调用 `fas consensus`（SPOA）
或 `fas multiz`（banded DP）。图遍历和 MSA 是正交的两个步骤，不应耦合。

pgr 的 MSA 质量可能优于 impg 的 per-bubble POA——`fas_multiz.rs` 实现了 banded DP 合并， 对 core
区段比纯 POA 更精确。

### 2.5 两段式 GFA，局部不合并回全局

- **两段式 GFA**（V4）：粗全局 GFA（V4a，提供"地图"）+ 区域精细 GFA（V4b，碱基级）。 比 minigraph
  （只有粗）和 impg（只有精）更完整。V4a 与 V4b 互不依赖。
- **局部 GFA 不合并回全局** — 保持"粗全局 GFA 是不可变投影"语义，不滑向"可变全局图" （等于重建
  minigraph `gfa_t` + augment）。需要全局精细视图走 V5 VCF/MAF（天然可 concat）。

### 2.6 V5 跳过 GFA 物化

BFS 等价类 + POA 直接产出 MSA → VCF，省去 query→GFA→VCF。这是 pgr 与 impg 的关键差异。

### 2.7 三层挑选问题

pgr 的"挑选"不是 impg 的"选哪些对跑比对"，而是分三层：

**第一层：从已有 MAF/PAF 挑选（查询层，无需新比对）** — impg 的传递闭包（`-x` BFS）： 若 A↔B、
B↔C 在同一区段有比对，则 A↔C 间接同源。所有 pairwise 比对当作图的边集，从目标区间出发 做 BFS，
自动发现所有直接和间接同源片段。**这一层不需要新比对**，只需把已有 MAF 转成 PAF 装入 区间树，
做查询层挑选。

**第二层：补充 pairwise 比对的挑选（align 层）** — 已有 MAF 只覆盖已跑过的对。以下场景 MAF 缺失
或不足：cohort 加入新基因组、已有 MAF 在某区段断开、某些 sample 对需要更精细的 region-level 重比对。

| 策略                | 来源                   | 适用                   | pgr 实现门槛           |
|---------------------|------------------------|------------------------|------------------------|
| 已有 PAF 覆盖度先验 | pgr 独有               | 已有部分 PAF 的 cohort | **推荐**，复用已有 PAF |
| `pgr lav lastz`     | pgr + Cactus 风格      | 特定 pair 需要新比对   | 已有（不含 `--self`）  |
| 系统发育树引导      | Cactus 风格            | 有 phylogeny           | 复用 `pgr nwk` 模块    |
| Mash KNN            | impg `--sparsify auto` | 无先验全选             | 需引入 mash crate      |

**已有 PAF 覆盖度先验策略**（pgr 推荐）：对每个 query_i 统计其在已有 PAF 上的覆盖区间集合 C_i； 对
query_i、query_j 计算 |C_i ∩ C_j| / |C_i ∪ C_j|（Jaccard）；选 Jaccard 高于阈值且尚未跑过 pairwise
的对补充比对。这样把 N² 降到"PAF 覆盖度共享的子集"。

> `pgr lav lastz --self` 是 Cactus 风格的**重复屏蔽**管道的一部分（碎片自比对检测基因组内重复），
> **不是**泛基因组比对工具。泛基因组 pairwise 用 `pgr lav lastz`（不含 `--self`）。
> 详见 [[cactus_lastz.md]] §5.6。

**第三层：region 级精细比对挑选** — 已有 MAF 是粗粒度的。某些 region（HLA、KIR、C4）需要更精细
pairwise，但全基因组精细比对代价高。从已有 PAF 的 gap/low-identity 区段筛选候选 region，对候选
region 跑 `pgr lav lastz`，合并回 PAF 网络。**这一层是第一层的补充**，不是泛基因组的核心路径，
按需开启。

### 2.8 paf query 的输出格式策略

**PAF 是默认输出，BED 为 `-o bed` 可选**：impg `query` 默认 `-o bed`，pgr V1 选择**PAF 为默认**
（含 CIGAR/gi/bi 完整比对记录），BED 通过 `-o bed` 可选——理由是 pgr 既有测试已断言 PAF 输出，且 PAF
对需要完整比对记录的场景更直接。BED 三列（`name start end`）是坐标投影的轻量产物，最 pipe 友好，用
`-o bed` 显式切换。

**FAS（block FASTA with shared coords）不输出**：FAS 格式的核心假设是所有序列共享一个统一坐标系
（通常以 reference 为锚），这在泛基因组场景不成立——PAF query 结果是各基因组**独立坐标系**下的
同源区段列表。

**fasta/maf 是可选附加，按依赖链后移**：impg 的 11 种输出按"是否需要序列文件"分两类——坐标类
（`bed`/`bedpe`/`paf`，不需 `-f`）是核心，序列类/MSA 类/图类（`fasta`/`maf`/`gfa` 等，需 `-f`）
是可选。pgr 按此分阶段：V1 坐标输出（PAF 默认，BED 可选）→ V2 未比对序列（需 `-f`）→V3 POA MSA（需
`-f`）。

---
## 3. 已实现能力（V1-V5）

| 阶段       | 命令                       | 实现                                       | 关键点                                                                                     |
|------------|----------------------------|--------------------------------------------|--------------------------------------------------------------------------------------------|
| **V1** ✅  | `pgr paf query` / `to-bed` | PAF 默认输出 + BED3 轻量坐标 + `-b` 批查   | PAF 含完整 CIGAR/gi/bi；BED 由独立子命令提供                                               |
| **V2** ✅  | `pgr paf to-maf`           | pairwise MAF（按 CIGAR 还原，需 `-f TSV`） | 不做 refine（上游 chain/net 已优化）；`-` 链 RC 处理                                       |
| **V3** ✅  | `pgr paf to-maf --msa`     | POA 多序列 MSA（需 `--transitive`）        | 复用 `libs/poa/`；target 第一条，queries 按 result 顺序                                    |
| **V4a** ✅ | `pgr paf graph -f refs.fa` | 粗全局 GFA（seqwish DSU 风格）             | `--min-var-len 100`；CIGAR 切分 → 段对 DSU → 节点序列 → 路径 + novel 段补全 → 边派生       |
| **V4b** ✅ | `pgr paf to-gfa`           | 区域精细 GFA（impg 风格）                  | unchop 默认开；`--crush` 可选 bubble 压缩；LN tag；多 region 独立                          |
| **V5** ✅  | `pgr paf to-vcf`           | POA MSA → VCF（SNP + INS/DEL）             | 复用 V3 `build_msa_entries`；三分支主循环；1bp anchor；indel 左对齐（`left_align_indels`） |
| **V6** ✅  | `pgr paf graph-report`     | 粗图拓扑报告（25 维度 TSV）                | 复用 V4a `PafGraph::build`；节点长度/覆盖分布 + 连通分量 + tips/self-loop + 路径长度分布   |

配套命令：`pgr maf to-paf`（MAF → PAF 转换）、`pgr paf index`（区间树索引，支持多文件合并
`build_multi`、`.paf.idx` 持久化、BGZF lazy CIGAR）。

### 3.1 增量增强（前三项已落地，第四项语义合并）

| 增强                   | 状态                                                                                    |
|------------------------|-----------------------------------------------------------------------------------------|
| indel VCF              | ✅ SNP + INS/DEL，1bp anchor，indel 左对齐（`left_align_indels`）                       |
| `-m/--max-depth`       | ✅ BFS 深度控制（默认 2，0=unlimited，5 子命令共享 `add_query_args`）                   |
| `--syntenic-filter`    | ✅ chain-level query 侧覆盖检查（装饰器，5 子命令共享 `add_query_args`）                |
| `--min-transitive-len` | ⚠️ 未单列选项；语义由 `--min-len` 覆盖（默认 10，impg `--min-transitive-len` 默认 101） |

### 3.2 V5 VCF 已知限制

- **1bp 锚定边界**：indel 位于 MSA 首列或前导列也 gap 时被跳过（VCF 要求 REF 非空）。
- **无 phasing**：GT 是单倍体（0/1/.）。
- **REF = target**：target 序列作为 REF（与 V3 MSA 的 target 选择一致）。

### 3.3 V4a 简化项（相对 seqwish）

无 disk-backed interval tree / SparseBitVec / lock-free DSU；路径方向恒 `+`（反向链段已翻转坐标）；
rGFA SN/SO/SR tag 已补全（S 行，origin 取 DSU 等价类中 `(seq_id, start)` 最小者）。

### 3.4 V4b `--crush` 边界

impg crush 8 阶段流程的最小子集（仅"邻居集合相同的节点合并"）。完整 15 种 ResolutionMethod 不做——
复杂 bubble 走外部工具（impg/smoothxg）。

---
## 4. 代码结构

### 4.1 模块组织

```
src/libs/paf/
├── mod.rs          # 模块导出
├── record.rs       # PafRecord — String 字段 + tags
├── parser.rs       # 纯文本 PAF 解析
├── cigar.rs        # CigarOp bit-packing + stats + identity
├── writer.rs       # PAF 行格式化
├── index.rs        # PafIndex + PafMetadata + SortedRanges + BFS
├── graph.rs        # V4a 粗全局 GFA 引擎（DSU 传递闭包）
└── persist.rs      # .paf.idx 磁盘持久化（bincode）

src/cmd_pgr/paf/
├── index.rs / query.rs / to_bed.rs / to_maf.rs / to_gfa.rs / to_vcf.rs / graph.rs
└── mod.rs
```

名字映射用 pgr 的 `IndexMap<String, u32>` 模式（与 `libs/loc.rs`、`libs/phylo/tree.rs` 一致），
不需要独立的 `SequenceIndex`。

### 4.2 PafRecord

string-based 设计（`src/libs/paf/record.rs`）：

```rust
pub struct PafRecord {
    pub query_name: String,   pub query_length: u32,
    pub query_start: u32,     pub query_end: u32,
    pub strand: char,
    pub target_name: String,  pub target_length: u32,
    pub target_start: u32,    pub target_end: u32,
    pub matches: u32,         pub block_length: u32,
    pub mapq: u8,
    pub tags: Vec<String>,    // gi:f:..., cg:Z:... 等
}
```

区间树节点使用独立紧凑 `PafMetadata`（u32 坐标 + `CigarStore` 引用），不存序列名。

### 4.3 PAF 解析与 BGZF 支持

双路径设计，依据输入文件类型自动 dispatch（`build_from_path`）：

- **Plain text / 普通 gzip**：走 `pgr::reader`（`src/libs/io.rs`）的 `flate2::read::MultiGzDecoder`，
  透明处理 plain gzip、multi-member gzip 和 BGZF。CIGAR 全量驻留内存（`CigarStore::Owned`）。
- **BGZF（lazy CIGAR）**：`is_bgzf()`（18 字节头检测：`1f 8b 08 04` + `BC` subfield）判定为 BGZF
  时，走 `build_lazy_bgzf()`，用 `noodles_bgzf::io::Reader<File>` 读取。每条 PAF 行的 BGZF virtual
  position（u64，高 33 位块偏移 + 低 16 位块内偏移）记录在 `CigarStore::Lazy(vpos)` 中；CIGAR
  字符串在查询时按需 seek + parse（`fetch_cigar`）。

用 `noodles_bgzf::IndexedReader` 而非 impg 的外部 `.gzi` 索引方案——自带 BGZF 块索引， seek
时直接定位块边界，无需二进制索引文件。

### 4.4 CIGAR 编解码

`CigarOp` bit-packing：3 位 op code + 29 位 length 压入单个 `u32`，支持 `=`/`X`/`I`/`D`/`M` 五种
op，单段最长 512Mbp。`parse_cigar("10=5X2I3D") → Vec<CigarOp>`，`format_cigar` 反向。

```rust
pub enum CigarStore {
    Owned(Vec<CigarOp>),
    Lazy(u64),              // BGZF virtual position
    LazyReversed(u64),      // mirror index entry, I/D swapped
}
```

`resolve_cigar(&CigarStore)` 是统一入口，Owned 直接 clone，Lazy 走 `fetch_cigar`。

identity 计算：对 CIGAR 做 fold 统计——`=`/`X` 计入 matches/mismatches，`I`/`D` 按事件计数（gi）
或按碱基计数（bi）。`gi` 评估同源性（对长 indel 宽容），`bi` 评估序列一致性（对长 indel 严格）。

PAF 输出：`write_paf_record` 输出 12 列 + 四个标准标签：`gi:f:`（gap-compressed identity）、`bi:f:`
（block identity）、`cg:Z:`（CIGAR string）、`an:Z:`（alignment name，可选）。

### 4.5 PafIndex 设计

```rust
pub struct PafIndex {
    pub names: IndexMap<String, u32>,
    pub(crate) trees: HashMap<u32, Arc<BasicCOITree<PafMetadata, u32>>>,
    pub(crate) reverse_trees: HashMap<u32, Arc<BasicCOITree<PafMetadata, u32>>>,
    pub(crate) lazy_source: Option<Mutex<bgzf::io::Reader<File>>>,
    pub(crate) lazy_source_path: Option<String>,
}
```

**关键设计**：

- **`IndexMap` 而非 impg 的 `SequenceIndex`**：与 pgr 既有 `loc.rs`、`phylo/tree.rs` 风格一致。
- **不用 `RwLock`**：V1 单线程构建+查询。V2 rayon 化时再包 `Arc<RwLock<>>`。
- **`Arc<BasicCOITree>`**：查询回调借用方便（非 disk-backed 需求）。
- **双向索引 `reverse_trees`**：`+` 链 record 在 `trees[target_id]` 和 `reverse_trees[query_id]`
  各插一条，mirror 条目 CIGAR 经 `reverse_cigar` 反转并交换 I/D。BFS 可双向传播。负链不建 mirror。
- **不需要 `ForestMap`**：纯内存，bincode 整体持久化（v3 格式含 `reverse_intervals`）。
- **不需要 `ImpgIndex` trait**：V1 只有一种 `PafIndex`，无 `MultiPafIndex`。
- **lazy CIGAR**：大 cohort 场景 CIGAR 全量驻留内存成本过高（4 万大肠杆菌 × 5K 比对 × 100bp CIGAR ≈
  数 GB），故 BGZF 输入走 `CigarStore::Lazy(vpos)`。

**构建**：`build_from_path` 检测 BGZF → `build_lazy_bgzf` 或 `build`。records 按 `target_id` 分组，
每组建一棵 `BasicCOITree`，串行（V2 可 rayon 化）。

**单跳查询**：`query(target_id, start, end)` → `tree.query_intersecting(range, callback)`。
回调计算 query 侧坐标投影。不返回 self-entry，不区分 normal/approximate 模式（默认不读序列=近似）。

**传递闭包 BFS**：`query_transitive_bfs(target_id, range, max_depth, ...)`。用
`HashMap<u32, SortedRanges>` 去重，每轮只把"未被已有区间覆盖的新增部分"入队。每轮对 `tid` 同时查询
`trees[tid]`（forward）和 `reverse_trees[tid]`（mirror），`Arc::clone` 取出后 `chain` 遍历。
`max_depth` 控制深度（默认 2，0=unlimited）。

**多文件索引**：`build_multi` 直接合并多 reader 的 records，统一 `IndexMap` + per-target 区间合并。
不需要 impg 的 `local_to_unified` 翻译表 / `RwLock` 缓存 / `.multi_impg` 缓存。

### 4.6 错误类型与依赖

`PafParseError` 枚举（`NotEnoughFields`/`InvalidInteger`/`InvalidStrand`/
`InvalidCigarFormat`/`InvalidFormat`/`IoError`），实现 `Display + Error`，在 `execute` 中用
`.map_err(|e| anyhow!("PAF: {}", e))` 桥接到 `anyhow::Result`。

唯一新增依赖是 `coitrees`（区间树）。pgr 已有 noodles-bgzf、rayon、serde/bincode 等基础设施。

---
## 5. 后续规划（V6-V8）

V1-V5 已完成索引→查询→图构建三层的最小可用闭环。后续按"图质量 → 规模扩展 → 应用层"递进。

### 5.0 近期打磨项（可选，V6 之前）

V1-V5 已可用，以下是审查中识别的零散收尾项。**均为可选**，按触发条件启动，不做预防性开发：

| 打磨项                      | 位置          | 触发条件                              | 工作量 |
|-----------------------------|---------------|---------------------------------------|--------|
| V4a rGFA tag（SN/SO/SR） ✅ | §3.3 简化项   | 需与 minigraph / odgi 工具链互操作时  | 小     |
| VCF 左对齐 ✅               | §3.2 已知限制 | 用户反馈 `bcftools norm` 后处理不够时 | 中     |
| `--min-tree-coverage`       | §6.4 Caf 维度 | 有 phylogeny 上下文且需按树分布过滤时 | 中     |

> `pgr paf` 基准测试（`benches/`）原列于此，因依赖 4 万大肠杆菌真实数据，已随 V7 规模扩展一并推迟。

**V4a rGFA tag** ✅：已在 `PafGraph::write_gfa` 的 S 行追加 `SN:Z`（源序列名）、`SO:i`（0-based
起始 偏移）、`SR:i:0`（rank 0 = primary）。origin 取 DSU 等价类中 `(seq_id, start)` 最小者（PAF
target 先 注册，故 target 优先），novel 节点 origin 取其填充时的 `(name, start)`。SR 暂恒为 0（pgr
路径方向 恒 `+`，无链翻转）。与 minigraph rGFA / odgi 工具链兼容。

**VCF 左对齐** ✅：在 `to_vcf.rs` 中实现 `left_align_indels` 辅助函数，对 INS/DEL 做左推。取 target
序列前 1000bp 前缀构建 `target_ext`，当锚点前碱基与所有非空 indel 序列末位碱基相同时，将锚点左移并
相应调整 indel 序列。POA MSA 通常已将 gap 左对齐到重复序列边界，`left_align_indels` 在此基础上做
二次规范化，确保不依赖 `bcftools norm` 后处理。

**基准测试**：`Cargo.toml` 仅有 `hier_benchmark`。V7 扩规模前应补 `paf_index_bench`（构建/查询
/BFS）和 `paf_graph_bench`（V4a DSU），为"选路径 A 哪几项优化"提供量化依据。基准测试依赖 4
万大肠杆菌 真实数据，目前不可用，随 V7 一并推迟。

### 5.0.1 测试改进（借鉴 impg/seqwish）

对照 impg `test_transitive_integrity.rs`（8 个传递闭包不变量）和 seqwish
`integration_test.rs` / `test/HLA/`（黄金回归），pgr paf 测试套件的覆盖存在以下盲区
（测试已按子命令拆分到独立文件，见文末文件结构表）：

| 改进项                              | 来源                          | pgr 现状              | 优先级 |
|-------------------------------------|-------------------------------|-----------------------|--------|
| 传递闭包不变量（4 场景）            | impg test 1/2/5/7             | ✅ 已落地（4 测试）   | 高     |
| GFA 路径拼写 round-trip             | impg `path_sequences` 校验    | ✅ 已落地（3 测试）   | 高     |
| indel 坐标在 query 层精度           | impg test 6                   | ✅ 已落地（2 测试）   | 中     |
| 黄金回归（固定输入 + md5 校验）     | seqwish HLA 30 基因座         | ❌ 无基础设施         | 低（暂不做） |

**传递闭包不变量** ✅（impg `test_transitive_integrity.rs`，每个用 2-4 行极小 PAF 构造精确拓扑）：

1. ✅ `command_paf_transitive_non_overlapping_regions_stay_separate` — A:0-100→B 和 A:500-600→C
   不互染（查询 A:0-100 不应找到 C）。pgr 现有测试只验证"能找到"，不验证"不该找到的没找到"。
2. ✅ `command_paf_transitive_coordinate_accuracy_subregion` — A:25-75 经传递闭包应投影到
   B:25-75，而非整段 B:0-100。pgr 现有 `command_paf_query_transitive` 只查整段（B:0-100），
   未验证子区间精度。
3. ✅ `command_paf_transitive_distant_regions_no_collapse` — 多 hop 不应把远距离区域错误连通
   （A:0-100→B→D:0-100 与 A:1000-1100→C→D:500-600 应保持分离）。pgr 无此测试。
4. ✅ `command_paf_transitive_multiple_alignments_to_same_target_stay_separate` — 同 target 的
   两个比对（A→B:0-100 和 A→B:500-600）应报告为两条独立结果，不合并。pgr 无此测试。

测试位置：`tests/cli_paf_query_bfs.rs`（4 个 `command_paf_transitive_*` 函数，均用 `to-bed` +
`--transitive` 子命令验证）。

**GFA 路径拼写 round-trip** ✅（impg `test_graph_poa.rs` / `test_graph_seqwish.rs`）：解析 GFA P 行
拼接路径序列，与输入 FASTA 逐字比对。这是图构建正确性的最强不变量（路径必须还原原序列）。
已落地 3 个测试（`command_paf_to_gfa_roundtrip_*`，位于 `tests/cli_paf_to_gfa.rs`）：
identical / SNP bubble / 2bp insertion bubble，
共用 `spell_gfa_paths` 辅助函数（解析 S/P 行，按 `+/-` orientation 拼写，`-` 时做 reverse-complement）。

**indel query 层精度** ✅（impg test 6）：验证 CIGAR 含 indel 时 query 投影坐标不漂移。已落地
2 个测试（`command_paf_to_bed_*_coordinate_accuracy`，位于 `tests/cli_paf_to_bed.rs`）：
insertion（`50=10I50=`，验证 insertion 边界
两侧的 B 子区间正确投影到 A）+ deletion（`50=10D50=`，验证 deletion 边界两侧的 B 子区间正确投影到 A）。
发现：query 正好在 indel 边界时，pgr 投影会包含相邻 indel 碱基（合理行为），用边界内移 10bp 的
子区间做精确 start 坐标验证。

**黄金回归**：seqwish `test/HLA/` 用 30 个 HLA 基因座（固定 fa + paf + 预计算 gfa.md5）做回归。
pgr 无固定测试数据集和 md5 校验基础设施，工作量大，暂不列入。

**大小写敏感性**（借鉴 seqwish `01_seqwish.t` 的 masked-sequence 测试）：
seqwish 对 lowercase/masked 输入做归一化，断言输出 GFA 的 md5 与 uppercase 版相同。
pgr 的 POA 不归一化（`poa/align.rs` 用 `seq_base == node_base` 严格比较），全小写输入会产生
不同拓扑（SNP bubble 场景：uppercase 4 segments vs lowercase 6 segments）。
因此 seqwish 的"md5 等价"不变量在 pgr 中不成立，改为验证 round-trip 不变量：小写输入下
每条 path 拼写出的序列仍等于输入小写序列（`command_paf_to_gfa_lowercase_roundtrip`，
位于 `tests/cli_paf_to_gfa.rs`）。这是 pgr 与 seqwish 的已知行为差异，若未来需要
case-insensitive 语义，应在 `to-gfa` 输入路径加 `to_ascii_uppercase()`，而非改 POA。

### 5.1 V6：图质量与归一化

| 能力                      | 来源                            | pgr 价值             | 优先级                             |
|---------------------------|---------------------------------|----------------------|------------------------------------|
| `pgr paf graph-report` ✅ | impg `graph report`（15+ 维度） | 图质量评估，无新算法 | 高                                 |
| `--gfaffix`（图归一化）   | impg `run_gfaffix`              | 重复 bubble 标准化   | 低（外部 binary）                  |
| 完整 crush 8 阶段         | impg `resolution.rs`            | 复杂 bubble 处理     | 低（V4b `--crush` 已覆盖主要场景） |

`graph-report` 已实现（25 维度 TSV：基础拓扑 + 节点长度/覆盖分布 + 连通分量 + tips/isolated/
self-loop + 路径长度分布），聚焦 impg 60 字段中可从 PafGraph 直接计算的核心子集，不依赖 povu/
flumble 等重型设施。完整 crush 不做。

**`--smooth` 不做**：impg `smooth.rs` 的核心是"对已有 GFA 重新做 POA 对齐"（smoothxg 块分解 + SPOA +
lace）。但 pgr 的 V4b `to-gfa` 已是 POA MSA 的直接产物，再套 smooth 的"重新 POA"是冗余的；对 V4a
粗图做 smooth 又与 V4b 功能重叠。经评估边际价值不明确，从计划中移除。

### 5.2 V7：规模扩展（4 万大肠杆菌级，待数据可用）

当前 V4a 在小 cohort（数十基因组）验证通过。扩到 [[ecoli-cohort.md]] 的 4 万大肠杆菌需要真实 cohort
数据（Mash 去冗余 + FastGA sparsify 产出的 PAF），**目前不可用，本阶段整体推迟**。数据就绪后两条
互补路径：

**路径 A：单图工程优化**（按收益/复杂度排序）

| 优化                                 | 来源                | 何时引入                       | 估计收益                     |
|--------------------------------------|---------------------|--------------------------------|------------------------------|
| spanning tree 剪枝                   | [[seqwish.md]] §3.1 | BFS 边遍历数成为瓶颈时         | N²→N-1 边遍历，BFS 提速 10×+ |
| lock-free DSU                        | [[seqwish.md]] §2.4 | V4a rayon 化时                 | 并行 union-find              |
| `--repeat-max` / `--min-repeat-dist` | [[seqwish.md]] §3.5 | 高拷贝重复区把图吹爆时         | 限制同序列在同节点的拷贝数   |
| disk-backed interval tree            | [[seqwish.md]] §2.3 | 全图超 RAM 时                  | 兜底可跑性                   |
| SparseBitVec 序列边界                | [[seqwish.md]] §2.2 | HPRC 规模（数百单倍型、Gb 级） | 边界索引内存 O(m)            |
| PosT 单 u64 编码                     | [[seqwish.md]] §2.1 | 反链投影成为热点时             | 单棵树同时表达正反链         |

**路径 B：分区分批构建**（impg partition+lace，[[impg.md]] §6.3）

`partition_size: Option<usize>` 开关式触发（同一个 `pgr paf graph` 命令切换直接构建 vs 分区构建）。
复用现有 `query_transitive_bfs`；lace 拼回需新增 `pgr paf lace` 子命令。

**判断标准**：数据就绪后先跑 4 万大肠杆菌基准测试（无优化版本）。单图能跑通但慢 → 选路径 A 的 1-2
项；单图超 RAM → 选路径 B。不做"预防性优化"。

### 5.3 V8：应用层 — genotyping（远期）

impg 能力栈顶端是 `Genotype` / `Infer`（[[impg.md]] §7）。pgr 路径：复用 V4a 粗全局 GFA + V5 VCF，
加 `pgr paf genotype`。推迟理由：genotyping 需要稳定的图构建（V4-V7）前置，且需要真实数据验证。
当前 pgr 核心价值在"图构建 + 区域查询"，应用层让位于 minigraph-cactus / vg 等成熟工具。

### 5.4 优先级路线图

```
当前 (V1-V6 graph-report ✅ + 增量增强 ✅* )
   │
   ├─ §5.0 打磨项（可选，按触发条件）: rGFA tag ✅ / VCF 左对齐 ✅ / --min-tree-coverage
   │
   ├─ §5.0.1 测试改进（借鉴 impg/seqwish）✅: 传递闭包不变量 / GFA round-trip / indel query 精度
   │
   ↓
V6 图质量 (graph-report ✅)
   ↓
V7 规模扩展 (4 万大肠杆菌基准 → 选 1-2 项 seqwish 优化) ← 待数据可用
   ↓
V8 应用层 (genotyping，远期)

* 增量增强：indel VCF / --max-depth / --syntenic-filter 已落地；
  --min-transitive-len 语义由 --min-len 覆盖，未单列选项（见 §3.1）
```

**§5.0.1 测试改进已全部完成**（9 个新测试，无需外部数据）：
传递闭包不变量 4 个 + GFA round-trip 3 个 + indel query 精度 2 个。
§5.0 打磨项中 `--min-tree-coverage` 需 phylogeny 上下文。V6 `graph-report` 已完成。V7 规模扩展
待 4 万大肠杆菌数据可用后启动。V8 视用户反馈与真实数据验证后再启动。

**测试文件结构重组**：原 `cli_paf_query.rs`（59 个测试，覆盖 query + to-maf/to-vcf/to-gfa/to-bed +
BFS + transitive 等多子命令）已按子命令拆分为 11 个独立文件，每个文件聚焦单一子命令或核心行为，
辅助函数（`write_bgzf_fa` / `spell_gfa_paths` / `revcomp`）按需在各文件内独立包含，无跨文件依赖：

- `cli_paf.rs`（12）— `paf index` 基础与持久化 roundtrip
- `cli_paf_bgzf.rs`（12）— BGZF/gzip 输入与懒加载 CIGAR
- `cli_paf_query.rs`（20）— `paf query` 核心参数与过滤
- `cli_paf_query_bfs.rs`（10）— BFS mirror 与 transitive 闭包不变量
- `cli_paf_graph.rs`（9）— `paf graph`（V4a 全局图）
- `cli_paf_graph_report.rs`（6）— `paf graph-report`（25 维度统计）
- `cli_paf_to_bed.rs`（4）— `paf to-bed` + indel 坐标精度
- `cli_paf_to_gfa.rs`（7）— `paf to-gfa`（V4b 局部图）+ path round-trip + lowercase round-trip
- `cli_paf_to_maf.rs`（13）— `paf to-maf`（pairwise + MSA）
- `cli_paf_to_vcf.rs`（6）— `paf to-vcf`（SNP + INS/DEL + 左对齐）
- `cli_paf_real.rs`（3）— 真实 multiz 数据回归

**跨平台临时文件**：所有 paf 测试已统一使用 `tempfile::TempDir` 替换硬编码 `/tmp/` 路径，
确保在 Windows 上可运行（`TempDir` 在 drop 时自动清理，删除了原 `fs::remove_file` 手动清理代码）。

---
## 6. 存量资产优势

通读 notes/ 下全部文档并分析 pgr 源码后，对 pgr 已有资产的认识持续深化。 以下发现显著降低了实现门槛。

### 6.1 `loc.rs` — pgr 的 IO 层比 impg 更成熟

分析了 `src/libs/loc.rs`（202 行）与 impg `paf.rs`（417 行）的对应关系。核心发现：

- **`Input` enum 比 impg 的 `PafHandle` 更强**：多了 `Buf` 变体（支持 stdin）， 且 `Bgzf` 变体使用
  `IndexedReader`（自带索引，seek 无需外部 `.gzi` 文件）
- **`read_offset()` 可直接替代 impg 的 `read_cigar_data()`**：同样是 seek+read+返回字节， pgr
  的实现更简洁（11 行 match + 2 行 I/O vs impg 的 46 行分支）
- **pgr 已有 BGZF 行读取能力**（`create_loc` 中对 `Input::Bgzf` 调用 `read_line`）， 只是需要抽象为
  `Input::read_line` 方法供 PAF 解析使用

**结论**：PAF 模块中最棘手的 IO 部分（多格式输入、BGZF 随机访问、CIGAR 懒加载） pgr 已经解决了 80%。
真正需要从零写的只有三样：区间树索引、PAF 行解析、CIGAR 编解码。

### 6.2 `IndexedReader` 自带索引能力，不需要 impg 的 GZI 机制

impg 的 `parse_paf_bgzf_with_gzi` 需要外部 `.gzi` 索引文件来做多线程解压，且需要
显式 `bgzf::VirtualPosition::from(offset)` 转换。pgr 的 `bgzf::io::IndexedReader`
在内部处理了这一切——调用者只需传字节偏移量。

这意味着 pgr 的 BGZF PAF 支持可以**跳过 impg 的模式 3**（GZI 两遍扫描）， 直接用 `IndexedReader`
做到同等性能。

### 6.3 pgr 已有的比对生成能力

pgr 有完整的 lastz 封装（7 套预设参数、并行执行），可以为特定 pair 生成 pairwise 比对。 `--self`
模式是重复屏蔽管道的一部分（碎片自比对），不是泛基因组比对工具。

### 6.4 Cactus Caf 的"退火-熔化"循环对 pgr 挑选机制的直接参考

`cactus.md` §8 详细分析了 Caf 模块（`caf.c`、`annealing.c`、`melting.c`）的迭代循环：

- **Annealing（加法）**：把两两比对捏合成 Pinch Graph 中的 Block。关键约束是
  `stCaf_annealBetweenAdjacencyComponents`——"只连接不同连通分量的序列对"，避免在
  同一连通区域形成复杂环
- **Melting（减法）**：按 Degree、Tree Coverage、Chain Length 进行多维过滤，
  `stCaf_getBlocksInChainsLessThanGivenLength` 丢弃破碎短链

对应的过滤维度可以搬到 pgr 的查询层：

| Caf 过滤维度   | pgr 对应参数           | 语义                       | 状态         |
|----------------|------------------------|----------------------------|--------------|
| Degree         | `--min-degree N`       | 过滤支持序列数 < N 的区间  | ✅ V1 已实现 |
| Tree Coverage  | `--min-tree-coverage`  | 过滤进化树上分布稀疏的区间 | 待实现       |
| Chain Length   | `--min-chain-length N` | 过滤总长 < N bp 的传递链   | ✅ V1 已实现 |
| Block End Trim | `--end-trim N`         | 切除比对边缘不可靠的 N bp  | ⏭️ 推迟      |

**但要警惕**：Caf 的 melting 在**图构建时**做（离线、全局视角），而 pgr 的挑选在**查询时**做
（在线、局部视角）。查询时无法做全图 Tree Coverage 计算。因此这些 Caf 过滤维度更适合作为传递闭包的
**后处理过滤**，而非 BFS 本身的中断条件。`--end-trim` 推迟——它需要 per-interval 修剪 CIGAR 两端，
与当前"区间整体投影"的输出模型不兼容，待 V2 引入序列输出时 一并处理。`--min-tree-coverage`
需要进化树上下文，留待后续阶段。

### 6.5 Minigraph-Cactus 分治策略对 pgr partition 的启示

`cactus.md` §3.1 详述了 Minigraph-Cactus 的五阶段流程：

```
Minigraph 骨架构建 → 图映射定位 → rgfa-split 切分 → 批量 Cactus 比对 → 合并
```

与 impg 的 Partition + Lace 模式（[[impg.md]] §6.3）对比：

| 维度     | Minigraph-Cactus               | impg Partition + Lace         |
|----------|--------------------------------|-------------------------------|
| 拆分依据 | Minigraph SV 图连通分量        | 传递闭包 + masking 去重       |
| 拆分粒度 | 染色体级（MB）                 | locus 级（KB-MB，窗口可控）   |
| 局部处理 | Cactus full pipeline (Caf+Bar) | 独立 GFA 构建 (seqwish/crush) |
| 合并方式 | HAL/VG join                    | lace（坐标驱动重新拼装）      |

**对 pgr 的启示**：如果 pgr 未来需要 partition（处理 > 100 基因组的 cohort）， 建议走
**Minigraph-Cactus 的"先粗后细"路线**：

- **粗拆分**：用已有的 Chain/Net syntenic 信息做染色体级拆分（类似 `cactus_graphmap_split.py` 的
  heuristic contig selection：regex + size + dropoff，见 `cactus.md` §2.4.2）
- **细拆分**：在每个大区块内，用传递闭包 BFS + masking 去重切分成 per-locus 批次
- **比对**：per-locus 跑 `pgr lav lastz`（不含 `--self`）生成局部 pairwise
- **合并**：per-locus PAF 汇总回全 cohort 区间树

这个路线复用了 pgr 的三项已有资产：Chain/Net syntenic 信息、`pgr lav lastz`、PAF 区间树。 但这是
**第二步或第三步的任务**，第一步不需要 partition。

### 6.6 pgr 已有的 MSA 资产（供后续阶段按需使用）

以下组件不在查询层使用，但在方向 D（图构建）或下游分析中可以直接调用：

| 组件               | 源码                          | 后续用途                                           |
|--------------------|-------------------------------|----------------------------------------------------|
| POA 引擎           | `libs/poa/poa.rs`             | 图构建阶段的 per-bubble 共识/比对                  |
| Banded DP          | `libs/fas_multiz.rs`          | partition 内多 pairwise 合并（比 impg POA 更精确） |
| `get_subs`         | `libs/alignment.rs:214`       | MSA 上的变体检测                                   |
| 裁剪函数           | `libs/alignment.rs:1351-1687` | BFS 结果边界清理                                   |
| crossbeam 并行管道 | `consensus.rs:250`            | `build_multi` 并行化                               |

但这些都是**独立的 CLI 命令或库函数**，通过 Unix pipe 组合，不与 `paf query` 耦合。

---
## 7. 暂不实现（明确边界）

V1-V5 图构建层已完成。以下功能仍明确排除，每条给了触发条件以防止 scope creep。

| 暂不实现                                | 理由                                       | 重新评估的触发条件                                                    |
|-----------------------------------------|--------------------------------------------|-----------------------------------------------------------------------|
| 补充 pairwise 比对（第二层）            | 已有 MAF/PAF 复用已足够                    | 传递闭包覆盖率不足                                                    |
| Mash KNN pair-selection                 | 小 cohort 有 MAF 先验时不需要              | 大 cohort 无先验（如 4 万 E. coli，见 [[ecoli-cohort.md]]）           |
| `pgr lav lastz --self` 全自比对         | 此 flag 用于重复屏蔽管道，非泛基因组比对   | 需要全新 cohort 的 pairwise 比对时评估 `pgr lav lastz`（不含 --self） |
| syng 免比对后端                         | [[impg.md]] §1.1 已明确不参考              | 永不                                                                  |
| partition / lace / refine               | 处理 >100 基因组的 cohort 时需要           | N > 50                                                                |
| stage DSL                               | 单命令不需要管道化                         | 出现三个以上 stage 串联                                               |
| 基因分型（genotype/infer）              | 能力栈顶端，依赖图构建层（[[impg.md]] §7） | 图构建层已就绪（V4-V5 ✅），需真实数据验证                            |
| 完整 crush 8 阶段 + 15 ResolutionMethod | V4b `--crush` 已覆盖主要场景               | 复杂 bubble 处理需求出现                                              |
| 全局精细 GFA 合并                       | 见 §2.5 决策                               | 永不                                                                  |
| minigraph `gfa_t` 对应实现              | pgr 走 seqwish DSU 路线                    | 永不                                                                  |
| 1ALN/TPA 格式支持                       | pgr 用 PAF/MAF                             | 永不                                                                  |

