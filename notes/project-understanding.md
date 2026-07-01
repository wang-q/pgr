# pgr 项目理解

本文档是我对 pgr (Practical Genome Refiner) 项目的整体理解，涵盖架构、设计哲学、代码模式、
当前能力与未来方向。写作时间：2026-06-27。

## 1. 项目定位

pgr 是一个**生物信息学 CLI 工具集**，定位是"基因组数据处理瑞士军刀"。它不追求成为某个领域的
一站式平台，而是在日常生信流程中充当格式转换、数据查询、快速分析的胶水层。

当前版本 0.2.0，作者 wang-q，MIT 协议，Rust 2021 edition。

### 1.1 与同类工具的差异

| 工具      | 定位                        | pgr 的区别                               |
|-----------|-----------------------------|------------------------------------------|
| UCSC kent | C 工具集，UCSC 格式权威实现 | Rust 实现，更安全的错误处理，无 panic    |
| samtools  | SAM/BAM/CRAM 专用           | pgr 不做 SAM/BAM，focus on UCSC 链       |
| bedtools  | BED/VCF/GFF 区间操作        | pgr 有 2bit/Chain/Net/Newick 等独特覆盖  |
| seqtk     | FASTA/FASTQ 轻量处理        | pgr 的 fa/fas 子命令更多更全             |
| pggb/odgi | 泛基因组图构建与分析        | pgr 目前是 pairwise 工具，泛基因组在路上 |

### 1.2 核心优势

- **格式覆盖广**：从序列 (FASTA/FASTQ/2bit) 到比对 (AXT/PSL/Chain/Net/MAF/LAV) 再到 系统发育 (Newick)
  和聚类，覆盖了生信日常的大多数格式。
- **UCSC Chain/Net 体系深耕**：`chain`/`net`/`axt`/`psl`/`lav`/`maf` 六个模块形成了 一套完整的
  pairwise 比对处理链，从 lastz 原始输出到 net 再到 multiz MAF。
- **Rust 实现的健壮性**：零 panic 策略，畸形输入返回友好错误而非崩溃。
- **Pipeline-friendly**：stdin/stdout 支持，可组合子命令。
- **已有 pairwise 比对基础设施成熟**：这是 pgr 最大的存量资产，也是走向泛基因组时 区别于 impg/pggb
  的独特起点。

## 2. 架构全景

### 2.1 双层结构

```
src/
├── pgr.rs              # 入口：clap 命令树定义 + dispatch
├── lib.rs              # 库入口：仅 re-export io 工具
├── cmd_pgr/            # 命令层：每个模块 = 一组子命令
│   ├── mod.rs          #   声明所有命令模块
│   ├── fa/             #   每个模块含 make_subcommand() + execute()
│   ├── fas/
│   ├── fq/
│   ├── twobit/         #   命令名 2bit
│   ├── gff/
│   ├── chain/
│   ├── net/
│   ├── axt/
│   ├── psl/
│   ├── lav/
│   ├── maf/
│   ├── paf/            #   PAF 泛基因组图操作
│   ├── clust/
│   ├── dist/
│   ├── mat/
│   ├── nwk/
│   ├── ms/
│   ├── pl/             #   Pipelines：编排外部工具
│   └── plot/           #   可视化输出 (TikZ/LaTeX)
└── libs/               # 核心库层：数据结构与算法
    ├── mod.rs
    ├── phylo/          #   系统发育树核心 (Arena 树结构)
    ├── poa/            #   偏序比对 (Partial Order Alignment)
    ├── chain/          #   Chain 算法 (连接、gap 计算、替换矩阵)
    │   └── net/        #     Net 格式处理 (class/filter/to-axt)
    ├── clust/          #   聚类算法 (hier/DBSCAN/MCL/k-medoids/NJ/UPGMA)
    │   ├── eval/       #     聚类评估
    │   └── tree_cut/   #     树切分方法
    ├── ms/             #   Hudson's ms 模拟器解析
    ├── fmt/            #   格式解析 (AXT/FAS/FA/FQ/LAV/MAF/PSL/2bit)
    ├── paf/            #   PAF 隐式图核心
    │   ├── index/      #     区间树索引 + BFS 传递闭包
    │   └── graph/      #     DSU 图构建 + GFA 输出
    └── ...             #   io, hash, hv, linalg, loc, nt, pairmat
```

**关键设计**：`cmd_pgr/` 管命令分发，`libs/` 管核心逻辑。命令层薄，逻辑在库层。这是好的分层。

### 2.2 命令分发模式

每增加一个子命令（如 `pgr fa size`），遵循固定模式：

1. 在 `src/cmd_pgr/<category>/<name>.rs` 实现 `make_subcommand()` 和 `execute()`
2. 在 `<category>/mod.rs` 中声明模块、在 `make_subcommand()` 中注册、在 `execute()` 中 dispatch
3. 在 `src/pgr.rs` 中注册顶层子命令并 dispatch

**两跳 dispatch**：`pgr.rs` (第一跳: `fa`) → `cmd_pgr/fa/mod.rs` (第二跳: `size`) →
`cmd_pgr/fa/size.rs` (实际执行)。这种模式在 clap 中常见，优点是每个叶子命令独立文件、易于维护；
缺点是新增命令需要改三处。

### 2.3 依赖策略

从 `Cargo.toml` 看，依赖选择明显偏好**成熟、稳定、轻量**的 crate：

- **CLI**: `clap` 4.x (非 derive 模式，手写 Command 树)
- **解析**: `nom` 8.x (Newick parser)、`regex`
- **并发**: `rayon` (数据并行)、`crossbeam`
- **生物信息学**: `noodles` (FASTA/FASTQ/BGZF/GFF)、`bio` (生物信息学基础类型)
- **数据结构**: `petgraph` (图)、`indexmap` (保序 HashMap)、`intspan` (区间集合)
- **哈希**: 多哈希支持 — `rapidhash`、`fxhash`、`murmurhash3`、`xxhash-rust`
- **输出**: `rust_xlsxwriter` (Excel)、`tera` (模板引擎)
- **外部工具编排**: `cmd_lib`、`which`

**没有引入的**：`tokio`（纯同步 CLI 工具，不需要异步运行时）、`sled`/`rocksdb`（不做嵌入式数据库）。
**已引入的**：`serde`（带 derive）+ `bincode`（用于 PAF 索引的 `.paf.idx` 持久化）。

### 2.4 构建配置

- `#![feature(portable_simd)]` — 使用了 nightly 的 portable SIMD
- `lto = true` — release 构建启用链接时优化
- 测试框架：`assert_cmd` + `predicates` 做 CLI 集成测试，`criterion` 做基准测试

## 3. 命令模块全景

按功能域分类：

### 3.1 序列 (Sequences)

| 模块     | 子命令数 | 核心能力                                              |
|----------|----------|-------------------------------------------------------|
| `fa`     | 18       | FASTA 全能操作：统计、筛选、切分、转换、mask、索引    |
| `fas`    | 20       | Block FA (多序列比对块)：统计、筛选、subset、变异检测 |
| `fq`     | 2        | FASTQ 交叉合并、转 FASTA                              |
| `twobit` | 5        | 2bit 二进制格式查询：range、sequence、masked 统计     |
| `gff`    | 1        | GFF 注释：rg (read group)                             |

**fa 和 fas 是序列模块的核心**，子命令最多、功能最全。`fas` 的 `multiz`、`variation`、 `refine`、
`to_vcf` 已经触及多序列比对和变异检测。

### 3.2 基因组比对 (Alignments)

| 模块    | 子命令数 | 核心能力                                             |
|---------|----------|------------------------------------------------------|
| `chain` | 6        | Chain 排序、过滤、split、stitch、反重复、转 net 准备 |
| `net`   | 6        | Net 分类、过滤、split、subset、syntenic、转 AXT      |
| `axt`   | 4        | AXT 排序、转 FAS/MAF/PSL                             |
| `psl`   | 8        | PSL 统计、直方图、lift、swap、转 chain、转 range     |
| `lav`   | 2        | LAV (lastz 原生输出) 转 PSL、lastz 调用封装          |
| `maf`   | 2        | MAF (multiple alignment format) 转 Block FA、转 PAF  |

**这是 pgr 最成熟的模块群**。完整覆盖了 UCSC 的 lastz → axtChain → chainAntiRepeat →
chainMergeSort → chainPreNet → chainNet → netSyntenic → netChainSubset → netToAxt →
axtToMaf 标准化流程中的大部分步骤。`chain`/`net` 模块在功能上可以替代 kent-tools 的 核心步骤
（虽然部分高级功能仍依赖外部工具）。

### 3.3 泛基因组 (Pangenome)

| 模块  | 子命令数 | 核心能力                                                                 |
|-------|----------|--------------------------------------------------------------------------|
| `paf` | 9        | PAF 隐式图：索引、查询、to-bed、to-fas、to-maf、graph、to-gfa、to-vcf、stat |

`paf` 模块是 pgr 走向泛基因组的核心载体。基于 PAF (Pairwise mApping Format) 的 all-vs-all
比对，构建隐式泛基因组图：

- **索引层**：`pgr paf index` 把 PAF 全量装入区间树，支持 `.paf.idx` 持久化
- **查询层**：`query` / `to-bed` / `to-maf` 按需投影目标区间，BFS 传递闭包找全同源片段
- **图构建层**：`graph` 粗全局 GFA（seqwish DSU 风格，零序列依赖拓扑模式）；`to-gfa` 区域精细 GFA；
  `to-vcf` POA MSA 导出变异；`stat` 图拓扑统计报告

详见 [[paf-pangenome.md]]。

### 3.4 聚类 (Clustering)

| 子命令     | 算法                            |
|------------|---------------------------------|
| `cc`       | 连通分量 (Connected Components) |
| `dbscan`   | DBSCAN 密度聚类                 |
| `hier`     | 层次聚类 (NN-chain 算法)        |
| `kmedoids` | K-medoids 划分聚类              |
| `mcl`      | Markov Cluster Algorithm        |
| `nj`       | Neighbor-Joining 建树           |
| `upgma`    | UPGMA 建树                      |

聚类模块算法覆盖广：基于密度的 (DBSCAN)、划分的 (k-medoids)、层次的 (hier/NJ/UPGMA)、 图的 (MCL)、
简单的 (CC)。评估子命令 (`eval`) 在设计阶段。

### 3.5 距离与矩阵 (Distance & Matrix)

| 模块   | 子命令数 | 核心能力                                    |
|--------|----------|---------------------------------------------|
| `dist` | 3        | 距离计算：hv (hypervariable)、seq、vector   |
| `mat`  | 6        | 矩阵操作：compare、format、subset、转换格式 |

`mat` 充当聚类流程的"数据预处理"环节：`mat to-pair` 把矩阵转成 pair 格式供聚类使用， `mat to-phylip`
转 PHYLIP 格式供外部工具。

### 3.6 系统发育 (Phylogeny)

| 模块  | 子命令数 | 核心能力                                            |
|-------|----------|-----------------------------------------------------|
| `nwk` | 18       | Newick 树全能操作：统计、比较、剪枝、reroot、可视化 |

`nwk` 模块功能非常丰富：树拓扑比较 (`cmp`、`topo`)、切分 (`prune`、`subtree`、`to-forest`)、
重标 (`rename`、`label`、`comment`)、重根 (`reroot`)、可视化 (`to-svg`、`to-dot`、`to-tex`)、
统计 (`stat`、`distance`、`support`)。底层 `libs/phylo/` 使用 Arena 树结构（非 Rc/RefCell 指针树），
参考了 `notes/design/phylo.md` 的设计讨论。

### 3.7 模拟、流程、可视化 (Simulation, Pipelines, Plot)

| 模块   | 子命令数 | 核心能力                                                |
|--------|----------|---------------------------------------------------------|
| `ms`   | 1        | Hudson's ms 模拟器输出转 DNA 序列                       |
| `pl`   | 7        | 集成流程：p2m、trf、ir、rept、ucsc、prefilter、condense |
| `plot` | 3        | TikZ/LaTeX 图：Venn、HH (hedgehog)、NRPS                |

`pl` (pipelines) 模块定位特殊——它**编排外部工具**（UCSC kent-tools、trf、FastK、Profex、
clustalw/muscle/mafft），充当工作流 glue。这与 `chain`/`net` 模块的纯 Rust 实现形成互补：能用 Rust
就自己实现，复杂/成熟的用外部工具。

## 4. 核心库层详解

### 4.1 `libs/phylo/` — 系统发育核心

- **Arena 树结构** (`node.rs`)：所有节点存储在 `Arena` 中，通过 `NodeId` 索引引用。 避免了 Rust
  中树结构的自引用问题（不用 `Rc<RefCell<>>`)。
- **Newick 解析** (`parser.rs`)：用 `nom` 手写解析器。
- **树比较** (`cmp.rs`)：Robinson-Foulds 距离等。
- **树算法** (`algo.rs`)：排序、reroot 等操作。

### 4.2 `libs/poa/` — 偏序比对

- 实现 Partial Order Alignment 算法（参考 SPOA）
- `graph.rs`：POA 图结构
- `align.rs`：序列到图的比对
- `consensus.rs`：从 POA 图提取一致性序列
- `msa.rs`：多序列比对接口

### 4.3 `libs/chain/` — Chain 核心逻辑

- `record.rs`：Chain 记录定义与解析
- `connect.rs`：Chain 连接（核心 chaining 算法）
- `gap_calc.rs`：Gap 计算（线性和仿射罚分）
- `sub_matrix.rs`：DNA 替换矩阵（如 HoxD55）
- `kdtree.rs`：高效前驱搜索的数据结构
- `anti_repeat.rs`：反重复处理
- `net/`：Net 格式处理子模块（builder/class/filter/finalize/reader/subset/syntenic/to-axt/types/writer）

### 4.4 `libs/clust/` — 聚类算法库

- `hier.rs`：NN-chain 层次聚类实现（参看 `docs/clust-hier.md`）
- `dbscan.rs`、`mcl.rs`、`k_medoids.rs`：各算法实现
- `nj.rs`、`upgma.rs`：建树算法
- `medoid.rs`：medoid 计算
- `feature.rs`：特征提取
- `format.rs`：格式处理
- `eval/`：聚类评估子模块（coordinates/distance/pairwise/partition）
- `tree_cut/`：树切分方法（clade/dynamic/hybrid/inconsistent/link/simple）

### 4.5 `libs/fmt/` — 格式解析

统一管理生物信息学格式的解析逻辑：

- `axt.rs`、`fas.rs`、`fa.rs`、`fq.rs`、`lav.rs`、`maf.rs`、`psl.rs`、`twobit.rs`、`vcf.rs`

最近重构过：原先在 `libs/` 根目录下的 `axt.rs`、`fas.rs`、`lav.rs`、`maf.rs` 等移入 `libs/fmt/`。

> **MAF 支持现状**（2026-07 确认）：

> - `maf.rs`（363 行）：完整的读写支持
>   - 读取：`MafComp`（s 行结构体）、`MafAli`（a 行 + components，含 `score=` 解析）、`next_maf_block()`（流式读取）、`parse_maf_block()`
>   - 写入：`MafWriter`（header + block 输出）
>   - 坐标转换：`MafComp::to_range()`（0-based → 1-based inclusive，含负链处理）
> - `cmd_pgr/maf/`：`to-fas` 和 `to-paf` 两个子命令

### 4.6 其他库

- `libs/io.rs`：I/O 辅助（`read_lines`、`reader`、`writer`）
- `libs/hash.rs`：哈希工具
- `libs/linalg.rs`：线性代数
- **`libs/loc.rs`**：FASTA 随机访问索引模块。`Input` enum（Buf/File/Bgzf）+
  `create_loc`（建 `.loc` 索引）+`read_offset`（seek+read）。
  **2026-06 发现：此 IO 抽象层可直接支撑 PAF 模块的 CIGAR 懒加载和 BGZF 随机访问，比 impg 的 `paf.rs` IO 层更成熟**。
  见 [[paf-pangenome.md]] §6.1。
- `libs/nt.rs`：核苷酸类型
- `libs/pairmat/`：pair 矩阵
- `libs/hv.rs`：hypervariable 区域
- `libs/chain/net/`：Net 格式处理（已移入 chain 子模块）
- `libs/fmt/twobit.rs`：2bit 格式读写
- `libs/fmt/psl.rs`：PSL 格式
- `libs/alignment/`：比对通用逻辑
- `libs/fas_multiz/`：Multiz 多序列比对处理（banded DP 合并）
- `libs/fas_xlsx.rs`：FASTQ 到 Excel 转换
- `libs/fasta/`：FASTA 处理工具（dedup/filter/stat）
- `libs/paf/`：PAF 隐式图核心（索引、查询、图构建、VCF 导出）
- `libs/ms/`：Hudson's ms 模拟器（解析器 + DNA 生成）
- `libs/plot/`：绘图工具（histogram/nrps/venn）
- `libs/lastz.rs`：lastz 调用封装
- `libs/par.rs`：并行辅助
- `libs/translate.rs`：序列翻译（六框翻译）

## 5. 设计模式与约定

### 5.1 命令模式

每个叶子命令文件遵循统一结构：

```rust
// make_subcommand: 定义 CLI 接口
pub fn make_subcommand() -> Command {
    Command::new("size")
        .about("Counts total bases and sequences in FASTA files")
        // .after_help(...)  可选
        .arg(arg!(-i --infile <FILE> "Input FASTA file"))
        .arg(arg!(-o --outfile <FILE> "Output file"))
}

// execute: 执行逻辑
pub fn execute(matches: &ArgMatches) -> anyhow::Result<()> {
    // 1. 提取参数
    // 2. 打开输入/输出
    // 3. 处理数据
    // 4. 写入结果
    Ok(())
}
```

- 用 `anyhow::Result<()>` 做返回值（CLAUDE.md 硬性要求）
- 不使用 clap derive 宏，手写 `Command` 构建
- `about` 用第三人称单数
- 输入参数统一命名：`infile`（单文件）或 `infiles`（多文件）

### 5.2 零 Panic 策略

所有用户输入（畸形数据、二进制文件、错误的命令行参数）必须返回友好错误，不能 panic。 这贯穿了 `nom`
解析器（返回 `Result`）、文件 I/O（`anyhow::Context` 附加错误上下文）、索引越界检查等所有层面。

### 5.3 测试约定

- 集成测试在 `tests/cli_<command>.rs`
- 测试数据在 `tests/<command>/`
- 使用 `PgrCmd` 辅助结构体（`tests/common/mod.rs`）
- 测试函数**不需要**返回 `anyhow::Result<()>`
- 必须用 `assert_cmd` 定位二进制

### 5.4 帮助文本规范

- `about`: 第三人称单数
- `after_help`: 用 raw string `r###"..."###`
    - Description → Notes (bullet `*`) → Examples (numbered `1.`)
- 标准 note (fa/fas): `* Supports both plain text and gzipped (.gz) files`
- 标准 note (twobit): `* 2bit files are binary and require random access (seeking)`

## 6. 项目现状评估

### 6.1 已完成的（成熟）

- **pairwise 比对全链路**：lastz → chain → net → axt → maf 的工具链在 `chain`/`net`/ `axt`/`psl`/
  `lav`/`maf` 六个模块中基本完整，且是纯 Rust 实现（不依赖 kent-tools）。
- **FASTA/FASTQ/2bit 处理**：`fa`(18 子命令) + `fas`(20 子命令) + `fq`(2) + `twobit`(5)，
  日常序列操作需求基本覆盖。
- **系统发育树操作**：`nwk`(18 子命令) 功能丰富，可视化 (SVG/DOT/TikZ) 也已有。
- **聚类算法**：7 种算法已实现。
- **距离/矩阵工具链**：`dist` → `mat` → `clust` 的数据流完整。

### 6.2 进行中的（活跃开发）

- **泛基因组方向**：已形成完整路线图，整合在 `notes/paf-pangenome.md`（路线决策 + 已实现能力 +
  代码结构 + 后续规划）。query-to-vcf 已全部完成。
  **2026-06 发现：`libs/loc.rs` 的 IO 抽象层可直接支撑 PAF 的 CIGAR 懒加载和 BGZF 访问， 实际实现量比最初估计少约 30%**。

- **`pl` 流程模块**：`ucsc`、`trf`、`rept`、`ir` 等 pipeline 在补充。

### 6.3 待补全的（TODO / 设计阶段）

- `pgr.rs` 末尾注释的 TODO：paralog、fas variation --indel、fas match、去完全包含序列
- `docs/clust-eval.md`：聚类评估（设计中）
- `notes/design/nwk-eval.md`：树结构多维度评估（设计中）

PAF 泛基因组方向（query-to-vcf）已全部完成，后续规划见 [[paf-pangenome.md]] §5（stat 规模扩展 / V7 图质量 /
V8 应用层）。

### 6.4 不做 / 不适合做的

- **SAM/BAM/CRAM 处理**：留给 samtools/pysam，pgr 聚焦 UCSC 链
- **Variant calling**：不是 pgr 的领域
- **Web 服务 / GUI**：pgr 是纯 CLI 工具
- **异步 I/O**：不需要 tokio，rayon 数据并行已足够
- **JSON/YAML 序列化**：目前不需要，输出以 TSV/FASTA/Newick 为主

## 7. 与周边项目的关系

### 7.1 UCSC kent-tools 的关系

pgr 的 `chain`/`net`/`axt`/`psl` 模块是 UCSC kent-tools 对应功能的**Rust 重实现**， 但并非完全替代：

- **pgr 自己实现的**：chain sort/split/stitch、net filter/split/subset/class、axt sort/to_fas、 psl
  stats/lift/swap
- **仍依赖 kent-tools 的**（通过 `pgr pl ucsc` 编排）：chainAntiRepeat、chainMergeSort、
  chainPreNet、chainNet、netSyntenic、netChainSubset 等复杂步骤

### 7.2 impg 的关系

参见 `notes/impg.md` 的详细分析。核心差异：

- impg 的"隐式图"用的是 PAF/wfmash 生态，pgr 用的是 Chain/Net 生态
- impg 需要 all-vs-all 比对（它不省比对，省的是图构建）
- pgr 已有 pairwise 比对基础设施，走向泛基因组时可以**复用已有比对而非重跑**

### 7.3 Cactus 的关系

参见 `notes/references/cactus.md` 和 `notes/references/cactus_lastz.md`。pgr 的 `pgr lav lastz` 子命令封装了 lastz 调用，
采用了 Cactus 风格的参数。Cactus 的 transitive alignment 机制对 pgr 泛基因组 方向有参考价值。

## 8. 关键风险与技术债

1. **nightly 依赖**：`#![feature(portable_simd)]` 是 nightly-only，锁死了编译器版本。 如果
   portable_simd 迟迟不稳定，可能成为包袱。
2. **maf 模块扩展**：`maf` 目前有 `to-fas` 和 `to-paf` 两个子命令。仍需 `filter`、`subset` 等扩展以支撑泛基因组管道。
3. **命令树深度嵌套**：三跳 dispatch (`pgr.rs` → `mod.rs` → `leaf.rs`) 在新增命令时容易遗漏
   某一层的注册。可以考虑宏简化。
4. **测试覆盖不均衡**：`fa`、`nwk` 等大模块有较全的集成测试，但 `chain`/`net`/`pl` 的覆盖 可能不足
   （`pl` 依赖外部工具，测试困难）。
5. **外部工具依赖的文档化**：`pl ucsc` 需要安装一整套 kent-tools，但错误提示不够友好。
6. **`fas` 模块职责过重**：20 个子命令塞在一个模块下，`fas multiz` 等复杂逻辑可能需要
   拆分为独立顶层命令。

