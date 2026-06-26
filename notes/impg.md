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

### 1.1 设计哲学: 隐式泛基因组图

传统泛基因组工具（如 pggb、Minigraph-Cactus）需要先物化一张 GFA 图， 再在图上做下游分析。`impg`
走第三条路：

- **比对即图**：把 PAF/1ALN/TPA 两两比对视为图的隐式描述 — 比对双方的坐标区间就是"边"，
  序列本身是"节点"。
- **按需投影**：给定一个目标区间 `seq:start-end`，在区间树上查找所有重叠的比对， 把目标坐标 lift
  到查询序列上，输出 BED/BEDPE/PAF/FASTA/GFA/VCF/MAF。
- **传递闭包**：`-x` 选项递归地把初次结果当作新的查询目标，找全所有同源片段， 类似 Cactus 的
  transitive alignment 但延迟到查询时才执行。
- **不物化图**：除非用户显式要求 `gfa`/`vcf` 输出，否则永远不会落到 GFA 文件。

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
`.rs` 文件、约 8.4 万行；其中 `main.rs` 1.57 万行、`resolution.rs` 1.72 万行、`syng.rs` 0.9 万行
三个超大文件占了近一半。`lib.rs` 导出 31 个模块，按职责分六组：

**核心索引（alignment 后端投影核心，4 文件）**

- `impg.rs` (3214 行) — 单文件 `.impg` 索引：`Impg`/`QueryMetadata`/`CigarOp`/`SortedRanges`，
  thread-local WFA aligner 缓存，CIGAR 懒加载。
- `impg_index.rs` (392 行) — `ImpgIndex` trait，统一单/多文件索引接口（§3.4）。
- `multi_impg.rs` (1093 行) — `MultiImpg`：协调多个 per-file `.impg` 子索引，`TreeLocation`
  全局↔本地 ID 翻译，staleness 检测。
- `forest_map.rs` (33 行) — `ForestMap`：`target_id → tree_offset` 反向索引，序列化用。

**比对格式统一抽象（5 文件）**

- `alignment_record.rs` (138 行) — `AlignmentRecord`：PAF/1ALN/TPA 三格式的统一 8 字段表示，
  `strand_and_data_offset` 高位存 strand、低位存文件偏移。
- `paf.rs` (416 行) — PAF 解析，支持 BGZF + GZI 索引。
- `onealn.rs` (894 行) — `.1aln` 解析（`onecode` 库），含 tracepoint + O(1) seek 元数据。
- `tpa_parser.rs` (236 行) — `.tpa` 解析（`tpa` 库），多线程 BGZF 解压。
- `seqidx.rs` (56 行) — 轻量 `SequenceIndex`：`name↔id↔len` 三映射，无 IO。

**GFA 图构建与 crush 算法（8 文件，泛基因组核心）**

- `graph.rs` (1207 行) — `SequenceMetadata`/`path_name()`（PanSN 坐标化）、`unchop_gfa`/`sort_gfa`
  （gfasort 集成）、`TerminalNRunClip`、SPOA 引擎构建。
- `resolution.rs` (17169 行) — **crush 算法实体**（§3.5），bubble-guided graph resolution，
  `ResolutionConfig`/`ResolutionMethod`/`MultiLevelWindowMode`/`ResolutionPolishMethod` 等
  大量枚举，`resolve_gfa_bubbles` 入口，POA/POASTA/abPOA 路由。
- `graph_pipeline.rs` (178 行) — `GraphPipelineSpec` 解析器（`stage,key=value:stage,...` DSL）。
- `graph_report.rs` (2688 行) — `GraphReport`/`GraphReportOptions`：white_space/sparse_coverage/
  repeat_context/path_jump/link_jump/topology 等 15+ 维度图质量报告，被所有 sweep 脚本消费。
- `smooth.rs` (2852 行) — smoothxg 风格块分解 + SPOA 平滑，`SmoothBlockSource` 三策略
  （PathOverlap/Flubble/NeighborMergePoasta）。
- `gfa_self_loops.rs` (1132 行) — blunt GFA self-loop 折叠 + `NormalizeSelfLoopsStats` 报告。
- `render_bundle.rs` (566 行) — `RenderManifest` + `RenderTranslationTables`：render 命令的
  可序列化图+坐标翻译包。
- `commands/syng2gfa.rs` (4652 行) — syng → GFA 物化，频率感知 syncmer 节点共享策略
  （`DEFAULT_GFA_MASK_TOP_FRACTION`/`DEFAULT_GFA_HIGH_FREQ_MIN_RUN` 等常量驱动）。

**命令实现（`commands/` 10 文件 + mod.rs）**

- `commands/mod.rs` (300 行) — `create_aligner`/`create_aligner_adaptive` 工厂（wfmash/fastga 后端，
  含 adaptive segment/sparsify/num_mappings 参数）。
- `commands/graph.rs` (1565 行) — `GraphBuildConfig` + `build_graph`/`induce_graph_from_alignment` /
  `run_graph_build`/`run_graph_build_poa`/`run_graph_build_pggb`/`run_graph_build_partitioned`。
- `commands/align.rs` (1321 行) — `AlignConfig` + sweepga sparsification 策略，PAF 输出。
- `commands/partition.rs` (1789 行) — cohort 窗口化切分，含 `rehome_singleton_slivers` （sliver
  重分配到 flank 邻居）。
- `commands/refine.rs` (948 行) — `RefineConfig`：收紧 loci 边界，支持 transitive BFS/DFS、
  blacklist、subset_filter。
- `commands/lace.rs` (2631 行) — 多格式（gzip/bgzf/zstd）GFA/VCF 拼接，路径名 `NAME:START-END`
  驱动。
- `commands/similarity.rs` (1053 行) — Jaccard/cosine/dice/estimated_identity 度量 + PCA/MDS。
- `commands/render.rs` (609 行) — `RenderConfig` + render bundle 输出。
- `commands/genotype.rs` (2886 行) — `SyngCosigtConfig`、`CandidateMode`（Spanning/Overlapping）、
  `parse_normalized_gfa`、`GraphContributionModel`。
- `commands/infer.rs` (1523 行) — `InferTarget`/`PartitionDiscoveryConfig`：跨区间/分区等位基因
  call + stitching。

**基因分型与投影（4 文件）**

- `genotyping.rs` (494 行) — 后端中立词汇：`FeatureSpace`（6 变体）/`EvidenceBackend`（6 变体） /
  `ScoringMethod`，为多后端扩展预留。
- `pack.rs` (344 行) — pack 二进制/TSV 格式（`IMPGPKB1` magic，1MB 块），`Coverage` 统计。
- `projection.rs` (114 行) — `ProjectionManifest` + `load`/`write_manifest`（projection bundle
  入口，format=`impg-projection` v1）。
- `projection/converter.rs` (646 行) — GAF→GFA projection 转换，
  `GfaProjectionOutputFormat`（ProjectionBundle/PackTsv）。

**序列索引与命名（5 文件）**

- `faidx.rs` (198 行) — `rust_htslib::faidx` 封装，LRU reader 缓存。
- `agc_index.rs` (264 行) — `ragc_core` AGC 归档支持，per-thread decompressor 池。
- `sequence_index.rs` (113 行) — `SequenceIndex` trait + `UnifiedSequenceIndex` 枚举 （Fasta 或
  Agc），按扩展名自动分发。
- `sequence_namespace.rs` (168 行) — `PanSn`/`PathIdentity`/`SourceSequenceRecord`/
  `SequenceNamespace`，`sample#haplotype#contig` 解析。
- `subset_filter.rs` (201 行) — `SubsetFilter`：exact/normalized/sample/sample+hap 四级匹配。

**syng 后端（6 文件，本文档不参考）**

- `syng.rs` (9072 行)、`syng_ffi.rs` (354 行)、`syng_graph.rs` (1606 行)、 `syng_graph_norm.rs` (307
  行)、`syng_parallel.rs` (184 行)、`syng_transitive.rs` (1956 行)。

### 1.4 关键依赖

- **`coitrees`** — cache-oblivious interval trees，区间查找的核心数据结构。
- **`gfa` + `handlegraph`** — GFA 读写与 handlegraph 抽象（用于 lace、crush 等图操作）。
- **`lib_wfa2`** — WFA 仿射比对，用于 BiWFA 边界精修。
- **`spoa_rs`、`poasta`** — POA MSA 引擎，用于 `gfa:poa` 与 crush 阶段。
- **`sweepga`、`seqwish`、`allwave`、`gfasort`、`bluntg`、`povu`** — 同一团队维护的泛基因组工具链，
  作为库直接嵌入（而非子进程调用）。
- **`ragc-core`** — AGC 序列归档格式支持。
- **`onecode`、`tpa`、`tracepoints`** — 1ALN / TPA 比对格式与 tracepoint 编解码。
- **`rayon`、`crossbeam-channel`、`indicatif`** — 并行、通道、进度条。
- **`noodles`、`rust-htslib`** — BIO 格式与 BAM/CRAM 支持。

## 2. main.rs — 命令分发与参数解析（重点）

`src/main.rs` 是 `impg` 二进制的入口，单文件约**613 KB / 1.5 万行+**， 承担了所有 clap 命令定义、
参数解析、stage 解析、命令分发与大量辅助逻辑。这是该项目最显著的工程特点（也是值得反思之处）。

### 2.1 文件顶层结构

文件大致按以下顺序组织（行号近似）：

| 区段                          | 行范围     | 内容                                                           |
|-------------------------------|------------|----------------------------------------------------------------|
| imports 与 `parse_*` 工具函数 | 1–250      | 度量后缀/syng/render 引擎解析                                  |
| GFA 简写解析                  | 250–4486   | `apply_gfa_output_engine_shorthand` 与各种 stage 解析器        |
| `GenotypeCommand` 子枚举      | 4486–4700  | `genotype` 下的 `Cos` (cosigt) 子命令                          |
| 顶层 `Args` enum              | 4706–6160  | 所有 20 个子命令的 clap 定义                                   |
| `main` / `run` 分发           | 6161–10250 | `match args { ... }` 巨型 dispatch                             |
| 辅助函数                      | 10250–end  | `build_graph_config`/`validate_*`/`initialize_threads_and_log` |

### 2.2 顶层 `Args` enum 与子命令清单

[main.rs#L4706](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L4706) 处的
`#[derive(Parser)] enum Args` 定义了所有子命令。下面按功能归类：

- **索引**
    - `Index` — 构建 `.impg` 索引（single 或 per-file 模式）；逻辑在 `main.rs` 内联 +   `impg.rs`/
      `multi_impg.rs`。
- **查询/投影**
    - `Query` — 把目标区间通过比对网络投影（核心命令）；路由到 `impg_index.rs` 的 trait 方法。
    - `Project` (别名 `projection`) — 投影命令；`projection.rs` + `projection/converter.rs`。
    - `Map` — 把短读段映射到 syng 索引（GAF/PAF/pack/proj 输出）；syng 后端。
- **分割/拼接**
    - `Partition` — 把 cohort 切成窗口化 loci；`commands/partition.rs`（含
      `rehome_singleton_slivers`）。
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
    - `Crush` — 解析 bounded bubble，支持多种 POA/POASTA/sweepga 路由； `resolution.rs`（§3.5，15
      种 `ResolutionMethod`）。
    - `Gfa2vcf` (别名 `gfa-to-vcf`, `povu`) — GFA → VCF；`main.rs` 内联 + povu 库。
    - `DescribeGraph` (别名 `describe-gfa`) — 输出图特征 Markdown/JSON/TSV 报告；
      `graph_report.rs`（15+ 维度）。
    - `Render` — 用 gfalook 渲染 1D 图；`commands/render.rs` + `render_bundle.rs`。
    - `Align` — 调用 wfmash/FastGA 跑比对；`commands/align.rs` +
      `commands/mod.rs::create_aligner_adaptive`。
- **syng 后端**
    - `Syng` — 从 FASTA/AGC 构建 syng 索引；`syng.rs` + `syng_parallel.rs`。
    - `SyngRepair` — 重建 `.syng.pstep`/`.syng.spos` 而不重读序列；`syng.rs`。
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

## 3. 核心数据结构 (Impg 隐式图)

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

### 3.5 `resolution.rs` — crush 算法实体（17169 行，全项目最大源文件）

`src/resolution.rs` 是 `crush` 子命令的实体实现，也是 impg 解决 bounded bubble 的核心。文件开头
明言："Bubble-guided graph resolution primitives... detect path-supported bubbles in a blunt GFA,
replace bounded single-entry/single-exit bubbles with exact path-preserving local graph induction,
and repeat until no eligible unseen bubbles remain." 关键设计：**不丢路径**（emitted paths are the
coordinate system）、**不做 lossy representative collapse**。

**入口函数** [resolution.rs#L997](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/resolution.rs#L997)：

```rust
pub fn resolve_gfa_bubbles(gfa: &str, config: &ResolutionConfig) -> io::Result<ResolvedGfa>
```

接收 GFA 1.0 字符串（仅 `S`/`L`/`P` 记录，link 必须是 blunt `0M`），解析为内部 `Graph`，调用
`resolve_graph_bubbles` 反复检测+替换+polish，直到无候选 bubble 或达到 `max_iterations`。返回
`ResolvedGfa { gfa: String, stats: ResolutionStats }`。

**`ResolutionConfig` 关键字段**
（[resolution.rs#L35](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/resolution.rs#L35)）：

- `max_iterations` — frontier 替换轮数上限（默认 1，`until-done` 解析为 `usize::MAX`）。
- `method: ResolutionMethod` — 替换算法选择（见下）。
- `auto_spoa_max_traversal_len` / `auto_poasta_max_traversal_len` — `method=auto` 时按中位遍历长度
  三档路由：`median < spoa` → sPOA；`spoa ≤ median < poasta` → POASTA；`≥ poasta` → sweepga。设为
  0 可禁用对应档。
- `auto_allwave_max_total_sequence` / `auto_allwave_max_traversals` — legacy 兼容字段，CLI 仍接受
  但 median 三档路由不使用。
- `motif_*` 系列 — `MotifLocal` 方法的 sparse offshoot 发现参数。

**[resolution.rs#L274](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/resolution.rs#L274)）**

这是 crush 算法最复杂的部分——同一套 bubble 替换框架支持 15 种 aligner/路由策略：

- `Auto` — 默认，按 median traversal length 三档路由 sPOA/POASTA/sweepga。
- `Poa` / `Poasta` / `Abpoa` — 直接 POA 替换（spoa_rs / poasta / abPOA）。
- `StarBiwfa` — 调试用：每条遍历对齐到 root，星形列图。path-preserving 但非质量默认。
- `Allwave` — sparse many-to-many BiWFA + AllWave + seqwish + SPOA polish。
- `Sweepga` — SweepGA/FastGA pair selection + seqwish + SPOA polish（默认 aligner=fastga）。
- `Wfmash` — 同 Sweepga 但 pin aligner=wfmash，用于与 PGGB 对比。
- `Hierarchical` — 按深度路由：level 0 → sweepga+seqwish，level ≥ 1 → POASTA。
- `ChainGreedy` — 贪心参考路径走，path-adjacent 连续 bubble 组成链，整链 POASTA。 （走独立的
  `resolve_graph_bubble_chains` 路径，绕过 POVU provenance tree。）
- `ChainPovu` — POVU 子树链分块 + smoothxg 风格局部 + bounded POASTA cleanup。
- `TopFlubbleSweepga` — 仅 level-0/root flubble，每个跑一次 all-vs-all SweepGA/seqwish。
- `IterativeMultiLevel` — 激进搜索：多 POVU 视图生成候选，多站点窗口先过 SweepGA/seqwish。
- `CoverageMultiBubble` — 覆盖驱动：iterative multi-level + outward residual + bp 加权覆盖报告。
- `MotifLocal` — motif 局部残差：从 path-step 支持度发现 sparse singleton/offshoot motif。

**`ResolutionPolishMethod`（3 变体）** — 替换后的 polish：`Poa`（POVU 重分解 + SPOA）、
`Poasta`（POVU 重分解 + POASTA）、`Smooth`（smoothxg 排序块平滑）。

**`MultiLevelWindowMode`（5 变体）** — `Largest`/`Parent`/`Sibling`/`Sliding`/`Local`，控制
`IterativeMultiLevel`/`CoverageMultiBubble` 的候选窗口生成策略。

**关键常量**
（[resolution.rs#L505-L569](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/resolution.rs#L505)）：
`DEFAULT_MAX_ITERATIONS=1`、`DEFAULT_MAX_TRAVERSAL_LEN=10_000`、
`DEFAULT_AUTO_SPOA_MAX_TRAVERSAL_LEN=1_000`、`DEFAULT_AUTO_POASTA_MAX_TRAVERSAL_LEN=10_000`、
`DEFAULT_REPLACEMENT_SEQWISH_MIN_MATCH_LEN=311`、`DEFAULT_MAX_PAIR_ALIGNMENTS=10_000`、
`DEFAULT_MAX_REPLACEMENT_PAF_BYTES=64MB` 等——这些数字直接对应 docs/crush-architecture-spec.md 的
Phase-2 阈值，是 C4 难用例调参的产物。

**内部数据结构**：`Graph`/`Segment`/`Path`/`Step`/`BubbleCandidate`/`BubbleNode`/`FlankContext`/
`TrimPlan`/`ReplacementPlan`/`CandidateFrontier`/`ChainCandidate` 等 20+ 个 struct/enum，构成完整的
bubble 检测→分块→替换→重写流水线。`DEBUG_REPLACEMENT_ID`/`DEBUG_APPLIED_FRONTIER_ID` 等原子计数器
为 `IMPG_CRUSH_DEBUG_DIR` 诊断提供唯一 ID（被 §7.6 的 `audit_poasta_replacement_cycles.py` 消费）。

**与 `pgr` 的相关性**：`resolution.rs` 是 impg 泛基因组部分最值得深读的单一文件。它把"bubble 检测

- 多 aligner 路由 + 迭代 polish + path 保真"完整实现了一遍，15 种 `ResolutionMethod` 的演进史
  本身就是 crush 算法的实验日志。`pgr` 若引入类似的 bubble 处理，可参考其 POVU provenance tree +
  median-length 三档路由的核心思路，但应避免 15 种 method 的复杂度膨胀——impg 自己的 C4 sweep 脚本
  （§7.6）显示大部分场景只用其中 3-4 种。

## 4. syng 免比对后端（略）

`impg` 还有一个 syng syncmer GBWT 免比对后端（`SyngIndex` / `SyngMatcher` / `SyncmerParams`，6
个 sidecar 文件 `.1khash`/`.1gbwt`/`.syng.names`/`.syng.pstep`/`.syng.spos`/`.syng.meta`，以及
`impg map` 的 GAF/PAF/pack/proj 输出）。需要时直接阅读 `impg-0.4.1/src/syng.rs` 与 `docs/` 下 syng
相关文档。

## 5. GFA 构建管道

`graph`、`query -o gfa`、`partition -o gfa` 三个命令共享同一套引擎实现，由 `--gfa-engine`
选择 ([lib.rs#L43](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/lib.rs#L43)的 `GfaEngine`
枚举)：

| Engine                   | Pipeline                                    | 用途                  |
|--------------------------|---------------------------------------------|-----------------------|
| `Pggb` (默认)            | sweepga + seqwish + smoothxg 平滑 + gfaffix | 平滑变异图            |
| `Seqwish`                | sweepga + seqwish + gfaffix                 | 原始（未平滑）图      |
| `Poa`                    | 单遍 SPOA                                   | 小区域、快速 MSA 输出 |
| `SyngNative`/`SyngLocal` | syng 锚点 + BiWFA + allwave + seqwish       | 近缘单倍型快速通道    |

### 5.1 stage 化管道

通过 `-o gfa:<stage1>:<stage2>:...` 简写，用户可在引擎前后插入 stage（见 §2.4）：

```
-o gfa:cut-n=100:pggb:crush,method=allwave:sort,pipeline=Ygs
       ^^^^^^^^^^ ^^^^ ^^^^^^^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^
       终端N裁剪  引擎  bubble 解析           最终排序
```

`build_graph_config` 与 `build_engine_opts` 把这些 stage 装配成 `EngineOpts` 结构，传给
`commands::graph::run` / `commands::partition::run` /`commands::syng2gfa::run`。

### 5.2 partitioned 模式与 `lace`

引擎名后追加 `:WINDOW`（如 `pggb:10000`）即进入分区模式：

1. 把目标区域按 `WINDOW` bp 切分。
2. 每个窗口独立构建 GFA（峰值内存受控）。
3. 最终用 `impg lace` 把 per-window GFA 拼回一张图，可选 `--fill-gaps` 用参考序列填充窗口间空隙。

`lace` 同时支持 GFA 与 VCF 输入，路径名必须遵循 `NAME:START-END` 约定（最后一个 `:` 是分隔符），
坐标驱动重新拼装。

### 5.3 GFA 管道源文件映射

GFA 构建管道分散在 8 个源文件中，职责边界如下：

- **`lib.rs`** — `GfaEngine` 枚举（5 变体）+ `EngineOpts`（已解析的引擎配置）+
  `SmoothPipelineConfig`（post-crush smoothxg pass 的配置，与 pggb 默认值一致：
  `target_poa_lengths=[700,1100]`、`max_node_length=100`）。
- **`graph_pipeline.rs`** — `GraphPipelineSpec` 解析器：把 `stage,key=value:stage,...` 字符串
  解析为类型化的 `Vec<GraphPipelineStage>`，只做语法校验，不决定可执行性。
- **`commands/graph.rs`** — 管道执行入口：`build_graph`/`induce_graph_from_alignment`（seqwish
  传递闭包）、`run_graph_build`/`run_graph_build_poa`/`run_graph_build_pggb`/
  `run_graph_build_partitioned`（分区模式）。`GraphBuildConfig` 含 20+ 字段（threads/frequency/
  min_aln_length/repeat_max/min_match_len/adaptive_min_match_len/sparse_factor/transclose_batch/
  disk_backed 等）。
- **`graph.rs`** — 图操作原语：`unchop_gfa`/`sort_gfa`（gfasort Ygs 集成）、 `build_spoa_engine`/
  `feed_sequences_to_graph`（SPOA）、`reverse_complement`、`terminal_n_clip_span`、
  `prepare_poa_graph_and_sequences`。
- **`resolution.rs`** — crush 阶段实体（§3.5），由 `parse_crush_stage` 装配 `ResolutionConfig`
  后调用 `resolve_gfa_bubbles`。
- **`smooth.rs`** — smoothxg 风格块分解 + SPOA 平滑，`SmoothConfig` 含
  `target_poa_lengths`/`max_node_length`/`poa_padding_fraction`/`scoring_params`/
  `block_source: SmoothBlockSource`（PathOverlap/Flubble/NeighborMergePoasta 三策略）。被
  `parse_smooth_stage` 装配，作为 post-crush stage 执行。
- **`gfa_self_loops.rs`** — `NormalizeSelfLoops` 命令实体：折叠 blunt GFA 路径局部的 self-loop
  重复单元，输出 `NormalizeSelfLoopsStats` 报告。
- **`commands/syng2gfa.rs`** — syng → GFA 物化（`SyngNative`/`SyngLocal` 引擎），频率感知 syncmer
  节点共享策略。`SyngGfaMode`（blunt/raw）控制输出形态：blunt 物化精确 source-spelling 0M paths；
  raw 输出 syng-native overlap 图。

**`graph_report.rs`**虽然不直接参与构建，但为所有 sweep 脚本提供评分依据： `GraphReportOptions`
含 15+ 阈值（`max_link_jump_frac`/`min_largest_component_frac`/`min_common_start_frac`/
`max_internal_tips`/`warn_duplicate_sequence_frac`/`min_white_space_gap_bp`/
`max_sparse_coverage_path_fraction` 等），`GraphReport` 输出 `status`/`failures`/`warnings` +
各维度测量值，序列化为 JSON 供 `c4-crush-cmaes.py` 等优化器消费。

### 5.4 `examples/` 诊断与实验工具（7 个独立可执行程序）

`examples/` 目录下有 7 个 `.rs` 文件，每个都是独立的 `cargo run --example <name>` 程序，
不是集成测试。它们主要服务于 crush/smooth 算法的**离线诊断与实验**，依赖 `impg::` 库 API，
不属于正式 CLI。

**path-preserving 校验类（2 个，crush 算法的端到端验证）**

- [`compare_gfa_paths.rs`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/examples/compare_gfa_paths.rs)
  (57 行) — 比较两个 GFA 的路径拼写是否完全一致。调用 `impg::resolution::path_sequences`，
  输出 `missing/extra/spelling_mismatches` 三类计数，退出码反映一致性。用法：
  `compare_gfa_paths <expected.gfa> <observed.gfa>`。
- [`validate_gfa_path_sources.rs`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/examples/validate_gfa_path_sources.rs)
  (104 行) — 验证 GFA 路径名 `NAME:START-END` 的拼写是否与源序列文件一致（含反向互补检查）。
  调用 `path_sequences` + `UnifiedSequenceIndex::fetch_sequence`，输出
  `forward/reverse_complement/unparsable/spelling_mismatches` 计数 + 首个 mismatch 上下文。用法：
  `validate_gfa_path_sources <graph.gfa> <sequence.fa|sequence.agc>...`。

**POVU/POA 诊断类（3 个，crush 算法内部结构探查）**

- [`povu_decomp_report.rs`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/examples/povu_decomp_report.rs)
  (74 行) — 用 `povu::NativeGfa` 库分解 flubbles，输出每个 site 的
  `id/parent_id/level/is_leaf/ref_start_step/ref_end_step/start/end` TSV，按
  ref_span 降序排列。用于理解 crush 算法的 POVU provenance tree 结构。用法：
  `povu_decomp_report <graph.gfa> [reference-name]`。

- [`neighbor_merge_existing_gfa.rs`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/examples/neighbor_merge_existing_gfa.rs)
  (141 行) — 对已有 GFA 跑 `SmoothBlockSource::NeighborMergePoasta` 平滑 + gfaffix + Ygs sort
  三步管道。这是 `smooth.rs` 的独立调用入口，用于在不跑完整 crush 的情况下对已有图做后处理。
  自动统计 PanSN haplotype 数，配置 `target_poa_lengths=[target_bp; iterations]`。用法：
  `neighbor_merge_existing_gfa <input.gfa> <output-prefix> [iterations=3] [target-bp=10000] [threads=32]`。

- [`poasta_order_driver.rs`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/examples/poasta_order_driver.rs)
  (~560 行，最大的 example) — **POASTA/abPOA 输入顺序实验驱动器**。测试 8+ 种序列顺序策略对 POA
  输出图质量的影响：

    - `longest_then_name` / `reverse_longest_then_name` — 长度基线
    - `medoid_first_then_longest` — k-mer Jaccard medoid 优先
    - `nearest_neighbor_guidetree` — 贪心最近邻引导树
    - `reference_chm13_first` — CHM13/GRCh38 优先（若存在）
    - `random_seed_N` — 确定性随机洗牌（splitmix64）
    - abPOA 变体（需 `--abpoa-bin`）每种顺序调用 `poasta_sequences_to_gfa_in_order` / `abpoa_sequences_to_gfa_in_order`， 再用
  `describe_gfa` 评估 25+ 指标（segments/links/cycles/coverage/white_space/depth 等），输出
  `orders.tsv` + `summary.tsv` + `skipped.tsv`，可选 `--render-svg` 调 gfalook 渲染。这是对应
  `docs/poasta-order-*.md` 文档的实验脚本，研究 POA 输入顺序敏感性。

**syng 探测类（2 个，本文档不参考）**

- `syng_anchor_probe.rs` (94 行) — 打印 `query_region_with_anchors` 返回的锚点位置。
- `syng_probe.rs` (118 行) — 对比 syng C 迭代器与 `query_region_with_anchors` 的位置一致性。

**与 `pgr` 的相关性**：`compare_gfa_paths` / `validate_gfa_path_sources` 是任何"路径保真图操作"都
需要的校验工具，`pgr` 若引入 GFA 操作可直接复用这两个 example 的思路（path 拼写 vs 源序列双向校验）。
`poasta_order_driver` 的实验框架（多策略 × 多指标 TSV 输出）也值得借鉴——它把"POA 顺序敏感性"
这个模糊问题变成了可排序的指标表。

### 5.5 `tests/` 测试体系（13 个集成测试 + 2 个验证脚本 + fixture 矩阵）

`impg` 没有传统的 `src/` 内单元测试，所有测试都在 `tests/` 目录下，共 13 个 `test_*.rs` 集成测试 文件
（合计 10866 行）+ 2 个 shell 验证脚本 + 一个中心化 fixture 矩阵。这与 `pgr` 的 `tests/cli_*.rs`
约定形成鲜明对比（详见 §6）。

**测试代码（`tests/*.rs`，13 文件，按规模降序）**

- `test_syng_integration.rs` (5302 行，最大) — syng CLI 端到端：`impg syng --agc` 建索引、round-trip
  查询、`impg partition -a <syng-prefix>`。`SYNG_LOCK` 全局 Mutex 串行化（C 库非线程安全）。
  本文档不参考。
- `test_crush_integration.rs` (1674 行) — crush 回归测试，对应 `docs/crush-audit.md` 的编号失败。
  真实数据：`c4_slice_1500_3000.gfa`（2942 segments/64k bp，C4A chr6:31891045-32123783，465
  haplotypes，seqwish-k=311）。已知 RED 检查用 `#[ignore]` 标记保持 CI 绿，文档注释说明复现命令。
  32 线程 `Once` 初始化 + `Mutex` 串行 C4 fragment 测试。
- `test_genotype_validation_suite.rs` (777 行) — syng pack/genotype/infer 的 truth-known 验证，
  断言精确 evidence vector 与基因型顺序，合成小面板。
- `test_transitive_integrity.rs` (767 行) — 传递闭包完整性：非重叠区域保持分离、坐标投影准确、
  双向查询对称、identity 过滤、传递链坐标准确。
- `test_local_compression_testbed.rs` (645 行) — 消费 `test_data/local_compression/manifest.json`
  矩阵（见下），13 个 `REQUIRED_CLASSES` × 12 个 `REQUIRED_METHODS` 笛卡尔积，每个 fixture 校验
  path spelling + 拓扑断言 + `allowed_ranges`（segment/link/depth/bubble/flubble/long_link 等
  min-max 区间）。这是 crush 算法最系统的回归网。
- `test_syng_startcount.rs` (351 行) — 隔离 `syngBWTpathStartNew` 的 startCount 自增行为。
- `test_pipeline_integration.rs` (268 行) — 全管道：index → partition → graph → lace，用酵母 chrV 7
  株系数据，`#[ignore]` 因依赖 wfmash + samtools。
- `test_genotype_gfa.rs` (226 行) — `impg genotype gfa` 的 `GraphContributionModel`/
  `GraphFeatureIdMode` 报告 section 提取与断言。
- `test_gfa_projection.rs` (206 行) — GAF→GFA projection 转换，`project_gaf_to_gfa` +
  `GfaProjectionOutputFormat`（ProjectionBundle/PackTsv），含 `tiny_graph` 内联 fixture。
- `test_agc_integration.rs` (186 行) — AGC vs FASTA 同内容校验，`AgcIndex` vs `FastaIndex` 交叉验证
  5 个 test case。
- `test_graph_seqwish.rs` (172 行) — seqwish 图诱导，内联 C4A/C4B 226bp 序列 + 单行 PAF
  （`cg:Z:65=1X160=`），验证 `chain_41`/`scaffold` 标签处理。
- `test_graph_poa.rs` (167 行) — `impg graph --gfa-engine poa` 端到端，临时 FASTA → GFA → path
  校验。
- `test_graph_output_crush.rs` (125 行) — `impg graph -o gfa:pggb:crush` 端到端，依赖 gfaffix
  binary。

**测试数据（`tests/test_data/`，三层 fixture）**

根目录小 fixture（`a.fa`/`b.fa`/`c.fa`/`ref.fa`/`ref2.fa` + `.fai`，30-80 字节级）用于基础索引查询；
`test.agc` (1KB) 用于 AGC 测试；`yeast.chrV.fa.gz` (1.2MB) 用于全管道测试。

`crush/` 子目录（6.9MB）—— crush 算法真实数据 fixture：

- `c4_slice_1500_3000.gfa` (7.2MB) — C4A 切片，最大的真实 GFA。
- `nested_bubbles_real.gfa` / `small_insertion.gfa` / `small_insertion_walks.gfa` — 小型 GFA。
- `top_flubble_seqwish_minrun.fa` + `.paf` — top-flubble seqwish 失败回归。
- `c4_fragments/` — 4 个 C4 子片段（`easy_shared_flank`/`bounded_multi_bubble`/
  `short_floor`/`duplicated_repeat`），每个 `.fa`+`.paf`，README.md 说明来源于
  `/home/erikg/impg/data/c4_top_flubble_fix_skipzero_20260527T185505Z/debug`。

`local_compression/` 子目录（216KB）—— **最有组织的测试数据集**。由
`scripts/local_compression_testbed.py write-fixtures` 生成，`manifest.json`
中心化管理（schema `local_compression_fixture_manifest_v1`）。13 个 fixture class
覆盖所有设计要求的局部压缩场景：`snp_bubble` / `short_indel` / `insertion_50_500bp` /
`alu_like_insertion` /`adjacent_bubbles_compress_together` / `bubble_split_by_fake_repeat_anchor`
/`repeated_motif_microtangle` / `duplicated_flank_requires_path_context` /
`tandem_copy_number_loop_cyclic` / `dispersed_repeat_glue_break_or_ignore` / `inversion_like` /
`nested_bubbles_top_level_right` / `nested_bubbles_top_level_wrong`。

每个 fixture 子目录含 4 件套：

- `input.fa` — 输入 FASTA（路径名遵循 `sample#hap#contig:START-END`）。
- `expected_paths.tsv` — 期望路径拼写（path-preserving 校验基准）。
- `metadata.json` — fixture 元数据，含 `expected_topology`（断言 ID + 描述）、
  `allowed_ranges`（segment_count/link_count/path_depth_median/node_depth_p95/
  white_space_proxy_bp/bubble_count/flubble_count/long_link_count/self_loop_count 等 min-max
  区间）、`long_link_policy`、`known_failure_mode`、`render_hints`。
- `notes.md` — 场景简述。

fixture 有 tier 标签：`ci`（快速 CI 子集）或 `local`（本地完整测试）。

**验证脚本（`tests/validation/`，2 个 bash）**

- `battery_syng_vs_paf.sh` (116 行) — syng vs PAF 批量对比：随机生成 N 个 2-5kb 查询区域， 每区域输出
  `syng_rows/paf_rows/common/syng_only/paf_only/mean_ds/max_abs_ds/mean_de/ max_abs_de/pct_start_5bp/pct_end_5bp`
  共 11 列 TSV。诊断工具，退出码恒 0。
- `compare_syng_vs_paf.sh` (110 行) — syng vs PAF 单次对比：行数/覆盖/首个 mismatch 上下文。

**测试基础设施约定**

- 二进制定位：所有测试用 `CARGO_BIN_EXE_impg` 环境变量（cargo 自动注入），回退到
  `target/{debug,release}/impg`，再回退到 `/home/erikg/impg/target/release/impg`（作者机器路径）。
- 外部依赖：gfaffix 必须存在（`CARGO_BIN_EXE_gfaffix`，测试会 copy 到 impg 旁边作 sibling）；
  wfmash/samtools 可选（`#[ignore]` 标记）。
- 串行化：syng 测试用 `LazyLock<Mutex<()>>` 全局锁（C 库非线程安全）；C4 fragment 测试用 `Once` +
  `Mutex` 32 线程池串行。
- 已知失败：`#[ignore]` 标记的 RED 检查保持 CI 绿，文档注释说明复现命令与对应 docs/ 审计条目。

**与 `pgr` 测试约定的对比**（详见 §6）：`pgr` 用 `tests/cli_<command>.rs` + `PgrCmd` 辅助结构体
（`tests/common/mod.rs`），测试数据放 `tests/<command>/`，强调 Zero Panic。impg 则把测试数据
中心化到 `manifest.json` + tier 分级 + `allowed_ranges` 区间断言——这套做法对"多方法 × 多场景"
的回归矩阵特别有效，`pgr` 若引入 crush/POA 类算法可借鉴。

## 6. 对比分析: impg vs pgr

`pgr` 的真正强项是**UCSC 体系的 pairwise 比对处理**（Chain/Net/MAF/AXT/PSL/LAV 全套， 见
[docs/chain.md](file:///Volumes/ExtHome/Scripts/pgr/docs/chain.md)）与**Block FA 多序列比对**
（`fas` 全套子命令 + `libs/poa/` 的 SPOA 移植 + `libs/fas_multiz.rs` 的 multiz 风格 banded DP 合并，
其 `FasMultizMode::Core` 即"多基因组共享 core 比对"）。pairwise 与 core 比对均已成熟。

**泛基因组图部分是 `pgr` 的空白**：[docs/gfa.md](file:///Volumes/ExtHome/Scripts/pgr/docs/gfa.md)
明确写"如果 `pgr` 未来涉及泛基因组操作"，是规划/知识背景文档而非实现；`src/cmd_pgr/` 下无 `gfa`
子命令，`src/libs/` 下无 GFA 模块。本节聚焦"泛基因组图"这一维度对比，作为 §7 启示的依据。

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

### 深度对比：区间投影 vs Chain lift

`pgr chain lift` 通过 Chain 的坐标映射把目标区间 lift 到查询序列，本质上是单条 Chain 的线性投影。
`impg query` 在 all-vs-all 比对网络上做区间树查找 + 传递闭包，等价于"在所有 Chain 的并集上做 BFS"。

- **优势**：impg 自动发现所有同源片段（包括间接通过第三序列的同源），pgr chain lift 需要用户手动选
  Chain。
- **代价**：impg 需要 all-vs-all 比对（O(n²) 内存与时间），pgr chain lift 直接用现成的 UCSC Chain
  （已经经过 syntenic 净化）。

## 7. 对 `pgr` 的启示（聚焦泛基因组部分）

`pgr` 的 pairwise 比对与多基因组 core 比对已成熟，无需从 impg 借鉴；且 `libs/poa/` 已有完整的
SPOA 移植（`Poa` struct，`add_sequence`/`consensus`/`msa`），被 `fas consensus` 与 `fas refine`
使用——这恰好是 impg crush 算法的 aligner 基础设施。下列启示聚焦**泛基因组部分**——即把已有比对
升级为 cohort 级隐式图、按需投影、bounded bubble 处理等。

1. **区间树 + 紧凑 CIGAR delta 的组合值得借鉴**：`pgr` 处理 PAF 时可直接复用 `coitrees` + `CigarOp`
   风格的 bit-packing，把全基因组 all-vs-all 比对装进内存。当前 `pgr` 的 PSL/Chain 处理是流式的，
   缺少随机访问能力——这是把 pairwise/core 比对升级为"全 cohort 可查询隐式图"的前提设施。
2. **trait 抽象单/多文件索引**：`ImpgIndex` trait + `MultiImpg` 是处理"单大文件 vs
   多小文件"两种部署模式的干净做法。`pgr` 若引入类似的索引层，可让命令代码与索引物理形态解耦。
3. **stage 化字符串 DSL 的两面性**：impg 的 `-o gfa:cut-n=100:pggb:crush:sort` 简写表达力强，
   但代价是 `main.rs` 膨胀到 60 万字符。`pgr`当前坚持 `pgr <format> <subcommand>`
   的多级结构更易维护，应保持。若未来需要类似的管道组合，可考虑专门的 pipeline 配置文件而非 CLI
   简写。
4. **避免 main.rs 巨型化**：impg 把 20 个子命令的 clap 定义与分发全塞进单文件是明显的反例。`pgr`
   的 `src/cmd_pgr/` 按格式/功能分组、每命令独立模块的结构更优，应继续坚持——`main.rs` 只做
   `ArgMatches` 分发，业务逻辑下沉到模块。
5. **thread-local 缓存模式**：impg 用 `thread_local!` 缓存 WFA aligner、 1aln/TPA 句柄、
   目标序列片段，避免重复分配。`pgr` 在并行处理 PSL/Chain 时可借鉴同样的模式。
6. **PAF `cg:Z:` 懒加载**：impg 不把 CIGAR 存进区间树节点，只存 `data_offset` + `data_bytes`，
   查询时按需读取。这对 `pgr` 处理大型 PAF 是直接可借鉴的内存优化。
7. **PanSN 命名约定**：impg 全程用 `sample#haplotype#contig` 命名（`#` 分隔），`pgr` 在 `pgr pl`
   流水线中若需要处理群体数据，可采用同一约定以与 pggb/impg/odgi 生态兼容。
8. **Zero Panic 与 AGENTS.md 的契合**：impg 源码中存在大量 `unwrap_or_else(|e| panic!(...))`（如
   `get_cigar_ops`、`get_target_sequence_cached`），违反了 `pgr` 的 Zero Panic 原则。`pgr`
   在借鉴其算法时应改为 `anyhow::Result` + `bail!`，把错误返回到调用方而非 panic。
9. **POA 基础设施可直接复用**：`pgr` 的 `libs/poa/`（`Poa` struct + `AlignmentParams` +
   `AlignmentType::Global`）已是 crush 算法 aligner 层的现成基础。impg crush 的核心是 "bubble 检测
   → aligner 路由（SPOA/POASTA/abPOA）→ polish"，`pgr` 若实现类似的 bubble resolution，可直接在现有
   `Poa` 之上扩展 `ResolutionMethod` 路由层，无需重新移植 aligner。`fas consensus` 已验证该 POA 在
   MSA 场景可用，迁移到 graph bubble 场景的门槛低于从零起步。

## 7.5. `impg-0.4.1/notes/` 目录

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

## 7.6. `impg-0.4.1/scripts/` 目录

`scripts/` 是 21 个独立可执行脚本（Python 13 个 + Bash 5 个 + R 1 个），几乎都是
**实验驱动器 (experiment driver)**：它们不属于 impg 二进制，而是封装"调用 impg 子命令 → 跑
sweep/参数搜索 → 收集时间/RSS/graph-report → 写 TSV/JSON → 上传 hypervolu.me 渲染"的工作流。
所有脚本都把重型 artifact 写到仓库外的 `/home/erikg/impg/data/<RUN_ID>/`，仅驱动器本身入库。
按用途分四类：

### (1) C4 难用例实验矩阵（10 个，主题占绝大多数）

围绕 GRCh38#0#chr6:31891045-32123783 (C4) 这个"超难"区域跑参数/算法矩阵，全部以 impg + graph-report +
`compare_gfa_paths` + gfalook 渲染为骨架：

- `c4-crush-cmaes.py` (64 KB) — CMA-ES 黑盒优化器，优化 crush 参数
- `c4-low-seqwish-k-sweep.py` (29 KB) — seqwish min-match k 值 sweep + POA polish
- `c4-highfreq-mask-crush-sweep.py` (28 KB) — high-frequency mask + crush 联合 sweep
- `c4-motif-local-polish.py` (37 KB) — 诊断 C4 残留 underaligned motif，局部 polish
- `c4_syng_tail_diagnosis.py` (7 KB) — syng 参数矩阵 ×C4，验证 FASTA 长度
- `run-c4-aggressive-motif-matrix.py` (27 KB) — C4 motif-window POA vs abPOA 对比
- `run-c4-hard-seqwish-k-sweep.sh` (17 KB) — hard 难度 seqwish k sweep，固定 baseline 对照
- `run-c4-k311-poa-threshold-sweep.sh` (21 KB) — k=311 下 POA 阈值 (500/1k/2k/5k/10k) sweep
- `run-c4-sweepga-abpoa10k.py` (25 KB) — SweepGA/seqwish seed + abPOA 10kbp motif polish
- `validate-flank-aware-crush-c4.py` (40 KB) — flank-aware crush 在 fixture 与 C4 上的验证

注：`c4-crush-cmaes.py` 明言 "optimizer is deliberately external to impg"——把候选生成与评分解耦，
评分逻辑不入主二进制。

### (2) crush 诊断与通用图 QC（3 个）

| 脚本                                   | 内容                                               |
|----------------------------------------|----------------------------------------------------|
| `audit_poasta_replacement_cycles.py`   | 消费 `IMPG_CRUSH_DEBUG_DIR`，输出 SCC/path TSV     |
| `graph-clean-qc.py` (14 KB)            | 拓扑 smoke test：无长程 link、tip 由端点解释       |
| `local_compression_testbed.py` (89 KB) | local compression fixtures（§8.10），写 scoreboard |

注：`graph-clean-qc.py` 明言 "not a biological truth check"；`local_compression_testbed.py`
是全目录最大脚本。

### (3) demo 与公共用法示例（4 个 bash + 1 个 R）

| 脚本                               | 内容                                                |
|------------------------------------|-----------------------------------------------------|
| `demo-gfa.sh` (15 KB)              | 验证 `query -o gfa` 与 `graph` 产出等价 GFA         |
| `demo-batching.sh` (14 KB)         | 验证 `--batch-bytes` 与 unbatched 产出等价          |
| `demo-partitioned.sh` (10 KB)      | 验证 full vs `pggb:10000` partitioned 等价          |
| `partitioning.sh` (3 KB)           | PAF+bedtools 窗口切分流水线（`makewindows` + mask） |
| `plot_partitioning_stats.R` (8 KB) | ggplot2 可视化 partition BED 统计                   |

### (4) 杂项工具（3 个）

| 脚本                           | 内容                                                   |
|--------------------------------|--------------------------------------------------------|
| `faln2html.py` (11 KB)         | FASTA alignment → HTML viewer（reactmsa/proseqviewer） |
| `procpaf.py` (3 KB)            | 用 `.fai` 把 subregion PAF 坐标 lift 回全染色体        |
| `hprcv2-syng-smoke.py` (10 KB) | syng 烟雾测试面板（多类窗口）                          |

### 共性观察

- **设计模式：driver 外置**。所有 sweep 脚本都把 impg 当黑盒 subprocess 调用，
  自身只做参数装配、artifact 落盘、TSV 汇总。这套思路对 `pgr` 未来的 benchmark/sweep
  工具直接可借鉴——把候选生成与评分解耦，评分逻辑不入主二进制。
- **artifact 路径强约定**：
  `/home/erikg/impg/data/<RUN_ID>/{graphs,reports,sorted,renders,logs,validation,debug}/`。这是
  single-developer 工作流的产物，`pgr` 若移植需参数化。
- **`compare_gfa_paths`**（`target/release/examples/`）是事实标准 baseline 比对工具，几乎每个 sweep
  都调用它做 path 序列保真度检查。
- **C4 占比失衡**：21 个脚本中 10 个专为 C4 一条区域服务，反映出 crush/C4 难用例是 impg
  开发后期的主战场，也提示该问题极具挑战性。
- **与 `pgr` 相关性**：除 `hprcv2-syng-smoke.py` 与 `c4_syng_tail_diagnosis.py` 属 syng
  主题不参考外，其余 19 个脚本的核心模式（driver 外置 + TSV 汇总 +path 保真检查 + graph-report
  评分）对 `pgr` 设计泛基因组图构建的 benchmark 体系有直接参考价值。

## 8. impg-0.4.1/docs 文档结构

`impg-0.4.1/docs/` 是项目的开发文档目录，规模庞大：
**111 个顶层 Markdown 文件 + 2 个子目录 (`designs/`、`evaluations/`)，合计约 3.5 MB / 34619 行**。
它不是面向最终用户的文档（那是 `README.md` 与 `--help`），而是开发过程中的设计笔记、实验报告、bug
诊断与审计记录。下面按实际内容主题分组梳理，每个文件给出真实的内容主旨（不只看文件名前缀）。

### 8.1 crush 算法 — 设计与规范（前瞻性设计文档，~21 个）

`crush` 是 impg 解决 bounded bubble 的核心算法：把比对/syng 产生的"bubble 化"图通过局部
condensation 压缩成适合下游分析/分型的紧凑图。这是 docs/中最大的主题，且算法本身经过多次重构，
所以设计文档非常多。

- **`crush-architecture-spec.md`** — crush 操作的**权威规范**。开篇明言："three days of
  agents working on crush have been flailing without a clear spec... No code change to crush
  should be proposed without referencing this document."定义 3 阶段算法：POVU flubble 检测 →
  按中位遍历长度路由 aligner → polish 直到收敛。目标是分钟级端到端，避免 all-to-all 比对。
- **`crush-design.md`** — crush 的最初设计稿，含 hot path 分析 （parse/render/validate
  是非对齐热路径）。
- **`crush-design-fix.md`** — 设计修复（针对早期实现的偏差）。
- **`crush-hierarchical.md`** — 层次化 crush：分层处理 bubble，level-0 用 SweepGA，深层用 POASTA。
- **`crush-level-descent.md`** — 层级下降算法（每轮重新 POVU 整图）。
- **`crush-true-level-descent.md`** — 真正的层级下降（修复 nested bubble 在 parent
  内部继续下降的问题）。
- **`crush-nested-bubble-test.md`** — nested bubble 回归测试
  `nested_bubble_level_descent_actually_descends`，固化"白色 ramp"可视化症状。
- **`crush-neighbor-merge-iterate.md`** — `:neighbor-merge-poasta` 阶段设计：把 path-adjacent
  bubble site 贪心合并成 reference-span 组，对每组跑 POASTA，重复 N 次（默认 3 次 10kb 迭代）。
- **`crush-per-bubble-isolation.md`** — 把每个 bubble 隔离成 FASTA 单独跑 SweepGA/seqwish 与
  wfmash+seqwish 对比，验证"应该完全压缩掉"的假设。
- **`crush-wider-context-bubbles.md`** — 给 aligner 加 wider context （flanking-aligner anchors），
  解决"mistakes propagate through hierarchy"问题。
- **`crush-chain-greedy-walk.md`** — greedy walk chain 设计。
- **`crush-chain-povu-tree-blocks.md`** — POVU tree blocks chain 设计。
- **`crush-flubble-guided-smoothing.md`** — flubble 引导的平滑。
- **`crush-identity-k-policy.md`** — identity-k 策略。
- **`crush-bail-removal.md`** — 移除 bail-out 逻辑（让所有 bubble 都尝试解）。
- **`flank-aware-crush-design.md`** — flank-aware crush 设计：短 indel/tandem
  motif/homopolymer/microtangle 的内部缺乏足够公共上下文，需要先确定可替换 target interval，再加
  occurrence-local 左右 flank 作为 resolver context。
- **`flank-aware-crush-quality-pass.md`** — flank-aware crush 的下游 WG 任务质量收紧。
- **`hierarchical-graph-resolution.md`** — 层次图解析总览。
- **`top-flubble-sweepga.md`** — `method=top-flubble-sweepga`：只对 POVU level-0 flubble 跑
  SweepGA，证明 SweepGA 能对大区域产生比对，但小区域仍欠对齐。
- **`local-graph-compression-testbed-design.md`** — 本地图压缩测试床设计： tiny 对抗性 fixture，
  用于快速暴露候选窗口选择、flubble 遍历、SweepGA/POA 等"是否会破坏 path"。
- **`local-graph-compression-testbed-quality-pass.md`** — 测试床下游任务质量收紧。
- **`remove-crush-replacement-quality-guards.md`** — 决定移除 crush replacement 的质量守卫：只保留
  correctness invariants（GFA 解析、path 名/顺序/拼写保留），compression ratio 等仅作诊断。

### 8.2 crush 算法 — 实验报告（~24 个）

每个实验通常假设 → 实验设置（binary commit、输出目录、PNG 链接）→ TL;DR → 详细结果。多数以 C4
GRCh38 chr6:31891045-32123783 为 canonical 测试区域。

- **`crush-experiment-synthesis.md`** — 7 个并行 crush 实验的**综合报告**（基于磁盘 GFA artifact 与
  `/usr/bin/time -v` 实测，无估算）。结论：aligner 选择（auto-routing）是 load-bearing knob。
- **`crush-aligner-speed-study.md`** — 关键发现：round-1 plan 2 (sPOA) 占 99.2% 时间（831s），
  POASTA 在相同输入上 9.93s（**84× 加速**），sPOA 被 POASTA 全面取代的起点。
- **`crush-exp-auto-k31.md`** / **`crush-exp-auto-k51.md`** — seqwish-k=31/51 对 auto-routing
  的影响。
- **`crush-exp-auto-fixed-filter.md`** — 固定 filter 实验。
- **`crush-exp-allwave-k31.md`** — AllWave 路由实验。
- **`crush-exp-sweepga-k31.md`** — 假设**被反驳**：bit-identical PAF 证明降低 seqwish-k 不能 rescue
  小 bubble。
- **`crush-exp-hybrid-sweepga-poasta.md`** — 2-tier 混合（POASTA+SweepGA，跳过 sPOA）。
- **`crush-exp-poasta-everywhere.md`** — 假设**确认**：POASTA 用于所有 bubble 是新最佳（51-200bp
  dup-extras 16→11）。
- **`crush-exp-min-run-5.md`** — 假设在**输入层被反驳**：min-run=5 与 min-run=3 mask 掉完全相同的 2
  个 syncmer 节点。
- **`crush-sweepga-everywhere-unfiltered.md`** — sweepga 全部 + no-filter。
- **`crush-sweepga-many-to-many-poasta-polish.md`** — many-to-many + POASTA polish。
- **`crush-smoothxg-on-output.md`** — smoothxg 后处理。
- **`crush-raw-fastga-seqwish.md`** — raw FastGA + seqwish。
- **`crush-gfaffix-run.md`** — gfaffix 后处理。
- **`crush-vs-pggb-comparison.md`** — 与 PGGB 控制对比：PGGB 13:38/64GB 产出 13288 seg / 234524
  bp；crush 当前最佳 36:53/101GB 产出 19836 seg / 553585 bp （**2.36× 更多序列**），并诊断 ≤200bp
  polyT/SNP bubble 未压缩的 3 个具体例子。
- **`crush-wfmash-replacement.md`** — 用 wfmash 替代 aligner（用户建议 "we should try WFMASH on the
  bubbles"）。
- **`eval-crush-smoothxg-on-2.md`** — smoothxg 评估 #2。
- **`c4-highfreq-mask-crush-sweep.md`** — 高频 mask 实验。
- **`coverage-driven-repeat-c4.md`** — coverage-driven repeat 检测。
- **`c4-hard-seqwish-k-sweep.md`** — seqwish-k sweep。
- **`c4-whole-region-sweepga-seed.md`** — 整个 C4 区域 SweepGA seed 图。
- **`local-syng-parameter-sweep.md`** — 本地 syng 参数扫描：测试 C4 残留碎片化是否由全局 HPRCv2
  syng 索引（k=63,s=8）引起；`gfa:syng-local` 提取 query-selected 序列重建区域 syng 索引。
- **`calibrate-spectrum-driven-20260530.md`** — spectrum 驱动校准（带日期快照）。

### 8.3 crush 算法 — 诊断与审计（~19 个）

回溯性文档，记录 bug 调查、性能瓶颈、行为审计。

- **`crush-audit.md`** — 真实 C4 数据行为审计：HEAD `90ba74f` 在 canonical 命令上 SIGTERM 超时（30
  min 无输出），对比 known-good baseline (`0af1a4c`之前，4:39 wall / 55GB RSS)。
- **`crush-spec-audit.md`** — spec 审计（实现是否匹配 `crush-architecture-spec.md`）。
- **`crush-quality-state.md`** — 质量状态迭代追踪：pggb-in-impg 参考图 vs
  syng→sweepga+seqwish-k311→crush auto 目标图的 gap 收敛过程。
- **`crush-verify-report.md`** — 5 个真实 GFA 输入 × 3 scale 的端到端验证，path 序列保留 ✓，
  whitespace 减少 18-66%。
- **`crush-perf-report.md`** — 性能报告（AMD EPYC 7713 × 2, 256 logical CPU）， 优化非对齐热路径。
- **`crush-big-bubble-diag.md`** — 大 bubble 诊断。
- **`crush-aligner-deep-diag.md`** — aligner 深度诊断。
- **`crush-aligner-failure-trace.md`** — aligner 失败追踪：用 `IMPG_CRUSH_DEBUG_DIR` 抓
  per-replacement 子图，无源码改动。
- **`crush-fragment-source-trace.md`** — fragment 来源追踪。
- **`crush-fix-routing.md`** — 路由修复。
- **`crush-fix-sweepga-short-filter.md`** — short filter 修复。
- **`crush-fixtures-redproof.md`** — fixtures 验证。
- **`crush-poasta-pass-through-audit.md`** — POASTA 替换 pass-through 审计：replacement 接受仅由
  path-validity gating，**无**图质量/compression-ratio/objective 检查守卫。
- **`crush-retry-on-poor-compression.md`** — 压缩不佳时重试。
- **`crush-scaffold-mass-zero.md`** — scaffold mass zero 问题。
- **`crush-trace.md`** — 通用追踪。
- **`crush-research-brief.md`** — 研究简报。
- **`crush-crush-handoff.md`** (即 `c4-crush-handoff.md`) — C4 crush 工作交接： 记录
  branch/PR/binary commit、known-good artifact、unrelated dirty state，便于下一轮从有用 artifact
  继续。

### 8.4 C4 (HLA) 难用例（~17 个）

C4 是 HLA 复合体区域（GRCh38#0#chr6:31891045-32123783，~232kb，HPRCv2 465 paths），impg
团队用作"难测试用例"——高重复、高 CNV、SV 富集。

- **blocker 系列**（编号 01-05b）— 阻塞 C4 正确处理的疑难问题排查：
    - `c4-blocker-01-poasta-scale.md` — POASTA scale blocker。
    - `c4-blocker-02-residual-routing.md` — 残留路由 blocker。
    - `c4-blocker-03-full-rerun-scoreboard.md` — 全量重跑记分板。
    - `c4-blocker-04-stop.md` — 停止点。
    - `c4-blocker-05b-complete-traversal-aggregation.md` — 完整遍历聚合： 候选只在 cluster 范围
      union 覆盖所有 graph path 时发出，bounded homologous node expansion 仅在 range-only union
      不足时使用。
- **CMA-ES 优化器**（3 篇）— 用 Python `cma` 包做黑盒参数优化：
    - `c4-cmaes-optimizer.md` — `scripts/c4-crush-cmaes.py` 外壳： 每个候选作为普通 `impg crush`
      命令运行，`impg graph-report` 后由 wrapper 单独判分；fallback 到确定性 seeded random sampler。
      `--mode crush-only`模式验证输出 GFA 的 path 名与拼写完全一致。
    - `c4-cmaes-results.md` — 2026-05-31 实测结果（53kb C4 graph，32 线程， 一次一个 trial）。
    - `c4-cmaes-target-shape-results.md` — 把 PGGB target spectrum 加入 objective：bp-weighted node
      copy-frequency distribution、ordered frequency distance (EMD) + total variation (TV)、node
      length distribution distance、excess total segment bp。
- **专项诊断/实验**：
    - `c4-fragment-regression-suite.md` — fragment 回归套件。
    - `c4-k311-1to1-noscaffold-shorttmp-20260606.md` — k311 1to1 no-scaffold 实验快照（带日期）。
    - `c4-node-path-coverage.md` — node/path coverage 分析。
    - `c4-self-loop-repeat-normalization.md` — self-loop repeat 归一化。
    - `c4-unified-highfreq-mask.md` — 统一高频 mask。
- **诊断与修复**：
    - `diagnose-and-fix-c4-compound.md` — C4 compound induction 诊断 （before/after artifact
      对比）。
    - `diagnose-residual-two-c4.md` — 残留两簇 C4 bubble 诊断（输入 16227 seg / 22882 link / 465
      path / 392671 seg bp / 100 trivial stringy）。
    - `diagnose-residual-underaligned-c4.md` — 欠对齐诊断。
    - `fix-c4-syng-20260530.md` — C4 syng 修复（带日期）。
    - `fix-top-flubble-sweepga.md` — top-flubble SweepGA 图诱导 bug 修复： 两个具体管道失败
      （min_match_len 大于所有 exact CIGAR run；零 PAF block 仍被接受为有效 replacement）。
- **集成与扩展**：
    - `integrate-c4-local.md` — C4 local 集成。
    - `iterative-multi-level-c4.md` — 迭代多层 C4。
    - `expand-multi-bubble-c4.md` — 扩展 multi-bubble。
    - `settle-local-replacement.md` — local replacement 收敛。
    - `evaluate-low-min-match-c4.md` — 低 min-match 评估。

### 8.5 syng 后端（~8 个，含 `designs/`）— 略

该主题文档包括 `designs/syng-integration.md`（集成架构总览）、 `syng-gfa-query.md`（local
GFA 查询配方，README 显式引用）、`syng-gfa-scaffold-filtering.md`、
`syng-parallel-construction.md`（并行构建 + 6 个 sidecar）、`syng-position-query-index.md`、
`syng-to-local-graph-translation.md`、`sequence-k-syng-filter.md`、`read-syncmer-index-design.md`。
需要时直接阅读 `impg-0.4.1/docs/` 下这些文件。

### 8.6 基因分型与 infer（~7 个）

- **`genotype-architecture.md`** — `impg genotype` 架构：围绕 graph-feature evidence 而非特定图表示。
  5 步模型（选 locus → 提取 candidate haplotype →表示为 graph feature 向量 → sample 表示为
  coverage/support 向量 → score ploidy-sized 组合）。`impg genotype cos` 用 COSIGT/LikeGT-style
  cosine similarity，`cosigt` 是 `cos` 的别名。
- **`genotype-evidence-audit.md`** — `impg genotype cos` 与 `impg infer` 的 evidence 路径审计：
  两个具体后端（syng syncmer-node 后端 + graph-node 后端），用于调试"显然正确的 assignment
  得分奇怪"。
- **`genotype-gfa-backend-design.md`** — Backend-Neutral GFA 基因分型 evidence 设计：从
  syng-syncmer-node 路径扩展到 backend-neutral GFA/variation-graph 分型，feature_space =
  `gfa-segment`，render-bundle translation feature ID。
- **`genotype-impute-debug-plan.md`** — impute 调试计划。
- **`genotype-validation-suite.md`** — 验证套件。
- **`infer-design.md`** — `impg infer` 的 stitching/mosaic 设计 （**README 显式引用**）：
  跨区间/分区输出等位基因 call。
- **`pangenome-genotyping-roadmap.md`** — 泛基因组分型路线图：
  IMPG 应成为泛基因组上的灵活分型与推断系统。长期目标管线
  `panel sequences/pangenome graph → implicit graph backend → sample evidence projection → local candidate subwalks → local genotype scoring → recombination/copying inference → inferred phased haplotype mosaics`。

### 8.7 图管道、DSL 与渲染（~3 个）

- **`graph-pipeline-dsl.md`** — `-o gfa:<stage>:<stage>` DSL 设计（对应 §2.4 的 stage 解析）。
- **`graph-quality-validation.md`** — 图质量验证。
- **`render-gbz-translation-design.md`** — Render/GBZ/Translation 设计：IMPG
  是隐式泛基因组上的翻译系统，scalable object 不是单张物化图，而是在源序列坐标、graph feature、
  evidence projection、inferred haplotype 之间不丢身份地移动的能力。Root namespace 是
  `source_sequence_id : [0, source_length)`。

### 8.8 外部工具封装审计（4 个）

同行评审系列：

- **`external-tool-wrapper-audit.md`** — 主审计文件。
- **`external-tool-wrapper-audit-aligner-review.md`** — aligner 封装审计 （wfmash/FastGA/SweepGA
  等）。
- **`external-tool-wrapper-audit-render-review.md`** — render 封装审计 （gfalook/odgi 等）。
- **`external-tool-wrapper-audit-peer-synthesis.md`** — 同行综合。

### 8.9 杂项单文件

- **`audit-eliminate-synthetic-local-ids.md`** — 合成 local ID 审计。

### 8.10 `evaluations/` 子目录（~35 个文件 + 2 个嵌套子目录）

测试运行结果与评估脚本输出。**重要特点**：多数文件以 `Task: <id>` + `Evaluator: agent-<N>` +
`Date:` 开头，是 agent 系统自动评分的产物，常见 `Overall score: 0.00 / 1.00` 与
`Rubric underspecified: true/false` 字段。

#### 8.10.1 C4 系列评估（~16 个）

- **`c4-autopoietic-recovery.md`** — 评分**0.00/1.00**，confidence 0.98。 要求作为 C4
  压缩恢复过程的稳定 supervisor。
- **`c4-refine-one.md`** — 评分**0.00/1.00**。两个 SPOA-only `impg crush` refinement。
- **`c4-lower-ceiling.md`** — Rubric underspecified: true。bounded iterative-multi-level C4 crush
  变体，sPOA 用于 <2kb，auto-POASTA ceiling 2.5kb/5kb，10-15kb 完全同源窗口路由到 SweepGA。
- **`c4-full-sweepga-seed-poa1kb.md`** — 全 C4 whole-region SweepGA seed + 1kb POA crush。
- **`c4-k311-poa-threshold-sweep.md`** — k311 + POA threshold 扫描 （500/1000/2000/5000/10000 bp）。
- **`c4-low-seqwish-k-sweep.md`** — 低 seqwish-k 扫描 + POA polish。
- **`c4-micro-tangle-spoa-audit-20260604.md`** — 微 tangle（百 bp 级）SPOA 审计：检查 POVU/crush
  是否将其作为 unique bubble 收集，测试直接 SPOA 修复。
- **`c4-spoa-true-series-validate.md`** — SPOA true series 验证 （100/250/500/1k/2k/5k/10k
  threshold 全部存在于 `series.tsv`）。
- **`c4-syng-tail-diagnosis-20260605.md`** — syng 3' 尾诊断：reproducer 重现报告的 tail，分析
  `syng2gfa` frequency-filter 后只过滤 3/7015 local syncmer node。
- **`compare-aggressive-c4.md`** — 激进对比。
- **`compare-c4-cap64.md`** — cap64 对比。
- **`diagnose-c4-underalignment.md`** — C4 欠对齐诊断。
- **`explore-c4-parent.md`** — parent-first 大 crush：`window-mode=outward` 让 multi-bubble
  resolver 生成 parent-sized 残留窗口。
- **`run-c4-sweepga-abpoa10k.md`** — SweepGA seed + abPOA 10kb polish。
- **`run-full-c4.md`** — 全 C4 motif-local POA condensation。
- **`run-larger-poasta.md`** — 更大 POASTA C4 internal crush。
- **`test-poasta-insertion-order.md`** — POASTA 插入顺序对 `Poasta_crush04` (span 6460bp, med
  6461bp, cov 47/465) 的影响。
- **`validate-flank-aware-crush-c4.md`** — flank-aware crush 在 C4 上的验证 （commit `b904c10`）。
- **`validate-current-c4.md`** — 当前 C4 local graph 质量验证（HEAD `e3c2494`）：
  **结论 "better than prior SYNG-derived C4 outputs on the specific repeat-artifact failure mode, but still not solved and not yet..."**。

- **`validate-current-c4-metrics.json`** / **`validate-current-c4-metrics.tsv`** —
  对应的指标数据文件。

#### 8.10.2 通用审计与评估

- **`add-pggb-frequency.md`** — pggb 频率过滤参数添加。
- **`audit-and-validate.md`** — 审计与验证。
- **`audit-poasta-replacement.md`** — POASTA 替换循环 C4 crush 审计 （`IMPG_CRUSH_DEBUG_DIR`
  instrumentation）。
- **`evaluate-abpoa-as.md`** — abPOA 评估。
- **`evaluate-grouped-multi.md`** — grouped multi-bubble 评估。
- **`evaluate-local-compression-testbed.md`** — 测试床评估 （**blocked by missing runner outputs**，
  follow-up `follow-up-produce`）。
- **`explain-and-improve.md`** — 评分**0.00/1.00**。解释与改进。

#### 8.10.3 实施记录（implement-\* 系列，agent 评分普遍低）

- **`implement-local-compression-fixtures.md`** — 评分**0.00/1.00**，confidence 0.98。
- **`implement-local-compression-testbed-runner.md`** — 评分**0.00/1.00**， confidence 0.99。
- **`implement-multi-bubble.md`** — 评分**0.00/1.00**，confidence 0.97。
- **`implement-occurrence-level.md`** — 评分**0.20/1.00**，confidence 0.78 （唯一非零分）。
- **`local-compression-testbed-fast-synthesis.md`** — testbed 综合报告，引用 fast-profile artifacts。

#### 8.10.4 `evaluations/local-compression-autopoietic/` — 自循环迭代测试

"autopoietic"（自生系统）— 自动反馈循环：每轮运行 → 分析候选 → 综合 → 进入下一轮。

- **`iter-1.md`** — iter 1：bounded SmoothXG-style sorted chunk window 应在 `nested_top_level_wrong`
  对抗 fixture 上分裂同长度独立变体簇，同时保持 exact path spelling。
- **`iter-2-candidate-analysis.md`** — iter 1 之后的候选算法对比决策报告 （无运行，无代码改动）。
- **`iter-2-synthesis.md`** — 综合 iter-2-candidate-analysis 与 iter-1。
- **`iter-2.md`** — iter 2：加入 metric-instrumented `chunk_window_sweepga_seqwish` 候选 +
  `path_replay_compression_ratio` 诊断 metric。
- **`next-metric-blind-spots.md`** — 基于 iter-1 fast scoreboard 的 metric 盲点分析。
- **`summary.md`** — 滚动状态表（含 iter 1/2/2b 三轮的 hypothesis/change/result/next recommendation
  列）。

#### 8.10.5 `evaluations/local-compression-testbed-fast/` — 测试床 fast

profile

结构与代码测试目录类似，含数据文件：

- **`report.md`** — 主报告，含 reproduction 命令 （`scripts/local_compression_testbed.py`）。
- **`artifact-index.md`** — 产物索引：`resolve-merge-for` merge 故意省略 bulky 生成 artifact，可从
  `wg/agent-626/continue-local-compression@9464fef` 恢复。
- **`validation.md`** — 验证：merge payload note + commands。
- **`fixture-validation.json`** — fixture 验证数据，列出 12 个 validated fixtures
  （`snp_bubble_3path`、`short_indel_3path`、`mid_insertion_200bp`、`alu_like_insertion`、
  `adjacent_bubbles_joint`、`fake_repeat_anchor_split`、`microtangle_repeat_motif`、
  `duplicated_flank_context`、`tandem_copy_loop_keep`、`dispersed_repeat_glue_break`、
  `inversion_like_case`、`nested_top_level_right`）。
- **`scoreboard.json`** — 计分板：每个 fixture × method（如 `local_syng_raw`） 的 `fixture_class`/
  `tier`/`method_family`/`method_parameters`/`command_line`/`output_gfa_path` 等字段。
- **`fixtures/`** — 测试 fixture 子目录（每 fixture 一份 input.fa/expected_paths.tsv/metadata.json +
  各 method 的 output.gfa）。

### 8.11 文档命名约定与特点

通读 111 个文件后观察到的约定：

1. **前缀主题化**：文件用 `主题-子主题.md` 命名（如 `crush-aligner-deep-diag.md`），便于
   `ls crush-*` 列出同主题文档。但**前缀不总等于主题**：`c4-crush-handoff.md` 前缀是 c4，内容是
   crush 交接；`pangenome-genotyping-roadmap.md` 前缀是 pangenome，主题属于 genotype。
2. **时间戳后缀**：少数文件以 `YYYYMMDD` 结尾 （`c4-k311-1to1-noscaffold-shorttmp-20260606.md`、
   `c4-micro-tangle-spoa-audit-20260604.md`、`calibrate-spectrum-driven-20260530.md`、
   `fix-c4-syng-20260530.md`），表示某天的运行快照。绝大多数文件不带时间戳，靠 git 历史记录演进。
3. **实验 vs 设计 vs 诊断三分类**（按内容而非仅按文件名）：
    - 设计/规范：`*-design.md`、`*-spec.md`、`*-architecture.md`、 `*-quality-pass.md` — 前瞻性设计
    - 实验结果：`*-exp-*.md`、`*-sweep.md`、`*-results.md`、`run-*.md`、 `compare-*.md` — 实验报告，
      多数含 hypothesis + TL;DR + 实测数据
    - 诊断/审计：`*-diag*.md`、`*-audit.md`、`*-trace.md`、`*-failure*.md`、 `diagnose-*.md`、
      `*-blocker-*.md` — 回溯性诊断
    - 实施/评估：`implement-*.md`、`evaluate-*.md`、`validate-*.md` — agent 系统产物，常含
      0.00/1.00 评分
4. **agent 系统产物特征**：`evaluations/` 下多数文件以 `Task: <id>` + `Evaluator: agent-<N>` +
   `Date:` + `Overall score: X.XX / 1.00` +`Confidence: 0.XX` + `Rubric underspecified: true/false`
   开头，是自动化 agent 评分系统的输出。`Overall score: 0.00` 频繁出现，反映评估严苛或任务未完成。
5. **无 README/index**：`docs/` 顶层无 `README.md` 或索引文件， 靠文件名前缀自组织。新读者需要
   `ls | sort` 后按前缀浏览。
6. **README 引用**：`README.md` 仅显式引用 `docs/syng-gfa-query.md`、 `docs/genotype-architecture.md`、
   `docs/infer-design.md` 三个文件，其余 100+文档是"暗物质"——只对开发者可见，不对最终用户暴露。
7. **artifact 路径风格**：文档大量引用 `/home/erikg/impg/data/<dir>/` 与 `data/...` 路径，以及
   `https://hypervolu.me/~erik/impg/<png>` PNG 链接，反映 single-developer (Erik Garrison) +
   agent-assisted 工作流。
8. **branch/worktree 引用**：多数实验文档记录 `Branch: wg/agent-<N>/<topic>` 与 binary commit hash，
   便于结果复现。

### 8.12 与 `pgr` 文档实践的对比

`pgr` 的 `docs/` 是**面向用户的设计文档**（每格式/每模块一篇，正文英文， AGENTS.md 第 10
行约束），数量受控（~25 篇），每篇对应一个 `pgr <command>`或核心算法。`impg` 的 `docs/` 是
**面向开发者的工作笔记**（英文为主、大量实验报告、无 README 索引、含 agent 评分产物），数量爆炸
（111 篇 + evaluations / 子树），反映其活跃开发节奏与"先写文档再写代码"+"agent-assisted"的工程文化。

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

