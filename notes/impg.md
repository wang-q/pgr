# impg 分析笔记

本文档旨在总结 `impg` (implicit pangenome graph) 项目的核心设计、命令分发架构与关键数据结构，重点剖析 `src/main.rs` 的组织方式，为 `pgr` 项目在泛基因组比对、区间投影与图构建管道方面提供参考。

`impg` 由 Andrea Guarracino、Bryce Kille、Erik Garrison 开发（v0.4.1，Rust 2021 edition），其核心思想是：**不显式构建泛基因组图**，而是把 all-vs-all 两两比对当作一张"隐式图"，通过区间树 (coitrees) 在比对网络上投影目标区间，按需提取同源序列。它同时集成了 syng syncmer GBWT 后端，可在不运行比对器的情况下做查询与分型。

## 1. 项目概览

### 1.1 设计哲学: 隐式泛基因组图

传统泛基因组工具（如 pggb、Minigraph-Cactus）需要先物化一张 GFA 图，再在图上做下游分析。`impg` 走第三条路：

*   **比对即图**：把 PAF/1ALN/TPA 两两比对视为图的隐式描述 — 比对双方的坐标区间就是"边"，序列本身是"节点"。
*   **按需投影**：给定一个目标区间 `seq:start-end`，在区间树上查找所有重叠的比对，把目标坐标 lift 到查询序列上，输出 BED/BEDPE/PAF/FASTA/GFA/VCF/MAF。
*   **传递闭包**：`-x` 选项递归地把初次结果当作新的查询目标，找全所有同源片段，类似 Cactus 的 transitive alignment 但延迟到查询时才执行。
*   **不物化图**：除非用户显式要求 `gfa`/`vcf` 输出，否则永远不会落到 GFA 文件。

### 1.2 双后端架构

`impg` 支持两种对齐数据来源，通过同一个 `-a` 参数透明切换：

| 后端 | 输入 | 工作方式 | 典型用途 |
|---|---|---|---|
| **Alignment 后端** | PAF / 1ALN / TPA 文件 | 用 coitrees 索引比对区间，按 CIGAR/tracepoint 投影 | 已有 all-vs-all 比对结果 |
| **syng 后端** | syng 索引前缀（`.1khash` / `.1gbwt` / `.syng.spos` 等 sidecar） | 用 syncmer GBWT 找共享锚点，链式锚点定义同源 | 免比对，直接从 FASTA/AGC 构建 |

`main.rs` 中的 [`detect_syng_prefix`](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L4604) 函数会根据 `-a` 参数的扩展名或磁盘上是否存在 `.1khash` 文件来自动判定走哪个后端，调用者无需显式声明。

### 1.3 源代码组织

`impg` 把库逻辑放在 `src/lib.rs`，把命令行入口与参数解析放在 `src/main.rs`。`lib.rs` 导出以下模块：

*   **核心索引**：
    *   `impg.rs` — 单文件 IMPG 索引，定义 `Impg` / `QueryMetadata` / `CigarOp` / `SortedRanges`。
    *   `impg_index.rs` — `ImpgIndex` trait，统一单文件与多文件索引接口。
    *   `multi_impg.rs` — `MultiImpg`，协调多个 per-file `.impg` 子索引，负责 ID 翻译与 staleness 检测。
*   **比对格式解析**：
    *   `paf.rs`、`onealn.rs`、`tpa_parser.rs`、`alignment_record.rs` — 三种比对格式的统一抽象。
*   **syng 后端**：
    *   `syng.rs` — 核心：`SyngIndex`、`SyngMatcher`、`SyncmerParams`、`SampledPositions`。
    *   `syng_ffi.rs` — 对 C 库 `syng` 的 FFI 绑定。
    *   `syng_graph.rs` / `syng_graph_norm.rs` / `syng_parallel.rs` / `syng_transitive.rs` — syng 图构建、归一化、并行、传递闭包。
*   **命令实现**：`commands/` 目录下每个子命令一个文件 — `align.rs`、`genotype.rs`、`graph.rs`、`infer.rs`、`lace.rs`、`partition.rs`、`refine.rs`、`render.rs`、`similarity.rs`、`syng2gfa.rs`。
*   **图构建管道**：`graph.rs`、`graph_pipeline.rs`、`graph_report.rs`、`smooth.rs`、`gfa_self_loops.rs`、`render_bundle.rs`。
*   **序列索引**：`faidx.rs`、`agc_index.rs`、`sequence_index.rs`、`seqidx.rs`、`sequence_namespace.rs`。
*   **杂项**：`genotyping.rs`、`pack.rs`、`projection.rs` + `projection/converter.rs`、`subset_filter.rs`、`forest_map.rs`。

### 1.4 关键依赖

*   **`coitrees`** — cache-oblivious interval trees，区间查找的核心数据结构。
*   **`gfa` + `handlegraph`** — GFA 读写与 handlegraph 抽象（用于 lace、crush 等图操作）。
*   **`lib_wfa2`** — WFA 仿射比对，用于 BiWFA 边界精修。
*   **`spoa_rs`、`poasta`** — POA MSA 引擎，用于 `gfa:poa` 与 crush 阶段。
*   **`sweepga`、`seqwish`、`allwave`、`gfasort`、`bluntg`、`povu`** — 同一团队维护的泛基因组工具链，作为库直接嵌入（而非子进程调用）。
*   **`ragc-core`** — AGC 序列归档格式支持。
*   **`onecode`、`tpa`、`tracepoints`** — 1ALN / TPA 比对格式与 tracepoint 编解码。
*   **`rayon`、`crossbeam-channel`、`indicatif`** — 并行、通道、进度条。
*   **`noodles`、`rust-htslib`** — BIO 格式与 BAM/CRAM 支持。

## 2. main.rs — 命令分发与参数解析（重点）

`src/main.rs` 是 `impg` 二进制的入口，单文件约 **613 KB / 1.5 万行+**，承担了所有 clap 命令定义、参数解析、stage 解析、命令分发与大量辅助逻辑。这是该项目最显著的工程特点（也是值得反思之处）。

### 2.1 文件顶层结构

文件大致按以下顺序组织（行号近似）：

| 区段 | 行范围 | 内容 |
|---|---|---|
| imports 与 `parse_*` 工具函数 | 1–250 | 度量后缀解析、syng 参数解析、render 引擎解析 |
| GFA 简写解析 | 250–4486 | `apply_gfa_output_engine_shorthand` 与各种 stage 解析器 |
| `GenotypeCommand` 子枚举 | 4486–4700 | `genotype` 下的 `Cos` (cosigt) 子命令 |
| 顶层 `Args` enum | 4706–6160 | 所有 20 个子命令的 clap 定义 |
| `main` / `run` 分发 | 6161–10250 | `match args { ... }` 巨型 dispatch |
| 辅助函数（验证、输出、初始化） | 10250–end | `build_graph_config`、`validate_*`、`initialize_threads_and_log` 等 |

### 2.2 顶层 `Args` enum 与子命令清单

[main.rs#L4706](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L4706) 处的 `#[derive(Parser)] enum Args` 定义了所有子命令。下表按功能归类：

| 类别 | 子命令 | 别名 | 职责 |
|---|---|---|---|
| **索引** | `Index` | — | 构建 `.impg` 索引（single 或 per-file 模式） |
| **查询/投影** | `Query` | — | 把目标区间通过比对网络投影（核心命令） |
| | `Project` | `projection` | 投影命令 |
| | `Map` | — | 把短读段映射到 syng 索引（GAF/PAF/pack/proj 输出） |
| **分割/拼接** | `Partition` | — | 把 cohort 切成窗口化 loci |
| | `Lace` | — | 把多个 per-window GFA/VCF 拼回一张图 |
| | `Refine` | — | 收紧 loci 边界以最大化 sample/haplotype 支持度 |
| **分析** | `Similarity` | — | 区间内成对相似度/距离 + PCA/MDS |
| | `Stats` | — | 汇总比对统计 |
| **图构建** | `Graph` | — | 从 FASTA 直接构建泛基因组图（不需要预比对） |
| | `NormalizeSelfLoops` | — | 折叠 blunt GFA 中路径局部的 self-loop 重复单元 |
| | `Crush` | — | 解析 bounded bubble，支持多种 POA/POASTA/sweepga 路由 |
| | `Gfa2vcf` | `gfa-to-vcf`, `povu` | GFA → VCF |
| | `DescribeGraph` | `describe-gfa` | 输出图特征 Markdown/JSON/TSV 报告 |
| | `Render` | — | 用 gfalook 渲染 1D 图 |
| | `Align` | — | 调用 wfmash/FastGA 跑比对 |
| **syng 后端** | `Syng` | — | 从 FASTA/AGC 构建 syng 索引 |
| | `SyngRepair` | — | 重建 `.syng.pstep`/`.syng.spos` 而不重读序列 |
| **基因分型** | `Genotype` | `gt` | 基因分型命名空间，含 `Cos` (cosigt) 子命令 |
| | `Infer` | — | 跨区间/分区输出等位基因 call，支持 stitching |

`Query` 命令的 `after_help` 文档（[main.rs#L4892](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L4892)）是整个项目最完整的输出格式说明，列出了 `auto/bed/bedpe/paf/gfa/vcf/maf/fasta/fasta+paf/fasta-aln/gbwt` 共 11 种输出与各自的约束。

### 2.3 命令分发模式 (`run` 函数)

[main.rs#L6168](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L6168) 处的 `fn run() -> io::Result<()>` 是命令分发中心。它的模式如下：

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

1.  `initialize_threads_and_log(&common)` — 设置 rayon 线程池与日志级别。
2.  解析/验证参数（`validate_output_format`、`require_merge_distance`、`validate_selection_mode`、`validate_region_size` 等）。
3.  通过 `engine_cli.parse_engine()` 解析 GFA 引擎简写（如 `pggb:10000`）。
4.  路由到 `commands::<sub>::run(...)` 中的实际实现。

`main.rs` 本身**不包含**子命令的业务逻辑（除了 `Index`/`Lace` 等较简单的命令），它只负责参数装配与调用 `commands::*` 模块。这种"瘦分发 + 胖模块"的边界是项目刻意维持的。

### 2.4 辅助函数分类

`main.rs` 顶部的约 4500 行是大量辅助函数，按职责可归为几类：

#### 参数解析 (`parse_*`)

```rust
fn parse_size(s: &str) -> Result<u64, String>            // 接受 k/m/g 后缀
fn parse_merge_distance(s: &str) -> Result<i32, String>  // i32 边界检查
fn parse_usize_size(s: &str) -> Result<usize, String>
fn parse_round_count(s: &str) -> Result<usize, String>   // "until-done" → usize::MAX
```

注意 `parse_merge_distance` 内部调用 `sweepga::parse_metric_number`，作者注释说明这是为了"让 impg 和 sweepga 用同样的方式解释用户输入的距离后缀"——这是与依赖库保持一致性的设计选择。

#### syng 参数与索引检测

*   `resolve_syng_syncmer_params` — 调和 legacy `--syncmer-k`/`--syncmer-w` 与新 `--smer-length`/`--syncmer-length` 两套互斥参数。
*   `detect_syng_prefix` — 通过后缀 (`.1khash`/`.1gbwt`/`.spos`/`.pstep`/`.names`/`.meta`) 或 sibling 文件存在性判断 `-a` 是否指向 syng 索引。
*   `resolve_syng_prefix` — 单文件场景下的封装。

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

[main.rs#L164](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L164) 的 `apply_gfa_output_engine_shorthand` 拆解这些冒号分隔的 stage，并分派给一系列 `parse_*_stage` 函数：

| Stage 函数 | 处理的内容 |
|---|---|
| `parse_terminal_n_clip_stage` | `cut-n=<bp>` 终端 N-run 裁剪 |
| `parse_syng_mask_stage` | syng 高频 syncmer mask (`top`, `max-occ`, `freq-run`, `freq-span`, `sequence-k`, `min-run` 等) |
| `parse_crush_stage` | bubble crush (`method`, `k-nearest`, `k-farthest`, `pair-trees`, `polish-rounds` 等) |
| `parse_smooth_stage` | smoothxg 风格平滑 (`target-poa-length`, `max-node-length`, `block-source`) |
| `parse_graph_sort_stage` | 最终 gfasort pipeline (默认 `Ygs`) |
| `parse_syng_assertion_params` | syng `k/s/seed` 参数断言 |

每个 stage 解析器返回 `Option<String>` 表示剩余的 stage 串，串成一条管道。这种"stage 化字符串 DSL"是 `impg` 的特色，但也使 `main.rs` 体积膨胀。

#### FASTA/FASTQ 流式读取

`read_fasta_records`、`read_fastq_records`、`stream_fasta_query_chunks`、`stream_fastq_query_chunks`、`stream_query_chunks_with` 等函数为 `Map` 命令提供流式读序支持，避免一次性把所有读段加载到内存。

#### syng map 输出

`emit_syng_map`、`emit_syng_map_gaf`、`emit_syng_map_paf`、`emit_syng_map_pack`、`emit_syng_map_pack_binary`、`emit_syng_map_projection` — 对应 `impg map` 的 GAF/PAF/pack/proj 四种输出。`pack` 是按 syng node ID 索引的 `u8` 计数向量，zstd 分块压缩以支持按 node ID 随机访问；`proj` 是 `sample.pack + reads.gaf.zst + manifest.json` 的目录 bundle。

#### 进程级 stdout 静默

[main.rs#L1883](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L1883) 的 `silence_stdout_for_process`（unix 分支）通过 `RawFd` 重定向屏蔽子进程的 stdout，用于 `gfaffix`、`odgi` 等外部工具的调用。这是 unix-only 代码，非 unix 平台有空实现（[main.rs#L1926](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L1926)）。

### 2.5 测试入口

[main.rs#L6151](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L6151) 的 `fn args_command_for_test() -> clap::Command` 暴露完整的命令定义给集成测试，用于断言 help 文本与参数互斥规则。

## 3. 核心数据结构 (Impg 隐式图)

### 3.1 `Impg` 结构

[impg.rs#L394](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/impg.rs#L394) 定义了单文件索引的核心结构：

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

其中 `TreeMap = FxHashMap<u32, Arc<BasicCOITree<QueryMetadata, u32>>>` — 每个 target 序列一棵区间树，节点 metadata 是 `QueryMetadata`。`RwLock` + `Arc` 允许查询时无锁共享。

### 3.2 `QueryMetadata` 与 `CigarOp` 紧凑编码

[impg.rs#L165](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/impg.rs#L165) 的 `QueryMetadata` 用 bit-packing 压缩 metadata：

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

`CigarOp` ([impg.rs#L74](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/impg.rs#L74)) 把 CIGAR 操作 + 长度压进单个 `u32`：高 3 位是 op (`=`/`X`/`I`/`D`/`M`)，低 29 位是长度。这种设计使区间树节点极紧凑，能装下全基因组规模的 all-vs-all 比对。

`Impg` 通过 `get_cigar_ops` 按需还原 CIGAR：
*   **PAF**：直接从原文件 `data_offset` 处读取 `cg:Z:` 标签的字节。
*   **1ALN/TPA**：从 tracepoint 解码，需要 `trace_spacing` 与（可能的）目标序列做 BiWFA 还原。

### 3.3 `SortedRanges` 与区间合并

[impg.rs#L243](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/impg.rs#L243) 的 `SortedRanges` 维护按起点排序的区间集合，支持基于 `min_distance` 的合并。`insert` 方法返回"未被现有区间覆盖的新增部分"，是传递闭包查询中"只把新发现区间加入下一轮 BFS/DFS"的关键。

### 3.4 `ImpgIndex` trait 与 `MultiImpg`

[impg_index.rs](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/impg_index.rs) 定义了 `ImpgIndex` trait：

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

`Impg`（单文件）与 `MultiImpg`（多文件）都实现这个 trait。`MultiImpg` 内部维护 `TreeLocation { index_idx, local_target_id }` 把全局 `target_id` 翻译到子索引的本地 ID。`main.rs` 中的命令代码只与 `&dyn ImpgIndex` 打交道，从而对单/多文件透明。

`MultiImpg` 还实现了 staleness 检测：当 `.impg` 索引比源比对文件旧时，警告并要求 `--force-reindex`。

## 4. syng 免比对后端

### 4.1 `SyncmerParams`

[syng.rs#L1464](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/syng.rs#L1464) 定义了 syncmer 参数：

```rust
pub struct SyncmerParams {
    pub k: u32,    // inner k-mer length (s-mer), 默认 8
    pub w: u32,    // window length; total syncmer length = w + k, 默认 55
    pub seed: u32, // hash seed, 默认 7
}
```

约定 `syncmer_length = k + w` 必须为奇数，默认 63。`impg syng` 命令接受 `--smer-length` (`s`) 与 `--syncmer-length` (`k+w`) 两套参数，`resolve_syng_syncmer_params` 处理 legacy `--syncmer-k`/`--syncmer-w` 与新参数的互斥关系。

### 4.2 `SyngIndex` 与 sidecar 文件

一个 syng 索引由 6 个 sidecar 文件组成，共享同一前缀：

| 后缀 | 内容 |
|---|---|
| `.1khash` | syncmer → node ID 字典 |
| `.1gbwt` | GBWT 索引（路径 → syncmer 出现） |
| `.syng.names` | path/序列名 ↔ ID 映射 |
| `.syng.pstep` | 采样 path-step 检查点（用于坐标定位） |
| `.syng.spos` | 采样 syncmer 出现位置 |
| `.syng.meta` | 参数元数据（`syncmer_k/w/seed`） |

`SyngIndex` ([syng.rs#L2249](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/syng.rs#L2249)) 在内存中持有这些 sidecar 的句柄。`SampledCheckpointIndex` 与 `SampledPositions` 实现"GBWT 出现 → 最近的采样 checkpoint → 通过 `.spos` 解析绝对坐标"的两级定位。`--position-sample-rate 256`（默认）表示每条 path 上每 256 个 syncmer-step 采样一次，外加 path 末端的终末 syncmer。

`SyngRepair` 命令可仅从 `.1gbwt` + `.1khash` 重建 `.pstep`/`.spos`，无需重读原始 FASTA。

### 4.3 syng 查询与 `map`

`SyngMatcher` ([syng.rs#L2119](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/syng.rs#L2119)) + `SyngMatcherWorker` 实现查询算法：

1.  在 query 区间上枚举 syncmer。
2.  按 `--syng-seed-drop-top-fraction` (默认 0.05%) 丢弃最高频的 query-local seed。
3.  用 GBWT walk seed（默认 5 个 syncmer，`--syng-seed-walk-anchors`）找候选 range。
4.  对候选 range 做 BiWFA 边界精修（需要 `--sequence-files`）或 `--syng-raw` 直通。
5.  链式锚点过滤：`--syng-min-chain-anchors`（自适应 cap）、`--syng-min-chain-fraction` (默认 0.5)。

`impg map` 把短读段投影到 syng 索引，输出格式：
*   `gaf`（默认）— per-read syncmer-node walk。
*   `paf` — 投影到基因组坐标。
*   `pack` — 紧凑 `u8` node 计数向量，zstd 分块压缩，按 node ID 随机访问。
*   `pack-tsv` — pack 的人类可读 TSV。
*   `proj` — `sample.pack + reads.gaf.zst + manifest.json` 的目录 bundle，携带 read-walk 证据，供 `impg infer` 的 stitching 使用。

## 5. GFA 构建管道

`graph`、`query -o gfa`、`partition -o gfa` 三个命令共享同一套引擎实现，由 `--gfa-engine` 选择 ([lib.rs#L43](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/lib.rs#L43) 的 `GfaEngine` 枚举)：

| Engine | Pipeline | 用途 |
|---|---|---|
| `Pggb` (默认) | sweepga + seqwish + smoothxg 平滑 + gfaffix | 平滑变异图 |
| `Seqwish` | sweepga + seqwish + gfaffix | 原始（未平滑）图 |
| `Poa` | 单遍 SPOA | 小区域、快速 MSA 输出 |
| `SyngNative` | syng 锚点 + BiWFA + allwave 稀疏化 + seqwish | 近缘单倍型快速通道，跳过外部 aligner |
| `SyngLocal` | 从 query-selected 序列重建本地 syng 图 | 区域 syncmer 参数扫描 |

### 5.1 stage 化管道

通过 `-o gfa:<stage1>:<stage2>:...` 简写，用户可在引擎前后插入 stage（见 §2.4）：

```
-o gfa:cut-n=100:pggb:crush,method=allwave:sort,pipeline=Ygs
       ^^^^^^^^^^ ^^^^ ^^^^^^^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^
       终端N裁剪  引擎  bubble 解析           最终排序
```

`build_graph_config` 与 `build_engine_opts` 把这些 stage 装配成 `EngineOpts` 结构，传给 `commands::graph::run` / `commands::partition::run` / `commands::syng2gfa::run`。

### 5.2 partitioned 模式与 `lace`

引擎名后追加 `:WINDOW`（如 `pggb:10000`）即进入分区模式：

1.  把目标区域按 `WINDOW` bp 切分。
2.  每个窗口独立构建 GFA（峰值内存受控）。
3.  最终用 `impg lace` 把 per-window GFA 拼回一张图，可选 `--fill-gaps` 用参考序列填充窗口间空隙。

`lace` 同时支持 GFA 与 VCF 输入，路径名必须遵循 `NAME:START-END` 约定（最后一个 `:` 是分隔符），坐标驱动重新拼装。

## 6. 对比分析: impg vs pgr

`pgr` 已有 UCSC 体系的 Chain/Net/MAF/AXT/PSL 处理与显式 GFA 操作（参考 [docs/cactus.md](file:///Volumes/ExtHome/Scripts/pgr/docs/cactus.md) 与 [docs/gfa.md](file:///Volumes/ExtHome/Scripts/pgr/docs/gfa.md)）。`impg` 与 `pgr` 在泛基因组处理上走的是互补路线。

| 维度 | `pgr` | `impg` |
|---|---|---|
| **图模型** | 显式 Chain/Net（线性投影）+ GFA | 隐式图（比对网络），按需物化 GFA |
| **核心数据结构** | Newick 树 / PSL / Chain / GFA | coitrees + 紧凑 CIGAR delta + syng GBWT |
| **比对输入** | AXT/MAF/PSL/Chain（多为 UCSC 风格） | PAF/1ALN/TPA（wfmash/minimap2 风格） |
| **查询模式** | 按 coordinate 直接读取 | 区间投影 + 传递闭包 |
| **免比对后端** | 无 | syng syncmer GBWT |
| **基因分型** | 无 | `genotype cos` + `infer` (cosigt 模型) |
| **图构建** | 操作已有 GFA | 内嵌 sweepga/seqwish/allwave 等完整管道 |
| **CLI 风格** | `pgr <format> <subcommand>` 多级 | `impg <command>` 单级，但每命令参数众多 |
| **代码组织** | `cmd_pgr/` 下按格式分组，命令实现独立 | `commands/` 下按命令分文件，但 `main.rs` 单文件巨大 |

### 深度对比：区间投影 vs Chain lift

`pgr chain lift` 通过 Chain 的坐标映射把目标区间 lift 到查询序列，本质上是单条 Chain 的线性投影。`impg query` 在 all-vs-all 比对网络上做区间树查找 + 传递闭包，等价于"在所有 Chain 的并集上做 BFS"。

*   **优势**：impg 自动发现所有同源片段（包括间接通过第三序列的同源），pgr chain lift 需要用户手动选 Chain。
*   **代价**：impg 需要 all-vs-all 比对（O(n²) 内存与时间），pgr chain lift 直接用现成的 UCSC Chain（已经经过 syntenic 净化）。

### 深度对比：syng vs pgr 的距离/聚类

`pgr dist` + `pgr clust` 走"显式计算距离矩阵 → 聚类"的路线。`impg similarity` 在比对网络或 syng 索引上直接出成对相似度，并支持 PCA/MDS。`impg partition` 的 `--selection-mode sample|haplotype`（PanSN 命名）展示了如何用样本/单倍型分组而非纯序列长度来驱动分区。

## 7. 对 `pgr` 的启示

1.  **区间树 + 紧凑 CIGAR delta 的组合值得借鉴**：`pgr` 处理 PAF 时可直接复用 `coitrees` + `CigarOp` 风格的 bit-packing，把全基因组 all-vs-all 比对装进内存。当前 `pgr` 的 PSL/Chain 处理是流式的，缺少随机访问能力。

2.  **trait 抽象单/多文件索引**：`ImpgIndex` trait + `MultiImpg` 是处理"单大文件 vs 多小文件"两种部署模式的干净做法。`pgr` 若引入类似的索引层，可让命令代码与索引物理形态解耦。

3.  **stage 化字符串 DSL 的两面性**：impg 的 `-o gfa:cut-n=100:pggb:crush:sort` 简写表达力强，但代价是 `main.rs` 膨胀到 60 万字符。`pgr` 当前坚持 `pgr <format> <subcommand>` 的多级结构更易维护，应保持。若未来需要类似的管道组合，可考虑专门的 pipeline 配置文件而非 CLI 简写。

4.  **避免 main.rs 巨型化**：impg 把 20 个子命令的 clap 定义与分发全塞进单文件是明显的反例。`pgr` 的 `src/cmd_pgr/` 按格式/功能分组、每命令独立模块的结构更优，应继续坚持——`main.rs` 只做 `ArgMatches` 分发，业务逻辑下沉到模块。

5.  **syng 免比对后端的思路**：syncmer GBWT + 采样 checkpoint 的两级坐标定位，是处理"大规模泛基因组快速查询"的有效方案。`pgr` 当前聚焦显式格式处理，若未来要做泛基因组快速查询，syng 的 sidecar 设计（6 个文件分工明确、`.meta` 自描述参数）是值得参考的工程模板。

6.  **thread-local 缓存模式**：impg 用 `thread_local!` 缓存 WFA aligner、1aln/TPA 句柄、目标序列片段，避免重复分配。`pgr` 在并行处理 PSL/Chain 时可借鉴同样的模式。

7.  **PAF `cg:Z:` 懒加载**：impg 不把 CIGAR 存进区间树节点，只存 `data_offset` + `data_bytes`，查询时按需读取。这对 `pgr` 处理大型 PAF 是直接可借鉴的内存优化。

8.  **PanSN 命名约定**：impg 全程用 `sample#haplotype#contig` 命名（`#` 分隔），`pgr` 在 `pgr pl` 流水线中若需要处理群体数据，可采用同一约定以与 pggb/impg/odgi 生态兼容。

9.  **Zero Panic 与 AGENTS.md 的契合**：impg 源码中存在大量 `unwrap_or_else(|e| panic!(...))`（如 `get_cigar_ops`、`get_target_sequence_cached`），违反了 `pgr` 的 Zero Panic 原则。`pgr` 在借鉴其算法时应改为 `anyhow::Result` + `bail!`，把错误返回到调用方而非 panic。
