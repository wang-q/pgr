# PAF 模块实现参考

本文档基于 impg-0.4.1 源码，梳理 pgr 需要实现的 PAF 相关组件的代码级设计。 这是纯实现参考，
不涉及路线决策（见 [[paf-route.md]]）或第一步行动计划（见 [[pairwise-selection.md]]）。

参考源码：`paf.rs`（417 行）、`alignment_record.rs`（138 行）、`seqidx.rs`（56 行）、 `main.rs` 的
`output_results_paf` 函数（`main.rs:11989-12101`）。

---
## 1. 模块结构

```
src/libs/paf/
├── mod.rs          # 模块导出
├── record.rs       # PafRecord struct — PAF 行在内存中的表示
├── parser.rs       # PAF 解析 — 纯文本 / BGZF / GZI 三种模式
├── cigar.rs        # CigarOp bit-packing + 字符串互转
├── writer.rs       # PAF 行格式化输出
├── lazy.rs         # CIGAR 懒加载（从源文件按偏移量读取）
└── index.rs        # PafIndex — 区间树索引 + 查询 + 传递闭包 BFS

src/libs/seqidx.rs  # SequenceIndex — 序列名↔ID 双向映射（paf 的前置依赖）
```

---
## 2. PafRecord — 核心数据结构

参考 impg 的 `AlignmentRecord`（`alignment_record.rs:12-21`）。

```rust
#[derive(Debug, Clone)]
pub struct PafRecord {
    pub query_id: u32,             // col 1: query sequence ID (from SequenceIndex)
    pub query_start: u32,          // col 3: query start (0-based)
    pub query_end: u32,            // col 4: query end
    pub target_id: u32,            // col 6: target sequence ID
    pub target_start: u32,         // col 8: target start
    pub target_end: u32,           // col 9: target end
    pub strand_and_offset: u64,    // MSB=strand, bits[62:0]=file offset to cg:Z:
    pub cigar_bytes: u16,          // byte length of CIGAR string (0 if absent)
    pub matches: u32,              // col 10: matching bases
    pub block_len: u32,            // col 11: alignment block length
    pub mapq: u8,                  // col 12: mapping quality (0-255)
}
// Total: 48 bytes
```

### 设计要点

**strand 编码在 MSB**（impg `alignment_record.rs:34`）。
`const STRAND_BIT: u64 = 0x8000000000000000`——MSB=0 为 Forward，MSB=1 为 Reverse。`strand()`/
`set_strand()` 封装位操作。

**CIGAR 懒加载**。`strand_and_offset` 的低 63 位是源 PAF 文件中 `cg:Z:` tag 的字节偏移量，查询时通过
`read_cigar_data()` 按需读取。区间树节点只存坐标和指针，不存 CIGAR 字符串——这是把全基因组 PAF
装入内存的关键优化（[[impg.md]] §3.2）。

**序列 ID 化**。`query_id`/`target_id` 是 `SequenceIndex` 中的 `u32`，区间树用整数 key 而非字符串，
大幅减少内存和比较开销。

**字段全部 `u32` 而非 `usize`**，跨平台一致，更紧凑。超过 4Gbp 的染色体可升级 `u64`。

### 与 impg 的差异

impg 的 `AlignmentRecord` 只有 8 个字段——`matches`/`block_len`/`mapq` 不存，需要时从原文件重读。
这是为了把区间树节点压到最小（cache 命中率优先）。pgr 第一步存 12 个完整字段，先用起来；
后续内存成瓶颈时再瘦身到 8 字段，内部改动不影响下游。

---
## 3. PAF 解析器

impg 的 `paf.rs` 支持三种读取模式。

### 模式 1：纯文本

```rust
pub fn parse_paf<R: BufRead>(
    reader: R, seq_index: &mut SequenceIndex,
) -> Result<Vec<PafRecord>, PafParseError> {
    let mut bytes_read: u64 = 0;
    let mut records = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let record = parse_paf_line(&line, bytes_read, seq_index)?;
        records.push(record);
        bytes_read += (line.len() + 1) as u64;
    }
    Ok(records)
}
```

`bytes_read` 跟踪每个 CIGAR tag 的文件偏移量，`SequenceIndex::get_or_insert_id` 在解析过程中动态构建
name→id 映射。`parse_paf_line`（impg `paf.rs:118-177`）做 12 列 tab split → 数字解析 → strand 判断
→ name→id 查找 → CIGAR tag 定位（`paf.rs:150-161`）。

### 模式 2：BGZF 压缩

与纯文本的关键差异：CIGAR 偏移量用**BGZF virtual position**而非字节偏移量。每行先记
`line_start_vpos`，读到 `cg:Z:` 后再记 `cigar_vpos`，seek 回行首再 forward scan。双 seek
模式确保跨越 BGZF 块边界时 CIGAR 定位正确（impg `paf.rs:199-270`）。

### 模式 3：BGZF + GZI 多线程

两遍扫描：先 `parse_paf` 得到 uncompressed byte offsets，再用 `gzi_index.query()` 转为
virtual positions。好处是多线程 BGZF 解压（`MultithreadedReader`），比单线程快 3-5x（impg
`paf.rs:274-302`）。

### 统一入口

`parse_paf_file`（impg `paf.rs:306-362`）做自动检测：扩展名 `.gz`/`.bgz` → 读 18 字节头判定 BGZF
（`is_bgzf`，`paf.rs:50-66`）→ 有 `.gzi` 走模式 3，无 `.gzi` 走模式 2，否则模式 1。普通 gzip
撞上检测会报错："Convert with: zcat file.paf.gz | bgzip > output.paf.gz"。

---
## 4. CIGAR 编解码

### CigarOp bit-packing

参考 impg `impg.rs:74`。3 位 op code + 29 位 length 压入单个 `u32`：

```rust
#[derive(Debug, Clone, Copy)]
pub struct CigarOp(u32);

impl CigarOp {
    pub fn new(op: char, len: u32) -> Self {
        let code = match op { '=' => 0, 'X' => 1, 'I' => 2, 'D' => 3, _ => 4 };
        CigarOp((len & 0x1FFF_FFFF) | (code << 29))
    }
    pub fn op(self) -> char { /* bits[31:29] → char */ }
    pub fn len(self) -> u32 { /* bits[28:0] → u32 */ }
}
```

支持 `=`、`X`、`I`、`D`、`M` 五种 op。单段最长 512Mbp（29 bits），基因组比对足够。

### 字符串互转

`parse_cigar("10=5X2I3D") → Vec<CigarOp>`，`format_cigar(ops) → "10=5X2I3D"`。

### CIGAR 懒加载

```rust
pub fn read_cigar_data(
    paf_path: &str, offset: u64, byte_len: usize,
) -> Result<Vec<CigarOp>, CigarReadError>;
```

根据扩展名 dispatch：纯文本走 `File::seek(Start(offset))`，BGZF 走
`bgzf::Reader::seek(virtual_position)`。读出来的字节就是 `cg:Z:` 的原始值，再调 `parse_cigar` 解析。
impg 内部维护 `thread_local!` CIGAR cache 避免重复 seek（[[impg.md]] §9.5），pgr 后续可借鉴。

---
## 5. SequenceIndex — 序列名↔ID 双向映射

参考 impg `seqidx.rs:1-56`。

```rust
#[derive(Clone, Debug, Default)]
pub struct SequenceIndex {
    name_to_id: FxHashMap<String, u32>,
    id_to_name: FxHashMap<u32, String>,
    id_to_len: FxHashMap<u32, u32>,
    next_id: u32,
}
```

六个方法：`get_or_insert_id(name, length)`（解析时动态构建）、`get_id`/`get_name`/`get_len`（O(1)
查找）、`len`/`is_empty`。

用 `FxHashMap` 而非标准 `HashMap`——FxHash 比 SipHash 快 2-3 倍，且 PAF 解析不受 HashDoS 威胁（key
来自自己的数据文件）。pgr 已有 `fxhash` 依赖，不需要新增。

---
## 6. PAF 输出格式化

参考 impg `main.rs:11989-12101`（`output_results_paf`）。

```rust
pub fn write_paf_record(
    out: &mut dyn Write,
    query_name: &str,  query_len: u32,  query_start: u32,  query_end: u32,
    target_name: &str, target_len: u32, target_start: u32, target_end: u32,
    strand: char, matches: u32, block_len: u32, mapq: u8,
    cigar: &[CigarOp], extra_tags: &[(&str, &str)],
) -> io::Result<()> {
    let cigar_str = format_cigar(cigar);
    write!(out, "{}\t{}\t...\t{}\t{}\t255", query_name, query_len, matches, block_len)?;
    let gi = gap_compressed_identity(cigar);
    let bi = block_identity(cigar);
    write!(out, "\tgi:f:{:.6}\tbi:f:{:.6}\tcg:Z:{}", gi, bi, cigar_str)?;
    for (key, val) in extra_tags { write!(out, "\t{}:{}", key, val)?; }
    writeln!(out)
}
```

四个标准标签：`gi:f:`（gap-compressed identity）、`bi:f:`（block identity）、`cg:Z:`（CIGAR string，
PAF 标准标签）、`an:Z:`（alignment name，可选）。

identity 计算（impg `main.rs:12042-12056`）：对 CIGAR 做 fold 统计——`=` 和 `X` 分别计入
matches/mismatches，`I` 和 `D` 按**事件**计数（gi）或按**碱基**计数（bi）。`gi` 评估同源性（对长
indel 宽容），`bi` 评估序列一致性（对长 indel 严格）。

---
## 7. 错误类型

```rust
#[derive(Debug)]
pub enum PafParseError {
    NotEnoughFields,               // < 12 tab-separated columns
    InvalidInteger(ParseIntError),
    InvalidStrand,
    InvalidCigarFormat,
    InvalidFormat(String),
    IoError(io::Error),
}
```

参考 impg `paf.rs:14-37`。实现 `Display` + `Error`，不做 `From` 自动转换。在 `execute` 中用
`.map_err(|e| anyhow!("PAF: {}", e))` 桥接到 `anyhow::Result`。

---
## 8. 实现优先级

### 第一期 — 支撑 `pgr maf to-paf` + `pgr paf index` + `pgr paf query`

| 优先级 | 组件                                               | 文件                 |
|--------|----------------------------------------------------|----------------------|
| P0     | `PafRecord` struct                                 | `libs/paf/record.rs` |
| P0     | `SequenceIndex`                                    | `libs/seqidx.rs`     |
| P0     | `parse_paf_line` + `parse_paf`（纯文本）           | `libs/paf/parser.rs` |
| P0     | `parse_cigar` + `format_cigar`                     | `libs/paf/cigar.rs`  |
| P0     | `write_paf_record`                                 | `libs/paf/writer.rs` |
| P0     | `PafIndex`：`build` + `query` + `query_transitive` | `libs/paf/index.rs`  |
| P1     | `CigarOp` bit-packing                              | `libs/paf/cigar.rs`  |

### 第二期 — 大 cohort 场景

| 优先级 | 组件                         | 说明                             |
|--------|------------------------------|----------------------------------|
| P1     | `parse_paf_bgzf`             | 模式 2，支持 `.paf.gz`           |
| P1     | `read_cigar_data`（懒加载）  | 索引层不存 CIGAR，查询时按需读取 |
| P2     | `parse_paf_file`（自动检测） | 统一入口，dispatch 三种模式      |
| P3     | `parse_paf_bgzf_with_gzi`    | 模式 3，多线程 BGZF 解压         |

### 第三期 — 查询层增强

| 优先级 | 组件                       | 说明                                               |
|--------|----------------------------|----------------------------------------------------|
| P2     | `PafRecord` 瘦身（8 字段） | 去掉 matches/block_len/mapq，回退到懒加载          |
| P2     | Caf 后处理过滤参数         | `--min-degree`、`--min-chain-length`、`--end-trim` |
| P3     | Chain/Net syntenic 过滤器  | `--syntenic-filter` 参数                           |

---
## 9. 新增依赖

唯一需新增的依赖是 `coitrees`（区间树）。pgr 已有 noodles-bgzf、fxhash、rayon 等所有其他基础设施。

如果不引入 `coitrees`，可选 `intspan::RangeMap`（pgr 已有，但语义是区间集合不是区间树）或
`BTreeMap` 手写（零依赖，但需要自己实现重叠查询）。第一期建议直接引入 `coitrees`。

---
## 10. 与 pgr 现有 IO 层的关系

pgr 的 `src/libs/loc.rs`（202 行）是 FASTA 随机访问索引模块，与 impg 的 `paf.rs` 在架构上高度平行。
两者都提供多格式输入抽象 + 偏移量 seek+read + BGZF 支持。

**`Input` enum 可以直接复用**。pgr 的 `Input { Buf, File, Bgzf }` 比 impg 的
`PafHandle { Plain, Compressed }` 更通用——多了 `Buf` 支持 stdin，且 `Bgzf` 变体用的是
`IndexedReader`（自带索引，seek 无需外部 `.gzi` 文件）。

**`read_offset()` 可以直接替代 impg 的 `read_cigar_data()`**。同样是 seek+read 返回字节，
pgr 的实现更简洁（11 行 match + 2 行 I/O vs impg 的 46 行分支）。CIGAR 懒加载可以直接写成对
`read_offset` 的薄封装：读字节 → UTF-8 解析 → `parse_cigar`。

**pgr 的 `IndexedReader` 比 impg 的基础 `Reader` 更强**。impg 的 `parse_paf_bgzf_with_gzi` 需要外部
`.gzi` 索引文件 + 显式 `VirtualPosition::from(offset)` 转换；pgr 的 `IndexedReader` 在内部处理
vpos，调用者只需传字节偏移量。这意味着 pgr 的 BGZF PAF 支持可以跳过 impg 的模式 3。

**需要小幅增强 `Input`**：加一个 `read_line` 方法统一三种变体的行读取（目前只在 `create_loc`
内部匹配），以及在 `Bgzf` 变体上暴露 `virtual_position()`（CIGAR 懒加载需要记 vpos）。

**完全缺失、必须新增的三样**：区间树（`coitrees`）、PAF 行解析、CIGAR 编解码。其余都是复用或薄封装。
汇总：

| loc.rs 组件               | PAF 角色                               |        状态        |
|---------------------------|----------------------------------------|:------------------:|
| `Input` enum              | 多格式输入抽象（plain/BGZF/stdin）     |      直接复用      |
| `read_offset()`           | CIGAR 懒加载的 seek+read 基础          |       薄封装       |
| `IndexedReader`           | BGZF seek（自带索引，无需外部 .gzi）   |      直接复用      |
| `reader_buf()`            | 纯文本 PAF 的 BufRead 创建             |      直接复用      |
| `Input::read_line`        | BGZF 行迭代读取                        |     需新增方法     |
| `virtual_position()` 暴露 | CIGAR 懒加载的 vpos 记录               |       需暴露       |
| 区间树                    | "chr1:1000-5000 与哪些 PAF 记录重叠？" | 需新增（coitrees） |
| PAF 行解析                | 12 列 tab split → PafRecord            |       需新增       |
| CIGAR 编解码              | CigarOp bit-packing + identity 计算    |       需新增       |

---
## 11. PAF 索引设计

本节基于 impg 源码 `impg.rs`（3214 行）、`impg_index.rs`（392 行）、`multi_impg.rs`（1093 行）
的详细分析。

### 11.1 impg 怎么做

**`Impg` struct 的 7 个字段**。`trees: RwLock<FxHashMap<u32, Arc<BasicCoitree<QueryMetadata>>>>`
是核心——`RwLock` 支持多线程并发读取，`Arc` 让同一棵树跨线程共享。`forest_map` 记录每个 target
的树在磁盘文件中偏移量（disk-backed 模式时按需加载）。`alignment_files` 存源 PAF 文件路径（CIGAR
懒加载时 seek 用）。`trace_spacing_cache` 是 1ALN/TPA 专用，pgr 不需要。

**构建流程**（`from_multi_alignment_records`，`impg.rs:1531-1648`）：输入是
`Vec<(Vec<AlignmentRecord>, String)>`——每个文件一组 records + 源文件路径。每条 record 被
flat_map 成 1-2 个 `Interval<QueryMetadata>` entry（bidirectional 模式下生成正反两条，
索引大小翻倍），然后按 `target_id` 分组，每组用 rayon 并行构建一棵 `BasicCoitree`，collect 成
`TreeMap`。双向索引的好处是"query_A 在 target_B 树上"和"query_B 在 target_A 树上"各有一条记录，
查询时不需要额外计算。代价是索引大小翻倍。pgr V1 不做双向。

**查询**（`query`，`impg.rs:1848-1924`）：拿 `target_id` → `get_or_load_tree` →
`tree.query(range, callback)`。回调对每个重叠 interval 做坐标投影（`project_overlapping_interval`），
产出 query 侧坐标 + CIGAR。结果集的第一个元素固定是查询区间自身（identity CIGAR）。CIGAR 只在
`store_cigar=true` 时计算，否则 `AdjustedInterval` 中 CIGAR 为空 Vec。提供两种投影模式：normal
（完整 CIGAR + 序列 I/O）和 approximate（快速无 I/O）。

**CIGAR 懒加载**（`get_cigar_ops`，`impg.rs:495-548`）：从 `alignment_files` 拿到源文件路径
→ PAF 格式走 `read_cigar_data`（seek 到 `cg:Z:` 偏移量 → 读字节 → `parse_cigar_to_delta`）；
1ALN/TPA 格式走 tracepoint 解码 + BiWFA（pgr 不需要）。还有 `populate_cigar_cache` 批量预加载 +
`query_with_cache` 复用——这个设计对 pgr V1 不适用，因为 V1 的 PafRecord 直接存 12 完整字段。

**多文件索引**（`MultiImpg`，`multi_impg.rs:108-133`）：核心是
`local_to_unified: Vec<Vec<u32>>`——每个子索引的 local target_id 到全局 unified_id 的翻译表。
`forest_map: FxHashMap<u32, Vec<TreeLocation>>` 记录一个 unified_id 可能在多个子索引中有树。
加载时只读 header（seq_index + forest_map），树数据 `get_or_load_tree` 时按需 lazy load。有
staleness 检测——`.impg` 比源 PAF 旧时警告。

**`ImpgIndex` trait**（`impg_index.rs:21-121`）：13 个方法签名统一了 `Impg` 和 `MultiImpg`
的接口——`seq_index`、`query`、`query_with_cache`、`populate_cigar_cache`、`query_transitive_dfs`、
`query_transitive_bfs`、`get_or_load_tree`、`target_ids`、`remove_cached_tree`、`num_targets`、
`sequence_files`。还有一个 `ImpgWrapper` enum 做 match dispatch。

### 11.2 pgr 怎么做

**数据结构**。pgr 的 `PafIndex` 比 impg 的 `Impg` 简单得多：

```rust
pub struct PafIndex {
    pub seq_index: SequenceIndex,           // name ↔ u32
    trees: FxHashMap<u32, Coitree<PafRecord>>,  // target_id → tree (owned, not Arc)
    source_files: Vec<PathBuf>,             // for CIGAR lazy loading
}
```

对比 impg，pgr 的数据结构有六处简化：

- **`trees` 不用 `RwLock` 包裹**。impg 用 `RwLock<TreeMap>` 是为了支持构建时独占写入、
  查询时共享读取的多线程模式。pgr V1 单线程构建+查询，不需要锁。V2 rayon 化时再包一层
  `Arc<RwLock<>>`。
- **树的持有不用 `Arc`**。impg 的每棵 `Coitree` 包在 `Arc` 里，因为 `disk_backed`
  模式下同一棵树可能被多个查询线程共享引用。pgr V1 树直接由 `FxHashMap` 拥有，无共享场景。
- **不需要 `ForestMap`**。impg 的 `ForestMap` 记录每个 target 的区间树在磁盘文件中的字节偏移量，用于
  disk-backed 按需加载（`load_tree_from_disk`，`impg.rs:1720`）。pgr V1 纯内存，不序列化到磁盘。V2
  加持久化时用 bincode 整体序列化 `PafIndex`，比 impg 的 per-tree 偏移量方案更简单。
- **不需要 `trace_spacing_cache`**。1ALN/TPA 格式专用，pgr 只支持 PAF。
- **不做 bidirectional 双索引**。impg 的 `from_multi_alignment_records` 对每条 A→B record
  额外生成一条 B→A entry，索引大小翻倍但免去反向查询的额外计算。pgr V1 不做——查询时若需要反向，
  单独扫一次即可。V3 评估后再决定是否加双向。
- **格式只支持 PAF**。impg 的 `Impg` 同时支持 PAF、1ALN、TPA 三种比对格式——`get_cigar_ops` 中 PAF
  走 `read_cigar_data`、1ALN/TPA 走 tracepoint 解码 + BiWFA 重建、`get_trace_spacing` 管理 1ALN 的
  lazy spacing 缓存。pgr 只需要 PAF，三种 code path 全部砍掉。

**构建流程**。输入和 impg 一样：`Vec<(Vec<PafRecord>, &str)>`——每个文件一组 records + 源文件路径。
遍历所有 records，按 `target_id` 分组为 `FxHashMap<u32, Vec<Interval<PafRecord>>>`，每组调用
`Coitree::new(&intervals)` 建一棵树。与 impg 构建流程有三点差异：

- `flat_map` 遍历用单线程 `for` loop，不用 rayon `par_iter`。V1 先跑通，rayon 化是 P1 优化。
- 只生成正向 entry（不做 bidirectional），每条 record 对应一个 `Interval`。
- 树构建用 `into_iter` 串行，而非 impg 的 `into_par_iter` 并行。同上，P1 优化。

**单跳查询**。`query(target_id, range_start, range_end)` 先 `trees.get(&target_id)`，
命中后 `tree.query_intersecting(range, callback)`。回调对每个重叠 interval 从 `PafRecord`
的坐标和查询区间重叠量计算 query 侧坐标。与 impg 查询有四点差异：

- **不做 self-entry**。impg 结果的第一个元素是查询区间自身（identity CIGAR）。pgr
  调用者知道查询区间，不需要冗余。
- **不区分 store_cigar/not**。impg 的 `AdjustedInterval` 中 CIGAR 可为空 Vec（`store_cigar=false`
  时）。pgr V1 的 PafRecord 存 12 完整字段，数据直接可用。
- **不区分 normal/approximate 两种投影模式**。impg 的 approximate 模式跳过序列 I/O 做快速投影。pgr
  V1 默认不读序列，近似模式是默认行为。
- **不做 CIGAR cache**。impg 有 `populate_cigar_cache` 批量预加载 + `query_with_cache` 复用缓存。
  pgr V1 数据已在 PafRecord 中，不需要。

**传递闭包 BFS**。`query_transitive(target_id, range, max_depth, merge_distance)` 只实现 BFS，不做
DFS。理由：BFS 按深度分层，语义直观，且 DFS 在 impg 中主要用于 partition（需要 `masked_regions`
预填充），pgr V1 不涉及。用 `FxHashMap<u32, SortedRanges>` 去重（参考 impg 的 `SortedRanges`），
每轮只把"未被已有区间覆盖的新增部分"入队。

**与 impg `ImpgIndex` trait 的对应**。pgr 不定义 trait——V1 只有一种 `PafIndex`（没有
`MultiPafIndex`），不需要 trait 抽象。V2 多文件索引时参考 impg 加 `PafIndex` trait + 两种 impl。

### 11.3 分级路线

**V1（最小原型，单文件全内存）**：`PafIndex` 纯内存，单文件 PAF → `Coitree` per target，单跳查询 +
BFS 传递闭包，`PafRecord` 12 完整字段（不懒加载 CIGAR），不持久化索引。新增依赖只有 `coitrees`。
交付物：`pgr maf to-paf` + `pgr paf index` + `pgr paf query`。

**V2（多文件 + 持久化 + 懒加载）**：多文件统一索引用 `local_to_unified` 翻译表（借鉴 MultiImpg），
`.paf.idx` 持久化（bincode 整体序列化），CIGAR 懒加载切换为 `libs/loc.rs::read_offset` 薄封装，
rayon 并行构建 + `Arc<Coitree>` 共享，`PafRecord` 瘦身到 8 字段。

**V3（大 cohort 优化，按需开启）**：bidirectional 双索引，`remove_cached_tree` 显式内存管理（大
cohort 时按 target 驱逐），staleness 检测（索引 ↔ PAF 源文件 mtime 比较），Caf 后处理过滤参数，
Chain/Net syntenic 过滤器。

### 11.4 impg 中 pgr 可以跳过的

以下九项加起来，pgr 比 impg 少约 40% 的索引相关代码量。

| impg 特性                   | 代码位置        | pgr 决策 | 理由                                           |
|-----------------------------|-----------------|:--------:|------------------------------------------------|
| `ForestMap` 磁盘偏移量      | `impg.rs:1631`  | V1 跳过  | V2 用 bincode 整体序列化，不走 per-tree 偏移量 |
| `trace_spacing_cache`       | `impg.rs:403`   | 永远跳过 | 1ALN/TPA 专用                                  |
| `bidirectional` 双索引      | `impg.rs:1543`  | V1 跳过  | 索引大小翻倍。V3 按需评估                      |
| `approximate_mode` 分支     | `impg.rs:1896`  | 功能融合 | pgr V1 默认不读序列，= impg 近似模式           |
| `RwLock<TreeMap>` 并发      | `impg.rs:395`   | V1 跳过  | 单线程。V2 加 `Arc<RwLock<>>`                  |
| `populate_cigar_cache`      | `impg.rs:1926`  | V1 跳过  | PafRecord 12 字段，数据已在内存                |
| `with_aligner` WFA/BiWFA    | `impg.rs:53`    | 永远跳过 | tracepoint→CIGAR 重建，pgr 不需要              |
| staleness 检测              | `multi_impg.rs` | V1 跳过  | V2 加 mtime 比较                               |
| `ImpgIndex` trait + wrapper | `impg_index.rs` | V1 跳过  | 单实现，V2 多文件时按需引入                    |

## 12. impg 源码关键位置速查

以下位置经精读 impg-0.4.1 源码核验，供实现时快速跳转对照。

| 组件 | 文件 | 行号 | 说明 |
|------|------|------|------|
| `CigarOp` bit-packing | `src/impg.rs` | 73 | 3 位 op + 29 位 len → u32 |
| `QueryMetadata` | `src/impg.rs` | 165 | 区间树节点元数据（8 字段，38 字节） |
| `TreeMap` | `src/impg.rs` | 226 | `FxHashMap<u32, Arc<BasicCOITree<...>>>` |
| `SortedRanges` | `src/impg.rs` | 242 | BFS 去重核心 |
| `Impg` struct | `src/impg.rs` | 394 | 7 字段 |
| CIGAR 懒加载 | `src/impg.rs` | 494 | `get_cigar_ops()` — PAF/1ALN/TPA 三路 dispatch |
| 坐标投影 | `src/impg.rs` | 1100 | `project_overlapping_interval()` |
| 索引构建 | `src/impg.rs` | 1549 | `from_multi_alignment_records()` — rayon 并行 |
| 单跳查询 | `src/impg.rs` | 1848 | `Impg::query()` |
| 传递闭包 BFS | `src/impg.rs` | 2291 | `Impg::query_transitive_bfs()` — 核心参考 |
| `ImpgIndex` trait | `src/impg.rs` | 2584 | 单/多文件统一抽象 |
| CIGAR 懒加载 seek | `src/paf.rs` | 68 | `read_cigar_data()` |
| PAF 行解析 | `src/paf.rs` | 118 | `parse_paf_line()` — 含 CIGAR 偏移量计算 |
| 统一解析入口 | `src/paf.rs` | 306 | `parse_paf_file()` — 3 模式自动检测 |
| 索引初始化 | `src/main.rs` | 11043 | `initialize_index()` |
| 查询执行 | `src/main.rs` | 11605 | `perform_query()` — 传递/非传递 dispatch |

## 13. wgatools 参考

[wgatools](https://github.com/wjwei-handsome/wgatools) (v1.1.0, Bioinformatics 2025)
是另一个 Rust 实现的全基因组比对工具集，支持 MAF/PAF/Chain/SAM 互转。
其 PAF/CIGAR 处理与 pgr 有高度重叠，以下是可借鉴的设计。

### 13.1 与 impg 的关键差异

| 维度 | impg | wgatools |
|------|------|----------|
| CIGAR 存储 | bit-packed `u32`（紧凑，4B/op） | `Cigar` struct（含 inv_* 倒位字段，信息丰富但内存更大） |
| PAF 解析 | 手写 `split('\t')` + 字段解析 | `csv::Reader` flexible 模式（更健壮，自动处理可变列数） |
| CIGAR 来源 | `cg:Z:` tag（PAF）或 tracepoint 解码（1ALN/TPA） | `cg:Z:` 或 `cs:Z:`（自动转换） |
| 格式统一 | `ImpgIndex` trait（索引层） | `AlignRecord` trait（记录层，统一 PAF/MAF/Chain 的字段访问） |
| MAF→PAF | 无（impg 直接从 PAF 启动） | `MAFRecord::convert2paf()` — 含 `query_name` 参数选 query |

### 13.2 pgr 可借鉴的设计

**`csv` crate flexible reader**（`parser/paf.rs:22-30`）：

wgatools 用 `csv::ReaderBuilder::flexible(true)` 解析 PAF，好处是自动处理
可变列数（12 列 + 可选 tags），同时过滤 `#` 注释行。pgr 手写 `split('\t')`
解析也可以，但如果后续需要支持 BGZF 流式读取，csv crate 的 Reader 更稳定。

**`AlignRecord` trait**（`parser/common.rs`）：

统一 PAF/MAF/Chain 三种格式的字段访问接口（`query_name`、`query_start`、
`target_end` 等）。pgr 目前只有 PAF，但后续做 Chain↔PAF 互转时需要类似抽象。
可以先不做 trait，但在 `PafRecord` 上预留类似方法签名。

**PAF 校验 `validate`**（`tools/validate.rs`）：

用 CIGAR 校验 `query_end` 和 `target_end` 是否与 CIGAR 推导值一致，
不一致时自动修正。pgr 的 `to-paf` 输出可以作为输入再验证一次，
或提供给用户做数据质量检查。

**`parse_maf_seq_to_trim`**（`parser/cigar.rs:155-199`）：

从 MAF 对齐串分析首尾 indel（用于裁剪链的块边界）。pgr 在做
MAF→Chain 转换时可以参考（后续）。

**`cs:Z:` → CIGAR 转换**（`parser/paf.rs:159-200`）：

部分工具（如 minimap2）输出 `cs:Z:` 而非 `cg:Z:`。wgatools 自动检测并转换。
pgr 后续如果要直接消费 minimap2 输出的 PAF，可以考虑支持。

### 13.3 wgatools CIGAR 统计维度（比 impg 更丰富）

`Cigar` struct（`parser/cigar.rs:16-29`）：

| 字段 | 含义 |
|------|------|
| `match_count` / `mismatch_count` | 匹配/错配碱基数 |
| `ins_event` / `ins_count` | 插入事件数 / 碱基数 |
| `del_event` / `del_count` | 删除事件数 / 碱基数 |
| `inv_ins_event` / `inv_ins_count` | 倒位区域的插入事件数 / 碱基数 |
| `inv_del_event` / `inv_del_count` | 倒位区域的删除事件数 / 碱基数 |
| `inv_event` | 倒位事件数 |

pgr V1 不需要 `inv_*` 字段（两序列 MAF 不涉及倒位），但 `match/mismatch/ins/del`
的事件数 vs 碱基数的区分在计算 gi/bi 时已经用到。

