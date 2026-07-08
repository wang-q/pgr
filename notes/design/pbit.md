# pbit 设计笔记

## 目标

借鉴 AGC (Assembled Genomes Compressor) v3.2.3 的 C++ 算法（LZ-diff、段级参考压缩、k-mer minimizer
参考选择，源码分析见 [agc-cpp.md](../references/agc-cpp.md)），设计 pgr 原生的群体基因组压缩格式
`pbit`（p = population 或 plus；2bit 记录 + delta），集成到 pgr 作为 `pgr pbit` 命令族。
**不兼容** C++ AGC 的 `.agc` 文件格式。

**关键约束**：

- **不兼容 `.agc` 文件格式**：设计原生 `pbit` 格式（扩展名 `.pbit`），仅借鉴 AGC 的算法（LZ-diff、
  段级参考压缩、k-mer minimizer 参考选择），不移植 CArchive 容器 / varint / 前缀编码
- 不引入 zstd 依赖，使用已有的 `flate2`（gzip/DEFLATE）替代
- **参考层直接复用 2bit 记录格式**（`dna_size + n_blocks + mask_blocks + packed_dna`），保留 N
  blocks 和 mask blocks（AGC 丢失了 mask），提取 `TwoBitFile::read_sequence` 核心逻辑为共享函数
  （见 §复用-2）
- 深度复用 pgr 现有基础设施（`libs/fmt/twobit`、`intspan::Range`、`lru::LruCache`、`libs/fmt/fa`）
- 遵循 pgr 分层原则：复杂逻辑放 `libs/pbit/`，`cmd_pgr/pbit/` 保持薄壳

## pgr 现有基础设施复用

### 1. 随机访问模式 — `libs/loc.rs`

pbit 的随机访问基于**文件偏移**（不再有 CArchive 的 stream/part 抽象）： `Decompressor` 解析
Header + Footer + Index 后，得到各参考段和 delta 的文件偏移，按需 seek + read。

| pbit 组件                                      | pgr Rust                                            | 说明             |
|------------------------------------------------|-----------------------------------------------------|------------------|
| `ref_groups[i].segment_offset`（u64 偏移）     | `loc::read_offset(reader, offset, size)`            | 按偏移读取数据块 |
| `delta_data_offset` + delta `packed_size`      | seek + `read_exact`                                 | delta 随机读取   |
| `IndexMap<String, ...>`（参考组/样本索引）     | `IndexMap<String, (u64, usize)>` (name→offset,size) | 名称→位置索引    |

**复用方式**：`Decompressor::read_sequence` 中 seek 到参考段偏移后调用 `read_2bit_record`， seek 到
delta 偏移后按 `packed_size` 读取——与 `loc::read_offset` 的"按偏移读取"模式一致。

参考：[fa/range.rs](file:///Volumes/ExtHome/Scripts/pgr/src/cmd_pgr/fa/range.rs) 的
`loc::open_indexed` + `loc::fetch_record` 模式。

### 2. 直接复用 2bit 记录格式 — `libs/fmt/twobit`（参考层底层）

pbit 的**参考层直接复用标准 2bit 记录**（`dna_size + n_blocks + mask_blocks + packed_dna`），
不是"借鉴模式"而是"调用同一函数"。这样 2bit 和 pbit 共享参考段的读写代码，pbit 仅在 2bit 记录
之上叠加 delta 层。

#### 2.1 提取共享函数 `read_2bit_record` / `write_2bit_record`

[twobit.rs:384-451](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs#L384) 的
`TwoBitFile::read_sequence` 内嵌了"从当前 reader 位置读取一个 2bit 记录"的逻辑（dna_size →n_blocks
→ mask_blocks → reserved → packed_dna → unpack）。将其核心抽取为模块级共享函数：

```rust
// libs/fmt/twobit.rs (新增 pub 函数)
/// Read a single 2bit record from the current reader position and return the
/// decoded DNA string with masks applied. Reused by TwoBitFile and pbit::Decompressor.
pub fn read_2bit_record<R: Read + Seek>(
    reader: &mut R,
    is_swapped: bool,
    start: Option<usize>,
    end: Option<usize>,
    no_mask: bool,
) -> Result<String>;

/// Write a single 2bit record (dna_size + n_blocks + mask_blocks + reserved +
/// packed_dna) to the writer. Reused by TwoBitWriter and pbit::Compressor.
pub fn write_2bit_record<W: Write>(
    writer: &mut W,
    dna: &str,
    do_mask: bool,
) -> Result<()>;
```

`TwoBitFile::read_sequence` 改为 `seek(offset) → read_2bit_record(...)` 的薄壳，
`TwoBitWriter::write` 改为循环调用 `write_2bit_record`。pbit 的 `Decompressor` 读参考段时也 调用
`read_2bit_record`，`Compressor` 写参考段时调用 `write_2bit_record`——**参考段字节级兼容 2bit 记录**。

#### 2.2 `SequenceReader` trait — 统一的随机访问接口

[libs/io.rs](file:///Volumes/ExtHome/Scripts/pgr/src/libs/io.rs) 定义了 pgr 的序列随机访问抽象：

```rust
pub trait SequenceReader {
    fn read_sequence(
        &mut self,
        name: &str,
        start: Option<usize>,
        end: Option<usize>,
    ) -> anyhow::Result<String>;
}
```

`TwoBitFile<R>` 在 [twobit.rs:500](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs#L500-L510)
实现了此 trait（委托给固有的 5 参数 `read_sequence`）。**pbit 的 `Decompressor<R>` 也实现此 trait**，
但语义限定为**读参考层序列**（参考 contig 名唯一），供未来 chain/net 按参考坐标查询 pbit
归档——与查询 2bit 完全等价。样本层提取（`getctg` 语义，遍历多样本输出多 FASTA）不走
`SequenceReader`，改调 `Decompressor::get_contig`（见 §2.6）。

#### 2.3 `TwoBitFile<R: Read + Seek>` 结构模式 → `Decompressor<R>`

twobit.rs 的 `TwoBitFile<R>` 结构
（[twobit.rs:277-284](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs#L277-L284)）展示了
pgr 处理二进制基因组格式的标准模式，pbit 读取器逐项镜像：

- **泛型边界** `R: Read + Seek`
    - `TwoBitFile<R>` — 2bit 读取器
    - `Decompressor<R>` — pbit 读取器（直接持有 `R`，无中间容器层）
- **三构造器模式**
  （[twobit.rs:286-360](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs#L286-L360)）：
    - `open<P: AsRef<Path>>(path)` → 文件支持的 `BufReader<File>`
    - `open_and_read<P: AsRef<Path>>(path)` → 内存 `Cursor<Vec<u8>>`（小文件/测试用）
    - `new(reader)` → 解析 header + 索引（2bit）/ 解析 header + ref_index + sample_index（pbit）
- **名称→位置索引**：
    - 2bit 单层：`sequence_offsets: HashMap<String, u64>`（name → 文件偏移，
      [twobit.rs:280](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs#L280)）
    - pbit 双层：参考段层 `ref_groups: Vec<RefGroupEntry>`（每个 group = 一个参考段，含
      `contig_name` + `segment_offset`）+ `contig_groups: IndexMap<String, Vec<u32>>`（contig 名 →
      ref_group_id 列表，供 `SequenceReader` 按 contig 名拼接多段）+ 样本层
      `collection: Collection`（sample → 各 contig 的段描述，供 `get_contig`/`get_sample` 遍历样本）
- **坐标提取 API**：
    - `TwoBitFile::read_sequence(name, start, end, no_mask)` — seek 到偏移，`read_2bit_record` 解包，
      应用 mask（[twobit.rs:384](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs#L384)）
    - `Decompressor::read_sequence(name, start, end)`（`SequenceReader`，读**参考层**）— 查
      `contig_groups[name]` 得 `Vec<ref_group_id>` → 遍历各 id: seek
      `ref_groups[id].segment_offset` → `read_2bit_record` → 拼接 → 按 `[start, end)` 切片
    - `Decompressor::get_contig(contig, start, end, strand, out)`（读**样本层**）— 遍历
      `collection.samples` → 对含此 contig 的样本：查 `SegmentDescs` → 每段按
      `ref_groups[ref_group_id].segment_offset` seek → `read_2bit_record` 读参考 → seek delta →
      flate2 解压 → LZ-diff 解码 →（if `is_rev_comp`: 反向互补）→ 拼接 → 切片 →（if
      `strand == "-"`: `rev_comp`）→ 写一条 FASTA

#### 2.4 `TwoBitWriter<W: Write>` 模式 → `Compressor<W>`

twobit.rs 的 `TwoBitWriter<W>`（
[twobit.rs:168-275](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs#L168-L275)）
展示了写入器模式，pbit 写入器对应：

- `TwoBitWriter<W: Write>` / `new(writer)` / `write(sequences, do_mask)` → 写 header + 索引 + 循环
  `write_2bit_record`
- `Compressor<W: Write + Seek>` / `new(writer, ...)` / `finish()` → 写 header → 写参考段
  （`write_2bit_record`）→ 写 delta 数据 → 写 sample_index → 写 footer（footer 在文件末尾，需
  `Seek`）。**参考段字节级与 2bit 记录一致**

#### 2.5 二进制解析工具函数

twobit.rs 的模块级 helper（`read_u32`、`read_u64`、`read_u32_vec`，
[twobit.rs:512-540](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs#L512-L540)）展示了
pgr 读取固定大小小端整数的模式。Phase 0 已将这些函数改为 `pub`（供未来复用），但 pbit 实际实现中
为避免冗余的 `is_swapped` 参数（pbit 统一小端，无需字节序切换），在 `format.rs` 内定义了独立的
`read_u32_le` / `read_u64_le` / `write_u32_le` / `write_u64_le`。记录级的
`read_2bit_record` / `write_2bit_record` 则被 pbit 直接复用（见 §复用-2.4）。

#### 2.6 命令层复用 — `pbit/range.rs` 与 `twobit/range.rs` 共享区间工具

[twobit/range.rs](file:///Volumes/ExtHome/Scripts/pgr/src/cmd_pgr/twobit/range.rs) 和
`pbit/range.rs` 共享**区间解析与输出工具**（`ranges_arg()` / `collect_ranges()` /`intspan::Range`
/ `nt::rev_comp` / `outfile_arg()` / `get_outfile()`），但提取逻辑不同——twobit 每个 range 读 1
条序列（`SequenceReader::read_sequence`），pbit 每个 range 读所有样本 （`get_contig`，`getctg`
语义，输出多 FASTA）：

```rust
// twobit/range.rs (现有)                     // pbit/range.rs (待实现)
let mut tb = TwoBitFile::open(infile)?;        let mut pbit = Decompressor::open(infile)?;
for el in ranges.iter() {                      for el in ranges.iter() {
    let rg = intspan::Range::from_str(el);         let rg = intspan::Range::from_str(el);
    let seq_id = rg.chr();                         let seq_id = rg.chr();
    if !tb.sequence_offsets.contains_key(seq_id)   if !pbit.contains_contig(seq_id)
        { warn; continue; }                            { warn; continue; }
    let (start, end) = ...;  // 1-based→0-based     let (start, end) = ...;  // 同
    let mut seq =                                  // 遍历所有样本，每个含此 contig 的样本写一条 FASTA；
        tb.read_sequence(seq_id, start, end, false)?;  // strand 由 get_contig 内部逐条 rev_comp
    if rg.strand() == "-" {                        pbit.get_contig(seq_id, start, end,
        seq = rev_comp(seq);                                              rg.strand(), &mut writer)?;
    }                                              }
    writeln!(writer, ">{}", rg)?;              }
    writeln!(writer, "{}", seq)?;
}                                              }
```

### 3. 坐标/区间基础设施

| 组件                                    | 位置                                                               | pbit 用途                               |
|-----------------------------------------|--------------------------------------------------------------------|-----------------------------------------|
| `intspan::Range`                        | crate                                                              | `range` 命令的区间解析（`chr1:1-1000`） |
| `args::ranges_arg()`                    | [args.rs](file:///Volumes/ExtHome/Scripts/pgr/src/cmd_pgr/args.rs) | CLI 区间参数                            |
| `args::collect_ranges()`                | args.rs                                                            | 收集区间列表                            |
| `args::infile_arg_required_with_help()` | args.rs                                                            | 输入文件参数                            |
| `args::outfile_arg()` / `get_outfile()` | args.rs                                                            | 输出文件参数                            |

**复用方式**：`pgr pbit range` 可支持 `chr1:1-1000` 形式的区间提取，复用 `ranges_arg` +
`collect_ranges`，与 `fa/range.rs`、`twobit/range.rs` 保持一致。

### 4. FASTA I/O — `libs/fmt/fa`

- `pgr::libs::fmt::fa::reader(infile)` — noodles-based FASTA reader（支持 gzip）
- `pgr::libs::fmt::fa::writer(outfile)` — FASTA writer
- `libs::fasta` — chunk/dedup/filter/stat 等操作

**复用方式**：AGC 的 `genome_io.cpp`（FASTA 读写）完全不需要移植。压缩端用 `fa::reader` 读输入，
解压端用 `fa::writer` 写输出。

### 5. LRU 缓存 — `lru::LruCache`

[fa/range.rs](file:///Volumes/ExtHome/Scripts/pgr/src/cmd_pgr/fa/range.rs#L77) 已用
`lru::LruCache<String, noodles_fasta::Record>` 缓存序列。

**复用方式**：AGC 解压器缓存解压的参考段和 delta pack，用 `LruCache<u32, String>`
（ref_group_id → 解码参考段；group = segment，单 u32 键即可）。

### 6. 压缩后端 — `flate2`（已有依赖）

AGC 用 zstd 的两处，均改为 flate2：

| AGC 用途                       | zstd API                                    | flate2 替代                                            |
|--------------------------------|---------------------------------------------|--------------------------------------------------------|
| 段级 delta 压缩（`CSegment`）  | `ZSTD_compressCCtx` / `ZSTD_decompressDCtx` | `flate2::write::GzEncoder` / `flate2::read::GzDecoder` |
| 元数据压缩（`CCollection_V3`） | 同上                                        | 同上                                                   |

**权衡**：gzip 压缩率略低于 zstd（通常差 5-15%），但无需新依赖。gzip 的压缩速度与 zstd level 3-6
相当，解压速度 gzip 略快。

**实现**：

```rust
use flate2::{write::GzEncoder, read::GzDecoder, Compression};

fn compress(data: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(data).unwrap();
    encoder.finish().unwrap()
}

fn decompress(data: &[u8]) -> Vec<u8> {
    let mut decoder = GzDecoder::new(data);
    let mut decoded = Vec::new();
    decoder.read_to_end(&mut decoded).unwrap();
    decoded
}
```

### 7. 其他可复用

| 组件                                   | 位置                                                          | pbit 用途                         |
|----------------------------------------|---------------------------------------------------------------|-----------------------------------|
| `libs::nt::rev_comp`                   | [libs/nt](file:///Volumes/ExtHome/Scripts/pgr/src/libs/nt.rs) | LZ-diff 的反向互补参考            |
| `indexmap::IndexMap`                   | crate                                                         | 保序 HashMap（pgr 统一模式）      |
| `minimizer-iter`                       | crate                                                         | 参考选择的 minimizer 计算         |
| `murmurhash3`                          | crate                                                         | LZ-diff 哈希表（与 C++ AGC 的 MurMur64 对齐） |
| `rayon`                                | crate                                                         | 并行压缩（多 contig 同时编码）    |

## Rust 架构建议

### 模块结构

```
src/libs/pbit/
├── mod.rs              # 模块导出
├── format.rs           # Header/Footer/Index 序列化、版本常量、二进制 helper
├── lz_diff.rs          # CLZDiffBase/V1/V2 → LzDiff (LZ-diff 算法)
├── segment.rs          # CSegment → Segment (段级 delta 压缩/解压)
├── collection.rs       # CCollection_V3 → Collection (sample/contig/segment 元数据)
├── compressor.rs       # CAGCCompressor → Compressor (写 header + 参考段 + delta + index + footer)
└── decompressor.rs     # CAGCDecompressor → Decompressor (读 header + index + 按需解压段)

src/cmd_pgr/pbit/
├── mod.rs              # 子命令注册
├── create.rs           # pgr pbit create
├── append.rs           # pgr pbit append
├── to_fa.rs            # pgr pbit to-fa
├── some.rs             # pgr pbit some
├── range.rs            # pgr pbit range
└── stat.rs             # pgr pbit stat
```

> **无 `archive.rs`**：pbit 不用 CArchive 多流容器，无需单独模块。Compressor/Decompressor 直接持有
> `W: Write + Seek` / `R: Read + Seek`，按文件偏移读写参考段（2bit 记录）和 delta 数据块。

### 核心数据结构

```rust
// format.rs — Header/Footer + binary helpers (independent read_u32_le / read_u64_le etc.;
//             reuses twobit's read_2bit_record / write_2bit_record at record level)
pub const PBIT_MAGIC: u32 = 0x54494250;  // 'PBIT' — native "2bit + delta" format
pub const PBIT_VERSION_MAJOR: u32 = 1;
pub const PBIT_VERSION_MINOR: u32 = 0;

/// File header (fixed 36 bytes, at file start).
pub struct PbitHeader {
    pub magic: u32,           // PBIT_MAGIC
    pub version: u32,         // major*1000 + minor
    pub segment_size: u32,
    pub kmer_len: u32,
    pub min_match_len: u32,   // LZ-diff min match length (CLI -l, default 18)
    pub ref_group_count: u32,
    pub sample_count: u32,
    pub ref_records_offset: u64,  // offset to Reference Records (usually 36)
}

/// Footer (fixed 24 bytes, at end of file). Located by seeking to file_size - 24.
pub struct PbitFooter {
    pub ref_index_offset: u64,    // offset to reference group index
    pub delta_data_offset: u64,   // offset to delta data section
    pub sample_index_offset: u64, // offset to sample index
}

/// Reference group index entry (in ref_index section). Each group = one reference
/// segment (one 2bit record), mirroring C++ AGC where group = segment.
pub struct RefGroupEntry {
    pub contig_name: String,  // contig this segment belongs to
    pub segment_offset: u64,  // file offset to one 2bit record
}

/// Delta entry (in delta data section).
pub struct DeltaEntry {
    pub is_rev_comp: bool,
    pub raw_length: u32,
    pub packed_data: Vec<u8>,       // flate2(LZ-diff encoded)
    // packed_size is not stored as a field; it's derived via `meta()` → packed_data.len()
}

// lz_diff.rs — LZ-diff encoder/decoder (V2), mirrors C++ CLZDiff_V2.
pub struct LzDiff {
    reference: Vec<u8>,       // 2-bit encoded (A=0,C=1,G=2,T=3,N=4, other=31)
    ht: HashTable,            // u16 or u32 positions (built lazily by prepare_index)
    min_match_len: u32,       // default 18
    key_len: u32,             // = min_match_len - hashing_step + 1 (= 15)
    key_mask: u64,            // 2 * key_len bits
    hashing_step: u32,        // 4 (USE_SPARSE_HT)
    ht_size: u64,             // hash table size
    ht_mask: u64,             // hash table mask (= ht_size - 1)
    index_ready: bool,        // lazy index build flag
    min_nrun_len: u32,        // 4
}

enum HashTable {
    Short(Vec<u16>),   // ref < 65535 * hashing_step
    Long(Vec<u32>),
}

impl LzDiff {
    pub fn new(min_match_len: u32) -> Self;
    /// 2-bit encode + store reference. Required for both encode and decode.
    pub fn prepare(&mut self, reference: &[u8]);
    /// Build hash table over the stored reference. Required only for encode
    /// (decode does not need the hash table). Split from `prepare` so that
    /// `Decompressor` can skip the expensive index build.
    pub fn prepare_index(&mut self);
    pub fn encode(&mut self, text: &[u8], encoded: &mut Vec<u8>);
    /// Decode using `self.reference` (stored by `prepare`). No external
    /// reference parameter — avoids redundancy with the internal field.
    pub fn decode(&self, encoded: &[u8], decoded: &mut Vec<u8>);
    pub fn estimate(&mut self, text: &[u8], bound: u32) -> usize;
}

// segment.rs — delta encoding/decoding per reference group (no archive dependency)
/// Segment-level delta compression/decompression.
/// Holds a single `LzDiff` (which stores the 2-bit encoded reference internally);
/// does NOT duplicate the reference as ASCII.
pub struct Segment {
    lz_diff: LzDiff,
}

impl Segment {
    /// Prepare reference (call before add/get). Delegates to `lz_diff.prepare`.
    pub fn prepare(&mut self, ref_dna: &[u8]);
    /// Build LZ-diff hash table (call before add, not needed for get).
    pub fn prepare_index(&mut self);
    /// Encode `seq` against the prepared reference, return uncompressed delta bytes.
    pub fn add(&mut self, seq: &[u8]) -> Vec<u8>;
    /// Decode `delta` (uncompressed) against the prepared reference, return raw sequence.
    pub fn get(&self, delta: &[u8]) -> anyhow::Result<Vec<u8>>;
}

// collection.rs — sample/contig/segment metadata (fixed u32 fields, flate2 compressed)
pub struct SegmentDesc {
    pub ref_group_id: u32,    // index into ref_groups
    pub delta_id: u32,        // index into ref_group's delta list
    // is_rev_comp / raw_length live in DeltaEntry (properties of the encoding,
    // shared by all segments pointing to the same delta).
}

pub struct ContigSegs {
    pub contig_name: String,
    pub segments: Vec<SegmentDesc>,
}

pub struct Collection {
    pub samples: IndexMap<String, Vec<ContigSegs>>,  // sample → contigs
    pub cmd_line: String,
}

impl Collection {
    pub fn register_sample_contig(&mut self, sample: &str, contig: &str);
    pub fn add_segment(&mut self, sample: &str, contig: &str,
                       ref_group_id: u32, delta_id: u32);
    pub fn get_contig_segments(&self, sample: &str, contig: &str) -> Option<&[SegmentDesc]>;
    pub fn list_samples(&self) -> Vec<&str>;
    pub fn list_contigs(&self, sample: &str) -> Vec<&str>;
    /// Serialize with fixed-size u32 LE fields + flate2 (no prefix coding).
    pub fn serialize(&self) -> Vec<u8>;
    pub fn deserialize(data: &[u8]) -> anyhow::Result<Self>;
}

// compressor.rs — holds W: Write + Seek directly (no archive wrapper)
pub struct Compressor<W: Write + Seek> {
    writer: W,
    header: PbitHeader,
    ref_groups: Vec<RefGroupEntry>,         // one per reference segment (group = segment)
    deltas: Vec<Vec<DeltaEntry>>,           // deltas[ref_group_id][delta_id]
    collection: Collection,
    segments: Vec<Segment>,                 // one per ref group (each Segment holds its own min_match_len)
    segment_size: usize,
    kmer_len: usize,
    // min_match_len is not stored as a field; it's passed into Segment::new() and
    // held there. The header still records it (PbitHeader.min_match_len).
}

impl Compressor<BufWriter<File>> {
    /// Create from output path + reference FASTA.
    /// Writes header (placeholder) + reference records (via write_2bit_record).
    /// The CLI `create` command calls this, then `append_sample` for each `-i`
    /// input FASTA, then `finish`.
    pub fn create<P: AsRef<Path>>(out_path: P, ref_fasta: &str, segment_size: usize,
                                  kmer_len: usize, min_match_len: u32) -> anyhow::Result<Self>;
    /// Open an existing `.pbit` for appending samples (powers `pgr pbit append`).
    /// Uses Decompressor to read existing header / ref_groups / deltas (with
    /// packed_data) / collection, reads each reference segment via
    /// `read_2bit_record` and rebuilds `Segment`s (prepare + prepare_index),
    /// positions the writer at `footer.ref_index_offset` (start of Reference
    /// Index). `finish` then rewrites Reference Index + Delta Data + Sample Index
    /// + Footer as a contiguous block. No reference FASTA needed — reuses
    /// reference records already embedded in the archive.
    pub fn open_for_append<P: AsRef<Path>>(in_path: P) -> anyhow::Result<Self>;
}

impl<W: Write + Seek> Compressor<W> {
    /// Append a sample: read FASTA → segment → k-mer minimizer reference selection
    /// → LZ-diff encode → flate2 compress → delta dedup → store. Sample name is
    /// derived from the FASTA basename by the caller (CLI layer).
    pub fn append_sample(&mut self, sample_name: &str, fasta_path: &str) -> anyhow::Result<()>;
    /// Finalize: write Reference Index → Delta Data → Sample Index → Footer →
    /// patch Header offsets. Called once after all samples are appended.
    pub fn finish(self) -> anyhow::Result<()>;
}

// decompressor.rs — holds R: Read + Seek directly (no archive wrapper)
pub struct Decompressor<R: Read + Seek> {
    reader: R,
    header: PbitHeader,
    footer: PbitFooter,
    ref_groups: Vec<RefGroupEntry>,
    contig_groups: IndexMap<String, Vec<u32>>,  // contig name → ref_group_id list (for SequenceReader)
    contig_set: HashSet<String>,                // all contig names in collection (for contains_contig)
    collection: Collection,
    ref_cache: LruCache<u32, String>,           // ref_group_id → decoded ref segment
    // Delta metadata loaded during `new` by scanning Delta Data headers (9 bytes each).
    // Stores is_rev_comp / raw_length / packed_size per delta so that `get_contig`
    // can compute segment coordinates without re-reading headers from disk.
    delta_meta: Vec<Vec<DeltaMeta>>,            // [ref_group_id][delta_id] → (is_rc, raw_len, packed_size)
    delta_offsets: Vec<Vec<u64>>,              // delta_offsets[ref_group_id][delta_id] → file offset of header
    delta_cache: LruCache<(u32, u32), Vec<u8>>,  // (ref_group_id, delta_id) → decoded raw seq
}

/// In-memory delta header (loaded during `new`, no packed_data).
pub struct DeltaMeta {
    pub is_rev_comp: bool,
    pub raw_length: u32,
    pub packed_size: u32,
}

impl Decompressor<BufReader<File>> {
    /// Open from file path (mirrors TwoBitFile::open).
    pub fn open<P: AsRef<Path>>(path: P) -> anyhow::Result<Self>;
}

impl Decompressor<std::io::Cursor<Vec<u8>>> {
    /// Open and read entire file into memory (mirrors TwoBitFile::open_and_read).
    pub fn open_and_read<P: AsRef<Path>>(path: P) -> anyhow::Result<Self>;
}

impl<R: Read + Seek> Decompressor<R> {
    /// Construct from an already-opened reader: parse header + footer + indexes
    /// (mirrors TwoBitFile::new parsing header + sequence_offsets). Builds
    /// `contig_groups`, `contig_set`, and scans Delta Data headers (9 bytes each)
    /// to populate `delta_meta` + `delta_offsets` without decompressing data.
    pub fn new(reader: R) -> anyhow::Result<Self>;

    /// Check if a contig name exists in any sample's collection (getctg/getset
    /// semantics). Mirrors tb.sequence_offsets.contains_key but queries the
    /// sample layer, not the reference layer. Built once in `new` as a
    /// `HashSet<String>` from `collection.samples[*][*].contig_name`.
    pub fn contains_contig(&self, name: &str) -> bool;
    pub fn list_samples(&self) -> Vec<&str>;
    pub fn list_contigs(&self, sample: Option<&str>) -> Vec<&str>;
    pub fn get_sample(&mut self, sample: &str, out: &mut impl Write) -> anyhow::Result<()>;
    /// Extract a contig from ALL samples (getctg semantics), optionally sliced to
    /// [start, end). Writes one FASTA entry per sample that has this contig.
    /// Uses `delta_meta` to compute segment coordinates and only decodes segments
    /// overlapping [start, end) (smart selection, like the reference layer).
    /// Output header: `>{sample_name} {contig}:{start}-{end}({strand})`.
    /// If strand is "-", each sequence is reverse-complemented before writing.
    pub fn get_contig(&mut self, contig: &str, start: Option<usize>, end: Option<usize>,
                      strand: &str, out: &mut impl Write) -> anyhow::Result<()>;
}

// SequenceReader reads the REFERENCE layer (reference contig names are unique),
// NOT sample sequences. This lets future chain/net modules query pbit by reference
// coordinate just like they query 2bit. Sample extraction goes through get_contig /
// get_sample (multi-sample output, doesn't fit SequenceReader's single-return contract).
//
// **Assumption**: reference contig names and sample contig names share the same
// namespace (e.g. both use "chr1"). This holds when reference and samples are from
// the same species assembly. `range` (getctg) queries the sample layer via
// `contains_contig`; `SequenceReader` queries the reference layer via
// `contig_groups` — both keyed by contig name, so the assumption must hold for
// a contig name to be resolvable in both layers.
impl<R: Read + Seek> crate::libs::io::SequenceReader for Decompressor<R> {
    fn read_sequence(
        &mut self,
        name: &str,
        start: Option<usize>,
        end: Option<usize>,
    ) -> anyhow::Result<String> {
        // name = reference contig name → contig_groups[name] = Vec<ref_group_id> →
        //   walk segments, accumulate lengths, skip those before `start`, read only
        //   segments overlapping [start, end) via read_2bit_record, concat the slice.
        // (Mirrors TwoBitFile::read_sequence: seek → read_2bit_record → slice,
        //  but spans multiple segments per contig via contig_groups without
        //  materializing the full contig.)
        ...
    }
}
```

### CLI 设计

遵循 twobit 的单字子命令 + 分组惯例（`2bit` 有 info/subset/transform 三组，见
[twobit/mod.rs](file:///Volumes/ExtHome/Scripts/pgr/src/cmd_pgr/twobit/mod.rs)），pbit 命令设计如下
（括注对应的 C++ AGC 命令，源码分析见 [agc-cpp.md](../references/agc-cpp.md)）：

- `create` → `create`（build 组，无 twobit 对应）
- `append` → `append`（build 组，无 twobit 对应）
- `getcol` → `to-fa`（transform 组，对应 `2bit to-fa`）
- `getset` → `some`（subset 组，对应 `2bit some`）
- `getctg` → `range`（subset 组，对应 `2bit range`）
- `info` + `listref` + `listset` + `listctg` → `stat`（info 组，对应 `2bit size`，用 flag 区分输出）

```
pgr pbit create    -i input.fa [-i input2.fa ...] -o out.pbit -r ref.fa [-s segment_size] [-k kmer_len] [-l min_match_len]
pgr pbit append    in.pbit -i input.fa [-i input2.fa ...] [-o out.pbit]  # 归档为位置参数输入，-o 省略时原地修改
pgr pbit to-fa     in.pbit -o out_dir/          # 提取所有样本为 FASTA，每样本一个文件 out_dir/{sample}.fa
pgr pbit some      in.pbit sample_list.txt [-o out.fa] [--invert]  # 按样本名列表提取，输出多 FASTA
pgr pbit range     in.pbit "chr1" "chr2:1-1000" [-o out.fa]  # 按 contig/区间提取（遍历所有样本，输出多 FASTA）
pgr pbit stat      in.pbit [--samples | --refs | --contigs [-s sample]]  # 统计/列表
```

> **参数说明**：
> - `-s segment_size`：分段大小（bp，默认 4096）
> - `-k kmer_len`：k-mer minimizer 参考选择的 k-mer 长度（默认 15）
> - `-l min_match_len`：LZ-diff 最小匹配长度（默认 18，对应 key_len=15）
> - 样本名从 `-i` 指定的 FASTA 文件 basename 派生（`libs::io::get_basename`），如
>   `path/to/sample1.fa` → 样本名 `sample1`
> - `append` 的归档文件为位置参数（输入），`-o` 可选（省略时原地修改 `in.pbit`，
>   指定时先复制再追加）；参考已嵌入归档，无需 `-r`
> - `some` / `range` 输出 FASTA header 格式：`>{sample_name} {contig}:{start}-{end}({strand})`
>   （`some` 无区间时为 `>{sample_name} {contig}`）

> **stat flag 语义**：
> - 默认（无 flag）：输出归档总览统计（ref_group_count、sample_count、segment_size、kmer_len、
>   min_match_len、各参考 contig 段数、各样本 contig 数）。
> - `--samples`：列出所有样本名（对应 C++ `listset`）。
> - `--refs`：列出参考 contig 名（来自 `ref_groups[*].contig_name`，去重保序；对应 C++ `listref`，
>   但 C++ 列"参考样本名"——pbit 无参考样本实体，参考由 contig 名标识）。
> - `--contigs [-s sample]`：无 `-s` 列出所有样本的 `sample\tcontig` 对（对应 C++ `listctg`）；
>   有 `-s` 仅列指定样本的 contig 名。

子命令分组（镜像 `2bit` mod.rs 结构）：

- build: `create` / `append`
- info: `stat`
- subset: `range` / `some`
- transform: `to-fa`

`range` 复用 `ranges_arg` 支持区间提取：`pgr pbit range in.pbit "chr1:1-1000" "chr2(-):500-1000"`。
区间解析与 `pgr 2bit range` 一致（共享 `ranges_arg` / `intspan::Range`），但提取语义不同——pbit
遍历所有样本输出多 FASTA（`getctg` 语义），2bit 读单序列（见 §复用-2.6）。

## 文件格式规范

pbit 格式（原生"2bit + delta"，扩展名 `.pbit`，区别于 C++ AGC 的 `.agc`）。所有整数使用
**固定大小 小端序**（u32 = 4 字节，u64 = 8 字节），不使用 varint / 前缀编码。字符串采用**长度前缀**
（u32 len + UTF-8 bytes），不使用 null 终止。

### 文件结构总览

```
┌─────────────────────────────────────┐  ← offset 0
│ Header (固定 36 字节)               │
├─────────────────────────────────────┤
│ Reference Records                   │  ← 参考段，每段为标准 2bit 记录
│   ref_group 0: seg 0, seg 1, ...    │     (dna_size + n_blocks + mask_blocks
│   ref_group 1: seg 0, seg 1, ...    │      + reserved + packed_dna)
│   ...                               │
├─────────────────────────────────────┤  ← footer.ref_index_offset
│ Reference Index                     │  ← 每个参考组的名称 + 段文件偏移列表
├─────────────────────────────────────┤  ← footer.delta_data_offset
│ Delta Data                          │  ← 每个参考组的 delta 列表
│   ref_group 0: delta 0, delta 1,... │     (is_rc + raw_len + packed_size +
│   ref_group 1: delta 0, delta 1,... │      flate2(LZ-diff encoded))
├─────────────────────────────────────┤  ← footer.sample_index_offset
│ Sample Index (collection)           │  ← flate2(序列化的 samples/contigs/segments)
├─────────────────────────────────────┤
│ Footer (固定 24 字节)               │  ← 三个 section 的偏移
└─────────────────────────────────────┘  ← EOF
```

### Header（36 字节，文件起始）

```
offset  size  field              说明
0       4     magic              0x54494250 ('PBIT', 小端)
4       4     version            major*1000 + minor (当前 1000)
8       4     segment_size       分段大小（bp，如 4096）
12      4     kmer_len           k-mer 长度（如 15）
16      4     min_match_len      LZ-diff 最小匹配长度（如 18，对应 -l 参数）
20      4     ref_group_count    参考段总数（每段 = 一个 group，非 contig 数）
24      4     sample_count       样本数（不含参考）
28      8     ref_records_offset Reference Records 起始偏移（通常 = 36）
```

### Reference Records（标准 2bit 记录，连续存储）

每个参考段为一个**标准 2bit 记录**（与
[twobit.rs:384](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs#L384)`read_sequence`
读取的格式字节级一致），由共享函数 `write_2bit_record` / `read_2bit_record`读写：

```
offset  size           field
0       4              dna_size                 序列长度
4       4              n_block_count            N-block 数
8       4*n_block_count n_starts                N-block 起始位置
        4*n_block_count n_sizes                 N-block 长度
        4              mask_block_count         Mask-block 数
        4*mask_block_count mask_starts          Mask-block 起始位置
        4*mask_block_count mask_sizes           Mask-block 长度
        4              reserved                 保留字段（0）
        (dna_size+3)/4 packed_dna               2-bit 打包 DNA
```

> **与 C++ AGC 的关键差异**：保留 mask blocks（C++ AGC 丢失 mask）；参考段不经过 zstd/flate2
> 二次压缩（2bit 打包已压缩 4 倍，delta 层负责进一步压缩）。

### Reference Index（参考组索引）

每个 ref_group = 一个参考段（一个 2bit 记录）：

```
offset  size           field
0       4              ref_group_count          (= header.ref_group_count)
for each ref_group:
  4      name_len       contig 名字节长度
  name_len name_bytes   contig 名（UTF-8，该段所属的参考 contig）
  8      segment_offset 该段的文件偏移（u64 LE，指向 Reference Records 中的一个 2bit 记录）
```

### Delta Data（delta 数据区）

每个参考组有一组唯一 delta（去重后的 LZ-diff 编码）。同一参考组的 delta 连续存储：

```
offset  size           field
0       4              ref_group_count
for each ref_group:
  4      delta_count    唯一 delta 数
  for each delta:
    1      is_rev_comp   0 = 正向, 1 = 反向互补
    4      raw_length    原始序列长度
    4      packed_size   flate2 压缩后字节数
    packed_size packed_data  flate2(LZ-diff encoded bytes)
```

> **随机访问**：delta 为变长记录（`packed_size` 各异），不在文件中存储 per-delta offset 表。
> `Decompressor::new` 顺序扫描一次 Delta Data 区，按每条 delta 头部的 `packed_size` 累加偏移，
> 构建内存中的 `delta_offsets: Vec<Vec<u64>>`（`[ref_group_id][delta_id]` → 文件偏移）。
> 扫描仅需读取每条 delta 的 9 字节头部（`is_rev_comp + raw_length + packed_size`），不解压数据。

### Sample Index（collection，flate2 压缩）

序列化为固定大小字段后用 flate2 压缩：

```
[未压缩的逻辑结构]
4      sample_count
for each sample:
  4      name_len
  name_len name_bytes     样本名
  4      contig_count
  for each contig:
    4      contig_name_len
    contig_name_len contig_name_bytes
    4      segment_count
    for each segment:
      4      ref_group_id    参考组索引
      4      delta_id        该参考组内的 delta 索引
      (is_rev_comp / raw_length 由 DeltaEntry 提供，此处不重复)
4      cmd_line_len
cmd_line_len cmd_line_bytes  命令行（记录创建参数）
```

### Footer（24 字节，文件末尾）

```
offset  size  field
0       8     ref_index_offset     Reference Index 起始偏移
8       8     delta_data_offset    Delta Data 起始偏移
16      8     sample_index_offset  Sample Index 起始偏移
```

> **设计要点**：

> - **无 magic number 在末尾**：Header 已含 magic，Footer 仅 24 字节偏移
> - **无变长整数**：所有字段固定大小，简化解析
> - **参考段 = 2bit 记录**：可直接用 `read_2bit_record` 读取，与 twobit.rs 共享代码
> - **delta 层独立压缩**：每个 delta 单独 flate2 压缩，支持随机访问单个样本/contig
> - **collection 整块压缩**：元数据通常较小，一次 flate2 压缩即可
> - **冗余计数字段**：`ref_group_count`/`sample_count` 在 Header、Reference Index、Delta Data、
>   Sample Index 中重复出现，这是有意为之——读取各 section 时可就地校验一致性（不匹配则报错），
>   避免损坏文件导致越界访问

## 参考资料

- AGC C++ 源码分析: [notes/references/agc-cpp.md](../references/agc-cpp.md)
- AGC GitHub: https://github.com/refresh-bio/agc
- AGC 论文: Deorowicz et al., "AGC: Assembly Genomes Compressor", Bioinformatics (2024)
- pgr loc 模块: [libs/loc.rs](file:///Volumes/ExtHome/Scripts/pgr/src/libs/loc.rs)
- pgr 2bit 模块: [libs/fmt/twobit](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs)
- pgr FASTA I/O: [libs/fmt/fa](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/fa.rs)
