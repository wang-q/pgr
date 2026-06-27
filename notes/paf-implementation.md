# PAF 模块实现参考

本文档基于 impg-0.4.1 源码，梳理 pgr 需要实现的 PAF 相关组件的代码级设计。
这是纯实现参考，不涉及路线决策（见 [[paf-route.md]]）或第一步行动计划（见 [[pairwise-selection.md]]）。

参考源码：`paf.rs`（417 行）、`alignment_record.rs`（138 行）、`seqidx.rs`（56 行）、
`main.rs` 的 `output_results_paf` 函数（`main.rs:11989-12101`）。

---

## 1. 模块结构

```
src/libs/paf/
├── mod.rs          # 模块导出
├── record.rs       # PafRecord struct — PAF 行在内存中的表示
├── parser.rs       # PAF 解析 — 纯文本 / BGZF / GZI 三种模式
├── cigar.rs        # CigarOp bit-packing + 字符串互转
├── writer.rs       # PAF 行格式化输出
└── lazy.rs         # CIGAR 懒加载（从源文件按偏移量读取）

src/libs/seqidx.rs  # SequenceIndex — 序列名↔ID 双向映射（paf 的前置依赖）
```

---

## 2. PafRecord — 核心数据结构

参考 impg 的 `AlignmentRecord`（`alignment_record.rs:12-21`）。

```rust
/// A single PAF alignment record in compact in-memory representation.
///
/// Stores all 12 mandatory PAF columns. CIGAR is NOT stored inline —
/// only file offset + byte length for lazy loading via `read_cigar_data()`.
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
// Total: 48 bytes — fits in a single 64-byte cache line
```

### 设计要点（直接借鉴 impg）

- **strand 编码在 MSB**（impg `alignment_record.rs:34`）：`const STRAND_BIT: u64 = 0x8000000000000000`。
  MSB=0 表示 Forward，MSB=1 表示 Reverse。`strand()`/`set_strand()` 方法封装位操作。
  参考 impg 实现：
  ```rust
  pub fn strand(&self) -> Strand {
      if (self.strand_and_offset & Self::STRAND_BIT) != 0 { Strand::Reverse }
      else { Strand::Forward }
  }
  pub fn set_strand(&mut self, strand: Strand) {
      match strand {
          Strand::Forward => self.strand_and_offset &= !Self::STRAND_BIT,
          Strand::Reverse => self.strand_and_offset |= Self::STRAND_BIT,
      }
  }
  ```

- **CIGAR 懒加载**（[[impg.md]] §3.2）：`strand_and_offset` 的低 63 位是源 PAF 文件中
  `cg:Z:` tag 的字节偏移量。查询时才调用 `read_cigar_data()` 从文件中读取——这是把
  全基因组 PAF 装入区间树的关键优化：区间树节点只存坐标和指针，不存 CIGAR 字符串。

- **序列 ID 化**：`query_id`/`target_id` 是 `SequenceIndex` 中的整数，区间树以 `u32` 为
  key（不用字符串），大幅减少内存和比较开销。

- **字段用 `u32` 而非 `usize`**：跨平台一致 + 更紧凑。对超过 4Gbp 的染色体可升级 `u64`。

### 与 impg `AlignmentRecord` 的差异及理由

| 维度 | impg `AlignmentRecord` | pgr `PafRecord` |
|------|------------------------|-----------------|
| 字段数 | 8 | 12 |
| matches/block_len/mapq | 不存（从原文件重读） | 解析时存储 |
| query_start/end 类型 | `usize` | `u32` |
| cigar_bytes 类型 | `usize` | `u16`（PAF 单行 CIGAR 不超过 64KB） |
| 用途 | 多格式通用（PAF/1ALN/TPA） | PAF 专用 |

**为什么不复用 impg 的 8 字段设计**：impg 把 `matches`/`block_len`/`mapq` 留在原文件、
需要时重读，是为了"单一记录最小化"——区间树节点越小，cache 命中率越高。但 pgr 的第一步
是"能用"，不是"优化到极致"。首次实现存 12 完整字段（48 bytes），后续若内存成瓶颈再
瘦身到 8 字段 + 懒加载。改动接口是 `PafRecord` 内部字段，不影响下游查询逻辑。

---

## 3. PAF 解析器 — 三种模式

impg 的 `paf.rs` 支持三种读取模式。

### 模式 1：纯文本 PAF（`parse_paf`，impg `paf.rs:179-194`）

```rust
/// Parse PAF lines from a buffered reader.
/// Each line = 12+ tab-separated fields. Builds SequenceIndex on the fly.
pub fn parse_paf<R: BufRead>(
    reader: R,
    seq_index: &mut SequenceIndex,
) -> Result<Vec<PafRecord>, PafParseError> {
    let mut bytes_read: u64 = 0;
    let mut records = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let record = parse_paf_line(&line, bytes_read, seq_index)?;
        records.push(record);
        bytes_read += (line.len() + 1) as u64; // include newline
    }
    Ok(records)
}
```

关键细节：
- `bytes_read` 跟踪每个 CIGAR tag 的字节偏移量（用于构建 `strand_and_offset`）
- `SequenceIndex::get_or_insert_id` 在解析过程中动态构建 name→id 映射
- `parse_paf_line`（impg `paf.rs:118-177`）解析单行：12 列 tab split → 数字解析 → strand 判断 →
  name→id 查找 → CIGAR tag 定位

**CIGAR tag 偏移量计算**（impg `paf.rs:150-161`）：
1. 从头扫描 tab-separated 字段
2. 遇到 `cg:Z:` 前缀的字段时，记录其在整个文件中的字节偏移量
3. 偏移量 = `bytes_before` + 前面所有字段的字节数（每个 tab 计 1 字节）
4. `cg:Z:` 前缀 5 字节 + 值 N 字节 → `cigar_offset += 5`，`cigar_bytes = N`

### 模式 2：BGZF 压缩 PAF（`parse_paf_bgzf`，impg `paf.rs:199-270`）

与纯文本的关键差异：
- CIGAR 偏移量用 **BGZF virtual position**（`reader.virtual_position()`）而非
  字节偏移量——因为 BGZF 块边界不可预测
- 每行先记 `line_start_vpos`，读到 `cg:Z:` 后再记 `cigar_vpos`
- 标准双 seek 模式：seek 回行首 → forward scan 到 `cg:Z:` → 记录 vpos → skip 剩余字节到行尾
- 确保 CIGAR 定位正确跨越 BGZF 块边界

```rust
pub fn parse_paf_bgzf<R: Read + Seek>(
    mut reader: bgzf::io::Reader<R>,
    seq_index: &mut SequenceIndex,
) -> Result<Vec<PafRecord>, PafParseError> {
    let mut records = Vec::new();
    let mut line_bytes = Vec::new();
    loop {
        let line_start_vpos = reader.virtual_position();
        line_bytes.clear();
        let bytes_read = reader.read_until(b'\n', &mut line_bytes)?;
        if bytes_read == 0 { break; }
        let line = std::str::from_utf8(&line_bytes[..line_bytes.len()-1])?;

        let mut record = parse_paf_line(line, 0, seq_index)?;
        let cigar_byte_offset = record.strand_and_offset & !AlignmentRecord::STRAND_BIT;

        // Seek back and forward-scan to CIGAR for accurate vpos
        reader.seek(line_start_vpos)?;
        if cigar_byte_offset > 0 {
            std::io::copy(&mut reader.by_ref().take(cigar_byte_offset), &mut std::io::sink())?;
        }
        let cigar_vpos = reader.virtual_position();
        // Update record with BGZF virtual position
        let strand_bit = record.strand_and_offset & AlignmentRecord::STRAND_BIT;
        record.strand_and_offset = u64::from(cigar_vpos) | strand_bit;
        records.push(record);
    }
    Ok(records)
}
```

### 模式 3：BGZF + GZI 多线程（`parse_paf_bgzf_with_gzi`，impg `paf.rs:274-302`）

两遍扫描：
1. 用 `parse_paf` 解析（得到 uncompressed byte offsets）
2. 用 `gzi_index.query()` 把每个 offset 转成 virtual position
3. 好处：多线程 BGZF 解压（`MultithreadedReader`），比单线程快 3-5x

### 统一入口（`parse_paf_file`，impg `paf.rs:306-362`）

```rust
pub fn parse_paf_file(
    paf_file: &str,
    file: File,
    threads: NonZeroUsize,
    seq_index: &mut SequenceIndex,
) -> io::Result<Vec<PafRecord>>
```

自动检测逻辑：
1. 扩展名 `.gz`/`.bgz` → 读 18 字节头判断是否 BGZF（`is_bgzf`）
2. 是 BGZF + 有 `.gzi` 索引 → 模式 3（多线程加速）
3. 是 BGZF 无 `.gzi` → 模式 2（单线程 BGZF reader）
4. 纯文本 → 模式 1

**BGZF 检测的 18 字节头部校验**（impg `paf.rs:50-66`）：
gzip magic（`0x1f 0x8b`）+ DEFLATE（`0x08`）+ FEXTRA（`0x04`）+
XLEN=6（`0x06 0x00`）+ BC subfield（`b'B' b'C'`）+ SLEN=2（`0x02 0x00`）。

普通 gzip 不支持 seek——若用户误用 `gzip` 而非 `bgzip`，必须报错并给出修复命令：
`"Convert with: zcat file.paf.gz | bgzip > output.paf.gz"`（impg `paf.rs:81`）。

---

## 4. CIGAR 编解码

### 4.1 CigarOp bit-packing（impg `impg.rs:74`，[[impg.md]] §3.2 详述）

```rust
/// Compact CIGAR operation: 3-bit op code + 29-bit length in one u32.
#[derive(Debug, Clone, Copy)]
pub struct CigarOp(u32);

impl CigarOp {
    pub fn new(op: char, len: u32) -> Self {
        let code = match op {
            '=' => 0, 'X' => 1, 'I' => 2, 'D' => 3, 'M' => 4, _ => 4,
        };
        CigarOp((len & 0x1FFF_FFFF) | (code << 29))
    }
    pub fn op(self) -> char {
        match self.0 >> 29 {
            0 => '=', 1 => 'X', 2 => 'I', 3 => 'D', _ => 'M',
        }
    }
    pub fn len(self) -> u32 { self.0 & 0x1FFF_FFFF }
}
```

支持 5 种 op：`=`（match）、`X`（mismatch）、`I`（insertion）、`D`（deletion）、
`M`（match/mismatch 不区分）。单个 CIGAR 段最长 512Mbp（29 bits），对基因组比对足够。

### 4.2 CIGAR 字符串 ↔ `Vec<CigarOp>` 转换

```rust
/// Parse CIGAR string "10=5X2I3D" into Vec<CigarOp>
pub fn parse_cigar(s: &str) -> Result<Vec<CigarOp>, CigarParseError>;

/// Format Vec<CigarOp> back to string "10=5X2I3D"
pub fn format_cigar(ops: &[CigarOp]) -> String;
```

### 4.3 CIGAR 懒加载（`read_cigar_data`，impg `paf.rs:68-114`）

```rust
/// Read CIGAR bytes from PAF file at given offset.
/// Dispatches on plain text (File::seek) vs BGZF (bgzf::Reader::seek).
pub fn read_cigar_data(
    paf_path: &str,
    offset: u64,
    byte_len: usize,
) -> Result<Vec<CigarOp>, CigarReadError>;
```

实现细节（impg `paf.rs:72-113`）：
- 根据文件扩展名判断用 `File::seek`（纯文本）还是 `bgzf::Reader::seek`（压缩）
- 纯文本：`file.seek(SeekFrom::Start(offset))` → `file.read_exact(&mut buf)`
- BGZF：`reader.seek(bgzf::VirtualPosition::from(offset))` → `reader.read_exact(&mut buf)`
- 读到的字节是 `cg:Z:` 的原始值，再调用 `parse_cigar` 解析

这是"索引层无 CIGAR、查询层按需读取"的物理基础。impg 内部维护一个
`thread_local!` CIGAR cache（[[impg.md]] §9.5），避免重复 seek——pgr 可后续借鉴。

---

## 5. SequenceIndex — 序列名与 ID 的双向映射

参考 impg `seqidx.rs:1-56`：

```rust
/// Bidirectional mapping between sequence names and compact integer IDs.
/// Used throughout the PAF index to replace string keys with u32.
#[derive(Clone, Debug, Default)]
pub struct SequenceIndex {
    name_to_id: FxHashMap<String, u32>,
    id_to_name: FxHashMap<u32, String>,
    id_to_len: FxHashMap<u32, u32>,
    next_id: u32,
}

impl SequenceIndex {
    pub fn new() -> Self { /* ... */ }

    /// Get existing ID or assign a new one. Optionally record sequence length.
    pub fn get_or_insert_id(&mut self, name: &str, length: Option<u32>) -> u32;

    pub fn get_id(&self, name: &str) -> Option<u32>;
    pub fn get_name(&self, id: u32) -> Option<&str>;
    pub fn get_len(&self, id: u32) -> Option<u32>;
    pub fn len(&self) -> usize;
    pub fn is_empty(&self) -> bool;
}
```

**设计要点**：
- 用 `FxHashMap`（`rustc_hash`）而非标准 `HashMap`——FxHash 比 SipHash 快 2-3 倍，
  且 PAF 解析不受 HashDoS 威胁（key 是序列名，来自自己的数据文件）
- `get_or_insert_id` 在解析时动态构建映射，避免事先扫描全文件
- `id_to_len` 记录序列长度（从 PAF 第 2/7 列提取），用于区间范围的边界校验
- pgr 已有 `fxhash` 依赖（`Cargo.toml:34`），可直接使用，不需要新依赖

---

## 6. PAF 输出格式化

参考 impg `main.rs:11989-12101`（`output_results_paf` 函数）。

```rust
/// Write a single PAF record. 12 mandatory columns + standard tags.
pub fn write_paf_record(
    out: &mut dyn Write,
    query_name: &str,  query_len: u32,  query_start: u32,  query_end: u32,
    target_name: &str, target_len: u32, target_start: u32, target_end: u32,
    strand: char,
    matches: u32, block_len: u32, mapq: u8,
    cigar: &[CigarOp],
    extra_tags: &[(&str, &str)],  // e.g., [("ms", "12345")]
) -> io::Result<()> {
    let cigar_str = format_cigar(cigar);
    write!(out, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        query_name, query_len, query_start, query_end,
        strand,
        target_name, target_len, target_start, target_end,
        matches, block_len, mapq,
    )?;
    let gi = gap_compressed_identity(cigar);
    let bi = block_identity(cigar);
    write!(out, "\tgi:f:{:.6}\tbi:f:{:.6}\tcg:Z:{}", gi, bi, cigar_str)?;
    for (key, val) in extra_tags {
        write!(out, "\t{}:{}", key, val)?;
    }
    writeln!(out)
}
```

### 标准输出标签

| Tag | 类型 | 含义 | pgr 是否需要 |
|-----|------|------|-------------|
| `gi:f:` | float | gap-compressed identity | ✅ 核心质量指标 |
| `bi:f:` | float | block identity（含 indel bp） | ✅ 核心质量指标 |
| `cg:Z:` | string | CIGAR string | ✅ PAF 标准标签 |
| `an:Z:` | string | alignment name（用于 trace back） | ⚠ 可选 |

### identity 计算公式（impg `main.rs:12042-12061`）

```
gap_compressed_identity = matches / (matches + mismatches + #indel_events)
block_identity         = matches / (matches + mismatches + indel_bp_total)
```

- `gi`（gap-compressed）：每个 indel **事件**计 1 个差异（不计 indel 长度），对长 indel 宽容
- `bi`（block）：每个 indel **碱基**都计入差异，对长 indel 严格
- 两者互补：`gi` 适合评估"同源性"（有无大段 SV），`bi` 适合评估"序列一致性"

具体的 CIGAR 统计逻辑（impg `main.rs:12042-12056`）：

```rust
let (matches, mismatches, insertions, inserted_bp, deletions, deleted_bp, block_len) =
    cigar.iter().fold((0,0,0,0,0,0,0), |(m, mm, i, i_bp, d, d_bp, bl), op| {
        let len = op.len();
        match op.op() {
            'M' => (m + len, mm, i, i_bp, d, d_bp, bl + len), // overestimate
            '=' => (m + len, mm, i, i_bp, d, d_bp, bl + len),
            'X' => (m, mm + len, i, i_bp, d, d_bp, bl + len),
            'I' => (m, mm, i + 1, i_bp + len, d, d_bp, bl + len),
            'D' => (m, mm, i, i_bp, d + 1, d_bp + len, bl + len),
            _   => (m, mm, i, i_bp, d, d_bp, bl),
        }
    });
```

---

## 7. 错误类型设计

参考 impg `paf.rs:14-37` 的 `ParseErr` enum：

```rust
#[derive(Debug)]
pub enum PafParseError {
    NotEnoughFields,               // < 12 tab-separated columns
    InvalidInteger(ParseIntError), // non-integer in numeric column
    InvalidStrand,                 // not '+' or '-'
    InvalidCigarFormat,            // malformed cg:Z: tag
    InvalidFormat(String),         // catch-all for semantic errors
    IoError(io::Error),            // wrapped I/O error
}
```

实现 `Display` + `Error`，不做 `From` 自动转换（避免隐式 error chain，与 impg 一致，
`ParseErr` 也没有实现 `From` trait）。在 `execute` 函数中用
`.map_err(|e| anyhow!("PAF parse: {}", e))` 桥接到 `anyhow::Result`。

---

## 8. 实现优先级

### 第一期 — 支撑 `pgr maf to-paf` + `pgr paf index` + `pgr paf query`

| 优先级 | 组件 | 文件 |
|--------|------|------|
| P0 | `PafRecord` struct | `libs/paf/record.rs` |
| P0 | `SequenceIndex` | `libs/seqidx.rs` |
| P0 | `parse_paf_line` + `parse_paf`（纯文本） | `libs/paf/parser.rs` |
| P0 | `parse_cigar` + `format_cigar` | `libs/paf/cigar.rs` |
| P0 | `write_paf_record` | `libs/paf/writer.rs` |
| P0 | `PafIndexBuilder`（建索引） | `libs/paf/index.rs` |
| P0 | `PafIndex::query`（区间投影） | `libs/paf/index.rs` |
| P0 | `PafIndex::query_transitive`（BFS） | `libs/paf/index.rs` |
| P1 | `CigarOp` bit-packing | `libs/paf/cigar.rs` |

### 第二期 — 大 cohort 场景

| 优先级 | 组件 | 说明 |
|--------|------|------|
| P1 | `parse_paf_bgzf` | 模式 2，支持 `.paf.gz` |
| P1 | `read_cigar_data`（懒加载） | 索引层不存 CIGAR，查询时按需读取 |
| P2 | `parse_paf_file`（自动检测） | 统一入口，dispatch 三种模式 |
| P3 | `parse_paf_bgzf_with_gzi` | 模式 3，多线程 BGZF 解压 |

### 第三期 — 查询层增强

| 优先级 | 组件 | 说明 |
|--------|------|------|
| P2 | `PafRecord` 瘦身（8 字段） | 去掉 matches/block_len/mapq，回退到懒加载 |
| P2 | Caf 后处理过滤参数 | `--min-degree`、`--min-chain-length`、`--end-trim` |
| P3 | Chain/Net syntenic 过滤器 | `--syntenic-filter` 参数 |

---

## 9. 新增依赖评估

| 依赖 | 用途 | pgr 现状 |
|------|------|----------|
| `noodles-bgzf` | BGZF 压缩 PAF 的读写 | ✅ 已有（`Cargo.toml:45`） |
| `fxhash` | `SequenceIndex` 的 HashMap hasher | ✅ 已有（`Cargo.toml:34`） |
| `coitrees` | 区间树——PAF 索引的物理基础 | ❌ 需新增 |

**`coitrees` 是唯一需新增的依赖**。pgr 已有 noodles-bgzf、fxhash、rayon 等所有其他基础设施。
`coitrees` 是轻量 crate（无额外依赖链），符合 CLAUDE.md 的"简洁优先"原则。

如果不引入 `coitrees`，可以用 `intspan`（pgr 已有，`Cargo.toml:40`）的 `RangeMap`
作为替代——但 `intspan` 是区间集合（存储占用区间），不是区间树（区间重叠查询），
语义不匹配。更合适的替代是用 `std::collections::BTreeMap` 手写区间索引——零依赖，
但需要自己实现区间重叠查询和平衡。第一期建议直接引入 `coitrees`，后续可评估手写替代。

---

## 10. 与 pgr 现有基础设施的关系：`paf.rs` vs `loc.rs`

pgr 已有的 `src/libs/loc.rs`（202 行）是 FASTA 随机访问索引模块，与 impg 的 `paf.rs`
（417 行）在架构上有显著的平行对应。本节分析两者关系，明确哪些可以复用、哪些需要新增。

### 10.1 架构对比

| 维度 | impg `paf.rs` | pgr `loc.rs` |
|------|--------------|-------------|
| 输入抽象 | `PafHandle { Plain(File), Compressed(bgzf::Reader) }` | `Input { Buf(Box<dyn BufRead>), File(File), Bgzf(bgzf::io::IndexedReader) }` |
| 随机访问 | `read_cigar_data(offset, byte_len)` | `read_offset(offset, size)` |
| BGZF seek | `bgzf::VirtualPosition::from(offset)` → `reader.seek(vpos)` | `bgzf::io::IndexedReader::seek(SeekFrom::Start(offset))` |
| 索引结构 | per-target coitrees 区间树 | `IndexMap<String, (offset, size)>` — 记录级 |
| 索引粒度 | 区间重叠查询 | 记录名查找 |
| 索引用途 | "哪些 PAF 记录与 chr1:1000-5000 重叠？" | "获取序列 chr1 的第 1000-5000 bp" |
| BufRead 支持 | ✅（`parse_paf`） | ✅（`Input::Buf`） |
| BGZF 支持 | ✅（`parse_paf_bgzf`） | ✅（`Input::Bgzf`） |
| GZI 多线程 | ✅（`parse_paf_bgzf_with_gzi`） | ❌（但 `IndexedReader` 自带索引能力） |

### 10.2 可以直接复用的

**`Input` enum**（`loc.rs:7-11`）：

```rust
pub enum Input {
    Buf(Box<dyn BufRead>),
    File(std::fs::File),
    Bgzf(bgzf::io::IndexedReader<std::fs::File>),
}
```

比 impg 的 `PafHandle` 更通用——多了 `Buf`（支持 stdin 和普通 BufRead），
且 `Bgzf` 变体使用 `IndexedReader`（自带 `.gzi` 索引），seek 能力比 impg 的
基础 `bgzf::io::Reader` 更强。pgr 的 PAF 模块可以直接复用 `Input` 来支持
纯文本、BGZF 压缩、stdin 三种输入模式，不需要像 impg 那样自己定义 `PafHandle`。

**`read_offset()`**（`loc.rs:185-201`）— CIGAR 懒加载的物理基础：

```rust
pub fn read_offset(reader: &mut Input, offset: u64, size: usize) -> anyhow::Result<Vec<u8>> {
    let mut data_buf = vec![0; size];
    match reader {
        Input::File(rdr) => {
            rdr.seek(SeekFrom::Start(offset))?;
            rdr.read_exact(&mut data_buf)?;
        }
        Input::Bgzf(rdr) => {
            rdr.seek(SeekFrom::Start(offset))?;  // IndexedReader handles vpos internally
            rdr.read_exact(&mut data_buf)?;
        }
        Input::Buf(_) => unreachable!(),
    }
    Ok(data_buf)
}
```

与 impg `read_cigar_data`（`paf.rs:68-114`）的功能完全对应——都是 seek + read + 返回字节。
关键差异：
- impg 用 `bgzf::VirtualPosition::from(offset)` 构造 vpos 显式 seek
- pgr 的 `IndexedReader` 在内部处理 vpos 转换——调用者只需传字节偏移量
- pgr 的实现更简洁：11 行匹配 + 2 行 seek/read vs impg 的 46 行分支

这意味着 pgr 的 `read_cigar_data` 可以直接写成对 `read_offset` 的薄封装：

```rust
pub fn read_cigar_data(input: &mut Input, paf_path: &str, offset: u64, byte_len: usize)
    -> Result<Vec<CigarOp>, CigarReadError>
{
    let bytes = libs::io::read_offset(input, offset, byte_len)?;
    let cigar_str = std::str::from_utf8(&bytes)?;
    parse_cigar(cigar_str)
}
```

**`reader_buf()`**（`loc.rs:66-78`）— 纯文本 PAF 解析的 reader 创建。

### 10.3 需要增强的

**`Input` enum 不支持 BGZF 行迭代读取**：

`loc.rs` 的 `create_loc()` 通过 `match &mut reader { Input::Bgzf(rdr) => rdr.read_line(...) }`
来做行读取（`loc.rs:30-34`）。但这个模式分散在 `create_loc` 函数内部，没有抽象为通用接口。
PAF 解析需要同样的能力（按行读取 BGZF 压缩文件）来支持 `parse_paf_bgzf`。

建议：在 `Input` 上实现一个 `read_line` 方法，统一三种变体的行读取：

```rust
impl Input {
    pub fn read_line(&mut self, buf: &mut String) -> io::Result<usize> {
        match self {
            Input::Buf(rdr) => rdr.read_line(buf),
            Input::Bgzf(rdr) => rdr.read_line(buf),
            Input::File(_) => unreachable!(), // File variant only used for seek
        }
    }
}
```

**`IndexMap` vs `FxHashMap`**：

`loc.rs` 用 `IndexMap<String, (u64, usize)>`（保序 HashMap），而 impg 用
`FxHashMap<String, u32>`（更快，不保序）。PAF 的 `SequenceIndex` 不需要保序，
应该用 `FxHashMap`（pgr 已有 fxhash 依赖）。

### 10.4 完全缺失、需要新增的

**区间树（interval tree）** — `loc.rs` 不提供区间重叠查询。

`loc.rs` 的索引是记录级（name→offset），不是区间级（interval→records）。PAF 索引
需要"给定 target `chr1:1000-5000`，返回所有与此区间有重叠的 PAF 记录"——这需要
区间树（如 `coitrees`）或等价物。这是唯一无法从 `loc.rs` 复用、必须新增的组件。

**PAF 行解析** — `loc.rs` 只解析 FASTA 格式（`>` 前缀的行）。

**CIGAR 编解码** — 完全空白。`CigarOp` bit-packing、`parse_cigar`/`format_cigar`、
identity 计算全部需要从零实现。

**PGZIP virtual position 跟踪** — `loc.rs` 的 `IndexedReader` 在内部处理 vpos，
不需要调用者关心。但 PAF 的 CIGAR 懒加载需要把 vpos 存入 `PafRecord.strand_and_offset`，
这意味着 BGZF 解析模式下需要显式获取 `reader.virtual_position()`。
pgr 当前没有在任何地方调用此方法——需要在 `Input::Bgzf` 分支上暴露。

### 10.5 分工结论

```
pgr 已有 (loc.rs)              直接复用           需增强             需新增
──────────────────────────────────────────────────────────────────────────
Input enum                     ✅ IO 抽象          read_line 方法      virtual_position 暴露
read_offset()                  ✅ CIGAR 懒加载
reader_buf()                   ✅ 纯文本 reader
IndexMap name→offset           -                   改为 FxHashMap      -
bgzf::io::IndexedReader        ✅ BGZF seek 能力   行迭代读取           -
─                               
区间树                          -                   -                  ✅ coitrees
PAF 行解析                       -                   -                  ✅ parse_paf_line
CIGAR 编解码                     -                   -                  ✅ CigarOp + parse/fmt
SequenceIndex                   -                   -                  ✅ SeqIndex (可借鉴 IndexMap 模式)
PAF 输出格式化                   -                   -                  ✅ write_paf_record
identity 计算                    -                   -                  ✅ gi/bi 公式
```

**核心结论**：pgr 的 `loc.rs` 已经解决了 PAF 模块中最棘手的 IO 问题——
多格式输入抽象 + BGZF 随机访问。这比 impg 的 `paf.rs` 在 IO 层更成熟（`IndexedReader` > 基础 `Reader`）。
PAF 模块真正需要从零写的只有三样：**区间树索引**、**PAF 行解析**、**CIGAR 编解码**。
其余都可以直接复用或薄封装。
