# pbit 设计笔记（移植 AGC 算法）

## 目标

借鉴 AGC (Assembled Genomes Compressor) v3.2.3 的 C++ 算法（LZ-diff、段级参考压缩、k-mer minimizer
参考选择），设计 pgr 原生的群体基因组压缩格式 `pbit`（p = population 或 plus；2bit 记录 + delta），
集成到 pgr 作为 `pgr pbit` 命令族。**不兼容** C++ AGC 的 `.agc` 文件格式。

**关键约束**：

- **不兼容 `.agc` 文件格式**：设计原生 `pbit` 格式（扩展名 `.pbit`），仅借鉴 AGC 的算法 （LZ-diff、
  段级参考压缩、k-mer minimizer 参考选择），不移植 CArchive 容器 / varint / 前缀编码
- 不引入 zstd 依赖，使用已有的 `flate2`（gzip/DEFLATE）替代
- **参考层直接复用 2bit 记录格式**（`dna_size + n_blocks + mask_blocks + packed_dna`）， 保留 N
  blocks 和 mask blocks（AGC 丢失了 mask），提取 `TwoBitFile::read_sequence` 核心逻辑为 共享函数
  （见 §复用-2）
- 深度复用 pgr 现有基础设施（`libs/fmt/twobit`、`intspan::Range`、`lru::LruCache`、 `libs/fmt/fa`）
- 遵循 pgr 分层原则：复杂逻辑放 `libs/pbit/`，`cmd_pgr/pbit/` 保持薄壳

## 源码分析（AGC C++）

> **采用范围说明**：本节描述的 AGC 文件格式（CArchive 多流容器、varint、前缀编码、footer 索引）
> **仅作参考**，移植时**不采用**——pbit 格式为原生"2bit + delta"（见 §文件格式规范）。本节
> 真正需要移植的是**算法**：LZ-diff（§LZ-diff 算法详解）、段级参考压缩流程（§压缩算法流程）、k-mer
> minimizer 参考选择。CArchive / collection 元数据 / 前缀编码等容器与编码细节**不移植**。

### 架构总览

```
agc-3.2.3/src/
├── app/                    # CLI 入口
│   ├── main.cpp            # main(), 调用 application.cpp
│   └── application.cpp     # 命令分发 (ketopt 参数解析)
├── core/                   # 核心算法
│   ├── agc_compressor.*    # 压缩器 (CAGCCompressor)
│   ├── agc_decompressor.*  # 解压器 (CAGCDecompressor)
│   ├── genome_io.*         # FASTA 读写 (CGenomeIO)
│   ├── kmer.h              # CKmer: 2-bit k-mer
│   └── hs.h                # 哈希集合
├── common/                 # 通用模块
│   ├── archive.*           # CArchive: 多流容器
│   ├── collection*.*       # CCollection_V1/V2/V3: 元数据
│   ├── segment.*           # CSegment: 段级压缩
│   ├── lz_diff.*           # CLZDiff_V1/V2: LZ-diff 算法
│   ├── agc_basic.*         # CAGCBasic: 压缩/解压基类
│   ├── agc_decompressor_lib.* # 解压库基类
│   ├── defs.h              # 类型定义
│   ├── io.h                # CInFile/COutFile
│   ├── utils.*             # 工具函数
│   └── queue.h             # 并发队列
├── lib-cxx/                # C++ API
└── py_agc_api/             # Python API
```

### CLI 命令（9 个）

C++ AGC 用**位置参数**（非 `-i`/`-r` flag），`-i` 在 create/append 中是"含 FASTA 文件名列表
的文件"（非直接输入 FASTA）：

| 命令      | 功能                        | 关键参数（位置参数 + flag）                                              |
|-----------|-----------------------------|--------------------------------------------------------------------------|
| `create`  | 从 FASTA 创建归档（含参考） | `<ref.fa> [<in1.fa>...]`, `-o out.agc`, `-s seg_size`, `-k kmer_len`, `-b pack_card`, `-l min_match_len` |
| `append`  | 向已有归档追加样本          | `<in.agc> [<in1.fa>...]`, `-o out.agc`                                   |
| `getcol`  | 提取所有样本                | `<in.agc>`, `-o out_dir/`, `-r`(不含参考), `-g`(gzip)                    |
| `getset`  | 提取指定样本                | `<in.agc> <sample_name>`, `-o out.fa`                                    |
| `getctg`  | 提取指定 contig             | `<in.agc> <contig_name>`, `-o out.fa`                                    |
| `listref` | 列出参考样本名              | `<in.agc>`                                                               |
| `listset` | 列出所有样本名              | `<in.agc>`                                                               |
| `listctg` | 列出样本和 contig 名        | `<in.agc>`                                                               |
| `info`    | 显示统计信息                | `<in.agc>`                                                               |

### 压缩算法流程

```
create 流程:
1. 读取参考 FASTA → 分段（~segment_size bp/段）
2. 每段作为一个 "group" 的参考 (group_id 递增, in_group_id=0)
   → C++: zstd 压缩存储到 archive stream "seg.{gid}.ref"
   → pbit: write_2bit_record 写标准 2bit 记录（保留 mask，不二次压缩）
3. 对每个输入样本:
   a. 读取 FASTA → 分段
   b. 对每段:
      - 计算 k-mer minimizer → 在所有 group 的参考中找最佳匹配
      - 若最佳匹配的参考段在反向上更好 → is_rev_comp=true, 反向互补
      - 用 LZ-diff 编码差异 → delta
      - 若 delta 与已有 delta 相同 → 复用 (in_group_id 指向已有)
      - 否则 → 新增到 group (in_group_id++)
   c. 每条 delta 单独 flate2 压缩（支持随机访问单样本，不批打包）
4. 存储元数据 (collection) → C++: archive stream; pbit: flate2 压缩到 Sample Index
5. C++: 存储 file_type_info → archive stream
6. 序列化 footer
```

> **pbit 与 C++ 的关键差异**：group = 参考段（与 C++ 一致），但 pbit 的 Reference Index 按
> `contig_name` 分组记录各段偏移（供 `SequenceReader` 按 contig 名拼接多段），C++ 无此需求（按
> k-mer splitter 对索引 group）。

### LZ-diff 算法详解

**核心思想**：LZ77 变体，在参考序列上建哈希表，用 (位置差, 长度) 编码匹配，未匹配部分为 literal。

**数据结构**：

- `reference`: 2-bit 编码的参考序列（A=0,C=1,G=2,T=3,N=4，其他=31）
- `ht16` / `ht32`: 开放寻址哈希表，存储参考中 k-mer 的位置
    - `short_ht_ver`: 参考长度 < 65535×hashing_step 时用 16-bit
    - `USE_SPARSE_HT`: 每 4 位取一个 key（hashing_step=4），减少表大小
    - `max_load_factor=0.7`, `max_no_tries=64`（线性探测上限）
- `key_len = min_match_len - hashing_step + 1`（默认 min_match_len=18 → key_len=15）
- `key_mask`: 2×key_len 位的掩码

**编码格式**（V2，当前版本）：

- **Literal**: `'A' + code`（单字节，code 是 0-20 的 2-bit 值）
- **特殊 literal `'!'`**: 表示 "与参考同位置相同"，解码时取 `reference[pred_pos]`
- **Match**: `<diff_pos>,<len-min_match_len>.` 或 `<diff_pos>.`（到序列末尾的匹配，len=~0u）
    - `diff_pos = ref_pos - pred_pos`（有符号，ASCII 十进制）
- **N-run**: `N_run_starter_code(30)` + `<len-min_Nrun_len>` + `N_code(4)`（≥4 个连续 N）

**V1 vs V2 差异**：

- V2 增加 "equal sequences" 优化（delta 为空）
- V2 增加 "match to end" 优化（len=~0u 省略长度字段）
- V2 增加 `'!'` back-reference literal
- V2 增加 `get_code_skip1` 快速扫描（在前一个 key 有效且当前有 literal 时，滑动窗口而非重新计算）

**编码流程**（`CLZDiff_V2::Encode`）：

```
i = 0, pred_pos = 0
while i + key_len < text_size:
    x = get_code(text + i)          # 提取 key_len 个 2-bit 碱基为 uint64
    if x == ~0u:                     # 含非 ACGT 字符
        Nrun_len = get_Nrun_len(...) # 检测 N-run
        if Nrun_len >= 4:
            encode_Nrun(Nrun_len)    # 编码 N-run
            i += Nrun_len
        else:
            encode_literal(text[i])  # 编码单个 literal
            i++, pred_pos++
        continue

    ht_pos = hash(x) & ht_mask
    find_best_match(ht_pos, text+i, ...)  # 在哈希表中找最佳匹配
    if no match:
        encode_literal(text[i])
        i++, pred_pos++
    else:
        if len_bck > 0:              # 回溯替换之前的 literal
            pop len_bck 个 literal
            match_pos -= len_bck
        if match_pos == pred_pos:    # 同位置匹配，标记 '!'
            替换匹配的 literal 为 '!'
        encode_match(match_pos, len, pred_pos)
        i += len, pred_pos = match_pos + len
```

### 文件格式

#### CArchive 多流容器

**footer-based 设计**（无 magic number）：

```
┌─────────────────────────────────────┐  ← 文件起始
│ Stream 0, Part 0 data              │
│ Stream 0, Part 1 data              │
│ ...                                │
│ Stream 1, Part 0 data              │
│ ...                                │
├─────────────────────────────────────┤
│ Footer                              │  ← footer_offset 处
│  ├─ no_streams (varint)            │
│  ├─ for each stream:               │
│  │   ├─ stream_name (null-term str) │
│  │   ├─ cur_id (varint)            │
│  │   ├─ raw_size (varint)          │
│  │   ├─ packed_size (varint)       │
│  │   ├─ packed_data_size (varint)  │
│  │   ├─ no_parts (varint)          │
│  │   └─ for each part:             │
│  │       ├─ offset (fixed uint64)  │
│  │       └─ size (fixed uint64)    │
├─────────────────────────────────────┤
│ footer_size (fixed uint64 LE)       │  ← 文件末 8 字节
└─────────────────────────────────────┘
```

**变长整数编码**（`CArchive::write<T>`）：

- 第 1 字节: 值的字节数 N
- 后续 N 字节: 值的大端表示

**固定整数**（`write_fixed<T>`）：8 字节小端 uint64

#### Stream 命名约定

- `file_type_info` — 归档元数据（producer, version, comment）
- `seg.{group_id}.ref` — 参考序列段（每 group 一个 stream，1 part）
- `seg.{group_id}.delta` — delta 编码段（每 group 一个 stream，多 part）
- `collection` — 样本/contig/segment 元数据

#### 元数据结构（collection.h）

```cpp
struct segment_desc_t {
    uint32_t group_id;       // 属于哪个参考组
    uint32_t in_group_id;    // 组内序号（用于定位 delta part）
    bool is_rev_comp;        // 是否反向互补
    uint32_t raw_length;     // 原始长度
};

// sample_desc_t = Vec<(contig_name, Vec<segment_desc_t>)>
// 每个 sample 包含多个 contig，每个 contig 包含多个 segment_desc
```

**元数据序列化**（V3）：

- 使用前缀编码变长 uint32（1-5 字节，类似 UTF-8 编码）
- 整个元数据块用 zstd 压缩后存入 `collection` stream

### 版本演进

| 版本 | file_version | LZ-diff | 元数据        | 说明                           |
|------|--------------|---------|---------------|--------------------------------|
| V1   | <2000        | V1      | collection_v1 | 初始版本                       |
| V2   | 2000-2999    | V1      | collection_v2 | 改进元数据                     |
| V3   | 3000+        | V2      | collection_v3 | 改进 LZ-diff + zstd 压缩元数据 |

当前版本: `AGC_FILE_MAJOR=3, AGC_FILE_MINOR=0` → 3000

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
pub fn read_2bit_record<R: Read>(
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
pgr 读取固定大小小端整数的模式。pbit 的 header / index / footer 直接复用这些函数（pbit 统一 小端，
`is_swapped` 参数固定 `false`）。

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
// format.rs — Header/Footer + binary helpers (reuse twobit's read_u32 / read_u64 / read_u32_vec)
pub const PBIT_MAGIC: u32 = 0x54494250;  // 'PBIT' — native "2bit + delta" format
pub const PBIT_VERSION_MAJOR: u32 = 1;
pub const PBIT_VERSION_MINOR: u32 = 0;

/// File header (fixed 32 bytes, at file start).
pub struct PbitHeader {
    pub magic: u32,           // PBIT_MAGIC
    pub version: u32,         // major*1000 + minor
    pub segment_size: u32,
    pub kmer_len: u32,
    pub ref_group_count: u32,
    pub sample_count: u32,
    pub ref_records_offset: u64,  // offset to Reference Records (usually 32)
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
    pub packed_size: u32,           // size after flate2
    pub packed_data: Vec<u8>,       // flate2(LZ-diff encoded)
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
    pub fn add(&mut self, seq: &[u8]) -> anyhow::Result<Vec<u8>>;
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
    segments: Vec<Segment>,                 // one per ref group
    segment_size: usize,
    kmer_len: usize,
    min_match_len: u32,                     // LZ-diff min match length (default 18)
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
[twobit/mod.rs](file:///Volumes/ExtHome/Scripts/pgr/src/cmd_pgr/twobit/mod.rs)），pbit 命令 设计如下
（括注对应的 C++ AGC 命令）：

- `create` → `create`（build 组，无 twobit 对应）
- `append` → `append`（build 组，无 twobit 对应）
- `getcol` → `to-fa`（transform 组，对应 `2bit to-fa`）
- `getset` → `some`（subset 组，对应 `2bit some`）
- `getctg` → `range`（subset 组，对应 `2bit range`）
- `info` + `listref` + `listset` + `listctg` → `stat`（info 组，对应 `2bit size`，用 flag 区分输出）

```
pgr pbit create    -i input.fa [-i input2.fa ...] -o out.pbit -r ref.fa [-s segment_size] [-k kmer_len] [-l min_match_len]
pgr pbit append    in.pbit -i input.fa [-i input2.fa ...] [-o out.pbit]  # 归档为位置参数输入，-o 省略时原地修改
pgr pbit to-fa     -i in.pbit -o out_dir/          # 提取所有样本为 FASTA，每样本一个文件 out_dir/{sample}.fa
pgr pbit some      -i in.pbit sample_list.txt [-o out.fa] [--invert]  # 按样本名列表提取，输出多 FASTA
pgr pbit range     -i in.pbit "chr1" "chr2:1-1000" [-o out.fa]  # 按 contig/区间提取（遍历所有样本，输出多 FASTA）
pgr pbit stat      -i in.pbit [--samples | --refs | --contigs [-s sample]]  # 统计/列表
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
> - `--contigs [-s sample]`：无 `-s` 列出所有样本的所有 contig 名（对应 C++ `listctg`）；
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
│ Header (固定 32 字节)               │
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

### Header（32 字节，文件起始）

```
offset  size  field              说明
0       4     magic              0x54494250 ('PBIT', 小端)
4       4     version            major*1000 + minor (当前 1000)
8       4     segment_size       分段大小（bp，如 4096）
12      4     kmer_len           k-mer 长度（如 15）
16      4     ref_group_count    参考段总数（每段 = 一个 group，非 contig 数）
20      4     sample_count       样本数（不含参考）
24      8     ref_records_offset Reference Records 起始偏移（通常 = 32）
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

## 实施阶段

### Phase 0: 共享函数提取（前置）

- [x] `libs/fmt/twobit.rs`: 提取 `read_2bit_record` / `write_2bit_record` 为 `pub` 模块级函数
- [x] `libs/fmt/twobit.rs`: `TwoBitFile::read_sequence` 改为 `seek(offset) → read_2bit_record` 薄壳
- [x] `libs/fmt/twobit.rs`: `TwoBitWriter::write` 改为循环调用 `write_2bit_record`
- [x] `libs/fmt/twobit.rs`: `read_u32` / `read_u64` / `read_u32_vec` 改为 `pub`（供 pbit 复用，
  pbit 统一小端，调用时 `is_swapped` 固定传 `false`；见 §复用-2.5）
- [x] `libs/fmt/twobit.rs`: 为 `read_2bit_record` / `write_2bit_record` 添加独立单元测试
  （含 mask/N-block/区间切片往返），不依赖 `TwoBitFile`/`TwoBitWriter` 薄壳
- [x] 验证：现有 `2bit` 命令的 `cargo test` 全部通过（行为不变）

### Phase 1: 文件格式 I/O

- [x] `format.rs`: `PBIT_MAGIC` / 版本常量、`PbitHeader` / `PbitFooter` 结构
- [x] `format.rs`: Header 读写（固定 32 字节）、Footer 读写（固定 24 字节）
- [x] `format.rs`: Reference Index 读写（`RefGroupEntry` 列表 + 段偏移）
- [x] `format.rs`: `DeltaEntry` 头部读写（`is_rev_comp` 1B + `raw_length` 4B + `packed_size` 4B，
  共 9 字节；`packed_data` 体由 Phase 5 填充，本阶段只定义头部格式，见 §文件格式规范-Delta Data）
- [x] `format.rs`: Sample Index 容器读写（`ref_group_count` / `delta_count` 计数字段 + 顺序结构，
  供 Phase 5 顺序扫描构建 `delta_offsets`）
- [x] `format.rs`: 字符串读写（u32 len + UTF-8 bytes）
- [x] 验证：能读写一个仅含 Header + 空 Reference Records + Footer 的空 `.pbit` 文件，往返一致

### Phase 2: LZ-diff 算法

- [x] `lz_diff.rs`: 2-bit 编码/解码（ACGT→0123, N→4）
- [x] `lz_diff.rs`: `prepare`（2-bit 编码 + 存储参考 + 填充 key_len 个 invalid_symbol）
- [x] `lz_diff.rs`: `prepare_index`（构建哈希表，16/32-bit 自适应，sparse 模式；仅 encode 需要，
  decode 跳过以节省开销）
- [x] `lz_diff.rs`: `encode`（V2：match/literal/N-run + '!' back-ref + end-match 优化）
- [x] `lz_diff.rs`: `decode`（V2：用 `self.reference` 解析编码 + 重建序列，不传外部 reference）
- [x] `lz_diff.rs`: `estimate`（仅计算编码大小，不生成输出）
- [x] 验证：编码后解码，序列与原始一致

### Phase 3: 段级 delta 压缩

- [ ] `segment.rs`: `Segment::prepare`（接收参考 DNA，委托 `LzDiff::prepare`）
- [ ] `segment.rs`: `Segment::prepare_index`（委托 `LzDiff::prepare_index`，encode 前调用）
- [ ] `segment.rs`: `Segment::add`（LZ-diff 编码 → 返回未压缩 delta bytes）
- [ ] `segment.rs`: `Segment::get`（LZ-diff 解码未压缩 delta → 返回原始序列）
- [ ] 验证：编码 N 条序列后能完整解码（flate2 压缩/解压由 `Compressor`/`Decompressor` 负责，
  `Segment` 仅处理未压缩 delta）

### Phase 4: 元数据

- [ ] `collection.rs`: `SegmentDesc` / `ContigSegs` / `Collection` 结构
- [ ] `collection.rs`: `register_sample_contig` / `add_segment`
- [ ] `collection.rs`: `serialize` / `deserialize`（固定 u32 LE 字段 + flate2）
- [ ] `collection.rs`: 查询方法（`get_contig_segments`, `list_samples`, `list_contigs`）
- [ ] 验证：元数据序列化/反序列化往返一致

### Phase 5: 压缩器/解压器

- [ ] `compressor.rs`: `Compressor<W: Write + Seek>` 泛型结构（直接持有 `W`，无 archive 包装）
- [ ] `compressor.rs`: `Compressor::create`（写 Header 占位 → 读参考 FASTA → 分段 → 每段
  `write_2bit_record` 写入 → 记录 `ref_groups[i].segment_offset` → 对每段创建 `Segment` 并调
  `prepare` + `prepare_index`，供后续 `append_sample` 编码使用）
- [ ] `compressor.rs`: `Compressor::open_for_append`（**从已有 `.pbit` 恢复 Compressor 状态**，供
  `pgr pbit append` 使用：用 `Decompressor` 读取已有归档的 Header / ref_groups / deltas（含
  packed_data）/ collection → 读取各参考段（`read_2bit_record`）→ 为每段创建 `Segment` 并调
  `prepare` + `prepare_index` → 重建 `Compressor` 的完整内存状态 → writer 定位到
  `footer.ref_index_offset`（Reference Index 起始处）。`finish` 时重写 Reference Index + Delta Data
  + Sample Index + Footer。参考 FASTA 不再需要，直接复用归档内已嵌入的参考段）
- [ ] `compressor.rs`: `append_sample`（读 FASTA → 分段 → k-mer minimizer 选择参考组 → **判断正向/
  反向匹配**：若反向更优则 `is_rev_comp=true` 并反向互补序列 → `Segment::add` 编码 → flate2 压缩
  → **delta 去重**：与同 ref_group 已有 delta 比对，相同则复用 `delta_id`，否则新增 → 写
  DeltaEntry → 记入 Collection）
- [ ] `compressor.rs`: k-mer minimizer 参考选择算法（含正向/反向匹配评分，取最优方向）
- [ ] `compressor.rs`: `finish`（写 Reference Index → 写 Delta Data → 写 Sample Index → 写 Footer →
  回填 Header 偏移）
- [ ] `decompressor.rs`: `Decompressor<R: Read + Seek>` 泛型结构（直接持有 `R`）
- [ ] `decompressor.rs`: 三构造器 `open` / `open_and_read` / `new`（镜像 twobit.rs:286-360）
- [ ] `decompressor.rs`: `new` 解析 Header + Footer + Reference Index + Sample Index，构建
  `contig_groups`（contig 名 → ref_group_id 列表）+ `contig_set`（`HashSet<String>`，供
  `contains_contig`）+ 顺序扫描 Delta Data 区构建 `delta_meta`（每条 delta 的
  `is_rev_comp`/`raw_length`/`packed_size`）和 `delta_offsets`（按 `packed_size` 累加偏移，
  仅读 9 字节头部，不解压数据）
- [ ] `decompressor.rs`: `contains_contig(name)`（查 `contig_set`，供 `range` 命令判存在）+
  `list_samples()`（委托 `Collection::list_samples`，供 `stat` 命令）+
  `list_contigs(sample)`（委托 `Collection::list_contigs`，供 `stat --contigs -s`）
- [ ] `decompressor.rs`: `impl SequenceReader for Decompressor<R>`（读**参考层**，镜像 twobit.rs:
  500-510 的 seek → read_2bit_record → slice 模式，但经 `contig_groups` 拼接多段；大 contig 按段
  累加长度（用 2bit 记录的 `dna_size`），只读包含 `[start, end)` 的段，避免拼接整条）
- [ ] `decompressor.rs`: `get_sample` / `get_contig(contig, start, end, strand, out)`（遍历 Collection
  样本 → 用 `delta_meta.raw_length` 累加计算各段坐标，仅解码包含 `[start, end)` 的段 → 每段 seek
  参考段 → `read_2bit_record` → 用 `delta_offsets` 定位 delta → seek → flate2 解压 → LZ-diff
  解码 → 拼接 → 切片 → if strand=="-" rev_comp → 写 FASTA）
- [ ] `decompressor.rs`: LRU 缓存参考段（`ref_group_id` → decoded DNA）+ delta 缓存
  （`(ref_group_id, delta_id)` → decoded raw seq，避免重复 flate2 解压 + LZ-diff 解码）
- [ ] 验证：压缩 → 解压 → FASTA 内容一致；`Decompressor` 的 `SequenceReader` 实现能读参考层 序列（供
  chain/net）；`get_contig` 能遍历样本输出多 FASTA；`open_for_append` 追加样本后 `to-fa` 能输出
  含新旧样本的完整 FASTA

### Phase 6: CLI 集成

- [ ] `cmd_pgr/pbit/mod.rs`: 子命令注册（分组：build/info/subset/transform，镜像 `2bit` mod.rs）
- [ ] `cmd_pgr/pbit/create.rs`: `pgr pbit create`（解析 `-r`/`-o`/`-s`/`-k`/`-l` →
  `Compressor::create` → 对每个 `-i` 输入用 `get_basename` 派生样本名 → `append_sample` → `finish`）
- [ ] `cmd_pgr/pbit/append.rs`: `pgr pbit append`（位置参数 = 输入归档，`-o` 可选 → 若指定 `-o` 则
  先复制归档；`Compressor::open_for_append` 打开 → 对每个 `-i` 输入 `append_sample` → `finish`）
- [ ] `cmd_pgr/pbit/to_fa.rs`: `pgr pbit to-fa`（对应 `2bit to-fa`，但输出为目录：
  用 `outdir_arg`，每样本一个文件 `{sample}.fa`；遍历 `list_samples` → `get_sample` 写各 contig）
- [ ] `cmd_pgr/pbit/some.rs`: `pgr pbit some`（对应 `2bit some`，复用 `fa_name_list_arg` + `invert_arg`）
- [ ] `cmd_pgr/pbit/range.rs`: `pgr pbit range` — 与 `twobit/range.rs` 共用区间解析工具 （`ranges_arg`
  / `collect_ranges` / `intspan::Range` / `nt::rev_comp`），但调 `Decompressor::get_contig`
  遍历样本输出多 FASTA（`getctg` 语义，见 §复用-2.6）
- [ ] `cmd_pgr/pbit/stat.rs`: `pgr pbit stat`（合并 C++ `info`/`listref`/`listset`/`listctg`， 用
  `--samples`/`--refs`/`--contigs` flag 区分）
- [ ] 在 `src/pgr.rs` 注册 `pbit` 子命令
- [ ] 验证：`cargo fmt && cargo clippy -- -D warnings && cargo test`

### Phase 7: 测试与基准

#### 7.1 单元测试（`#[cfg(test)] mod tests`，各模块内嵌）

- [ ] `lz_diff.rs`:
  - 编解码往返：随机序列（ACGT-only、含 N、含小写）→ `prepare` + `prepare_index` → `encode` → `decode` → 比对
  - 边界：空序列、1 bp 序列、纯 N 序列、超长序列（> segment_size）
  - `prepare` 不调 `prepare_index` 时 `decode` 仍正常（验证 decode 不依赖哈希表）
  - `min_match_len` 不同值（15/18/21）下编解码正确
- [ ] `segment.rs`:
  - `prepare` + `add` → `get` 往返；多条序列 add 后各自 get 一致
  - 同参考不同样本（相似 vs 差异大）的 delta 大小合理（相似 < 差异）
- [ ] `format.rs`（Header/Footer/Index I/O）:
  - 空 archive 往返：Header（占位偏移）+ 空 Reference Records + Footer（零偏移）
  - 最小 archive 往返：1 ref_group / 0 sample / 0 delta
  - 多 ref_group + 多 sample + 多 delta 往返：偏移回填正确，`Decompressor::new` 能完整解析
  - `delta_meta` 扫描：`is_rev_comp`/`raw_length`/`packed_size` 与写入时一致
- [ ] `collection.rs`:
  - `add_sample` / `add_segment` 后 `list_samples` / `list_contigs` / `get_segments` 返回正确
  - `serde` 序列化/反序列化往返（`bincode` 或自定义二进制，与 Phase 1 一致）

#### 7.2 集成测试（`tests/cli_pbit.rs`，使用 `PgrCmd` 辅助）

遵循 pgr 惯例：单文件 `tests/cli_pbit.rs`，测试数据放 `tests/pbit/`，用 `PgrCmd::new().args(&[...]).run()`。

- [ ] 测试数据准备 — 优先复用现有材料，仅新建必需的派生样本：
  - **复用** [tests/pgr/pseudocat.fa](file:///Volumes/ExtHome/Scripts/pgr/tests/pgr/pseudocat.fa)
    作为主参考（1 contig `cat`，18803 bp，含小写 mask，无 N；> `segment_size`，测多段拼接）
  - **复用** [tests/2bit/expected/testMask.fa](file:///Volumes/ExtHome/Scripts/pgr/tests/2bit/expected/testMask.fa)
    + [tests/2bit/expected/testN.fa](file:///Volumes/ExtHome/Scripts/pgr/tests/2bit/expected/testN.fa)
    作为补充参考（小规模 mask/N 处理往返测试）
  - **复用** [tests/index/final.contigs.fa](file:///Volumes/ExtHome/Scripts/pgr/tests/index/final.contigs.fa)
    作为多 contig 参考（72 contigs，测 Collection 多 contig 路径）
  - **新建** `tests/pbit/sample_cat_snp.fa`：从 `pseudocat.fa` 派生（保持 contig 名 `cat`），
    引入 ~1% SNP + 少量短 indel — 测 delta 压缩有效性（delta 应远小于原始序列）
  - **新建** `tests/pbit/sample_cat_div.fa`：从 `pseudocat.fa` 派生（保持 contig 名 `cat`），
    大段替换/缺失 — 测低相似度退化行为（delta 接近原始大小，不应 panic）
  - **复用** [tests/pgr/pseudopig.fa](file:///Volumes/ExtHome/Scripts/pgr/tests/pgr/pseudopig.fa)
    作为不匹配样本（contig 名 `pig1`/`pig2` ≠ `cat`）— 测不匹配 contig 处理
    （参考 `toy_ex/c.fa`，确认 pbit 是跳过还是单独存储）
  - **新建** `tests/pbit/sample_list.txt`：含 `sample_cat_snp`/`sample_cat_div`，每行一个样本名
- [ ] `test_pbit_create_basic`：`pgr pbit create -r pseudocat.fa -i sample_cat_snp.fa -o out.pbit` → 文件存在、非空
- [ ] `test_pbit_create_multiple_samples`：`-i sample_cat_snp.fa -i sample_cat_div.fa` → `stat --samples` 列出两个样本
- [ ] `test_pbit_stat_overview`：`stat`（无 flag）→ stdout 含 `ref_group_count`/`sample_count`/`segment_size`/`kmer_len`/`min_match_len`
- [ ] `test_pbit_stat_refs`：`stat --refs` → 列出参考 contig 名（`cat`）
- [ ] `test_pbit_stat_samples`：`stat --samples` → 列出 `sample_cat_snp`/`sample_cat_div`
- [ ] `test_pbit_stat_contigs`：`stat --contigs` → 列出所有样本的 contig；`stat --contigs -s sample_cat_snp` → 仅该样本的 contig
- [ ] `test_pbit_to_fa_roundtrip`：`create` → `to-fa -o out_dir/` → 各 `out_dir/{sample}.fa` 与原始样本 FASTA 内容一致（按 contig 比对）
- [ ] `test_pbit_range_full_contig`：`range out.pbit "cat"` → 每个含 `cat` 的样本输出一条 FASTA，序列与 `to-fa` 提取的 `cat` 一致
- [ ] `test_pbit_range_slice`：`range out.pbit "cat:1-100"` → 输出切片序列正确（与 `pgr fa range` 对原始 FASTA 的切片比对）
- [ ] `test_pbit_range_neg_strand`：`range out.pbit "cat(-):1-100"` → 输出为正向切片的反向互补
- [ ] `test_pbit_range_multi_ranges`：多个区间参数 → 每个区间各样本各一条 FASTA
- [ ] `test_pbit_range_multicontig`：用 `final.contigs.fa` 作参考，`range "k81_130" "k81_88:1-50"` → 多 contig 区间提取正确
- [ ] `test_pbit_some_basic`：`some out.pbit sample_list.txt -o out.fa` → 仅含列表中样本的序列
- [ ] `test_pbit_some_invert`：`some out.pbit sample_list.txt --invert -o out.fa` → 仅含列表外样本的序列
- [ ] `test_pbit_append`：`create -r pseudocat.fa -i sample_cat_snp.fa -o out.pbit` → `append out.pbit -i sample_cat_div.fa` → `stat --samples` 列出两个样本
- [ ] `test_pbit_append_overwrite`：`append out.pbit -i sample_cat_div.fa`（省略 `-o`）→ 原地修改，`stat --samples` 正确
- [ ] `test_pbit_create_custom_params`：`-s 1024 -k 10 -l 15` → `stat` 输出对应参数值
- [ ] `test_pbit_empty_contig`：参考含空 contig（长度 0）→ create 不 panic，range 返回空序列
- [ ] `test_pbit_single_sample`：仅 1 个样本 → to-fa/range/some 均正常
- [ ] `test_pbit_identical_samples`：两个相同样本 → delta 去重生效（`stat` 显示 delta 数 ≤ ref_group 数）
- [ ] `test_pbit_no_match_contig`：`create -r pseudocat.fa -i pseudopig.fa` → 不 panic，
  `stat --contigs -s pseudopig` 行为符合设计（跳过或单独存储，取决于 pbit 规格）
- [ ] `test_pbit_mask_roundtrip`：用 `testMask.fa` 作参考 + 样本 → `to-fa` 提取后 mask 大小写保留一致
- [ ] `test_pbit_n_roundtrip`：用 `testN.fa` 作参考 + 样本 → `to-fa` 提取后 N 位置一致

#### 7.3 属性测试 / 随机往返

- [ ] `tests/cli_pbit.rs::test_pbit_random_roundtrip`：用 `rand` 生成随机参考（多 contig，含 N）+ 随机样本（SNP/indel 变异）→ create → to-fa → 比对，重复 N 次
- [ ] 大 contig 分段往返：参考 contig 长度 = `segment_size * 3 + 1`（跨 4 段）→ range 提取各段边界区间 → 序列正确

#### 7.4 性能基准（`benches/pbit_benchmark.rs`，`criterion`）

- [ ] 压缩速度：参考 1 Mb + N 个样本（N=1/10/100），测 `create` 耗时
- [ ] 解压速度：`to-fa` 全量提取 vs `range` 单 contig 提取 vs `range` 区间提取
- [ ] 压缩率：`.pbit` 文件大小 vs 原始 FASTA vs `.2bit`（直接存储参考）→ 三者对比
- [ ] `delta_cache` 命中率：重复 `range` 同一 contig → 第二次应命中缓存
- [ ] 验证：`cargo bench` 能跑通，输出合理数据

## 当前状态

- 规划中 — 尚未开始实现

## 参考资料

- AGC 源码: `agc-3.2.3/src/`
- AGC GitHub: https://github.com/refresh-bio/agc
- AGC 论文: Deorowicz et al., "AGC: Assembly Genomes Compressor", Bioinformatics (2024)
- pgr loc 模块: [libs/loc.rs](file:///Volumes/ExtHome/Scripts/pgr/src/libs/loc.rs)
- pgr 2bit 模块: [libs/fmt/twobit](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/twobit.rs)
- pgr FASTA I/O: [libs/fmt/fa](file:///Volumes/ExtHome/Scripts/pgr/src/libs/fmt/fa.rs)
- pgr fa/range 命令: [cmd_pgr/fa/range.rs](file:///Volumes/ExtHome/Scripts/pgr/src/cmd_pgr/fa/range.rs)
- pgr 2bit/range 命令:
  [cmd_pgr/twobit/range.rs](file:///Volumes/ExtHome/Scripts/pgr/src/cmd_pgr/twobit/range.rs)
- pgr 参数构建器: [cmd_pgr/args.rs](file:///Volumes/ExtHome/Scripts/pgr/src/cmd_pgr/args.rs)
- SPOA 移植笔记（风格参考）:
  [notes/design/spoa_port.md](file:///Volumes/ExtHome/Scripts/pgr/notes/design/spoa_port.md)

