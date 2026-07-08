# AGC (Assembled Genomes Compressor) C++ 源码分析

> 本文档记录 AGC v3.2.3 C++ 源码的分析，作为 pgr `pbit` 格式设计的算法来源参考。
> pbit 仅借鉴 AGC 的算法（LZ-diff、段级参考压缩、k-mer minimizer 参考选择），
> **不兼容** `.agc` 文件格式。pbit 格式规范见 [pbit 设计笔记](../design/pbit.md)。

> **采用范围说明**：本文档描述的 AGC 文件格式（CArchive 多流容器、varint、前缀编码、footer 索引）
> **仅作参考**，移植时**不采用**——pbit 格式为原生"2bit + delta"（见
> [pbit.md §文件格式规范](../design/pbit.md#文件格式规范)）。真正需要移植的是**算法**：
> LZ-diff（§LZ-diff 算法详解）、段级参考压缩流程（§压缩算法流程）、k-mer minimizer 参考选择。
> CArchive / collection 元数据 / 前缀编码等容器与编码细节**不移植**。

## 架构总览

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

## CLI 命令（9 个）

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

## 压缩算法流程

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

## LZ-diff 算法详解

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

## 文件格式

### CArchive 多流容器

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

### Stream 命名约定

- `file_type_info` — 归档元数据（producer, version, comment）
- `seg.{group_id}.ref` — 参考序列段（每 group 一个 stream，1 part）
- `seg.{group_id}.delta` — delta 编码段（每 group 一个 stream，多 part）
- `collection` — 样本/contig/segment 元数据

### 元数据结构（collection.h）

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

## 版本演进

| 版本 | file_version | LZ-diff | 元数据        | 说明                           |
|------|--------------|---------|---------------|--------------------------------|
| V1   | <2000        | V1      | collection_v1 | 初始版本                       |
| V2   | 2000-2999    | V1      | collection_v2 | 改进元数据                     |
| V3   | 3000+        | V2      | collection_v3 | 改进 LZ-diff + zstd 压缩元数据 |

当前版本: `AGC_FILE_MAJOR=3, AGC_FILE_MINOR=0` → 3000

## 参考资料

- AGC 源码: `agc-3.2.3/src/`
- AGC GitHub: https://github.com/refresh-bio/agc
- AGC 论文: Deorowicz et al., "AGC: Assembly Genomes Compressor", Bioinformatics (2024)
- pbit 设计笔记: [notes/design/pbit.md](../design/pbit.md)
