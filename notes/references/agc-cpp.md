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
| `listctg` | 列出样本和 contig 名        | `<in.agc> <sample_name>`                                                 |
| `info`    | 显示统计信息                | `<in.agc>`                                                               |

## 压缩算法流程

```
create 流程:
1. 读取参考 FASTA → 分段（~segment_size bp/段）
2. 每段作为一个 "group" 的参考 (group_id 递增, in_group_id=0)
   → C++: ZSTD 压缩存储到 archive stream（V1/V2: `seg-{gid}-ref`, V3: `x{base64(gid)}r`）
     * 先检测序列周期性 best_frac：若 < 0.5 用 bytes2tuples 转换 + ZSTD(level 13, marker=1)，
       否则 plain ZSTD(level 19, marker=0)；压缩后追加 1 字节 marker 供解码判断
   → pbit: write_2bit_record 写标准 2bit 记录（保留 mask，不二次压缩）
3. 对每个输入样本:
   a. 读取 FASTA → 分段
   b. 对每段:
      - 计算 k-mer minimizer → 在所有 group 的参考中找最佳匹配
      - 若最佳匹配的参考段在反向上更好 → is_rev_comp=true, 反向互补
      - 用 LZ-diff 编码差异 → delta
      - 若 delta 为空（与参考相同）→ 复用参考 (in_group_id=0)
      - 若 delta 与已有 delta 相同 → 复用 (in_group_id 指向已有)
      - 否则 → 新增到 group (in_group_id++)
   c. 每 `contigs_in_pack`（`-b` 参数）条 delta 用 0xff 分隔符拼接为一个 part，
      ZSTD 压缩（level 17）后存入 delta stream（支持按 part 随机访问）
4. 存储元数据 (collection) → C++: archive stream; pbit: flate2 压缩到 Sample Index
5. C++: 存储 file_type_info → archive stream
6. 序列化 footer
```

> **pbit 与 C++ 的关键差异**：group = 参考段（与 C++ 一致），但 pbit 的 Reference Index 按
> `contig_name` 分组记录各段偏移（供 `SequenceReader` 按 contig 名拼接多段），C++ 无此需求（按
> k-mer splitter 对索引 group）。

### 段级复用：Identity（delta 为空）与 delta 去重（segment.cpp）

`CSegment::add`（`segment.cpp`）是段级压缩的核心，也是 pbit v1010 Identity 优化的来源：

```
if (no_seqs == 0):
    lz_diff->Prepare(s);  store_in_archive(s);   # 首条 = 该 group 的参考（id 0）
else:
    delta = lz_diff->Encode(s)                    # 对参考做 LZ-diff
    # IMPROVED_LZ_ENCODING（V1 条件启用；V2 无条件启用）
    if delta.empty():      return 0               # ★ 与参考完全相同 → in_group_id=0 = 参考段
    if delta 已在 v_lzp 中: return 已有序号        # 相同 delta 复用
    else: v_lzp.push(delta); return 新序号        # 新增 delta（id 从 1 起）
```

语义要点：
- **Identity 不是单独的类型，而是 `in_group_id = 0`**：`get(id_seq)` 对 `id_seq == 0`
  直接返回参考序列（`segment.cpp:get` 的 `if (id_seq == 0) ctg = ref_seq`），不经过
  LZ-diff 解码。样本段与参考**整段相同**时（delta 为空）零载荷指向参考。
- **delta 去重与 Identity 独立**：相同 delta 复用返回的是 `v_lzp` 中的已有序号（≥1），
  只有空 delta 才指向 0。pbit v1010 用 `DeltaEncoding::Identity` 显式表达空载荷，
  并按 `(encoding, is_rev_comp, raw_length)` 去重共享条目（pbit 的 delta 表按参考组
  存储，无需"id 0 = 参考"的隐式约定）。
- **局限性**：Identity 仅覆盖"段 == 整条参考"（AGC 段通常 ~2 kb）；段与参考**部分相同**
  仍走 LZ-diff match 指令。pbit v1007 的 CIGAR 混合编码把"参考任意子区间"纳入
  CIGAR，因此 pbit v1010 的 Identity 比 AGC 更强：**任意参考区间**全等即零载荷
  （`ref_start/ref_end` 由段描述符携带）。

### delta part 打包（contigs_in_pack，`-b` 参数）

delta 不逐条独立压缩，而是攒 `contigs_in_pack`（`-b pack_cardinality`，默认 128）条，
用 `0xff` 分隔符拼接成一个 part，整块 ZSTD（level 17）写入 delta stream
（`store_in_archive` / `get` 按 `part_id = id_seq / contigs_in_pack` 随机访问，解压后
再按 `0xff` 切出 `seq_in_part_id`）。raw 特殊段（`add_raw`，如空段占位）同机制。

设计动机：单条 delta 通常几十~几百字节，独立压缩要摊薄压缩器头/窗口开销；part 打包
把多条小 delta 合并压缩，压缩率更高，随机访问粒度退化为 part（一次解压 128 条）。

**pbit 借鉴建议（记录，未实现）**：pbit 目前每条 delta 独立 flate2（10 字节头 +
gzip），Identity 后纯 `=` 段零载荷、碎链 PAF 恢复区已整块 flate2；若未来 delta 单条
过小且数量大，可借鉴 part 打包（多条 delta 拼接后统一压缩），代价是随机访问粒度变大
（当前按 delta 独立 seek 简单直接）。与"碎链 cg 位打包"同属细粒度压缩优化，用户已裁定
暂缓，先记录不实现。

## LZ-diff 算法详解

**核心思想**：LZ77 变体，在参考序列上建哈希表，用 (位置差, 长度) 编码匹配，未匹配部分为 literal。

**数据结构**：

- `reference`: 2-bit 编码的参考序列。FASTA 转换表 `cnv_num`（`agc_basic.h`）：A=0,C=1,G=2,T=3,N=4，其他 IUPAC(R/Y/S/W/K/M/B/D/H/V/U)=5-14，无效字符=30。LZ-diff 内 `invalid_symbol=31` 仅用作参考序列尾部 padding 哨兵（`prepare_gen` 追加 key_len 个 31）
- `ht16` / `ht32`: 开放寻址哈希表，存储参考中 k-mer 的位置
    - 哈希函数为 **MurMur64Hash**（`lz_diff.cpp` 内 `MurMur64Hash mmh`，`ht_pos = mmh(x) & ht_mask`），线性探测最多 `max_no_tries=64` 槽；表大小取 2 的幂（`while (ht_size & (ht_size-1)) ht_size &= ht_size-1; ht_size <<= 1`），负载因子 0.7
    - `short_ht_ver`: 参考长度 / hashing_step < 65535 时用 16-bit
    - `USE_SPARSE_HT`: 每 4 位取一个 key（hashing_step=4），减少表大小
- `key_len = min_match_len - hashing_step + 1`（默认 min_match_len=18 → key_len=15）
- `key_mask`: 2×key_len 位的掩码

**编码格式**（V2，当前版本）：

- **Literal**: `'A' + code`（单字节，code 0-20：0-3 为 ACGT，4 为 N，5-14 为其他 IUPAC 码；`is_literal` 判定 `'A'`..`'U'`）
- **特殊 literal `'!'`**: 表示 "与参考同位置相同"，解码时取 `reference[pred_pos]`
- **Match**: `<diff_pos>,<len-min_match_len>.` 或 `<diff_pos>.`（到序列末尾的匹配，len=~0u）
    - `diff_pos = ref_pos - pred_pos`（有符号，ASCII 十进制）
- **N-run**: `N_run_starter_code(30)` + `<len-min_Nrun_len>` + `N_code(4)`（≥4 个连续 N）

**V1 vs V2 差异**：

- "equal sequences" 优化（delta 为空）：V1 由 `IMPROVED_LZ_ENCODING`（`defs.h` 默认定义）条件启用，V2 无条件启用
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
        if match_pos == pred_pos:    # 同位置匹配
            回看之前已编码的 literal，若值 == reference[match_pos - i] 则替换为 '!'
        if 到序列末尾:               # match-to-end 优化
            encode_match(match_pos, ~0u, pred_pos)   # len=~0u 省略长度字段
        else:
            encode_match(match_pos, len, pred_pos)
        i += len, pred_pos = match_pos + len
```

## 文件格式

### CArchive 多流容器

**footer-based 设计**（无 magic number）：

```
┌─────────────────────────────────────┐  ← 文件起始
│ Stream 0, Part 0: metadata(varint) + raw_data
│ Stream 0, Part 1: metadata(varint) + raw_data
│ ...                                │
│ Stream 1, Part 0: ...              │
│ ...                                │
├─────────────────────────────────────┤
│ Footer                              │  ← file_size - 8 - footer_size 处
│  ├─ no_streams (varint)            │
│  ├─ for each stream:               │
│  │   ├─ stream_name (null-term str) │
│  │   ├─ no_parts (varint)          │  ← 读入 cur_id，即 parts 数
│  │   ├─ raw_size (varint)          │
│  │   └─ for each part:             │
│  │       ├─ offset (varint)        │  ← 非 fixed，用 write()
│  │       └─ size (varint)          │  ← 非 fixed，用 write()
│  （注：packed_size / packed_data_size 仅内存字段，不写入 footer）
├─────────────────────────────────────┤
│ footer_size (fixed uint64 LE)       │  ← 文件末 8 字节
└─────────────────────────────────────┘
```

> **数据区每个 part** = `metadata`(varint) + `raw_data`。metadata 由写入方指定
> （如 delta part 存 raw_size=解压后字节数，0 表示未压缩）。

**变长整数编码**（`CArchive::write<T>`）：

- 第 1 字节: 值的字节数 N
- 后续 N 字节: 值的大端表示

**固定整数**（`write_fixed<T>`）：8 字节小端 uint64（`COutFile::WriteUInt` 逐字节低位先行）

### Stream 命名约定

- `file_type_info` — 归档元数据（producer, version, comment）
- `params` — 压缩参数（kmer_length, min_match_len, pack_cardinality, segment_size）
- 参考序列段（每 group 一个 stream，1 part）：
    - V1/V2: `seg-{group_id}-ref`
    - V3: `x{base64(group_id)}r`（base64 编码 group_id，紧凑命名）
- delta 编码段（每 group 一个 stream，多 part）：
    - V1/V2: `seg-{group_id}-delta`
    - V3: `x{base64(group_id)}d`
- collection 元数据（按版本拆分为多个 stream，无单一 `collection` 流）：
    - V1: `collection-desc`
    - V2: `collection-main` + `collection-details`
    - V3: `collection-samples` + `collection-contigs` + `collection-details`

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
- 元数据拆分为 samples / contigs / details 三部分，分别 ZSTD 压缩后存入对应的 `collection-*` stream（见上文 Stream 命名约定）

## 版本演进

| 版本 | file_version | LZ-diff | 元数据        | 说明                           |
|------|--------------|---------|---------------|--------------------------------|
| V1   | <2000        | V1      | collection_v1 | 初始版本                       |
| V2   | 2000-2999    | V1      | collection_v2 | 改进元数据                     |
| V3   | 3000+        | V2      | collection_v3 | 改进 LZ-diff + zstd 压缩元数据 |

当前版本: `AGC_FILE_MAJOR=3, AGC_FILE_MINOR=0` → 3000

> **源码版本怪癖**：目录名为 `agc-3.2.3`，但 `defs.h` 内 `AGC_VER_MAJOR/MINOR/BUGFIX = 3/2/2`
> （build `20260326.1`），即源码自报 v3.2.2；各 `.cpp/.h` 头注释版本号 `3.2`（2024-11-21）。
> 记笔记/对照 pgr 移植时以 `defs.h` 为准。

> **CLI 补充**：`create` 除表格所列外还有 `-a`(adaptive)、`-c`(concatenated genomes)、
> `-d`(不存 cmd-line)、`-f`(fall-back minimizers 比例)、`-t`(线程)。其中 `-a` 自适应压缩
> 与 `-f` 与参考选择相关，pgr 的 Reference Index 设计可参考其"minimizer 找参考 + 局部
> fall-back"思路（细节在 application.cpp / agc_compressor.cpp）。

## 参考资料

- AGC 源码: `agc-3.2.3/src/`
- AGC GitHub: https://github.com/refresh-bio/agc
- AGC 论文: Deorowicz et al., "AGC: Assembly Genomes Compressor", Bioinformatics (2024)
- pbit 设计笔记: [notes/design/pbit.md](../design/pbit.md)
