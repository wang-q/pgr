# 支持的文件格式 (Supported File Formats)

本文档描述 `pgr` 支持的基因组比对文件格式及其实现细节。

## LAV 格式 (LAV Format)

LAV (Local Alignment View) 是 BLASTZ 等比对工具使用的格式。`pgr` 在 `src/libs/lav.rs` 中实现了 LAV 文件的解析。

### 结构 (Structure)

LAV 文件由一系列的“节” (stanza) 组成，以 `>:` 开头（但在实际文件中通常以 `#:lav` 开头，节以 `{}` 包裹）。主要包含以下几种节：

*   `s` (size): 定义序列的大小。
*   `h` (header): 定义序列的名称和可能的反向互补标记 `(reverse)`。
*   `d` (data): 包含比对数据的矩阵等信息（`pgr` 目前主要关注比对块）。
*   `a` (alignment): 具体的比对块，包含得分 (`s`), 块坐标 (`b`), 或具体的区间列表 (`l`)。

### 坐标系 (Coordinate System)

LAV 格式使用 **1-based 闭区间** (1-based, fully closed) 坐标系统。

*   序列索引: 从 1 开始。
*   区间表示: `[start, end]`，即包含 `start` 和 `end` 位置的碱基。
*   `a` 节中的坐标:
    *   `s` (score): 比对得分。
    *   `b` (block): `b <begin_t> <begin_q>`。表示比对块在目标序列和查询序列的起始位置（1-based）。
    *   `e` (end): `e <end_t> <end_q>`。表示比对块的结束位置（1-based）。
    *   `l` (location): `l <begin_t> <begin_q> <end_t> <end_q> <percent_identity>`。详细的局部比对片段，包含起始和结束位置。

### `pgr` 实现 (`libs/lav.rs`) 与 UCSC `lav` 库的对比

`pgr` 的 LAV 解析实现 (`src/libs/lav.rs`) 旨在与 UCSC 的 `lav.c` 和 `lavToPsl.c` 保持高度兼容，以确保转换结果的一致性。以下是详细的对比和实现细节：

*   坐标系转换 (Coordinate System Conversion)
    *   LAV 文件: 使用 1-based 闭区间坐标。
    *   内部表示: `pgr` 与 UCSC 一致，在读取时将坐标减 1 (`val - 1`)，转换为 0-based 半开区间 (half-open interval) `[start, end)`。这对于后续转换为 PSL 格式至关重要。

*   得分调整 (Score Adjustment)
    *   逻辑: UCSC 的 `lavToPsl` 在读取 LAV 得分时会执行 `score = score - 1` 操作。
    *   实现: `pgr` 复刻了这一逻辑，确保生成的 PSL 得分与 UCSC 工具一致。

*   序列名称处理 (Sequence Name Processing - `justChrom`)
    *   逻辑: UCSC 在处理 LAV 头部 (`h` stanza) 时，通常会调用 `justChrom` 函数，该函数会移除文件路径（只保留文件名），并移除 `.nib` 或 `.fa` 等特定前缀/后缀，甚至处理 `nib:` 格式的路径。
    *   实现: `pgr` 实现了类似的 `parse_header_word` 逻辑：
        1.  移除可能的引号。
        2.  移除 `>` 前缀。
        3.  移除路径目录部分，只保留文件名。
        4.  检查并处理 `(reverse)` 标记以确定链方向。

*   序列大小解析 (Size Parsing)
    *   格式: `s` 节通常包含 3 个字段：`filename start length`。
    *   实现: `pgr` 严格按照 UCSC 逻辑，解析第三个字段作为序列全长 (`tSize`/`qSize`)。

*   边缘碎片移除 (Frayed Ends Removal)
    *   逻辑: 比对过程中可能会产生长度为 0 的块（虽然少见，但在某些 LAV 文件中存在）。UCSC 有 `removeFrayedEnds` 函数来清理这些块。
    *   实现: `pgr` 实现了 `remove_frayed_ends` 函数，自动检测并移除比对块列表首尾长度为 0 的块，防止生成无效的 PSL 记录。

*   忽略的节 (Ignored Stanzas)
    *   `pgr` 目前主要关注核心比对信息，因此会解析并忽略 `d` (data), `m` (matrix), `x` (gap penalties) 等节，这不会影响核心坐标和序列信息的提取。

## PSL 格式 (PSL Format)

PSL 是 UCSC 使用的标准比对格式。`pgr` 在 `src/libs/psl.rs` 中提供了一个 `Psl` 结构体，匹配 21 列标准。

### 结构与字段 (Structure and Fields)

PSL 文件是制表符分隔的文本文件，通常包含 21 列：

1.  `matches`: 匹配的碱基数（不包括重复序列）。
2.  `misMatches`: 不匹配的碱基数。
3.  `repMatches`: 重复序列中的匹配碱基数。
4.  `nCount`: N 的数量。
5.  `qNumInsert`: 查询序列中的插入数（gap 数量）。
6.  `qBaseInsert`: 查询序列中的插入碱基总数。
7.  `tNumInsert`: 目标序列中的插入数。
8.  `tBaseInsert`: 目标序列中的插入碱基总数。
9.  `strand`: 链的方向 (`+` 或 `-` 表示查询链；对于蛋白质比对等情况可能为 `++` 等）。
10. `qName`: 查询序列名称。
11. `qSize`: 查询序列大小。
12. `qStart`: 查询序列比对起始位置 (0-based inclusive)。
13. `qEnd`: 查询序列比对结束位置 (0-based exclusive)。
14. `tName`: 目标序列名称。
15. `tSize`: 目标序列大小。
16. `tStart`: 目标序列比对起始位置 (0-based inclusive)。
17. `tEnd`: 目标序列比对结束位置 (0-based exclusive)。
18. `blockCount`: 比对块的数量。
19. `blockSizes`: 每个比对块的大小（逗号分隔）。
20. `qStarts`: 每个比对块在查询序列中的起始位置（逗号分隔）。
21. `tStarts`: 每个比对块在目标序列中的起始位置（逗号分隔）。

### 坐标系说明 (Coordinate System)

*   **0-based, half-open**: PSL 格式标准使用 0-based 起始，半开区间 `[start, end)`。
*   负链坐标 (Negative Strand Coordinates):
    *   查询序列 (Query): 如果 strand 是 `-`，`qStart` 和 `qEnd` 以及 `qStarts` 数组中的坐标是相对于查询序列**反向互补**后的正向坐标，还是相对于原始序列的坐标？
    *   UCSC 标准:
        *   `qStart`, `qEnd`: 在 PSL 文件中，如果是 `-` 链，这两个值通常是相对于 **反向互补后** 的序列末尾计数的（即从原始序列末尾向前的坐标），或者说是“反转”的坐标。
        *   但是，在 `psl.c` 的 `pslLoad` 中，如果是 `-` 链，它通常会把坐标转换回正链坐标吗？不，`pslLoad` 只是直接加载文件内容。
        *   关键点: UCSC 的 `psl` 结构体在内存中通常存储的是文件中的原始值。
        *   对于 `-` 链，文件中的 `qStart` 是从查询序列末尾算起的吗？
            *   UCSC FAQ 说明: "If the strand is '-', then the qStart and qEnd fields are coordinates in the reverse complement of the query sequence."
            *   这意味着：`qStart` 是从反向互补序列的 0 开始算的。
            *   换算回原始序列索引：`original_index = qSize - 1 - rc_index` (对于点坐标)。对于区间 `[s, e)`，`orig_start = qSize - rc_end`, `orig_end = qSize - rc_start`。
    *   `pgr` 实现：目前 `Psl` 结构体只是数据的容器，直接存储文件中的值，不进行自动坐标转换，这与 UCSC 的 `psl` 结构体行为一致。

### `pgr` 实现 (`libs/psl.rs`) 与 UCSC `psl` 库的对比

`pgr` 的 PSL 实现 (`src/libs/psl.rs`) 旨在与 UCSC 的 `psl.c`/`psl.h` 兼容。

*   整数类型与字段对应 (Integer Types)
    *   UCSC: `match`, `misMatch` 等为 `unsigned` (32-bit)；`BaseInsert` 为 `int` (signed)。
    *   pgr: `match_count`, `mismatch_count` 等显式使用 `u32`；`BaseInsert` 使用 `i32`。
    *   备注: 名称略有调整以符合 Rust 风格，类型语义保持一致。

*   坐标字段 (Coordinate Fields)
    *   UCSC: `qStart`, `qEnd` 等为 `unsigned`。
    *   pgr: 使用 `i32`。
    *   备注: `i32` 理论范围减半（20亿），但在绝大多数基因组数据下无影响。

*   数组字段 (Array Fields)
    *   UCSC: `blockSizes`, `qStarts`, `tStarts` (unsigned *)。
    *   pgr: `block_sizes`, `q_starts`, `t_starts` (Vec<u32>)。
    *   备注: 动态数组实现一致。

*   Strand
    *   UCSC: `char strand[4]` (通常是 2 字符 + null)。
    *   pgr: `String`。
    *   备注: 语义一致。

*   读写支持 (Read/Write Support)
    *   UCSC: 支持完整的 `pslLoad` 和 `pslWrite`。
    *   pgr: 实现了 `from_str` (对应 `pslLoad`) 和 `write_to` (对应 `pslTabOut`)。支持从字符串解析和写入到 writer。
    *   输出格式: `write_to` 生成标准的制表符分隔格式，数组字段以逗号分隔，与 UCSC 严格一致。

*   功能函数 (Functional Parity)
    *   `pslRc`: `pgr` 实现了 `rc` 方法。**注意**: 由于 `pgr` 的 `Psl` 结构体目前不包含序列字段 (`qSequence`/`tSequence`)，因此 `rc` 只处理坐标和 Strand 的反向互补，不处理序列本身。
    *   `pslScore`: `pgr` 实现了 `score` 方法，计算逻辑与 UCSC 一致 (考虑了 protein 乘数)。
    *   `pslIsProtein`: `pgr` 实现了 `is_protein` 方法，通过检查 block 坐标判定是否为蛋白质比对。
    *   `pslFromAlign`: `pgr` 实现了 `from_align` 构造函数，支持从序列比对构建 PSL 记录（包含 trimIndel 逻辑）。
    *   缺失功能: 尚未实现 `pslCheck`, `pslRecalcMatchCounts` 等校验和重算功能。

*   结构差异 (Structure Differences)
    *   序列字段: UCSC `psl` 结构体包含 `qSequence` 和 `tSequence` 字段存储具体序列；`pgr` 目前未实现这些字段，仅存储坐标和统计信息。

## AXT 格式 (AXT Format)

AXT 是 UCSC 定义的一种用于表示成对基因组比对的简单文本格式，通常由 blastz 等工具生成。`pgr` 在 `src/libs/axt.rs` 中实现了 AXT 文件的读写。

### 结构 (Structure)

AXT 文件由一系列的比对块 (block) 组成，块与块之间通常由空行分隔。每个块包含三行核心信息：

1.  **Summary Line** (Header): 包含比对的元数据。
    *   格式: `id tName tStart tEnd qName qStart qEnd qStrand score`
    *   示例: `0 chr19 3001012 3001075 chr11 70568380 70568443 - 3500`
    *   字段说明:
        *   `id`: 比对编号。
        *   `tName`: 目标序列名称。
        *   `tStart`, `tEnd`: 目标序列比对区域的起始和结束位置。
        *   `qName`: 查询序列名称。
        *   `qStart`, `qEnd`: 查询序列比对区域的起始和结束位置。
        *   `qStrand`: 查询序列链方向 (`+` 或 `-`)。目标序列假设总为 `+`。
        *   `score`: 比对得分。

2.  **Primary Sequence Line**: 目标序列的比对碱基（包含 gap `-`）。
3.  **Aligning Sequence Line**: 查询序列的比对碱基（包含 gap `-`）。

### 坐标系 (Coordinate System)

*   **1-based, fully closed**: AXT 文件标准使用 1-based 起始，闭区间 `[start, end]`。
*   **负链处理**:
    *   如果 `qStrand` 是 `-`，`qStart` 和 `qEnd` 是相对于查询序列**反向互补**后的坐标。
*   **内部表示**:
    *   `pgr` 在读取时将起始坐标减 1 (`val - 1`)，转换为 **0-based 半开区间** (0-based, half-open) `[start, end)` 存储在 `Axt` 结构体中。
    *   写入时会自动加 1 还原为 1-based 格式。

### `pgr` 实现 (`libs/axt.rs`)

`pgr` 提供了 `Axt` 结构体和 `AxtReader` 迭代器。

*   **解析** (Parsing):
    *   严格解析 9 个字段的 Header 行。
    *   验证两行序列长度是否一致。
    *   自动处理空行和以 `#` 开头的注释行。
*   **输出** (Writing):
    *   `write_axt` 函数将内存中的 0-based 坐标转换回 1-based 写入文件。
    *   每个块输出后会追加一个空行，符合标准格式规范。

## 2bit 格式 (2bit Format)

2bit 是 UCSC 开发的一种高效的基因组序列二进制存储格式。`pgr` 在 `src/libs/twobit.rs` 中实现了 2bit 文件的读取。

### 结构 (Structure)

2bit 文件主要由文件头、索引和序列数据组成：

*   Header:
    *   `magic`: 标识文件类型的魔数 (`0x1A412743`)。
    *   `version`: 版本号。版本 0 使用 32 位偏移量，版本 1 使用 64 位偏移量（支持 >4GB 文件）。
    *   `seqCount`: 序列数量。
    *   `reserved`: 保留字段。

*   Index:
    *   包含序列名称（长度 + 字符串）和该序列在文件中的偏移量。

*   Data (Sequence Records):
    *   `dnaSize`: 序列总长度 (bp)。
    *   `nBlockCount` / `nStarts` / `nSizes`: 硬屏蔽 (Hard Masking) 区域，表示 'N' 的位置。
    *   `maskBlockCount` / `maskStarts` / `maskSizes`: 软屏蔽 (Soft Masking) 区域，表示小写字母（通常用于重复序列）的位置。
    *   `packedDna`: 实际的 DNA 数据，每字节存储 4 个碱基。

### 数据编码 (Data Encoding)

*   Packed DNA: 采用 2-bit 编码，T=00, C=01, A=10, G=11。
    *   每个字节包含 4 个碱基，高位在前 (Big-Endian bit order within byte)。
    *   最后一个字节可能包含填充位。

### `pgr` 实现 (`libs/twobit.rs`) 与 UCSC `twoBit` 库的对比

`pgr` 的实现旨在提供高效的随机访问读取能力。

*   版本支持 (Version Support)
    *   UCSC: 支持版本 0 (32-bit) 和版本 1 (64-bit)。
    *   pgr: **同时支持版本 0 (32-bit) 和版本 1 (64-bit)**。写入时默认使用版本 1 以支持大文件。

*   字节序处理 (Endianness)
    *   UCSC: 通过检查魔数自动检测字节序并进行 swap。
    *   pgr: 实现了相同的逻辑，支持在不同字节序机器上读取 2bit 文件。

*   屏蔽处理 (Masking Handling)
    *   UCSC: `twoBitToFa` 工具通过命令行参数控制是否应用屏蔽。
    *   pgr: `read_sequence` 方法提供 `no_mask` 参数。
        *   N-blocks: 总是应用，将对应的碱基替换为 'N'。
        *   Mask-blocks: 如果 `no_mask` 为 `false` (默认)，将对应的碱基转换为小写；否则保持大写。

*   随机访问 (Random Access)
    *   UCSC: 使用 `twoBitOpen` 加载索引，支持 `twoBitReadSeqFrag` 读取片段。
    *   pgr: `open` 时预加载所有序列的偏移量到 `HashMap`，并**保留原始序列顺序**。`read_sequence` 支持 `start`/`end` 参数读取任意片段，只解码必要的字节块，高效处理大文件。

## chain 格式 (Chain Format)

Chain 格式用于描述成对序列比对，允许两个序列同时存在 gap。它是 UCSC Genome Browser 工具集（如 `liftOver`）的核心格式。

### 结构 (Structure)

Chain 文件由一系列的链 (chain) 组成，每条链以 header 行开始，后跟多行数据：

1.  **Header Line**:
    *   格式: `chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id`
    *   示例: `chain 4900 chr19 59128983 + 561234 561300 chr11 1234567 + 23456 23522 1`
    *   字段说明:
        *   `chain`: 固定关键字。
        *   `score`: 比对得分。
        *   `tName`, `tSize`, `tStrand`: 目标序列名称、大小、链方向（总是 `+`）。
        *   `tStart`, `tEnd`: 目标序列比对区域的起始和结束位置。
        *   `qName`, `qSize`, `qStrand`: 查询序列名称、大小、链方向 (`+` 或 `-`)。
        *   `qStart`, `qEnd`: 查询序列比对区域的起始和结束位置。
        *   `id`: 链的唯一标识符。

2.  **Data Lines** (Alignment Blocks):
    *   格式: `size dt dq`
    *   字段说明:
        *   `size`: 无 gap 比对块的长度。
        *   `dt`: 目标序列中，当前块结束到下一个块开始之间的 gap 长度。
        *   `dq`: 查询序列中，当前块结束到下一个块开始之间的 gap 长度。
    *   **注意**: 每一条链的最后一行只包含 `size`，表示该链结束。

### 坐标系 (Coordinate System)

*   **0-based, half-open**: Chain 格式使用 0-based 起始，半开区间 `[start, end)`。这与 Python 切片或 BED 格式一致。
*   **负链处理**:
    *   如果 `qStrand` 是 `-`，`qStart` 和 `qEnd` 是基于**反向互补链**的坐标系。这意味着坐标是从序列末尾向前计算的。
    *   `pgr` 在内部处理时会保留这种表示，但在转换为 Block（绝对坐标）时会根据链方向正确计算。

### `pgr` 实现 (`libs/chain/record.rs`)

`pgr` 提供了完整的 Chain 格式读写支持，核心实现在 `src/libs/chain/record.rs` (原 `chaining` 模块)。

*   **数据结构**:
    *   `ChainHeader`: 对应 header 行的所有字段。
    *   `ChainData`: 对应数据行的 `size`, `dt`, `dq`。
    *   `Chain`: 包含一个 header 和一个 `ChainData` 向量。

*   **坐标转换**:
    *   `to_blocks()`: 将差异编码 (`size`/`dt`/`dq`) 转换为绝对坐标的 `Block` 列表，方便后续处理（如 `net` 转换）。
    *   `from_blocks()`: 逆向操作，将排序好的 `Block` 列表转换为 Chain 格式的数据行。

*   **算法移植**:
    *   `pgr` 的 `chaining` 模块（`src/libs/chain`）实现了 UCSC `axtChain` 的核心算法，包括 KD-tree 索引构建、动态规划打分以及重叠区域微调 (`trim_overlaps`)，确保生成的 chain 文件与 UCSC 标准工具兼容。

## Net 格式 (Net Format)

Net 格式用于以分层结构表示最佳比对链，通常由 `chainNet` 工具从 chain 文件生成。这种格式对于 synteny 映射和 `liftOver` 链选择至关重要。

### 结构 (Structure)

Net 文件是一种分层结构，展示了在目标基因组上如何用不同层次的比对链来填充空隙。

*   **Hierarchy**: `net` -> `fill` -> `gap` -> `fill` -> `gap` ...
    *   `net`: 染色体级别的容器。
    *   `fill`: 使用特定的 Chain 来填充某个区域。
    *   `gap`: `fill` 内部未被比对覆盖的空隙，可能会被下一层级的 `fill` 填充。

### 字段说明 (Fields)

1.  **Net Line**:
    *   格式: `net tName tSize`
    *   示例: `net chr1 248956422`
    *   说明: 定义目标序列名称和大小。

2.  **Fill Line**:
    *   格式: `fill tStart tLength qName qStrand qStart qLength id chainId score ali qDup type ...`
    *   示例: `  fill 12345 100 chr2 + 67890 100 id 1 score 5000 ali 95`
    *   关键字段:
        *   `tStart`, `tLength`: 目标序列上的填充区域起始和长度。
        *   `qName`, `qStrand`: 查询序列名称和方向。
        *   `qStart`, `qLength`: 查询序列上的对应区域。
        *   `chainId`: 来源 Chain 的 ID。
        *   `score`: 该区域的得分。
        *   `ali`: 对齐的碱基数。

3.  **Gap Line**:
    *   格式: `gap tStart tLength qName qStrand qStart qLength`
    *   示例: `    gap 12400 20 chr2 + 67945 20`
    *   说明: 定义当前 `fill` 内部的一个空隙，该空隙可能会包含子级 `fill`。

### 坐标系 (Coordinate System)

*   **0-based, half-open**: 使用 `start` 和 `length` 表示。
*   **负链处理**:
    *   目标序列 (Target) 始终为正链 (`+`)。
    *   如果 `qStrand` 为 `-`，则 `qStart` 是基于反向互补链的坐标。

### `pgr` 实现 (`libs/net.rs`)

`pgr` 在 `src/libs/net.rs` 中实现了 Net 格式的读写和操作。

*   **数据结构**:
    *   `Chrom`: 对应 `net` 行，包含根 `Gap` 和空间索引。
    *   `Fill`: 对应 `fill` 行，包含比对信息和子级 `Gap` 列表。
    *   `Gap`: 对应 `gap` 行，包含子级 `Fill` 列表。
    *   `NetNode`: 枚举类型，用于统一处理 `Fill` 和 `Gap` 节点。

*   **与 UCSC 对比**:
    *   结构完全对应 UCSC `chainNet.c` 中的 `chrom`, `fill`, `gap` 结构体。
    *   `pgr` 实现了 `read_nets` 解析器，能够还原文件的层级结构。
    *   实现了 `ChainNet` 类，支持将 Chain 添加到 Net 中，复刻了 UCSC `chainNet` 的核心逻辑。
    *   **参数对应**: `add_chain` 方法支持 `min_space` (默认 25) 和 `min_fill` (通常为 min_space/2) 参数，与 UCSC `chainNet` 工具的行为一致，用于控制填充空隙的最小阈值。

## loc 格式 (Location Index)

`.loc` 是 `pgr` 为解决传统 FASTA 索引（如 `.fai`）局限性而设计的索引格式。

*   **格式无关 (Format Agnostic)**: 与 `samtools faidx` 依赖固定行宽不同，`.loc` 格式能够索引**任何**合法的 FASTA 文件，无论其行宽是否固定或一致。这使得它能够可靠地处理草图组装（draft assemblies）和格式混乱的文件。
*   **统一架构 (Unified Architecture)**: 对纯文本和 BGZF 压缩数据采用一致的索引策略，屏蔽了底层压缩细节，实现了对压缩数据的无缝随机访问。

该格式没有对应的 C 实现。

### 结构 (Structure)

`.loc` 文件是制表符分隔的文本文件，每一行对应 FASTA 文件中的一条序列记录。

*   **列定义**:
    1.  `name`: 序列名称（FASTA header 中第一个空格前的字符串）。
    2.  `offset`: 记录在文件中的起始字节偏移量（包含 `>` 字符）。
    3.  `size`: 记录的总字节大小（包含 header、序列和换行符）。

### 示例 (Example)

```text
k81_130	0	272
k81_88	272	428
k81_68	700	241
...
```

### `pgr` 实现 (`libs/loc.rs`)

`pgr` 在 `src/libs/loc.rs` 中实现了 `.loc` 索引的创建和读取。

*   **创建索引**: `create_loc` 函数遍历 FASTA 文件，记录每条序列的偏移量和大小。
*   **读取记录**: `fetch_record` 函数利用索引中的 `offset` 和 `size`，通过 `seek` 直接定位并读取完整记录，无需解析整个文件。
*   **BGZF 支持**: 借助 `noodles_bgzf`，该格式同样适用于 BGZF 压缩的 FASTA 文件，实现压缩数据的随机访问。
