# PSL 格式 (PSL Format)

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

