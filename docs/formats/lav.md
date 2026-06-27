# LAV 格式 (LAV Format)

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

