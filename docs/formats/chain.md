# chain 格式 (Chain Format)

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

