# AXT 格式 (AXT Format)

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

