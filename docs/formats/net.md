# Net 格式 (Net Format)

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

