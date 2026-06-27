# loc 格式 (Location Index)

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

