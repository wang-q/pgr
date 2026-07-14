# 2bit 格式 (2bit Format)

2bit 是 UCSC 开发的一种高效的基因组序列二进制存储格式。`pgr` 在 `src/libs/fmt/twobit.rs` 中实现了 2bit 文件的读取。

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

### `pgr` 实现 (`libs/fmt/twobit.rs`) 与 UCSC `twoBit` 库的对比

`pgr` 的实现旨在提供高效的随机访问读取能力。

*   版本支持 (Version Support)
    *   UCSC: 支持版本 0 (32-bit) 和版本 1 (64-bit)。
    *   pgr: **同时支持版本 0 (32-bit) 和版本 1 (64-bit)**。写入时默认使用版本 1 以支持大文件。

* 字节序处理 (Endianness)
    * UCSC: 通过检查魔数自动检测字节序并进行 swap。
    * pgr: 写入时始终使用 little-endian；读取时通过 magic 检测自动处理跨平台字节序差异。

*   屏蔽处理 (Masking Handling)
    *   UCSC: `twoBitToFa` 工具通过命令行参数控制是否应用屏蔽。
    *   pgr: `read_sequence` 方法提供 `no_mask` 参数。
        *   N-blocks: 总是应用，将对应的碱基替换为 'N'。
        *   Mask-blocks: 如果 `no_mask` 为 `false` (默认)，将对应的碱基转换为小写；否则保持大写。

*   随机访问 (Random Access)
    *   UCSC: 使用 `twoBitOpen` 加载索引，支持 `twoBitReadSeqFrag` 读取片段。
    *   pgr: `open` 时预加载所有序列的偏移量到 `HashMap`，并**保留原始序列顺序**。`read_sequence` 支持 `start`/`end` 参数读取任意片段，只解码必要的字节块，高效处理大文件。

