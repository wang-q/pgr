# Cactus Lastz Repeat Masking 模块详解

本文档详细解析 `c:\Users\wangq\Scripts\pgr\cactus-master\preprocessor\lastzRepeatMasking` 目录下的代码逻辑及其在 Cactus 预处理流程中的作用。

## 1. 概述

`LastzRepeatMasking` 是 Cactus 预处理阶段的一个关键模块，用于通过序列比对（通常是自身比对或与近缘物种比对）来识别和屏蔽重复序列。

与传统的 `RepeatMasker`（基于库）不同，Cactus 的这个模块基于“过度比对”（Over-alignment）的原理：如果查询序列的某个区域能比对到目标序列的很多地方（覆盖深度过高），则该区域被认为是重复序列。

## 2. 目录结构与核心组件

该功能模块的代码分布在两个位置：

1.  **编排层 (`src/cactus/preprocessor/lastzRepeatMasking/`)**:
    *   `cactus_lastzRepeatMask.py`: Python (Toil Job)，主控程序，负责调度整个流程。

2.  **工具层 (`preprocessor/lastzRepeatMasking/`)**:
    包含实际执行任务的脚本和 C 程序（通常会被编译/安装到 bin 目录供主控程序调用）：
    *   `cactus_fasta_fragments.py`: **切片工具**。将查询序列切割成重叠的片段。
    *   `cactus_covered_intervals.c`: **深度计算工具** (C Source)。编译后为 `cactus_covered_intervals`。
    *   `cactus_fasta_softmask_intervals.py`: **屏蔽工具**。应用屏蔽区间。
    *   `Makefile`: 构建脚本。

| 文件名 | 类型 | 功能描述 |
| :--- | :--- | :--- |
| `cactus_lastzRepeatMask.py` | Python Job | **调度**。定义了 `LastzRepeatMaskJob`。 |
| `cactus_fasta_fragments.py` | Python Script | **切片**。生成重叠片段。 |
| `cactus_covered_intervals.c` | C Program | **深度计算**。滑动窗口算法核心。 |
| `cactus_fasta_softmask_intervals.py` | Python Script | **屏蔽**。修改 FASTA 序列。 |

## 3. 工作流程 (Workflow)

整个流程在 `LastzRepeatMaskJob.run()` 方法中定义，分为三个主要步骤：

### 步骤 1: 序列切片 (Fragmentation)
*   **方法**: `getFragments`
*   **工具**: `cactus_fasta_fragments.py`
*   **逻辑**:
    *   为了提高灵敏度并允许并行化，将查询序列（Query）切割成短片段。
    *   **参数**:
        *   `--fragment`: 片段长度（默认 200bp）。
        *   `--step`: 步长（默认 `fragment / 2`，即 50% 重叠）。
    *   **思考：为什么选择 200bp？**
        *   **并行粒度与调度**: 避免任务过小（导致调度开销过大）或过大（导致局部高深度截断）。
        *   **内存效率**: 200bp 的 Bitmap 占用内存极低，支持大规模并发。
        *   **边界效应**: 50% 的重叠确保任何 <100bp 的重复元件至少在一个片段中是完整的，避免因切分点导致的检测丢失。
            *   **注意**: 这里的覆盖度是指**还原到原始 Query 全长序列**后的深度。由于相邻片段有 50% 的重叠（例如 Fragment 1 覆盖 0-200bp，Fragment 2 覆盖 100-300bp），原始 Query 的中间区域（如 100-200bp）会被这两个片段分别比对一次。因此，在 `maskCoveredIntervals` 统计全局深度时，这些重叠区域的基准深度就是 2。程序通过 `scale_period = 2` 来校正这个基准值。
            *   **警告**: 这个 50% 的重叠比例 (`step = fragment // 2`) 和深度校正因子 (`scale_period = 2`) 是在代码中**硬编码**的。Cactus 不支持自定义重叠比例。如果强行修改为 80% 重叠（即 5 倍覆盖），但 `scale_period` 仍为 2，会导致严重的**过度屏蔽 (Over-masking)**。
        *   **种子灵敏度**: 足够容纳多个 Lastz 种子（Seed），保证比对灵敏度。
    *   **输出**: 一个包含所有片段的 FASTA 文件。片段名称通常编码了原始坐标（用于后续还原）。

### 步骤 2: 序列比对 (Alignment)
*   **方法**: `alignFastaFragments`
*   **工具**: `lastz` (CPU) 或 `run_kegalign` (GPU)
    *   **注**: 虽然 `FastGA` 等工具速度更快，但 Cactus 目前的代码硬编码了 `lastz`/`run_kegalign` 的调用接口及特定的输出格式（CIGAR/General），暂不支持直接替换。
*   **逻辑**:
    *   将上一步生成的片段（Query Fragments）与目标序列（Target）进行比对。
    *   **注意**: 目标序列通常是该物种所有分块（Chunks）合并后的全长序列，以确保能检测到全基因组范围内的重复。对于小基因组或未分块的情况，它就是完整的原始染色体。
    *   **Query 与 Target 的选择规模**:
        *   **Query**: 每次 Job 处理 **1 个** 基因组分块（Chunk）。
            *   **分块大小**: 由 `chunkSize` 参数控制（默认通常较大，如 10MB+）。
            *   **定义**: 这是指**原始基因组序列**的长度（例如 10MB）。
            *   **注意**: 由于后续的切片（Fragmentation）步骤有 50% 的重叠，实际生成的 Query Fragments 总碱基数会是这个 `chunkSize` 的约 **2 倍**。
        *   **Target**: 取决于 `proportionToSample` 参数（默认 1.0）。
            *   如果是 1.0 (100%)，Target 就是**全基因组**（所有 Chunks 合并）。
            *   如果 < 1.0 (例如 0.2)，Target 则是从全基因组的所有**互斥分块（Non-overlapping Chunks）**中，选取包含当前 Query Chunk 在内、以及其他随机选取的合计约 20% 的分块。
            *   **文件生成**: 程序会将这些**互不重叠**的 Target Chunks **串联 (concatenate)** 成一个巨大的**包含多条序列的 FASTA 文件 (Multi-FASTA)** 供 Lastz 使用。（注：由于 Target Chunks 之间在基因组上是物理分割且互斥的，因此直接串联不会引入人为的重复序列）。
    *   **Lastz 参数与调用细节**:
        *   **Target 输入修饰符 `[multiple]`**: 
            *   默认情况下，Lastz 倾向于将 Target 处理为单条序列。
            *   为了支持上述的 Multi-FASTA Target，Cactus 在调用 Lastz 时会显式添加 `[multiple]` 修饰符（例如 `target.fa[multiple]`）。
            *   这强制 Lastz 将输入文件识别为包含多条独立序列的文件，从而正确处理跨序列的比对。
        *   **序列名解析 `[nameparse=darkspace]`**:
            *   指定 Lastz 仅使用 FASTA 标题行中**第一个空格之前**的字符串作为序列名称（即截断描述信息）。
            *   例如，对于 `>seq1 description info`，Lastz 会将其识别为 `seq1`。
            *   这确保了输出结果（如 CIGAR 或 General 格式）中的序列名简洁且不含空格，避免破坏列格式解析。
        *   `--querydepth=keep,nowarn:N`: **深度截断**。如果查询序列某位置的比对深度超过 N，Lastz 通常会截断或报错。这里设置为 `keep` (保留所有比对) 但 `nowarn` (不报错)，但通过 `:N` 实际上是告诉 Lastz 关注高深度区域。在 Cactus 语境下，结合后续的 `cactus_covered_intervals`，这是为了确保能捕获高拷贝重复，同时防止极度重复（如着丝粒）导致输出文件过大或运行时间过长。
        *   `--format=general:name1,zstart1,end1,name2,zstart2+,end2+`: **输出格式**。
            *   `name1`: Query 片段名称（包含原始坐标信息）。
            *   `zstart1, end1`: Query 上的匹配区间（0-based）。
            *   `name2`: Target 序列名称。
            *   `zstart2+, end2+`: Target 上的匹配区间（正义链坐标）。
            这种简洁的格式只保留了计算覆盖度所需的坐标信息，丢弃了序列本身，极大地减小了 I/O 开销。
        *   `--markend`: 标记输出结束，用于完整性检查。
    *   **GPU 支持**: 如果启用了 GPU，会调用 `run_kegalign`（SegAlign 的一部分）进行加速。

### 步骤 3: 覆盖度计算与屏蔽 (Masking)
*   **方法**: `maskCoveredIntervals`
*   **工具**: `cactus_covered_intervals` + `cactus_fasta_softmask_intervals.py`
*   **逻辑**:
    1.  **计算区间**: `cactus_covered_intervals` 读取比对结果。
        *   **坐标系选择**: 核心目标是屏蔽 **Query** 序列。因此，尽管比对涉及 Target，但深度统计是累加在 Query 坐标轴上的。Target 只是作为“证据”来证明 Query 某些区域在基因组中多次出现。
        *   **坐标还原**: 利用片段名称中的偏移量信息（如 `seq1_1000` 表示偏移 1000bp），将片段内部坐标还原到原始全长 Query 坐标系。
        *   使用滑动窗口算法计算每个碱基的覆盖深度。
        *   如果深度 > 阈值 `M`（由 `period` 参数决定），则输出该区间。

### 思考：为什么不使用 Target 坐标系？

你可能会问，既然 Target 通常是全基因组，统计 Target 上的堆积深度岂不是更直观？这里有几个关键考量：

1.  **局部性与独立性 (Key)**:
    *   **场景**: Cactus 将 Query 切分成无数小片段（如 Fragment A）并行处理。
    *   **Query 视角**: Fragment A 中的一个转座子会比对到 Target 的 100 个位置。在当前 Job 中，我们立刻看到 Fragment A 的该位置深度为 100 -> **判定为重复**。
    *   **Target 视角**: Target 上的这 100 个位置，每个位置只被 Fragment A 覆盖了 **1次**。
        *   要看到 Target 上的高深度，必须汇总 **所有** Query 片段的比对结果（All-vs-All）。
        *   但在分布式的 Job 中，我们只有当前 Fragment A 的信息。
    *   **结论**: 使用 Query 坐标系允许我们在**单个 Job 内部**，仅凭局部信息就能完成重复序列的鉴定，无需全局汇总。

2.  **并行与内存**: 
    *   在内存中维护一个小片段（如 200kb）的覆盖度数组（Bitmap）非常廉价且高效。
    *   如果要统计 Target（全基因组）的深度，每个小 Job 都需要维护巨大的全基因组数组，或者需要一个昂贵的后处理步骤来合并所有 Job 的结果。

3.  **直接性**: 
    *   最终目的是修改 Query 的 FASTA 文件（Soft-masking）。直接计算 Query 上的屏蔽区间，可以直接应用。
    2.  **应用屏蔽**: `cactus_fasta_softmask_intervals.py` 读取原始 Query FASTA 和上一步生成的区间文件。
        *   将指定区间内的序列转换为小写（Soft-masking）。
        *   也可以配置为硬屏蔽（Hard-masking）或反向操作（Unmask）。

## 4. 关键参数 (`RepeatMaskOptions`)

在 `cactus_lastzRepeatMask.py` 中，`RepeatMaskOptions` 类定义了控制流程的关键参数：

*   `fragment` (int, default 200): 切片大小。
*   `minPeriod` (int, default 10): 最小重复周期/覆盖度阈值。
*   `proportionSampled` (float, default 1.0): 采样比例，会影响最终的覆盖度阈值计算 (`period = proportionSampled * minPeriod`)。
*   `lastzOpts` (str): 传递给 `lastz` 的额外参数。
*   `gpu` (int): 指定使用的 GPU 索引。
*   `unmaskInput` / `unmaskOutput`: 控制输入输出的屏蔽状态。

## 5. C 代码细节 (`cactus_covered_intervals.c`)

这是该模块中唯一的 C 代码，负责核心的深度统计算法。

*   **输入**: 标准输入读取 `lastz` 的 `general` 格式输出。
*   **数据结构**: 使用位图（Bitmap）或字节数组作为滑动窗口，记录当前窗口内每个位置的覆盖次数。
*   **特殊处理**:
    *   `--queryoffsets`: 这是一个关键标志。开启后，程序会解析 Query 序列名（如 `>seq1_1000`），将其视为 `seq1` 从 1000bp 开始的片段，从而正确累加深度到全长序列上。

## 6. 总结

`LastzRepeatMasking` 展示了 Cactus "分而治之" 的设计哲学：
1.  **Python** 负责工作流编排和文件管理（Toil）。
2.  **C/C++** (Lastz, cactus_covered_intervals) 负责计算密集型的核心任务。
3.  通过文件流（Pipes/Files）进行组件间通信。

对于 `pgr` 项目，如果需要实现类似的高灵敏度重复序列屏蔽，该模块的架构（切片 -> 比对 -> 深度阈值 -> 屏蔽）是非常值得参考的范例。

## 7. PGR 工具链替代方案设计

为了使用 Rust 高效替代现有的 Python 脚本，计划在 `pgr` 中新增 `fa window` 子命令。

### 7.1 `pgr fa window` vs `cactus_fasta_fragments.py`

我们实现了 `pgr fa window` 来替代 Cactus 的切分逻辑，并增强了内存控制和灵活性。

*   **窗口切分与覆盖度**:
    *   通过 `--len` (窗口长度) 和 `--step` (步长) 控制切分。
    *   Cactus 默认行为: `-l 200 -s 100` (2x 覆盖度，50% 重叠)。
    *   `pgr` 允许任意组合，如 `-l 200 -s 200` (1x 覆盖度，无重叠) 或 `-l 200 -s 10` (20x 覆盖度)。
*   **大文件处理与内存优化**:
    *   **Cactus**: 先生成巨大的完整片段文件，再进行 shuffle 和 split，内存消耗极大。
    *   **pgr**: 引入 `--chunk N` 参数，在生成阶段直接切分输出文件。
        *   `--chunk N --shuffle`: 仅在内存中缓冲 N 条记录，洗牌后写入并清空缓冲区。大幅降低内存峰值。
        *   `--chunk N` (无 shuffle): 流式处理，内存占用极低。
*   **全 N 过滤**:
    *   自动跳过仅包含 N 的窗口，减少无效计算。
*   **1-based 坐标**:
    *   输出 Header 格式 `>name:start-end` 采用 1-based 闭区间，符合人类阅读习惯及下游工具（如 Samtools）标准。

| 功能点 | `cactus_fasta_fragments.py` | `pgr fa window` (设计目标) |
| :--- | :--- | :--- |
| **切片方式** | 滑动窗口 (Fragment/Step) | 滑动窗口 (Length/Step) |
| **输入** | STDIN 流式 | File/STDIN 流式 |
| **Header 格式** | `>name_start` (Origin-1 default) | `>name:start-end` (1-based Range 风格) |
| **过滤** | 跳过全 N | 跳过全 N |
| **随机化** | `--shuffle` (内存密集) | 暂不直接支持内存 Shuffle，建议后续通过 `shuf` 管道处理 |

### 7.2 用法

详细用法请参考 `pgr fa window --help`。

### 7.3 实现细节
*   **流式处理**: 类似于 `pgr fa size`，逐条读取 Record，不需要 `.loc` 索引文件，适合处理巨大文件流。
*   **内存优化**: 仅持有当前 Record 的 Sequence，不加载整个文件。

