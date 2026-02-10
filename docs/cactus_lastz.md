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
*   **工具**: `lastz` (CPU)
    *   **注**: 虽然 `FastGA` 等工具速度更快，但 Cactus 目前的代码硬编码了 `lastz` 的调用接口及特定的输出格式（CIGAR/General），暂不支持直接替换。
    *   **PGR 实现**: 在 `pgr` 中，我们在 `c:\Users\wangq\Scripts\pgr\src\cmd_pgr\lav` 下实现了 `lastz` 的包装器。为了兼容后续的链式比对（Chaining）流程，我们将输出格式调整为 **LAV**。
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
            *   **注意**: `[multiple]` 是 Lastz 命令行语法的一部分，**不是文件名的一部分**。磁盘上实际存在的文件仍然是 `target.fa`。Lastz 在读取该文件时会应用此修饰符。
            *   这强制 Lastz 将输入文件识别为包含多条独立序列的文件，从而正确处理跨序列的比对。
        *   **序列名解析 `[nameparse=darkspace]`**:
            *   指定 Lastz 仅使用 FASTA 标题行中**第一个空格之前**的字符串作为序列名称（即截断描述信息）。
            *   例如，对于 `>seq1 description info`，Lastz 会将其识别为 `seq1`。
            *   这确保了输出结果（如 CIGAR 或 General 格式）中的序列名简洁且不含空格，避免破坏列格式解析。
        *   `--querydepth=keep,nowarn:N`: **深度截断**。如果查询序列某位置的比对深度超过 N，Lastz 通常会截断或报错。这里设置为 `keep` (保留所有比对) 但 `nowarn` (不报错)，但通过 `:N` 实际上是告诉 Lastz 关注高深度区域。在 Cactus 语境下，结合后续的 `cactus_covered_intervals`，这是为了确保能捕获高拷贝重复，同时防止极度重复（如着丝粒）导致输出文件过大或运行时间过长。
        *   `--format=lav`: **输出格式**。
            *   为了支持后续的 Chaining 操作（将片段比对串联成完整比对），我们使用 LAV 格式。
            *   LAV 格式保留了详细的比对信息，可以被 `axtChain` 或类似的工具处理。
        *   `--markend`: 标记输出结束，用于完整性检查。

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
        *   `--chunk N --shuffle --seed <S>`: 仅在内存中缓冲 N 条记录，洗牌后写入并清空缓冲区。大幅降低内存峰值。支持确定性随机种子。
        *   `--chunk N` (无 shuffle): 流式处理，内存占用极低。
*   **全 N 过滤**:
    *   自动跳过仅包含 N 的窗口，减少无效计算。
*   **1-based 坐标**:
    *   输出 Header 格式 `>name:start-end` 采用 1-based 闭区间，符合人类阅读习惯及下游工具（如 Samtools）标准。

| 功能点 | `cactus_fasta_fragments.py` | `pgr fa window` (已实现) |
| :--- | :--- | :--- |
| **切片方式** | 滑动窗口 (Fragment/Step) | 滑动窗口 (Length/Step) |
| **输入** | STDIN 流式 | File/STDIN 流式 |
| **Header 格式** | `>name_start` (Origin-1 default) | `>name:start-end` (1-based Range 风格) |
| **过滤** | 跳过全 N | 跳过全 N |
| **随机化** | `--shuffle` (内存密集) | `--shuffle` + `--chunk` (低内存) + `--seed` (可复现) |

### 7.2 用法

详细用法请参考 `pgr fa window --help`。

### 7.3 实现细节
*   **流式处理**: 类似于 `pgr fa size`，逐条读取 Record，不需要 `.loc` 索引文件，适合处理巨大文件流。
*   **内存优化**: 仅持有当前 Record 的 Sequence，不加载整个文件。

### 7.4 Target 构建策略 (PGR Design)

Cactus 采用复杂的 "分块-采样-物理合并" 策略来构建 Target，旨在平衡大基因组的覆盖度与计算资源。对于 `pgr`，我们采用更简化的策略：

*   **小基因组 (Bacteria/Fungi/Insects)**:
    *   直接将**整个基因组文件**作为 Target。
    *   简单高效，无需拆分。
*   **大基因组 (Human/Plants)**:
    *   **按染色体拆分**: 每条染色体（或 Scaffold）作为一个独立的 Target 处理单元。
    *   **逻辑**: 染色体是自然的生物学单元，且单条染色体（如 Human Chr1 ~250Mb）通常已包含足够的重复序列样本，足以进行有效的 RepeatMasking。
    *   **注意 (All-vs-All)**: 无论是哪种拆分策略，为了同时检测**重复序列 (Transposons)** 和 **重复基因 (Gene Duplication/Paralogs)**，必须执行 **All-vs-All** 比对（即 Query 的每一部分都要与 Target 的每一部分比对）。
        *   在按染色体拆分的情况下，这意味着如果我想找 Chr1 上的片段在 Chr2 上的同源拷贝，我需要确保比对空间覆盖了 Chr1 vs Chr2。
        *   但在 Cactus 的 RepeatMasking 阶段，它主要关注**高拷贝重复**（深度 > 50 或更多）。这种高拷贝元件通常在单条染色体内部就有足够多的拷贝数来触发阈值。因此，对于 RepeatMasking 目的，**Self-Alignment (Chr1 vs Chr1)** 通常是足够的。
        *   如果您的目标是包含**低拷贝的重复基因**（例如只有 2-3 个拷贝的旁系同源基因），那么确实需要全基因组范围的 All-vs-All。这会显著增加计算量（N^2 复杂度）。
    *   **优势**:
        *   **天然并行**: 每条染色体一个任务。
        *   **避免边界问题**: 染色体内部连续。
        *   **实现简单**: 无需复杂的随机采样和重组逻辑。

### 7.5 `pgr pl lastz` vs `lastz` (Wrapper)

为了简化上述复杂的 Lastz 调用流程，我们实现了 `pgr pl lastz` 命令：

*   **功能**:
    *   **自动合并 Target**: 如果提供多个 Target 文件（如分染色体的文件），自动合并为临时文件。
    *   **参数预设**:
        *   自动添加 `[multiple][nameparse=darkspace]` 修饰符。
        *   自动计算 `--querydepth=keep,nowarn:N` (根据 `--period`)。
        *   自动设置 `--format=general...` 输出格式。
    *   **参数预设** (`--preset`):
        *   支持 `near`/`primates` (近缘, `set01`), `medium`/`mammals` (中等, `set03`), `distant`/`vertebrates` (远缘, `set06`) 等预设参数集。
        *   支持 `--lastz-args "..."` 透传其他高级参数。
*   **用法示例**:
    ```bash
    # 使用预设参数进行哺乳动物间比对
    pgr pl lastz query.fa target.fa --preset mammals > alignment.cigar
    
    # 手动传递额外参数
    pgr pl lastz query.fa target.fa --lastz-args "--hspthresh=3000" > alignment.cigar
    ```

## 8. 附录：Lastz 命令构建实现参考

以下是 `LastzRepeatMaskJob` 中构建比对命令的详细规范，可供实现参考：

1.  **输入准备**:
    *   **Query**: 上一步生成的切片文件（Fragments）。
    *   **Target**: 将选中的多个 Target Chunks 文件**物理合并 (cat)** 为一个临时文件。
2.  **构建 Lastz 命令行**:
    *   **输入文件修饰符**:
        *   **Target**: `filename[multiple][nameparse=darkspace]`
            *   `[multiple]`: 必选。告诉 Lastz 这是一个包含多条序列的文件（即使合并后只有一条，加上也无妨）。
            *   `[nameparse=darkspace]`: 必选。只取标题行第一个空格前的 ID。
            *   `[unmask]`: 可选。如果 `unmaskInput=True`，则添加此项，忽略输入序列中的软屏蔽（小写字母）。
        *   **Query**: `filename[nameparse=darkspace]`
            *   `[nameparse=darkspace]`: 必选。
            *   `[unmask]`: 可选。同上。
    *   **核心参数**:
        *   `--querydepth=keep,nowarn:<N>`:
            *   `N = period + 3`。其中 `period` 是屏蔽阈值（通常为 10 或根据采样率调整）。
            *   `+3` 是一个经验性的 "Fudge Factor"（修正因子），确保能召回足够多的比对供后续统计。
            *   `keep`: 不丢弃高深度区域的比对。
            *   `nowarn`: 超过深度时不输出警告。
        *   `--format=general:name1,zstart1,end1,name2,zstart2+,end2+`:
            *   输出自定义的表格格式，无 Header。
            *   `zstart`: 0-based start。
            *   `end`: 1-based end (open interval)。
            *   `+`: 强制 Target 坐标始终为正义链（Lastz 默认如果比对到反义链，Start 会大于 End 或使用负坐标，这里强制标准化）。
        *   `--markend`: 在输出文件末尾写入一行标记，防止因进程崩溃导致的文件截断未被发现。
3.  **输出**:
    *   生成 `.cigar` 文件（尽管扩展名是 cigar，实际内容是上述 `general` 格式）。

## 9. 附录：Python 代码详解 (`cactus_lastzRepeatMask.py`)

文件路径: `cactus-master/src/cactus/preprocessor/lastzRepeatMasking/cactus_lastzRepeatMask.py`

### `RepeatMaskOptions` 类
*   **功能**: 数据类，用于存储重复序列屏蔽的配置选项。
*   **关键属性**:
    *   `fragment`: 切片大小（默认 200）。如果为奇数会自动加 1 保证偶数，以便能被 2 整除。
    *   `minPeriod`: 最小重复周期。
    *   `proportionSampled`: 采样比例。
    *   `period`: 实际使用的深度阈值，计算公式 `max(1, round(proportionSampled * minPeriod))`。

### `LastzRepeatMaskJob` 类 (继承自 `RoundedJob`)
这是 Toil Job 的具体实现，负责执行实际的屏蔽任务。

#### `__init__`
*   **功能**: 初始化 Job，计算资源需求。
*   **资源计算**:
    *   **Memory**: 根据 Target 大小动态计算。
    *   **Disk**: 预留 4 倍于输入文件大小的空间。

#### `getFragments(self, fileStore, queryFile)`
*   **功能**: 调用 `cactus_fasta_fragments.py` 对 Query 进行切片。
*   **输入**: 原始 Query FASTA 文件。
*   **输出**: 包含重叠片段的 FASTA 文件路径。
*   **关键调用**:
    ```python
    cactus_call(..., parameters=["cactus_fasta_fragments.py", 
               "--fragment=%s", "--step=%s", "--origin=zero"])
    ```
    *   `--step` 被硬编码为 `fragment // 2` (50% 重叠)。

#### `alignFastaFragments(self, fileStore, targetFiles, fragments)`
*   **功能**: 执行核心的比对步骤。
*   **流程**:
    1.  **合并 Target**: 将所有 Target Chunks 合并为一个临时文件 (`catFiles`)。
    2.  **构建修饰符**: 为 Target 添加 `[multiple][nameparse=darkspace]`，为 Fragments 添加 `[nameparse=darkspace]`。
    3.  **构建命令**: 组装 `lastz` 命令，包含 `--querydepth=keep,nowarn:N` 和 `--format=general`。
    4.  **执行**: 调用 `cactus_call` 运行比对。
*   **输出**: CIGAR (General format) 比对结果文件。

#### `maskCoveredIntervals(self, fileStore, queryFile, alignment)`
*   **功能**: 根据比对结果计算高深度区间并应用屏蔽。
*   **流程**:
    1.  **计算区间**: 调用 C 程序 `cactus_covered_intervals`。
        *   参数 `M`: 深度阈值，计算为 `period * 2` (因为 50% 重叠导致基准深度翻倍)。
        *   参数 `--queryoffsets`: 启用 Query 坐标还原。
    2.  **应用屏蔽**: 调用 `cactus_fasta_softmask_intervals.py`。
        *   读取原始 Query 和上一步生成的区间文件。
        *   生成最终的 Soft-masked FASTA 文件。

#### `run(self, fileStore)`
*   **功能**: Job 的主入口点，串联上述步骤。
*   **流程**:
    1.  从 FileStore 读取 Query 和 Target 文件到本地临时目录。
    2.  调用 `getFragments` 切片。
    3.  调用 `alignFastaFragments` 比对。
    4.  调用 `maskCoveredIntervals` 屏蔽。
    5.  将最终结果写回 FileStore。

