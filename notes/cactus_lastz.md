# Cactus Lastz Repeat Masking 模块详解

本文档详细解析 `src/cactus/preprocessor/lastzRepeatMasking` 目录下的代码逻辑及其在 Cactus 预处理流程中的作用，并介绍了 `pgr` 项目如何通过 Rust 实现高效替代方案。

## 1. 概述与核心概念

`LastzRepeatMasking` 是 Cactus 预处理阶段的关键模块，其核心思想与传统的 `RepeatMasker`（基于重复序列库）完全不同。它基于**“过度比对”（Over-alignment）**的原理：

> 如果基因组中某段序列（Query）能在自身或其他基因组（Target）中找到大量高相似度的比对位置（即覆盖深度极高），那么这段序列很可能就是重复序列（如转座子、简单重复序列）。

### 形象化理解：过度比对 (Over-alignment)

想象一下我们把一条染色体（Query）切成无数小片段，拿其中一个片段去基因组里“搜寻”：

*   **单拷贝基因**: 这个片段只能在基因组里找到 **1个** 完美匹配的位置（它自己）。
*   **重复元件 (TE)**: 这个片段（比如它是某个 Alu 元件的一部分）能在基因组里找到 **1000个** 匹配位置。
    *   **结论**: 深度为 1000 -> **这是一个重复序列**。
    *   **操作**: 我们把这个片段对应的区域标记为“需要屏蔽”。

### 为什么需要这个流程？

虽然 `RepeatMasker` 等基于库的工具很强大，但它们依赖于已知的重复序列库（RepBase/Dfam）。对于很多非模式物种，我们可能没有完善的重复序列库。

Cactus 的这种**从头（De novo）**检测方法不依赖外部数据库，仅靠序列自身的重复性就能工作，因此对新测序的物种非常有效。

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

整个流程在 `LastzRepeatMaskJob.run()` 方法中定义。为了方便理解，我们可以将其抽象为以下数据流：

```mermaid
graph LR
    A[Query Genome] -->|切片| B(Fragments)
    B -->|Lastz 比对| C(Alignments)
    D[Target Genome] -->|作为参考| C
    C -->|深度统计| E(High Depth Regions)
    E -->|Soft Masking| F[Masked Genome]
```

流程分为三个主要步骤：

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

## 4. 关键参数 (`RepeatMaskOptions`)

在 `cactus_lastzRepeatMask.py` 中，`RepeatMaskOptions` 类定义了控制流程的关键参数：

*   `fragment` (int, default 200): 切片大小。如果为奇数会自动加 1 保证偶数，以便能被 2 整除。
*   `minPeriod` (int, default 10): 最小重复周期/覆盖度阈值。
*   `proportionSampled` (float, default 1.0): 采样比例。实际使用的深度阈值 `period` 计算公式为 `max(1, round(proportionSampled * minPeriod))`。
*   `lastzOpts` (str): 传递给 `lastz` 的额外参数。
*   `unmaskInput` / `unmaskOutput`: 控制输入输出的屏蔽状态。

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

3.  **工具特性 (Lastz Optimization)**:
    *   Lastz 的 `--querydepth=keep,nowarn:N` 参数。
    *   **机制**: 当 Query 序列的某位置比对深度超过阈值 N 时，Lastz 会**强制停止**该 Query 区域的后续搜索。
    *   **性能影响**: 若不设置此参数，运行速度会**大大降低**（Lastz 会试图寻找所有拷贝，导致指数级变慢）。
    *   **结果**: 这意味着对于高拷贝重复序列，我们**只能找到前 N 个 Target 位置**，而遗漏了剩余所有的拷贝位置（即 Target 上的结果是不完整的）。
    *   **结论**: 如果试图从 Target 视角统计深度，会因为大量丢失比对而无法正确识别重复。反之，Query 视角只要达到阈值 N，就已有充分证据将其判定为“重复”并进行屏蔽，这与 Lastz 的截断策略完美契合。

## 5. PGR 工具链替代方案设计

PGR (Phylogenetics in Rust) 项目提供了一套高效的 Rust 工具链，旨在替代上述复杂的 Python/C 混合流程。我们的目标是：**更快、更省内存、更易于部署和使用**。

与 Cactus 的复杂调度不同，PGR 采用模块化的 CLI 工具设计，用户可以通过标准的 Shell 管道或脚本灵活组合。

以下是各个步骤的替代方案详解：

### 5.1 `pgr fa window` vs `cactus_fasta_fragments.py`

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

### 5.2 用法

详细用法请参考 `pgr fa window --help`。

### 5.3 实现细节
*   **流式处理**: 类似于 `pgr fa size`，逐条读取 Record，不需要 `.loc` 索引文件，适合处理巨大文件流。
*   **内存优化**: 仅持有当前 Record 的 Sequence，不加载整个文件。

### 5.4 Target 构建策略 (PGR Design)

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

### 5.5 `pgr lav lastz` vs `lastz` (Wrapper)

为了简化复杂的 Lastz 调用流程，我们实现了 `pgr lav lastz` 命令：

*   **功能**:
    *   **目录/文件递归**: 支持对输入目录进行递归搜索，自动识别 `.fa` 和 `.fa.gz` 文件。
    *   **参数预设**:
        *   自动添加 `[multiple][nameparse=darkspace]` 修饰符。
        *   自动添加 `--querydepth` (默认 13)。
        *   自动设置 `--format=lav` 输出格式（专为 PGR Chaining 流程优化）。
    *   **预设参数集** (`--preset`):
        *   支持 `set01` 到 `set07` (来自 UCSC/Cactus 的经典参数集)。
        *   包含参数配置和打分矩阵 (Matrix)。
        *   使用 `--show-preset` 可查看详细参数。
    *   **自定义参数**:
        *   支持 `--lastz-args "..."` 透传其他高级参数（覆盖预设）。
    *   **自比对优化**:
        *   `--self` 参数：当 Target 和 Query 为同一文件时，自动传递 `--self` 给 Lastz，避免冗余计算。
    *   **并行计算**:
        *   使用 Rust `rayon` 库实现多线程并行执行 (Task Parallelism)，替代了 Perl 版本中的 MCE 模块。

*   **用法示例**:
    ```bash
    # 使用预设参数进行比对 (Human vs Chimp)
    pgr lav lastz target.fa query.fa --preset set01 -o lastz_out
    
    # 目录递归与并行执行
    pgr lav lastz target_dir/ query_dir/ --preset set03 --parallel 8
    
    # 自比对
    pgr lav lastz genome.fa genome.fa --self --preset set01
    ```

### 5.6 实战指南：从头构建 RepeatMasking 流程


为了替代 Cactus 的 `cactus_covered_intervals` 及其复杂流程，我们可以使用 PGR 工具链构建一个清晰的 Shell 脚本。

**场景假设**: 你有一个新组装的基因组 `genome.fa`，想要对其进行重复序列屏蔽（Soft-masking）。

**完整脚本示例**:

```bash
#!/bin/bash
set -e

# 1. 准备工作
INPUT_FA="genome.fa"
WORK_DIR="masking_work"
mkdir -p $WORK_DIR

# 获取染色体大小 (用于后续坐标还原)
pgr fa size $INPUT_FA > $WORK_DIR/genome.sizes

# 2. 切片 (Fragmentation)
# 将基因组切成 200bp 的片段，步长 100bp (50% 重叠)
# 这一步对应 Cactus 的 cactus_fasta_fragments.py
pgr fa window $INPUT_FA --len 200 --step 100 --out $WORK_DIR/fragments.fa

# 3. 比对 (Alignment)
# 将切片后的片段比对回原始基因组 (Self-Alignment)
# 使用 preset set01 (高灵敏度)，输出为 LAV 格式
# 这一步对应 lastz 的执行
pgr lav lastz $INPUT_FA $WORK_DIR/fragments.fa --preset set01 --self -o $WORK_DIR/fragments.lav

# 4. 格式转换与坐标还原
# 4.1 LAV -> PSL
pgr lav to-psl $WORK_DIR/fragments.lav -o $WORK_DIR/fragments.psl

# 4.2 Lift Coordinates
# 将 Fragment 内部坐标 (如 seq1:100-300 的第 10bp) 还原为全基因组坐标 (seq1 的第 110bp)
pgr psl lift $WORK_DIR/fragments.psl --q-sizes $WORK_DIR/genome.sizes -o $WORK_DIR/lifted.psl

# 5. 深度计算与区间提取
# 5.1 PSL -> Range (.rg)
# 提取比对到的区间，自动处理正负链
pgr psl to-range $WORK_DIR/lifted.psl > $WORK_DIR/query_coverage.rg

# 5.2 计算深度 (Depth > 10)
# 任何深度超过 10 的区域被视为重复序列
# 对应 cactus_covered_intervals
spanr coverage $WORK_DIR/query_coverage.rg -m 10 > $WORK_DIR/mask_regions.json

# 6. 应用屏蔽 (Masking)
# 将识别出的区域在原基因组中标记为小写 (Soft-masking)
pgr fa mask $INPUT_FA $WORK_DIR/mask_regions.json -o masked_genome.fa

echo "Masking completed: masked_genome.fa"
```

这个流程清晰地展示了数据是如何流动的：
`Genome` -> `Fragments` -> `Alignments (LAV/PSL)` -> `Coordinates (Lifted)` -> `Depths (Ranges)` -> `Mask (JSON)` -> `Masked Genome`。

