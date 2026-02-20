# Multiz 分析笔记

本文档分析位于 `multiz/` 目录下的源码，该项目是 UCSC `multiz-tba` 软件包的核心组件之一。

## 1. 项目概览

*   **工具名称**: `multiz`
*   **版本**: 11.2
*   **功能**: 对两个多序列比对文件（MAF 格式）进行比对合并。
*   **核心假设**: 两个输入 MAF 文件的第一行（Top Row）必须是相同的参考序列（Reference Sequence）。`multiz` 利用这个共同的参考序列作为锚点，将其他的序列比对合并进来。

## 2. 核心算法

基于 `mz_yama.c` 和 `multiz.c` 的源码分析。

### 2.1 动态规划 (Dynamic Programming)
`multiz` 使用动态规划算法来比对两个 profile（或两个 alignment block 集合）。

*   **打分机制**: Sum-of-pairs (SP) substitution scores。
*   **Gap Penalty**: 采用 "Quasi-natural gap costs" (Altschul 1989) 和仿射 Gap 罚分 (Affine Gap Scores)。
    *   End-gaps (首尾 Gap) 通常不收取 Gap open penalty。
*   **实现**: 基于 Myers & Miller (1989) 的网格模型。

### 2.2 带状优化 (Banded DP)
为了提高效率，`multiz` 并不总是计算完整的 DP 矩阵，而是采用带状（Banded）动态规划。
*   **Radius (R)**: `multiz` 接受一个参数 `R` (默认为 30)，定义了 DP 搜索的半径。
*   **LB/RB**: 代码中使用 `LB` (Left Boundary) 和 `RB` (Right Boundary) 数组来限制每一行 DP 的计算范围，只在参考序列对齐位置附近的 `R` 范围内进行搜索。

## 3. 源码结构

| 文件 | 描述 |
| :--- | :--- |
| `multiz.c` | **主程序入口**。负责处理命令行参数，读取输入 MAF，调度比对逻辑。 |
| `mz_yama.c` | **核心算法**。实现了 `yama` 函数，即上述的动态规划比对引擎。 |
| `mz_preyama.c` | 预处理模块，负责确定比对的边界和重叠区域，为 `yama` 准备数据。 |
| `maf.c` / `maf.h` | MAF 格式的读写库。定义了 `mafAli` 结构体。 |
| `mz_scores.c` | 定义了打分矩阵和 Gap 罚分参数。 |
| `util.c` | 通用工具函数（内存分配、错误处理等）。 |

## 4. 使用方法

```bash
multiz file1.maf file2.maf v [out1 out2]
```

*   `file1.maf`, `file2.maf`: 两个输入的 MAF 文件，第一行必须是相同的参考序列。
*   `v`: 版本/模式控制。
    *   `0`: 参考序列在两个文件中都不固定（允许微调）。
    *   `1`: 参考序列在第一个文件中固定。
*   `out1`, `out2` (可选): 分别用于收集 `file1.maf` 和 `file2.maf` 中未能比对上的（Unused）Block。
    *   `out1`: 收集 `file1.maf` 中未被使用的 Block。
    *   `out2`: 收集 `file2.maf` 中未被使用的 Block。
*   **参数控制**:
    *   `R=30`: 调整 DP 半径。
    *   `M=1`: 最小输出宽度。

## 5. 对 `pgr` 项目的启示

1.  **Profile Alignment**: `multiz` 本质上是一个 Profile-Profile Aligner。在 `pgr` 的 `fas consensus` 或图比对中，如果涉及多序列合并，可以参考其 Sum-of-pairs 的处理方式。
2.  **锚点策略**: 利用共同参考序列作为锚点是处理大规模比对的有效手段（这也是 `pgr fas join` 的逻辑基础）。
3.  **MAF 处理**: `multiz` 的 `maf.c` 提供了一个轻量级的 MAF 解析参考实现。**更新**: `pgr` 现已实现完整的 MAF 读写功能（见第 6 节），采用了更符合 Rust 特性的设计，但字段定义上仍与 MAF 规范保持一致。

## 6. pgr MAF 实现现状 (2025-02)

目前 `pgr` 已在 `src/libs/fmt/` 下实现了完整的 MAF 读写支持。

### 6.1 模块分布与功能

*   **读取 (Reader)**: 位于 `src/libs/fmt/fas.rs`
    *   **核心函数**: `next_maf_block`, `parse_maf_block`。
    *   **特性**: 
        *   支持 `a` (alignment) 和 `s` (sequence) 行解析。
        *   **坐标转换**: 内置 `to_range()` 方法，将 MAF 的 0-based 坐标转换为 1-based inclusive 格式（如 `chr:start-end`）。
        *   **负链处理**: 自动处理负链坐标，将其转换为相对于正链的坐标范围。
*   **写入 (Writer)**: 位于 `src/libs/fmt/maf.rs`
    *   **核心结构**: `MafWriter`。
    *   **特性**: 支持输出标准 MAF 头信息 (`##maf`) 和对齐块，自动处理列宽对齐。

### 6.2 数据结构对比
 
 目前读写使用两套略有不同的结构，开发时需注意转换：
 
 | 特性 | 读取 (fas.rs) | 写入 (maf.rs) |
 | :--- | :--- | :--- |
 | **结构体** | `MafEntry` | `MafComp` |
 | **序列存储** | `Vec<u8>` | `String` |
 | **数值类型** | `u64` | `usize` |
 
 ### 6.3 MAF 格式实现对比 (pgr vs multiz vs UCSC)
 
 通过对比 `multiz/maf.c` (mini-maf) 与 `chainnet/src/lib/maf.c` (UCSC完整版)，以及 `pgr` 的实现，可以发现：
 
 1.  **代码同源性**:
     *   `multiz` 中的 `maf.c` (header 注明 "version 12") 明确标注了 "Stolen from Jim Kent & seriously abused"。它是一个精简版 (mini-maf)，移除了大量依赖 (如 `linefile.h`, `common.h` 等)，直接使用标准 C 库函数。
     *   UCSC `chainnet` 中的 `maf.c` 是完整版，依赖于 Kent Source 庞大的基础设施库 (`common.h`, `linefile.h`, `hash.h` 等)。
 
 2.  **功能差异**:
     *   **UCSC 完整版**:
         *   支持 `i` (synteny breaks), `q` (quality), `e` (empty/bridging), `r` (region definition) 等多种扩展行。
         *   拥有复杂的内存管理 (`AllocVar`, `slAddHead` 等 Kent 库特有宏)。
         *   包含大量辅助函数，如 `mafSubset` (切片), `mafFlipStrand` (反向互补), `mafScoreMultiz` (打分) 等。
     *   **multiz 精简版**:
         *   仅保留核心的 `a` (alignment) 和 `s` (sequence) 行解析，足以支持比对算法。
         *   内存管理简化为 `ckalloc`。
         *   去除了大部分与比对算法无关的辅助功能。
     *   **pgr 实现**:
         *   **解析策略**: 类似于 `multiz`，专注于核心的 `a` 和 `s` 行，忽略非标准行（但解析器通常具有鲁棒性，能跳过未知行）。
         *   **坐标系统**: 与 UCSC/multiz 保持一致 (0-based start, 1-based size)，但内部提供了向 `1-based inclusive` (GFF/VCF 风格) 转换的接口，适应现代分析需求。
         *   **内存模型**: 使用 Rust 的所有权模型 (`String`, `Vec`) 替代 C 指针，杜绝了内存泄漏风险，且无需手动管理 `free`。
 
 ### 6.4 总结
 
 `pgr` 的 MAF 模块已经成熟，具备处理 UCSC MAF 格式的能力，并集成了坐标标准化功能，适合作为后续基因组比对分析的基础组件。

## 7. pgr 与 multiz 的异同分析

`pgr` 采用 **"Stitch + Refine" (拼接+精炼)** 的分步策略，与 `multiz` 的 **"Integrated Alignment" (一体化比对)** 形成鲜明对比。

### 7.1 pgr 的分步工作流

#### 第一步：机械拼接 (Stitch)
*   **工具**: `pgr fas join` (或流程脚本 `pgr pl p2m`)
*   **实现细节**:
    1.  **锚定**: 用户指定一个参考物种（Reference Target）。
    2.  **交集计算**: `pgr pl p2m` 流程会调用 `fas cover` 和 `spanr intersect` 计算所有物种共有的基因组区域。
    3.  **切片**: `fas slice` 根据交集范围从原始 Pairwise 文件中提取序列。
    4.  **堆叠**: `fas join` 以参考序列的坐标范围为 Key，将来自不同文件的 Block 机械地堆叠在一起。
*   **特点**:
    *   **极速**: 仅涉及 I/O 和坐标计算。
    *   **局限**: 不处理 Gap。它假设所有输入文件在同一坐标下的参考序列片段是完全一致的。若不同 Pairwise Alignment 中参考序列的 Gap 状态不一致，直接 Join 会导致非参考序列错位。它没有 `multiz` 的 Yama 动态规划引擎来解决冲突。

#### 第二步：重新比对 (Refine)
*   **工具**: `pgr fas refine`
*   **逻辑**: 对拼接后的 Block 进行多序列比对 (MSA)。
*   **作用**: 弥补第一步的缺陷。调用 `mafft`, `muscle` 或内置 POA 引擎，修正由于机械堆叠导致的 Gap 错位，生成最终的高质量比对。

### 7.2 策略深度对比

#### multiz (UCSC)
*   **模式**: **一体化动态规划 (Integrated DP)**。
*   **核心算法**: Yama (Sum-of-pairs + Gap Costs)。
*   **处理逻辑**: 在合并过程中实时解决 Gap 冲突，插入新的 Gap 以对齐所有序列。
*   **输出目标**: **Union/Mesh** (保留所有可能的比对区域)。

#### pgr (p2m + join + refine)
*   **模式**: **拼接后修正 (Post-hoc Refinement)**。
*   **核心算法**: Set Ops (Intersection) -> Stacking -> MSA (Refinement)。
*   **处理逻辑**: 先根据坐标“硬”合并，再通过 MSA 工具“软”调整。
*   **输出目标**: **Core/Intersection** (通常仅关注严格的交集核心区)。

### 7.3 结论
*   `multiz` 适合构建复杂的、包含大量 Indel 和重排的全基因组比对 (WGA)。
*   `pgr` 流程适合快速构建**核心基因组 (Core Genome)** 或处理**基于无 Gap 参考骨架**的数据。通过 `refine` 步骤，`pgr` 也能产出高质量的比对，但其依赖于“共同核心”的存在。

### 7.4 在 pgr 中是否需要 Yama DP

结合以上分析，可以给出一个比较实用的结论：

*   对于当前 `pgr` 的主用例（构建严格交集的 core genome，比对区域由 `cover`/`slice` 控制，再由 `fas refine` 精修），**可以不实现 multiz 的 Yama 动态规划引擎**。
    *   这一路线本质上假设：在共同核心区域内，各个 pairwise/MAF 中的参考序列 gap 模式差异不大，或者差异可以在后续局部 MSA 中被吸收。
    *   代价是：更偏向“交集/核心”，不会像 multiz 那样追求覆盖所有可能的比对区域（union/mesh）。
*   Yama DP 主要解决的是两类问题：
    1.  **参考序列 gap 冲突**：不同 MAF 中参考序列的 gap pattern 不一致时，通过 DP 在合并过程中实时插入/调整 gap，使所有序列在同一参考坐标系下保持一致。
    2.  **block 级冲突与 unused block 判定**：通过 sum-of-pairs 全局评分决定哪些 block 被合并，哪些落到 `out1`/`out2`。
    在 `pgr` 当前的“Set Ops + Refine” 模式下，这两类问题分别由“交集窗口的选取”和“后续 MSA 的局部重比对”粗略处理，并不追求与 multiz 完全一致的行为。
*   因此，**在如下前提下，舍弃 DP 是可接受的工程权衡**：
    *   只关注严格交集区域，对边缘和稀有对齐不敏感。
    *   参考骨架事先经过统一处理（同一版本、类似 masking/裁剪策略），大型重排通过 `chain/net` 等流程已先行解决。
    *   `fas refine` 作用在规模适中的窗口上，用于修正机械堆叠引入的局部错位，而不负责重新定义 block 边界。
*   若未来需要在 `pgr` 中支持接近 multiz 行为的 “Union/Mesh WGA” 模式，更自然的路径是：**在 `fas` 层实现 multiz 类功能**，而不是在 MAF 解析层再实现一次 `multiz`：
    *   上游已经可以通过 `pgr axt/maf to-fas` 等命令，将 pairwise 或 MAF 转换为 block FA (`.fas`)。
    *   在 `.fas` 层进行 profile 合并，可以直接对齐到 `pgr fas` 现有生态（`cover`/`slice`/`join`/`refine`/`stat` 等），避免重复造 MAF 级别的轮子。

## 8. 在 fas 层实现 multiz 类功能的设想

在 `pgr` 现有设计中，`fas` 是“块级多序列比对”的核心抽象：每个 block 表示一段参考坐标下的多物种比对。若要实现 multiz 类功能，更自然的选择是围绕 `fas` 做 profile 合并，而不是在 MAF 文本层面复刻 `multiz`。

### 8.1 目标与输入/输出

*   **目标**：
    *   给定多个 block FA 文件（例如多个 pairwise-derived `.fas` 或不同 pipeline 生成的 `.fas`），在共享参考物种的坐标系下，将它们合并为一个“union/mesh 风格”的 block FA。
    *   和现有 `p2m + join + refine` 所产出的 **Core/Intersection** 结果互为补充：一个偏交集（core），一个偏并集（union/mesh）。
*   **输入**：
    *   `k` 个 `.fas` 文件（`k >= 2`），它们的 block 中均包含同名的参考序列（例如 `ref`）。
    *   可选的：一个核心交集区域（来自 `fas cover` + `spanr`），用于限制计算范围。
*   **输出**：
    *   一个新的 `.fas` 文件，包含合并后的多序列比对 block：
        *   在交集区域内，行为应与当前 `p2m + join + refine` 相兼容。
        *   在边缘/非完全交集区域内，会尽量保留来自不同输入的对齐（union 行为）。

### 8.2 相对于 multiz 的主要差异

*   **工作层级不同**：
    *   multiz 直接操作 MAF，对齐的是两个 MAF profile。
    *   `pgr` 中的 fas-multiz 将直接操作 block FA，对齐的是若干 `.fas` profile。
*   **上游数据准备不同**：
    *   在 `pgr` 中，pairwise MAF/AXT 等通常已经通过 `pgr maf/axt to-fas` 等步骤规整为更统一的 block FA 表达。
    *   这意味着 fas-multiz 可以假设输入已经经过一次“标准化”，不需要自己处理复杂的 MAF 语法和扩展行。
*   **与现有命令的关系**：
    *   fas-multiz 更像是 `pgr fas join` 的“智能版/DP 版”：
        *   `fas join`：根据参考坐标“机械堆叠”，不处理 gap 冲突。
        *   fas-multiz：在堆叠时引入 profile–profile 的 DP/启发式，解决参考 gap 冲突和 block 选择问题。
    *   输出仍然是 `.fas`，可以直接接上 `fas refine`, `fas stat`, `fas to-vcf` 等命令。

### 8.3 可能的数据流设计（草案）

1.  **标准化输入**：
    *   所有 upstream 比对结果（MAF/AXT 等）先通过现有命令统一转为 `.fas`。
    *   如有需要，可加一步 `fas normalize`（对序列名、物种名、参考 ID 做统一）。
2.  **block 级别的配对与聚类**：
    *   按参考物种与坐标对 block 做分组，将“位置相近”的 block 视为候选合并单元。
    *   这一层可以重用 `fas cover` / `spanr` 得到的区间信息。
3.  **profile 合并（multiz-like）**：
    *   对每个候选区间内的多个 block profile，执行简化版的 profile–profile DP 或其他启发式：
        *   在参考坐标附近采用带状 DP（Radius R），解决不同 `.fas` 之间参考 gap 的不一致。
        *   根据 sum-of-pairs 打分决定保留哪些列/序列，以及如何插入额外 gap。
    *   输出合并后的单个 block（或少数几个 block）。
4.  **后处理与 refine**：
    *   输出的 `.fas` 可以再交给 `pgr fas refine` 做局部 MSA，以获得更“平滑”的 alignment（尤其是在非参考序列上）。

### 8.4 与现有 core 流程的互补关系

*   `p2m + join + refine`：
    *   假设参考骨架在各数据源中基本一致。
    *   倾向于“只相信大家都同意的部分”（严格交集），适合构建 core genome。
*   fas-multiz：
    *   允许不同数据源在边缘和 gap pattern 上存在一定差异，通过 profile 合并策略尽量“合在一起”。
    *   输出更偏 union/mesh，适合探索 union pan-genome 或 WGA 风格的结果。

在实现层面，fas-multiz 可以作为一个新的子命令（例如 `pgr fas multiz` 或 `pgr fas merge-mesh`），并明确声明它与 `p2m + join + refine` 的适用场景不同：前者追求覆盖度（union），后者继续服务于一致性（intersection）。
