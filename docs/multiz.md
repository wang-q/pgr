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

结合以上分析，可以给出一个比较实用的结论（并补充当前实现现状）：

*   对于最初设计的 `pgr` 主用例（构建严格交集的 core genome，比对区域由 `cover`/`slice` 控制，再由 `fas refine` 精修），**可以不直接复刻 multiz 的完整 Yama 动态规划引擎**。
    *   这一路线本质上假设：在共同核心区域内，各个 pairwise/MAF 中的参考序列 gap 模式差异不大，或者差异可以在后续局部 MSA 中被吸收。
    *   代价是：更偏向“交集/核心”，不会像 multiz 那样追求覆盖所有可能的比对区域（union/mesh）。
*   在此基础上，`pgr` 目前在 `fas` 层引入了一个**简化版的带状 DP 引擎**（见第 8 节 `libs::fas_multiz`）：
    *   在参考坐标网格上做带状 DP，用两个 profile 的物种交集上的 sum-of-pairs 打分（base–base 使用替换矩阵，base–gap 使用统一 gap 罚分），用于解决不同输入之间参考 gap 模式和列选择的冲突。
    *   目前主要针对两个输入 `.fas` 的窗口合并场景，作为 `fas join` 的“智能版”补充，而不是完整的 yama 复刻。
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

### 7.5 pgr 中 multiz 前置链路：LASTZ 与链化

在 UCSC 的典型 WGA 流程中，`multiz` 位于“pairwise 比对 + 链化 + net + mafFromNet”之后，只消费已经整理好的 MAF。`pgr` 目前在这一前置链路上，也已经有相当完整的 Rust 封装，主要对应到：

*   `pgr lav lastz`：LASTZ 前端
    *   位置：`src/cmd_pgr/lav/lastz.rs`。
    *   作用：包装 `lastz`，生成 LAV 格式输出，参数设计对齐 Cactus/UCSC 风格。
    *   特点：
        *   内置 UCSC 风格 preset（`set01`..`set07`），包括常见物种组合（Hg vs Pan/Mm/Bos/DanRer 等），每个 preset 绑定一套参数串和一个 4x4 替换矩阵（通过临时文件写给 `Q=` 选项）。
        *   自动加上 `--format=lav`、`--markend`、`--ambiguous=iupac`、`--querydepth=keep,nowarn:N` 等选项，行为贴近 Cactus repeat-masking 里的 LASTZ 调用约定。
        *   支持单文件和目录递归：对 target/query 目录做递归扫描（`.fa`/`.fa.gz`），生成笛卡尔积 job 列表。
        *   使用 `rayon` 并行跑多个 lastz 进程，并为每个 target–query 组合生成类似 `[t]vs[q].lav` 的输出文件名（带冲突规避逻辑）。
    *   对 multiz 的意义：
        *   对应于“blastz/lastz pairwise 比对”这一步，为后续链化、net 和 multiz/fas-multiz 提供高质量的成对比对基础。

*   `pgr psl chain`：PSL 链化（axtChain 风格）
    *   位置：`src/cmd_pgr/psl/chain.rs`，调用 `libs::chain` 中的 DP 引擎。
    *   作用：把 PSL 对齐 block 链成较长的 syntenic chain，逻辑类似 UCSC 的 `axtChain`/`chainNet` 里的链化步骤。
    *   打分与 gap 模型：
        *   使用 `SubMatrix` 作为替换矩阵，默认 Identity（匹配 +100 / 不匹配 -100），也可通过 `--score-scheme` 选择 HoxD55 或读取 LASTZ 格式打分文件。
        *   gap 成本由 `GapCalc` 驱动：
            *   线性模式：`--linear-gap loose|medium`，对应 Kent 源码中针对远缘/近缘物种的 quasi-natural gap 曲线。
            *   仿射模式：`--gap-open` + `--gap-extend` 显式指定 open/extend，内部通过 `GapCalc::affine` 生成 gap 曲线。
        *   链化 DP 中的评分公式与 UCSC axtChain 一致：`Score = BlockScore + max(PrevScore - GapCost)`。
    *   结构与实现：
        *   依据 `(t_name, q_name, strand)` 分组，内部使用 KD-tree 等结构（见 `libs::chain`）加速前驱 block 搜索。
        *   允许在有 2bit 序列的情况下，用 `ScoreContext` 和 `calc_block_score` 精确重算每个 block 的序列得分，而不是只依赖 PSL 自带分数。
    *   对 multiz 的意义：
        *   在 `pgr` 里，这一步提供了与 UCSC 链化阶段等价的“整理过的 syntenic 对齐骨架”，可以作为（通过 AXT/MAF/FA 转换后）fas-multiz 的上游输入。

综上，`pgr lav lastz` + `pgr psl chain` 组合，大致覆盖了 UCSC 链路中 “blastz/lastz 比对 + axtChain 链化” 这两步。它们提供了 multiz/fas-multiz 所需的 pairwise 对齐基础，而 `libs::fas_multiz` 则承担了更上游的 profile 合并角色：在已经有 syntenic 对齐骨架的前提下，对多个 `.fas` profile 做带状 DP 合并，构建 union/mesh 风格的多序列比对。

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

### 8.5 命令行接口草案

*   子命令名称（示例）：
    *   `pgr fas multiz`
*   核心参数（示例）：
*   `-r, --ref <NAME>`：参考物种/序列名称，必须在所有输入 `.fas` 中存在。
*   `<infiles>...`：位置参数，输入的 block FA 文件，数量 `>= 2`，行为与 `pl p2m` 一致。
    *   `--radius <INT>`：带状 DP 半径 `R`，类似 multiz 中的 `R`，控制参考坐标附近的搜索宽度。
    *   `--min-width <INT>`：最小输出 block 宽度，对标 multiz 的 `M`。
    *   `-o, --out <FILE>`：输出 `.fas` 文件名。
    *   `--score-matrix <FILE>`：可选，指定替换矩阵（默认可复用 `libs/chain/sub_matrix.rs` 中已有配置）。
    *   `--mode <core|union>`：模式切换：
        *   `core`：在交集区域内行为尽量贴近 `p2m + join`，只对 gap 冲突做最小修复。
        *   `union`：尽量保留所有输入的对齐信息，生成 mesh 风格结果。

### 8.6 约束与实现注意事项

*   **参考骨架一致性**：
    *   要求所有输入 `.fas` 的参考序列来自同一基因组版本，且建议事先经过相同的 masking/裁剪流程。
*   **窗口化处理**：
    *   实现时应采用窗口化策略（例如按固定长度或按 block 切分），避免在超长区间上运行大规模 profile DP。
    *   每个窗口内的 profile 合并结果可以再交给 `fas refine` 做一次本地 MSA。
*   **打分与带状 DP**：
    *   可以重用 `pgr` 中现有的打分矩阵和 gap 参数（如 `libs/chain` 相关代码），避免在 `fas` 层重新定义一套 scoring。
    *   带状 DP 的半径 `R` 和最小宽度 `M` 建议和 multiz 保持同一数量级，以便结果直观可控。
*   **失败与降级策略**：
    *   当某个窗口内 profile DP 无法找到合理路径（打分过低或冲突过多）时，可以退回到简单的 `fas join` 行为，或干脆将该窗口标记为“未合并”，交给上游/下游流程决定如何处理。

### 8.7 与现有模块的集成点

*   **输入准备**：依赖现有的 `pgr axt/maf to-fas` 和 `fas` 系列命令，将所有上游结果规整为块级 `.fas`。
*   **区间计算**：复用 `fas cover` 和 `spanr` 的区间逻辑，定义候选合并窗口。
*   **比对与 refine**：在新实现的 fas-multiz 中完成 profile 合并后，调用现有 `pgr fas refine` 作为可选的精修步骤。
*   **下游分析**：输出 `.fas` 可以继续被 `fas stat`, `fas to-vcf`, `fas split` 等命令消费，与当前 `p2m + join + refine` 的结果处于同一生态。

### 8.8 libs 实现草案

> 2026-02 更新：本节给出的设计已经在 `libs::fas_multiz` 中基本落地实现（包括 `FasMultizMode`/`FasMultizConfig`/`Window` 以及 `merge_window`、`merge_fas_files`、自动窗口推导等），但仍保留“草案”形式以便对比 multiz 原始设想与当前实现。当前实现细节与局限见 8.10 小节。

为了方便后续在 Rust 中实现 fas-multiz，这里给出一个初步的 libs 级别设计。

*   **模块位置**：
    *   新增 `src/libs/fas_multiz.rs`，在 `src/libs/mod.rs` 中通过 `pub mod fas_multiz;` 暴露。
*   **依赖复用**：
    *   解析 `.fas`：复用 `libs::fmt::fas` 中的 `FasEntry`、`FasBlock`、`next_fas_block` 等。
    *   区间坐标：继续使用 `intspan::Range`。
    *   打分与碱基类型：复用 `libs::nt::NT_VAL` 以及 `libs::chain` 中已有的替换矩阵和 gap 参数。
    *   简单统计/评估：如有需要，可调用 `libs::alignment::alignment_stat` 做 sanity check。
*   **核心类型**：
    *   合并模式：
        ```rust
        pub enum FasMultizMode {
            Core,
            Union,
        }
        ```
    *   配置结构：
        ```rust
        pub struct FasMultizConfig {
            pub ref_name: String,
            pub radius: usize,
            pub min_width: usize,
            pub mode: FasMultizMode,
        }
        ```
    *   窗口定义：
        ```rust
        pub struct Window {
            pub chr: String,
            pub start: u64,
            pub end: u64,
        }
        ```
*   **对外 API 草图**：
    *   文件级合并（供 CLI 使用）：
        ```rust
        pub fn merge_fas_files(
            ref_name: &str,
            infiles: &[impl AsRef<Path>],
            windows: &[Window],
            cfg: &FasMultizConfig,
        ) -> anyhow::Result<Vec<FasBlock>>;
        ```
        *   读入多个 `.fas` 文件，根据给定窗口把 block 分组，对每个窗口调用 `merge_window`，最终返回按参考坐标排序的一组 `FasBlock`。
    *   单窗口合并（算法核心）：
        ```rust
        pub fn merge_window(
            ref_name: &str,
            window: &Window,
            blocks_per_input: &[Vec<FasBlock>],
            cfg: &FasMultizConfig,
        ) -> Option<FasBlock>;
        ```
        *   输入是某个窗口内来自多个文件的 block 集合，输出是一个合并后的 block（或在无法合理合并时返回 `None`）。
*   **merge_window 内部步骤概述**：
*   将每个输入中参考物种的 `FasEntry` 映射到统一的参考坐标网格上，得到多条略有差异的参考轨迹。
*   在参考轨迹之间执行带状 profile 对齐：DP 网格仍然只在参考坐标上展开，但每个对角单元的得分由两个 profile 的物种交集上的 sum-of-pairs 决定（共享物种的 base–base 使用 `libs::chain::SubMatrix::hoxd55` 的替换分数并做适当缩放，base–gap 使用统一的 gap 罚分，gap–gap 不计分），对两个及以上输入采用 progressive 带状 DP。
*   按照合并后的参考轨迹，对每个输入的非参考序列进行重采样：在缺失列处插入 gap，在 Union 模式下允许在参考 gap 位置引入新列，在 Core 模式下则尽量丢弃不一致列；非参考物种在 DP 打分中参与 sum-of-pairs，但坐标仍然沿参考轨迹重采样。
*   将重采样后的各物种序列按列拼接，构造新的 `FasBlock`，并为参考 entry 生成合适的 `Range`（可以取窗口的 Range 或交集 Range）。
*   如果在某个窗口内 profile 对齐得分过低或冲突过多，则返回 `None`，由调用者决定使用简单 `fas join` 还是跳过该窗口。

### 8.9 与 multiz-multiz 源码的异同

这里的 fas-multiz 方案是从 `multiz-multiz` 的源码和算法抽象出来的一个 “pgr 版本”，既保留了一些核心思想，也刻意做了简化和调整。

*   **共同点（继承 multiz 的部分）**：
    *   都是以参考物种坐标为主轴，在参考坐标上定义窗口/段落，再在每个窗口内做 profile 合并。
    *   都采用带状 DP（或类似思想）限制搜索空间，在参考附近做局部优化，而不是在全空间做 MSA。
    *   在 union/mesh 场景下，都试图尽可能保留不同输入中的真实比对关系，只在必要时删除或压缩冲突列。
    *   支持“核心交集 + 扩展区域”的思路：核心部分倾向于各输入一致，边缘部分允许有差异并通过 DP 协调。
*   **差异（pgr 有意做的调整）**：
    *   **工作层级不同**：
        *   multiz-multiz 在 MAF 层操作，直接对齐两个 MAF profile。
        *   fas-multiz 在 block FA 层操作，输入是多个 `.fas` 文件，链路由 `pgr axt/maf to-fas` 标准化过，因此语法和元数据更简单。
    *   **DP 引擎复杂度不同**：
        *   multiz-multiz 的 `yama` 部分实现了一套完整的 profile–profile DP 引擎，并考虑了较丰富的状态与回溯（sum-of-pairs + 替换矩阵 + 仿射 gap 模型）。
        *   fas-multiz 只在参考轨迹上做简化的带状对齐：采用 progressive pairwise DP 和可配置的 match/mismatch/gap 标量打分，重点解决参考 gap 冲突和列选择问题，不追求完全复刻 yama 的所有细节状态。
    *   **实现位置与职责边界不同**：
        *   multiz-multiz 是一个专门服务于 MAF/多序列比对构建的独立 C 项目。
        *   fas-multiz 被设计为 `pgr` 的一个 libs 模块（`libs::fas_multiz`），与现有 `fas cover/slice/join/refine` 等命令协作，而不是独立的 pipeline。
    *   **输入准备和预处理链路不同**：
        *   multiz-multiz 直接消费上游链路输出的 MAF（如 blastz/last 等）。
        *   在 pgr 中，上游的 pairwise 结果通常已经通过若干步骤转换、规范成 `.fas`，fas-multiz 可以假设这些输入已经做过一次清洗/规整。
    *   **目标偏好与使用场景不同**：
        *   multiz-multiz 更偏“通用 WGA 引擎”，追求在大范围基因组上做 mesh 式对齐。
        *   fas-multiz 明确被设计成 pgr 的一个“union/mesh complement”：在 core/intersection 流程之外，提供一个额外的 union 视角，并保持与现有 `p2m + join + refine` 在交集区域内尽量兼容。

### 8.10 当前 fas-multiz 实现状态（2026-02）

> 本节描述的是当前 `pgr` 仓库中已经落地的 `libs::fas_multiz` 实现，用于对照前文的 multiz 设想。实现仍然是“轻量级 fas-multiz”，未来可以继续向更完整的 profile–profile DP 演进。

**实现位置与对外 API**

*   模块位置：`src/libs/fas_multiz.rs`，通过 `pub mod fas_multiz;` 暴露为 `pgr::libs::fas_multiz`。
*   核心类型：与 8.8 草案基本一致，并在配置中加入了 DP 打分参数：
    *   `FasMultizMode { Core, Union }`
    *   `FasMultizConfig { ref_name, radius, min_width, mode, match_score, mismatch_score, gap_score }`
    *   `Window { chr, start, end }`
*   对外函数：
    *   `merge_window(ref_name, window, blocks_per_input, cfg) -> Option<FasBlock>`
    *   `merge_fas_files(ref_name, infiles, windows, cfg) -> Result<Vec<FasBlock>>`
    *   `merge_fas_files_auto_windows(ref_name, infiles, cfg) -> Result<Vec<FasBlock>>`

**窗口推导与 Core/Union 语义**

*   `merge_fas_files` 需要调用方显式给出 `windows`，行为与草案一致。
*   `merge_fas_files_auto_windows` 会：
    *   从所有输入 `.fas` 中提取参考物种 `ref_name` 的 `Range`，按 `radius` 向两侧扩展。
    *   按染色体合并重叠区间，再按 `min_width` 过滤过短窗口。
    *   按 `cfg.mode` 过滤窗口：
        *   `Core`：只保留“在所有输入中都有参考覆盖”的窗口（严格交集）。
        *   `Union`：只要有任意一个输入在该窗口有参考覆盖即可保留（并集风格）。

**窗口内合并逻辑（带状 DP 合并）**

*   一般情况（任意输入个数）：
*   对于给定窗口，先从每个输入文件中选出在窗口内与参考重叠的 block，组成 `blocks`。
*   若 `blocks` 为空，或（在 Core 模式下）某些输入找不到参考 block，则直接返回 `None`。
*   若 `blocks.len() >= 2`，先尝试 progressive 带状 DP 合并：
*   使用内部函数 `merge_blocks_with_dp`，按顺序对 `blocks` 做两两 DP 合并。
*   每一步都要求参与合并的参考 entry 在去掉 `'-'` 后的序列完全相同（ungapped equal），否则这一轮 DP 失败。
*   在参考坐标网格上调用 `banded_align_refs`：
*   只在 diagonal ± `radius` 的带内做 DP。
*   对每个对角单元，使用两个 profile 的物种交集上的 sum-of-pairs 打分：共享物种的 base–base 组合通过 `libs::chain::SubMatrix::hoxd55` 获取替换分数并按固定比例缩放到与 `match_score` 相近的量级；base–gap 组合收取由 `gap_model` 推导出的统一 gap 罚分（`constant` 模式直接使用 `gap_score`，`medium`/`loose` 模式则从 `GapCalc` 的 quasi-natural 曲线抽取一个与替换分数量级相当的 gap 罚分），gap–gap 不计分；横向/纵向移动同样只收取一次该 gap 罚分。
*   将 DP 生成的参考轨迹映射到所有物种：
*   对每一列，优先从前一个累积结果（或第一个输入）的对应位置取碱基，不存在时再从当前输入取；两边都缺失则填 `'-'`。
*   `Core` 模式下只合并在当前累积结果和新输入中都存在的物种；`Union` 模式下允许物种只存在于其中一边。
*   在 Core 模式下，任一步 DP 失败都会导致整个 progressive 合并失败，随后回退到“保守合并”逻辑。
*   在 Union 模式下，如果某一步 DP 失败，则跳过该输入，继续尝试将后续输入与当前累积结果进行 DP 合并；成功的部分会被保留，无法对齐的输入则在该窗口中被忽略。
*   progressive DP 完成后（无论是否跳过了一些输入），若至少完成了一次成功的 DP 合并，则直接返回最终累积的 block。
*   如果 progressive DP 入口阶段就失败（例如前两条参考轨迹 ungapped 不同，或带宽内找不到合理路径），则自动回退到“保守合并”逻辑：
*   要求所有候选 block 的参考 entry 完全相同（包含 gap），否则返回 `None`。
*   `Core` 模式下只保留在所有输入中都存在的物种；`Union` 模式下保留物种并集。
*   参考物种的 `Range` 继承自模板 block；其他物种继承其来自的原始 block。

**当前实现的局限与后续扩展方向**

*   目前 DP 采用 progressive pairwise 策略，对多个输入的合并顺序敏感，尚未实现真正意义上的多维 profile–profile sum-of-pairs 动态规划。
*   DP 网格仍然只在参考行的坐标上展开，非参考物种没有各自独立的坐标轴，它们通过物种交集上的 profile–profile sum-of-pairs 打分参与评分，但不改变 DP 网格结构，与 multiz/yama 中更完整的多维 DP 仍有差距。
*   替换分数已经复用 `libs::chain::SubMatrix` 做 base–base 的 sum-of-pairs 打分：默认使用 `hoxd55`，也支持通过 `--score-matrix` 读取 LASTZ 格式文件或预设名称（例如 `hoxd55`），并通过简单缩放与当前 `match_score` 的量级对齐。gap 支持三类模型：`constant`、`medium`/`loose`、以及显式仿射：
    *   `constant`：直接使用 `gap_score` 作为统一线性 gap 罚分。
    *   `medium`/`loose`：从 `GapCalc::medium`/`GapCalc::loose` 的 quasi-natural 曲线中取 `len=1,2` 两点，反推出一组近似的仿射参数 `(open, extend)`，再按 HoxD55 的打分尺度和 `match_score` 做线性缩放，在带状 DP 中用“open + extend × length”的形式累积 gap 罚分，从而实现长度依赖的 quasi-natural 近似。
    *   显式仿射：当通过 `--gap-open`/`--gap-extend` 提供 open/extend 时，fas-multiz 在 DP 中直接使用这一组仿射参数（同样按 `match_score` 缩放）进行三状态的仿射 gap 计分。
*   已提供 CLI 子命令 `pgr fas multiz`（见 8.5 小节），支持 `--mode core|union`、`--radius`、`--min-width`、`--gap-model`、`--gap-open`、`--gap-extend` 以及 `--score-matrix` 等参数；gap 配置风格与 `pgr psl chain` 保持一致，而替换矩阵也不再局限于内置的 HoxD55，可与链化阶段共享同一套 matrix 配置；`libs::fas_multiz` 仍作为底层引擎，便于在 pipeline 或其他子命令中复用。
*   在 gap 行为上，对端部 gap（leading/trailing gap）增加了简单的“首尾特化”规则：在带状 DP 回溯得到参考物种之间的对齐路径后，会自动裁剪掉首尾连续的单侧 gap 列（即仅一侧为碱基、另一侧为 gap 的前缀/后缀列），使这些端部 gap 在行为上视为 free end gaps，而中间区域仍按标准仿射 gap 计分；若需要更复杂的端部 gap 放宽或偏置策略，可以在这一基础上继续细化。
*   在上述基础上，仍可以在后续逐步接近 multiz 的完整行为，例如：
    *   将 progressive DP 升级为真正的多输入 profile–profile sum-of-pairs 动态规划。
*   在 DP 失败时更智能地选择降级策略（退回 `fas join`、标记窗口未合并等）。
