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

### 6.3 总结
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
