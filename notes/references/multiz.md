# multiz 源码分析

> 整理于 2026-02，源自对 `multiz/` 目录源码的通读。目的：理解 multiz 的 profile–profile DP 算法，
> 为 pgr fas-multiz 设计提供参考。设计部分见 [[fas-multiz.md]]。

本文档分析位于 `multiz/` 目录下的源码，该项目是 UCSC `multiz-tba` 软件包的核心组件之一。

## 1. multiz 概览

### 1.1 工具与输入输出

*   **工具名称**: `multiz`（本文基于 11.2 版本）。
*   **功能**: 对两个多序列比对文件（MAF 格式）进行比对合并，生成新的多序列比对。
*   **核心假设**: 两个输入 MAF 文件的第一行（Top Row）必须是相同的参考序列（Reference Sequence）。`multiz` 利用这个共同的参考序列作为锚点，将其他物种的比对统一到同一坐标系中。

**命令行接口（简化版）**：

```bash
multiz file1.maf file2.maf v [out1 out2]
```

*   `file1.maf`, `file2.maf`: 两个输入的 MAF 文件，Top Row 为同一参考。
*   `v`: 版本/模式控制。
    *   `0`: 参考在两个文件中都可微调。
    *   `1`: 参考在第一个文件中固定，第二个文件的参考可相对其滑动。
*   `out1`, `out2`（可选）:
    *   `out1`: 收集 `file1.maf` 中未被使用的 block。
    *   `out2`: 收集 `file2.maf` 中未被使用的 block。
    *   若未提供这两个文件，未被使用的 block 将不会额外输出。

常用参数（与 DP 行为直接相关）：

*   `R=30`: 带状 DP 的半径（band radius），控制在参考坐标附近的搜索宽度。
*   `M=1`: 最小输出宽度，小于该宽度的结果不会输出。

### 1.2 核心算法与 Gap 模型

基于 `mz_yama.c` 和 `multiz.c` 的源码，可以把 multiz 的核心算法概括为"参考锚定下的 profile–profile 动态规划"：

*   **Sum-of-pairs 打分**:
    *   对 profile–profile 的每一列，按所有物种对的替换矩阵分数求和（sum-of-pairs）。
*   **Gap 成本**:
    *   使用 quasi-natural gap costs（Altschul 1989）和仿射 gap（affine gap）相结合的模型。
    *   对内部 gap（序列中间）收取 open + extend × length 的成本。
    *   对 end-gaps（首尾 gap）通常不收取 gap open 罚分，更偏向 free-end gap 的行为。
*   **DP 实现**:
    *   采用类似 Myers & Miller (1989) 网格模型的实现方式，在参考坐标上做带状 DP。
    *   Yama 引擎在 DP 过程中实时解决参考 gap 冲突、插入必要的列，以维持所有物种的一致对齐。

### 1.3 带状优化与边界控制

为了避免在全矩阵上做 O(L²) 的 DP，multiz 引入了带状优化：

*   **半径 R**: 通过 `R` 参数限定允许的偏移范围，只在参考两条轨迹的"对齐对角线"附近计算 DP。
*   **LB/RB 边界数组**:
    *   每一行 DP 上维护 `LB`（Left Boundary）和 `RB`（Right Boundary），限制当前行能访问的列范围。
    *   在参考上本来就不可能匹配到的位置不会进入 DP，大幅减少计算量。

从 pgr 的角度看，这套机制和 `libs::fas_multiz` 里的"radius + 带状 DP"是一致的，只是 multiz 在 MAF/profile 层实现，而 fas-multiz 在 `.fas`/block 层实现。

### 1.4 源码结构（bird's eye）

| 文件 | 描述 |
| :--- | :--- |
| `multiz.c` | 主程序入口：命令行解析、输入 MAF 读取、调度比对逻辑。 |
| `mz_yama.c` | 核心 Yama 算法：profile–profile DP 引擎。 |
| `mz_preyama.c` | 预处理模块：确定比对边界和重叠区域，为 Yama 准备数据。 |
| `maf.c` / `maf.h` | MAF 读写库，定义 `mafAli` 等结构。 |
| `mz_scores.c` | 打分矩阵和 gap 罚分参数的定义。 |
| `util.c` | 通用工具函数（内存管理、错误处理等）。 |

## 2. 对 pgr 项目的启示

1.  **Profile Alignment**: `multiz` 本质上是一个 Profile-Profile Aligner。在 `pgr` 的 `fas consensus` 或图比对中，如果涉及多序列合并，可以参考其 Sum-of-pairs 的处理方式。
2.  **锚点策略**: 利用共同参考序列作为锚点是处理大规模比对的有效手段（这也是 `pgr fas join` 的逻辑基础）。
3.  **MAF 处理**: `multiz` 的 `maf.c` 提供了一个轻量级的 MAF 解析参考实现。**更新**: `pgr` 现已实现完整的 MAF 读写功能（见第 3 节），采用了更符合 Rust 特性的设计，但字段定义上仍与 MAF 规范保持一致。

## 3. pgr MAF 实现现状 (2025-02)

目前 `pgr` 已在 `src/libs/fmt/` 下实现了完整的 MAF 读写支持。

### 3.1 模块分布与功能

*   **读写统一**：位于 `src/libs/fmt/maf.rs`
    *   **读取 (Reader)**:
        *   **核心函数**: `next_maf_block`, `parse_maf_block`。
        *   **特性**:
            *   支持 `a` (alignment) 和 `s` (sequence) 行解析，`a` 行 `score=` 字段已解析到 `MafAli.score`。
            *   **坐标转换**: `MafComp::to_range()` 将 MAF 的 0-based 坐标转换为 1-based inclusive 格式（如 `chr:start-end`）。
            *   **负链处理**: 自动处理负链坐标，将其转换为相对于正链的坐标范围。
    *   **写入 (Writer)**:
        *   **核心结构**: `MafWriter`。
        *   **特性**: 支持输出标准 MAF 头信息 (`##maf`) 和对齐块，自动处理列宽对齐。

### 3.2 数据结构

 读写使用同一套结构体：

 | 结构体 | 用途 | 关键字段 |
 | :--- | :--- | :--- |
 | `MafComp` | `s` 行（一条序列组件） | `src`、`start`、`size`、`strand`、`src_size`、`text`（均为 `String`/`usize`） |
 | `MafAli` | `a` 行 + block 内所有 `MafComp` | `score: Option<f64>`、`components: Vec<MafComp>` |

 ### 3.3 MAF 格式实现对比 (pgr vs multiz vs UCSC)

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

 ### 3.4 总结

 `pgr` 的 MAF 模块已经成熟，具备处理 UCSC MAF 格式的能力，并集成了坐标标准化功能，适合作为后续基因组比对分析的基础组件。
