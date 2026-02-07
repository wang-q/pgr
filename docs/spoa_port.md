# Spoa 到 Rust 的移植计划

## 1. 目标
将 [Spoa](https://github.com/rvaser/spoa) (SIMD POA) 库从 C++ 移植到 Rust，并将其集成到 `pgr` 项目中。目标是提供 Rust 原生的高性能偏序比对 (POA) 和一致性序列生成能力。

**注意**: 根据性能评估，外部 `spoa` 二进制程序（C++ SIMD 实现）目前仍比纯 Rust 标量实现更快。因此，我们将 **保留** 对外部 `spoa` 的支持，并将其作为 `pgr fas consensus` 的一个选项（或默认选项），同时提供 Rust 原生实现作为可移植的替代方案。

## 2. 源码分析 (Spoa C++)
Spoa 经过高度优化，利用 SIMD (SSE4.1, AVX2) 进行动态规划。

*   **核心组件**:
    *   `Graph`: 有向无环图 (DAG)，表示偏序比对。存储节点 (碱基) 和加权边。
    *   `AlignmentEngine`: 比对算法的抽象基类。
    *   `AlignmentEngine` 子类: 处理特定的空缺罚分 (线性、仿射、凸函数) 和比对模式 (局部、全局、半全局)。
    *   `SIMD Implementation`: 使用内联函数 (intrinsics) 对 DP 矩阵进行并行处理 (垂直或对角线并行化)。

## 3. Rust 架构建议

### 3.1 模块结构
位置: `src/libs/poa/`

```text
src/libs/poa/
├── mod.rs          # 公共 API
├── poa.rs          # Poa 结构体封装
├── graph.rs        # Graph 类型定义 (基于 petgraph), NodeData
├── align.rs        # AlignmentEngine trait 及标量实现 (SISD)
├── consensus.rs    # 一致性序列生成逻辑
└── msa.rs          # 多序列比对 (MSA) 生成逻辑
```

### 3.2 数据结构

**Graph (`graph.rs`)**
*   **实现**: 复用 `petgraph::graph::DiGraph<NodeData, EdgeData>`。
*   `NodeData`: 存储碱基 (`u8`) 和 `aligned_to` (用于失配分支)。
*   `EdgeData`: 存储边权重 (`u32`)。
*   利用 `petgraph::algo::toposort` 进行拓扑排序。

**Alignment (`align.rs`)**
*   `AlignmentType`: 支持 Global (NW), Local (SW), SemiGlobal (OV)。
*   `AlignmentParams`: 仿射空缺罚分 (Match, Mismatch, GapOpen, GapExtend)。

### 3.3 SIMD 策略
本项目暂不追求 SIMD 优化，仅实现单线程标量版本。

**策略**:
1.  仅实现 **标量 (SISD)** 版本以确保正确性。
2.  性能优化仅限于算法层面和 Rust 语言层面的优化，不涉及 SIMD 指令集。

## 4. 实施阶段

### 第一阶段：核心数据结构与标量比对 (已完成)
*   [x] 基于 `petgraph` 定义 `Graph` 类型及 `NodeData`/`EdgeData`。
*   [x] 实现 `add_alignment` 逻辑（将序列对齐添加到图中）。
*   [x] 利用 `petgraph::algo::toposort` 验证拓扑排序功能。
*   [x] 实现一个基本的标量 `AlignmentEngine`，支持：
    *   [x] 全局 (NW)、局部 (SW)、半全局 (OV) 模式。
    *   [x] 仿射空缺罚分 (最常用)。
*   [x] 使用简单的测试用例对照 Spoa 进行验证。

### 第二阶段：一致性序列生成与集成 (已完成)
*   [x] 实现 `generate_consensus` (重束算法 / heaviest bundle algorithm)。
*   [x] 实现 `generate_msa` (从图生成多序列比对)。
*   [x] 封装 `Poa` 结构体以便于调用。
*   [x] **集成到 `pgr fas consensus`**:
    *   [x] 修改 `src/libs/alignment.rs`，增加 `get_consensus_poa_builtin`。
    *   [x] 修改 `src/cmd_pgr/fas/consensus.rs`，增加 `--engine <spoa|builtin>` 参数。
    *   [x] 默认行为：默认为 `builtin`，但允许用户指定 `spoa`。
    *   [x] 确保支持并行处理 (`--parallel`)。
*   [x] **集成到 `pgr fas refine`**:
    *   [x] 修改 `src/libs/alignment.rs`，`align_seqs` 支持 `builtin` (原 `poa`) 和 `spoa`。
    *   [x] 修改 `src/cmd_pgr/fas/refine.rs`，增加 `--msa <builtin|spoa>` 选项。
*   [x] 清理:
    *   [x] 删除临时的 `src/cmd_pgr/poa/` 模块及子命令。
    *   [x] 更新相关测试。

### 第三阶段（已移除）：SIMD 优化
*   本次移植不包含 SIMD 优化部分。

## 5. 当前状态

*   核心 POA 库 (`src/libs/poa`) 已完全实现，支持：
    *   全局、局部、半全局比对。
    *   仿射空缺罚分。
    *   **一致性序列生成** (Consensus)。
    *   **多序列比对生成** (MSA)。
*   `pgr fas consensus` 已集成双引擎支持 (`builtin` 和 `spoa`)。
*   `pgr fas refine` 已支持 `builtin` (默认) 和 `spoa` 选项，生成多序列比对。也支持 `clustalw`, `mafft`, `muscle`。
*   `builtin` 引擎的输出已验证与 `spoa` 一致。

## 6. 使用说明

### 一致性序列生成 (Consensus)

`pgr fas consensus` 命令用于生成一致性序列。

*   **输入**: Block FA 格式 (`.fas`)，支持 gzip。
*   **引擎选择**:
    *   `--engine builtin` (默认): 使用内置 Rust 实现，无需外部依赖。
    *   `--engine spoa`: 调用外部 `spoa` 命令（需在 PATH 中）。
*   **比对参数**:
    *   `--match <int>`: 匹配分 (默认: 5)
    *   `--mismatch <int>`: 失配罚分 (默认: -4)
    *   `--gap-open <int>`: Gap 开启罚分 (默认: -8)
    *   `--gap-extend <int>`: Gap 延伸罚分 (默认: -6)

### 多序列比对 (MSA)

`pgr fas refine` 命令用于重新比对 Block FA 文件中的序列。

*   **使用内置 POA**: `pgr fas refine input.fas --msa builtin` (默认)
*   **使用外部 Spoa**: `pgr fas refine input.fas --msa spoa`
*   **使用外部 ClustalW**: `pgr fas refine input.fas --msa clustalw`
*   **使用外部 Mafft**: `pgr fas refine input.fas --msa mafft`
*   **使用外部 Muscle**: `pgr fas refine input.fas --msa muscle`
*   **并行**: 支持 `--parallel <N>` 多线程加速。
*   **Outgroup**: 支持 `--outgroup` 选项正确处理外群。

示例:
```bash
# 使用内置引擎 (默认)，全局比对
pgr fas consensus tests/fas/example.fas

# 使用外部 spoa 引擎，局部比对
pgr fas consensus tests/fas/example.fas --engine spoa --algorithm local

# 自定义打分矩阵
pgr fas consensus tests/fas/example.fas -m 2 -n -3 -g -5 -e -1

# 并行处理
pgr fas consensus input.fas -p 4 -o output.fas
```

## 7. 参考资料
*   [Spoa GitHub](https://github.com/rvaser/spoa)
*   [Rust SIMD Guide](https://rust-lang.github.io/portable-simd/)
