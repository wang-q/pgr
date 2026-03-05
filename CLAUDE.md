# CLAUDE.md

此文件为 AI 助手在处理本仓库代码时提供指南与上下文。

## 项目概览

**当前状态**: 活跃开发中 | **主要语言**: Rust

**语言约定**: 为了便于指导，本文件 (`CLAUDE.md`) 使用中文编写，且**与用户交流时请使用中文**。但项目代码中的**所有文档注释 (doc comments)**、**行内注释**以及**提交信息**必须使用**英文**。

**目录约定**: 任何被 `.gitignore` 完全忽略的目录，均仅作为参考资料，**不是本项目的一部分**。

`pgr` (Practical Genome Refiner) 是一个多功能的基因组数据处理工具集。它旨在提供一套高效、易用的命令行工具，用于处理各种生物信息学格式（FASTA, FASTQ, AXT, MAF, PSL, Chain, Net, Newick 等）以及执行常见的基因组分析任务。

## 构建命令

### 构建

```bash
# 开发构建
cargo build

# 发布构建 (高性能)
cargo build --release
```

### 测试

```bash
# 运行所有测试
cargo test
```

## 架构

### 源代码组织

- **`src/pgr.rs`** - 主程序入口，负责命令行解析和分发。
    - 使用 `clap` 进行参数解析。
    - 在 `main` 函数中注册所有子命令模块。
- **`src/lib.rs`** - 库入口，导出模块。
- **`src/cmd_pgr/`** - 命令实现模块。按功能/格式分组：
    - **Sequences**: `fa` (FASTA), `fq` (FASTQ), `fas` (Block FA), `twobit` (2bit), `gff`.
    - **Alignments**: `axt`, `chain`, `net`, `maf`, `psl`, `lav`.
    - **Analysis**: `clust` (Clustering), `dist` (Distance), `mat` (Matrix), `nwk` (Phylogeny/Newick).
    - **Misc**: `ms` (Simulation), `pl` (Pipelines).
- **`src/libs/`** - 共享工具库和核心逻辑。
  - **`phylo/`** - 系统发育分析核心库。
    - **`tree/`**: 树结构定义、遍历、I/O、统计 (`stat.rs`)、切分 (`cut.rs`)。
    - **`algo/`**: 树操作算法（排序、reroot 等）。
  - **`poa/`** - 偏序比对 (Partial Order Alignment) 实现。
  - **`chain/`** - 基因组比对链 (Chain/Net) 处理逻辑。
  - **`clust/`** - 聚类算法实现。
    - **`hier.rs`**: 层次聚类 (NN-chain 算法)。
    - **`dbscan.rs`, `mcl.rs`, `k_medoids.rs`**: 其他聚类算法。
  - **`io.rs`** - I/O 辅助函数。
  - **`alignment.rs`**, **`fas_multiz.rs`**, **`psl.rs`** 等 - 特定格式处理逻辑。

## 关键设计文档

### Analysis (Phylogeny, Clustering, Matrix)
- **`docs/phylo.md`**: 系统发育树核心数据结构设计 (Arena vs Pointer)。
- **`docs/nwk-cut.md`**: 树切分与 Scan 模式。
- **`docs/nwk-eval.md`**: 树结构多维度评估（设计中）。
- **`docs/clust.md`**: 聚类模块总体规划。
- **`docs/clust-hier.md`**: 层次聚类算法与实现细节。
- **`docs/clust-eval.md`**: 通用聚类评估（设计中）。
- **`docs/dist.md`**: 距离计算与度量。
- **`docs/mat.md`**: 距离矩阵操作与转换。

### Algorithms & Ports
- **`docs/spoa_port.md`**: SPOA (SIMD POA) 移植笔记。
- **`docs/ms2dna_port.md`**: MS 模拟器移植笔记。
- **`docs/multiz.md`**: Multiz 多序列比对格式处理。

### General
- **`docs/axt.md`**: AXT 格式操作与转换。
- **`docs/chain.md`**: Chain 格式高级处理。
- **`docs/fa.md`**: FASTA 格式全能操作。
- **`docs/fas.md`**: Block FASTA (FAS) 格式操作与转换。
- **`docs/formats.md`**: 支持的文件格式概览。
- **`docs/fq.md`**: FASTQ 格式操作。
- **`docs/gff.md`**: GFF 格式操作。
- **`docs/lav.md`**: LAV 格式操作与转换。
- **`docs/maf.md`**: MAF 格式操作与转换。
- **`docs/ms.md`**: Hudson's ms 模拟器数据处理。
- **`docs/net.md`**: Net 格式操作与转换。
- **`docs/nwk.md`**: Newick 系统发育树操作与可视化。
- **`docs/pl.md`**: 集成分析流程 (Pipelines)。
- **`docs/plot.md`**: 各种生物数据可视化工具。
- **`docs/psl.md`**: PSL 格式操作与转换。
- **`docs/twobit.md`**: 2bit 格式操作。

## 命令结构 (Command Structure)

每个命令在 `src/cmd_pgr/` 下作为一个独立的模块实现，通常包含两个公开函数：

1.  **`make_subcommand`**: 定义命令行接口。
    -   返回 `clap::Command`。
    -   使用 `.about(...)` 设置简短描述 (第三人称单数)。
    -   推荐使用 `.after_help(...)` 提供详细帮助信息。
2.  **`execute`**: 命令执行逻辑。
    -   接收 `&clap::ArgMatches`。
    -   返回 `anyhow::Result<()>`。

### 关键依赖

- **`noodles`**: 处理 FASTA, FASTQ, BGZF, GFF 等标准格式。
- **`clap`**: 命令行参数解析。
- **`anyhow`**: 错误处理。
- **`rayon`**: 并行计算。
- **`nom`**: 文本解析 (Newick 等)。
- **`regex`**: 正则表达式。

## 开发工作流

### 添加新命令

1.  在 `src/cmd_pgr/` 下相应的类别目录中创建新文件 (或新建目录)。
2.  在 `src/cmd_pgr/mod.rs` (或子目录的 `mod.rs`) 中声明该模块。
3.  在 `src/pgr.rs` 中注册该子命令。
4.  实现 `make_subcommand` 和 `execute`。
5.  添加测试文件 `tests/cli_<command>.rs`。

### 测试约定

- 集成测试位于 `tests/` 目录下，文件命名为 `cli_<command>.rs`。
- 测试数据通常放在 `tests/data/<command>/` 目录下。
- **推荐使用 `PgrCmd` 辅助结构体**（定义在 `tests/common/mod.rs`）来编写集成测试，以简化子进程调用和断言。
- 测试函数**不需要**返回 `anyhow::Result<()>`，也不需要以 `Ok(())` 结尾。直接在函数体中执行断言即可。
- 必须使用 `assert_cmd` 来定位二进制文件，以兼容自定义构建目录。
- **稳定性原则 (Zero Panic)**: 任何用户输入（包括畸形数据、二进制文件）都不应导致程序 Panic。必须捕获所有错误并返回友好的错误信息。
- **基准测试**: 性能敏感的变更必须伴随 `benches/` 下的基准测试结果（使用 `criterion`）。

## 代码规范

- 使用 `cargo fmt` 格式化代码。
- 使用 `cargo clippy` 检查潜在问题。
- 优先使用标准库和项目中已引入的 crate。
- 保持代码简洁，注重性能。

## 帮助文本规范 (Help Text Style Guide)

- **About**: 第三人称单数动词 (e.g., "Calculates...", "Converts...").
- **Args**:
    - Input: `infile` / `infiles`.
    - Output: `outfile` (`-o`).
- **Description**: 简明扼要，解释命令的核心功能。
- **Examples**: 提供典型的使用示例。
