# CLAUDE.md

此文件为 AI 助手在处理本仓库代码时提供指南与上下文。

## 项目概览

**当前状态**: 活跃开发中 | **主要语言**: Rust

**语言约定**: 为了便于指导，本文件 (`CLAUDE.md`) 使用中文编写，且**与用户交流时请使用中文**。但项目代码中的**所有文档注释 (doc comments)**、**行内注释**以及**提交信息**必须使用**英文**。

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
  - **`phylo/`** - 系统发育分析核心库 (Tree, Node, Algo)。
  - **`poa/`** - 偏序比对 (Partial Order Alignment) 实现。
  - **`chain/`** - 基因组比对链 (Chain/Net) 处理逻辑。
  - **`clust/`** - 聚类算法 (DBSCAN, MCL, K-Medoids)。
  - **`io.rs`** - I/O 辅助函数。
  - **`alignment.rs`**, **`fas_multiz.rs`**, **`psl.rs`** 等 - 特定格式处理逻辑。

### 命令结构 (Command Structure)

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

- 集成测试位于 `tests/` 目录下，命名格式为 `cli_<command>.rs`。
- 测试数据位于 `tests/` 下对应的子目录中 (e.g., `tests/axt/`, `tests/newick/`)。
- 推荐使用 `assert_cmd` 进行 CLI 测试。
- **Newick 测试数据**: 集中在 `tests/newick/`，使用 `.nwk` 扩展名。

## 代码规范

- 使用 `cargo fmt` 格式化代码。
- 使用 `cargo clippy` 检查潜在问题。
- 优先使用标准库和项目中已引入的 crate。
- 保持代码简洁，注重性能。

### Hash 算法选择

- **HashMap**: 使用 `ahash` (默认) 或 `fxhash` 用于高性能内存哈希。
- **One-shot Hashing**: 使用 `rapidhash` 用于非加密的高速哈希 (如去重)。
- **Stable/Crypto**: 使用 `SipHash` 或 `SHA-256` 当需要稳定或安全哈希时。

## 帮助文本规范 (Help Text Style Guide)

- **About**: 第三人称单数动词 (e.g., "Calculates...", "Converts...").
- **Args**:
    - Input: `infile` / `infiles`.
    - Output: `outfile` (`-o`).
- **Description**: 简明扼要，解释命令的核心功能。
- **Examples**: 提供典型的使用示例。
