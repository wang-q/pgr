# CLAUDE.md

此文件是我（AI 助手）在本仓库工作时的行为准则。所有规则都是硬性要求，除非用户明确覆盖。

## 语言规则

- **与用户交流**: 中文
- 本文件 (`CLAUDE.md`) : 使用中文编写
- **代码注释 (doc comments `///` `//!` 和行内 `//`)**: 英文
- **Git 提交信息**: 英文
- **文档正文** (如 `docs/*.md` 中的说明文字): 英文
- **Notes** (如 `notes/*.md` 中的供我自己看的说明文字): 中文

## 代码风格

**权衡：** 这些准则偏向谨慎而非速度。对于简单任务，自行判断。

### 编码前先思考

**不要假设。不要隐藏困惑。呈现权衡。**

在实现之前：
- 明确陈述你的假设。如果不确定，就问。
- 如果存在多种理解，把它们都列出来 — 不要默默地选一个。
- 如果有更简单的方法，说出来。必要时提出反对意见。
- 如果有不明白的地方，停下来。指出困惑之处。问。

### 简洁优先

**用最少的代码解决问题。不做任何推测性设计。**

- 不添加未被要求的功能。
- 不添加未被要求的"灵活性"或"可配置性"。
- 不为不可能发生的场景写错误处理。
- 如果你写了 200 行但其实 50 行就够了，重写它。

问自己："资深工程师会觉得这过于复杂吗？" 如果是，简化。

### 分层原则

**复杂逻辑放 `libs/`，`cmd_pgr/` 保持薄壳。**

- `src/libs/` 是复杂逻辑、算法、格式 I/O、共享工具的归宿。
- `src/cmd_pgr/` 仅负责：CLI 参数解析、参数转换、调用 `libs`、输出格式化。
- 单命令专用的复杂逻辑也放 `libs/`，即使当前只有一个消费者。
- 命令文件中内联的算法/业务逻辑应回迁 `libs/`。

判断标准：涉及算法、数据结构、复杂流程控制的代码属 `libs/`；只是 `clap` 参数 → 调用 → 打印的代码属 `cmd_pgr/`。

反例：在 `cmd_pgr/foo.rs` 里实现距离计算函数 → 应迁到 `libs/`。
正例：`cmd_pgr/foo.rs` 只做 `let args = parse(matches); let result = libs::foo::run(args); println!("{result}")`。

> 注："三次相似代码"原则针对的是重复代码的抽象提取，与本节的代码分层无关。

### 精准修改

**只改必须改的。只清理自己造成的混乱。**

编辑现有代码时：
- 不要"改进"相邻的代码、注释或格式。
- 不要重构没有坏的东西。
- 匹配现有风格，即使你不会这样写。
- 如果你注意到无关的死代码，提出来 — 不要删除它。

当你的修改产生了孤立代码时：
- 删除因你的修改而变得未使用的 import/变量/函数。
- 不要删除之前就存在的死代码，除非被要求。

检验标准：每一行改动都应该能追溯到用户的请求。

### 目标驱动执行

**定义成功标准。循环直到验证通过。**

将任务转化为可验证的目标：
- "添加验证" → "为无效输入写测试，然后让它们通过"
- "修复 bug" → "写一个能复现它的测试，然后让它通过"
- "重构 X" → "确保重构前后测试都通过"

对于多步骤任务，陈述简要计划：
```
1. [步骤] → 验证: [检查]
2. [步骤] → 验证: [检查]
3. [步骤] → 验证: [检查]
```

强有力的成功标准让你可以独立循环。薄弱的标准（"让它能用"）需要不断澄清。

### 必须遵守

- 每个 PR / commit 跑 `cargo fmt` 和 `cargo clippy -- -D warnings`，clean 之后再提交
- 公共 API (pub fn / pub struct / pub trait) 必须写 doc comment (英文，一行即可)
- 不写冗余注释 — 如果函数名和类型签名已经说明了行为，不要画蛇添足
- 用 `anyhow::Result<T>` 做函数返回值，`anyhow::bail!` / `anyhow::anyhow!` 构造错误

### 禁止

- 不要引入新依赖，除非用户明确要求
- 不要为了"可能"的未来需求写抽象 — 三次相似代码出现之后再考虑提取
- 不要写半成品实现 — stub / TODO 必须有明确的后续任务链接
- 不要用 `unsafe`，除非有充分理由且用户同意
- 不要写超过一行的 doc comment，除非是 trait 定义或复杂不变量
- 不要反向兼容的 shim（rename `_vars`、re-export 旧类型等）

## 项目概览

**当前状态**: 活跃开发中 | **主要语言**: Rust

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
    - **Sequences**: `fa` (FASTA), `fq` (FASTQ), `fas` (Block FA), `twobit` (命令 `2bit`), `gff`.
    - **Alignments**: `axt`, `chain`, `net`, `maf`, `paf`, `psl`, `lav`.
    - **Analysis**: `clust` (Clustering), `dist` (Distance), `mat` (Matrix), `nwk` (Phylogeny/Newick), `plot`.
    - **Misc**: `ms` (Simulation), `pl` (Pipelines).
- **`src/libs/`** - 共享工具库和核心逻辑。
  - **`fmt/`** - 格式 I/O: `fa`, `fas`, `fq`, `axt`, `maf`, `psl`, `lav`, `twobit`, `vcf`.
  - **`phylo/`** - 系统发育分析核心库。
    - **`node.rs`/`parser.rs`/`error.rs`**: 树节点定义、Newick 解析、错误类型（位于 `phylo/` 根级）。
    - **`cmp.rs`**: 树拓扑比较 (`TreeComparison` trait，Robinson-Foulds 等)。
    - **`tree/`**: 树算法与操作。
      - `ops.rs`/`algo.rs`: 节点操作 (add/remove/reroot/prune 等) 与算法。
      - `traversal.rs`/`query.rs`: 遍历 (pre/post/level-order) 与查询 (LCA/路径/距离/单系性)。
      - `stat.rs`/`balance.rs`/`distance.rs`/`support.rs`: 统计、平衡性指标、距离、支持值。
      - `io/`: 格式 I/O — Newick/DOT/SVG/Forest。
  - **`poa/`** - 偏序比对 (Partial Order Alignment) 实现。
  - **`chain/`** - 基因组比对链 (Chain/Net) 处理逻辑。
  - **`clust/`** - 聚类算法实现。
    - **`hier.rs`**: 层次聚类 (NN-chain 算法)。
    - **`dbscan.rs`, `mcl.rs`, `k_medoids.rs`**: 其他聚类算法。
    - **`nj.rs`, `upgma.rs`**: 建树算法 (Neighbor-Joining, UPGMA)。
    - **`tree_cut/`**: 树切分算法。
    - **`eval/`**: 聚类评估指标。
  - **`paf/`** - PAF 处理: 记录读写、查询、图构建、索引。
  - **`fasta/`** - FASTA 操作 (chunk/dedup/filter/stat)。
  - **`pairmat/`** - 配对距离矩阵。
  - **`io.rs`** - I/O 辅助函数。
  - **`alignment/`**, **`fas_multiz/`**, **`ms/`**, **`plot/`** 等 - 特定格式处理逻辑。

## 关键设计文档

- **`docs/`**: 用户面向命令文档（英文），每个 `pgr <command>` 对应一个 `docs/<command>.md`（子命令采用 `<command>-<subcommand>.md` 形式，如 `clust-cut.md`）；`docs/formats/` 为格式规范参考。注：`2bit` 命令的文档为 `docs/twobit.md`（历史命名）。
- **`notes/`**: 开发者面向笔记（中文）：`notes/design/`（设计稿/移植笔记）、`notes/references/`（外部工具源码分析）、`notes/` 根（项目理解/场景规划）。
- **`notes/project-understanding.md`**: 项目整体理解（架构、命令模块、核心库、现状评估、设计模式），含各文档的索引与定位，需要时查阅。

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
- **`coitrees`**: 区间树索引 (PAF/loc 模块核心)。
- **`intspan`**: 区间集合数据结构。
- **`petgraph`**: 图结构 (chain/paf 图构建)。
- **`indexmap`**: 保序 HashMap (名称→id 映射统一模式)。
- **`serde` + `bincode`**: PAF 索引 `.paf.idx` 持久化。

## 开发工作流

### 添加新命令

1.  在 `src/cmd_pgr/` 下相应的类别目录中创建新文件 (或新建目录)。
2.  在 `src/cmd_pgr/mod.rs` (或子目录的 `mod.rs`) 中声明该模块。
3.  在 `src/pgr.rs` 中注册该子命令。
4.  实现 `make_subcommand` 和 `execute`。
5.  添加测试文件 `tests/cli_<command>.rs`。

### 测试约定

- 集成测试位于 `tests/` 目录下，文件命名为 `cli_<command>.rs`。
- 测试数据通常放在 `tests/<command>/` 目录下。
- **推荐使用 `PgrCmd` 辅助结构体**（定义在 `tests/common/mod.rs`）来编写集成测试，以简化子进程调用和断言。
- 测试函数**不需要**返回 `anyhow::Result<()>`，也不需要以 `Ok(())` 结尾。直接在函数体中执行断言即可。
- 必须使用 `assert_cmd` 来定位二进制文件，以兼容自定义构建目录。
- **稳定性原则 (Zero Panic)**: 任何用户输入（包括畸形数据、二进制文件）都不应导致程序 Panic。必须捕获所有错误并返回友好的错误信息。
- **基准测试**: 性能敏感的变更必须伴随 `benches/` 下的基准测试结果（使用 `criterion`）。

### 辅助命令

**Changelog 生成**（以最新 tag 为起点）:
```bash
git tag | sort -V | tail -1          # find latest tag
git log v0.2.0..HEAD > gitlog.txt
git diff v0.2.0 HEAD -- "*.rs" "*.md" > gitdiff.txt
```

**Code coverage**:
```bash
rustup component add llvm-tools
cargo install cargo-llvm-cov
cargo llvm-cov
```

**WSL 构建**（避免 Windows 文件系统性能问题）:
```bash
mkdir -p /tmp/cargo
export CARGO_TARGET_DIR=/tmp/cargo
cargo build
```

## 帮助文本规范 (Help Text Style Guide)

- **`about`**: Third-person singular (e.g., "Counts...", "Calculates...").
- **`after_help`**: Uses raw string `r###"..."###`.
    - **Description**: Detailed explanation.
    - **Notes**: Bullet points starting with `*`.
        - Standard note for `fa`/`fas`: `* Supports both plain text and gzipped (.gz) files`
        - Standard note for `fa`/`fas`: `* Reads from stdin if input file is 'stdin'`
        - Standard note for `2bit`: `* 2bit files are binary and require random access (seeking)`
        - Standard note for `2bit`: `* Does not support stdin or gzipped inputs`
    - **Examples**: Numbered list (`1.`, `2.`) with code blocks indented by 3 spaces.
- **Arguments**:
    - **Input**: `infiles` (multiple) or `infile` (single).
        - Help: `Input [FASTA|block FA|2bit] file(s) to process`.
    - **Output**: `outfile` (`-o`, `--outfile`).
        - Help: `Output filename. [stdout] for screen`.
- **Terminology**:
    - `pgr fa` -> "FASTA"
    - `pgr fas` -> "block FA"
    - `pgr 2bit` -> "2bit"
