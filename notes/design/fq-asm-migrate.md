# fq / asm 命令组迁移到 anchr 方案

> 2026-08-12 决策定稿（用户方向）：reads 处理（`fq`）与组装（`asm`）业务
> 迁到 anchr（组装器定位），FASTA/FASTQ **读取与 Phred 编码基础留在 pgr**
> （与 FA 读取不迁移同理）。anchr 依赖 pgr crate（基础库），跨项目依赖
> 可接受（均自有项目）。

## 1. 决策背景

pgr 定位"通用基因组数据处理工具集"，anchr 定位"染色体级组装流程编排器"。
`fq`（reads 处理）与 `asm`（组装）业务逻辑主要服务 anchr 流程，归位
anchr 更内聚，同时给 pgr 瘦身（约 1.5 万行代码 + 21 个测试文件）。

关键前提：**格式读取与 Phred 编码是基础层，不属于任何命令的"业务"**——
FASTA/FASTQ 读入（`libs/fmt/`）与 FA 读取一样永远留在 pgr；Phred 编码
转换/检测（`libs/fq/qual` + `detect_quality_base`）同样留 pgr（被
`kmer qhist` 与 `asm` 使用）。

## 2. 分层原则

| 层 | 内容 | 归属 |
| :--- | :--- | :--- |
| 基础层 | 格式 I/O（FASTA/FASTQ 读入 `libs/fmt/`）、Phred 编码（qual）、k-mer 基础、PAF 解析、io/ds/loc/sys | **留 pgr** |
| 业务层 | `fq` 清洗/合并/纠错/归一/采样/分块算法、`asm` 组装算法（unitig/contig/OLC/map） | **迁 anchr** |
| 命令壳 | `cmd_pgr/fq/`、`cmd_pgr/asm/` 的 clap 解析与执行 | **迁 anchr** |

## 3. 迁移清单（逐文件）

### 3.1 留在 pgr（基础层）

| 文件 | 说明 |
| :--- | :--- |
| `libs/fmt/seq.rs` / `fa.rs` / `fq.rs` | FASTA/FASTQ 读入与写出（SeqReader 等） |
| `libs/fq/qual.rs` | `base_to_number` / `to_phred` / `from_phred`（Phred 编码转换） |
| `libs/fq/trim.rs` 的 `detect_quality_base` | 质量编码检测（**从 trim.rs 抽到基础层**，如 `libs/fq/qual.rs` 或独立模块） |
| `libs/fq/pairs.rs` | `PairReader`（双端读入，属读取能力） |
| `libs/kmer/`、`libs/paf/`、`libs/io.rs`、`libs/ds/`、`libs/loc.rs`、`libs/sys*` | 被 fq/asm 依赖的基础库 |

### 3.2 迁到 anchr（业务层 + 命令壳）

| 文件 | 说明 |
| :--- | :--- |
| `libs/fq/`（除 qual.rs/pairs.rs/trim.rs 的 detect_quality_base） | trim/trim_adapter/merge/norm/sample/split/clump/bbnet/overlap 业务 |
| `libs/asm/`（assemble.rs/tadpole.rs/mod.rs） | tadpole 组装、unitig 逻辑 |
| `libs/olc/`（consensus.rs/layout.rs/overlap.rs/mod.rs） | OLC 三阶段 |
| `cmd_pgr/fq/`（16 文件） | 命令壳 |
| `cmd_pgr/asm/`（9 文件） | 命令壳 |
| `tests/cli_fq*.rs`（17 文件）+ `tests/cli_asm*.rs`（4 文件） | 集成测试随迁 |
| `docs/fq*.md`、`docs/asm*.md` | 用户文档随迁 |

### 3.3 依赖关系（迁移后）

```
anchr（fq/asm 业务 + 命令壳）
  └─ 依赖 pgr crate（libs::fmt / fq::qual / kmer / paf / io / ds / loc / sys）
pgr（格式/基础 + 比对/分析/索引/遮蔽等其余命令组）
```

`fq` 与 `asm` 在 anchr 内互相依赖（`fq/merge` 用 `asm::tadpole`、
`asm` 用 `fq::qual`），同一 crate 内无碍。pgr 侧 `kmer qhist` 的
`detect_quality_base` 依赖由基础层满足，不随迁。

## 4. 成本与风险

1. **pgr libs pub 化**：目前 `lib.rs` 不导出 fq/asm 依赖的基础模块；
   anchr 依赖 pgr crate 需将 `fmt`/`fq`（qual 部分）/`kmer`/`paf`/`io`/
   `ds` 等设为 `pub`，并承担库 API 版本维护（pgr 从"二进制 + 内部库"
   变为"对外库 + 二进制"）；
2. **`detect_quality_base` 抽取**：从 `libs/fq/trim.rs` 移到基础层，
   确认无循环依赖（trim 业务迁走后，qual 基础不依赖业务）；
3. **版本同步**：pgr 基础库演进需与 anchr 协调（Cargo 依赖版本、breaking
   变更沟通）；可用 workspace 或版本约束缓解；
4. **测试/文档/CI 迁移**：21 个测试文件 + fq/asm 文档 + CI 步骤随迁；
5. **CHANGELOG/help 一致性**：pgr 的 `tests/cli_consistency.rs`（after_help
   命令集合核对）需同步移除 fq/asm 条目；anchr 侧新增对应核对。

## 5. 实施步骤（分阶段）

采用**双轨迁移**策略（用户定稿）：pgr 全程保持完整可用，迁移增量推进、
可回退；只有在两个目录核对确认后才从 pgr 删除。

### 阶段 1：pgr libs pub 化（不删除任何代码）

目标：把 anchr 迁移需要的 pgr 基础模块设为公开 API，pgr 二进制与命令
行为完全不变。

- `lib.rs` 展开 `pub mod`：`fmt`（seq/fa/fq）、`fq`（qual/pairs 基础）、
  `kmer`、`paf`、`io`、`ds`、`loc`、`sys` 等被 fq/asm 依赖的模块；
- `detect_quality_base` 从 `libs/fq/trim.rs` 抽到基础层（如
  `libs/fq/qual.rs`），`kmer qhist` 改引用新位置；
- 只做 pub 化与抽取，**不删代码、不改命令注册**；`cargo test` 全绿
  （行为不变）；
- pgr 版本作为 anchr 的依赖（path 依赖开发期，定稿后切 git/版本）。

**2026-08-12 已完成**：

- `detect_quality_base` + `PHRED33`/`PHRED64` 从 `libs/fq/trim.rs` 抽到
  `libs/fq/qual.rs`（基础层），`kmer qhist`/`fq s-filter`/`trim` 内部引用
  更新，4 个 detect 测试随迁（qual.rs 内）；
- `kmer::base_codes`、`kmer::count::count_keys` 由 `pub(crate)` 改 `pub`
  （fq norm 迁移后 anchr 需要调用）；
- 新增 `tests/migrate_api.rs`（integration test 以外部 crate 视角编译，
  `use` 全部解析 = anchr 依赖 pgr crate 可用全部基础符号）；
- 未删除任何代码、未改命令注册；`cargo test` 1755 全绿，fmt/clippy 干净。

### 阶段 2：anchr 逐步迁移（双轨）

在 anchr crate 内逐命令/逐模块移植 fq/asm（业务 libs + 命令壳 + 测试 +
docs），依赖 pgr crate 的基础模块：

- 迁移顺序建议：先 `fq` 基础命令（trim/sample/filter/split/merge）→
  `fq` 纠错（ec-kmer/ec-overlap/extend）→ `asm`（unitig/contig）→
  `asm` OLC（ovlp/layout/cns/olc）→ `asm map`；
- 每迁一个命令：代码进 anchr（`use pgr::libs::...` 引用基础），对应测试
  随迁并全绿；pgr 侧代码**原样保留**；
- anchr 与 pgr 的 fq/asm 命令在此阶段**同时存在**（双轨），anchr 流程可
  逐步切到 anchr 侧命令验证。

### 阶段 3：核对两个目录

对每个迁移完成的命令，核对 pgr 与 anchr 的实现一致：

- 同一测试数据跑两边，输出**逐字节一致**（golden 对照，参照
  `anchr-trim-replace.md` / `anchr-merge-replace.md` 的既有对照流程）；
- 统计/统计输出一致（如覆盖量、unitig 数、PSL 行数）；
- anchr 侧集成测试全绿 + pgr 侧对应测试仍全绿（双轨期）。

核对不通过的命令**不进入删除阶段**，回 anchr 侧修或评估差异。

### 阶段 4：从 pgr 去除（仅限核对确认完成的）

对核对通过的命令/模块，从 pgr 删除：

- 移除 `cmd_pgr/fq`、`cmd_pgr/asm` 对应命令文件 + `pgr.rs` 注册 +
  `cmd_pgr/*/mod.rs` 声明；
- 移除 `libs/fq`（业务部分）、`libs/asm`、`libs/olc` 中已迁移且无 pgr
  内消费者的代码；
- 移除随迁测试文件、docs；更新 `cli_consistency.rs`（after_help 命令集合）、
  project-understanding 索引、todo、CHANGELOG；
- 逐批删除、逐批 `cargo test`，保持剩余命令全绿。

**2026-08-13 已完成**：

- 命令壳：`cmd_pgr/fq`（15 命令）、`cmd_pgr/asm`（8 命令）删除，`pgr.rs`
  注册/分发/after_help 移除，`cmd_pgr/mod.rs` 移除声明；
- 业务 libs：`libs/asm`、`libs/olc`、`libs/map` 删除，`libs/fq` 仅留
  qual.rs/pairs.rs/mod.rs（基础层）；
- 测试/文档：`tests/cli_fq*.rs`/`cli_asm*.rs`（21 个）删除，docs
  `fq.md`/`asm.md` 删除，`benches/fq_assemble_benchmark.rs`/
  `asm_map_benchmark.rs` 及 Cargo.toml bench 定义删除；
  `tests/bbtools/`（Lambda 数据）**保留**——kmer 命令组测试仍在用
  （`cli_kmer.rs` 5 个真实数据测试），非 fq/asm 专属；
- 索引：project-understanding §3.1/§4.5/§6.1/§10/§11 移除 fq/asm 内容；
  todo/khmer 引用改指向 anchr；CHANGELOG 记录；docs/sam.md、dist.md、
  usage_examples.md、sam to_rg/ihist doc comment 的 `pgr fq/asm` 示例改
  `anchr`；
- 参考笔记随迁：references 9 个（bcalm/canu/celera/cutadapt/fairy/
  metaMDBG/quorum/sickle/skesa）+ design 6 个（anchr-trim/merge-replace、
  fq-assemble/asm-map/fq-index/olc）从 pgr 移除；
- 验证：`cargo test` 全绿，`rg libs::(asm|olc|map)|cmd_pgr::(fq|asm)`
  无残留，fmt/clippy 干净；阶段 1 pub 化基础层（fmt/fq::qual/pairs/kmer/
  paf/io/ds/loc/sys）保留供 anchr 依赖。

### 阶段 5：收尾

- 更新 `references/anchr.md` 边界划分（fq/asm 业务已在 anchr）、todo 状态；
- 确认 pgr 剩余命令组不依赖已删 libs（`rg libs::fq|asm|olc` 无残留）；
- anchr 流程端到端冒烟（anchr fq/asm + pgr 基础 + 模板）。

## 6. 影响面小结

- pgr 减少：约 1.5 万行（fq/asm libs + cmd）+ 21 测试 + fq/asm 文档；
- pgr 保留：格式读入/写出（fa/fas/fq 基础）、Phred 编码、k-mer、PAF、
  比对（align）、分析（dist/kmer/rept/runlist）、索引（pgi/pbit）、
  泛基因组（paf）、模拟/流程/可视化；
- anchr 增加：reads 处理 + 组装全套命令，依赖 pgr 基础库。

## 7. 阶段 2 操作手册（anchr 侧迁移）

> 用户选择"只写文档"路线：pgr 侧代码不动（除阶段 1 已完成的 pub 化），
> 本手册供在 anchr 仓库（`~/Scripts/anchr`）照单操作。

### 7.1 anchr crate 准备

1. **Cargo.toml** 加依赖：
   ```toml
   [dependencies]
   pgr = { path = "../pgr" }   # 开发期 path；定稿后切 git/版本
   ```
2. **`src/lib.rs`** 加模块（与现有 libs 并列）：
   ```rust
   pub mod asm;
   pub mod fq;
   pub mod olc;
   ```
3. **`src/cmd/mod.rs`** 加子命令模块（建议子目录，保持 pgr 组织方式）：
   ```rust
   pub mod fq;
   pub mod asm;
   ```
   并在 `anchr.rs`（clap dispatch）注册 `fq`/`asm` 顶层命令组。

### 7.2 文件迁移映射

**库（libs）**——从 `pgr/src/libs/` 复制到 `anchr/src/libs/`：

| 源（pgr） | 目标（anchr） | 说明 |
| :--- | :--- | :--- |
| `fq/bbnet.rs` | `fq/bbnet.rs` | 业务 |
| `fq/clump.rs` | `fq/clump.rs` | 业务 |
| `fq/merge.rs` | `fq/merge.rs` | 业务（依赖 `asm::tadpole`） |
| `fq/norm.rs` | `fq/norm.rs` | 业务（依赖 `kmer`，pgr 基础） |
| `fq/overlap.rs` | `fq/overlap.rs` | 业务 |
| `fq/sample.rs` | `fq/sample.rs` | 业务 |
| `fq/split.rs` | `fq/split.rs` | 业务 |
| `fq/trim.rs` | `fq/trim.rs` | 业务（`detect_quality_base` 已抽走，无 qual 定义） |
| `fq/trim_adapter.rs` | `fq/trim_adapter.rs` | 业务 |
| `asm/assemble.rs` | `asm/assemble.rs` | 业务 |
| `asm/tadpole.rs` | `asm/tadpole.rs` | 业务（依赖 `pgr::fq::qual`） |
| `olc/consensus.rs` | `olc/consensus.rs` | 业务 |
| `olc/layout.rs` | `olc/layout.rs` | 业务 |
| `olc/overlap.rs` | `olc/overlap.rs` | 业务 |
| `fq/mod.rs` / `asm/mod.rs` / `olc/mod.rs` | 同名 mod.rs | 去掉 `pairs`/`qual` 声明（留在 pgr） |

**命令（cmd）**——从 `pgr/src/cmd_pgr/` 复制到 `anchr/src/cmd/`：

| 源（pgr） | 目标（anchr） | 命令 |
| :--- | :--- | :--- |
| `fq/clean.rs` 等 15 个（clean/clump/ec_kmer/ec_overlap/extend/filter/interleave/merge/norm/range/s_filter/sample/split/to_fa/trim_qual） | `cmd/fq/*.rs` | `fq` 组 |
| `asm/cns.rs` 等 8 个（cns/contig/layout/map/olc/ovlp/unitig + common.rs） | `cmd/asm/*.rs` | `asm` 组 |

**测试**——从 `pgr/tests/` 复制：`cli_fq*.rs`（17 个）+ `cli_asm*.rs`（4 个）；
libs 内嵌 `#[cfg(test)]`（qual.rs 的 detect 测试留 pgr，其余随 libs 走）。

### 7.3 引用改写规则

1. **pgr 基础库**：`crate::libs::fmt/...`、`crate::libs::kmer/...`、
   `crate::libs::nt::rev_comp`、`crate::libs::ds::radix_sort::...`、
   `crate::libs::par::ordered_map`、`crate::libs::sys::mem_cap`、
   `crate::libs::io::...` → 改 `pgr::libs::...`（阶段 1 已 pub，可直接引用）；
2. **fq/asm/olc 内部**（迁移后同 crate）：`crate::libs::fq/asm/olc::...`
   保持 `crate::libs::...` 不变（anchr 内同路径）；
3. **`crate::cmd_pgr::args::*`**（get_outfile/outfile_arg/infiles_arg/
   ensure_outfile_distinct/parse_parallel_auto/collect_ranges 等）：
   anchr 复制这批 clap helper 到 `src/utils.rs` 或新建 `src/cmd/args.rs`，
   引用改为 `crate::cmd::args::...`；
4. **`crate::cmd_pgr::kmer::qhist::{qual_thresh_arg, bits_arg}`**
   （`fq/s_filter.rs` 用）：两个小 helper 复制到 anchr（内联或放 cmd/fq 公共处）；
5. **`use pgr::libs::fq::qual`**：`asm/tadpole.rs`、`fq/*` 中引用
   `crate::libs::fq::qual::{from_phred, to_phred}` → `pgr::libs::fq::qual::...`
   （qual 留在 pgr 基础层）；
6. **`use pgr::libs::fq::pairs::PairReader`**：`fq` 命令/`fmt` 依赖
   → `pgr::libs::fq::pairs::PairReader`（pairs 留 pgr）。

### 7.4 迁移批次与顺序（依赖序）

1. **批 1：`asm` libs**（tadpole/assemble）+ `olc` libs（consensus/layout/
   overlap）——只依赖 pgr 基础（fmt/kmer/paf/qual），无 fq 业务依赖；
2. **批 2：`fq` libs 基础**（trim/trim_adapter/sample/split/overlap）——
   只依赖 pgr 基础（fmt/qual；overlap.rs 无 crate::libs 依赖；trim 已无
   qual 定义）；
3. **批 3：`fq` libs 进阶**（merge→依赖批 1 tadpole；norm/clump/bbnet）——
   norm 依赖 pgr `kmer`（base_codes/count_keys 已 pub）；
4. **批 4：`fq` cmd**（15 个）+ `asm` cmd（8 个）——cmd 壳 + 复制 args
   helper；
5. **批 5：测试随迁**（cli_fq*/cli_asm* + libs 内嵌）——每批迁移时对应
   测试一并搬。

### 7.5 每批验证（双轨 golden）

- 同一测试数据分别跑 pgr 与 anchr 的对应命令，输出 `diff` **逐字节一致**
  （参照 `anchr-trim-replace.md` / `anchr-merge-replace.md` 的既有对照
  流程：golden 文件 + 统计对比）；
- anchr 侧 `cargo test` 对应测试全绿；pgr 侧对应测试**仍全绿**（双轨期
  不删 pgr 代码）；
- 统计一致：覆盖量 / unitig 数 / PSL 行数 / 直方图 bin 等。

### 7.6 阶段 3 核对清单（删除前逐命令确认）

对每个迁移命令：

- [ ] pgr / anchr 同输入输出逐字节一致（golden diff 为空）；
- [ ] anchr 侧测试全绿，pgr 侧测试仍全绿；
- [ ] 边界统计一致（长度/数量/覆盖）；
- [ ] 确认 pgr 内无残留消费者（`rg libs::fq|asm|olc` 仅剩基础层与
  已删除项）。

核对通过后按阶段 4 从 pgr 删除。

## 8. 参考笔记迁移清单（fq/asm 专属文档随迁）

很多 `notes/` 参考文档只为 fq/asm 提供背景（外部工具源码分析、移植记录），
与 pgr 其余部分无关；fq/asm 迁到 anchr 后这些笔记应一并迁走（作为 anchr
项目的笔记），pgr 只留通用/其他模块的文档。

### 8.1 迁到 anchr（fq/asm 专属，9 个 references）

| 文档 | 服务对象 |
| :--- | :--- |
| `bcalm.md` | `asm unitig` 移植来源 |
| `canu.md` | `asm olc` 参考（OLC 设计） |
| `celera.md` | `asm olc` 参考 |
| `cutadapt.md` | `fq trim-qual`/`clean` 参考 |
| `fairy.md` | `fq norm` 大数据方案调研（k-mer 采样路线） |
| `metaMDBG.md` | OLC v1 素材（`asm olc`） |
| `quorum.md` | read 纠错参考（`fq ec-kmer` + anchr `quorum` 命令） |
| `sickle.md` | `fq trim-qual` 算法来源 |
| `skesa.md` | OLC v1 素材（`asm olc`） |

### 8.2 迁到 anchr（design 移植记录，6 个）

| 文档 | 服务对象 |
| :--- | :--- |
| `anchr-trim-replace.md` | `fq` trim 系列移植（BBTools 对照） |
| `anchr-merge-replace.md` | `fq` merge/ec 系列移植（BBTools 对照） |
| `fq-assemble.md` | `asm` contig/unitig/olc |
| `asm-map.md` | `asm map` |
| `fq-index.md` | `fq range`（FASTQ `.loc` 索引） |
| `olc.md` | `asm ovlp/layout/cns/olc` |

### 8.3 迁到 anchr（audit，1 个）

| 文档 | 说明 |
| :--- | :--- |
| `audit-fq.md` | fq 命令族审计记录 |

### 8.4 留在 pgr（28 个 references）

`agc-cpp`（pbit）、`alnfill`（align fill/rest）、`app-egaz`（align/
chainnet）、`biser`（SD）、`cactus`/`cactus_lastz`（泛基因组）、`fastga`
（pgi）、`fastk`（k-mer）、`gfa`（格式）、`hv`（距离）、`impg`（泛基因
组）、`kaks`、`khmer`（k-mer）、`merqury-fk`（k-mer）、`minigraph`（泛
基因组）、`mosdepth`（rg coverage）、`multiz`（fas_multiz）、`ntsynt`、
`pangenome-tools`、`repeatmasker`（rept）、`ropebwt3`、`seqwish`、
`smoothxg`、`spoa`（poa）、`syng`（syncmer）、`ucsc`（chain/net）、
`wgatools`（paf）、`anchr.md`（边界划分档案，留 pgr 作迁移记录）。

### 8.5 迁移后 pgr 侧索引更新

- `project-understanding.md` §10/§11/§12：移除迁走的 design/references
  条目（8.1/8.2/8.3），保留 §7 本迁移方案档案；
- `todo.md` 中引用迁走文档的条目改指向 anchr（如 anchr 模板替换的 golden
  对照流程引用 `anchr-trim-replace.md` → 在 anchr 侧）；
- 基准（benchmarks/）均为 pgr 通用（kmer-throughput 属 kmer 命令，留）；
  `fq-asm-migrate.md` 本身作为迁移档案留在 pgr，迁移完成后归档。
