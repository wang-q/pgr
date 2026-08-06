# 复刻 RepeatMasker 的 RMBlast 路径：`pgr rept rmblast` 设计

> 2026-08-07。目标：按 RepeatMasker 4.2.4 的 `-lib` 搜索路径**逐参数复刻**，
> 调用外部 `makeblastdb` / `rmblastn`，只输出 runlist 区间（直接喂
> `pgr fa mask`），不做 family/class 注释与 ProcessRepeats 后处理。
> TRF 阶段已有 `pgr rept trf`，本命令不重复（对应 RM 的 `-nolow` 视角）。

> **状态（2026-08-07 更新）**：`pgr rept rmblast` 已实现落地——
> `src/libs/rmblast.rs`（矩阵/参数/tab 解析 + 单测）、
> `src/libs/pl/repeat.rs::run_rmblast_repeat_pipeline`、
> `src/cmd_pgr/rept/rmblast.rs`。两处复刻修正后与 RM 默认输出几乎完全一致：
> ① minscore 直接用 cutoff=225（早期误读了 `runTestStage` 的 7.5% 折扣，
> 正式 `runStage` 不打折）；② 默认 `--frag 60000` 复刻 RM 分片 + 每片段
> GC 选矩阵。10 株 × TnCentral：RM 原始 IS 被我们覆盖 100.0%、我们被 RM
> 覆盖 99.7%，RM-strict 100%（残差是 RM ProcessRepeats 的元件端点裁剪，
> 见 §2.8）。

## 1. 复刻对象：RM 的 `-lib` 单阶段路径

源码依据（`RepeatMasker/` 4.2.4）：

* `runSearchStages`（RepeatMasker:3537）：用户给 `-lib` 时**只跑一个搜索阶段**，
  `general_search_parameters`（RepeatMasker:3639），不涉及物种配方；
* `runStage`（RepeatMasker:2560 起）：完成 GC→矩阵、minscore 折扣、参数传递；
* `NCBIBlastSearchEngine.pm` `getParameters()`（:439 起）：把配方翻译成 rmblastn
  命令行；
* `search()`（RepeatMasker:1896）：组装并调用，含失败重试（带宽/word_size 调整）。

## 2. 逐参数对应表（默认档）

配方 `general_search_parameters` → rmblastn 参数：

| 配方字段 | 值 | rmblastn 参数 | 最终值 | 备注 |
| :--- | :--- | :--- | :--- | :--- |
| minscore | 225 | `-min_raw_gapped_score` | **225** | 正式 `runStage` 直接用配方值；7.5% 折扣只在未使用的 `runTestStage` 里（早期误读，已修正） |
| minmatch | [8,9,11,13] | `-word_size` | 9 | 默认档取 index 1；`-s`=8 / `-q`=11 / `-qq`=13 |
| gap_initValue | −30 | `-gapopen` | 24 | abs(−30−(−6)) |
| ins_gap_extValue | −6 | `-gapextend` | 6 | abs(−6)；del_gap_ext −5 仅 crossmatch 用 |
| bandwidth | 14 | xdrop 三件套 | 450 / 225 / 112 | 正带宽走 MaskerAid 默认：ungap=ms×2、gap_final=ms、gap=int(ms/2)，ms=225 |
| matrix | 20p##g.matrix | `-matrix` + `BLASTMAT` | 20p<GC>.matrix | GC 选择见 §3；BLASTMAT=Matrices/ncbi/nt |
| masklevel | 101 | `-mask_level` | 101 | runStage 硬编码 101（非配方 90） |
| raw=0 | — | `-complexity_adjust` | 开 | |
| — | — | `-dust` | no | 低复杂度交给 masklevel，不是 blast dust |
| — | — | `-num_alignments` | 9999999 | 不做默认 500 条截断 |
| — | — | `-num_threads` | 4 | 引擎默认；RM 主程序不设（`-pa` 控制 batch 并行） |
| — | — | `-db` / `-query` | 库 / 基因组 batch | query=基因组，subject=库 |
| — | — | `-outfmt` | `6 score perc_sub perc_query_gap perc_db_gap qseqid qstart qend qlen sstrand sseqid sstart send slen kdiv cpg_kdiv transi transv cpg_sites` | 18 列 tab；RM 因 `setGenerateAlignments(1)` 再附 qseq/sseq，对 hit 集无影响，我们不输出 |

建库：`makeblastdb -dbtype nucl -in <lib> -out <db>`（RepeatMasker:6549）。

分片（仅复刻时）：fragmentSize=60000、overlapLen=2000（RepeatMasker:629-638），
`-frag` 可覆盖且必须 ≥ 2×overlap；hit 坐标按 batch 偏移映射回基因组。

## 3. GC → 矩阵选择（runStage / chooseMatrices）

* `GC_frac`：默认 43；`-gccalc` 或 batch 为单序列且 >2000 bp 时用 batch 平均 GC
  （RepeatMasker:2561）；
* `chooseMatrices`（RepeatMasker:4229）：≤36→35g、≤38→37g、≤40→39g、≤42→41g、
  ≤44→43g、≤46→45g …（每 2% 一档，到 53g 封顶）；
* 最终 `-matrix 20p<GC>.matrix`，文件在 `Matrices/ncbi/nt/20p35g.matrix …
  20p53g.matrix`（共 10 个）；`BLASTMAT` 指向该目录（NCBIBlastSearchEngine.pm:716）。

## 4. 命令与管道设计

### 4.1 命令

`pgr rept rmblast <repeats> <genome> [options]`，输出 runlist JSON
（与 e-kmer / e-align / trf 同构）。

命名理由：**与 `trf` 保持一致**（工具名即子命令名，无 e/s 前缀）；库仍作
第一个位置参数（同 e-kmer / e-align）。RM 参数面（cutoff/word_size/matrix/
frag）与 e-align 的 PGI 参数面完全正交，塞进 `--engine rmblast` 会让
e-align 的 CLI 变脏，故单独命令。

选项（初始）：

| 选项 | 默认 | 说明 |
| :--- | :--- | :--- |
| `--cutoff` | 225 | 复刻 `-cutoff`（直接传给 rmblastn，不打折） |
| `--word-size` | 9 | 或 `--speed slow\|default\|quick\|rush` → 8/9/11/13 |
| `--matrix-gc` | auto | auto=按染色体 GC 计算（=RM 单序列 batch）；`43` 复刻多序列 batch 默认；或显式数值 |
| `--gccalc` | off | 复刻 `-gccalc`（batch 平均 GC，仅配合 `--frag` 有意义） |
| `--frag` | 60000 | 复刻 RM 默认分片（fragmentLen/overlap=2000，SimpleBatcher 算法）；0=整条染色体 |
| `--min-len` / `--fill-fragment` | 0 / 0 | 0=完全忠实 RM 原始 hit；设 >0 与 e-align 可比 |
| `--parallel` | 8 | 每个 rmblastn `-num_threads 4`（复刻 RM 4 核/batch），并行数 ≤8 |
| `--matrix-gc` 的矩阵文件 | 内置 | 仿 lastz.rs 模式内置 20p 矩阵（见 §5.3） |
| `--rmblast-dir` | PATH | makeblastdb/rmblastn 所在目录 |

### 4.2 管道

```
genome.fa(.gz) + repeats.fa(.gz)
  ├─ gz 解压 / 拆染色体（复用 fa split name + name_map 处理点号名）
  ├─ 可选 RM 式分片 60kb/2kb（记录 batch 偏移）
  ├─ makeblastdb -dbtype nucl -in <lib> -out <db>
  ├─ 每染色体/分片：rmblastn <§2 全部参数> -query <chr> -db <db>
  │     -matrix 20p<GC>.matrix  （BLASTMAT=<临时矩阵目录>，内置矩阵，见 §5.3）
  ├─ 解析 18 列 tab：qseqid/qstart/qend（+偏移 → 基因组坐标），sstrand 忽略
  ├─ run_repeat_runlist_pipeline：cover → excise(min-len) → fill(fill-fragment)
  └─ runlist JSON（空结果输出 {}）
```

实现位置：薄壳 `src/cmd_pgr/rept/rmblast.rs`；管道/解析逻辑在
`src/libs/pl/repeat.rs` 新增 `run_rmblast_repeat_pipeline`（与
`run_align_repeat_pipeline` 并列），tab 解析器独立函数 + 单测。

## 5. 一致性决策（复刻 vs 简化）

1. **复刻**：makeblastdb 命令、全部引擎参数、outfmt 列、`-num_alignments`、
   `-mask_level 101`、`-complexity_adjust`、minscore=225（不打折）、
   GC→矩阵映射。
2. **分片默认复刻**（`--frag 60000`、overlap 2000，SimpleBatcher 算法，
   每片段独立 batch + 片段 GC 选矩阵）：这是 RM 默认行为；实测与 RM .out
   对拍后达到 100% 覆盖（§2.8）。`--frag 0` 提供整染色体模式。
3. **矩阵来源（仿 lastz.rs）**：pgr 对 lastz 的得分矩阵处理（
   `src/libs/lastz.rs`）是现成先例——把 UCSC 矩阵内置为 `pub const &str`，
   运行时写临时文件、`Q=<path>` 引用、`NamedTempFile` 保活到进程结束。
   rmblast 照搬该模式：把 10 个 20p 矩阵（`20p35g.matrix` …
   `20p53g.matrix`，文本、各几 KB、OSL-2.1 许可）内置为常量，运行时把选中的
   矩阵写入临时目录（文件名保持 `20p<GC>.matrix`），`std::process::Command`
   `.env("BLASTMAT", <tempdir>)` 后以 `-matrix 20p<GC>.matrix` 调用 rmblastn
   （rmblastn 只按名查矩阵，必须经 BLASTMAT，不能像 lastz 那样直接给路径）；
   临时目录在调用期间保活。**好处**：与 lastz 一致、零外部路径依赖、结果可复现。
4. **不做 ProcessRepeats**：cycleReJoin（碎片连链）、边界精修、edge effect
   移除、family/class、K2P、报告——全部不做；区间侧只做 cover+excise+fill。
5. **TRF 阶段不包含**：`pgr rept trf` 已有；组合用法 = RM 默认
   （散在 + 串联）的完整遮蔽。

预期与 RM 的差异（验证时量化）：无碎片连链 → 同一元件可能多段（合并后影响
小）；无边界精修 → 端点略毛糙；无 edge 修复 → 仅分片模式下有边界假片段。

## 6. 验证

1. 单元测试：18 列 tab 解析（正/负链、batch 偏移）、GC→矩阵映射、
   minscore 折扣、`--cutoff` 覆盖。
2. 集成测试（`tests/cli_rept.rs` 风格）：合成数据（400 bp 重复插入两次）检出；
   `makeblastdb`/`rmblastn` 不在 PATH 时跳过（同 trf 测试）。
3. 对拍（10 株 E. coli + TnCentral，沿用 `scripts/rm-gold-compare.sh` 流程）：
   * rmblast vs RM 原始 IS（`rm-is`）：期望高重合，差异 ≈ 无后处理损失；
   * rmblast（`--min-len 50`）vs RM-strict：验证阈值对齐；
   * 与 e-align PGI(k31) / LASTZ 并列，补齐方法梯队。
   指标：bp、片段数、`runlist statop` 双向覆盖、diff 平均长度。
4. 结论写回 `notes/ecoli-repeats.md`（新增 §2.8）与 `docs/rept.md`。

## 7. 实施步骤

1. `src/libs/pl/repeat.rs`：tab 解析 + GC/矩阵 + minscore 折扣 + 管道
   → 验证：`cargo test`
2. `src/cmd_pgr/rept/rmblast.rs` + `mod.rs` 注册 + `--help`
   → 验证：合成数据 e2e
3. 矩阵内置（方案 B）→ 验证：无 RM 目录依赖下跑通
4. 10 株对拍 → 写回 notes/docs
5. `cargo fmt` + `cargo clippy -- -D warnings` clean

## 8. 风险

* RMBlast 版本：2.14.1 支持 tab 格式与全部参数（已验证）；旧版 2.13 才引入
  tab，更老版本不可用——文档注明最低 2.13。
* 逐字节一致不可达：ProcessRepeats 的后处理（连链/精修）不做，目标定为
  **hit 集一致**（同一搜索参数下的 rmblastn 输出一致），区间层允许差异。
* 性能：rmblastn 比 pgi 慢一个量级；细菌规模无压力，真核需 `--frag` + 并行。
* 矩阵许可：OSL-2.1 允许复制，内置前确认文件头无额外限制。
