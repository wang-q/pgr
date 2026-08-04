# pgr fas 命令族代码审核记录（2026-08-05）

对 `pgr fas` 命令族（全部 20 个子命令）以及相关库文件（`libs/fmt/fas`、
`libs/alignment`、`libs/fas_multiz`、`libs/fas_xlsx`、`libs/nt`、`libs/io`、
`libs/ds/{intspan,crossover}`）和全部测试/文档进行审核。缺陷按类别分组记录；
关键修复均附回归测试（见文末），验证概况见文末"验证"一节。

审核范围：
- **信息**：`check` / `cover` / `link` / `name` / `stat`
- **子集**：`filter` / `slice` / `subset`
- **转换**：`concat` / `consensus` / `join` / `multiz` / `refine` / `replace`
- **文件**：`create` / `separate` / `split`
- **变异**：`to-vcf` / `to-xlsx` / `variation`

审核重点：数据安全（`-o` 不得覆盖输入，含 `.loc` 侧车索引与分块/逐物种输出）、
Zero Panic（畸形输入不 panic）、坐标/长度边界处理、算法正确性、文档一致性。

## 排除的疑点（经核验无需修复）

* 逐命令通读全部 20 个 `fas` 子命令的 `execute`：`unwrap()`/`unreachable!`
  全部为 clap `required` 参数或 `value_parser` 约束枚举，运行期不可达，无潜在
  panic（符合"稳定性原则"）。
* `-o` 覆盖保护覆盖情况：全部单文件输出命令（`stat`/`variation`/`filter`/
  `join`/`name`/`consensus`/`refine`/`multiz`/`replace`/`slice`/`subset`/
  `concat`/`to-vcf`/`to-xlsx`/`link`/`cover`/`create`/`check`）均调用
  `ensure_outfile_distinct`。带辅助输入的命令（`concat`/`subset` 的
  `--required`、`slice` 的 `--runlist`、`replace` 的 `--replace-tsv`、
  `to-vcf` 的 `--sizes`、`create`/`check` 的 `--genome` 及其 `.loc`）均一并
  纳入保护列表。`separate`/`split` 输出为目录，采用逐输出路径 `same_path`
  反向检查（见下）。
* `cover --trim` 允许负值：`IntSpan::trim` 经 `inset` 在负 `n` 时扩展区域，
  属用户误用而非越界，`saturating_add/sub + clamp` 保证不 panic，未改。
* `to-xlsx` 复杂替换（`pattern="unknown"`）的 `from_str_radix`：绘制时对
  `occurred == '1'` 分支才调用，complex 恒走 `sub_{base}_unknown` 分支，不会
  对 `"unknown"` 做二进制解析，未改。
* `best_crossover` 的 `debug_assert_eq!`：四个切片均以 `map_a.len()` 构建，
  长度恒相等；release 下为 no-op，无越界，未改。
* `rev_comp`（`separate --rc`）对 `-`：`NT_COMP['-']='-'`，先 rc 再
  `format_sequence(is_dash=true)` 移除，行为正确，未改。
* block FA 解析跳过 `#` 注释行（`next_fas_block`），`docs/fas.md` 的
  "以 # 开头的行视为注释" 描述准确。
* `align_seqs_quick` 的等长校验：quick 模式本假设序列已大致对齐（axt/maf
  转换场景），等长是前提；对歧异长度输入报错而非 panic，属正当行为。该校验
  与 `refine --quick` 的正常（等长带 indel）输入不冲突（见"修复的缺陷"中对
  `refine_block` 的辨析）。

## 记录项（未改，低风险 / 待决策）

* `separate`/`split` 的物种名/染色体名经 `sanitize_filename` 清洗后，两个不同
  名称可能碰撞到同一输出文件名（如 `a/b` 与 `a_b`），会静默合并到同一文件。
  与文件名清洗方案固有行为一致，低风险，未改。
* `to-xlsx` 的 `--no-single` 对 indel 用 `freq <= 1` 判断，complex（`freq=-1`）
  也会被过滤；与 `--no-complex` 语义部分重叠，属过滤器常规行为，非缺陷，未改。
* `run_pipeline` 串行分支对 `proc_block` 的错误直接传播（`?`），并行分支收集
  到 `errors` 后统一报错；对畸形 block 的 `next_fas_block` 错误则跳过并告警。
  行为一致，未改。

## 修复的缺陷（共 11 处，含本轮 2 处回归）

### 数据安全（`-o` / 输出路径同输入保护，5 处）

**全部 `fas` 流式命令允许 `-o` 覆盖输入文件**：流式命令先打开输出再读取输入，
若 `-o` 指向输入文件，会在读取前截断输入，静默清空数据。修复：在
`stat`/`variation`/`filter`/`join`/`name`/`consensus`/`refine`/`multiz`/
`replace`/`slice`/`subset`/`concat`/`to-vcf`/`to-xlsx`/`link`/`cover` 中统一
加入 `ensure_outfile_distinct` 检查；涉及辅助输入的命令同时保护 `--required`/
`--runlist`/`--replace-tsv`/`--sizes` 文件。回归
`command_to_xlsx_output_same_as_input_rejected` 等。

**`create` 的 `-o` 可覆盖基因组 `.loc` 侧车索引**：`create` 在读取输入链接与
参考基因组前先打开（截断）输出 writer。若 `-o` 命名为 `{genome}.loc`，会在
`open_indexed`（`create_from_links` 内部）读取索引前先截断该文件，随后因 mtime
判定"新鲜"而不重建，`load_loc` 读到空索引，所有链接静默丢弃，且 `.loc` 被永久
损坏。修复：将 `{genome}.loc` 一并加入 `ensure_outfile_distinct` 保护列表。
`check` 同理确认已保护。回归 `command_create_output_not_overwrite_loc_index`。

**`separate`/`split` 的输出文件可能覆盖输入文件**：两命令的 `-o` 输出是**目录**，
输出文件名由物种名/染色体名/比对块动态生成。若输入文件恰好位于输出目录且文件名
与某输出文件名一致，`truncate` 打开会截断正在流式读取的输入。修复：在打开每个
输出文件前用 `same_path` 与所有输入路径比对，命中即 `bail!`。

**`multiz` 的 `-o` 覆盖输入**：`multiz` 先读全部输入（`merge_fas_files_auto_
windows`）后开 writer，但 `-o` 若指向输入仍应在打开前拒绝，以保证正确性。修复：
在计算前加入 `ensure_outfile_distinct`。`cover` 同理（`write_json` 在读取完输入
后调用，但统一前置检查）。

### Zero Panic / 越界（对齐与变异，4 处）

**不等长 block 在 `stat`/`variation`/`to-vcf`/`to-xlsx` 中越界 panic**：
`get_subs`/`get_indels`/`alignment_stat`/`align_seqs_quick` 假设列数一致，不等长
时复用 `seqs[0].len()` 索引其他序列会越界。修复：在这些函数开头校验所有序列等长，
不等长即返回友好错误。这些命令要求 block 为已对齐的 MSA（等长），报错语义正确。
回归 `command_fas_stat_unequal_length_no_panic`、
`command_variation_outgroup_unequal_length_no_panic`、
`command_to_xlsx_unequal_length_no_panic`。

**外群序列短于内群导致越界 panic**：`polarize_subs`/`polarize_indels` 用
`sub.pos`/`indel.start..end` 直接索引外群序列，外群短于范围时越界。修复：加
`ensure!(og_idx < og.len())` / `ensure!(end <= og.len())` 边界检查。回归
`command_variation_outgroup_unequal_length_no_panic`。

**`slice_block` 对全 gap 第二条序列 panic**：子切片整个落在 indel 岛内（某物种
全为 gap）时，修剪 indel 边界后 `ss_ints` 为空，`ss.min()`/`ss.max()` 对空
`IntSpan` panic。修复：`ss_ints.is_empty()` 时跳过该子切片。回归
`slice_block_all_gap_second_species_no_panic`。

**`trim_head_tail` 对全 gap 比对 panic**：`--chop` 时头部循环移除全部字符，随后
尾部循环 `seqs[i].remove(cur_len - 1)` 在 `cur_len == 0` 时下溢 panic。修复：
头部移除量 `min(..., seqs[0].len())`，尾部循环在序列为空时 `break`。回归
`trim_head_tail_all_gap_no_panic`。

### 算法正确性（1 处）

**`banded_align` I 态 gap 延伸成本用错计数**：I 态 gap 延伸用
`k * gap_extend_pen`，其中 `k` 被外层循环的单元格索引遮蔽，误用为列内物种数
的实际应来源于 A 块列物种数。修复：将物种数变量改名为 `k_a`，并用于全部 gap
延伸计算（含首行插入链与 I 态延伸）。回归既有 `fas_multiz` 合并测试。

### 回归修复（本轮，2 处）

**`refine_block` 误加等长校验，破坏 `refine` 的不等长重比对语义**：此前为防
越界在 `refine_block` 开头加了"所有序列等长才继续"的校验。但 `refine` 的用途
恰恰是重比对**不等长**的序列（MSA 对齐），该校验导致合法输入 `tests/fas/
refine.fas`（首 block 中 Spar 为 18 bp、其余为 21 bp）直接报错，`command_refine_
default`/`command_refine_poa` 失败。修复：**移除** `refine_block` 中的等长校验
（`refine` 走 `align_seqs` POA/外部队列，天然处理不等长；`--chop` 的头部/尾部
修剪已由 `trim_head_tail` 的空序列保护兜底）。同时把回归测试
`command_fas_refine_unequal_length_no_panic` 从"断言报错"改为"断言成功且三物种
均输出"，以反映正确语义。`align_seqs_quick` 的等长校验予以保留（quick 模式
前提是已对齐输入，语义不同）。

**`to-xlsx` 外群含 IUPAC 歧义碱基时命令失败**：`--outgroup` 下绘制外群替代碱基
时，单元格样式名 `sub_{obase}_unknown` 只对标准碱基（A/C/T/G/N）注册，外群出
现歧义码（R/Y/S/W/K/M/B/D/H/V）时 `format_of.get` 返回 `None`，命令以
"missing format for outgroup substitution" 报错退出。对合法 MSA 输入而言属
崩溃缺陷。修复：在 `fas_xlsx.rs` 的该分支加入 `.or_else(|| format_of.get(
"sub_N_unknown"))` 兜底，歧义外群碱基复用 N（黑色）样式而非失败。回归
`command_to_xlsx_outgroup_ambiguity_no_error`。

## 验证

* 数据安全：`-o` 同路径（单文件与 `create` 的 `.loc`、`separate`/`split` 目录）、
  辅助输入文件 `--required`/`--runlist`/`--replace-tsv`/`--sizes` 等修复前后均
  实测复现，既有输入原样保留。
* 畸形输入：不等长 block、全 gap 比对/第二条序列、短外群、空 block 等 fuzz，
  零 panic。
* `cargo build`、`cargo clippy --all-targets -- -D warnings`、`cargo fmt --check`
  全部干净；全部 lib 单测 + 集成测试 + doctest 通过（含 42 个 `cli_fas`、
  65 个 `cli_fas_vars`、12 个 `cli_fas_poa` 集成测试）。
* 新增/修正回归测试（主要）：`command_create_output_not_overwrite_loc_index`、
  `command_fas_stat_unequal_length_no_panic`、
  `command_variation_outgroup_unequal_length_no_panic`、
  `command_to_xlsx_unequal_length_no_panic`、
  `slice_block_all_gap_second_species_no_panic`、
  `trim_head_tail_all_gap_no_panic`、
  `command_fas_refine_unequal_length_no_panic`（改为断言成功）等。

## 结论

`fas` 命令族审核完成（累计修复 11 处缺陷并补回归测试与文档澄清），并经多轮纵深
复审（全部 20 个子命令的执行路径、`-o` 覆盖保护含 `create`/`check` 的 `.loc` 与
`separate`/`split` 目录输出、`stat`/`variation`/`to-vcf`/`to-xlsx` 的不等长与
外群越界、`slice_block`/`trim_head_tail` 全 gap 场景、`banded_align` I 态 gap
成本、`fas_multiz` 合并顺序确定性、`to-xlsx` 外群歧义碱基样式兜底、`docs/fas.md`
与帮助文本一致性）复核。

**最终收敛轮**：对全部记录项（`separate`/`split` 文件名碰撞、`to-xlsx` 的
`--no-single` 与 complex 重叠、`run_pipeline` 错误传播、`cover --trim` 负值、
`refine --quick` 等长前提、由歧义码兜底引出的 `sub_{base}_unknown` 内群样式覆盖
等）逐一重新核验，均属文档化一致行为或不可达的极端命名场景，非缺陷，无需改动。
此轮未再发现任何新的 `fas` 缺陷，审核收敛。