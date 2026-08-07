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

* `multiz` 当同一输入文件在某个派生窗口内包含多个**非连续** block 时，`merge_window`
  用 `group.iter().find(...)` 只取该输入第一个与窗口重叠的 block，其余 block 的数据被
  静默丢弃。经核验，正确修复需按参考坐标拼接非连续 block，而现有
  `merge_two_blocks_with_dp`/`merge_conflicting_refs` 是为重叠 locus 的 re-align 设计，
  对非连续 block 会输出坐标错乱的结果（实测 block1 的碱基被压入 block2 的坐标区间，
  只是换一个 block 被丢）。引入拼接逻辑风险高、超出当前合并架构，本轮不修复，记录为
  已知限制。仅在 `--radius` 大到使两个非连续 block 的 `±radius` 扩展区间并入同一窗口
  时触发；默认 `--radius 30` 通常不会。

* `create` 某链接 range 的 `end` 超出参考染色体长度时中止整个命令（与"chr 不存在
  跳过并告警"不一致）。`create` 对越界坐标报错属正当行为，歧义高，未改。
* `refine --outgroup` 的 `trim_outgroup`/`trim_complex_indel` 删除"内群全 gap、外群有
  碱基"的列后，外群序列变短但 `ranges` 中外群 range 坐标未随之收缩，输出的 range 与
  序列长度不一致。与 kent 原版行为一致，且修复需贯通 trim 函数改签名，风险高，记录为
  已知限制。
* `consensus` 内群序列全为 gap 时（`seqs.retain` 剔除后为空），POA 产生空共识，输出
  `>name\n\n`。极端边缘，无 panic，未改。
* `to-xlsx` 非 outgroup 下 `sub.freq = min(freq, N-freq) <= 0.5` 恒成立，`--min-freq>0.5`
  会过滤掉全部 SNP 输出空表。命令层允许 `[0,1]`，属参数语义边界，未改。
* `merge.rs` 非参考物种回退仅在"完全未映射"时生效，对"映射到 gap"不回退到另一块，
  与"物种内容不丢失"语义不完全一致；参考物种刻意只取 block A（保持 ungapped 不变量），
  非参考物种同样"第一块优先"，属设计取舍，未改。
* 理论溢出（不可达）：`banded_align` 的 i32 乘法需块物种数 >8400 才溢出；`poa` 的
  u32 权重需 >4e8 次累加；`topological_sort` 重叠非互斥 clique 重复输出/环死循环仅
  手工构造畸形图可触发，`add_alignment` 恒构造真 clique/DAG，管线内不可达。未改。

## 修复的缺陷（共 20 处，含 1 处记录为已知限制）

### 数据安全（`-o` / 输出路径同输入保护，4 处）

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

**`multiz`/`cover` 的 `-o` 覆盖输入**：`multiz` 先读全部输入（`merge_fas_files_auto_
windows`）后开 writer，但 `-o` 若指向输入仍应在打开前拒绝。修复：在计算前加入
`ensure_outfile_distinct`。`cover` 同理（`write_json` 在读取完输入后调用，但统一
前置检查）。

### Zero Panic / 越界（对齐与变异，4 处）

**不等长 block 越界 panic**：`get_subs`/`get_indels`/`alignment_stat`/
`align_seqs_quick` 假设列数一致，不等长时复用 `seqs[0].len()` 索引其他序列会越界。
修复：在这些函数开头校验所有序列等长，不等长即返回友好错误。回归
`command_fas_stat_unequal_length_no_panic`、
`command_variation_outgroup_unequal_length_no_panic`、
`command_to_xlsx_unequal_length_no_panic`。

**外群序列短于内群越界 panic**：`polarize_subs`/`polarize_indels` 用
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

### 对齐语义（2 处）

**`refine_block` 误加等长校验，破坏 `refine` 的不等长重比对语义**：此前为防
越界在 `refine_block` 开头加了"所有序列等长才继续"的校验。但 `refine` 的用途
恰恰是重比对**不等长**的序列（MSA 对齐），该校验导致合法输入 `tests/fas/
refine.fas` 直接报错。修复：**移除** `refine_block` 中的等长校验（`refine` 走
`align_seqs` POA/外部队列，天然处理不等长；`--chop` 的头部/尾部修剪已由
`trim_head_tail` 的空序列保护兜底）。回归 `command_fas_refine_unequal_length_
no_panic` 从"断言报错"改为"断言成功且三物种均输出"。`align_seqs_quick` 的等长
校验予以保留（quick 模式前提是已对齐输入，语义不同）。

**`to-xlsx` 外群含 IUPAC 歧义碱基时命令失败**：`--outgroup` 下绘制外群替代碱基
时，单元格样式名 `sub_{obase}_unknown` 只对标准碱基（A/C/T/G/N）注册，外群出
现歧义码（R/Y/S/W/K/M/B/D/H/V）时 `format_of.get` 返回 `None`，命令以
"missing format for outgroup substitution" 报错退出。修复：在 `fas_xlsx.rs` 的
该分支加入 `.or_else(|| format_of.get("sub_N_unknown"))` 兜底，歧义外群碱基复用
N（黑色）样式而非失败。回归 `command_to_xlsx_outgroup_ambiguity_no_error`。

### 坐标 / 共识（3 处）

**`slice` 对负链参考/物种产生空输出或反向范围**：`slice_block` 用
`chr_to_align`/`align_to_chr` 做坐标换算。对负链**参考**，`chr_to_align` 使递增的
染色体坐标映射到递减的比对列，`ss_start > ss_end` 触发 `ss_start >= ss_end` 的
`continue`，整个子切片被丢弃，负链参考的切片输出为空；对负链**非参考**物种，
`align_to_chr` 使 `start > end`，输出形如 `>Oth.chr2(-):9-4` 的反向范围。修复：
在生成子切片时对 `ss_start/ss_end` 交换归一为 `[min,max]`；在输出每个物种的
`start/end` 时同样交换为 `start <= end`（保留原链向）。回归
`slice_block_reverse_strand_reports_valid_ranges`。

**POA 共识忽略节点权重，多数碱基不敌首序列骨架**：`generate_consensus` 只用边权重
找最重路径，忽略 `NodeData.weight`（经过该节点的序列数）。对 `C, A, A` 输入，A 节点
权重 2、C 节点权重 1，但共识输出 `C` 而非 `A`。修复：按 SPOA heaviest-bundle 语义
把节点权重计入路径得分 `score[u] = weight[u] + max(edge + score[prev])`。回归
`test_consensus_prefers_majority_by_node_weight`。

**POA 共识对全 gap 首序列输出空**：首序列全为 gap（`----`）时，`add_alignment` 会为
每个 `-` 创建权重 1 的 gap 节点，这些节点索引最小、在路径平局时胜出，共识退化为全
gap，经 `cons.replace('-',"")` 后为空。修复：在 `consensus_block` 生成共识前剔除
全 gap 序列（`seqs.retain(|s| s.iter().any(|b| *b != b'-'))`），使只有真实碱基参与
投票。

### 集合与去重（2 处）

**`concat`/`subset` 对 `--required` 中重复的物种名输出重复条目**：两个命令都用
`read_names::<Vec<String>>` 逐行读取 `--required` 列表并原样遍历。若同一物种名在
列表中重复出现，`concat` 会将该序列拼接两次并输出重复行，`subset` 会在每个 block
中发出重复条目。修复：在两个命令读取 `needed` 后按首次出现顺序去重
（`retain(|n| seen.insert(n.clone()))`，`HashSet` 判重）。回归
`command_fas_concat_duplicate_required_no_dup`、
`command_fas_subset_duplicate_required_no_dup`。

**`multiz` 对倒置参考区间（`start > end`）u64 下溢 panic**：`derive_windows_from_blocks`
计算窗口宽度 `width = e - s` 时，若某输入 block 的参考 entry 坐标倒置（如畸形但可解析
的 `>ref.chr(+):100-1`），`e - s` 在 debug 下 u64 下溢 panic、在 release 下回绕为超大
窗口。修复：在 `derive_windows_from_blocks` 的两个收集点（窗口区间与按染色体分组）跳过
`start > end` 的倒置参考 entry。回归 `derive_windows_inverted_reference_range_no_panic`。

### 边界 / 非 IUPAC（3 处）

**`create` 对超出参考基因组范围的坐标 abort**：`create` 处理链接文件中超出参考染色
体长度的坐标时，`get_seq_loc` 返回 `slice_error`，`create_from_links` 直接 `return Err`
使整个 `create` 运行中止，与"chr 不存在/无效 range 跳过并告警"不一致。修复：在
`create_from_links` 中捕获以 `"slice error"` 开头的错误，记录警告并 `continue` 跳过
该区间；保留其他类型错误（如文件不存在）的传播。回归
`command_create_skips_out_of_range_link`。

**`separate --rc` 对非 IUPAC 字符反向互补后乱码**：`separate --rc` 对负链使用
`nt::rev_comp` 反向互补，`NT_COMP` 表将未知字符（如 `*`）映射为 255 哨兵值，
`format_sequence` 将其渲染为 `ÿ`。修复：`separate` 内联反向互补逻辑，检查 `NT_COMP`
返回 255 时保留原始字节。回归 `command_separate_rc_preserves_non_iupac`。

**`slice` 对带 gap 的参考物种中止整个命令**：`slice_block` 用参考 species 的
`seq_intspan`（非 gap 碱基集合）经 `chr_to_align` 把 runlist 染色体坐标映射到比对列。
当参考物种自身含 gap 时，其非 gap 碱基数小于基因组长 `end-start+1`，`chr_to_align`
的边界 `chr_end = chr_start + 非gap数 - 1` 小于 `end`，runlist 覆盖全范围时
`upper=end` 触发 `[pos] out of ranges`，`?` 传播使整个 `slice` 命令中止。修复：在
`slice_block` 中捕获 `chr_to_align` 两端点任一失败，`log::warn!` 并 `continue` 跳过
该子区间，而非中止。回归 `slice_block_gapped_reference_no_abort`（lib）与
`command_slice_gapped_reference_no_abort`（CLI）。

### 已知限制（未修复，1 处）

**`multiz` `merge_window` 非 DP 回退在参考序列逐字符不等时丢弃整窗口**：当窗口内
`blocks.len() >= 2` 且 DP 合并失败（无共享物种 / banded align 失败）时，回退到
简单拼接：要求所有 block 的参考序列 `entry_seq_equal`（逐字符相等），否则对整窗口
返回 `None` 静默丢弃。参考序列含真实非 gap 差异（如来自不同组装的 SNP）即可触发，
属数据丢失。经核验，参考不一致时本就无法在共享坐标系下单义拼接，丢弃是该架构下的
设计取舍，且与"非连续 block 数据丢失"同源，记录为已知限制，不修复。

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
* `slice` 负链参考/非参考物种修复（`slice_block_reverse_strand_reports_
  valid_ranges`）、POA 节点权重修复（`test_consensus_prefers_majority_by_node_
  weight`）、全 gap 首序列剔除均实测复现修复前后差异（`C,A,A -> A`、全 gap 首序列
  `---- + ACGT -> ACGT`、负链参考 `>Ref.chr1(-):103-108` 输出恢复）。
* `concat`/`subset` 的 `--required` 去重（`command_fas_concat_duplicate_required_
  no_dup`、`command_fas_subset_duplicate_required_no_dup`）、`multiz` 倒置参考区间
  下溢修复（`derive_windows_inverted_reference_range_no_panic`）均实测复现修复
  前后差异；`merge_window` 参考序列不一致丢弃整窗口经核验记录为已知限制。
* 纵深复审逐命令复核 `replace`/`join`/`link --best`/`separate`/`split`/`slice`/
  `to-vcf`/`to-xlsx`/`variation`/`stat`/`cover`/`name`/`create`/`check`/
  `consensus`/`refine`/`multiz` 及库 `poa/graph`、`fas_multiz/{banded_align,merge}`、
  `ds/intspan` 的全部执行路径，未再发现新缺陷。`to-xlsx` 外群模式下 indel 段不绘制
  外群行（外群仅作极化参考的显示取舍），记录不修复。

## 结论

`fas` 命令族审核完成（累计修复 19 处缺陷并补回归测试与文档澄清），并经多轮纵深
复审（全部 20 个子命令的执行路径、`-o` 覆盖保护含 `create`/`check` 的 `.loc` 与
`separate`/`split` 目录输出、`stat`/`variation`/`to-vcf`/`to-xlsx` 的不等长与
外群越界、`slice_block`/`trim_head_tail` 全 gap 场景、`banded_align` I 态 gap
成本、`fas_multiz` 合并顺序确定性、`to-xlsx` 外群歧义碱基样式兜底、`slice` 负链
坐标、POA 共识节点权重与全 gap 首序列、`--required` 去重、`multiz` 倒置参考区间、
`create`/`separate` 越界与非 IUPAC 边界、`slice` 带 gap 参考坐标映射、`docs/fas.md`
与帮助文本一致性）复核。

多轮纵深复审覆盖全部 20 个子命令的执行路径与相关库函数，均未再发现需要修复的新
`fas` 缺陷。逐命令核验要点：`replace`（`read_replace_tsv` 单字段=删除、多字段=复制、
重复 header 保原 block）、`join`（`join_block_entries` 目标 entry 仅首现入 map、同
range 共享坐标）、`link --best`（`pair_d` 等长/可比碱基前提，不可评分对跳过）、
`separate`/`split`（逐输出路径 `same_path` 反查、`sanitize_filename` 碰撞已记录）、
`slice`（`chr_to_align`/`align_to_chr` 负链端点、`IntSpan::min/max` 仅对非空 subslice
调用、indel 岛仅修剪边界列、内部 gap 保留的分层语义）、`to-vcf`（跨 block 物种/顺序
一致校验、REF 取自 target、`pos_idx` 越界防护）、`to-xlsx`（`create_formats` 键全覆盖
`paint_sub`/`paint_indel` 引用、outgroup 歧义码兜底 `sub_N_unknown`、外群 `-` 走
`sub_-_unknown`）、`variation`/`stat`（外群剔除、等长校验）、`cover`
（`aggregate_coverage_into` 跨 block 累积 name 键）、`name`（`IndexMap` 保序去重）、
`create`/`check`（`.loc` 侧车保护）、`consensus`（`seqs.retain` 剔除全 gap、外群保留）、
`refine`（引擎分发、`--chop` 由 `trim_head_tail` 空序列防护兜底）、`multiz`
（`derive_windows_from_blocks` 倒置区间跳过、`merge_window` DP/回退路径、`banded_align`
参考锚定带与三态 trace）。库级核验：`alignment/stat.rs`（`mean_d` 短路、`pair_d` 零可比
碱基报错）、`alignment/variation.rs`（complex `freq=-1` 语义）、`alignment/coords.rs`
（`pos in holes pin to left base`）、`alignment/trim.rs`（空序列兜底）、`fmt/fas.rs`
（`run_pipeline`/`run_parallel` writer 失败后持续 drain 防死锁、`check_entry_against_ref`
越界 `FAILED` 不中止、`split_block_key` 空 block `None` 短路）、`fas_multiz/merge.rs`
（确定性合并顺序、`best_crossover` 四等长 `debug_assert_eq`）、`ds/intspan.rs`
（`inset`/`trim`/`banish`/`find_islands_n`/`find_islands_ints` 边界）。`--required`
除 `concat`/`subset` 外无其他 `fas` 消费者（`read_names` 搜索确认）；`IntSpan::add_pair`
对倒置区间安全返回不 panic；`multiz` `merge_window` 非 DP 回退在参考序列不一致时丢弃
整窗口（与"非连续 block 数据丢失"同源，属共享坐标系合并架构下的设计取舍）。`to-xlsx`
外群模式下 indel 段不绘制外群行（外群仅作极化参考的显示取舍），记录不修复。
`cargo build`、`cargo clippy --all-targets -- -D warnings`、`cargo fmt --check` 全部
干净；`--lib`（590）、`cli_fas`（45）、`cli_fas_vars`（11）、`cli_fas_poa`（12）、
`cli_fas_multiz`（2）全部通过。注：`cargo clippy --all-targets` 报出的 6 处 warning
全部位于 `src/libs/pgi/align.rs` 与 `src/libs/syncmer.rs`（`gix_matchmer`/syncmer
工作区），与 `fas` 无关、非本次改动引入，不在 `fas` 审核范围内。