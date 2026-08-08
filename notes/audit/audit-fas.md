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
* **POA 回溯 `curr_j - 1` 下溢疑点**（`poa/align.rs` E 态）：E 态分支访问
  `e[curr_i][curr_j-1]` 前无显式 `curr_j > 0` 守卫。经核验，进入 E 态时恒有
  `curr_j >= 1`：所有模式（global/local/semi）下 `e[*][0]` 均为 `neg_inf`，而任何
  转移进 E 态（M 态 `e[u][j-1]==target`、F 态 `e[u][j]==target`、E 态延伸
  `e[i][j-1]+extend==target`）在 `curr_j==1→0` 处都要求 `e[*][0]` 等于一个真实分值，
  恒不可能；E 态 else 分支无条件 `curr_j -= 1` 只会转移到 F 态（不越界）。用 3000 余
  个随机图 × 3 模式 × 随机序列（含经 `add_alignment` 构建的真管线拓扑）fuzz，无一
  panic。判定为不可达，不加守卫（避免推测性防御代码）。
* **`write_vcf_block` 的 `pos_idx`**：`checked_sub(1)` 与 `pos_idx >= seqs[0].len()`
  双重校验，越界返回友好错误，无 panic。
* **`merge.rs` `cut`/`map_a`/`map_b` 越界疑点**：`banded_align_refs` 在返回前强制
  `map_a.len() == map_b.len() && !map_a.is_empty()`（banded_align.rs L425），且裁剪后
  同区间切片保持等长；`col()` 用 `seq.get(idx)` 带边界，无越界。
* **`slice.rs` `ss_start`/`seq_len_i32`**：`i32::try_from` + `ss_start < 1 || ss_end >
  seq_len_i32` 边界检查，`start_idx = ss_start-1`、`end_idx = ss_end` 均在界内。
* **`variation.rs` `uniq_indel_seqs` `min_by_key`**：`indel_seqs` 恒非空（至少 1 条
  序列），`unique` 后非空，`.ok_or_else` 为防御性代码，不 panic。
* **`align_seqs_quick` 负 pad / 下溢**：`pad`/`fill` 经 `i32::try_from` 校验；`pad >
  align_len` 时 `align_len - pad` 为负，`add_pair(负, align_len)` 非倒置区间被接受，但
  随后与 `IntSpan::from_pair(1, align_len)` 相交，负坐标被裁掉，最终 `lower>=1`、切片
  越界前有 `intersect` 兜底，无 panic。
* **`consensus_block` 保留空序列 / `refine_block` 外群校验**：`seqs.retain` 剔除非全
  gap 序列，空序列（len 0）在 POA 中不产生节点；`refine` 外群要求 `n >= 3` 为显式
  校验，属正当行为。
* **`Range::from_str` 大坐标溢出**：生产路径 `decode`→`match_at`→`tail_match`→
  `parse_i32` 对溢出数字串返回 `None`（不匹配，不 panic）；带 `.unwrap()` 的
  `from_str_regex` 仅用于单测对照（坐标受限），不在发布二进制可达路径。

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
* **`multiz` `merge_window` 非 DP 回退在参考序列逐字符不等时丢弃整窗口**：当窗口内
  `blocks.len() >= 2` 且 DP 合并失败（无共享物种 / banded align 失败）时，回退到
  简单拼接：要求所有 block 的参考序列 `entry_seq_equal`（逐字符相等），否则对整窗口
  返回 `None` 静默丢弃。参考序列含真实非 gap 差异（如来自不同组装的 SNP）即可触发，
  属数据丢失。经核验，参考不一致时本就无法在共享坐标系下单义拼接，丢弃是该架构下的
  设计取舍，且与"非连续 block 数据丢失"同源，记录为已知限制，不修复。
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

## 修复的缺陷（共 32 处）

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

### Zero Panic / 越界（9 处）

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

**外部对齐器返回越界序列 id 时 panic**：`align_seqs`（`refine --engine` 非 builtin 走
外部 MSA 程序）解析外部对齐器输出时，用记录头 `>N` 的 `N` 直接索引 `out_seqs`。若外部
对齐器返回额外/重编号的记录（id 超出输入序列索引范围），`out_seqs[idx]` 越界 panic。
修复：在写入前校验 `idx >= out_seqs.len()` 即返回友好错误，拒绝写入错标输出。

**`multiz` 窗口边界 off-by-one**：`ref_overlaps_window` 原用 `start < window.end &&
end > window.start` 判断 block 参考区间与窗口是否重叠。两区间均为 1-based 且含端点，
一个恰好落在窗口边界上的单碱基 block 会被判为不重叠而丢弃。修复：改为含端点的
`start <= window.end && end >= window.start`。回归
`ref_overlaps_window_includes_boundary_positions`。

**`IntSpan` runlist 解析空 token 静默注入坐标 0**：`IntSpan::try_from` 解析 runlist 时，
空 token（前导逗号 `,1`、连续逗号 `1,,3`）会退化为隐式 `(0,0)` 对，静默把坐标 0 注入
集合。`slice`/`cover` 等命令的 `--runlist` 走此解析。修复：遇到逗号开头的空 run 即返回
"Number format error: empty run" 错误；尾随逗号（`1,`）仍作无害换行终止符。回归
`try_from_rejects_empty_run_tokens`。

**`align_to_chr` 对空序列 intspan panic**：`alignment/coords.rs` 的 `align_to_chr`
在 `ints.is_empty()` 时会调用 `ints.min()`/`ints.max()`，二者对空 `IntSpan` 直接
`panic!`。经核验，`fas` 现有全部调用方（`slice_block`、`write_variations`、
`write_vcf_block`）当前均不可达该路径：`write_variations`/`write_vcf_block` 的
`t_ints_seq` 仅在目标全 gap 时为空，而 `get_subs` 只在"所有序列该列均为真实碱基
（`NT_VAL<=3`）且存在差异"时产出 substitution，目标全 gap 时无任何 substitution 需
换算坐标；`slice_block` 中任一物种全 gap 会使 `indel_ints` 覆盖全列，子切片被
`ss_ints.is_empty()` 完全剔除。但 `align_to_chr` 是公开库函数，与其对称的
`chr_to_align`（对空参考自然返回错误）行为不一致，属潜在 panic 隐患。修复：在函数
入口对空 intspan 返回友好错误。回归 `align_to_chr_empty_intspan_errors`。

**`multiz` 对倒置参考区间（`start > end`）u64 下溢 panic**：`derive_windows_from_blocks`
计算窗口宽度 `width = e - s` 时，若某输入 block 的参考 entry 坐标倒置（如畸形但可解析
的 `>ref.chr(+):100-1`），`e - s` 在 debug 下 u64 下溢 panic、在 release 下回绕为超大
窗口。修复：在 `derive_windows_from_blocks` 的两个收集点（窗口区间与按染色体分组）跳过
`start > end` 的倒置参考 entry。回归 `derive_windows_inverted_reference_range_no_panic`。

### 算术溢出（2 处）

**`to-xlsx` 序列数 >32 时颜色索引溢出**：`paint_indel` 用 `u32::from_str_radix(&occurred,2)`
、`paint_sub` 用 `u32::from_str_radix(&pattern,2)` 把与内群序列等长的二进制串解析为 u32 以
取色索引。当 block 内群序列数 >32 时二进制串超过 u32 位宽，`from_str_radix` 返回溢出错误，
命令失败。修复：改为按字节折叠模 `color_loop` 计算索引（`(acc*2 + (b=='1')) % color_loop`），
对短串结果一致、对长串无溢出。

**`to-xlsx` `--wrap` 大值导致 u16 溢出**：`export_to_xlsx` 用 `opt.wrap + 3` 遍历设置列宽，
当 `--wrap` 取接近 u16 上限（如 65535）时 `wrap + 3` 溢出回绕，列宽设置错乱。修复：改为
`saturating_add(3)`；`paint_indel` 的 `col_cursor + col_taken` 同样改 saturating，`paint_sub`
的 `col_cursor += consumed` 亦防溢出。

### 算法正确性（3 处）

**`banded_align` I 态 gap 延伸成本用错计数**：I 态 gap 延伸用
`k * gap_extend_pen`，其中 `k` 被外层循环的单元格索引遮蔽，误用为列内物种数
的实际应来源于 A 块列物种数。修复：将物种数变量改名为 `k_a`，并用于全部 gap
延伸计算（含首行插入链与 I 态延伸）。回归既有 `fas_multiz` 合并测试。

**POA 共识忽略节点权重，多数碱基不敌首序列骨架**：`generate_consensus` 只用边权重
找最重路径，忽略 `NodeData.weight`（经过该节点的序列数）。对 `C, A, A` 输入，A 节点
权重 2、C 节点权重 1，但共识输出 `C` 而非 `A`。修复：按 SPOA heaviest-bundle 语义
把节点权重计入路径得分 `score[u] = weight[u] + max(edge + score[prev])`。回归
`test_consensus_prefers_majority_by_node_weight`。

**`to-xlsx` 列游标推进错误导致单元格重叠/数据损坏**：`export_to_xlsx` 每绘制一个变异后
固定 `opt.col_cursor += 1`。但 `paint_indel` 的多碱基 indel 实际占用 `indel.length.min(3)`
列（合并单元格），固定 +1 会让下一个变异绘制在 indel 已占用的列上，覆盖其扩展列。修复：
`paint_sub`/`paint_indel` 返回实际占用列数（`paint_indel` 返回 `col_taken`），游标按其
推进。

### 对齐语义（1 处）

**`refine_block` 误加等长校验，破坏 `refine` 的不等长重比对语义**：此前为防
越界在 `refine_block` 开头加了"所有序列等长才继续"的校验。但 `refine` 的用途
恰恰是重比对**不等长**的序列（MSA 对齐），该校验导致合法输入 `tests/fas/
refine.fas` 直接报错。修复：**移除** `refine_block` 中的等长校验（`refine` 走
`align_seqs` POA/外部队列，天然处理不等长；`--chop` 的头部/尾部修剪已由
`trim_head_tail` 的空序列保护兜底）。回归 `command_fas_refine_unequal_length_
no_panic` 从"断言报错"改为"断言成功且三物种均输出"。`align_seqs_quick` 的等长
校验予以保留（quick 模式前提是已对齐输入，语义不同）。

### 坐标 / 共识（2 处）

**`slice` 对负链参考/物种产生空输出或反向范围**：`slice_block` 用
`chr_to_align`/`align_to_chr` 做坐标换算。对负链**参考**，`chr_to_align` 使递增的
染色体坐标映射到递减的比对列，`ss_start > ss_end` 触发 `ss_start >= ss_end` 的
`continue`，整个子切片被丢弃，负链参考的切片输出为空；对负链**非参考**物种，
`align_to_chr` 使 `start > end`，输出形如 `>Oth.chr2(-):9-4` 的反向范围。修复：
在生成子切片时对 `ss_start/ss_end` 交换归一为 `[min,max]`；在输出每个物种的
`start/end` 时同样交换为 `start <= end`（保留原链向）。回归
`slice_block_reverse_strand_reports_valid_ranges`。

**POA 共识对全 gap 首序列输出空**：首序列全为 gap（`----`）时，`add_alignment` 会为
每个 `-` 创建权重 1 的 gap 节点，这些节点索引最小、在路径平局时胜出，共识退化为全
gap，经 `cons.replace('-',"")` 后为空。修复：在 `consensus_block` 生成共识前剔除
全 gap 序列（`seqs.retain(|s| s.iter().any(|b| *b != b'-'))`），使只有真实碱基参与
投票。

### 集合与去重（1 处）

**`concat`/`subset` 对 `--required` 中重复的物种名输出重复条目**：两个命令都用
`read_names::<Vec<String>>` 逐行读取 `--required` 列表并原样遍历。若同一物种名在
列表中重复出现，`concat` 会将该序列拼接两次并输出重复行，`subset` 会在每个 block
中发出重复条目。修复：在两个命令读取 `needed` 后按首次出现顺序去重
（`retain(|n| seen.insert(n.clone()))`，`HashSet` 判重）。回归
`command_fas_concat_duplicate_required_no_dup`、
`command_fas_subset_duplicate_required_no_dup`。

### 多序列合并（multiz，2 处）

**`multiz` 渐进合并失败时静默丢弃物种**：`merge_blocks_with_dp` 对 `merge_two_blocks_with_dp`
返回 `None` 的中间块直接丢弃（不保留其物种），与"保留各输入物种并集"的承诺不符。修复：改为
`match` 捕获 `None` 分支，`log::warn!` 提示被丢弃的 block 及其物种。

**`multiz` 覆盖率过滤误删单碱基参考块**：`derive_windows_from_blocks` 原有的 `DupeTree`
覆盖率过滤中，`DupeTree::add` 忽略 `start == end` 的零宽单碱基区间，导致单碱基参考 block
派生的窗口被过滤掉。修复：移除该覆盖率过滤（每个窗口本就从至少一个输入的参考区间派生出，
覆盖恒成立）。

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

### to-xlsx 显示 / 样式（1 处）

**`to-xlsx` 外群含 IUPAC 歧义碱基时命令失败**：`--outgroup` 下绘制外群替代碱基
时，单元格样式名 `sub_{obase}_unknown` 只对标准碱基（A/C/T/G/N）注册，外群出
现歧义码（R/Y/S/W/K/M/B/D/H/V）时 `format_of.get` 返回 `None`，命令以
"missing format for outgroup substitution" 报错退出。修复：在 `fas_xlsx.rs` 的
该分支加入 `.or_else(|| format_of.get("sub_N_unknown"))` 兜底，歧义外群碱基复用
N（黑色）样式而非失败。回归 `command_to_xlsx_outgroup_ambiguity_no_error`。

### 文档一致性（4 处）

**`consensus` 输出格式描述与实现不符**：`docs/fas.md` 原称"每个 block 的首条序列变为一致性
序列，其余序列保留"，但 `consensus_block` 实际用共识序列替换所有内群序列，仅 `--outgroup`
时保留末条外群。已更新为准确描述。

**`calling 变异` 措辞**：`docs/fas.md` 变异组描述中 `从比对中 calling 变异` 中英混杂，改为
`从比对中检测变异（calling variants）`。

**`collect_subs`/`collect_indels` 外群语义描述误导**：`variation.rs` 原 doc 称"`seqs`
最后一条序列被当作外群"。实际约定是：调用方需把同一序列既作为 `seqs` 末元素传入、又作为
`outgroup` 参数传入，末元素被排除出内群。已更新 doc 澄清该约定。

**`trim_head_tail` doc 描述与实现不符**：`trim.rs` 原 doc 称"返回每序列被删除碱基对应的
(head, tail) 元组向量"，实际函数就地修改序列与其 range 并返回 `String`（complex 区域）。
已更新 doc。

## 验证

* 数据安全：`-o` 同路径（单文件与 `create` 的 `.loc`、`separate`/`split` 目录）、
  辅助输入文件 `--required`/`--runlist`/`--replace-tsv`/`--sizes` 等修复前后均
  实测复现，既有输入原样保留。
* 畸形输入：不等长 block、全 gap 比对/第二条序列、短外群、空 block 等 fuzz，
  零 panic。
* `cargo build`、`cargo clippy --all-targets -- -D warnings`、`cargo fmt --check`
  对 `fas` 相关文件全部干净（`pbit` 工作区的 fmt diff 与 `pgi/align.rs`/
  `syncmer.rs` 的 clippy warning 为并行任务独立工作区/既有问题，非 `fas`、非本次
  改动引入）；全部 lib 单测 + 集成测试 + doctest 通过（`--lib` 36 个 fas、
  `cli_fas` 48、`cli_fas_vars` 15、`cli_fas_poa` 12、`cli_fas_multiz` 2）。
* 新增/修正回归测试（主要）：`command_create_output_not_overwrite_loc_index`、
  `command_fas_stat_unequal_length_no_panic`、
  `command_variation_outgroup_unequal_length_no_panic`、
  `command_to_xlsx_unequal_length_no_panic`、
  `slice_block_all_gap_second_species_no_panic`、
  `trim_head_tail_all_gap_no_panic`、
  `command_fas_refine_unequal_length_no_panic`（改为断言成功）、
  `command_to_xlsx_many_sequences_no_overflow`、
  `command_to_xlsx_outgroup_ambiguity_no_error`、
  `ref_overlaps_window_includes_boundary_positions`、
  `try_from_rejects_empty_run_tokens`、`align_to_chr_empty_intspan_errors` 等。
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

## 2026-08-09 Review Fixes（第七轮）

第七轮以全新视角复核全部 20 个子命令与相关库（`fmt/fas.rs`、`alignment/*`、
`fas_multiz/*`、`fas_xlsx.rs`、`poa/*`、`ds/intspan.rs`），并逐项比对
`docs/fas.md` 与 CLI 帮助文本。本轮修复 3 处缺陷：

### 数据完整性（subset）

**`subset` 块内重复物种名静默丢序列**：`subset` 用 `HashMap<name, &FasEntry>` 为每个
block 建名→序列索引，`collect()` 对同名重复条目保留**最后一次**，另一个重复序列被静默
丢弃；而 `concat`（用 `.position()` 取首个）与 `replace`（重复 header 告警）对同类畸形
输入均有明确处理，行为不一致。修复：改为 `entry.entry(name).or_insert(e)` **首现优先**
，与 `concat` 的首匹配语义对齐，不再静默丢序列。回归
`command_fas_subset_duplicate_species_in_block_keeps_first`。

### to-xlsx 版式（paint_indel 换行边界）

**`paint_indel` 换行 off-by-one 与小节双重计数**：`paint_indel` 内部换行判断
`col_cursor + col_taken > wrap` 要求 indel 结束列 ≤ `wrap-1`，使恰好结束在第 `wrap`
列的多碱基 indel 被无谓换入新小节（而 `paint_sub` 可占用第 `wrap` 列，语义不一致）；
且内部换行后调用方又无条件执行 `col_cursor + consumed > wrap` 的后置换行，当
`wrap ≤ col_taken`（即 `--wrap 1/2/3`）时 `sec_cursor` 被连续递增两次，留下空小节并
把名称写进多余小节。修复：将换行预检上移到调用方统一处理（`col_cursor + width >
wrap + 1` 才换行，允许 `width` 列恰结束于第 `wrap` 列），移除 `paint_indel` 内部换行，
消除双重计数。回归 `command_to_xlsx_indel_fits_in_wrapped_section`。

### 多序列合并（multiz merge_conflicting_refs 数据丢失）

**`merge_conflicting_refs` 对仅存在于单块的物种静默截断**：两输入 block 参考序列真实
不一致（`ungapped_equal` 为假）走 `merge_conflicting_refs` 在最佳 crossover `cut` 处拼接
时，原实现对**每个**物种都按 `pos < cut` 用 `map_a`+`group[0]`、`pos >= cut` 用
`map_b`+`group[1]`。对仅存在于 block A 的物种（`group[1]` 为 `None`），其 `pos >= cut`
半段全部塌成 `-`；对仅存在于 block B 的物种同理左半段全 `-`。即单块物种的一半序列被
静默丢弃。对比非冲突路径 `merge_two_blocks_with_dp`（两块回退 `group[0]`→`group[1]`，
单块物种经 `map_a`/`map_b` 全程携带）行为不一致。修复：按物种**在哪些块中存在**选择
映射——共享物种（含参考）保持 `cut` 拼接（左随 A、右随 B，与拼接后的参考一致），仅
存在于 A 的物种用 `map_a` 全程、仅存在于 B 的物种用 `map_b` 全程，不再截断。回归
`merge_window_conflicting_refs_keeps_single_block_species`。

### 记录不修复（低风险 / 设计取舍）

* `separate`/`split` 的 `sanitize_filename` 名称碰撞 → 已在前文记录为已知限制。
* `merge_conflicting_refs` 合并后参考序列与其声明 range 的 ungapped 长度可能不一致：
  仅在"同一窗口两个输入 block 参考序列真实不一致"（矛盾输入）时触发，与既有的
  "参考不一致丢弃整窗口"同源，属共享坐标系合并架构下的设计取舍，记录不修复。
* `cover --name` 在"输入无块"（报错）与"有块但无该物种"（输出空 `{}`）两种情形行为
  不一致：均为极边缘信息性问题，无数据损坏，未改。
* `link --best` 函数 doc 措辞 "best-to-best bilateral" 与实现（最近邻+去重）略有出入，
  但 CLI 帮助文本已准确描述为 "nearest-neighbor bilateral links (deduplicated)"，行为与
  文档一致，仅函数注释措辞，未改。

### 验证

`cargo test --test cli_fas --test cli_fas_vars --test cli_fas_poa --test cli_fas_vcf
--test cli_fas_multiz` 全部通过（`cli_fas` 49、`cli_fas_vars` 16、`cli_fas_vcf` 6 等），
`cargo clippy --all-targets -- -D warnings` 与 `cargo fmt --check` 干净。

## 2026-08-09 Review Fixes（第八轮：纵深复核确认无新缺陷）

第八轮对第七轮修复的 `merge_conflicting_refs` 单块物种截断问题做回归验证（新增
`merge_window_conflicting_refs_keeps_single_block_species`，通过），并以全新视角完整
重读全部 20 个 `fas` 子命令（`check`/`concat`/`consensus`/`cover`/`create`/`filter`/
`join`/`link`/`multiz`/`name`/`refine`/`replace`/`separate`/`slice`/`split`/`stat`/
`subset`/`to_vcf`/`to_xlsx`/`variation`）与相关库（`fmt/fas.rs`、`alignment/{coords,
variation,trim,slice,stat}`、`fas_multiz/{merge,banded_align,windows,mod}`、
`fas_xlsx.rs`、`poa/*`、`ds/intspan.rs`），并比照 `docs/fas.md` 与各命令帮助文本。

第八轮纵深补充复核（逐文件精读，未发现新缺陷）：
- `fas_xlsx.rs`：`vars` 以位置为键，sub（全真实碱基列）与 indel（含 gap 列）的起点按
  构造不相交，键冲突不可达；`col_taken = length.min(3)` 与预检
  `col_cursor+width > wrap+1` 一致，`paint_indel` 内不再换行、由调用方统一推进，游标
  `saturating_add`；outgroup 行 `pos_row+seq_count+1` 与 `paint_name` 的 entries 行
  对齐，`sec_height` 多出的一行仅为间隔空白，非缺陷。
- `alignment/variation.rs`：`get_subs`/`get_indels` 等长校验 + `bail!`；`polarize_subs`
  `og_idx < og.len()`、`polarize_indels` `end <= og.len()` 边界；`freq` 取少数等位/极化
  派生计数语义正确；`vcf_alt_bases` 去重排除参考碱基。doctest 全对。
- `alignment/slice.rs`：`ss_start < 1 || ss_end > seq_len_i32` 在切片前守卫，杜绝
  `(ss_start-1) as usize` 下溢；全 gap 物种使子切片退化为纯 indel 岛被 `ss_ints.is_empty()`
  剔除，永不触发空 intspan 的 `align_to_chr`；负链端点交换归一正确。
- `alignment/trim.rs`：`trim_head_tail` 头部移除量 `min(max, seqs[0].len())`、尾部
  `cur_len==0` break，全 gap 不 panic；`head/tail_indel_ints` 空判守卫
  `min()/max()`；`replace_range(lower-1..upper)` 的 lower≥1。
- `fmt/fas.rs` `consensus_block`/`refine_block`、`alignment/msa.rs`、`poa/consensus.rs`：
  全 gap 内群 `retain` 剔除、空 POA 图 `generate_consensus` 返回空、节点权重计入
  路径得分、`refine_block` 外群 `n>=3` 校验、`chop`/`pad`/`fill` 的 `try_from` 与
  `chop*2` 理论溢出（需 chop≥2^63，不可达）。

验证：`cargo test --lib fas`（37 通过，含新增回归）、`cli_fas`（49）、`cli_fas_multiz`
（2）、`cli_fas_vars`（16）、`cli_fas_poa`、`cli_fas_vcf` 全部通过；`cargo clippy
--all-targets -- -D warnings` 与 `cargo fmt --check` 干净。第八轮（含纵深补充复核）
未发现需修复的新缺陷，`fas` 审核收敛。

## 结论

`fas` 命令族审核完成（累计修复 33 处缺陷并补回归测试与文档澄清），并经多轮纵深
复审（全部 20 个子命令的执行路径、`-o` 覆盖保护含 `create`/`check` 的 `.loc` 与
`separate`/`split` 目录输出、`stat`/`variation`/`to-vcf`/`to-xlsx` 的不等长与
外群越界、`slice_block`/`trim_head_tail` 全 gap 场景、`banded_align` I 态 gap
成本、`fas_multiz` 合并顺序确定性、`to-xlsx` 列游标/颜色索引/外群歧义碱基样式兜底、
`slice` 负链坐标、POA 共识节点权重与全 gap 首序列、`--required` 去重、`multiz`
倒置参考区间/窗口边界/渐进合并丢物种/覆盖率过滤、`create`/`separate` 越界与非
IUPAC 边界、`slice` 带 gap 参考坐标映射、外部对齐器输出越界 id、`IntSpan` runlist
空 token、`to-xlsx` `--wrap` 算术溢出、`align_to_chr` 空 intspan、`docs/fas.md`
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

多轮纵深复审逐项核验的候选疑点（POA 回溯 `curr_j-1` 下溢、`write_vcf_block` 的
`pos_idx`、`merge.rs` `cut`/`map_a`/`map_b` 越界、`slice.rs` `ss_start`、`variation.rs`
`uniq_indel_seqs`、`align_seqs_quick` 负 pad、`consensus_block` 空序列、`Range::from_str`
大坐标溢出）全部为已守卫或不可达，未产生新修复。自第四轮起连续两轮（第五、六轮）未
发现需修复的新缺陷，`fas` 审核收敛。
