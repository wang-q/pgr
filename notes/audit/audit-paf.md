# paf 命令族代码审核记录（2026-08-05）

对 `pgr paf` 命令族（`index` / `query` / `to-bed` / `to-fas` / `to-maf` /
`to-vcf` / `to-gfa` / `graph` / `stat`，共 9 个子命令）及其库文件
（`libs/paf`：`parser`、`record`、`cigar`、`index`（`builder`/`query`/`bfs`/
`persist`）、`fasta`、`query`、`msa_build`、`to_maf`、`to_fas`、`vcf`、
`poa_compact`、`graph`（`builder`/`dsu`/`segment`/`gfa`/`report`）、
`maf_import`）和全部测试/文档进行审核。缺陷按类别分组记录；关键修复均附
回归测试（见文末），验证概况见文末"验证"一节。

## 排除的疑点（经核验无需修复）

* `query` 的 subset / syntenic 过滤中使用 `idx.id_to_name(*qid).unwrap_or("")`
  并非 panic——`unwrap_or` 带默认值，`qid` 无效时返回空串而非崩溃。`output_paf`
  同理用 `unwrap_or("?")`。均安全。
* `index/query.rs::project` 空 CIGAR 分支：`len <= 0` 时返回 `None`，不会返回
  空/反向区间；带 CIGAR 分支在 `found && q_min < q_max && t_min < t_max` 才
  `Some`，仅覆盖 I（插入）或仅覆盖 D（删除）的区间因 query 或 target 侧未
  累积到有效区间而正确返回 `None`。边界处理完备。
* `parse_region` / `load_bed_regions` 对 `start`/`end` 均做了非负与 `end > start`
  校验，非法输入返回友好错误，不 panic。
* `graph` / `stat` 的 `--min-var-len` 在 execute 层 `anyhow::ensure!(min_var_len > 0)`
  校验；S 行在 topology-only 模式（无 `-f`）输出 `S <id> * LN:i:<len>`
  （无 SN/SO），有 FASTA 时输出 `SN:Z`/`SO:i`/`SR:i:0`，与 `docs/paf.md` 描述一致。
* 全量扫描 `cmd_pgr/paf` 与 `libs/paf` 的 `unwrap()`/`unreachable!`：全部为
  clap `required` 参数、`value_parser` 约束后的枚举，或防御性分支（如
  `left_align_indels` 前对 `t_aln_pos <= 0` 的 guard），运行期不可达，无
  Zero-Panic 隐患。
* `output_paf` 的 `matches`/`block_len`/`gi`/`bi`/`cg` 取完整源 CIGAR 而非投影
  子区间，属文档已声明的已知限制（`libs/paf/query.rs` 头注释），非 bug。
* `--merge-distance > 0` 时强制要求 `--fasta-tsv`（`anyhow::bail!`），与
  `docs/paf.md` 一致。
* `to-vcf` 的 GT 编码（0=REF、1..=N=ALT、`.`=gap/非 ACGT）与 `docs/paf.md` 和
  `after_help` 描述一致；indel 左对齐规则（锚点前参考碱基 == 每个非空 indel
  序列末位时左移）与文档一致。

## 记录项（未改，低风险 / 待决策）

* `output_paf`/`to-bed` 输出的 `query_length`/`target_length`（PAF 第 2/7 列）
  恒为 `0`，因 `PafIndex` 不保留每序列总长。若要填充需改索引格式持久化
  `src_size`。属已知限制（代码注释已声明），跨格式变更，未改。
* `to-vcf` 的部分删除样本（仅删除 DEL 区一部分）回退为 best-effort 等位基因
  且不做完全左对齐（`docs/paf.md` 与 `after_help` 已声明）。行为正确，未改。
* `query --transitive --merge-distance` 依赖 `--fasta-tsv` 重算 CIGAR；若用户
  未提供 `-f` 会硬报错。属有意约束，未改。

## 修复的缺陷（共 13 处）

### 坐标投影正确性（Zero-Panic / 精度，3 处）

**空 CIGAR 投影左越界时 over-count 长度**（`index/query.rs::project`）。当查询
区间起点落在记录 target 区间左侧（如实例 `[0,50)` 对记录 target `[20,80)`），
空 CIGAR 分支原按查询区间原始宽度投影，导致投影到 query 上过长（50bp 而非
实际重叠 30bp）。修复：将查询跨度 clamp 到实际重叠 `[ts,te)` 与记录 target
extent 的交集，再按重叠计算投影长度。回归 `test_project_empty_cigar_left_
overhang_clamps_len`、`test_project_empty_cigar_left_overhang_minus_strand`。

**空 CIGAR 且 `-` 链的子区间投影坐标错误**。`-` 链 PAF 记录 `query_start/end`
是正向坐标，但比对以反向读取 query，target 偏移 `off` 应映射到 query 的
*末端* 侧：`[query_end - off - len, query_end - off)`。空 CIGAR 分支原先未镜像
（与带 CIGAR 的 `rc_to_forward` 不一致），导致 `-` 链投影到错误的正向坐标。
修复：空 CIGAR 分支对 `-` 链同样做末端镜像。回归
`test_project_empty_cigar_minus_strand_subinterval`。

**`maf to-paf` 反向链区间越界可能 usize 下溢**（`maf_import.rs`）。`-` 链 MAF
中 `start + size > src_size` 时，`reverse_range_pair(start, start+size, ...)`
的 `src_size - end` 会 usize 下溢（debug panic）。修复：调用前
`checked_add` + 越界 `bail!` 友好报错。回归
`test_minus_strand_interval_exceeding_src_size_rejected`。

### 数据安全（`-o` 覆盖保护，5 处）

`paf` 家族各输出命令此前未调用 `ensure_outfile_distinct`，`-o` 指向输入文件
会把 PAF/TSV 输入覆盖成输出（exit 0，静默数据丢失）。修复：在 `index`（保护
全部 infiles）、`query` / `to-bed`（保护 infile）、`to-fas` / `to-maf` /
`to-vcf` / `to-gfa`（保护 infile + fasta_tsv）、`graph` / `stat`（保护 infile
+ fasta_tsv）中就地补齐。回归
`command_paf_index_output_same_as_input_rejected`、
`command_paf_query_output_same_as_input_rejected`、
`command_paf_to_bed_output_same_as_input_rejected`、
`command_paf_graph_output_same_as_input_rejected`，断言输入文件未被改动。

### 参数校验（CLI，3 处）

**POA 评分参数接受非法值**（`args.rs`）。`--match` 可为负数/0、罚分参数可为
正数，POA 引擎会静默产生逆向/无效打分。修复：新增 `parse_match_score`
（`--match` 必须 > 0，正奖励）与 `parse_poa_penalty`（`--mismatch` /
`--gap-open` / `--gap-extend` 必须 <= 0，罚分），在 clap 解析期拒绝。回归
`command_paf_poa_score_params_validated`。

**查询过滤参数接受负值被静默当作"关闭"**（`args.rs`）。`--min-dist` /
`--min-output-len` / `--merge-distance` / `--min-chain-length` 为负数时被查询
管线当作 0/关闭处理，属用户误用且无提示。修复：新增 `parse_non_negative_i32`
在解析期拒绝 `>= 0` 之外的值。回归 `command_paf_negative_query_filters_rejected`
（覆盖空格分隔与 `=` 两种语法）。

**`--min-identity` 越界未被校验**（`args.rs`）。此前用 `clap::value_parser!(f64)`
接受任意浮点，越 `[0.0,1.0]` 的值照常参与过滤。修复：新增 `parse_min_identity`
在解析期拒绝越界值。回归 `command_paf_query_min_identity_out_of_range_rejected`。

### 变体正确性 / 对齐（2 处）

**VCF 中 DEL 在 INS 之后被静默丢弃**（`vcf.rs`）。当连续的 target gap（INS）后
紧跟 target 非 gap 但有 query gap（DEL）时，旧代码取 `msa[0][col_start-1]`
（恰为 INS 的 target gap）作锚点并丢弃该 DEL。修复：提取 MSA 变异走查为
`emit_msa_variants`，DEL 锚点取 INS 之前最后的 target 非 gap 碱基（经
`t_aln_pos`），并补单元测试。回归 `test_emit_msa_variants_del_after_ins`
（手工构造 `A--CGT`/`ATT-GT` 的 INS 后接 DEL 的 MSA）。

**MSA 构建对相接片段不合并，导致 per-name 去重丢碱基**（`msa_build.rs`）。
`build_msa_entries` 的合并判据用 `qs < last.2`，当 target 删除把一个连续 query
区域切成两个恰好相接的片段（`qs == last.2`）时不会合并，随后 per-name 去重会
静默丢弃第二片的碱基。修复：判据改 `qs <= last.2`。回归
`test_build_msa_entries_merges_touching_fragments`。

### 输出坐标 / 图节点序列（2 处）

**`to-fas` pairwise `-` 链头坐标用正向区间而非反向坐标**（`to_fas.rs`）。pairwise
FAS 头原先用 `q_start_fwd`/`q_end_fwd`（正向 query 区间），对 `-` 链显示的
是反向互补序列，头坐标应为其反向坐标（`src_size - qe`），与 `to-maf`/MSA 输出
一致。修复：改用 `q_start_maf` 推导显示坐标。回归
`command_paf_to_fas_reverse_strand_partial_query_uses_rev_coords`。

**`graph` 节点序列未按代表段正向取向存储**（`graph/builder.rs`）。当 DSU 组件
的代表段序列缺失、而由 `-` 取向的成员段填充节点时，原代码直接存成员段正向区
域，导致该 `-` 步长与存储序列不一致。修复：由 `-` 取向段填充时对区域做反向互
补，使节点序列保持在代表段正向取向。回归
`test_node_sequence_revcomp_when_filling_from_minus_segment`。

## 文档修复

* `docs/paf.md`：引言"查询类子命令"列表补上遗漏的 `to-fas`（原先只列
  `query` / `to-bed` / `to-maf` / `to-vcf` / `to-gfa`，与"通用 query 选项"一节
  及首段功能概览不一致）。

## 验证

* `cargo fmt` 与 `cargo clippy --all-targets -- -D warnings` 干净。
* `cargo test`：569 个单元测试 + 71 个 doctest + 全部集成测试通过，零失败。
  paf 相关：`cli_paf`（19）、`cli_paf_to_fas`（11）、`cli_paf_to_maf`（15）、
  `cli_paf_to_vcf`（6）、`libs::paf::*`（含持久化、索引、BFS、MSA、VCF 等）
  全部通过。
* 新增回归测试（主要）：`test_project_empty_cigar_left_overhang_clamps_len`、
  `test_project_empty_cigar_left_overhang_minus_strand`、
  `test_project_empty_cigar_minus_strand_subinterval`、
  `test_minus_strand_interval_exceeding_src_size_rejected`、
  `test_negative_score_preserved_in_ms_tag`、
  `test_build_msa_entries_merges_touching_fragments`、
  `test_emit_msa_variants_del_after_ins`、
  `test_node_sequence_revcomp_when_filling_from_minus_segment`、
  `command_paf_poa_score_params_validated`、
  `command_paf_negative_query_filters_rejected`、
  `command_paf_query_min_identity_out_of_range_rejected`，以及
  `command_paf_index/query/to_bed/graph_output_same_as_input_rejected` 的
  `-o` 同输入保护各用例（断言输入文件未被改动）。
* 复审轮（2026-08-05）复跑 `cargo test`、`cargo clippy --all-targets -- -D
  warnings`、`cargo fmt --check`，均干净；对 `project`、`query` 过滤、`vcf`
  左对齐、`graph` 节点取向等关键路径做了防御性逐行复核，未再发现新问题。

## 结论

`paf` 命令族审核完成（累计修复 13 处缺陷、补回归测试与文档澄清），并经多轮
纵深复核（`index/project` 坐标投影、BFS 传递遍历、`query` 过滤、`vcf` 变异
调用、`msa_build` 合并、`graph` 节点取向、`maf_import` 反向链、`-o` 覆盖保护、
POA/查询参数校验）收敛，未再发现新问题。