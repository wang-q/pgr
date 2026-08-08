# runlist / rg 命令族代码审核记录（2026-08-04）

对 `pgr runlist` 与 `pgr rg` 两个命令族及相关库文件（`libs/runlist`、
`libs/ds/intspan`、`libs/ds/range`、`libs/fmt/gff`、`libs/io`，以及迁移
spanr 调用的 `libs/pl/repeat`、`cmd_pgr/pl/p2m`、`cmd_pgr/rept/trf`）和全部
测试/文档进行审核。以下仅保留有借鉴意义的结论；逐轮验证过程已精简。

## 与外部参考实现的语义一致性核对

runlist 家族对照 spanr 0.6.7 源码、rg 家族对照 rgr 源码逐条核对并逐字节对拍。
有意差异（均已记录）：
- `stat`/`statop` 输出用 TSV（spanr 用 CSV），字段顺序与数值不变。
- `runlist combine --op intersect`：后续集合缺失的染色体按空集折叠。
- `rg sort` 保留重复行并稳定排序；`rg merge` 簇内按输入顺序输出。
- `rg span` 极端坐标饱和 + clamp 到可表示域（rgr 裸算术回绕/panic）。
- 畸形输入：rgr/spanr 部分输入直接 panic，pgr 统一友好跳过（Zero-Panic）。

## 排除的疑点（安全不变量，经核验无需修复）

- coitrees 会内部排序，`RgIndex`/`rg merge` 按行序喂入不排序区间不会错乱。
- `rg count` 的 COITree 区间为**闭区间** `[start, end]`，仅端点相接
  （`chr1:10-20` 对 `chr1:20-20`）仍计为重叠。
- `depth_runs` 扫描线：事件合并、`run_depth` 归属、`pos>s` 关段、尾部开口
  区间收尾正确；`by_level` 键为精确深度串，深度受 `>= min_depth` 过滤不会为负。
- `Range` 扫描器与参考正则一致：最左匹配、贪婪分组、`end=0` 缺省、溢出即无效、
  回退取首个空白 token、Unicode/UTF-8 字宽。
- `rg_merge_mapping` 的 DSU/COITree 聚类正确、`part == merged` 自映射跳过、
  HashSet 判重、f32 比例判据。
- `IntSpan` 边界：坐标受 `POS_INF-1` 约束，`upper+1`/`find_pos(upper+1)`/
  `edge-1` 均不溢出；`covered` 两次二分 + i64 累加 + 饱和；`holes` 直接取空隙。
- `rg runlist`/`rg prop`/`rg count` 的 `length = end - start + 1` 在
  `usable_range` 约束下（start ≥ 1、end ≤ POS_INF-1）不溢出 i32。

## 带点 contig 名截断 bug 的处置结论

spanr 时代 `chr:start-end` 按 `.` 截断 contig 名（`NC_000913.1` → `"1"`），当时
以 `c1..cN` 名字映射规避；解析器迁入后重新评估：**不修改解析语义**。理由：
`name.chr(strand)` 是 `.rg` 物种前缀约定，且 `Range` 拆分同时服务 FASTA 头与
`psl to-rg`；四个 rept 管线的实际影响已由 `c1..cN` 映射修复验证。若需原生支持
带点名，正解是新增严格解析模式（新特性，不建议按 bug 修）。

## 记录项（未改，低风险 / 待决策）

- `__single__` 哨兵键与同名单集合理论上碰撞（与 spanr 一致）。
- `runlist split` 的键含 `/` 或 `..` 时输出可写到 outdir 之外；outdir 为输入
  所在目录且某键等于输入 basename 时，输出会覆盖输入文件（先读全量再写，计算
  不受影响，仅磁盘输入被覆盖）。与 spanr 一致，窄边角，未加 `ensure_outfile_
  distinct`（split 输出为目录内多文件，逐路径检查成本高于收益）。
- `runlist merge`（未加 `--all`）两个输入首段 stem 相同时，后者静默覆盖前者。
  属键控方案固有行为（`--all` 用完整 stem 可规避）。
- `-o stdin` 会被 `pgr::writer` 当作字面文件名创建（`writer` 只对 `stdout`
  哨兵特判）。属全局约定（输出用 `stdout`、输入用 `stdin`），非本家族特有。
- **库级观察（范围外）**：`IntSpan::find_islands_n` 的
  `self.find_pos(val + 1, 0)` 在 `val == i32::MAX` 时 `val + 1` 溢出（与
  `contains` 已用 `checked_add` 修复的模式不一致）。该函数仅被
  `libs/alignment/slice.rs`（`fas slice`）调用，且调用处坐标受序列长度（i32
  检查）约束，当前不可达；若未来在 `val` 无上界的路径使用需先加 `checked_add`。
  属 alignment 子系统，留待后续审计。

## 已知限制（有意保留）

- `stat`/`statop` 对 0 长度染色体输出 inf/NaN（与 spanr 一致）。

## 修复的缺陷（根因模式）

### Zero-Panic / 溢出（IntSpan 与解析器）

- **解析器/集合运算多处整数溢出与越界 panic**（runlist 数字串 i32/i64 累加溢出、
  反转区间 `add_pair`、坐标上限 `> POS_INF-1`、`.rg` 行 start>end、非法 JSON 值、
  `genome size ≤ 0`、span trim/pad 极端 `-n`、excise/fill/cardinality 全幅跨度、
  `rg span` 8 种 op × 近上限坐标、`holes` 补集 i32::MIN、`contains` 的 `n+1`、
  `spans()` 的 i32::MIN 上边、`at` 的 `abs(i32::MIN)`、`index` 的 span_len、
  `at`/`index`/`slice` 对空集/越界直接 panic）。修复：统一改
  `checked_add`/i64/`unsigned_abs`/`saturating`，解析器拒绝 `> POS_INF-1`、
  `lower > upper` 报错，`at`/`index`/`slice` 返回 `Option`、`add_pair` 跳过反转
  区间，`holes` 直接取相邻 span 空隙（不经补集），`cardinality` 饱和到 i32::MAX。
  回归 `extreme_ops_do_not_overflow`、`holes_fill_on_i32_min_coordinates_do_not_
  overflow`、`contains_wide_domain_does_not_overflow`、
  `invert_and_complement_on_i32_min_set_do_not_overflow`、
  `indexing_wide_spans_do_not_overflow`、`reversed_pairs_are_skipped_not_panicked`、
  `invalid_index_arguments_return_none`、`invalid_runlists_are_ignored_not_panicked`。
- **`pgr runlist` 裸调用 panic**（无子命令时 unreachable）。修复：
  `subcommand_required` + `arg_required_else_help`，exit 2。回归
  `bare_runlist_shows_help`。
- **空集 `is_neg_inf`/`is_pos_inf` unwrap panic**。修复：`is_some_and`。

### 输入校验 / 静默错误

- **multi 文件被当作 single 传参时静默变空集**（`statop`/`compare`）。修复：非
  字符串值报错；`json_to_sets` 混合形态静默丢数据 → 两分支对异形值均报错。
- **`Range` end 坐标溢出被静默折叠为 start**（`chr1:5-99999999999` → 点区间）。
  修复：`parse_i32` 返回 Option，溢出使整行无效。回归
  `overflow_end_is_invalid_not_start`。
- **`IntSpan::inset` clamp 下界把 i32::MIN 静默改写**。修复：下界改 i32::MIN。

### 数据安全（`-o` 同输入保护 / 输出截断）

- **流式命令 `-o` 与输入同路径时先截断输入**；**`-o` 同输入检查被 symlink/
  hardlink 别名绕过**。修复：`same_path` 用 canonicalize + dev/inode 比较。
- **五个流式命令在后续输入打开失败时截断输出**。修复：先打开/读取全部输入再
  创建 writer。回归 `command_rg_output_preserved_on_missing_input`、
  `command_runlist_convert_output_preserved_on_missing_input`。
- **目录作为输入通过打开探针、读取失败后输出被截断**。修复：`libs/io::reader`
  打开前拒绝目录。回归 `command_rg_output_preserved_on_directory_input`。
- **`ensure_outfile_distinct` 对屏幕哨兵 `stdout` 误判 / 输入侧流哨兵 `stdin`
  被同路径检查误拒**。修复：跳过字面 `stdout` 输出与字面 `stdin` 输入。
- **`-o` 覆盖保护缺口**：`rg runlist`/`rg prop` 未含 runlist.json 参考文件；
  `rg cover`/`rg merge`/`rg coverage` 未调用 `ensure_outfile_distinct`；
  `runlist merge/compare/some/combine/span` 五个 JSON 输出命令缺保护。修复：就地
  补齐（`rg sort` 原地排序安全合理，不加入检查）。回归
  `command_rg_output_same_as_input_rejected`/`command_runlist_output_same_as_
  input_rejected` 各用例。

### 性能

- **`rg merge` 去重 O(n²) 退化**（100k 区间 37.8 s → HashSet 0.28 s）。基准
  `benches/rg_merge_benchmark.rs`。

### 解析一致性 / 与 rgr 数值一致性

- **`Range` 扫描器对非 ASCII 字符名静默失配**（Unicode `\w`）。修复：字符类按
  UTF-8 字符判定，数字仍只收 ASCII。回归 `command_rg_unicode_contig_names_parsed`。
- **`rg merge` 覆盖度用 f64 判据，0.8 边界与 rgr 的 f32 分叉**。修复：比值改用
  f32 算术。回归 `command_rg_merge_exact_threshold_parity`。

### 迁移遗留 / 行为一致性

- **rept/p2m 测试仍以 spanr 缺失为跳过条件**（管道已全部内建）→ 移除守卫；
  用户文档仍把 spanr 列为依赖 → 改为 trf/FastK/Profex 与 `pgr runlist
  stat/statop` 示例；src 注释残留 spanr 现态描述 → 改为 "runlist 解析器/管道"。
- **`#` 注释行在各 rg 子命令间不一致**。修复：所有 rg 命令统一跳过 `#` 行。
- **`rg merge` 自映射 parity**：`part == merged` 时输出无意义映射行 → 跳过。
- **`rg span` op/mode 校验推迟到逐行**（空输入静默成功）→ 读取前校验
  `shift/flank + both`。
- **`IntSpan::write_json` 未显式 flush**（写盘失败静默 0 退出）→ 末尾 `flush?`。
- **`runlist split -o stdout` 丢弃键名**（无法区分同值不同键）→ stdout 模式输出
  `key\tvalue`。
- **`runlist merge --all` 帮助文本与代码/文档不一致**（--all 用完整 stem 作 key）
  → 帮助文本修正。
- **删除未使用且有 panic 隐患的 `gff_to_set`**。

## 结论

`rg`/`runlist` 两个命令族审核完成（累计修复 48 处缺陷），经纵深复核收敛；与
spanr/rgr 逐字节对拍、畸形输入 fuzz 零 panic，`cargo fmt`/`clippy` 干净。
