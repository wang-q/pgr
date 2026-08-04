# runlist / rg 命令族代码审核记录（2026-08-04）

对 `pgr runlist` 与 `pgr rg` 两个命令族及相关库文件（`libs/runlist`、
`libs/ds/intspan`、`libs/ds/range`、`libs/fmt/gff`、`libs/io`，以及迁移
spanr 调用的 `libs/pl/repeat`、`cmd_pgr/pl/p2m`、`cmd_pgr/rept/trf`）和
全部测试/文档进行审核。缺陷按类别分组记录；关键修复均附回归测试（见
文末），验证概况见文末"验证"一节。

## 与外部参考实现的语义一致性核对

runlist 家族对照 spanr 0.6.7 源码、rg 家族对照 rgr 源码逐条核对并逐字节
对拍，语义一致。有意差异（均已记录）：

* `stat`/`statop` 输出用 TSV（spanr 用 CSV），字段顺序与数值不变。
* `runlist combine --op intersect`：后续集合缺失的染色体按空集折叠。
* `rg sort` 保留重复行并稳定排序；`rg merge` 簇内按输入顺序输出。
* `rg span` 极端坐标饱和 + clamp 到可表示域（rgr 裸算术回绕/panic）。
* 畸形输入：rgr/spanr 部分输入直接 panic，pgr 统一友好跳过
  （Zero-Panic）。

## 排除的疑点（经核验无需修复）

* 行首/尾空白与 `|Species=Yeast` 等后缀：`Range::from_str` 最左匹配天然
  容忍，各命令一致。
* coitrees 会内部排序，`RgIndex`/`rg merge` 按行序喂入不排序区间不会错乱。
* `stat`/`statop` 的 `--all` 表头字段数在 multi/single 形态下一致。
* `rg span` 饱和在 Range 层保留、CLI 层 `clamp_to_domain`，职责分开。
* `rg span` pad 越过坐标 1 输出空行、`1:-100` 回退整行作 chr：与
  vendored crate/rgr 一致，未改。

## 记录项（未改，低风险 / 待决策）

* `__single__` 哨兵键与同名单集合理论上碰撞（与 spanr 一致）。
* `runlist split` 的键含 `/` 或 `..` 时输出可写到 outdir 之外（与 spanr
  一致）。
* `combine --op intersect` 的空集折叠语义（见上），不追 spanr 逐字节。

## 已知限制（有意保留）

* `stat`/`statop` 对 0 长度染色体输出 inf/NaN（与 spanr 一致）。

## 带点 contig 名截断 bug 的处置结论

spanr 时代 `chr:start-end` 按 `.` 截断 contig 名（`NC_000913.1` → `"1"`），
当时以 `c1..cN` 名字映射规避；解析器迁入后重新评估：**不修改解析语义**。
理由：`name.chr(strand)` 是 `.rg` 物种前缀约定（`S288c.I(-):190-200` →
`I`），且 `Range` 拆分同时服务 FASTA 头与 `psl to-range`；四个 rept 管线的
实际影响已由 `c1..cN` 映射修复验证。若需原生支持带点名，正解是新增严格
解析模式（新特性，不建议按 bug 修）。

## 修复的缺陷（共 43 处）

### 崩溃 / 越界 / 溢出（Zero Panic，19 处）

**解析器尾部 `-` 越界 panic**：`runlist_to_ranges` 遇 `-` 时越界取字节
   （`"1-"`）。修复：先查长度再判 `upper_is_neg`。
**解析器超大数字 i32 溢出**：digits 用 i32 累加，`"99999999999"`
   debug panic。修复：i64 累加 + i32 范围检查。
**反转区间在 `add_pair` panic**（`"5-3"`、`"1--1"` 等）。修复：
   `runlist_to_ranges` 对 `lower > upper` 报 `Bad order`，`valid` 返回 false。
**坐标上限溢出（JSON/.rg/GFF 三入口）**：`upper + 1` 超 i32。修复：
   解析器拒绝 `> POS_INF - 1`；`rg_to_set`/`rg_to_intervals` 跳过越界行；
   `gff runlist` 报 "coordinates out of range"。
**`.rg` 行 start > end panic**。修复：`rg_to_set`/`rg_to_intervals`
   跳过（不改 `is_valid`）。回归 `rg_to_set_skips_reversed_ranges`。
**非法 runlist JSON 值 panic**（`json_to_set` 直接 `IntSpan::from`）。
   修复：`json_to_set`/`json_to_sets` 返回 Result、先 `valid` 校验；
   `io::read_runlist` 同步。
**genome size ≤ 0 或超上限 panic**。修复：`genome_set` 返回 Result。
**span trim/pad 极端 `-n` 溢出 panic**。修复：`saturating_add/sub` +
   clamp。回归 `extreme_ops_do_not_overflow`。
**excise/fill/cardinality 全幅跨度溢出**。修复：i64 计算，
    `cardinality` 饱和到 i32::MAX。
**gff runlist start > end 记录 panic**。修复：命令层报错。
**`pgr runlist` 裸调用 panic**（unreachable）。修复：
    `subcommand_required` + `arg_required_else_help`，exit 2。回归
    `bare_runlist_shows_help`。
**空集 `is_neg_inf`/`is_pos_inf` unwrap panic**。修复：`is_some_and`。
    回归 `infinity_predicates_on_empty_set`。
**`rg span` Range 运算加减溢出**（8 种 op × 近上限坐标 + 极端 `-n`）。
    修复：saturating 算术；`shift_3p` 去掉取负。
**`rg span` pad 路径 `-n=-2147483648` 取负溢出**。修复：
    `saturating_neg()`（span.rs 与 `IntSpan::pad`）。
**解析器 i64 累加溢出**（19+ 位数字串）。修复：每位累加后检查越界
    提前退出。
**`rg span` 合法结果超出可表示坐标域**。修复：CLI 层 `clamp_to_domain`
    夹回 `1..=POS_INF - 1`。
**`IntSpan::holes` 对 `i32::MIN` 集合补集溢出**（`spans()` 边下溢）。
    修复：`holes` 直接取相邻 span 空隙，不经过补集。回归
    `holes_fill_on_i32_min_coordinates_do_not_overflow` + CLI 用例。
**`IntSpan::from`/`add_runlist`/`remove_runlist` 非法输入 panic**（原为
    外部 API 兼容保留）。修复：`from` 非法输入返回空集，两个 runlist 方法
    前置 `valid` 忽略。回归 `invalid_runlists_are_ignored_not_panicked`。
**IntSpan 其余 API 极端输入溢出**：`contains` 的 `n+1`、`spans()` 的
    `i32::MIN` 上边、`at` 的 `abs(i32::MIN)`、`at_pos`/`at_neg`/`index`
    的 span_len。修复：`checked_add`/i64/`unsigned_abs`，`spans()` 跳过
    退化 span，`index` 结果饱和。回归
    `contains_wide_domain_does_not_overflow`、
    `invert_and_complement_on_i32_min_set_do_not_overflow`、
    `indexing_wide_spans_do_not_overflow`。
**IntSpan 最后的参数校验 panic**：`at`/`index`/`slice` 对空集、索引
    0、越界、元素不存在等直接 panic，`add_pair` 对反转区间 panic
    （"Bad order"）。修复：`at`/`index`/`slice` 改为返回 `Option`
    （非法输入为 `None`，`alignment/coords.rs` 两个调用方接
    `ok_or_else` 报错），`add_pair` 对 `lower > upper` 跳过（与
    `from_pairs` 语义一致）。回归
    `reversed_pairs_are_skipped_not_panicked`、
    `invalid_index_arguments_return_none`。

### 输入校验 / 静默错误（3 处）

**multi 文件被当作 single 传参时静默变空集**（`statop`/`compare`）。
   修复：非字符串值报 "runlist value for ... is not a string"。
**`json_to_sets` 混合形态静默丢数据**。修复：flat/multi 两分支对异形
    值均报错。
**删除未使用且有 panic 隐患的 `gff_to_set`**。

### 外部工具与参数 / CLI / 文档（1 处）

**`span -n` 缺帮助文本**。补 "Number of bases to trim or pad; length
    threshold for excise/fill"。

### 迁移遗留 / 行为一致性（5 处）

**rept/p2m 测试仍以 spanr 缺失为跳过条件**（管道已全部内建）。修复：
    移除 spanr 守卫。
**用户文档仍把 spanr 列为依赖**。修复：改为 trf/FastK/Profex 与
    `pgr runlist stat/statop` 示例。
**src 注释残留 spanr 现态描述**。修复：改为 "runlist 解析器/管道"
    （历史出处保留）。
**`#` 注释行在各 rg 子命令间不一致**。修复：所有 rg 命令统一跳过
    `#` 行，文档补充说明。
**`rg merge` 自映射 parity**：part 等于合并串时输出无意义映射行。
    修复：`part == merged` 时跳过。

### 输入校验 / 静默错误（2 处）

**`Range` end 坐标溢出被静默折叠为 start**（`chr1:5-99999999999` →
    点区间）。修复：`parse_i32` 返回 Option，溢出使整行无效。回归
    `overflow_end_is_invalid_not_start` + `command_rg_overflow_end_skipped`。
**`IntSpan::inset` clamp 下界把 i32::MIN 静默改写**。修复：下界改
    i32::MIN。回归 `inset_identity_at_i32_min`。

### 参数校验 / 数据安全 / 错误传播（3 处）

**`rg span` op/mode 校验推迟到逐行**（空输入静默成功）。修复：读取前
    校验 `shift/flank + both`。回归
    `command_rg_span_invalid_mode_checked_before_input`。
**流式命令 `-o` 与输入同路径时先截断输入**。修复：`same_path` 检查
    拒绝（后升级为 canonicalize + dev/inode 覆盖别名）。回归
    `command_rg_output_same_as_input_rejected`、
    `command_runlist_convert_output_same_as_input_rejected`。
**`IntSpan::write_json` 未显式 flush**（写盘失败静默 0 退出）。修复：
    末尾 `writer.flush()?`。

### 性能（1 处）

**`rg merge` 去重 O(n²) 退化**（100k 区间 37.8 s）。修复：HashSet 判重
    （0.28 s）。基准 `benches/rg_merge_benchmark.rs`（10k ≈ 2.7 ms、
    50k ≈ 18 ms）。回归 `command_rg_merge_dedups_identical_lines`。

### 解析一致性 / 数据安全（2 处）

**`Range` 扫描器对非 ASCII 字符名静默失配**（Unicode `\w`）。修复：
    字符类按 UTF-8 字符判定，数字仍只收 ASCII。回归
    `regex_and_manual_decoders_agree` + `command_rg_unicode_contig_names_parsed`。
**`-o` 同输入检查被 symlink/hardlink 别名绕过**。修复：canonicalize +
    dev/inode 比较。回归 `command_rg_output_alias_of_input_rejected`。

### 与 rgr 的数值一致性（1 处）

**`rg merge` 覆盖度用 f64 判据，0.8 边界与 rgr 的 f32 分叉**。修复：
    比值改用 f32 算术。回归 `command_rg_merge_exact_threshold_parity`。

### 数据安全 / 溢出（2 处）

**五个流式命令在后续输入打开失败时截断输出**。修复：先打开/读取全部
    输入再创建 writer。回归
    `command_rg_output_preserved_on_missing_input`、
    `command_runlist_convert_output_preserved_on_missing_input`。
**`IntSpan::covered` 近全幅查询 per-span 累计 i32 溢出**。修复：i64
    相减。回归 `covered_wide_domain_does_not_overflow`。

### 数据安全 / 参数校验（2 处）

**`ensure_outfile_distinct` 对屏幕哨兵 `stdout` 误判**。修复：
    `outfile == "stdout"` 跳过检查。回归
    `command_rg_stdout_named_input_allowed`。
**输入侧流哨兵 `stdin` 被同路径检查误拒**（`stdout` 哨兵误判的镜像）。
    修复：跳过字面 `stdin` 输入。回归
    `command_rg_stdin_sentinel_output_allowed`。

### 数据安全（1 处）

**目录作为输入通过打开探针、读取失败后输出被截断**。修复：
    `libs/io::reader` 打开前拒绝目录。回归
    `command_rg_output_preserved_on_directory_input`。

### 文档修复

* `notes/design/runlist.md`：子命令数 12 → 10（cover/coverage 迁出、
  gff 归位）。
* `notes/project-understanding.md`：补 `rg`/`runlist` 家族与
  `libs/runlist` 目录、§9 区间链路、rgr-tva-audit 状态。
* `notes/design/rgr-tva-audit.md`："待定稿"/"不再移植" 更新为已实现。

## 验证

* 差分对拍：与 spanr/rgr 逐字节 2000+ trial、与朴素实现 3000+ trial、
  `Range` 扫描器 vs 参考正则 20 万 trial（8 种 Unicode 文字系统），全部
  一致。
* 畸形输入 fuzz：累计 5000+ trial（JSON/.rg/GFF/sizes/二进制、超大数字、
  反转区间、极端 `-n`、近上限坐标、缺失/目录输入），零 panic。
* 数据安全：`-o` 同路径、缺失/目录输入、`stdout`/`stdin` 哨兵等修复
  前后均实测复现，既有输出原样保留。
* 性能：`rg merge` 100k 37.8 s → 0.28 s；criterion 10k ≈ 2.7 ms、
  50k ≈ 18 ms；`rg runlist` 20 万行 × 5 万 span 0.29 s。
* `cargo fmt`、`cargo clippy --all-targets -- -D warnings` 干净；全部测试
  + doctest 通过；docs 示例逐一执行通过。
* 新增回归测试（主要）：runlist 阶段 14 个函数；rg 阶段
  `command_rg_span_extreme_no_panic`、`command_rg_comments_skipped`、
  `extreme_ops_do_not_overflow`、`overflow_end_is_invalid_not_start`、
  `inset_identity_at_i32_min`、`command_rg_overflow_end_skipped`、
  `regex_and_manual_decoders_agree`、`command_rg_unicode_contig_names_parsed`、
  `command_rg_output_alias_of_input_rejected`、
  `covered_wide_domain_does_not_overflow`、
  `command_rg_output_preserved_on_missing_input`、
  `command_runlist_convert_output_preserved_on_missing_input`、
  `command_rg_stdout_named_input_allowed`、
  `command_rg_output_preserved_on_directory_input`、
  `command_rg_stdin_sentinel_output_allowed`、
  `invalid_runlists_are_ignored_not_panicked`、
  `contains_wide_domain_does_not_overflow`、
  `invert_and_complement_on_i32_min_set_do_not_overflow`、
  `indexing_wide_spans_do_not_overflow`、`reversed_pairs_are_skipped_
  not_panicked`、`invalid_index_arguments_return_none` 等。

## 后续补充（2026-08-04 复核）

`rg runlist` 再次核对 rgr 源码：`superset` 语义确认为 `set.superset(range)`
（range 完全落在 runlist 内，`size == length`），与 rgr 逐字节一致。补充：

* 文档：`rg runlist` 的 `after_help` 与 `docs/rg.md` 补 `--op superset`
  示例并澄清其"range 完全包含在 runlist 内"的语义（原描述仅说
  "contained"，未点名 `superset`，且 op 名易与"range 是 runlist 超集"混淆）。
* 测试：原 `command_rg_runlist` 只覆盖精确包含与完全不重叠，无法区分
  `superset` 的正确语义（range 内含于 runlist）与另一种解释（range 包含某
  span）。新增 `command_rg_runlist_superset_partial_overlap`，用
  `chr1:5-25`（含 span 但非内含）锁死"range 在 runlist 内"语义，并断言
  overlap / non-overlap 的包含/部分重叠/不相交三态输出。
