# runlist / rg 命令族代码审核记录（2026-08-04）

对新增的 `pgr runlist` 与 `pgr rg` 两个命令族及相关的库文件进行多轮深入
审核。范围：`cmd_pgr/runlist/` 10 个子命令（combine/compare/convert/genome/
merge/some/span/split/stat/statop）、`cmd_pgr/rg/` 8 个子
命令（cover/coverage/count/merge/prop/runlist/sort/span）、
`cmd_pgr/gff/runlist`、`libs/runlist`、迁入的 `libs/ds/intspan` 与手写扫描
器版 `libs/ds/range`、`libs/fmt/gff`、`libs/io::read_runlist`，以及迁移了
spanr 调用的 `libs/pl/repeat`、`cmd_pgr/pl/p2m`、`cmd_pgr/rept/trf` 和全部
测试/文档。

轮次：第 1–5 轮审 runlist 家族（#1-19），第 6 轮审 rg 家族并复核共享库
（#20-25）。每轮发现问题后修复并进入下一轮复核；第 5 轮与第 6 轮均经全量
重读 + 随机畸形输入 fuzz + 差分对拍未再发现新问题后收束。最终
`cargo fmt`、`cargo clippy --all-targets -- -D warnings` 干净，全部测试
二进制 + doctest 通过（1178 个断言）。

第 7 轮（本文件定稿后追加）：全量重读 `rg` / `runlist` 家族与
`libs/ds/range.rs`、`libs/ds/intspan.rs`、`libs/runlist/mod.rs`，用暴力
实现差分对拍（集合运算、深度扫描线、`covered`、`rg merge` 聚类、Range
极端运算、CLI 近上限坐标 200 trial）与 300 trial × 20 条命令的畸形输入
fuzz 复核，新发现 2 处缺陷（#26、#27）并修复；修复后重跑全部验证。

第 8 轮（复核轮）：对 #26/#27 修复做极端坐标差分对拍（[i32::MIN,
POS_INF-1] 域 3000 trial 的 union/intersect/diff/xor/trim/pad，含无限臂
与饱和语义）与 120 trial × 18 条命令的全新种子 fuzz，另手工核验
`stat` / `statop` 的 `--all` 表头与行列对齐——未再发现新问题，收束。

第 9 轮（终审轮）：换全新角度复核——5000 组随机集合的结构不变量
（edges 严格递增 / 偶数 / span 不相邻）+ Display 往返解析一致、150 组
CLI 往返（`rg cover` → `runlist convert` → `rg cover` 并集一致、
`runlist span trim/pad -n 0` 恒等）、5000 组 `Range` 解析-重编码
往返——未发现新问题，确认收束条件成立。

第 10 轮（复核轮）：换全新角度复核——CLI 参数组合的校验时机、输出文件
与输入文件同路径时的数据安全、JSON 写出的错误传播、帮助文本完整性，并
重跑四组差分对拍（rg count/prop/runlist/coverage 250 trial、runlist
compare/combine/span/convert/some/split/genome 200 trial、rg merge
300 trial、集合运算与 covered 150 trial）与两轮新种子畸形输入 fuzz
（120 trial × 2）——发现 3 处缺陷（#28-#30）并修复；修复后重跑全部验证。

第 11 轮（性能复核轮）：放大输入规模检查各命令复杂度，发现 `rg merge`
去重为 O(n²)（`Vec::contains`），100k 条单染色体区间耗时 37.8 s
（debug）；修复后同输入 0.28 s（约 135 倍）。新增 criterion 基准
`benches/rg_merge_benchmark.rs`（disjoint / clustered × 10k / 50k，
release 下 10k ≈ 2.7–2.9 ms、50k ≈ 17.8–18.4 ms，近线性），并补充
去重语义回归测试。

第 12–15 轮（复核轮，未发现新缺陷）：

* 第 12 轮：把此前只手工核对的 `runlist stat` / `statop` 纳入差分对拍
  （150 trial，覆盖 single/multi、`--all`、`--base`、四种 op，逐字段
  一致）；multi 形态 `compare` 四种 op 120 trial 一致；500 区间链式
  重叠验证 `rg merge` 传递聚类（全部映射到同一 merged）；将第 10 轮
  五处同路径检查提取为共享 helper `ensure_outfile_distinct`
  （`cmd_pgr/args.rs`，消除 5 份重复）。
* 第 13 轮：对 `rg span` 全部 Range 运算（trim / trim_5p / trim_3p /
  shift_5p / shift_3p / flank_5p / flank_3p / excise + CLI 层
  `clamp_to_domain`）做逐公式差分对拍，三个种子共 700 trial（含近上限
  坐标、负 `n`、i32 极值、name/strand 前缀），全部一致。
* 第 14 轮：复核全部改动 diff；裸命令 `pgr rg` / `pgr runlist` 与未知
  子命令均 exit 2 并显示帮助；全量测试 + doctest 通过。
* 第 15 轮：`docs/rg.md` 与 `docs/runlist.md` 的全部示例命令逐一执行
  验证（含 stdin、`--longest`、`--full`、`--all`、`--suffix`、
  `--base` 等形态），全部按文档工作。

## 修复的缺陷（共 31 处）

修复按发现顺序全局编号（#1-19 为 runlist 阶段、#20-25 为 rg 阶段、
#26-27 为第 7 轮），
按类别分组。

### 崩溃 / 越界 / 溢出（Zero Panic，16 处）

1. **IntSpan runlist 解析器尾部 `-` 越界 panic**：`runlist_to_ranges` 在
   遇到 `-` 时无条件 `bytes.get(idx + i + 1).unwrap()`，输入 `"1-"` 直接
   越界 panic（`valid("1-")` 也 panic）。修复：先检查 `idx + i + 1 < len`
   再判断 `upper_is_neg`。
2. **解析器超大数字 i32 溢出**：digits 累加用 `i32`，`"99999999999"`
   在 debug 构建 panic（release 静默回绕）。修复：改为 i64 累加并做
   i32 范围检查，越界返回 "Number format error: out of range"。
3. **反转区间在 `add_pair` panic**：`"5-3"`、`"1-0"`、`"1--1"` 解析通过
   `is_valid` 后进入 `add_pair` 触发 "Bad order" panic。修复：
   `runlist_to_ranges` 对 `lower > upper` 直接返回 `Bad order` 错误，
   `IntSpan::valid` 对这些输入返回 false（`add_runlist("1--1")` 的 panic
   消息仍含 "Bad order: 1,-1"，保留原测试契约）。
4. **坐标上限溢出（JSON / .rg / GFF 三入口）**：`add_pair` 把
   `upper + 1` 存为边，坐标 ≥ 2147483646（`POS_INF - 1` 之上）时
   `find_pos(upper + 1)` 溢出 panic（如 runlist `"1-2147483646"`、
   `.rg` 行 `chr1:2147483647-2147483647`、GFF 记录 end=2147483647）。
   修复：解析器拒绝 `upper_val > POS_INF - 1`；`rg_to_set` /
   `rg_to_intervals` 跳过越界坐标；`gff runlist` 命令层报
   "coordinates out of range"。回归测试覆盖三个入口。
5. **`.rg` 行 start > end panic**：`chr1:10-5` 通过 `Range::is_valid`
   （仅查 start != 0）后 `add_pair(10, 5)` panic。修复：
   `rg_to_set` / `rg_to_intervals` 跳过 `start > end` 的行（`Range` 的
   `is_valid` 语义被其他命令依赖，如 `2bit range` 对反转区间有专门报错，
   故不改 `is_valid` 本身）。回归测试
   `rg_to_set_skips_reversed_ranges`。
6. **非法 runlist JSON 值 panic**：`json_to_set` 直接用
   `IntSpan::from(s)`，`"abc"`、`"1-"` 等触发 `add_runlist` 的
   `unwrap()` panic（`pgr runlist span` 等命令均可复现）。修复：
   `json_to_set` / `json_to_sets` 改为返回 `anyhow::Result`，先
   `IntSpan::valid` 校验再构造，非法值报 "invalid runlist for ..."；
   `libs/io::read_runlist`（`fa mask` / `fas slice` 共用）同步校验。
8. **genome 染色体 size ≤ 0 或超上限 panic**：`chr.sizes` 中出现 0 / 负数
   / ≥ 2147483646 时 `IntSpan::from_pair(1, v)` panic。修复：
   `genome_set` 返回 `Result`，分别报 "invalid chromosome size" /
   "out of range"。
9. **span trim/pad 极端 `-n` 溢出 panic**：`inset` 的 `lower += n` /
   `upper -= n` 对 `-n 2147483647` 等溢出（debug panic）。修复：
   `saturating_add/sub` + `clamp(NEG_INF, POS_INF - 1)`，结果截到可表示
   坐标范围。回归测试 `extreme_ops_do_not_overflow` 与 CLI 用例。
10. **excise/fill/cardinality 全幅跨度溢出 panic**：对近全幅 i32 区间
    （如 `-2147483647-2147483645`），`upper - lower + 1` 溢出 i32。
    修复：`span_len` 用 i64 计算并比较；`cardinality` 用 i64 累加并饱和到
    `i32::MAX`。
11. **gff runlist start > end 记录 panic**：`record.start > record.end` 时
    `add_pair` panic。修复：命令层报 "invalid GFF record: start X > end Y"。
12. **`pgr runlist` 裸调用 panic**：不带子命令时 `execute` 落入
    `unreachable!`（exit 101），而 fa/fas/gff 等命令组都用
    `subcommand_required`。修复：加 `subcommand_required(true)` +
    `arg_required_else_help(true)`，裸调用显示帮助并 exit 2。回归测试
    `bare_runlist_shows_help`。
13. **`is_neg_inf` / `is_pos_inf` 空集 unwrap panic**：空集上
    `edges.front()/back().unwrap()` panic，`is_infinite()` 因此对空集 panic。
    修复：`is_some_and`，空集返回 false。回归测试
    `infinity_predicates_on_empty_set`。
20. **`rg span` 的 Range 操作加减溢出 panic**：`trim` / `trim_5p` /
    `trim_3p` / `shift_5p` / `shift_3p` / `flank_5p` / `flank_3p` 对
    `start + n`、`end - n`、`end + n + 1`、`start - n - 1` 等用裸算术，
    对合法的近上限坐标（如 `chr1:2147483645-2147483645`）配极端 `-n`
    （如 2147483647）在 debug 构建直接 panic（release 静默回绕）。复现：
    `pgr rg span mx.rg -n 2147483646`。修复：全部改为 saturating 算术
    （与 `IntSpan::inset` 的既有修复一致），`check` 继续把越序结果归为
    非法 (0,0)；`shift_3p` 重写为直接计算，去掉 `-n` 取负。
21. **`rg span` pad 路径 `-number` 取负溢出**：`--op pad -n=-2147483648`
    在 `-number` 处 panic（span.rs 与 `IntSpan::pad` 各一处）。修复：
    `number.saturating_neg()`；`IntSpan::pad` 同样改 `n.saturating_neg()`。
22. **runlist 解析器 i64 累加溢出**：`runlist_to_ranges` 用 i64 累加数字，
    超过 19 位的纯数字串（如 `"99999999999999999999"`）在
    `lower * radix` 处溢出 panic，影响所有 runlist JSON 入口（`rg prop` /
    `rg runlist` / `runlist span/compare/combine/convert/stat/statop`，
    以及 `fa mask` / `fas slice` 共用的 `io::read_runlist`）。修复：每位
    累加后检查是否低于 i32::MIN（累加单调递减，提前退出安全），越界返回
    "out of range" 错误。
23. **`rg span` 合法结果超出可表示坐标域**：saturating 后可能出现
    `chr1:2147483646-2147483647` 这类"合法但超出 POS_INF - 1"的输出，
    下游 rg 命令会静默丢弃。修复：CLI 层 `clamp_to_domain` 把合法结果夹回
    `1..=POS_INF - 1`，越序则归为非法行。

### 输入校验 / 静默错误（3 处）

7. **multi 文件被当作 single 传参时静默变空集**：`statop` / `compare` 的
   infile2 为多层 JSON 时，`json_to_set` 的 `filter_map` 把对象值全部静默
   跳过，统计结果全 0 而无提示。修复：非字符串值报
   "runlist value for ... is not a string"；`json_to_sets` 的 multi 分支对
   非对象值同样报错（混合形态不再静默丢数据）。
14. **`json_to_sets` 混合形态静默丢数据**：`{"a":"1-5","b":{...}}` 等
    multi/flat 混杂输入按首值判定形态后，另一形态的值被静默丢弃。修复：
    flat 分支的非字符串值、multi 分支的非对象值均报错。
15. **删除未使用且有 panic 隐患的 `gff_to_set`**：`libs/runlist` 中该函数
    除单元测试外无调用者（`gff runlist` 命令走 `libs/fmt/gff` 的 noodles
    解析），且对 start > end 记录会 panic。直接删除函数及其测试，避免误导
    后续调用者。

### 外部工具与参数 / CLI / 文档（1 处）

19. **`span -n` 缺帮助文本**：`--number` 无 `.help()`，帮助页空白。补：
    "Number of bases to trim or pad; length threshold for excise/fill"。

### 迁移遗留 / 行为一致性（5 处）

16. **rept / p2m 测试仍以 spanr 缺失为跳过条件**：`tests/cli_rept.rs` 的
    e-align/s-kmer/e-kmer/trf/s-align 用例与 `tests/cli_pl.rs` 的 p2m 用例
    在 `spanr` 不在 PATH 时静默跳过，但管道已全部内建、不再需要 spanr。
    修复：移除 spanr 守卫（e_align 端到端测试恢复真正执行）。
17. **用户文档仍把 spanr 列为依赖**：README.md、docs/usage_examples.md、
    docs/pl.md、docs/rept.md 中的依赖清单与 stat/statop 示例仍写 spanr。
    修复：改为 `trf` / `FastK` / `Profex`（p2m 无外部依赖）与
    `pgr runlist stat/statop` 示例。
18. **src 注释残留 spanr 现态描述**：`libs/pl/repeat.rs`、`cmd_pgr/rept/
    trf.rs`、`libs/pl/mod.rs` 的注释把当前管道说成 "spanr cover 截断点号
    名"、"spanr 管道" 等。修复：改为 "runlist 解析器" / "runlist 管道"
    （历史出处说明保留）。
24. **`#` 注释行在各 rg 子命令间不一致**：`cover`/`coverage`（走
    `rg_to_set`/`rg_to_intervals`）显式跳过 `#` 行，而 `count`/`span`/
    `sort`/`prop`/`runlist`/`merge` 只靠解析失败跳过，`# chr1:1-10`
    会被当成数据（count 计入、sort 排进正文）。修复：所有 rg 命令统一
    跳过 `trim_start().starts_with('#')` 的行（与 rgr 默认行为及
    cover/coverage 一致）；`docs/rg.md` 与 `rg sort` 帮助文本补充说明；
    新增全家族回归测试。
25. **`rg merge` 自映射 parity**：当某个 part 恰好等于合并串
    （`chr1(+):min-max`）时，rgr 跳过该 `part→merged` 自映射，pgr 会
    输出无意义的映射行。修复：`rg_merge_mapping` 在
    `parts[i].0 == merged` 时跳过。

### 输入校验 / 静默错误（第 7 轮，2 处）

26. **`Range::from_str` 的 end 坐标溢出被静默折叠为 start**：
    `tail_match` 对超过 i32 的数字串返回 0（意图是"无效"），但 `decode`
    把 `end == 0` 当作"缺失"并默认成 start，导致
    `chr1:5-99999999999` 被解析为合法的点区间 `chr1:5`——`rg count /
    span / sort / prop / runlist` 全部基于错误坐标计算而输出原行，
    静默产生错误结果（`rg span` 还会把错误结果写进输出）。修复：
    `parse_i32` 改为返回 `Option<i32>`（溢出返回 `None`），
    `tail_match` 传播 `None` 使整行不匹配（落入无效行分支）；字面量
    `end = 0`、缺失 end、裸 start 的默认行为不变。回归测试
    `overflow_end_is_invalid_not_start` 与 CLI 用例
    `command_rg_overflow_end_skipped`。
27. **`IntSpan::inset` 的 clamp 下界把可表示坐标 i32::MIN 静默改写**：
    `inset` 把饱和结果 `clamp(NEG_INF, POS_INF - 1)`，而解析器接受
    `-2147483648`（`test_valid` 显式断言），于是 `runlist span --op
    trim/pad -n 0`（恒等操作）会把 `{-2147483648}` 变成
    `{-2147483647}`。修复：clamp 下界改为 `i32::MIN`（与解析器域一致），
    只保留上界 `POS_INF - 1` 的截断。回归测试
    `inset_identity_at_i32_min`。

### 参数校验 / 数据安全 / 错误传播（第 10 轮，3 处）

28. **`rg span` 的 op/mode 组合校验推迟到逐行处理**：`--op shift/flank
    -m both` 的报错发生在处理到第一条有效行时；输入为空或全注释/无效行
    时命令静默成功（exit 0），有有效行时才失败（原版 rgr 这里直接
    `unreachable!` panic）。修复：读取任何输入前先校验
    `matches!(op, "shift" | "flank") && mode == "both"` 并提前报错，
    逐行 match 的兜底分支改为 `unreachable!`。回归测试
    `command_rg_span_invalid_mode_checked_before_input`（empty 与
    comments-only 输入 × shift/flank 两条命令）。
29. **流式命令 `-o` 与输入同路径时先截断输入文件，静默数据丢失**：
    `rg span` / `rg prop` / `rg runlist` 在读取输入行之前创建 writer，
    `rg count` 建立索引后创建 writer 再读 target，`runlist convert` 在
    读取 JSON 之前创建 writer——`-o 同输入` 会把输入清空后继续执行
    （`rg sort` 因先缓冲不受影响；先读后写的 cover/coverage/merge 及
    runlist 其余子命令天然支持 in-place）。修复：`libs/io` 新增
    `same_path`（`std::path::absolute` 词法归一化后比较，无文件系统访问），
    五个命令在读取前对所有位置参数做同路径检查并报
    "output file ... is also an input file"。回归测试
    `command_rg_output_same_as_input_rejected`（span/prop/runlist/count
    四例 + 断言输入文件未被改动）与
    `command_runlist_convert_output_same_as_input_rejected`。
30. **`IntSpan::write_json` 未显式 flush**：写盘失败时错误被
    `PgrWriter::drop` 降级为 stderr warning 且进程以 0 退出，与其余命令
    显式 `writer.flush()?` 传播错误的行为不一致（影响 `rg cover` /
    `rg coverage` / `runlist span` / `combine` / `genome` / `merge` /
    `some` / `gff runlist` 的 JSON 输出路径）。修复：`write_json` 末尾
    增加 `writer.flush()?`。

文档小修：`rg span -n` 的帮助文本补充 excise 阈值说明（"length threshold
for excise"，与 `runlist span -n` 一致）。

### 性能（第 11 轮，1 处）

31. **`rg merge` 去重 O(n²) 退化**：`rg_merge_mapping` 用
    `Vec::contains` 逐条去重，单染色体区间数 n 的去重为 O(n²)；
    100k 条 disjoint 区间（debug 构建）耗时 37.8 s。修复：并行维护
    `HashSet<(line, start, end)>` 做 O(1) 判重，`parts` 仍按输入顺序
    保留以维持稳定输出；同一输入修复后 0.28 s。新增基准
    `benches/rg_merge_benchmark.rs`（release：10k ≈ 2.7 ms、50k ≈
    18 ms，disjoint 与 clustered 两种形态），新增回归测试
    `command_rg_merge_dedups_identical_lines`（重复行只保留一份、
    不产生自簇，第三条重叠区间正常入簇）。

## 与外部参考实现的语义一致性核对

逐条对照原始 spanr 0.6.7 源码（docs.rs）确认 runlist 家族未改语义：

* `statop` 的 `c2 = s2_size / s2_length`、`ratio = c2 / c1`、`all` 行汇总
  公式逐字段一致（含 f32/f64 格式化的差异保留）。
* `stat` 的逐染色体行 + `all` 行、`--all` 只保留全基因组行一致。
* `compare` 的 `chrs` 取并集、缺失染色体补空、多 others 顺序折叠一致。
* `merge` 的 stem 命名（`--all` = 完整 stem）一致。
* `coverage` 的 `[start, end+1)` 半开区间 + 深度过滤语义一致；扫描线
  实现与 rust-lapper `depth()` 输出 diff 验证一致（基准见
  `notes/benchmarks/interval-overlap.md`）。
* `span` 的 trim/pad/excise/fill 与迁入的 intspan crate 语义一致。

逐条对照原始 rgr（intspan 项目 `cmd_rgr`）确认 rg 家族未改语义：

* `rg count` 的 inclusive 端点判交与 rgr 的 `[start, end+1)` 半开转换等价。
* `rg sort` 的排序键 (chr, start, strand) 与 rgr `sort_by_cached_key`
  一致（rgr 用 BTreeMap 去重，pgr 保留重复行并稳定排序，为有意差异）。
* `rg runlist` 的 superset = "runlist ⊇ range" 与 rgr
  `set[chr].superset(&intspan)` 一致。
* `rg merge` 的 reciprocal coverage 判据与 `part == merged` 自映射跳过
  与 rgr 一致（rgr 用 petgraph 连通分量，pgr 用 coitrees + union-find，
  聚类结果等价）。

## 排除的疑点（经核验无需修复）

* **行首/尾空白与后缀文本**：`Range::from_str` 的非锚定最左匹配天然容忍
  前导/尾随空白与 `|Species=Yeast` 等后缀，各 rg 命令行为一致，无需额外
  trim。
* **coitrees 输入顺序**：`BasicCOITree::new` 内部会按 (start, end) 排序，
  `RgIndex` / `rg merge` 按行序喂入不排序的区间不会错乱。
* **stat / statop 的 `--all` 表头字段数**：multi 与 single 两种形态下
  header 与数据行列数均一致（含 `--all` 去掉 `chr,` / `all,` 后），
  spanr parity 测试通过。
* **`rg span` 饱和结果的取舍**：Range 层保留饱和语义（与 vendored crate
  对齐），CLI 层 `clamp_to_domain` 保证输出可回读（见 #23），两处职责
  分开，不在 Range 层重复夹紧。
* **`rg span` pad 使 start 越过 1 时输出空行**：`Range::check` 把负坐标
  夹到 0（0 即"无效"哨兵），pad 处于坐标 1 的区间会得到无效行——与
  vendored crate / rgr 行为一致，未改（改动会破坏 spanr parity）。
* **`Range::from_str` 对 `1:-100` 等以 `-` 开头坐标的 fallback**：`:`
  后不是数字时不匹配，回退为整行作为 chr（`fa_headers` 固定语义），
  修复 #26 不影响该路径。

## 记录项（未改，低风险 / 待决策）

* `__single__` 哨兵键与"恰好命名为 `__single__` 的单集合"理论上碰撞，
  与原始 spanr 的 `__single` 设计一致，未处理。
* `runlist split` 的 `<key><suffix>` 文件名直接由 JSON 顶层键拼出：键含
  `/` 或 `..` 时文件会写到 outdir 之外（或嵌套目录）。与 spanr 行为一致，
  本地工具由用户自持输入，未加防御。

## 已知限制（有意保留）

* `IntSpan::from` / `add_runlist` 保持外部 crate 的 panic API（API 兼容），
  所有命令行入口已前置 `IntSpan::valid` 校验，CLI 不再可达。
* `stat` / `statop` 对 0 长度染色体输出 inf/NaN，与 spanr 行为一致。
* `banish` / `distance` / `at` / `index` / `slice` / `contains` 等 vendored
  API 的极端输入溢出未修（runlist/rg 命令不经过这些路径）。

## 带点 contig 名截断 bug 的处置结论（对应 repeat-masking.md 记录）

`notes/design/repeat-masking.md` 记录了 spanr 时代"带点 contig 名被
`chr:start-end` 解析按 `.` 截断"的 bug（`NC_000913.1` → `"1"`，多 contig
时 `chr1.1` / `chr2.1` key 冲突），当时以 `c1..cN` 名字映射规避并验证。
解析器迁入本地后重新评估：**不修改解析语义**，维持 `c1..cN` 映射。

理由：

1. `name.chr(strand):start-end` 是 `.rg` 格式的物种前缀约定，兼容测试
   明确断言 `S288c.I(-):190-200` → 染色体 `I` 且输出不含 `S288c`
   （`command_cover` / `command_coverage`）；`docs/runlist.md` 承诺与
   外部 spanr 输出一致。全局去掉 `.` 拆分是破坏性变更。
2. `Range::from_str` 的 `name.chr` 拆分同时服务 FASTA 头解析、`psl
   to-range` 等消费者，`range.rs` 的差分测试把 `a.b.c:1-2` → name `b`、
   chr `c` 作为固定语义（与正则解析器对拍）。
3. bug 记录中的实际影响（四个 rept 管线的 runlist key 被截断导致
   `fa mask` 失配）已由 `c1..cN` 映射完全修复并验证（`NC_000913.1` 与
   Fusarium 44 条 scaffold 均完整），且映射代码在管道内、成本极低。

若未来希望带点 contig 名在 `pgr rg cover/coverage` 输入中原生可用，
正解是新增 rg 专用的严格解析模式（或 `--no-species-prefix` 开关），
属于新特性，需与 spanr 兼容目标权衡，不建议按 bug 修复处理。

## 验证

### 交叉 / 差分验证

* `rg count` / `rg coverage` / `rg prop` 与暴力实现随机对拍 80 trial 一致；
  `rg merge` 的合并覆盖区间与簇成员并集 300 trial 一致（含
  `part == merged` 自映射跳过的解释）。
* 与 spanr / rgr 的语义一致性核对见上（输出逐字段对照 + fixture 兼容
  测试）。

### 鲁棒性（Zero Panic）

* 第 1–5 轮：1600+ 条随机畸形输入（JSON / .rg / GFF / sizes / 二进制）对
  全部 runlist 命令及 `fa mask` / `gff runlist` fuzz，无 panic。
* 第 6 轮：rg 家族两轮 fuzz 共 1000+ trial（随机畸形行 + 近上限坐标 +
  极端 `-n`，每 trial 覆盖 16~19 条命令调用）、runlist 家族 400 trial
  （随机 JSON / 尺寸文件 + 极端 `-n`）、收尾轮 500 trial 新种子——零 panic。
* 第 7 轮：暴力差分对拍（set ops 3000 trial、depth 3000 trial、
  covered 3000 trial、rg merge 120 trial、Range 极端运算、CLI 近上限
  坐标 200 trial）全一致；300 trial × ~20 条命令的随机畸形输入 fuzz
  零 panic。修复 #26/#27 后全部单元/CLI 测试与 52 个测试目标通过。
* 第 10 轮：修复 #28-#30 后重跑四组差分对拍（250 + 200 + 300 + 150
  trial，全部一致）与两轮全新种子 fuzz（各 120 trial，含极端 `-n`、
  畸形 JSON / 尺寸文件 / `.rg` 文本），零 panic；全量测试 + doctest
  通过，`cargo fmt`、`cargo clippy --all-targets -- -D warnings` 干净。
* 第 11 轮：`rg merge` 100k 区间从 37.8 s 降至 0.28 s（debug 实测）；
  criterion 基准（release）10k ≈ 2.7–2.9 ms、50k ≈ 17.8–18.4 ms。
* 第 12–15 轮：stat/statop 差分 150 trial、multi compare 120 trial、
  Range 运算差分 700 trial 全部一致；docs 示例逐条执行通过；全量测试
  + doctest 通过；`cargo fmt`、`cargo clippy --all-targets -- -D
  warnings` 干净。

### 最终状态

* `cargo fmt`、`cargo clippy --all-targets -- -D warnings` 干净；全部测试
  二进制 + doctest 通过（1178 个断言）。
* 新增回归测试：runlist 阶段 14 个函数（含 test_valid / genome / cover 等
  既有用例扩展）；rg 阶段 `command_rg_span_extreme_no_panic`（8 种 op/mode
  × 极端 `-n` + `pad i32::MIN`）、`command_rg_comments_skipped`（7 个子
  命令）、`extreme_ops_do_not_overflow`（Range 层）、
  `command_runlist_span_invalid_runlist_errors` 增补 20 位数字串、
  `command_runlist_span_extreme_ops_do_not_panic` 增补 `pad i32::MIN`。
  第 7 轮增补：`overflow_end_is_invalid_not_start`（Range 层）、
  `inset_identity_at_i32_min`（IntSpan 层）、
  `command_rg_overflow_end_skipped`（cover/span/count/sort 四条命令）。

## 提交状态

* 第 1 阶段的大部分改动已在提交 `a20d2dc`，其余（裸命令帮助、`span -n`
  帮助、`gff runlist` 坐标上限、`.rg` max_coord 守卫及回归测试）已在
  `4dd6f8c` 落地。
* 第 2 阶段（rg 家族深审）的代码与测试改动已在 `aac5c50` 提交；本文档的
  合并、重命名与结构调整（与 `audit-sd-rept-align.md` 统一骨架）暂存待
  提交。
* 第 7 轮的 #26/#27 修复与回归测试、第 8 轮核验记录在工作区（本环境
  `.git` 只读，无法创建提交；用户侧可自行 `git commit`）。
