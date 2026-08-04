# runlist / rg 命令族代码审核记录（2026-08-04）

对新增的 `pgr runlist` 与 `pgr rg` 两个命令族及相关的库文件进行多轮深入
审核。范围：`cmd_pgr/runlist/` 12 个子命令（combine/compare/convert/cover/
coverage/genome/merge/some/span/split/stat/statop）、`cmd_pgr/rg/` 8 个子
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

## 修复的缺陷（共 25 处）

修复按发现顺序全局编号（#1-19 为 runlist 阶段、#20-25 为 rg 阶段），
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

## 记录项（未改，低风险 / 待决策）

* `__single__` 哨兵键与"恰好命名为 `__single__` 的单集合"理论上碰撞，
  与原始 spanr 的 `__single` 设计一致，未处理。

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

### 最终状态

* `cargo fmt`、`cargo clippy --all-targets -- -D warnings` 干净；全部测试
  二进制 + doctest 通过（1178 个断言）。
* 新增回归测试：runlist 阶段 14 个函数（含 test_valid / genome / cover 等
  既有用例扩展）；rg 阶段 `command_rg_span_extreme_no_panic`（8 种 op/mode
  × 极端 `-n` + `pad i32::MIN`）、`command_rg_comments_skipped`（7 个子
  命令）、`extreme_ops_do_not_overflow`（Range 层）、
  `command_runlist_span_invalid_runlist_errors` 增补 20 位数字串、
  `command_runlist_span_extreme_ops_do_not_panic` 增补 `pad i32::MIN`。

## 提交状态

* 第 1 阶段的大部分改动已在提交 `a20d2dc`，其余（裸命令帮助、`span -n`
  帮助、`gff runlist` 坐标上限、`.rg` max_coord 守卫及回归测试）已在
  `4dd6f8c` 落地。
* 第 2 阶段（rg 家族深审）的代码与测试改动已在 `aac5c50` 提交；本文档的
  合并、重命名与结构调整（与 `audit-sd-rept-align.md` 统一骨架）暂存待
  提交。
