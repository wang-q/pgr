# runlist 命令族代码审核记录（2026-08-04）

对新增的 `pgr runlist` 命令族及相关的库文件进行多轮深入审核。范围：
`cmd_pgr/runlist/` 12 个子命令（combine/compare/convert/cover/coverage/
genome/merge/some/span/split/stat/statop）、`cmd_pgr/gff/runlist`、
`libs/runlist`、迁入的 `libs/ds/intspan` 与手写扫描器版 `libs/ds/range`、
`libs/fmt/gff`、`libs/io::read_runlist`，以及迁移了 spanr 调用的
`libs/pl/repeat`、`cmd_pgr/pl/p2m`、`cmd_pgr/rept/trf` 和全部测试/文档。
每轮发现问题后修复并进入下一轮复核；最后一轮（第 5 轮）经全量重读、
1600+ 条随机畸形输入 fuzz 与边界探测未再发现新问题后收束。
最终 fmt / clippy（--all-targets）干净，51 个测试二进制 + 71 个 doctest
全绿。

## 修复的缺陷（19 处）

### 崩溃 / 溢出（13 处，Zero Panic）

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
7. **multi 文件被当作 single 传参时静默变空集**：`statop` / `compare` 的
   infile2 为多层 JSON 时，`json_to_set` 的 `filter_map` 把对象值全部静默
   跳过，统计结果全 0 而无提示。修复：非字符串值报
   "runlist value for ... is not a string"；`json_to_sets` 的 multi 分支对
   非对象值同样报错（混合形态不再静默丢数据）。
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

### 输入校验 / 静默错误（2 处）

14. **`json_to_sets` 混合形态静默丢数据**：`{"a":"1-5","b":{...}}` 等
    multi/flat 混杂输入按首值判定形态后，另一形态的值被静默丢弃。修复：
    flat 分支的非字符串值、multi 分支的非对象值均报错。
15. **删除未使用且有 panic 隐患的 `gff_to_set`**：`libs/runlist` 中该函数
    除单元测试外无调用者（`gff runlist` 命令走 `libs/fmt/gff` 的 noodles
    解析），且对 start > end 记录会 panic。直接删除函数及其测试，避免误导
    后续调用者。

### 迁移遗留（4 处）

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
19. **`span -n` 缺帮助文本**：`--number` 无 `.help()`，帮助页空白。补：
    "Number of bases to trim or pad; length threshold for excise/fill"。

## 与 spanr / intspan 的语义一致性核对

逐条对照原始 spanr 0.6.7 源码（docs.rs）确认未改语义：

* `statop` 的 `c2 = s2_size / s2_length`、`ratio = c2 / c1`、`all` 行汇总
  公式逐字段一致（含 f32/f64 格式化的差异保留）。
* `stat` 的逐染色体行 + `all` 行、`--all` 只保留全基因组行一致。
* `compare` 的 `chrs` 取并集、缺失染色体补空、多 others 顺序折叠一致。
* `merge` 的 stem 命名（`--all` = 完整 stem）一致。
* `coverage` 的 `[start, end+1)` 半开区间 + 深度过滤语义一致；扫描线
  实现与 rust-lapper `depth()` 输出 diff 验证一致（基准见
  `notes/benchmarks/interval-overlap.md`）。
* `span` 的 trim/pad/excise/fill 与迁入的 intspan crate 语义一致。

## 验证

* `cargo fmt`、`cargo clippy --all-targets -- -D warnings` 干净。
* 51 个测试二进制全部通过（含 71 个 doctest）；新增回归测试 14 个函数
  （另有对 test_valid / genome / cover 等既有用例的扩展），覆盖反转区间、
  非法 runlist、越界坐标、极端 span 参数、multi 误用、裸命令、GFF 畸形
  记录等。
* 1600+ 条随机畸形输入（JSON / .rg / GFF / sizes / 二进制）对全部
  runlist 命令及 `fa mask` / `gff runlist` fuzz，无 panic。

## 已知限制（有意保留）

* `IntSpan::from` / `add_runlist` 保持外部 crate 的 panic API（API 兼容），
  所有命令行入口已前置 `IntSpan::valid` 校验，CLI 不再可达。
* `__single__` 哨兵键与"恰好命名为 `__single__` 的单集合"理论上碰撞，
  与原始 spanr 的 `__single` 设计一致，未处理。
* `stat` / `statop` 对 0 长度染色体输出 inf/NaN，与 spanr 行为一致。
* `banish` / `distance` / `at` / `index` / `slice` 等 vendored API 的极端
  输入溢出未修（runlist 命令不经过这些路径）。

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

## 提交状态

大部分改动（解析器含坐标上限拒绝、JSON 校验、gff 反转记录、inset/excise/
fill/cardinality 溢出、文档与测试守卫）已在提交 `a20d2dc`。工作区仍有未
提交改动（`.git` 目录只读，待环境或用户提交）：`pgr runlist` 裸命令修复
（`subcommand_required`）、`span -n` 帮助、`gff runlist` 坐标上限收紧到
`POS_INF - 1`、`.rg` 入口的 max_coord 守卫，以及对应回归测试。
