# rg 命令族代码审核记录（2026-08-04）

对新增的 `pgr rg` 命令族（cover/coverage/count/merge/prop/runlist/sort/span）
做多轮深入审核，并复核 runlist 家族与共享底层（`libs/runlist`、
`libs/ds/intspan`、`libs/ds/range`）。上一轮审核记录见
`notes/audit/audit-runlist.md`（runlist 家族 5 轮、19 处修复）；本轮在其
基础上对 rg 家族做全新深审，并对共享库新增路径重新 fuzz。每轮发现问题后
修复并进入下一轮复核；最后一轮全量重读 + 三组 fuzz + 差分对拍未再发现
新问题后收束。

## 修复的缺陷（6 处）

### 崩溃 / 溢出（4 处，Zero Panic）

1. **`rg span` 的 Range 操作加减溢出 panic**：`trim` / `trim_5p` /
   `trim_3p` / `shift_5p` / `shift_3p` / `flank_5p` / `flank_3p` 对
   `start + n`、`end - n`、`end + n + 1`、`start - n - 1` 等用裸算术，
   对合法的近上限坐标（如 `chr1:2147483645-2147483645`）配极端 `-n`
   （如 2147483647）在 debug 构建直接 panic（release 静默回绕）。复现：
   `pgr rg span mx.rg -n 2147483646`。修复：全部改为 saturating 算术
   （与 `IntSpan::inset` 的既有修复一致），`check` 继续把越序结果归为
   非法 (0,0)；`shift_3p` 重写为直接计算，去掉 `-n` 取负。
2. **`rg span` pad 路径 `-number` 取负溢出**：`--op pad -n=-2147483648`
   在 `-number` 处 panic（span.rs 与 `IntSpan::pad` 各一处）。修复：
   `number.saturating_neg()`；`IntSpan::pad` 同样改 `n.saturating_neg()`。
3. **runlist 解析器 i64 累加溢出**：`runlist_to_ranges` 用 i64 累加数字，
   超过 19 位的纯数字串（如 `"99999999999999999999"`）在
   `lower * radix` 处溢出 panic，影响所有 runlist JSON 入口（`rg prop` /
   `rg runlist` / `runlist span/compare/combine/convert/stat/statop`，
   以及 `fa mask` / `fas slice` 共用的 `io::read_runlist`）。修复：每位
   累加后检查是否低于 i32::MIN（累加单调递减，提前退出安全），越界返回
   "out of range" 错误。
4. **`rg span` 合法结果超出可表示坐标域**：saturating 后可能出现
   `chr1:2147483646-2147483647` 这类"合法但超出 POS_INF - 1"的输出，
   下游 rg 命令会静默丢弃。修复：CLI 层 `clamp_to_domain` 把合法结果夹回
   `1..=POS_INF - 1`，越序则归为非法行。

### 行为一致性 / parity（2 处）

5. **`#` 注释行在各 rg 子命令间不一致**：`cover`/`coverage`（走
   `rg_to_set`/`rg_to_intervals`）显式跳过 `#` 行，而 `count`/`span`/
   `sort`/`prop`/`runlist`/`merge` 只靠解析失败跳过，`# chr1:1-10`
   会被当成数据（count 计入、sort 排进正文）。修复：所有 rg 命令统一
   跳过 `trim_start().starts_with('#')` 的行（与 rgr 默认行为及
   cover/coverage 一致）；`docs/rg.md` 与 `rg sort` 帮助文本补充说明；
   新增全家族回归测试。
6. **`rg merge` 自映射 parity**：当某个 part 恰好等于合并串
   （`chr1(+):min-max`）时，rgr 跳过该 `part→merged` 自映射，pgr 会
   输出无意义的映射行。修复：`rg_merge_mapping` 在
   `parts[i].0 == merged` 时跳过。

## 验证

* `cargo fmt`、`cargo clippy --all-targets -- -D warnings` 干净；
  全部 52 个测试二进制 + doctest 通过（共 1178 个断言）。
* 新增回归测试：`command_rg_span_extreme_no_panic`（8 种 op/mode ×
  极端 `-n` + `pad i32::MIN`）、`command_rg_comments_skipped`（7 个子
  命令）、`extreme_ops_do_not_overflow`（Range 层）、
  `command_runlist_span_invalid_runlist_errors` 增补 20 位数字串、
  `command_runlist_span_extreme_ops_do_not_panic` 增补 `pad i32::MIN`。
* fuzz：rg 家族两轮共 1000+ trial（随机畸形行 + 近上限坐标 + 极端
  `-n`，每 trial 覆盖 16~19 条命令调用）、runlist 家族 400 trial
  （随机 JSON/尺寸文件 + 极端 `-n`）——零 panic。
* 差分对拍：`rg count` / `rg coverage` / `rg prop` 与暴力实现随机对拍
  80 trial 一致；`rg merge` 的合并覆盖区间与簇成员并集 300 trial 一致。

## 已知限制（沿用上一轮，未变）

* `IntSpan::from` / `add_runlist` 保持外部 crate 的 panic API（API 兼容），
  所有命令行入口已前置校验，CLI 不可达。
* `__single__` 哨兵键与"恰好命名为 `__single__` 的单集合"理论上碰撞，
  与原始 spanr 设计一致。
* `stat` / `statop` 对 0 长度染色体输出 inf/NaN，与 spanr 一致。
* `banish` / `distance` / `at` / `index` / `slice` 等 vendored API 的
  极端输入溢出未修（runlist/rg 命令不经过这些路径）。
* `Range::contains(i32::MAX)` 等 vendored API 的 `n + 1` 溢出未修
  （新命令不调用）。
