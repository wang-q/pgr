# IntSpan 集合运算与批量构建基准（criterion）

> 目的：把 `intersect`/`union`/`diff`/`xor` 线性化与 `from_pairs` 批量构建
> 的收益固化为库级 criterion 基准，直接对比新旧实现（旧实现用公开 API
> 在 bench 内重建）。2026-08-04 实测；同日复测（修复轮后）数值基本一致。

## 基准文件

`benches/intspan_setops_benchmark.rs`：

* `setops`：5k / 20k span 的两条 runlist（100 Mb 染色体，长度 100–2000），
  四个集合运算 × 新旧实现。
* `construction`：10k / 100k 个随机 1 bp span（互不重叠的 adversarial
  输入），`from_pairs`（排序构建）vs `add_pair` 循环（逐个插入）。
* `covered`：2k / 5k span 的 runlist + 2,000 条查询，比较
  `IntSpan::covered`（二分重叠段）、`intersect+cardinality`（线性基线）、
  以及"把 span 抽成 Vec 后用 `partition_point`"（曾用 SpanIndex 的做法，
  已合并回 IntSpan）。

```bash
cargo bench --bench intspan_setops_benchmark
```

## 结果（median，release，measurement-time 2s；2026-08-04 复测）

### setops

| op | n=5k 旧 | n=5k 新 | 加速 | n=20k 旧 | n=20k 新 | 加速 |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| intersect | 673.1 µs | 30.6 µs | **22×** | 11.42 ms | 102.7 µs | **111×** |
| union | 1.147 ms | 27.5 µs | **42×** | 15.26 ms | 96.4 µs* | **158×** |
| diff | 1.548 ms | 28.9 µs | **54×** | 17.01 ms | 80.6 µs | **211×** |
| xor | 2.084 ms | 84.8 µs | **25×** | 32.86 ms | 444.9 µs | **74×** |

\* union 20k 整轮跑时为 128.1 µs（受前序重基准的热噪声影响），单独复测
为 96.4 µs（91.1–101.8 µs），与初测 87.9 µs 同量级，取复测值。

加速随 n 放大（O(n·m) → O(n+m)）；20k 时 74–221×，远高于 CLI 基准的
5–8×（CLI 时间被 runlist JSON 加载占据）。

### construction（随机 1 bp 稀疏区间，adversarial）

| n | `from_pairs` | `add_pair` 循环 | 加速 |
| :--- | ---: | ---: | ---: |
| 10k | 125.0 µs | 2.241 ms | **18×** |
| 100k | 1.738 ms | 185.9 ms | **107×** |

`add_pair` 逐个插入无序区间是 O(n²)（VecDeque 中间搬移），`from_pairs`
排序 + 单遍合并为 O(n log n)，差距随 n 放大。

### covered（2,000 条查询，2000 / 5000 span）

| 实现 | 2000 span | 5000 span |
| :--- | ---: | ---: |
| `IntSpan::covered`（as_slices 快路径） | 34.0 µs | 39.2 µs |
| `partition_point` on Vec（SpanIndex 式） | 27.0 µs | 35.2 µs |
| `intersect` + `cardinality`（线性基线） | 8.89 ms | 22.25 ms |

covered 每查询 ~13–19 ns，比线性 intersect 快两个数量级，比"抽成 Vec 再
partition_point"慢 ~1.3×（VecDeque 索引/闭包开销）；代价是省掉了 SpanIndex
那份重复的 span 数据。`covered()` 用 `as_slices()` 快路径（append 构建的
VecDeque 是单切片）后较初版（37.9/42.1 µs）提升 ~11%，剩余常数差距可接受，
不再追。

## 结论

库级基准证实：集合运算线性化（20k span 时 ~100–220×）与批量构建
（100k 时 ~105×）都是数量级收益，且随规模增长；这与 CLI 基准
（bench-rg-prop.md、../design/runlist.md）互补——CLI 看端到端，criterion
看纯库层。旧的 O(n·m)/O(n²) 实现保留在 bench 内作基线，便于后续回归
监控。

## 附：模块第二轮系统性审视（2026-08-04）

针对"IntSpan 是大量工作的基础"再做一轮盘点：

* **修复 `banish` 溢出**：`end - start + 1` 与平移坐标在极端参数
  （如 `banish(i32::MIN, i32::MAX)`）下溢出 panic。改 i64 计算 +
  clamp，`start > end` 返回原集合（无可 ban 内容）。回归测试
  `banish_extreme_args_do_not_overflow`。
* **核实安全的点**（无需改）：
  * `list_to_ranges` 的 `vec[end-1] + 1`：`vec[end-1] = i32::MAX` 时
    `end < len` 不可能成立（去重后无更大元素），短路保证不溢出。
  * `to_vec`/`elements`：仓库内无调用者（公共 API，vendored 语义）。
  * `at`/`index`/`slice`：O(n) 扫描，仅 alignment 小集合使用。
  * `merge`/`subtract`/`add_pair` 增量 O(n) 搬移：大集合路径已全部走
    线性/排序构建（union/diff/from_pairs/rg_files_to_set），剩余调用点
    小规模或有序追加。
  * edges 的 `±1` 运算都在坐标契约（≤ POS_INF-1）内，无越界。
* **基准覆盖**：setops / construction / covered 三组 criterion 已就位，
  未来结构改动先量化再决定。
