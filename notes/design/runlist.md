# runlist 命令族（spanr 迁移）

## 背景

外部 `spanr`（intspan repo 的 CLI）是 pgr 重复/比对管线的区间运算依赖：
`rept e-kmer/s-kmer/e-align/s-align/trf` 与 `pl p2m` 通过 PATH 调用
`spanr cover / coverage / span / compare / merge`。为减少外部依赖，迁为
pgr 内部命令族 `pgr runlist`（命名避免与旧 spanr 冲突，符合 pgr 格式优先
命名习惯）。

## 结构

* `libs/runlist/`：核心逻辑（rg 解析、深度扫描线、span/compare/merge、
  JSON 读写），复用 `libs/ds::IntSpan` 与已迁入的 `set2json` 等辅助。
* `cmd_pgr/runlist/`：五个薄壳子命令（cover/coverage/span/compare/merge），
  与 spanr CLI 参数兼容（含 `--detailed`、`--op`、`--all` 等）。
* 管线改为进程内调用：`repeat.rs`（cover→fill→excise→fill、coverage）、
  `p2m.rs`（compare→span excise→merge）、`trf.rs`（cover）。

## coverage 实现与性能

`coverage` 用**扫描线**（按 start/end 事件排序后单遍累计深度），
O(n log n)，不需要区间树——纯深度聚合场景下这是标准最优做法（coitrees
适合重叠枚举而非逐点深度）。详细模式（`-d`）在同一遍扫描里按精确深度
分桶输出。

1M 随机区间（chr 100 Mb）实测：

| 实现 | 耗时 |
| --- | ---: |
| pgr runlist coverage（扫描线） | 2.5 s |
| 外部 spanr coverage（lapper depth） | 26.8 s |

输出逐字节一致（100k 与 1M 两组都 diff 通过）。

受控对比（release，1M 区间、仅深度计算、不含解析）：

| 实现 | 耗时 |
| --- | ---: |
| pgr 扫描线 `depth_at_least` | 0.20 s |
| rust-lapper `Lapper::new` + `depth()` | 22.8 s（build 0.08 s + depth 22.8 s） |

覆盖碱基数完全一致（100,046,409）。rust-lapper 确实"适合"覆盖度场景
（spanr 就是用它的 `depth()`），但实测其 `depth()` 在密集区间上比扫描线
慢约两个数量级；扫描线实现更简单且无额外依赖。`coverage depth` 基准组
（benches/interval_overlap_benchmark.rs）在 10k/100k 上保留对比，release
中位数：

| n | pgr 扫描线 | lapper `depth()` |
| --- | ---: | ---: |
| 10k | 1.3 ms | 469 ms |
| 100k | 15.6 ms | 3.68 s |

## 兼容性

* 与 spanr 相同的语义：`.rg` 行为 1-based inclusive；`Range` 解析的
  非锚定最左匹配/回退行为一致；dotted 名字仍由管线预先映射为安全名
  （`c1..cN`）再还原。
* 外部 spanr 仍保留（anchr 等其它项目在用）；pgr 不再依赖它。
