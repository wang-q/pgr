# runlist 命令族（spanr 迁移）

## 背景

外部 `spanr`（intspan repo 的 CLI）是 pgr 重复/比对管线的区间运算依赖：
`rept e-kmer/s-kmer/e-align/s-align/trf` 与 `pl p2m` 通过 PATH 调用
`spanr cover / coverage / span / compare / merge`。为减少外部依赖，迁为
pgr 内部命令族 `pgr runlist`（命名避免与旧 spanr 冲突，符合 pgr 格式优先
命名习惯）。

## 结构

* `libs/runlist/`：核心逻辑（rg 解析、深度扫描线、span/compare/merge、
  combine/convert/genome/gff/some/split/stat/statop、JSON 读写），复用
  `libs/ds::IntSpan` 与已迁入的 `set2json` 等辅助。
* `cmd_pgr/runlist/`：**10 个子命令**（combine/compare/convert/genome/
  merge/some/span/split/stat/statop），与 spanr CLI 参数兼容（含 `--op`、
  `--all`、`--longest` 等），逐条与外部 spanr 输出 diff 验证一致
  （stat/statop 输出分隔符为有意差异：pgr 用 TSV、spanr 用 CSV）。
  原 spanr 的 `cover`/`coverage`（.rg 输入）迁出到 `pgr rg` 家族，`gff`
  子命令归位到 `pgr gff runlist`（GFF 输入转换归 GFF 命令管），参数与
  行为不变。
* 管线改为进程内调用：`repeat.rs`（cover→fill→excise→fill、coverage）、
  `p2m.rs`（compare→span excise→merge）、`trf.rs`（cover）。

统计类命令（stat/statop）的差异：spanr 在染色体缺失于 sizes 或输入非法
（如 statop 的 infile2 为多层 JSON）时直接 panic；pgr 改为友好报错或按
空集处理（statop 的 `s2`/`set_op` 缺失染色体按 0 计，Zero-Panic）。
输出用 TSV（tab 分隔）而非 spanr 的 CSV，字段顺序与数值公式一致。

## 测试迁移

外部 intspan 的 spanr 测试套件（`tests/cli_spanr.rs`，23 个用例）与 16 个
真实夹具（`tests/spanr/`）已迁移：夹具复制到 `tests/runlist/`，用例改写为
`pgr runlist` 调用（`tests/cli_runlist_compat.rs`），全部一次通过，覆盖
全部子命令及各 op/flag 组合。gff 相关两个用例（含 --tag 与 merge 工作流）
随命令移入 `tests/cli_gff.rs`，夹具移到 `tests/gff/`。另有自造的
`tests/cli_runlist.rs`（13 个用例）与 `libs/runlist` 单元测试兜底边界行为。

## coverage 实现与性能

`coverage` 用**扫描线**（按 start/end 事件排序后单遍累计深度），
O(n log n)，不需要区间树——纯深度聚合场景下这是标准最优做法（coitrees
适合重叠枚举而非逐点深度）。详细模式（`-d`）在同一遍扫描里按精确深度
分桶输出。

1M 随机区间（chr 100 Mb）实测：

| 实现 | 耗时 |
| --- | ---: |
| pgr rg coverage（扫描线） | 2.5 s |
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

---

## 证据附录：区间结构三方案基准：IntSpan vs coitrees vs rust-lapper

> 目的：回答"判断某个点是否落在基因组区间内"以及"枚举与窗口重叠的区间"
> 时，三种候选结构的取舍。数据生成镜像 rust-lapper 官方基准
> （`rust-lapper-master/benches/lapper_benchmark.rs`）：100 Mb 染色体、
> 区间长度 500..80 kb，另加一条覆盖 90% 染色体的超长区间作为病态用例。
> 2026-08-04 复测，数值与初测基本一致（±5%）。

## 三种结构（语义不同）

* `libs::ds::IntSpan`：合并后的 runlist 集合——只回答"点是否被覆盖"
  （`contains` 二分），丢弃区间身份；构造即做区间合并。
* `coitrees::BasicCOITree`：COITree 区间索引（PAF/pbit 已在用），枚举所有
  重叠区间，最坏情况有保证。
* `rust-lapper`：按起点排序 + 二分 + 最长区间补偿；`find`/`count`（BITS），
  外部 intspan 项目的 spanr coverage / rgr count 在用。

## 环境与执行

- `criterion 0.5.1`，`--measurement-time 2 --warm-up-time 1`，样本 100
- `rust-lapper 1.3.0`（dev-dependency，仅基准用）、`coitrees 0.4.0`
- 随机数据固定种子（`benches/interval_overlap_benchmark.rs` 内 SEED）

```bash
cargo bench --offline --bench interval_overlap_benchmark
```

## 复测结果（median，n = 1k / 10k / 100k 区间，2026-08-04）

### 构造

| 方案 | 1k | 10k | 100k |
| --- | ---: | ---: | ---: |
| intspan add_pair（合并） | 0.0475 ms | 1.078 ms | 2.422 ms |
| coitrees new | 0.0284 ms | 0.325 ms | 6.262 ms |
| lapper new | 0.0161 ms | 0.519 ms | 6.870 ms |

### 点成员查询（n 个随机点；用户核心场景）

| 方案 | 1k | 10k | 100k |
| --- | ---: | ---: | ---: |
| intspan contains（普通/超长） | 0.0099 / 0.0077 ms | 0.0846 / 0.0500 ms | 0.197 / 0.200 ms |
| coitree point query（普通/超长） | 0.0163 / 0.0339 ms | 0.688 / 0.914 ms | 28.460 / 32.512 ms |
| lapper count（普通/超长） | 0.0146 / 0.0141 ms | 0.807 / 0.815 ms | 12.093 / 11.972 ms |

### 区间重叠枚举（n 个 2 kb 查询窗口）

| 方案 | 1k | 10k | 100k |
| --- | ---: | ---: | ---: |
| coitree query（普通/超长） | 0.0173 / 0.0354 ms | 0.688 / 0.903 ms | 28.445 / 32.407 ms |
| lapper find（普通/超长） | 0.0081 / 0.204 ms | 0.611 / 18.374 ms | 22.287 / **1583.1 ms** |

## 结论

1. **纯点成员查询用 IntSpan**：100k 时比 coitree 快 ~140x、比 lapper 快
   ~60x（合并后 span 数远小于区间数，`contains` 为二分）；代价是只有
   "是否覆盖"而没有区间身份。超长区间对它无影响。
2. **重叠枚举（需要区间身份）**：普通数据 lapper 最快（~1.2x），但
   rust-lapper README 自己警告的病态场景确实发生——一条覆盖 90% 染色体的
   区间使 lapper find 在 100k 时从 22 ms 恶化到 **1.57 s**（~71x）；
   coitrees 查询时间几乎不受影响。
3. **构造**：三者量级相当；100k 时 intspan 反而最快（add_pair 顺带合并）。

建议：pgr 的点成员 / 掩码路径维持 IntSpan；需要枚举重叠区间的生产代码
（PAF/pbit 索引）维持 coitrees（最坏情况有保证）；若某数据集确认不存在
超长区间，可考虑用 rust-lapper 换查询速度。

---

## 证据附录：Range 字符串解析基准：正则 vs 手写扫描器

> 目的：`libs::ds::Range` 的 `from_str` 用正则解析区间字符串（如
> `S288c.I(-):27070-29557`），对正则性能存疑，验证手写解析器能否等价替代
> 并提速。2026-08-04 复测。

## 方案（已切换）

* **生产路径**：`Range::from_str` → `decode`，已切为手写逐字节扫描器；
  复刻正则的全部语义：非锚定最左匹配、贪婪 name/chr/strand、`start` 缺
  `-end` 时 `end = start`、显式 `end = 0` 视为缺失（`c:911_0` → 911-911）、
  无匹配时回退为第一个空白 token；唯一差异：数字溢出 i32 时返回 0 而不是
  panic（修复了正则路径的 `parse::<i32>().unwrap()` panic）。
* **正则保留为文档**：原正则原文与语义说明在 `src/libs/ds/range.rs` 的
  模块文档里；测试模块保留正则解码器作为对拍 oracle，基准文件里也留了
  一份作为对比基线。

等价性由 `regex_and_manual_decoders_agree` 保证：固定语料（含
`foo I:1-100`、`a.b.c:1-2`、`1:-100` 等边界）+ 2 万条随机 fuzz 逐字段对拍
一致。

## 执行

```bash
cargo bench --offline --bench range_parse_benchmark
```

语料为 17 条真实格式混合（普通 / 带链向 / 带物种前缀 / 单坐标 /
斜杠下划线 contig / 回退用例）。

## 复测结果（median，2026-08-04）

| 方案 | 17 条语料整体 |
| --- | ---: |
| 正则（基线） | 5.610 µs |
| 手写 `from_str`（生产） | 1.153 µs |

手写版快约 **4.9 倍**（单条约 68 ns vs 330 ns）。三次 criterion 运行
（7.651/5.066/5.610 µs vs 0.859/1.015/1.153 µs）中位数有波动，加速比稳定
在 **5–9 倍**区间，已作为生产实现。
