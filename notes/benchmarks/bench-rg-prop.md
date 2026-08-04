# `rg prop` 命令行基准：pgr vs rgr（同为 IntSpan 交集算法）

> 目的：对比 `pgr rg prop`（`libs/runlist::range_prop`）与外部 `rgr prop`
> 的耗时与内存。与 `rg count` 不同，两侧用的是**同一套 IntSpan 交集算法**
> （pgr 已迁入 intspan crate 代码），预期接近持平；基准用于确认差距量级。
> 2026-08-04 实测。

## 环境与版本

* pgr：本仓库 release 构建（v0.4.0，含 `rg prop`）
* rgr：`~/.cbp/bin/rgr` 0.8.6（release）
* 机器：本机（hyperfine 3.x，`/usr/bin/time -v` 量内存）

## 数据（合成，种子 20260805）

* runlist `rl.json`：8 条染色体（chr1..chr8）× 100 Mb，每染色体 25,000 个
  随机区间（长度 100–2,000）合并后共 **154k spans**（约覆盖 30% 基因组，
  JSON 2.7 MB）
* `target.100k.rg`：100,000 条随机查询区间（与 count 基准同一文件）

## 复现命令

```bash
pgr rg prop rl.json target.100k.rg -o /dev/null
rgr  prop rl.json target.100k.rg -o /dev/null

hyperfine --warmup 1 --runs 3 \
  'pgr rg prop rl.json target.100k.rg -o /dev/null' \
  'rgr  prop rl.json target.100k.rg -o /dev/null'
```

## 结果（100k target，3 次取均值）

| 实现 | 时间 | RSS |
| :--- | ---: | ---: |
| pgr `rg prop` | 6.113 ± 0.040 s | 15.8 MB |
| rgr `prop` | 5.773 ± 0.011 s | 9.5 MB |

rgr 约快 **1.06×**（pgr 慢 ~6%）；内存 pgr 约 1.7×。

补充验证（8 次运行）：pgr 6.208 ± 0.021 s，rgr 5.844 ± 0.037 s——差距
稳定（~6%），非噪声。

## 正确性验证

20k target 上两者输出 `sort` 后 `diff` 为空（20,000 行逐行一致）。

## 分析

1. **算法同源，性能持平符合预期**：`prop` 的核心是每条 target 与 runlist
   `IntSpan` 求交集并数基数；pgr 迁入的就是 intspan 的同一份代码，因此
   没有 count（coitrees vs lapper）那种结构性的性能差。
2. **耗时主因是逐查询的全量交集**：每染色体 ~19k spans，100k 次
   `intersect`（O(spans)）共约 20 亿次操作，构成 6 s 的主体；解析/输出
   占比很小。
3. **~6% 差距的可能来源**：pgr 的 runlist 加载走 `IntSpan::valid` +
   `IntSpan::from` **两次解析**（`json_to_set` 先校验再构造），rgr 只解析
   一次；154k spans 的重复解析可解释大部分差距。次要来源是 pgr 的
   `usable_range`/`Range` 守卫开销。
4. **优化方向（未做）**：`prop` 可用"按染色体排序 span 数组 + 二分定位
   重叠段"把单查询从 O(spans) 降到 O(log n + k)（k 为实际重叠段数），对
   稀疏覆盖数据可再快一个量级；或对目标也建前缀和。属后续优化，本次仅
   记录现状。

## 为什么 pgr 没有像 count 那样领先？（2026-08-04 复核）

* **count 的 3.4× 来自索引结构代差**：count 是区间查询——pgr 用 coitrees
  （查询有界），rgr 用 rust-lapper；数据结构不同，差距是结构性的。
* **prop 两边跑同一份代码**：`intersect`（complement→merge→invert，
  O(spans)）在 intspan 0.8.6↔0.8.7 与 pgr vendored 版本间 diff 为空，
  无索引可优化，持平是必然。
* **定位实验**：
  * 小 runlist（1 span/染色体）：pgr 50 ms vs rgr 90 ms——pgr 反而快
    1.8×（手写 `Range` 扫描器 vs rgr 0.8.6 的正则解析，逐行解析优势）。
  * 空 target（纯加载）：pgr 16 ms vs rgr 11 ms——加载差 5 ms，不是主因。
  * 差距随 runlist 规模出现（小规模 pgr 赢、大规模 rgr 赢 6%）→ 差距在
    交集热路径的常数因子，而该路径两侧代码相同。
* **~6% 差距的最可能解释**：编译产物差异。rgr 0.8.6 二进制在 GitHub
  Actions 上构建（`strings` 可见 `/home/runner/work/intspan/...` 路径），
  优化参数未知（可能 LTO / target-cpu）；同一份代码在不同工具链下差
  ±5–6% 属正常范围。pgr 侧的小开销（`json_to_set` 双重解析、跨模块
  `usable_range`/`range_prop`、`cardinality` i64 版）只影响加载/常数，
  不足以解释 340 ms 量级的差距。
* **结论**：prop 持平（略慢 6%）是"算法同源 + 编译差异"的结果，不是 pgr
  代码回归；若要在 prop 上反超，唯一正道是换算法（§分析 4）。

## 附：消除"双重解析"（try_from）的两次尝试（2026-08-04）

`json_to_set` / `read_runlist` 对每个 runlist 字符串先 `IntSpan::valid`
再 `IntSpan::from`，解析两遍（154k spans 约多 6 ms 加载）。

* **第一次尝试（intersect 版热路径）**：加 `IntSpan::try_from` 后加载
  16.2→10.1 ms，但 prop 从 6.2 s 恶化到 7.1 s（+14%，稳定）——查询循环
  未动，纯属 LTO 全程序优化被新增 pub 函数扰动（无 LTO 构建 9.8 s）。
  当时热路径是 memmove 密集的 intersect，14% 代价大，故还原。
* **第二次尝试（二分版热路径）**：换算法后热循环只剩 ~35 ms，LTO 扰动
  的影响面大幅缩小；实测 try_from 版 35.1 ms vs 双解析版 48.7 ms——
  除 ~6 ms 加载收益外，LTO 决策这次还偏向有利（查询循环更快）。
* **结论**：保留 `try_from`。代价与收益取决于热路径形态：intersect 时代
  微优化不值得（14% 回归），二分时代净收益明显（~28%）。给热路径函数加
  `#[inline]` 稳定 LTO 决策仍可作为后续保险。

## 附：IntSpan 集合运算线性化（intersect/union/diff/xor，2026-08-04）

prop 的火焰图暴露出 `intersect` 内部 `complement+merge+invert` 链的代价：
`merge` 对每个区间调 `add_pair`，在大集合上做 O(n) VecDeque 搬移，最坏
O(n·m)。据此把三个核心集合运算重写为线性双指针合并（O(n + m)，结果
直接拼 edges，不再经 `add_pair`）：

* `intersect`：两指针走重叠区间；
* `union`：两指针归并 + 合并重叠/相邻 span；
* `diff`：逐 self span 减去 other 的重叠段；
* `xor`：`union(...).diff(&intersect(...))`，随三者变快。

`merge`/`subtract`/`add_pair` 作为增量 API 保留（`rg_to_set` 等按序插入
场景仍是 O(1) 摊还）。正确性：新旧实现差分测试（含随机 200×200 集合对
与空集/无穷集边界）逐字节一致；全量测试通过。

`pgr runlist compare`（两个 154k-span runlist，8 次取均值）：

| op | 旧实现 | 新实现 |
| :--- | ---: | ---: |
| intersect | 183.7 ms | **23.6 ms**（~7.8×） |
| union | 215.6 ms | **40.7 ms**（~5.3×） |
| diff | 222.0 ms | **34.0 ms**（~6.5×） |
| xor | 308.7 ms | **51.7 ms**（~6.0×） |

剩余时间主要被 runlist JSON 加载（~16–25 ms）占据；集合运算本身已降到
毫秒级。`runlist span` / `venn` / alignment trim 等 `intersect` 消费者
同步受益。

## 附：IntSpan 模块系统性审视（2026-08-04）

`intersect` 的发现促使对整个 `libs/ds/intspan` 做了一轮系统性审查，重点
是同类 O(n·m) / O(n²) 模式：

### 已修复

1. **批量构建 O(n²)（`add_pair` 逐个插入 unsorted 区间）**：`rg_to_set`
   原逐行 `add_pair`，1M 随机稀疏区间 cover 需 2.5 s（有序输入 0.1 s，
   差 25×）。新增 `IntSpan::from_pairs`（排序 + 单遍合并，O(n log n)），
   `rg_to_set` 改为收集 pairs 后批量构建；再抽 `rg_files_to_set` 供
   cover/trf/repeat 多文件共用（原多文件 `merge` 大集合又是 O(n·m)）：
   * cover 单文件 1M 稀疏：2.5 s → 181 ms（~14×）
   * cover 两文件 1M 稀疏：8.4 s → 299 ms（~28×）
2. **`find_islands_ints` O(n·m) → O(n+m)**：原对每个 self span 做一次
   `intersect`；改为双指针收集与 other 重叠的整段 island（alignment
   trim/variation 使用）。
3. **`distance` O(n·m) → O(n+m)**：双指针求最近相邻边界间隙；顺带消除
   原 `(lower - upper).abs()` 的 i32 溢出风险（i64 中间量 + 饱和）。

以上均带新旧实现差分测试（含空集/无穷集边界），全量测试通过。

### 审查后确认可接受

* `merge` / `subtract` / `add_pair` 保持增量 O(n) 搬移语义：重构后大集合
  路径已避开（`union`/`diff`/`from_pairs`/`rg_files_to_set` 均为线性或
  排序构建）；剩余调用点是 alignment 逐记录的小规模合并、`depth_runs`
  扫描线的有序插入、GFF 有序文件。
* `add_ranges` 对无序输入仍是 O(n²)，但调用方（`add_runlist`/`merge`）
  传的是有序 span，安全。
* `at`/`index`/`slice` 为 O(n) 扫描 + 文档化 panic，量级可接受。
* `to_vec`/`elements` 对超大 span 有内存爆炸风险（vendored 语义，未用）。

### 结论

IntSpan 的"昂贵"集中在这类**经 `add_pair` 在 VecDeque 中间做 O(n) 搬移**
的路径；凡能换成排序构建或双指针线性合并的地方都已替换。模块现在的主要
热点是解析与 IO，集合运算均为线性或 O(n log n)。

## 结论

`pgr rg prop` 与 `rgr prop` 性能基本持平（rgr 快 ~6%），内存多 ~1.7×，
输出逐行一致。若 prop 进入高频路径，优先做"二分重叠段"优化并消除
runlist 双重解析，可预期反超；当前量级（100k target 约 6 s）对交互使用
可接受。

## 更新（2026-08-04）：二分重叠段优化已实现，pgr 反超 ~120×

火焰图（perf + inferno，25556 采样）确认 6 s 的绝对主体是 `intersect`
链中的 `add_pair`——在 ~19k 段 runlist 上做 O(n) VecDeque 搬移
（`add_pair` 总耗时 66.7%，其中 `remove` 64.3% / `memmove` 21.6%）。
据此实现 `libs/runlist::SpanIndex`：利用 IntSpan span 有序且互不相交的
不变量，把重叠段定位从线性扫描换成两次 `partition_point` 二分
（重叠段必为连续区间），单查询 O(spans) → O(log n + k)。

复测（8 次取均值）：

| 实现 | 时间 | RSS |
| :--- | ---: | ---: |
| pgr `rg prop`（二分 + try_from） | **35.1 ± 1.0 ms** | — |
| pgr `rg prop`（二分，双解析） | 48.7 ± 1.1 ms | 15.6 MB |
| rgr `prop` | 5.820 ± 0.022 s | 9.5 MB |

* pgr 自优化前（6.2 s）提升 ~177×；相对 rgr 快 **~166×**。
* 正确性：20k target 输出与优化前（已与 rgr 逐行一致）sort 后 diff 为空。
* 顺带消除了每查询的 `IntSpan::new` + `add_pair` + `cardinality` 分配，
  并保留 `try_from`（runlist 单次解析，见下）。
* 教训：基准里"持平/慢 6%"的结论是**算法同源**的必然，而换算法后是
  数量级反超——性能差异主要看数据结构，不是解析/常数。
