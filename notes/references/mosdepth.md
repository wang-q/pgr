# mosdepth: 快速 BAM/CRAM 深度计算器

> 整理于 2026-08，源自对 `mosdepth-master`（v0.3.14，Brent Pedersen，MIT）源码的分析。
> 目的：理解 BAM 深度/覆盖度计算的主流实现（CIGAR 粒度事件差分 + 逐碱基数组），
> 与 pgr 的 `rg coverage`（.rg 区间扫描线）对比，评估可借鉴的语义与边界处理。
> 源码：`mosdepth.nim`（989 行）+ `depthstat.nim`（统计）+ `int2str.nim`（整数格式化）。

## 1. 简介

mosdepth 是面向 WGS / exome / targeted sequencing 的 BAM/CRAM 深度计算工具，
按染色体单遍流式处理，输出 per-base 深度、窗口/区域平均深度、阈值统计、
深度分布和汇总表。核心实现是**差分数组（difference array）+ 前缀和（cumsum）**：
把每条 read 的覆盖区间转成 `start: +1 / end: -1` 事件写入按位置索引的数组，
一次 cumsum 得到逐碱基深度——与 pgr `depth_runs` 的扫描线是同一思想，只是
mosdepth 用染色体长度数组直接索引（O(chrom) 内存），pgr 用事件排序（O(n) 内存）。

## 2. 输出与命令行

`mosdepth [options] <prefix> <BAM-or-CRAM>`，输出（BGZF + CSI 索引，`tbx` UCSC
preset）：

| 文件 | 内容 | 触发 |
| :--- | :--- | :--- |
| `*.per-base.bed.gz` | 逐碱基深度，游程压缩 `chrom start stop depth` | 默认（`-n` 关闭） |
| `*.regions.bed.gz` | 每个窗口/BED 区域的 mean（或 `-m` median）深度 | `--by <window\|bed>` |
| `*.quantized.bed.gz` | 相邻碱基合并进同一深度 bin（如 `1:10:20`） | `--quantize` |
| `*.thresholds.bed.gz` | 每个区域内深度 ≥ 各阈值的碱基数 | `--thresholds` + `--by` |
| `*.mosdepth.global.dist.txt` | 每染色体/全基因组"深度 ≥ X 的碱基占比"累积分布 | 总是 |
| `*.mosdepth.summary.txt` | 每染色体/总计的 length/bases/mean/min/max | 总是 |
| `*.per-base.d4` | d4 格式逐碱基深度（`--d4`，编译期开关） | 可选 |

关键过滤：`-F/--flag` 默认 **1796**（0x704 = unmapped | secondary | QC-fail |
duplicate），`-Q/--mapq`、`-l/-u/--min-frag-len/--max-frag-len`（按 `abs(isize)`）、
`-R/--read-groups`、`-c/--chrom`（支持 `chr:start-end`）。CRAM 需要 `-f/--fasta`
或 `REF_PATH`（`check_cram_has_ref` 缺引用直接 quit(1)）。

## 3. 深度计算核心

### 3.1 事件生成：CIGAR 粒度（`gen_start_ends` / `inc_coverage`）

对每条 read 的 CIGAR 逐 op 扫描，只处理 `consumes.reference` 的 op：

* **只有消耗 reference 的 op 才推进参考位置**（`if not con.reference: continue`；
  insertion/soft-clip 等只消耗 query 的 op 直接跳过，其长度不进入参考坐标；
  match `M`/`=`/`X` 同时消耗 query 与 reference）；
* **deletion / intron（N）消耗 reference 但不消耗 query，推进参考位置但该段不计深度**——
  即 `-split` 语义：每个对齐块贡献独立的 `(pos, +1)` 与 `(last_stop, -1)`，两个块之间的
  deletion/N 区域深度为 0（如块 `[10,15)`、`[20,25)` 中间 5 bp deletion：事件
  `+1@10, -1@15, +1@20, -1@25` 归位后 `[15,20)` 深度 0）；
* 相邻参考块之间只有 deletion/N 时**不重复开段**（`pos == last_stop` 检查），
  块的起止事件精确对齐 CIGAR 的 M/= /X 段。

事件直接写入差分数组：`arr[p.pos] += p.value`，单条 read 的多个块产生多个
`+1/-1` 事件。`to_coverage` 对全数组 `cumsum()` 得到逐碱基深度。

### 3.2 三种计数模式

| 模式 | 行为 | 用途 |
| :--- | :--- | :--- |
| 默认 | CIGAR 粒度 + **mate 重叠校正**（见 3.3） | 精确深度 |
| `-x --fast-mode` | 只看 `rec.start/rec.stop`（read 外边界），不看内部 CIGAR、不校正 mate | 快速、绝大多数场景推荐 |
| `-a --fragment-mode` | 只统计 proper pair 的 read1，覆盖整个 fragment（`start ~ start+abs(isize)`） | fragment 深度 |

### 3.3 mate 重叠校正

proper pair 且两条 read 在同染色体、有重叠时（`rec.stop > matepos`），若不校正，
重叠区会被双端各计一次。mosdepth 用 `seen` 表（key = qname）暂存位置靠前的
第一条 read；第二条到达时，把两条 read 的 CIGAR 事件**合并排序**，累积
`pair_depth`，凡深度达到 2 的区间（即重叠段）做 `dec(last_pos) / inc(pos)`，
把重叠区从差分数组中扣掉（注释中示例：`pair.depth` 累积到 2 时开始递减）。
单 CIGAR 的常见情况有快速路径（`n_cigar == 1` 时直接
`dec(arr[rec.start]) / inc(arr[mate.stop])`）。

### 3.4 输出生成

* **per-base**：`gen_depths` 遍历逐碱基数组，深度变化处切分，输出
  `(start, stop, depth)` 游程；整数转字符串用自带的 `fastIntToStr`（Milo Yip
  itoa 查表法，输出是主要性能瓶颈）。
* **regions**：窗口（`window_gen`，半开 `[start, start+window)`）或 BED 区域
  （`region_gen`，按染色体预读、边消费边删）；mean 直接求和除以长度，
  `-m` median 用 `CountStat` 直方图（65536 桶，超出 65535 的深度计入末桶，
  高深度时中位数是近似值）。
* **thresholds**：`write_thresholds` 对每个 region 逐碱基统计"深度 ≥ 各阈值"的碱基数；
  阈值经 `threshold_args` 排好序，内层 `if v < t: break` 提前跳出（`v < t` 后更大的阈值也
  必然不满足，避免扫完全部阈值）；无数据的染色体（`tid==-2`）直接输出全 0 行。
* **quantized**：`gen_quantized` 把相邻、同 bin 的碱基合并，bin 由
  `--quantize 1:10:20` 定义（`:inf` 表示开区间；环境变量 `MOSDEPTH_Q*` 可
  覆盖输出标签）。
* **distribution**：直方图数组从 512 动态增长（`inc` 在 `v >= len` 时
  `set_len(v+10)`），深度 > 400000（`MAX_COVERAGE`）截断到 399990（只影响
  分布，不影响 per-base/quantized）；输出时反序 cumsum，得到
  "深度 ≥ X 的碱基占比"；先跳过深度索引 > 300 的 0 计数桶（`irev > 300 and
  v == 0`），再跳过累计占比 < 8e-5 的行。
* **summary**：length / bases（= 深度总和，即总比对碱基数）/ mean / min / max；
  空染色体 min 记 0；数值精度由 `MOSDEPTH_PRECISION` 控制（默认 2）。

## 4. 性能设计

* **单遍 + 数组复用**：`arr.init(tlen+1)` 只分配一次，后续染色体
  `set_len` + `zeroMem` 复用；深度数组为 `seq[int32]`（约 4 B/碱基，1 Gb
  染色体约 4 GB 内存）。
* **输出 I/O**：BGZF 压缩级别 1；BGZI + CSI 分级（`get_min_levels` 按最长
  染色体定级，1<<14 起 ×8）；`fastIntToStr` 避免 `$int` 的通用格式化开销。
* **线程**：`-t/--threads` 只用于 htslib BAM 解压，深度计算单线程。
* **CRAM 优化**：按需裁剪 required fields（非 fast-mode 才要
  QNAME/RNEXT/PNEXT），关闭 MD 解码。
* 构建要求 `--mm:refc`（README 明确"必须"，影响性能与正确性）；
  `GC_disableMarkAndSweep`。

## 5. 与 pgr 的对比

### 5.1 相同点

* **差分/扫描线思想**：`start +1 / end -1` 事件累积 → 深度。pgr
  `libs/runlist::depth_runs` 对 `(start, end)` 半开区间做同样的事（事件排序后
  单遍累积），二者数学上等价。
* **游程输出**：mosdepth per-base 的 `(start, stop, depth)` 与 pgr runlist
  JSON 的 `"start-end"` span 都是对深度恒定段做压缩；mosdepth 输出 0 深度段
  （如 `0 80 1` 后跟 `80 16569 0`），pgr runlist 只输出覆盖段（语义更紧凑）。
* **quantized ≈ `pgr rg coverage -d`**：都按深度分组输出，mosdepth 按
  bin 合并、pgr 按精确深度。
* **thresholds/regions ≈ `pgr runlist stat/statop`**：区域覆盖率统计。

### 5.2 差异

| 维度 | mosdepth | pgr rg coverage |
| :--- | :--- | :--- |
| 输入 | BAM/CRAM（需索引），CIGAR 粒度 | `.rg` 区间（来自 PSL block 等），block 粒度 |
| 内存 | O(染色体长度) 差分数组 | O(区间数) 事件向量（稀疏友好） |
| mate 校正 | 有（默认模式） | 无（上游 PSL 已给出块） |
| 过滤 | flag/mapq/插入长度/read group | 无（上游完成） |
| 空染色体 | 输出 0 深度区间 | 不出现该染色体 |
| 数值上限 | int32 深度、分布截断 400000 | IntSpan 坐标上限 POS_INF-1 |

pgr 的 rept/s-align 管道本质上等价于 mosdepth 的 **fast-mode**（PSL block 的
target 区间 → `.rg`，只计块、不校正 mate），因此 pgr 侧无需实现 CIGAR 解析与
mate 校正；若未来要支持"直接读 BAM 算深度"（如 `pgr rg coverage` 直连
samtools 输出），本笔记的 CIGAR 事件语义（deletion/N 计 0、ins 不推进）是
对齐基准。

## 6. 源码中值得注意的健壮性问题（对 pgr 的参考意义）

* **BED 行 `start > end` 直接 abort**：`bed_line_to_region` 用 `doAssert
  s <= e`（Nim 的 doAssert 在 release 也生效），畸形 BED 使进程崩溃，而不是
  报错/跳过。对照 pgr 的 Zero-Panic 原则（`.rg` 反转区间跳过、runlist 非法值
  报错），这是两种工程取向的差异。
* **`get_tid` 顺序假设**：快速路径 `tgts[last_tid + 1]` 无边界检查，依赖
  BAM 目标按 tid 顺序、查询严格递增；`-c` 限制单染色体时安全，乱序 region
  查询理论上有越界风险。
* **fragment-mode 越界**：`arr[fragment_start + abs(rec.isize)] -= 1` 未校验
  fragment 末端是否超出染色体长度（Nim 默认边界检查会 IndexDefect，release
  可关）。
* **中位数近似**：`CountStat` 65536 桶截断，>65535 深度的区域 median 偏小。
* **`gen_depths` 的 `offset/istop` 参数**：全代码只在 per-base 输出处以默认值
  `gen_depths(arr)`（`offset=0, istop=0`）调用，区域用途的偏移/截断参数从未
  生效，属死参数；尾部三路 `yield` 是历史遗留（正常全染色体输出由
  `last_i < stop` 分支闭合到 `len(arr)-1`，另两路实际不可达）。
* **`gen_quantized` 末两位 off-by-one**：循环 `for pos in 0..<(arr.high-1)` 只
  比较到 `len-3` 位置处的 bin，最后两个位置不参与比较，被并入上一段的 bin
  （尾部 `yield (last_pos, len(arr)-1, ...)`）；当染色体末尾两碱基落入不同 bin
  时 quantized 输出会错标（功能测试的 MT 用例末尾深度稳定，未暴露该问题）。
* **数值**：`arr` 为 `int32`，深度 > 2^31 会溢出（实际不可能）；分布截断
  `MAX_COVERAGE` 只影响直方图。

## 7. 测试

* `functional-tests.sh`：ssshtest 框架 + 真实小型 BAM（`tests/ovl.bam`、
  `tests/overlapping-pairs.bam`、`tests/nanopore.bam`），黄金断言覆盖 mate
  重叠校正（默认 vs fast-mode 的 MT 深度）、CRAM 缺引用报错、缺失染色体、
  乱序 BED、区域越界、插入长度过滤等；构建带
  `--boundChecks:on -d:useSysAssert -d:useGcAssert`（`tests/funcs.nim` 同参数）。
* `tests/funcs.nim`（由 `tests/all.nim` 导入执行）：`depthstat` min、quantize
  参数解析、`linear_search`/`make_lookup`、threshold 解析、region 解析、BED
  index preset 等单测；CountStat 中位数单测在 `depthstat.nim` 自己的
  `isMainModule` 块里（需直接 `nim c -r depthstat.nim` 才会执行）。
