# khmer（master）：固定内存 k-mer 计数 + 数字归一化（源码分析）

> 2026-08 整理，纯源码分析（`khmer-master/`，2019-03-13 快照，3.0/oxli 重构
> 后的 master 未发布版）。khmer 是 DIB Lab 的 k-mer 工具集，代表作品是
> **数字归一化（diginorm）**：用 Count-Min Sketch 计数表做流式 read 归一化。
> 对应 pgr 需求：为 `pgr fq norm` 的"大数据量近似哈希路径"找一个现成的
> 参考实现（此前只调研过 bbnorm `bits=16`/prefilter，没有源码级对照）。

## 1. 概况

- **定位**：k-mer 计数（countgraph）、数字归一化（normalize-by-median）、
  组装图（hashgraph）等；脚本面向宏基因组/转录组 read 降冗余。
- **语言**：C++ 核心（`src/oxli/` + `include/oxli/`）+ Cython 绑定
  （`khmer/_oxli/`）+ 薄 Python 脚本（`scripts/`）。
- **数据流**：一遍流式扫描，内存**固定**（由 `-N`/`-x`/`-M` 决定），
  号称 "never runs out of memory"（与 bbnorm 官方文档同款卖点）。
- **与 bbnorm 的本质区别**：khmer 是**在线 diginorm**——边读边计数边判定，
  结果依赖输入顺序；bbnorm 是两遍式（先全量建表再过滤，结果确定）。

## 2. 架构

| 层次 | 位置 | 内容 |
|---|---|---|
| 哈希 | `src/oxli/kmer_hash.cc` | 2-bit 滚动哈希、Murmur、cyclic 三种 |
| 存储 | `include/oxli/storage.hh` | BitStorage / NibbleStorage / ByteStorage / QFStorage |
| 表 | `src/oxli/hashtable.cc` | Counttable 建表、consume、median_at_least |
| 绑定 | `khmer/_oxli/graphs.pyx` | Countgraph / SmallCountgraph / QFCounttable |
| 脚本 | `scripts/normalize-by-median.py` | diginorm 流程（422 行，逻辑很薄） |

## 3. k-mer 哈希（kmer_hash.cc / kmer_hash.hh）

- **默认 2-bit 编码**：A=0、T=1、C=2、G=3；非 ACGT 一律按 3（G）处理。
  前向 `f` 左移 2 位累积，反向互补 `r` 反向累积；
  canonical = `min(f, r)`（`uniqify_rc`）。u64 支持 k ≤ 32。
  **source quirk**：`twobit_repr/twobit_comp` 是宏，快速分支注释
  `"NOTE: Assumes data is already sanitized as it should be by parsers"`——
  它**假定输入已清洗**（否则 N/小写会走默认分支映射为 3/G）。这正是
  `clean_input_reads` 必须先把 `N→A` 的原因；若在 pgr 直接复用此编码，
  必须自行保证清洗或处理非 ACGT。
- **滚动迭代** `KmerIterator::next()`：
  `f = (f << 2 | twobit(ch)) & bitmask`，
  `r = (r >> 2) | (twobit_comp(ch) << (2k-2))`。与 pgr 现有滚动方式同构。
- **可选哈希**：MurmurHash3_x64_128（seed 0，正反链各一次，
  canonical = `f ^ r`，自互补特判）；cyclic hash（fwd+rev 求和）。
  表内寻址一律 `khash % tablesize`，tablesize 取**素数**
  （`get_n_primes_near_x`），避免取模热点。

## 4. 计数存储（storage.hh）：Count-Min Sketch

`ByteStorage` 即标准 CMS：

- `n_tables` 张 byte 表，每表 `tablesize` 个 **u8 计数器**（0-255）；
  `add(khash)` 对每张表 `khash % tablesize` 处原子 +1；
  `get_count` 取各表最小值（标准 CMS 的 min 语义）。
- **饱和处理**：所有表都满（255）时，若 `_use_bigcount` 则改走
  `_bigcounts`（`unordered_map<u64, u16>`，上限 65535）。
  **默认关闭**——`Storage()` 构造把 `_use_bigcount` 置 false，
  normalize-by-median 的 argparse 没有 `bigcount` 属性，不会开启；
  只有 `load-into-counting.py`/`abundance-dist.py` 显式 `-b` 默认 True。
- 其他存储：`BitStorage`（nodegraph，1-bit 存在性）、
  `NibbleStorage`（SmallCountgraph，4-bit，max 15）、`QFStorage`（计数商过滤器）。
- 每字节桶数 `_buckets_per_byte`：countgraph=1、smallcountgraph=2、nodegraph=8。
  **内存 = n_tables × tablesize bytes**（countgraph 情形）。

## 5. 参数与内存（khmer_args.py）

- 默认：`k=32`、`-N/--n_tables=4`、`-x/--max-tablesize=1e6`（太小，脚本会警告）、
  `-M/--max-memory-usage` 上限、`--small-count` 切 4-bit 表。
- 自动配置公式（`estimate_optimal_with_K_and_f`）：
  fp ≈ `(1 - e^{-n/H})^Z`（n = unique k-mer 数，H = 单表 size，Z = 表数）；
  给定内存时 `Z = ln2 · mem/n`，给定 fp 时反推 H。
- `-C/--cutoff` 默认 20，即 median k-mer 覆盖度阈值；范围 [0, 256)。

## 6. normalize-by-median 语义（脚本 + hashtable.cc）

```
batch = ReadBundle(read0, read1)          # 双端一对（或单端一条）
if not all(read.median_at_least(cutoff) for read in batch):
    for read in batch:                    # 任一 read 低于 cutoff → 整批保留
        countgraph.consume(read)          #   且保留的 read 才计数进表
        yield read
```

- **成对策略**：`coverages_at_least = all(...)`——只要批内**任一条** median
  低于 cutoff，两条都保留（与 bbnorm 的 keepboth 同思路）。
- **median_at_least 短路优化**（hashtable.cc）：不做全排序。设
  `min_req = ceil((len-k+1)/2)`，统计 k-mer 计数 ≥ cutoff 的个数，
  一旦 ≥ min_req 即判 `median ≥ cutoff`（前 min_req 个先扫，失败才继续）。
- **清洗**：`clean_input_reads` 把序列大写、`N → A`；短于 k 的 read 在
  `broken_paired_reader(min_length=k)` 处被过滤。
- **顺序相关**：被丢弃的 read 不计数，表随保留 read 在线演化——同一个文件
  换一种顺序处理，keep/drop 集合可能不同。这是与 pgr norm（bbnorm 语义）
  最根本的差异，移植时不能照搬整条流程，只能借用数据结构。

## 7. 对 pgr 的启示

1. **CMS 骨架可作近似路径的最小参考**：`n_tables × tablesize` 的 byte 表 +
   min 查询，实现量很小；pgr 移植时可用 splitmix64/murmur 等现成哈希 +
   mask/取模，不必抄素数表。
2. **bigcount 的教训（重要）**：高覆盖度数据（我们讨论过的 1000×）下 u8
   计数会饱和，而 bigcount 是 unordered_map，内存随**超饱和 key 数**增长，
   一样会失控。近似路线要么接受饱和（低 cutoff 判定足够），要么不提供
   精确大计数——bbnorm 的 `bits=16` 也是同理（固定 2-byte 计数、无溢出表）。
3. **median_at_least 可移植**：若 pgr norm 以后提供 khmer 式 median 语义的
   近似模式，这个短路算法直接照搬；当前 bbnorm 移植用的是
   truedepth/depthAL 分位数判定，不冲突。
4. **1TB 决策不变**：CMS 路线内存固定但精度随装填率下降（且 k=32 时
   khmer 的 2-bit 编码本身不占大内存，占内存的是计数表）；精确路线仍走
   pgr 自己的 `.pkt` + sort-merge。khmer 只补全"近似"一侧的细节。

## 8. 与 pgr 现有设施对照

| 项 | pgr 现有 | khmer |
|---|---|---|
| 精确计数 | `KmerTable`（u128+u32）、`count.rs` .pkt 排序表 | 无（近似） |
| 近似计数 | 无（bbnorm 移植为精确表语义，`bits=16` 待评估） | ByteStorage CMS（u8/u16 饱和） |
| 判定 | truedepth/depthAL 分位数 + toss（bbnorm） | median ≥ cutoff（在线） |
| 大表查询 | sort-merge 归并 / 二分 | `%` 素数表 |

## 9. 结论

khmer 的价值集中在两点：**固定内存 CMS 的最小实现范式**和
**median_at_least 短路判定**。它的整体流程（在线 diginorm）与 pgr
要保留的 bbnorm 语义（两遍式）不兼容，不应整体移植。若最终选择近似路线，
建议按 bbnorm 的 `bits=16`/prefilter 语义实现（保持与 BBTools 输出一致），
khmer 的存储结构和判定优化作旁证参考。
