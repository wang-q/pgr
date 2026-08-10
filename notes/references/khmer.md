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
| 哈希 | `src/oxli/kmer_hash.cc` | 2-bit 滚动哈希 `_hash`、MurmurHash3 `_hash_murmur`、cyclic `_hash_cyclic` 三种 |
| 存储 | `include/oxli/storage.hh` | BitStorage / NibbleStorage / ByteStorage / QFStorage |
| 表 | `src/oxli/hashtable.cc` + `include/oxli/hashtable.hh` | `Hashtable` 基类（2-bit）→ `Hashgraph` 派生 `Countgraph`/`SmallCountgraph`/`Nodegraph`；`MurmurHashtable` 派生 `Counttable`/`SmallCounttable`/`Nodetable`/`QFCounttable`；`CyclicHashtable` 派生 `CyclicCounttable`。提供建表、`consume_string`、`median_at_least` |
| 绑定 | `khmer/_oxli/graphs.pyx` | `Countgraph`/`SmallCountgraph`/`Nodegraph`（← `CpCountgraph`/`CpSmallCountgraph`/`CpNodegraph` → oxli::Hashgraph，**2-bit**）；`Counttable`/`SmallCounttable`/`Nodetable`/`QFCounttable`（← oxli::MurmurHashtable，**MurmurHash3**） |
| 脚本 | `scripts/normalize-by-median.py` | diginorm 流程（422 行，逻辑很薄） |

> **哈希选择的关键差异（§3 详述）**：diginorm 实际用的 `Countgraph` 继承
> `Hashgraph`→`Hashtable`（基类 `new_kmer_iterator` 返回 `TwoBitKmerHashIterator`、
> `hash_dna` 走 `_hash`），因此走 **2-bit 编码**；而独立表类
> `Counttable`/`SmallCounttable`/`Nodetable`/`QFCounttable` 全部继承
> `MurmurHashtable`，走 **MurmurHash3**。两种哈希的 canonical 化方式也不同
> （2-bit 取 `min(f,r)`，murmur 取 `f^r`）。

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
  canonical = `f ^ r`，自互补特判 `rev==kmer` 时直接返回 `f` 不再异或）；cyclic
  hash（fwd+rev 求和）。**注意**：在 3.0/oxli master 里"可选"已成部分类
  的**默认**——`Counttable`/`SmallCounttable`/`Nodetable`/`QFCounttable` 都继承
  `MurmurHashtable`（用 `_hash_murmur`），只有 `CyclicCounttable` 用 cyclic；
  2-bit `_hash` 仅保留给基类 `Hashtable`（即 `Countgraph`/`Nodegraph` 等 graph
  类）。因此 **canonical 化有两条路**：graph 类 2-bit 走 `min(f,r)`，
  murmur 表类走 `f^r`（不保证等于 2-bit 的 min 语义），移植对照时须分清。
- 表内寻址一律 `khash % tablesize`，tablesize 取**素数**
  （`get_n_primes_near_x(n, x)`：从 `x-1` 起向下数奇数、逐个判素数，
  返回**不高于 x 的 n 个素数**；`x==1` 时直接返回 `{1}`（非素数），
  x 过小可能不足 n 个）。素数表规避 `%` 与 2 的幂取模的碰撞聚集热点。
- **约束**：CLI 层强制 `k ≤ 32`（超出报错退出，`create_countgraph`/
  `create_nodegraph`）、`n_tables ≤ 20`（默认需 `--force` 才能越过）——
  k 上限来自 u64 2-bit 编码，n_tables 上限是经验值。

## 4. 计数存储（storage.hh）：Count-Min Sketch

`ByteStorage` 即标准 CMS：

- `n_tables` 张 byte 表，每表 `tablesize` 个 **u8 计数器**（0-255）；
  `add(khash)` 对每张表 `khash % tablesize` 处原子 +1；
  `get_count` 取各表最小值（标准 CMS 的 min 语义）。
- **饱和处理**：所有表都满（255）时，若 `_use_bigcount` 则改走
  `_bigcounts`（`KmerCountMap = unordered_map<HashIntoType,u16>`，
  即 `BoundedCounterType = unsigned short`，上限 `MAX_BIGCOUNT=65535`，
  首个超过的 key 记 `_max_count+1 = 256`）。
  **默认关闭**——基类 `Storage()` 构造把 `_supports_bigcount=false`
  `_use_bigcount=false`（`ByteStorage` 构造才置 `_supports_bigcount=true`），
  normalize-by-median 的 argparse 没有 `bigcount` 属性，不会开启；
  仅 `load-into-counting.py`/`abundance-dist.py` 加了
  `-b/--no-bigcount, dest='bigcount', default=True`（即 **bigcount 默认开，
  `-b` 关掉**），`create_countgraph` 里 `if hasattr(args,'bigcount')` 才调用
  `set_use_bigcount`。
  两个并发细节：`add` 写 bigcount 前用 `__sync_bool_compare_and_swap`
  自旋锁保护；`get_count` 仅在 `min_count == max_count`（已饱和）时才去
  `_bigcounts` 里查，未饱和时零开销。另注：`add` 里多线程可能让 u8 计数
  略超出 `_max_count`（代码注释明示可接受的小 slop）。
- **BitStorage（Nodegraph）= Bloom filter**：每表 `tablesize` 个 bit，
  每表分配 `tablesize/8+1` 字节。`test_and_set_bits` 用
  `__sync_fetch_and_or` 原子置位，`_occupied_bins` 只在**表 0** 置新位时
  +1（作为全局代理），`_n_unique_kmers` 首次发现新 k-mer 原子 +1；
  `get_count` 需**所有表**对应位都置 1 才返回 1（多表 bit AND，
  即标准多哈希 Bloom filter）。`Nodegraph::update_from` 用按字节 OR 合并
  （Bloom 可 union），并借 `__builtin_popcountll(me^tmp)` 统计新增 occupied
  位。**1-bit 存在性 → nodegraph 只有 0/1，无计数**。
- **NibbleStorage（SmallCountgraph）= 4-bit CMS**：每表 `tablesize/2+1`
  字节，每字节两个 nibble。寻址：`_table_index=(k%tablesize)/2`，
  `_mask`/`_shift` 由 `(k%tablesize)%2` 决定高/低半字节（240/15、4/0）。
  `_max_count=15`，计数到 15 即饱和停加（防溢出）；每表配一个
  `std::mutex`（固定 32 个 mutex 池，故断言 `n_tables ≤ 32`）。
- **QFStorage（QFCounttable）= 计数商过滤器 CQF**：单表
  `qf_init(&cf, 1ULL<<size, size+8, 0)`（size 必须 2 的幂，来自 Cython 校验；
  第 2 参为槽数 `2^size`，第 3 参 `size+8` 为 key 位数，末参 value 位数=0 未用），
  `add`/`get_count` 用 `khash % cf.range` 寻址，计数用 `qf_count_key_value`。
  khmer 3.0 试验特性（底层 `third-party/cqf/gqf.h` 不在本快照内，为外部依赖），
  每槽约 1.3 字节，槽被占满后会**停止接受 `add`**（内存不能像 CMS 那样
  预先严格固定），`get_tablesizes` 返回 `cf.xnslots`。
- **磁盘存储格式**（`SAVED_SIGNATURE="OXLI"`，`SAVED_FORMAT_VERSION=4`）：
  头部固定 `OXLI`(4B) + version(1B) + `ht_type`(1B)，类型常量
  `SAVED_COUNTING_HT=1`、`SAVED_HASHBITS=2`、`SAVED_SMALLCOUNT=7`、
  `SAVED_QFCOUNT=8`。随后因类型而异——ByteStorage/Countgraph 文件：
  `use_bigcount`(1B) + `ksize`(4B) + `n_tables`(1B) + `n_occupied`(8B)，
  然后每表写 `tablesize`(8B)+该表原始字节（Nibble/Bit 分别为
  `tablesize/2+1`、`tablesize/8+1` 字节/表），末尾写 `n_counts`(8B) +
  逐条 bigcount 条目（`kmer` u64 + `count` u16）；`BitStorage`/`NibbleStorage`
  文件则无 `use_bigcount` 段。`.gz` 用 zlib `gzread/gzwrite` 读写，按扩展名
  自动分派（`ByteStorageFile` 看文件名末段是否 `gz`）。
- 每字节桶数 `_buckets_per_byte`（`khmer/__init__.py` 字典）：
  countgraph=1、smallcountgraph=2、nodegraph=8、qfcounttable=1/1.26。
  注意 `calculate_graphsize` 给出的 `tablesize` 单位是**桶**（entries），
  实际内存 = `n_tables × tablesize / _buckets_per_byte` bytes；countgraph
  因 `_buckets_per_byte=1` 才简化为 `n_tables × tablesize`。

## 5. 参数与内存（khmer_args.py）

- 默认：`k=32`、`-N/--n_tables=4`、`-x/--max-tablesize=1e6`（太小，脚本会警告）、
  `-M/--max-memory-usage` 上限、`--small-count` 切 4-bit 表
  （SmallCountgraph/NibbleStorage，max 15）。
  `-M` 用 `memory_setting` 解析：支持裸数字/科学计数/`K|M|G|T` 后缀，
  **十进制 1000 进制、不带尾 `B`**。
- 自动配置公式（`estimate_optimal_with_K_and_M` / `_and_f`，二者共用
  fp 估计 `fp ≈ (1 - e^{-n/H})^Z`，n = unique k-mer 数，H = 单表 size，
  Z = 表数）：给定内存时 `Z = ln2 · mem/n`（`estimate_optimal_with_K_and_M`），
  给定目标 fp 时反推 H（`estimate_optimal_with_K_and_f`）。
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
  `min_req = 0.5 + float(len-k+1)/2`（对整数即 `ceil((len-k+1)/2)`），
  统计 k-mer 计数 ≥ cutoff 的个数，
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
5. **canonical 化与 pgr `.pkt` 同族**：diginorm 的 `Countgraph` 用 2-bit，
   khmer `uniqify_rc = min(f, r)` 取
   正/反向互补两条 **packed 2-bit 数值**的较小者；pgr `.pkt` 的 canonical
   key 取"正/反向互补 2-bit 编码中**字典序**较小者"。在固定字母表下
   （左侧碱基是 packed 整数的高位），序列字典序 == packed 数值序，二者
   是**同一 canonical 化逻辑**（都取 fwd/rev 编码较小的那条）。pgr 的
   `.pkt` 与 khmer 在"选择哪条链"上结论一致，可互为一致性参照；区别只在
   pgr 精确排序去重、khmer 用哈希进概率表。
   **告诫**：该同族关系**仅对 2-bit/graph 类成立**；khmer 的 murmur 表类
   （`Counttable` 等）canonical 走 `f ^ r`，与 pgr `.pkt` 的字典序取小
   **不是同一逻辑**，若拿 murmur 类表作对照需单独核对。
6. **当前 pgr 状态核实（2026-08 审计）**：`pgr fq norm` 现走**精确
   canonical KmerTable** + bbnorm 逐 read 判定（truedepth/depthAL 分位数
   + toss），计数表非近似（见 `notes/audit/audit-fq.md`、
   `notes/design/anchr-trim-replace.md` §M6）；`pgr kmer table` 走精确
   `.pkt` 排序表。因此 khmer 的 CMS/median 判定只作为"未来若新增近似路径"
   的参考，当前精确路线下不直接落地——与 §9 结论一致。

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
