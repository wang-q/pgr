# bcalm (BCALM 2)：紧凑 de Bruijn 图构建（源码分析）

> 2026-08 整理，纯源码分析（`bcalm/`，版本 `v2.2.3`）。BCALM 2 是
> Rayan Chikhi 等人的紧凑 de Bruijn 图（compacted de Bruijn graph, cdBG）构建
> 工具，ISMB 2016 / Bioinformatics 32(12): i201–i208。**后续更新**：用户补全了
> `gatb-core` 子模块（`git submodule` 已检出到 `gatb-core/gatb-core/`），核心算法
> 现可读，见 §2.2 与 §3。仓库内容构成：① `src/` 薄封装（bcalm_1.cpp 委托给 GATB
> 的 `GraphUnitigsTemplate`）；② `gatb-core/.../src/gatb/bcalm2/`——**真实算法**
> （`bcalm_algo.cpp` 分桶压缩 + `bglue_algo.cpp` 全局拼接 + `ograph.cpp` 桶内
> 压缩 + `unionFind.hpp` 无锁并查集）；③ `bidirected-graphs-in-bcalm2.md` 双向图
> 形式化定义；④ 若干脚本（convertToGFA/pufferize/split_unitigs/abundance_stats/
> unitigEvaluator）；⑤ `thirdparty/` 并发原语；⑥ `debruijn/impl/LinkTigs.cpp`
> 计算 unitig 间 `L:` 边。

## 1. 概况

- **定位**：从测序 reads（FASTA/FASTQ，可 gzip，多文件）构建**紧凑 de Bruijn 图**
  ——输出图中所有 unitig（非分支路径的序列），并报告 unitig 之间的连接边。
- **核心思想**（论文，与源码对应）：先用 dsk 做**外部排序的 k-mer 计数**
  （磁盘友好，产出 `.h5`），再按 **minimizer** 把 k-mer 分桶（超级桶 = DSK
  partition），桶内用 `graph3` 做**局部压缩**产出部分 unitig（写 `*glue*` 文件），
  最后用 **union-find（UF）把各桶产出的单元片段 "glue" 拼成完整 unitig**，并
  用 `LinkTigs` 重算 unitig 间 `L:` 边。完整三段流水线见 §3。
- **双链语义**：所有 k-mer 转 canonical 表示（k-mer 与其反向互补视为同一对象，
  RC 后只出现一次）；unitig 方向不保证跨 run 一致。
- **输入输出**：`-in`（单文件或多文件列表 `ls -1 *.fastq > list_reads`）；
  输出默认 `<prefix>.unitigs.fa`；GFA 需用 `scripts/convertToGFA.py` 后处理。

## 2. 仓库结构

```
bcalm/
├── src/            # 薄封装：main.cpp + bcalm_1.cpp/hpp，全部委托给 GATB
├── gatb-core/
│   └── gatb-core/src/gatb/
│       ├── bcalm2/                      # ★ BCALM2 真实算法
│       │   ├── bcalm_algo.{hpp,cpp}     #   分桶 + 桶内压缩（bcalm2）
│       │   ├── bglue_algo.{hpp,cpp}     #   UF 全局拼接（bglue）
│       │   ├── ograph.{h,cpp}           #   graph3 桶内压缩
│       │   ├── unionFind.hpp            #   无锁并行并查集
│       │   └── ThreadPool.h / logging  #   线程池 / 日志
│       └── debruijn/impl/
│           ├── GraphUnitigs.{hpp,cpp}   #   GATB 建图调度
│           ├── UnitigsConstructionAlgorithm.{hpp,cpp}
│           └── LinkTigs.cpp             #   计算 unitig 间 L: 边
├── bidirected-graphs-in-bcalm2/
│   └── bidirected-graphs-in-bcalm2.md   # ★ 双向图形式化定义（核心文档）
├── scripts/        # convertToGFA / pufferize / split_unitigs / abundance_stats / unitigEvaluator
├── thirdparty/     # ThreadPool / concurrentqueue / lockbasedqueue / lockstdqueue / lockstdvector
├── example/        # run-tiny、circular_unitigs_unittests、pufferize、uf
└── test/           # simple_test.sh + minitip.fa
```

### 2.1 `src/` 薄封装（`bcalm_1.cpp`）

`bcalm_1` 继承 GATB 的 `Tool`。构造里直接取 GATB `GraphUnitigsTemplate<32>`
的参数解析器（`getOptionsParser`），并隐藏了与 cdBG 无关的计数选项
（`STR_KMER_ABUNDANCE_MIN_THRESHOLD`、`STR_HISTOGRAM_MAX`、
`STR_SOLIDITY_KIND`、`STR_URI_SOLID_KMERS`），把 `STR_REPARTITION_TYPE` 与
`STR_MINIMIZER_TYPE` 默认值设为 `"1"`（frequency-based minimizer）。`execute()`
按 k-mer 大小调用 `Integer::apply` 分派到对应位宽的 `GraphUnitigsTemplate<span>`
模板实例，然后 `GraphType::create()` 建图，最后删掉 `.h5` 中间文件。注释明说：
`where did all the code go? now it's mostly in ../gatb-core/.../gatb/bcalm/`。

> `main.cpp`：`bcalm -v`/`--version` 打印版本；其余走 `bcalm_1().run()`，
> GATB `Exception` 捕获后打印 `EXCEPTION: <msg>` 并返回非零。

### 2.2 调度链（`GraphUnitigs.cpp` → `UnitigsConstructionAlgorithm.cpp`）

`GraphUnitigsTemplate<span>::create()` 读参数并决定做哪几步：`do_bcalm`（bcalm2
分桶压缩）、`do_bglue`（UF 拼接）、`do_unitigs`/`do_links`（LinkTigs 算边）。
`UnitigsConstructionAlgorithm::execute()` 依次调
`bcalm2<>()` → `bglue<>()` → `link_tigs<>()`。可选 `-skip-bcalm`/`-skip-bglue`/
`-redo-links`/`-skip-links`（`pufferize.py` 末尾就用到这些来只重算链接）。

## 3. 三段流水线（真实算法，`gatb-core/.../bcalm2/`）

BCALM 2 把 cdBG 构建拆成三段：**分桶压缩（bcalm2）→ UF 全局拼接（bglue）→
算边（LinkTigs）**。核心数据结构：`graph3`（桶内压缩）、`unionFind`（无锁并查集）、
`markedSeq`（待拼序列的端点信息）、`BooPHF`（MPHF）。

### 3.1 `bcalm2`：分桶 + 桶内压缩（`bcalm_algo.cpp`）

主循环 `for each superbucket (= DSK partition) p`，词表："partition/super-bucket"
= 一个 DSK partition；"bucket" = 一个 minimizer 对应的桶。

1. **扩展超级桶**（`InsertIntoQueues`）：对 partition 内每个 solid k-mer，若丰度
   ≥ 阈值则保留；用 `modelK1`（k-1 长 minimizer 模型）求其**左/右 minimizer**；
   凡 minimizer 属于本 partition `p` 的，把 `(minimizer, kmer, abundance, leftmin,
   rightmin)` 元组推入**本线程的 flat bucket queue**。
2. **traveller k-mers**：当一个 k-mer 的左右 minimizer 落在**不同 partition**（跨
   桶），把它以 ASCII 落到 `prefix.doubledKmers.<p>` 文件（注释里存丰度），等轮到
   `repart(max_minimizer)` 那个 partition 时再读回——这是"跨超级桶的边"，只能串
   到后续 partition 处理。代码注释解释了为何必须按 minimizer 顺序迭代 partition。
3. **按 minimizer 排序**：每线程用 `std::sort` 按元组首元素（minimizer）排序，并
   记录 `start_minimizers[thread][minimizer]` 与每 minimizer 的 k-mer 计数。
4. **逐 bucket 压缩**（`ThreadPool` 并行）：对每个实际出现的 minimizer，构造
   `graph3<SPAN>(k-1, minimizer, minSize, nb)`，把该 minimizer 的所有 k-mer
   `addtuple` 进去，调 `debruijn()` 压缩（§3.3），再把每个 unitig 重新算左右
   (k-1)-mer minimizer，标记 `lmark/rmark`（首/尾是否与当前 minimizer 不同），连同
   每个 k-mer 丰度向量写入**本线程的 glue 文件** `prefix.glue.<thread>`。
5. 每 partition 结束 flush glue 文件；最后写 `prefix.glue`（非空 glue 文件清单）。

> 线程安全细节：用 `ThreadGroup::findThreadInfo` 找线程 id 决定写哪个 queue/glue
> 文件，避免加锁；`flat_bucket_queues` 每线程一份。`BINSEQ` 宏可切换 string↔二进制
> 序列存储（默认关，作者说 "graph4 is not ready"）。
>
> 设计点：**丰度阈值过滤发生在 InsertIntoQueues**（`operator()` 里 `abundance <
> threshold` 就 return），而不是 dsk 计数阶段——所以同一份 `.h5` 计数文件可以
> 用更高的 `-abundance` 复用重跑，无需重新数 k-mer（`bcalm_algo.cpp` 注释明说）。
> 对 pgr 意义：若 `KmerTable` 做精确计数，把阈值过滤推迟到消费侧，可一次计数
> 多次调阈值。

### 3.2 `bglue`：UF 全局拼接（`bglue_algo.cpp`）

把各 glue 文件里带 `lmark/rmark` 的部分 unitig 按端点 k-mer 拼接成完整 unitig：

1. **构建端点哈希集合**（`prepare_uf`，3 趟）：对每个 glue 序列，凡 `lmark`/`rmark`
   为真的端点 (k-)mer，用 `ModelCanon` 编码后经 `Hasher_T` 哈希，按 `hash % nb_passes
   == pass` 分到 3 个趟文件 `prefix.glue.hashes.<pass>`（每趟内排序 + 去重 + k 路归并）。
2. **构造 BooPHF MPHF**（`boomphf::mphf`）：把全部端点哈希建成最小完美哈希，映射到
   `[0, nb_uf_keys)`。
3. **建无锁并查集** `unionFind ufkmers(nb_uf_keys)`：对每段 glue 序列，凡左右都标记
   （`lmark && rmark`），把 `uf[mphf(ks)]` 与 `uf[mphf(ke)]` `union_` 起来——**同属
   一个 unitig 的端点 k-mer 并到同一等价类**。UF 结果镜像成 `ufkmers_vector`
   （32 位类号）。
4. **按 UF 类分片**：把 glue 序列按 `ufclass % nbGluePartitions` 写进
   `prefix.gluePartition.<i>`（`nbGluePartitions` = min(2000, max_open_files/2)，防文件
   句柄超限）；无标记（`!found_class`）的序列直接作为最终 unitig 输出。
5. **逐片拼链**（线程池）：`determine_order_sequences` 用 `kmerIndex`（端点 k-mer →
   序列索引集合）从**链的端点**（无 `lmark` 或 `rmark` 的那端）出发，沿 `rmark`
   逐个找后继拼接成 `chain`（必要时 `revcomp` + 高位置 1 标记反链）；剩余未拼上的
   是**环状 unitig**（`while nb_chained < size` 兜底，`expect_circular=true` 检测环
   并在环处断开）。`glue_sequences` 按链顺序"砍掉后续序列首 k-mer"拼出最终序列。
6. **写头**：`make_header` 输出 `LN:i:<len>` + `KC:i:<sum>`/`km:f:<mean>`（或
   `-all-abundance-counts` 的 `ab:Z:<...>`），写 `BufferedFasta`（自实现带互斥锁的
   缓冲 FASTA 输出，`needs_consecutive_ids` 时补连续 ID）。

> 端点标记语义：glue 序列的 comment 前两字符是 `lmark rmark`（'1'/'0'）；凡
> `lmark`/`rmark` 为 1 说明该端与桶 minimizer 不一致、**必然还有另一段要接**——
> 这正是 bglue 判断"是否需拼接"的依据；两端都没标记的序列是完整 unitig，直接输出。

### 3.3 `graph3`：桶内压缩（`ograph.cpp`）

`graph3` 在**一个 minimizer 桶内**构建并压缩 de Bruijn 图。输入是若干 k-mer
（注释说"初始是 k-mer，之后成为 unitig"），输出该桶的 unitig：

- **索引**：`addtuple(seq, leftmin, rightmin, abund)` 把每条序列放入 `unitigs[]`
  数组，并把其**首/尾 (k-1)-mer 的 canonical 形式**分别登记进 `left`/`right`
  两个 `kmerIndice` 列表（带 `SEQ_LEFT`/`SEQ_RIGHT` 位置标记与序列索引）。
  左右 minimizer 与桶 minimizer 匹配的那侧才登记（`indexed_left/right`）。
- **压缩**（`debruijn()`）：把 `left`/`right` 按 (k-1)-mer 排序（末尾推 -1 哨兵
  免越界检查），双指针扫描；当 left 与 right 的 (k-1)-mer 相等时，若无其他
  (k-1)-mer 竞争（`go` 保持 true），调 `compaction(iL,iR,kmmer)` 把两条序列
  拼接。`compaction` 依据 4 种 overlap 方向（正/反、首/尾匹配）选择拼接方式，
  并把被吸收的一端改写为数字索引（`unitigs[i]=to_string(iR)`，即"压缩痕迹"，
  `isNumber` 判断），用 `compact_abundances` 合并丰度向量。
- **输出**（`output(i)`）：`isNumber(unitigs[i][0])` 为真表示已被吸收，跳过；
  否则是最终 unitig。`pre_tip_cleaning`（尖端清理，SPAdes tip 长度约定
  `<3*(k+1)`）被作者**默认关闭**（注释说 `indexed_*`/`connected_*` 未正确处理，
  且收益不大）。

> 细节：`reverseinplace` 用位运算（`^=4`、`^=17`）原地翻转并互补；k-mer 编码用
> GATB `LargeInt`（`beg2int128`/`end2int128`/`rcb`），2-bit/碱基。`debruijn()` 里
> 注释明确："先标记要压缩的对，再统一压缩，以保证所有 unitig 的 connection 信息
> 正确"——但当前 `to_compact` 方案被注释掉，改为扫描时即时压缩。

### 3.4 `unionFind` 无锁并行并查集（`unionFind.hpp`）

基于 Wenzel Jakob 的无锁并查集（Anderson & Woll 论文），**路径压缩 + 按秩合并**，
用 CAS（`compare_exchange_weak`）并发安全。每个元素存 64 位：低 32 位 parent、高
32 位 rank。`union_` 循环里 CAS 失败就重试。`normalize()`（把类号规范为最小 id）
需 3× 内存，默认关；`printStats` 输出集合数/均值/最大秩与内存。

> 注意：UF 的键是**端点 k-mer 的 64 位哈希**（`uf_hashes_t`），不是 k-mer 本身——
> 注释明说"有碰撞，但应该可以接受"；元素数受 `UINT32_MAX` 限制（超限报错退出）。

### 3.5 `LinkTigs`：计算 unitig 间 `L:` 边（`debruijn/impl/LinkTigs.cpp`）

bglue 输出的 unitig **不带连接边**；`link_tigs` 重新扫描所有 unitig，用**端点
(k-1)-mer 哈希表**找重叠，产出带 `L:` 的最终 FASTA：

- 内存受限，分 **8 趟**（`nb_passes`），每趟只处理属于该趟的 unitig 端点
  （`is_in_pass`）；每趟建 `utigs_links_map`（(k-1)-mer → 端点 `ExtremityInfo` 列表），
  再为每个 unitig 端点查 `in_links`/`out_links`，写 `<prefix>.links.<pass>`。
- 方向判定：比较 unitig 端点 (k-1)-mer 的 canonical 值与其实际序列是否同向
  （`beginInSameOrientation`）；`ExtremityInfo` 打包 (unitig id, rc 标记, 端点类型)。
  **回文特例**：当 (k-1)-mer 为回文（长度偶数时可能）用 `nevermindInOrientation`
  放宽方向。
- 输出两种格式：`edge_km_representation` 为真时输出 `J:0:<utig>:<rc>`（k-mer 级），
  否则输出 `L:-:<utig>:<sign>`；`write_final_output` 对 8 个 links 文件做 k 路归并，
  把每条 unitig 的链接追加到头（可 `renumber_unitigs` 重排 ID）。
- 最终 FASTA 头：`<id> LN:i:<len> KC:i:<sum> km:f:<mean> L:<sign>:<to>:<sign> ...`
  （`L:` 语义见 §5.1）。

## 4. 双向图模型（`bidirected-graphs-in-bcalm2.md`，核心）

> 2025 作者更新：文中这套双向图形式化**不是作者推荐的通用定义**——他们后来改用
> "vertex sides 之间的无向边" 定义（见 bioRxiv 2022.01.20.477068 Section 2）。
> 但本文档精确描述了 **BCALM 2 到底怎么表示 cdBG**，仍是理解其输出的权威来源。

### 4.1 边是 5 元组，有镜像约束

边 `e = (from, to, fromSign, toSign, label)`，`fromSign/toSign ∈ {+,-}`。给定两节点
x、y 共有 8 种边型，但按**镜像**分为 4 种连接类型（镜像边成对存在）：

| 镜像类型 | 边 | 镜像边 |
|---|---|---|
| 1 | (x,y,+,+) | (y,x,-,-) |
| 2 | (x,y,-,-) | (y,x,+,+) |
| 3 | (x,y,+,-) | (y,x,+,-) |
| 4 | (x,y,-,+) | (y,x,-,+) |

**自环例外**（e.from = e.to）：类型 3/4 的镜像边完全相同，只能保留一条
（self-mirror）。作者吐槽"这个特例耗费了人类数千小时的调试时间"。奇数长度的
字符串不可能等于自己的反向互补，因此 self-mirror 的 overlap 长度必为偶数。

### 4.2 在 DNA 上的解读（spelling rule）

节点代表一对字符串（label 及其反向互补 rc）。符号 + 用 label，− 用 rc(label)：

| fromSign | toSign | overlap |
|---|---|---|
| + | + | label 后缀 = 另一 label 前缀 |
| + | − | label 后缀 = rc(另一 label) 前缀 |
| − | + | rc(label) 后缀 = 另一 label 前缀 |
| − | − | rc(label) 后缀 = rc(另一 label) 前缀 |

两个字符串可按 overlap 拼接成更长的字符串（spelling rule）。两个方向相反的
overlap 恰好构成类型 4 的镜像对——所以 overlap 天然满足双向图的镜像约束。

### 4.3 双向 de Bruijn 图与 unitig/walk

- **节点**：所有去重后的 canonical k-mer（k-mer 与 rc 只占一个节点）。
- **边**：所有长度为 k-1 的 overlap。
- **walk**：边序列满足 `e_i.to = e_{i+1}.from` 且内部顶点处 `e_i.toSign =
  e_{i+1}.fromSign`（方向上的"sign 连续"）。单个顶点也是 walk。walk 按 spelling
  rule 拼出字符串，其**镜像 walk** 拼出反向互补。
- **unitig**：单顶点或满足条件的路径——内部顶点只与两条边（及各自镜像）关联；
  两端点是"唯一的出/入边"。**最大 unitig** 无法再向两端扩展。
- **compact 图**：把每个最大 unitig 及其镜像合并成一个顶点，边表示所有 k-1
  overlap。性质（作者未给形式证明）：最大 unitig 是图的一个顶点分解；walk 要么
  完整经过 unitig 顶点，要么只以 unitig 的前缀/后缀开头或结尾。

## 5. 输出格式

### 5.1 FASTA（默认，`<prefix>.unitigs.fa`）

```
><id> LN:i:<len> KC:i:<total_ab> km:f:<mean_ab> L:<+/->:<other_id>:<+/-> [..]
```

- `LN`：unitig 长度；`KC`/`km`：unitig 内 k-mer 总丰度 / 平均丰度。
- 每个连接边一个 `L:` 条目，记在**边的 from 节点**头上。`L:+:y:+` 表示
  forward-forward 出边；`L:+:y:-` 表示 forward-reverse。入边以"当前节点的 RC 的
  出边"形式编码，如 `L:-:x:+`。
- **BCALM 2 会记录所有边**（不只记每条镜像的一条）。
- `-all-abundance-counts` 时输出 `ab:Z:<ab_0> <ab_1> ... <ab_(len-k)>`（unitig 内
  每个 k-mer 的丰度向量）。

### 5.2 GFA（`scripts/convertToGFA.py`）

用法 `python convertToGFA.py in.fa out.gfa k`：
- 写 `H VN:Z:1.0 ks:i:<k>` 头；每个 unitig 一个 `S` Segment；每个 `L:` 一条 `L` Link，
  overlap 长度固定为 `k-1`（`<k-1>M`）。
- `-s/--single-directed`：只输出两节点间**一条**边（`name < b[2]`；同名时
  `-(b[1]==b[3]=='-')`），避免输出整个反对称图；否则输出全部 `L:`。
- 旧版 `MA=` 标签兼容处理为 `MA:f:`。

## 6. 参数语义（README）

| 参数 | 含义 |
|---|---|
| `-in` | 输入（单个 FASTA/FASTQ，或文件列表）；`-out` 输出前缀（默认 = 输入 basename） |
| `-kmer-size` | k-mer 长度（节点长度） |
| `-abundance-min` | 丰度阈值 X：seen **严格小于** X 次的 k-mer 被过滤（典型去测序错误） |
| `-minimizer-size` | minimizer 长度（分桶粒度，示例用 5-8） |
| `-minimizer-type` | 0 = 字典序 minimizer；1 = frequency-based（默认 1） |

**更大的 k**：源码编译时用 `cmake -DKSIZE_LIST="32 64 ... 320"` 指定 k 的
倍数（只能 32 的倍数，必须含 32），运行时可用到列表最大值；中间值用更小的
模板实例加速。`KSIZE_LIST="32 320"` 时 k>32 的速度等同 k=320。

**中间产物**：`.h5`（或 `_gatb/`，k-mer 计数）与 `*glue*`（待拼接压缩序列），
执行后安全删除；真正输出只有 `.unitigs.fa`。

## 7. 脚本与工具

- **`pufferize.py` / `split_unitigs.py`**：在参考基因组**端点 k-mer**（集合 B/E）
  处把 unitig 切开，使每个 unitig 首 k-mer ∈ B、尾 k-mer ∈ E，适配 pufferfish /
  使 unitig 与参考端点对齐。两者逻辑相同，前者输出 GFA + 重建 `P` path、后者输出
  FASTA。切割规则：遇"起点/终点 k-mer"在其处（含/不含该 k-mer）切分。
- **`abundance_stats.py`**：读 unitigs FASTA，统计每个 `km:f:` 平均丰度对应的
  unitig 个数与总长度直方图。
- **`unitigEvaluator.cpp`**（来自 BRAW）：评估"reads/参考 k-mer 是否都在 unitigs
  中"——按 k-mer 哈希分 `2^n` 趟 + 1024 桶哈希表 + OpenMP，输出 TP/FP/FN、
  错误 k-mer 率（*10000）与缺失 k-mer 率。用 OpenMP lock 数组做桶级同步。
- **`thirdparty/`**：`ThreadPool.h` 是 progschj 线程池改版（任务回调带
  `thread_id`）；`concurrentqueue.h` 是 moodycamel 无锁队列；另有若干 lock 型
  队列/vector。这些是 GATB 压缩/glue 阶段并发的底层原语。
- **example/`circular_unitigs_unittests`**：专门针对"BCALM2 不拼装环状 contig"的
  长期 issue——`test1.fa`（完美环状 unitig）、`test2.fa`（polyA 尾）、`test3.fa`
  （加"随机垃圾"使 k-mer 落入同桶）的回归。

## 8. 与 pgr 的关联

- **pgr 现状**：`libs/kmer`（`KmerTable`：canonical 2-bit u128 key、精确计数、
  radix sort、rayon 并行）已具备精确 k-mer 计数能力；`libs/chain`/`libs/paf` 有
  图与区间数据结构。目前**没有** cdBG / unitig 压缩功能。
- **可借鉴点**：
  1. **双向图表示法**（§4）：BCALM 2 的 `L:+:-` 边编码 + spelling rule + unitig
     定义，是表示"双链 de Bruijn 图 + 压缩"的清晰模板——若 pgr 未来做
     `pgr kmer` 层面的 unitig 或图输出，可直接照搬其边语义。
  2. **unitig 的数学定义**：内部顶点 degree ≤ 2（含镜像）、端点唯一出/入边、
     最大 unitig 是顶点分解——这是判断"某条压缩路径是否合法 unitig"的判据。
  3. **三段流水线分工**（§3）：分桶压缩 → UF 拼接 → 重算边，边界清晰；每段可用
     不同并行策略（桶内 `graph3` 线程池、UF 无锁并查集、LinkTigs 分趟磁盘哈希）。
  4. **桶内压缩算法 `graph3`**（§3.3）：把 k-mer 序列按首/尾 (k-1)-mer 建左右索引、
     排序后双指针合并 overlap——这是"给定一批 k-mer 求其 unitig"的直接可移植实现
     （pgr 的 `KmerTable` 精确计数恰好提供这批 k-mer）。
  5. **无锁并行并查集 `unionFind`**（§3.4）：路径压缩 + 按秩合并 + CAS，低 32 位
     parent / 高 32 位 rank 打包，可复用到 pgr 需要"按等价类合并"的并行场景。
  6. **端点哈希 + BooPHF MPHF + UF**（§3.2）：用哈希而非原 k-mer 做并查集键以省
     内存（容忍碰撞），并查集结果镜像成 32 位类号再按类分片拼链。
  7. **FASTA 头嵌边信息**（`L:`）+ 用独立脚本转 GFA（`convertToGFA.py`）的分工
     模式：bcalm 本体只输出带边注释的 FASTA，GFA 交给后处理。pgr 若加 GFA 输出
     可参考此分层。
  8. **minimizer 分桶 + traveller k-mer 落盘**（§3.1）：跨桶边落盘延迟处理，避免
     内存膨胀——内存-磁盘权衡的典型设计。
- **pgr 的差异优势**：`KmerTable` 是精确计数（无近似/无哈希碰撞），内存中建表；
  BCALM 2 依赖外部排序（dsk）+ minimizer 分桶以在低内存跑超大基因组。对 pgr
  的定位（单机、精确、内存友好场景），unitig 压缩更可能在 `fas`/`kmer` 层面的
  小图场景落地——此时可复用 `graph3` 的桶内压缩思路（§3.3），而无需完整移植
  dsk 外部排序与 traveller 落盘这些大内存场景才需要的机制。

## 9. 局限

- 依赖 GATB（大型 C++ 框架，C++11、OpenMP），构建链较重；`gatb-core` 是 git
  子模块，需 `--recursive` 或 `submodule update` 拉取。
- 环状 unitig 的拼装曾是长期未解 issue（见 example/circular_unitigs_unittests）：
  bglue 用 `determine_order_sequences` 的 `expect_circular` 兜底处理环，拼装后
  用"砍掉最后一个碱基"的方式避免自闭合（`glue_sequences` 里 `res_seq.substr(1)`）。
- UF 键是 k-mer 哈希（可碰撞）；UF 类号 32 位，序列总数理论上限
  `2^(32-1)`（`bglue_algo.cpp` 文件头注释），超过会出错。
- **`link_tigs` 不支持 k<4**（`LinkTigs.cpp` 检查 `kmerSize < 4` 直接报错退出，
  作者说"tiny dBG 请用 Python 构造"）；`is_in_pass` 的分趟哈希也不支持超过
  8 趟（`normalized_smallmer` 只处理 4 个碱基）。
- unitig 方向跨 run 不保证稳定；输出是 FASTA（丢质量信息）。
- 代码中多处 `#if 0` / `#ifdef` 遗留实验路径（`BINSEQ`、`pre_tip_cleaning`、
  `to_compact`、旧版 unordered_map 并查集），作者自己也标注了多处
  "could be further optimized"、"didn't test it yet"。
- 单线程旧版是独立项目（Malfoy/bcalm，即 BCALM 1），本仓库是 BCALM 2 多线程版。

---

*参考来源: 本项目源码 `bcalm/`（src/ + gatb-core/gatb-core/src/gatb/bcalm2/
bcalm_algo.cpp、bglue_algo.cpp、ograph.cpp、unionFind.hpp + debruijn/impl/LinkTigs.cpp、
UnitigsConstructionAlgorithm.cpp + bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md
+ scripts/ + README.md + CMakeLists.txt + example/ + thirdparty/）*
