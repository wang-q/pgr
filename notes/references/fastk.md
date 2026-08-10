# FastK: 高保真 k-mer 计数器

> 整理于 2026-06，更新于 2026-08-11，源自对仓库内 `FASTK-master/` 源码的分析
> （= FASTK-1.2 快照，2025-09-13 下载，README 标注 Current: April 18, 2021；
> **不是上游当前 master**，见 [[../design/kmer.md|kmer.md]] §2.3 版本核对）。
> 目的：理解 FastK 的 Super-mer + Minimizer k-mer 计数策略，为 pgr 的
> `rept s-kmer`/`e-kmer` 与 `kmer` 命令组（原生 `libs/kmer`）提供参考。
> **pgr 已不再依赖外部 FastK/Profex**（2026-08-09 迁移完成），本笔记的
> 借鉴价值在于：理解其算法机制（minimizer 分桶、两阶段排序、profile 编码），
> 以及对照 pgr 原生实现取舍了什么（见 §4）。

## 1. 简介
FastK 是一个专为处理高质量 DNA 组装数据（如 Illumina 或 PacBio HiFi 数据）而优化的 k-mer 计数工具。它采用了一种新颖的基于 Minimizer 的分发方案，能够利用磁盘存储来处理任意大小的数据集。

其核心优势在于能够直接生成序列的 k-mer 计数档案（Profile），并且在处理低错误率（1% 或更低）数据时，通过两阶段的“先 Super-mer 后加权 k-mer”排序策略，实现了极高的速度。

## 2. 核心概念 (Key Concepts)

### 1. K-mer 与 Canonical K-mer
*   **K-mer**: 长度为 k 的 DNA 序列片段（FastK 支持 k >= 5，默认为 40）。
    注：CLI 层对 k 只做 `ARG_POSITIVE` 校验（k > 0），**没有显式的 k >= 5 下限**；
    "5" 是 split.c 的 seed minimizer 长度（`#define MIN_LEN 5`）。但 `Determine_Scheme`
    末尾有硬校验 `if (KMER < PAD_LEN) Clean_Exit`（`PAD_LEN = MIN_LEN + PAD >= 5`），
    因此实际运行时 k 必须 ≥ PAD_LEN（至少 5），太小会直接报错退出。
*   **Canonical K-mer (规范 K-mer)**: 为了处理测序读取方向未知的问题，FastK 将一个 k-mer 及其反向互补序列（Watson-Crick complement）视为同一个 k-mer。在两者中，字典序较小的那个被称为“规范 k-mer”。FastK 的统计表（Table）只记录规范 k-mer。

### 2. Super-mer (超 k-mer)
*   背景: 在传统的 k-mer 计数中，如果一条读长（Read）有 150bp，K=40，那么它包含 111 个 k-mer。如果直接把这 111 个 k-mer 全部拆散单独存下来排序，数据量会膨胀得非常大（111 * 40 bytes），而且丢失了它们原本相邻的信息。
*   定义: Super-mer 是原序列中一段连续的、共享同一个 Minimizer 的子序列。
    *   当 FastK 扫描序列时，只要相邻的 k-mer 的 Minimizer 没有变（或者新出现的 Minimizer 依然是当前这一段里最小的），就把它们连在一起，形成一个 Super-mer。
    *   直到 Minimizer 发生变化，不得不划分到另一个桶时，才切断。
*   例子:
    *   假设序列是 `ABCDE`，K=3。
    *   包含的 k-mer 有: `ABC`, `BCD`, `CDE`。
    *   如果这三个 k-mer 算出来的 Minimizer 都是 `B`。
    *   传统做法: 存 3 个条目 `ABC`, `BCD`, `CDE`。
    *   FastK 做法: 存 1 个条目 `ABCDE` (Super-mer)。
*   优势:
    *   压缩数据: 极大减少了需要写入磁盘和排序的条目数量（通常减少 10-50 倍）。
       （注：本 checkout 的 README/代码中无此数字出处，系估计值，可能来自论文或其他版本。）
    *   提升速度: 排序 1 个长条目比排序 10 个短条目要快得多。在最后阶段，程序再从 Super-mer 中还原出具体的 k-mer 进行计数。

### 3. Minimizer (最小标识符)
*   定义: Minimizer 是 k-mer 序列内部的一个特定的、较短的子序列（m-mer，m < k）。通常选择字典序最小的那个 m-mer。
*   分发策略: FastK 根据 k-mer 中包含的 Minimizer 来决定将其放入哪个存储桶（Bucket）。
*   核心作用:
    *   保证同类归并: 如果两个 k-mer 的序列完全相同，它们必然拥有相同的 Minimizer。因此，它们一定会被分发到同一个桶中。这就像把所有姓“李”的书都扔进“L”号箱子，不管这书是从哪里来的，只要是同一本书，它一定在“L”箱子里。
    *   并行独立性: 因为所有相同的 k-mer 都在同一个桶里，我们在统计“L”箱子时，完全不需要去问“Z”箱子有没有漏掉的。这意味着不同的桶可以完全独立地由不同线程并行处理，互不干扰。
    *   内存控制: 无论数据集多大（比如 1TB），我们都可以通过增加桶的数量（比如分成 1000 个桶），让每个桶只有 1GB，从而可以轻松读入内存进行快速排序。

**实现细节（比"固定 m-mer"更精巧）**：FastK 的 minimizer 不是简单的固定长度 m-mer，而是
一个 **core prefix trie（核心前缀树）+ 自适应 padding** 的方案（`Determine_Scheme`，
`split.c`）：

*   以深度 5 的全 4 叉树为种子（`MIN_LEN=5`，即 5-mer，共 `MIN_TOT=4^5=1024` 个叶节点）。
*   先在输入前 1 Gbp 上（`Get_First_Block(io, 1000000000)`）统计每个 core prefix 的出现
    频次，按目标分块数 `npieces=2*NPARTS` 动态给每个叶节点加 `PAD` 位（每次不够再 +2，
    `PAD += 2`），把过热的 minimizer 继续往下细分，直到每个 core prefix 的频次都
    ≤ 阈值 `kthresh = ktot/npieces`（或提升 < 2% 时停止）。
    `PAD_LEN = MIN_LEN + PAD`；**`MAX_SUPER = KMER - PAD_LEN + 1`**（代码里
    `MAX_SUPER = KMER - PAD_L1`，`PAD_L1 = PAD_LEN-1`），即超 k-mer 的最大长度由
    padding 反推——分发时 `force = (p-m >= MAX_SUPER)` 强制切断超长的 super-mer。
*   迭代停止条件（三选一）：① `o == Min_States`（无需再分裂，trie 已"core"）；
    ② `PAD > 0 && last_max < 1.02*max_count`（细化收益 < 2%，此时**重新计算**
    `NPARTS = ktot/max_count+1` 并 `NPARTS /= 2`，桶数按实际数据回退调整）；
    ③ `PAD_LEN >= KMER-1`（padding 已接近 k）。
*   `assign_pieces` 再把 core prefix 尽可能均匀地装箱到 `NPARTS` 个桶（先按频次降序，
    加权随机贪心）。`Min_Part[]` 数组即编码了这棵 trie：`>=0` 为桶号，`<0` 为
    进入子树的状态偏移。
*   因此同一个 minimizer 可能对应不同桶；分发时沿 trie 逐位下钻定位桶（`b<0` 时
    `b = Min_Part[((mc >> y) & 0x3) - b]`）。
*   **碱基编码是按频率排序的**（`Tran[]`）：对训练集统计 A/C/G/T 频次，把最频繁的映射
    为 0、次之为 1……这样 minimizer 值小的区域天然承载更多数据，有助于前缀树平衡。
    这不同于固定 A=0/C=1/G=2/T=3 的编码。计数阶段（`count.c`）的 canonical 比较用的是
    另一套固定 `Dran/Fran`（a=0,c=1,g=2,t=3），两套编码职责不同：`Tran` 只管分桶路由，
    `Dran/Fran` 管 k-mer 的 2-bit 编码与反向互补。

**source quirk**：`refine_tree` 里 `PAD += 2`（仅 `FIND_MTHRESH` 未定义时 +1），且
`last_max < 1.02*max_count` 时停止细化并可能改变 `NPARTS`——桶数在训练阶段可能被
实际数据规模回退调整。末尾有硬校验 `if (KMER < PAD_LEN) Clean_Exit`（"K-mer must be
at least PAD_LEN"），k 太小直接报错退出。

> **注意：与 `pgr` (Sketching) 的区别**
> FastK 中的 Minimizer 与 `pgr/src/libs/hash.rs` 中用于 MinHash Sketch 的 Minimizer 虽然原理相同，但用途截然不同：
> *   FastK (无损路由): 针对每一个 k-mer，找到其内部的 m-mer 标签，用于决定去哪个桶。不丢弃任何数据，目的是全量统计。
> *   PGR Sketch (有损采样): 在长序列的滑动窗口中选出一个 k-mer 代表该窗口。丢弃绝大部分数据，目的是生成稀疏指纹（Signature）用于快速比对。

### 4. 输出文件格式
FastK 生成以下几种核心文件（均为二进制；最终结果写入 `-N` 指定的目录，中间
分桶/排序文件落在 `-P` 指定的 sort 目录）：
*   **直方图 (.hist)**: 记录了每种频次（1次, 2次...）出现的 k-mer 数量。最大计数限制为 32,767。
    *   二进制布局（`count.c` 写、`libfastk.c::Load_Histogram` 读）：
        固定头 `int32 k | int32 low(=1) | int32 high(=32767) | int64 ilowcnt | int64 max_inst`
        + `int64 hist[1..=32767]`（共 32767 个频次 bin）。`max_inst` 是 ≥ high 的
        k-mer 实例总数（封顶累计）。pgr `libs/kmer/hist.rs` **字节级复刻**了此布局
        （见 §4 借鉴）。
*   **K-mer 表 (.ktab)**: 一个排序的列表，包含数据集中所有（或满足特定阈值）的规范 k-mer 及其计数。
    *   结构 = **stub 文件 + 分片**（`table.c::Merge_Tables`）：
        `<root>.ktab` stub 含 `int32 k | int32 NTHREADS | int32 cutoff | int32 IDX_BYTES`
        + 前缀索引 `int64 pindex[2^(8*IDX_BYTES)]`（偏移表）；数据在隐藏分片
        `.<root>.ktab.#`，每条记录 = k-mer 字节 + `uint16 count`。这就是
        `libfastk.c::_Kmer_Stream` 的"前缀压缩 + 索引"结构（`ibyte` 前缀 →
        `index` 偏移表 + `inverse_index`），供 O(1) 前缀定位 + 桶内小范围查找。
    *   `IDX_BYTES` 按总条目数自适应（`count.c`：>2^26 用 3、≥2^16 用 2、否则 1）。
*   **档案 (.prof)**: (可选) 数据集中每条序列的 k-mer 计数概况。这是一个压缩格式，不仅记录了 k-mer 的出现，还按原序列顺序保留了位置信息。
    *   结构 = stub + 分片（`merge.c::Merge_Profiles`）：`<root>.prof` stub 含
        `int32 k | int32 ITHREADS`；数据在 `.<root>.pidx.N`（A-file，read id →
        偏移索引）+ `.<root>.prof.N`（D-file，RLE 压缩的逐 super-mer profile）。
    *   profile 编码（`merge.c`）用 **delta（前向差分）+ RLE**：每条 read 首值
        绝对编码，后续值与前一值做差，差为 0 则 RLE 计数（≤63），小差（±31）
        单字节 0x40 编码、大差 0x80 两字节；README 称实际约 **4.7 bits/base**。

## 2.5 CLI 参数与内存/并行设计（`FastK.c`）

**用法**：`FastK [-k<int(40)>] [-t[<int(1)>]] [-p[:<table>[.ktab]]] [-c] [-bc<int>]
[-v] [-N<path_name>] [-P<dir($TMPDIR)>] [-M<int(12)>] [-T<int(4)>] <source>...`

| 参数 | 含义 | 默认 |
| :--- | :--- | :--- |
| `-k` | k-mer 长度 | 40 |
| `-t[<level>]` | 建 k-mer 表，只保留计数 ≥ level 的 k-mer | 无参数=cutoff 1 全量 |
| `-p[:<table>]` | 生成 profile；带表则生成**相对 profile**（值=表内 count，缺省 0） | 关 |
| `-c` | 先把输入做 **homopolymer 压缩**（连续同碱基折叠为单个）再计数 | 关 |
| `-bc` | 忽略每条 read 开头的 N 个碱基（如 barcode） | 0 |
| `-v` | 详细输出各阶段统计 | 关 |
| `-N` | 输出目录 + 根名前缀 | 取自首个输入 |
| `-P` | 块级排序的外部文件目录 | `$TMPDIR` 或 `/tmp` |
| `-M` | 排序阶段可用内存（GB） | 12 |
| `-T` | 线程数 | 4 |

**内存/并行设计**：
*   **NPARTS（分桶数）自适应**：先用前 1 Gbp 块估计总 k-mer 规模
    `gsize = (block->totlen - k*nreads) * ratio * rsize`，再
    `NPARTS = ceil(gsize / SORT_MEMORY)`（`FastK.c`）。每个桶的排序数据量由
    `-M` 控制，保证能进内存；磁盘只做中转。
*   **两级文件扇出**：Phase 1 把 super-mer 写到 `NPARTS × ITHREADS` 个 `.T*` 文件
    （每个线程 t 对每个桶 n 一个），再按桶并行处理；`RLIMIT_NOFILE` 会校验
    需同时打开的 `(NPARTS+3)*NTHREADS` 个文件描述符。
*   **`-p:<table>` 要求 k 一致**：`FastK.c` 校验 `PRO_TABLE->kmer != KMER` 直接报错；
    且 `-p:<table>` 时 `-t` 被忽略、默认直方图也不再输出（相对 profile 不做直方图）。
*   全部阶段用 pthread 并行；`-T` 同时是桶内排序线程数与 profile 分片数
  （`NTHREADS`），输入读取另用 `ITHREADS`。

## 3. 处理流程 (Processing Workflow)

FastK 的内部处理逻辑主要分为四个阶段（Phase）：

### 第一阶段：分片与分发 (Phase 1: Partitioning)
*   输入扫描: 程序首先扫描输入数据集的前 1GB 数据。
*   方案确定: 基于这部分数据，计算 Minimizer 分布，确定如何将 Super-mer 均衡地分发到临时桶中。
*   全量分发: 扫描整个数据集，计算 Super-mer，并根据 Minimizer 方案将它们写入磁盘上的不同临时文件（Buckets）。
*   *代码对应*: `split.c`

### 第二阶段：排序与计数 (Phase 2: Sorting & Counting)
*   并行处理: 对每个临时桶（partition）并行执行操作。
*   两级排序（`count.c::Sorting`）:
    1.  **Super-mer 排序**：`supermer_list_thread` 解包每个 `.T*` 文件的位压缩
        super-mer 到数组，`Supermer_Sort` 排序（LSD/MSD 基数排序）。
    2.  **加权 k-mer 排序**：`kmer_list_thread` 对排序后的 super-mer 计数，
        把每个唯一 super-mer 内**所有 k-mer 以其出现次数加权**（同一 super-mer
        出现 ct 次，则其内每个 k-mer 权重 = ct，写入 `uint16`）放入数组，
        `Weighted_Kmer_Sort` 再排序。
*   **canonical 判定**：`kmer_list_thread` 对每个 k-mer 的正向与反向互补
  （用 `Comp` 表逐 2-bit 互补实现反转）逐字节比较，取字典序较小者
  （`kb<hb` 取正向，否则取反向）。
*   **count 上限 32767**：`ct >= 0x8000` 时溢出部分记入 `overflow`，count 封顶
  `0x7fff`（与 .hist 的 high 一致）。
*   **统计生成**：`Weighted_Kmer_Sort` 过程中累积 k-mer 频次直方图，最后写
  `<root>.hist`（`PRO_TABLE==NULL` 时才输出）。
*   按需输出：`-t` → `table_write_thread` 写 `L*` 分片；`-p` → 逆序倒排两次
  （`cmer_list`/`cmer_merge` + LSD sort + `profile_list`/`profile_write`）生成
  `P*` profile 分片。
*   *代码对应*: `count.c`, `LSDsort.c`, `MSDsort.c`

### 第三阶段：表格合并 (Phase 3: Table Merging)
*   合并: 将第二阶段各线程/桶生成的排序好的 k-mer 片段，按字典序合并成一个单一的、全局有序的 `.ktab` 文件。
*   实现（`table.c::Merge_Tables`）：对每个线程，把 NPARTS 个 `L*` 分片做
  **k-way 堆合并**（小根堆，按 `KMER_BYTES` 前缀比较），逐条写出；同时统计各
  `IDX_BYTES` 前缀的出现次数累积成前缀偏移索引，连同 stub 头一起回填
  `.ktab`。合并后删除 `L*` 分片。
*   *代码对应*: `table.c`

### 第四阶段：档案合并 (Phase 4: Profile Merging)
*   合并: (仅当启用 `-p` 选项时) 将分布在各处的 Profile 片段合并成最终的 `.prof` 文件。
*   实现（`merge.c::Merge_Profiles`）：每个线程把 NPARTS×NPANELS 个 `P*` 分片
  按 super-mer id 归并（分块 `PAN_SIZE=1024*NPARTS`），对同一 read 的相邻
  super-mer profile 做 **delta + RLE 压缩**，写出 D-file（profile 数据）与
  A-file（read id → 偏移索引）。
*   压缩: Profile 数据采用 delta+RLE 编码压缩（约 4.7 bits/base），以节省空间。
*   *代码对应*: `merge.c`

## 4. 对 pgr 的启示

> **现状（2026-08-09 起）**：pgr 已用原生 `src/libs/kmer/` 替换外部 FastK/Profex，
> `pgr rept s-kmer` / `e-kmer` 不再依赖外部工具；并新增 `pgr kmer` 命令组
> （table/profile/hist/gc/qhist/qcheck/gsize）。完整设计与取舍见
> [[../design/kmer.md|kmer.md]]（§9 功能对照、§10 三种格式）。

1. **Super-mer + Minimizer 策略（未复刻）**：FastK 通过 Super-mer（共享同一
   Minimizer 的连续 k-mer 段）将数据量压缩 10-50 倍，再按 Minimizer 分桶并行排序。
   这套"先聚合后排序 + 磁盘分桶"是为 **TB 级数据**设计的。pgr 原生实现
   **主动放弃**了 super-mer 与磁盘分桶（`kmer.md` §3.2）：直接全量收集 canonical
   u128 key → 全局 `radix_sort_u128_par` 排序 → 一趟分组计数。细菌/真菌级
   （~5 Mb–50 Mb）内存可承受（u128+u32 ≈ 20 B/唯一 k-mer，50 Mb 真菌 ≈ 1 GB），
   用完即释放。**判断**：pgr 场景（单基因组/单库，≤ 数百 Mb）用不上 TB 级
   分桶，简化是对的；若未来目标是 Gb 级基因组，才需在 `build_table` 内加分块。

2. **Minimizer 的两种用途**：FastK 的 Minimizer 是**无损路由**（每个 k-mer 都被
   分配到对应桶，不丢弃数据，用于全量统计）；pgr 的 `src/libs/hash.rs` 中用于
   MinHash Sketch 的 Minimizer 是**有损采样**（在滑动窗口中选出代表 k-mer，丢弃
   大部分数据，用于生成稀疏指纹）。两者原理相同但用途截然不同，`pgr dist seq`
   用的是后者（详见 §2 核心概念中的对比注释）。**Minimizer 的"自适应 padding
   前缀树 + 频率排序碱基映射"（`Determine_Scheme`）pgr 未复刻**——那是为了把
   数据均衡散到 NPARTS 个磁盘桶；pgr 的 `KmerTable` 走单块内存排序，无此需求。

3. **两阶段排序策略（借鉴了"排序合并"这一思路）**：FastK 的"先 Super-mer 排序，
   再 k-mer 加权排序"在低错误率数据上极快。pgr 原生实现**采用了其"排序合并
   替代哈希/查表"的核心哲学**：`profile.rs` 曾用 `partition_point` 二分查表，
   profiling 实测 `table_profiles` 占 78.5%（73 MB 表 DRAM 随机访问 cache-miss
   主导），改用 FastK 式的**收集全部窗口 key → 排序 → 与 `table.keys` 线性归并**，
   self_profiles 5.2×、relative_profiles 5.4×、`rept s-kmer` 整命令 3.4× 提速
   （`kmer.md` §3.5，基准见 `notes/benchmarks/bench-profile-hotspots.md`）。

4. **直方图 `.hist` 字节兼容（已复刻）**：`count.c` 写 `.hist` 的固定布局
   （`k | low=1 | high=32767 | ilowcnt | max_inst | hist[32767]`）被 pgr
   `libs/kmer/hist.rs` 逐字节复刻，外部 Histex/KatGC/GenomeScope 可直接读 pgr
   输出（实测与 FastK 自产逐行一致）。这是"单一小文件、外部有消费者"格式做
   兼容的范例（成本 ~50 行）。

5. **`.ktab`/`.prof` 格式不做兼容（决策）**：两者都是"stub + 隐藏分片"的多文件
   布局（`.ktab` 含前缀索引、`.prof` 含 delta+RLE 与 read-id 索引），pgr 分别用
   单文件 `.pkt`（紧凑 packed key + u32 count）和 `.pkp`（header + raw u16）替代。
   原因：pgr 无外部消费者、不喜欢分片；`KmerTable` 查表也不需要 `.ktab` 的
   前缀压缩索引（内存排序数组 + 二分/归并即可）。**借鉴点**：`.ktab` 的
   "前缀压缩 + 偏移索引"结构与 pgr `pgi`/PAF 前缀索引思想同源，未来若需大表
   O(1) 前缀定位可参考，当前无此需求。

6. **外部工具依赖已消除**：`pgr rept s-kmer` / `e-kmer` 及 `pgr kmer` 全链路
   不再调用 FastK/Profex/Histex；`--keep-index` 的旧 `.ktab` 缓存升级后自动重建为
   `.pkt`。理解 FastK 内部机制仍有助于在跨工具对照（如与 KMC3/Jellyfish 语义
   核对）时定位差异。

7. **`-c` homopolymer 压缩（未做）**：PacBio/HiFi 专用（homopolymer 错误率是
   其他错误的 5 倍，`-c` 折叠后 hoco k-mer 错误率降 5 倍）。pgr 无此选项；
   若未来要支持 HiFi 组装质量评估（MerquryFK 场景）可低成本补上（`kmer.md` §9
   列为剩余缺口优先级 2）。

---
*参考来源: [FastK GitHub Repository](https://github.com/thegenemyers/FASTK)*
