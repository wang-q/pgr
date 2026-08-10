# syng: Syncmer 图与序列集合的高效表示

> 整理于 2026-08，源自对 `syng-main/` 目录源码的分析。目的：理解 Richard Durbin 的 syncmer（特别是 closed syncmer）实现，为 pgr 用 syncmer 替代当前 minimizer 提供算法与移植参考。本文聚焦 syncmer 部分，GBWT/图构建等内容不在 pgr 当前需求范围内，仅作背景介绍。

## 1. 简介

`syng`（**Syn**cmer **G**raph）是 Richard Durbin（剑桥）开发的序列图工具，核心思想是：把任意 DNA 序列表示为一组 *syncmer* 上的 path。syncmer 是所有 k-mer 的一个子集，**保证对任何 DNA 序列提供稀疏但完整的覆盖**（平均深度约 2×）。这与 minimizer 的"有损采样"形成关键对比。

- **syncmer 定义**：本文采用 Edgar (2021) 的 *closed syncmer*——一个长度为 `w+k` 的窗口被称为 syncmer，当且仅当其内部最小的 s-mer（长度为 `k` 的子串）出现在窗口的**首端或末端**。详见 [Edgar, 2021](https://peerj.com/articles/10805/)。
- **文件产物**：`.1khash`（syncmer 序列集合）、`.1path`（序列表示为 kmer 索引路径）、`.1gbwt`（GBWT 编码的图+路径）、`.1seq`（重建序列）。均基于 [ONEcode](https://github.com/thegenemyers/ONEcode) 二进制格式。
- **工具套件与最新进展**：仓库由 Makefile 构建出 `syng`（主程序）、`syngmap`（读段→图映射）、`syngpath2gbwt`（`.1path`→`.1gbwt`）、`syngstat`（ONEcode 统计）、`k31type`，外加 ONEcode 自带的 `ONEview` 与 `SEQUENCE_UTILITIES`（`seqconvert`/`seqstat`/`seqextract` 互转 `.1seq`）。2026-03 起 GBWT 用 `rskip.[ch]`（**行程长度编码跳表**，核心操作 O(logN)）重写，README 明确推荐**两步法**：先用 `syng` 产出 `.1path`，再用 `syngpath2gbwt` 生成 `.1gbwt`。各命令细节见 §3.4。
- **性能参考**：20× PacBio HiFi 的 935Mb 慈鲷基因组（~19Gbp）→ 1.05GB `.1khash` + 493Mb `.1gbwt`（(1023,32)-syncmer），MacBook Pro 上 62 秒。

> **范围说明**：pgr 不需要复现 syng 的图/GBWT/ONEcode 部分。我们只需要 `seqhash.[ch]` 里的 **syncmer 迭代器算法**，用来替换 `src/libs/hash.rs` 中基于 `minimizer_iter` crate 的 minimizer 采样。

## 2. 核心概念 (Key Concepts)

### 2.1 参数约定（syng 与 pgr 的差异，极易混淆）

syng 的参数命名与 pgr/minimizer 习惯**相反**，移植时必须注意：

| 符号 | syng 含义 | pgr/minimizer 习惯 | 本文澄清 |
| :--- | :--- | :--- | :--- |
| `k` | **s-mer 长度**（用于哈希的小子串） | k-mer 长度（被采样的串） | syng 的 `k` 是"哈希粒度" |
| `w` | 窗口内 s-mer 个数的参数（**注意**：syng 主程序经 `+1` 后实际窗口含 `w+1` 个 s-mer，见下行） | 窗口内 k-mer 个数 | 语义一致（pgr 按 w 个） |
| `w + k - 1` | **syncmer 跨度**（w 个 s-mer 覆盖的碱基数） | — | 裸 `seqhash.c` 迭代器按 w 个 s-mer（跨度 `w+k-1`）；syng 主程序（`syng.c:432`）创建 Seqhash 时用 `w+1` 补偿（注释 "need the +1 here, awkwardly"），实际窗口含 `w+1` 个 s-mer、跨度 `w+k`，与 README 一致；pgr 移植按 w 个 s-mer（`w+k-1`） |
| `seed` | 哈希函数种子 | — | 一致 |

默认参数 `k=8, w=55`：syng 实际输出跨度 `63`（`w+1` 个 s-mer，`w+k`，`.1khash` 存储长度也是
`w+k`），pgr 移植按 w 个 s-mer、跨度 `62`（`w+k-1`）。README 中 `(1023,32)-syncmer` 指
`(w,k) = (1023,32)`。**下文如无特别说明，"syncmer 跨度"=`w+k-1`，"s-mer 长度"=`k`，窗口含
`w` 个 s-mer（pgr 移植语义）。**

> **最小序列长度差异**：裸 `syncmerIterator` 要求 `len >= w+k`（`seqhash.c:189`
> `if (len < sh->w + sh->k)`），比"w 个 s-mer 构成一个完整窗口"的理论最小值 `w+k-1` 多 1——
> 这会漏掉长度恰为 `w+k-1`（仅一个完整窗口）的序列中的 syncmer；syng 主程序因 `w+1` 补偿
> 实际要求 `len >= w+k+1`。pgr 实现按理论最小值判定（`closed_syncmers_stream` 在
> `hashes.len() < window` 时返回空），即单窗口序列若满足端点最小也会被采到。这是 pgr 比 syng
> 更宽松的边界处理。

### 2.2 Closed Syncmer 的定义与性质

设窗口含 `w` 个长度为 `k` 的 s-mer（位置 `0..w-1`），窗口跨度 `L = w + k - 1` 个碱基。设 `h_i` 为第 `i` 个 s-mer 的规范哈希（canonical，见 §3.1）。

- **closed syncmer 判定**：窗口为 syncmer ⟺ 窗口最小哈希值出现在位置 `0`（首位）或 `w-1`（末位）。syng 与 pgr 都用"最小值是否出现在端点"的**值判定**（`hash==min` / `x<=min`，见 `seqhash.c:201,206-209`），而非 argmin 位置。两者差异在**输出**：syng 的 `syncmerNext` 返回窗口首端 s-mer 的哈希（`hash[iStart]`），但 syng 主流程（`syng.c:73`）调用时 kmer 参数传 `0` 并不使用它——图路径节点身份来自位置 `pos` 处的 **`w+k` 长窗口序列**（`kmerHashFindThreadSafe`/`syncmerAdd` 按 canonical 存进 `KmerHash`，`.1khash` 里就是这些窗口序列）；s-mer 哈希输出仅在 `seqhash.c` 的 TEST 里打印。pgr 输出最小 s-mer 的 canonical 哈希（使序列与反向互补产生同一哈希集合，Mash/Jaccard 所需，见 `syncmer.rs` 的 `closed_syncmers_stream`/`syncmer_dna`）。并列时 syng 滚动循环先查末端（`x<=min`），pgr 先查首端（`==`）；但因 pgr 取的是最小值哈希，并列两端哈希相同，tie-break 不影响集合的链对称性。
- **密度保证**：相邻 syncmer 间隔有上界、无大 gap（密度数值见下条）。但**不保证序列首尾被覆盖**（syng 用 X/Y ends 补首尾；pgr 的 sketch 用途不依赖首尾覆盖）。这是 syncmer 相对 minimizer 的核心优势——采样位置由"端点最小"的几何约束决定，对 indel/重排局部化扰动。
- **密度**：平均约 `2/(w+1)` 的 s-mer 是 syncmer 的端点，对应 syncmer 在序列上的平均深度约 `2×`（每个位置平均被 2 个 syncmer 覆盖）。注：这是"窗口含 w 个 s-mer"的理论值；syng 实际用 `w+1` 个 s-mer，密度约 `2/(w+2)`，pgr 移植（w 个）为 `2/(w+1)`。
- **"closed" vs "open"**：closed 要求最小 s-mer 在两端；open syncmer 只要求在某个固定偏移。syng 只实现 closed。

### 2.3 三种采样器对比（同一基础设施）

`seqhash.[ch]` 在同一套 `Seqhash` + 滚动哈希基础设施上提供了四种迭代器，对比 pgr 场景：

| 采样器 | 判定规则 | 是否保证覆盖 | 密度 | pgr 现状 |
| :--- | :--- | :--- | :--- | :--- |
| **kmer** | 全部 k-mer | 是（满覆盖） | 1× | — |
| **minimizer** | 窗口内 hash 最小的 k-mer | 否（窗口间会跳） | ~1/w | `src/libs/hash.rs` 现用 |
| **mod-minimizer** | `hash % w == 0` 的 k-mer | 否 | 1/w | pgr 的 `"mod"` 选项 |
| **closed syncmer** | 最小 s-mer 在窗口首/末 | **是** | ~2/(w+1) | **本次移植目标** |

> 关键差异：minimizer 在每个窗口选最小值后"跳跃"到新窗口，可能跳过中间区域；syncmer 通过"端点最小"的几何约束保证每个位置都在某个 syncmer 窗口内，因此**对序列变更和 indel 更鲁棒**（插入/缺失只局部影响 syncmer 集合，而 minimizer 的跳跃会放大扰动）。这对 `pgr dist seq` 的 Mash/Jaccard 距离稳定性有直接价值。

## 3. syncmer 算法详解（移植核心）

### 3.1 哈希设计：乘加移位 + 滚动 + 规范化

`Seqhash` 结构体（`seqhash.h:15`）的核心字段：

```c
typedef struct {
  int seed ;            // 随机种子
  int k ;               // s-mer 长度（哈希粒度，< 32）
  int w ;               // 窗口内 s-mer 个数
  U64 mask ;            // 2*k 位的掩码
  int shift1, shift2 ;
  U64 factor1, factor2 ;// 乘加移位哈希的奇数因子
  U64 patternRC[4] ;    // 每个 base 的反向互补移位模式
} Seqhash ;
```

**关键设计**（`seqhash.c:17` `seqhashCreate`）：

1. **2-bit 编码**：碱基以 `0,1,2,3`（a,c,g,t）编码，一个 k-mer 占 `2k` 位（`k<32` 保证 fit 在 `u64`）。
2. **乘加移位哈希**（`seqhash.h:73`）：
   ```c
   static inline U64 kHash(Seqhash *sh, U64 k) {
     return ((k * sh->factor1) >> sh->shift1);
   }
   ```
   - `factor1` 是奇数（`| 0x01`），由 `srandom(seed)` 生成，`shift1 = 64 - 2*k`（`seqhashCreate` 同时算 `factor2`/`shift2`，但 `kHash` 只用 `factor1`/`shift1`）。
   - **seed 是实验调优来的**：`seqhash.c:301-365` 的 `#ifdef H_EXPLORE` 工具暴力遍历 seed
     `1..999999`（`for (h = 1; h < 1000000; ++h)`），对一段高度重复序列（poly-a/c/g/t +
     12 种二周期重复，每段 32nt × 16 段共 512bp）用全 k-mer 迭代器采样（每 32 个取一个），
     记录采样哈希的最小值并**最大化该最小值**（分布最均匀、无低值塌缩）的 seed。syng 默认
     `seed=7`（`syncmerset.c:17 syncmerParamsDefault`、`syng.c:432` 的 `seqhashCreate` 调用均用 7），很可能
     即由此工具搜出。pgr 的 `hash_factor` 用 splitmix64 由 seed 生成 factor
     （[syncmer.rs:159](file:///home/wangq/Scripts/pgr/src/libs/syncmer.rs#L159)），雪崩性质良好、
     任意 seed 都均匀，故无需此调优步骤。
   - 这是一种 fast universal hashing，比 Murmur/Fx 更轻量，且对短 k-mer 足够均匀。最终实现的哈希选择见 §5.1.3（DNA 用此 2-bit 乘加移位，蛋白用 `RapidHash` 字节哈希）。
3. **正反向同步滚动**（`seqhash.c:67` `advanceHashRC`）：
   ```c
   si->h    = ((si->h << 2) & sh->mask) | *s;           // 正向滚动
   si->hRC  = (si->hRC >> 2) | sh->patternRC[(int)*s];  // 反向互补滚动
   ```
   - `patternRC[i] = (3-i) << 2*(k-1)`：把互补碱基放到高位。
   - 每加入一个碱基，正向哈希左移 + 新碱基入低位；反向哈希右移 + 互补碱补充入高位。O(1) 滚动。
4. **规范化（canonical）**（`seqhash.c:57` `hashRC`）：取 `min(hashF, hashR)`，并记录 `isForward`。这等价于 pgr `seq_sketch` 里的 `.canonical()`。

### 3.2 closed syncmer 迭代器算法

这是本次移植的核心，源码在 `seqhash.c:186-250`，约 65 行 C。算法用一个长度 `w` 的**环形缓冲区** `hash[]` 存储当前窗口的 s-mer 哈希，`iStart` 指向窗口起点。

**初始化**（`syncmerIterator`）：

1. 计算前 `w` 个 s-mer 的哈希填入 `hash[0..w-1]`，求 `min`。
2. 若 `hash[0] == min` 或 `hash[w-1] == min`（最小 s-mer 已在端点），当前窗口即 syncmer，返回。
3. 否则向前滑动：每读入一个新 s-mer `x` 放到 `hash[iStart]`（覆盖最旧值，`iStart` 循环递增）：
   - 若 `x <= min` → `x` 是新窗口**末端**的最小 s-mer → syncmer，更新 `min = x`，返回；
   - 若 `hash[iStart]`（新窗口**起点**的值）`== min` → syncmer，返回。

**步进**（`syncmerNext`）：

1. 输出当前 syncmer 的 k-mer（`hash[iStart]`）、位置、链方向。
2. 若刚输出的位置就是 `min` 持有者，重置该槽为 `U64MAX` 并重新线性扫描 `hash[]` 求 `min`（注释提到可换堆，但 `w` 不大时线性扫描够用）。
3. 同初始化的逻辑向前滑动找下一个 syncmer。

**正确性直觉**：窗口每次前进一格，新窗口 = 旧窗口去掉首 s-mer、加入尾 s-mer。新窗口成为 syncmer ⟺ 新 min 在首或尾。分两种情况：
- 新加入的尾 s-mer `x` 是新 min → 端点（尾）；
- 旧窗口的 min 不在被移除的首位时，它仍在窗口内；若它在**新的首位**（即旧窗口的第二位）→ 端点（首）。

代码用 `x <= min`（而非 `<`）处理并列，保证边界情形的稳定性。

### 3.3 SyncmerSet / KmerHash 存储（参考用，pgr 暂不需要）

`syncmerset.[ch]` 在 syncmer 迭代器之上构建去重集合，`kmerhash.[ch]` 是底层哈希表。pgr 当前用 `rapidhash::RapidHashSet<u64>` 做去重，无需移植这套 ONEcode 持久化机制。但有两点设计值得记录：

- **1-based 索引 + 负号表示反向**（`kmerhash.h:19-22`）：`kmerHashAdd` 返回 `index > 0` 表示新增正向，`index < 0` 表示命中反向互补。pgr 的 `MinimizerInfo.strand` 字段可复用此思路。
- **canonical 定向**（`kmerhash.c:57` `isCanonical`）：k-mer 存储时统一取向（kmer < revcomp(kmer)），比较时只需正向比对。这与 pgr `seq_sketch` 的 `.canonical()` 一致。
- **同聚物过滤**：`syng.c` 建立 SyncmerSet 时预先把 4 条同聚物（poly-a/c/g/t，长度 `w+k`）插入
  `KmerHash`（索引 `1/2/-2/-1`），主流程跳过 `|sync| ≤ 2` 的命中（注释 "don't record
  poly-A/C/G/T"）——`.1khash` 因此不含同聚物 syncmer。
- **长度必须为奇数**：`kmerHashCreate` 要求存储长度 `len` 为奇数（`if (!(len & 0x01)) die`），故
  `w+k` 必须为奇数（默认 63，示例 (1023,32)→1055、(63,8)→71 均满足）。
- **KmerHash 哈希表机制**（`kmerhash.c`）：开地址哈希——`table[loc]` 存 1-based 索引，冲突探测步长
  `delta = hashDelta(pack, plen, dim)` 由 64 位字 X-OR 折叠而成，末尾 `v |= 1` **保证为奇数**（与
  `2^dim` 互质，故能遍历整个表）；条目按 2-bit 压缩进 `pack[]`（`plen=(len+31)>>5` 个 u64 每项，
  `psize=size*0.3`，`dim` 初值 20，`max` 逼近 `psize` 时 `doubleTable` 翻倍扩容并重排）。查询线程安全版
  `kmerHashFindThreadSafe` 用调用方提供的 `uBuf` 缓冲，供 syng/syngmap 的并行抽取直接使用。
- **SyncmerSet 计数与位置字段**（`syncmerset.h:27`）：`count`（跨所有输入的累计频次，I64）、
  `thisCount`（当前输入内频次，char 1..127）、`maxCount`（任意单个输入内最大频次，char）、`loc`
  （位置，I64，正负表方向）。`.1khash` 按 schema 以 `S`（2-bit DNA）、`C`（counts）、`M`（maxCount）、
  `L`（locations）行分块（128Mb/块）写出；`syncmerUpdateMaxCount` 在每轮输入结束后把 `thisCount`
  并入 `maxCount`。
  > **勘误**：`syncmerset.h:34` 注释把 `loc` 写成 "the first location this syncmer was seen"，但
  > `syng.c` 实际用**蓄水池采样**（`rand() % count == 0` 时更新，`syng.c:551-554,561-564`），即等概率抽到的
  > 代表性位置，并非字面 "第一个"。

### 3.4 syng 主程序流程与各命令作用

**主程序 `syng`**（`syng.c`）是一个流式、多线程（pthread，默认 `-T 8`，每线程 100Mb 序列缓冲）管道：
读入 FASTA[.gz]/FASTQ[.gz]/BAM/CRAM/SAM/`.1seq`/`.1path`/`.1gbwt` 或 fofn 文件列表，用
`syncmerIterator`/`syncmerNext` 抽 syncmer（调用时 kmer 参数传 `0`，节点身份由
`kmerHashFindThreadSafe` 查 `pos` 处 `w+k` 长窗口得到），最后按 `outType`（SEQ/PATH/GBWT）写出。
每个输入序列文件结束时报 "yielding N syncs with M extra syncmers"，全程结束报 "average X coverage"
（= 总实例数 / 去重 syncmer 数，即平均深度）。

`syng` 主要选项（`syng.c:291-321` 的 usage）：

| 选项 | 作用 |
| :--- | :--- |
| `-w` / `-k` / `-seed` | syncmer 参数，默认 55 / 8 / 7 |
| `-T` | 线程数，默认 8 |
| `-o <prefix>` | 输出前缀，默认 `syngOut`，作用于其后所有 `write*` |
| `-readK <file>` | 从既有 `.1khash` 继续（用于增量构建） |
| `-zeroK` / `-limitK <min> <max>` | 清零 kmer 计数 / 按计数过滤 kmer 集（`max=0` 无上界） |
| `-histK` | 输出 kmer 计数的二次直方图（`qhist`，含 mean/median/N50） |
| `-noAddK` | 不新增 syncmer，未命中置 0 |
| `-writeK` / `-writeKfa` | 写 `.1khash` / 写 gzip fasta（`.1kmer.fa.gz`） |
| `-writeNewK <prefix>` | 只写新增 syncmer；隐含 `-noAddK` |
| `-writePath` / `-writeGBWT` / `-writeSeq` | 写 `.1path` / `.1gbwt` / `.1seq` |
| `-outputEnds` / `-noEnds` | 写/不写路径两端非 syncmer 的 `X`/`Y` DNA 段（默认写） |
| `-noNames` | 不把序列名写入 path/gbwt（read 集用） |

**其余命令**：
- `syngmap`（`syngmap.c`）：`syngmap <.1khash> <.1gbwt> <query>`，把查询读段以 **MEM**（maximal exact
  match）映射到图，输出 `.1map`。`M` 行记 mem 的 start/end/count，`U` 行记唯一比对
  （file/path/offset，负 offset 表反向），`X` 行记图中缺失的 syncmer（附其序列），`F` 行记被过滤序列。
  过滤器：`-filterG <nG>`（连续 G 超 nG 的坏 Illumina 读）、`-filterQ <QT>`（平均质量低于 QT）、
  `-filterIllumina`（等价 `-filterG 60 -filterQ 20`）、`-outputIds`。MEM 定位用 GBWT 前向 + 回退搜索
  （`syngBWTmatchStart`/`syngBWTmatchNext`），并借助唯一 syncmer（`count==1`）的 `loc` 回溯 `syngBWTlocFind`。
- `syngpath2gbwt`（`syngpath2gbwt.c`）：`XX.1path YY.1gbwt`，把显式路径列表转成隐式 GBWT（每条 path 正反向
  各加入一次），是生成 `.1gbwt` 的两步法第二步。
- `syngstat`（`syngstat.c`）：对 ONEcode 文件统计；对 `gbwt` 报告顶点/边/序列数并调 `syngBWTstat`。
- `k31type`（`k31type.c`）、`ONEview`（ONEcode 自带的 ONE 文件查看器）。

## 4. 与 pgr 当前 minimizer 实现的对比

pgr 现有实现在 [src/libs/hash.rs](file:///home/wangq/Scripts/pgr/src/libs/hash.rs)，关键 API：

| pgr API（现有，保留） | 作用 | 对应 syncmer 新增项 |
| :--- | :--- | :--- |
| `seq_mins(seq, hasher, k, w) -> RapidHashSet<u64>` | 生成 minimizer 哈希集合（用于 Jaccard/Mash） | `seq_syncmer_set`（同签名风格） |
| `seq_sketch(seq, seq_id, k, w, soft_mask, filter) -> Vec<MinimizerInfo>` | 带 positional/strand 的 sketch（用于 mapping） | `syncmer_dna` / `syncmer_protein`（保留 `MinimizerInfo`） |
| `load_minimizers(infile, hasher, k, w, is_merge) -> Vec<MinimizerEntry>` | 从 FASTA 加载 | `load_syncmers`（同模式） |
| `set_distances` / `mash_distance` / `mash_to_sim` | 距离度量 | **不变**（与采样器无关，共用） |

调用点（已接入；`--sampler` 在 `dist hv` 的 `execute` 中分流到 syncmer 或 minimizer 路径；`dist seq` 已于 2026-08-08 删除，其 syncmer 路径随 `pgi build` 保留、minimizer 路径由 `dist mini` 继承）：

- [src/cmd_pgr/dist/mini.rs](file:///home/wangq/Scripts/pgr/src/cmd_pgr/dist/mini.rs) — minimizer 草图 → `pgr dist mini`（原 `dist seq` 的 minimizer 模式）
- [src/cmd_pgr/dist/hv.rs](file:///home/wangq/Scripts/pgr/src/cmd_pgr/dist/hv.rs) — `--sampler syncmer` → `load_hv_from_fasta_syncmer`；否则 `load_hv_from_fasta` → `pgr dist hv`
- [src/cmd_pgr/pgi/build.rs](file:///home/wangq/Scripts/pgr/src/cmd_pgr/pgi/build.rs) — `pgr pgi build` 用 closed syncmer 作为稀疏 k-mer 索引种子：默认 `--smer 8 --window 5`（与 FastGA GIX 参数一致）、`-k 40`，**锚定在窗口最小 s-mer 端点**（syng 式 closed syncmer），而非 GIX 的 window-start match-mer。这是 syncmer 算法在 pgr 中除 `dist hv` 采样外的第二处落地（索引/比对侧），说明同一 `seq_syncmer_set`/`syncmer_dna` 核心可同时服务"草图距离"与"稀疏种子"两类用途
- `set_distances` / `calc_distances` / `mash_distance`（不变，与采样器无关）
- `--kmer`/`--window` 默认值由 `args::resolve_kmer_window` 统一分流（两命令共用）

**算法层面差异**：

- pgr 的 `JumpingMinimizer`（`hash.rs:43`）先对全文所有 k-mer 预算哈希（`hash_kmers`），再做"跳跃式"选最小——O(n) 内存且语义是经典 minimizer。
- pgr 的另一条路径用 `minimizer_iter` crate（`hash.rs:111`、`seq_sketch`），已是滚动窗口式。
- syng 的 syncmer 迭代器是 O(w) 内存的滚动式（环形缓冲区，无需预算全部哈希）。**已优化（2026-08-03 及后续）**：`dna_canonical_hashes` 改为流式迭代器（`syncmer.rs:187`），`closed_syncmers_stream` 采用**分块前缀/后缀最小值**的滑动窗口最小方案（simd-minimizers 思路的标量版，`syncmer.rs:92-97` 注释明确 "no monotonic deque needed"）：用大小为 `window` 的块前缀/后缀最小值求窗口最小，另配一个环形缓冲 `ring` 保留最近 `window` 个原始条目（哈希 + extra）以做"端点最小"判定——DNA/蛋白路径均不再预算全部 s-mer 哈希（内存 O(window)）。`closed_syncmers_from_hashes` 现在是流式核心的薄包装，语义与输出完全一致（strand-symmetry / bounded-gap / density 测试全绿，pgi build 的 single-pass 对照测试亦通过）。

## 5. 对 pgr 的启示与实现计划

> **状态：已实现** — `src/libs/syncmer.rs` 已落地（核心算法 + DNA/蛋白双轨），`pgr dist seq` 与 `pgr dist hv` 均已接入 `--sampler syncmer`。以下为设计稿与实现记录。
>
> **定位**：不全面替换 minimizer，而是新增 syncmer 作为采样器选项，与现有 mod-minimizer/minimizer 长期并存（见 §5.2 采样器矩阵）。目标是为 `pgr dist seq` / `dist hv` 提供对 indel/重排更稳定的备选采样器，而非强制切换。
>
> **与原计划的分歧**：原计划（§5.3.1）拟用 `--sampler mod-minimizer|syncmer|minimizer` 三选一。实际实现中 mod-minimizer 仍保留在 `--hasher mod` 下（非破坏性，沿用既有 CLI），`--sampler` 只承载 `minimizer|syncmer` 两选一。

### 5.1 核心移植要点

1. **算法体量极小**：closed syncmer 迭代器核心约 65 行 C（`seqhash.c:186-250`），Rust 移植后预计 80–120 行，是本次工作的全部硬核内容。**不要过度设计**——不需要堆优化、不需要 ONEcode 持久化、不需要 GBWT。
2. **参数命名陷阱**：移植时务必在 Rust 侧用清晰命名避免 syng 的 `k`/`w` 歧义。建议：
   ```rust
   pub struct SyncmerParams {
       pub smer: usize,   // syng 的 k，哈希粒度（小 k-mer 长度）
       pub window: usize, // syng 的 w，窗口内 s-mer 个数
       pub seed: u64,
   }
   // syncmer 长度 = smer + window
   ```
3. **哈希函数选择**：syng 用乘加移位（`k * factor1 >> shift1`）。最终实现按双轨选择——DNA 路径采用 syng 式 2-bit packed 乘加移位（`k_hash = x * factor >> (64-2k)`，[syncmer.rs:198](file:///home/wangq/Scripts/pgr/src/libs/syncmer.rs#L198)，`factor` 用 splitmix64 由 seed 生成，与 syng 的 `libc random()` 值不同但同样均匀），蛋白路径复用 `RapidHash` 作用于 s-mer 字节串。`--hasher` 在 syncmer 路径被**整体忽略**：DNA 用 2-bit 自带哈希，蛋白也固定用 `RapidHash`（`seq_syncmer_set` 硬编码，`dist hv` 的 syncmer 分支 `hv.rs:137-145` 不传 `--hasher`，`load_hv_from_fasta_syncmer` 亦无 hasher 参数）——`--hasher`（rapid/fx/murmur）只作用于 minimizer/frachash 路径。
4. **canonical 处理**：syng 同时维护 `h` 与 `hRC` 取 min。pgr 的 `seq_sketch` 已通过 `.canonical()` 做了等价事；移植 syncmer 时需在迭代器内部完成（因为判定"端点最小"必须用 canonical 哈希），不能依赖外部 crate 的后处理。

5. **氨基酸适配（硬约束）**：pgr 当前 minimizer 同时服务 DNA 和蛋白（[dist/mini.rs](file:///home/wangq/Scripts/pgr/src/cmd_pgr/dist/mini.rs)，蛋白 `-k 7 -w 2`、DNA `-k 21 -w 5`），靠的是字节串哈希（`rapid`/`fx`/`murmur`）对任意字母表工作。但 syng 的 syncmer 实现 **DNA 强绑定**：2-bit 编码（仅 4 碱基）、`patternRC` 反向互补、canonical 三处都假设 DNA。**蛋白没有反向互补链概念**，因此蛋白 syncmer 反而更简单——去掉 canonical 即可。移植必须双轨：DNA 路径保留 canonical（链无关性对 `pgr dist hv` 的距离稳定性必要），蛋白路径用字节哈希、不做 canonical。

### 5.2 建议的模块结构

遵循 pgr 的分层原则（复杂逻辑放 `libs/`），新建独立模块，**不动 `hash.rs` 里的距离度量部分**。核心是把"窗口端点最小"的几何判定与"如何对 s-mer 哈希"解耦——核心算法只写一次，DNA/蛋白各自包装：

```text
src/libs/
├── hash.rs        # 保留：SetDistances, mash_distance, mash_to_sim, MinimizerEntry, Hasher trait
└── syncmer.rs     # 新增：核心算法 + DNA/蛋白双轨包装
```

API 草案：

```rust
pub struct SyncmerParams { pub smer: usize; pub window: usize; pub seed: u64; }

// 字母表无关核心：给定 s-mer 哈希流，返回最小哈希落在窗口首/末的位置。
// 即 closed syncmer 的几何判定，~30 行，不关心 DNA 还是蛋白。
fn closed_syncmers_from_hashes(hashes: &[u64], w: usize) -> Vec<usize>;

// DNA 路径：2-bit 滚动哈希 + canonical（忠实移植 syng seqhash.c）。
// 返回 strand，保证两条链产生相同 sketch。
pub fn syncmer_dna(seq: &[u8], params: &SyncmerParams) -> Vec<(u64, usize, bool)>;

// 蛋白路径：复用 hash.rs 的 Hasher trait 作用于 s-mer 字节串，无 canonical。
// 行为与当前 minimizer 的字节哈希一致，drop-in 替换。
pub fn syncmer_protein(seq: &[u8], params: &SyncmerParams, h: impl Hasher)
    -> Vec<(u64, usize)>;

// 便捷 dispatch：按 is_protein 在 syncmer_dna/syncmer_protein 间选择。
// cmd_pgr 层按 --sampler 统一分发，三套采样器长期并存：
//   --sampler mod-minimizer → hash.rs 现有 mod 路径（DNA canonical，保留）
//   --sampler syncmer   → 本模块 syncmer_dna（DNA）或 syncmer_protein（蛋白）
//   --sampler minimizer → hash.rs 现有 rapid/fx/murmur 路径（DNA+蛋白，保留）
pub fn seq_syncmer_set(seq: &[u8], params: &SyncmerParams, is_protein: bool)
    -> rapidhash::RapidHashSet<u64>;
```

### 5.3 切换策略

1. **长期并存（不强制替换，已实现）**：`dist seq` 与 `dist hv` 均已增加 `--sampler minimizer|syncmer` 选项（mod-minimizer 仍走 `--hasher mod`，见上方分歧说明）。三套采样器长期共存，默认值可后续按实证表现调整，但不以"替换 minimizer"为目标——mod-minimizer 与 minimizer 作为回退与对照保留。
2. **测试基准（已实现）**：在 `syncmer.rs` 中用随机化属性测试断言"有界间隔"性质——相邻 syncmer 端点位置之差 ≤ `2(w-1)`（w≥2；即任意 `2(w-1)` 个连续 s-mer 位置至少含一个 syncmer 端点）。这是 minimizer 不具备的可验证不变量。
   > **勘误**：原设计稿曾写"任何长度 ≥ `smer+window` 的子串至少命中一个 syncmer"（等价于间隔 ≤ `w+1`），实测**错误**——随机序列下间隔可达 `w+2` 及以上。**紧致上界为 `2(w-1)`**（穷举 w∈{3,4,5} 的全排列与 4-字母表序列、随机 w≤32 均验证，且可达）：取非端点连续段内的最小哈希位置 `p*`，因 `p*` 非端点，其左右两个 w-窗口内必各有更小值；而 `p*` 已是段内最小，更小值只能在段外，故段两端 `a<b` 满足 `b-w+1 ≤ p* ≤ a+w-1`，即 `b-a ≤ 2(w-1)`。原稿"两侧 w 范围"应为 `w-1`（窗口内除端点自身外的 `w-1` 个位置），故上界是 `2(w-1)` 而非 `2w`。注意序列首尾**不保证**被覆盖（§2.2），仅保证内部间隔有界。
3. **参数等价性**：切换到 `--sampler syncmer` 时 `--kmer`/`--window` 语义会变（见 §2.1 表），需在文档与 CLI 帮助中显式说明，避免用户沿用旧参数得到不同密度。

### 5.4 预期收益与风险

- **收益**：序列距离对 indel/重排更稳定（syncmer 集合的局部扰动性）；采样密度有理论保证，便于跨样本对齐比较；保留 mod-minimizer/minimizer 作为回退，迁移风险可控。
- **风险**：syncmer 平均深度 ~2×，sketch 集合比 minimizer 大约一倍，Jaccard/Mash 距离的尺度会变化——需重新校准 `pgr dist seq` 在 `--sampler syncmer` 下的阈值与 `--kmer` 默认值。**默认值已校准**（`resolve_kmer_window`：DNA smer=8/window=55 即 syng 默认；蛋白 smer=7/window=5，k=7 使随机碰撞概率可忽略且与 minimizer 蛋白惯例一致，w=5 给 ~33% 密度适配短序列）。当前无距离阈值参数，距离尺度变化待实证评估。
- **硬约束**：必须同时支持 DNA 和蛋白（见 §5.1 第 5 点）。蛋白走字节哈希无 canonical 的独立 syncmer 路径，不能为了忠实移植 syng 的 2-bit+canonical 而漏掉蛋白。验证标准：蛋白 syncmer 集合与同序列的 minimizer 集合在密度量级上一致（不要求完全相同，因采样器语义不同），且 `pgr dist seq` 蛋白用法（`-k 7 -w 2`）在 `--sampler syncmer` 下可用。

---

*参考来源: [syng GitHub](https://github.com/richarddurbin/syng) | [Edgar, 2021 — Syncmers are more sensitive than minimizers](https://peerj.com/articles/10805/) | [ONEcode](https://github.com/thegenemyers/ONEcode)*
