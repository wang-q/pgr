# Ropebwt3 源码分析

> 整理于 2026-08，源自对 `ropebwt3-3.10/` 目录源码（Heng Li，V3.10-r281，2025-11）及
> README/NEWS 的通读。目的：理解基于 FM-index 的**无 pairwise 比对**泛基因组索引路线
> （RLE-BWT + SMEM + BWA-SW），与 pgr 的"pairwise 比对 + PAF 隐式图"路线对照，
> 评估可借鉴的算法。

## 1. Ropebwt3 概览

- **工具定位**：为高度冗余的序列集（泛基因组、高覆盖 reads）构建 FM-index，并查询
  SMEM（supermaximal exact matches）或做局部比对（修订版 BWA-SW）。
- **作者/版本**：Heng Li（lh3）；首版 3.0（2024-06），当前 3.10-r281（2025-11）。
- **关键数字**：
  - 7.3 Tb 常见细菌基因组 → **30 GB RLE-BWT**（无损压缩）；
  - 152 M. tuberculosis 基因组 / 472 个人类长读组装（18 GB 下载）的预建索引已发布。
- **路线本质**：不依赖任何 pairwise 比对——直接对**序列集合**建 FM-index，
  查询（SMEM/局部比对）直接在索引上进行。与 pgr 的"外部比对器产出 PAF →
  PAF 隐式图"是完全不同的架构。

## 2. 整体架构与数据流

```
FASTA/FASTQ 序列集
  │  build（两种算法：ropebwt2 动态 / libsais 广义 SA）
  ▼
BWT（FMR 动态 rope 或 FMD 静态 RLE）
  │  FMD + ssa（采样后缀数组）+ len.gz（序列名/长度）
  ▼
FM-index（rld rank + 采样 SA）
  ├─ mem   → SMEM（Gagie 算法 / 原始算法），BED/覆盖统计
  ├─ sw    → BWA-SW 局部比对 / end-to-end，PAF（含 rh/cs 等 tag）
  ├─ hapdiv→ 101-mer 单倍型多样性
  └─ get   → 按索引取回序列
```

- 完整索引 = 三个文件：`.fmd`（RLE-BWT，支持 rank）、`.fmd.ssa`（采样后缀数组，
  供 sw 报告坐标）、`.fmd.len.gz`（序列名/长度，供 PAF 输出）。
- **双链对称性假设**：默认输入第 i 条序列作为 BWT 的第 2i 条、其反向互补作为第
  2i+1 条；`mem`/`sw`/`hapdiv` 在搜索前会用 `rb3_fmi_is_symmetric` 校验 BWT 双侧链
  计数对称，否则报错退出。查询端由 `sid` 的奇偶区分正反链（`pos_stranded`）。
- 两种 BWT 格式：
  - **FMR**（ropebwt2/fermi 动态格式）：可增量追加/合并（`-i` 续建），BWT 构建
    用；同一 BWT 不保证同一 FMR。文件头 magic = `"RB\2"`（mrope.c `mr_dump`，
    `mrope.c:152-159`），随后 6 条 rope（ACGTN$ 各一）逐条 dump。
  - **FMD**（静态格式）：结构简单、加载快、内存小、可 mmap，查询用。文件头
    magic = `"RLD\3"`（rld0.c `rld_dump`，`rld0.c:222-243`）；`rld_restore` 按
    内存加载，`rld_restore_mmap` 按 `mmap` 零拷贝映射（`rld0.c:322-341`）。
  - 可互转（`build -i in.fmd -bo out.fmr` / 反向），转换走 `rb3_enc_fmr2fmd` /
    `rb3_enc_fmd2fmr`（fm-index.c，前者遍历 FMR 重编码 RLE，后者按 `cnt[]` 累积
    计数切分每个 rope 段）。
  - 其余二进制格式：SSA 采样后缀数组 magic = `"SSA\1"`（ssa.c:204）、BRE 文件
    magic = `"BRE\1"`（rld0.c:279）——解析时都用 magic 判断格式，跨格式复用同一
    `rb3_fmi_restore`（fm-index.h:123-133），先试 FMD 失败再试 FMR。

## 3. 核心算法

### 3.1 BWT 构建（build.c，两种算法）

1. **ropebwt2 动态算法**（`-2`/`-s`/`-r` 选项）：
   - 把序列**反转**后逐条插入 multi-rope（`mr_insert_multi`）；
   - rope（rope.c）是 B+-tree 风格的块结构（`ROPE_DEF_MAX_NODES=64`、
     `ROPE_DEF_BLOCK_LEN=512`，可用 `-n`/`-l` 覆盖；`max_nodes` 强制偶数、
     `block_len` 强制 8 对齐），内部节点 `rpnode_t` 用位域压缩（`l:54, n:9,
     is_bottom:1`），每个节点保存 6 个 marginal counts（ACGTN$）；
   - **叶子块**即一段 RLE 编码的 BWT 子串，采用 rle.h 的 "43+3" codec（详见 3.2）；
   - 插入 = 在 BWT 中做 backward 插入，沿树**自上而下"搜索 + 分裂"单趟**完成：
     节点满（`n==max_nodes`）或叶子 `n_runs + RLE_MIN_SPACE(18) > block_len`
     时 `split_node` 分裂，动态维护 BWT；`rpcache_t` 缓存上次插入位置加速连续插入。
   - `-r` 用 RCLO（reverse-complement + 计数排序优化），适合短读。
   - **免释放内存池**（rope.c `mempool_t`，`rope.c:13-49`）：节点与叶子各一个
     bump 分配器，按 `MP_CHUNK_SIZE=0x100000`（1MB）分块、只分配从不 free——
     避免索引构建期百万级 `rpnode_t` 的小对象 malloc/free 碎片与开销；树深上界
     `ROPE_MAX_DEPTH=80`（rope.h:7）。这是 pgr 大规模索引构建时可借鉴的
     "固定大小池 + 无释放"模式。
2. **libsais 批量算法**（默认）：
   - 按 batch 读入序列（含正反链），用 **libsais 构建广义后缀数组**（libsais.c，
     外部库；`sais-ss.c` 的 `rb3_build_sais` 按 batch 长度选择：`len + sais_extra_len >=
     INT32_MAX` 时用 64 位 `libsais64_gsa`，否则用 32 位 `libsais_gsa` 省内存；
     "libsais16x64" 并不存在）→ 由 SA 生成部分 BWT；
   - batch 大小默认 **7G** 符号（`-m`，`rb3_parse_num` 支持 K/M/G 后缀）；
   - **SA→BWT 直接构造**（`rb3_build_sais32/64` 的亮点）：不必对 SA 再排序，而是对每个
     后缀取其**前驱字符**填回原位——`SA[i] = T[SA[i]==0 ? len-1 : SA[i]-1]`（位置 0 的
     后缀取最后一个字符做循环哨兵），一趟把广义 SA 就地转成 BWT，省一次重排；
   - 多批次时把部分 BWT **增量合并**到已有 rope（`rb3_fmi_merge_plain`）；
   - `-p` 可让部分 SA 构建与 BWT 合并并行（3.9+，`worker_pipeline` 两步流水线，SA 用
     `sais_threads`、合并用 `n_threads - sais_threads` 线程）；注意仅在**单输入文件**时
     启用（build.c 判断 `argc - o.ind == 1`），多文件时回到逐文件串行流程。

输出格式：`-b` FMR（rope dump）、`-d` FMD（rld 编码）、`-T` TREE、`-e` BRE
（block-run encoding）、默认 PLAIN 文本 BWT。

> **外部构建路径（grlBWT，README:140-145）**：当输入为单文件时，可用
> `fa2line` 把序列（正反链）转成每行一条 → 交给外部 **grlBWT** 构建（可能更快
> 但需工作磁盘空间）→ `grl2plain` 出纯文本 BWT → `plain2fmd` 转回 FMD。
> 这体现 ropebwt3 的"构建后端可插拔"设计：`plain2fmd`（main.c:299-331）只是把
> 纯文本 BWT 逐符号 `rld_enc` 成 FMD，任何外部 SA/BWT 工具的产物都能接入。

### 3.2 RLE-BWT 存储（rld0.c + rle.h）

FMD 的核心数据结构 `rld_t`（rld0.h），其 RLE 编码分两层：

- **rope 叶子层的 "43+3" codec（rle.h）**：BWT 连续相同字符压缩为 (symbol, run)，
  每个 run 低 3 位存符号、高位存长度，长度按阈值分档**变长**：`<2^4` 占 1 字节、
  `<2^8` 占 2、`<2^19` 占 4、否则占 8 字节（`rle_enc1`/`rle_dec1`）——这是
  FMR/rope 在内存与磁盘上的 RLE 格式；
- **FMD 的 delta 编码（rld0.c）**：把 run 流再压缩——run 长度用 **delta 编码**
  （`rld_delta_enc1`，宽度 = `2⌈log₂⌈log₂l⌉⌉ + 1 + ⌈log₂l⌉` 的自适应变长码），连同
  symbol（`abits=3` 位，`asize=6`）打包成连续 bit 流写入 64-bit 字；
  - 数据按 **`2^23` 字（`RLD_LSIZE`）分块**（`z[]`），每块内再细分为
    `ssize=2^sbits` 字（默认 `bbits=3` → 8 字）的 small block，每个 small block
    头部存**该块内的累计符号计数**（按块内总数自动选 16/32/64-bit 三种宽度，
    即 `rld_block_type` 0/1/2 与 `offset0[]`），便于随机 seek 与 skip；
  - **rank 支持**：`rld_rank_index` 建 `frame` 索引（每 `2^ibits` 个符号采样一次
    块级前缀计数），配合 `cnt`（累计）/`mcnt`（边际）计数，`rld_rank*`/`rld_extend`
    近似 O(1)——这是 FM-index 查询（occ/LF）的基础；
- alphabet：DNA 模式 `asize=6`（ACGTN$，`RB3_ASIZE`），`_DNA_ONLY` 时启用
  `rld_dec0_fast_dna` 快速解码热路径。

### 3.3 FM-index 查询（fm-index.c）

- `rb3_fmi_t` = rld（FMD）+ mrope（FMR）+ `rb3_ssa_t`（采样 SA）+ `rb3_sid_t`
  （序列名索引）+ `acc`（累计计数表）。
- **backward search / extend**：`rb3_fmd_extend(f, ik, ok, is_back)` 用 rank
  扩展区间（LF 映射）；`rb3_fmd_set_intv` 初始化单字符区间。
- **采样后缀数组**（`rb3_ssa_t`）：每 `1<<ss` 个位置存一个 SA 值，低位 `ms` bits
  存序列 ID、高位存序列内偏移；`r2i` 在 backward search 到哨兵时给出序列 ID。
- 3.8 优化：m 个高度相似基因组时，定位一个区间内位置集合的期望时间 O(s/m)
  （s = 采样率）。

### 3.4 SMEM（search.c）

- **SMEM 定义**：在查询序列上不被任何更长 MEM 包含的 maximal exact match。
- **原始算法** `rb3_fmd_smem1`（`--old-mem` 切回）：对每个起点 x 先**前向**扩展区间
  （`rb3_fmd_extend`），记录所有 size ≥ `min_occ` 的区间；再**后向**扩展，区间无法延伸时
  输出一个 MEM（长度 ≥ `min_len`）。每位置尝试一次，汇总得全部 SMEM。
- **Gagie 算法** `rb3_fmd_smem1_TG`（3.2 起默认，更快；fm-index.c:483-518）：
  对每个起点 `x` 从 `x + min_len - 1` 处的单字符区间起步，**先反向**扩展到
  `min_len`（中途 `size < min_occ` 即判定无 MEM、返回下一起点）；再**前向**扩展
  直到无法延伸，输出这一个长 MEM（`info = x<<32 | j`）；随后**再反向**扩展一次
  求下一个起点。全程只输出长 MEM，省去原始算法对每个位置尝试、收集大量短 MEM
  的重复扫描。`check_long` 参数切换"只探测存在性"模式，供 `-j` 的
  `rb3_fmd_smem_present` 门控预检（fm-index.c:530-538）复用同一热路径。
- `mem` 输出 BED（query name, start, end, #hits），默认不输出位置（`-p` 可选，
  3.8 起输出半随机子集）；`--gap`/`--cov` 报告未覆盖/覆盖长度。
- **歧义碱基**：query 非 ACGT 字符经 `rb3_nt6_table` 统一编码为 N（code 5，合法 BWT
  符号），避免 mem 段错误（3.4 修复）；sw 的 query DAWG 构建中再把 N 转成 A。

### 3.5 BWA-SW 局部比对（bwa-sw.c，修订版）

- **DAWG**（dawg.c）：对查询序列构建有向无环词图，BWA-SW 在 query 的 DAWG 上跑 DP，
  而不是 query 后缀数组。局部模式用 `rb3_bwtl_gen` 先建 query 的 BWT 再 `rb3_dawg_gen`
  构造完整 DAWG；end-to-end 模式则用**线性 DAWG**（`rb3_dawg_gen_linear`，query 本身即
  一条路径，适合 `hapdiv`/`-e` 的整条查询）；
- **query 轻量 BWT**（`rb3_bwtl_gen`，dawg.c:28-76）：对短查询用 libsais 建 SA，
  BWT 用 **2-bit 打包**（每 4 基一字节、每 16 基 32-bit 字），rank 每 16 基存一次
  前缀计数（`occ[]`）；`rb3_bwtl_rank1a`（dawg.c:78-89）用预计算的
  `bwtl_cnt_table[256]` **查表统计单字内 2-bit 符号数**（一次处理 4 字节），
  是 SW 里"小字符串上反复 rank"的紧凑 SIMD 式优化——pgr 若做 query 端小索引，
  这种 2-bit 打包 + 查表 rank 是省内存且快的模板；
- **候选集**：`sw_cell_t` 携带 H/E/F 三个 DP 状态与 **SA 双区间**（`lo, lo_rc, hi-lo`）；
  `sw_candset_t` 以 `(lo,hi)` 为键去重/合并同一区间的多个来源，`sw_update_candset`
  裁剪，`sw_heap_insert1` 用堆保留 top `n_best` 候选；
- **rank 缓存**：`rb3_r2cache_init` 用哈希表缓存 `(k,l)` 的 rank 结果（默认容量
  `0x10000`，`-C` 可调），DP 中 `rb3_fmd_extend_cached` 复用，省去大量重复 occ 计算；
- **矩阵布局**：`cell` 按 `n_node × n_best` 铺开（每 DAWG 节点至多保留 `n_best` 个
  cell，`row[i].a` 指向各节点行）——`-N`（`n_best`，默认 25）即 DP 的"带宽"；
- **F 状态**：`sw_track_F` 用 `fpar`（`rb3_u128_t` 数组，存前驱区间）跟踪 F 数组
  （前向 gap 延伸）的 `F_from_off`，供回溯；`sw_cell_dedup` 删除 out-of-band cells
  （3.4 提速 20%）；
- **回溯**：`sw_backtrack1` / `sw_backtrack` 从最优位置回溯路径，输出 CIGAR
  （`cs` tag，3.7+）；
- **评分**：默认 BLASTN 风格参数（`rb3_swopt_init`：`match=1, mis=3, gap_open=5,
  gap_ext=2, min_sc=30`，`-A/-B/-O/-E/-m` 可覆盖）；`-N` 设每个 DAWG 节点保留的候选
  hit 数（`n_best`，默认 25），`-k` 要求比对末端有 k-mer 精确匹配（`end_len`，默认 11），
  `-j` 设启动比对所需的最小 MEM 长度（`min_mem_len`，默认 0；>0 时先跑
  `rb3_fmd_smem_present` 做**长 MEM 门控预检**，query 不含足够长 MEM 就直接跳过 SW，
  是省算力的廉价前置过滤）；`-e` end-to-end 模式
  输出相似单倍型（`--all-e2e` 紧凑输出）。
- **性能提示**：局部比对比 SMEM 慢数十倍，不用于高通量 reads。

> **`sw` 的 PAF 输出 tag**（`search.c` `write_paf`）：12 列之外固定追加 `AS:i:`（得分）、
> `qh:i:`（query 端命中数 n_qoff）、`rh:i:`（ref 端区间长度 hi-lo）、`cg:Z:`（CIGAR）、
> `cs:Z:`（cs 串）；`--seq` 时再追加 `rs:Z:`（参考序列）；多位置命中（n_pos>1）时追加
> `ap:Z:`/`aq:Z:`（有/无序列名时，逐位置 `name,strand,pos;`）。其中 `cg`（CIGAR）与 `cs`
> （cs 字符串）是**两个不同的 tag**，分别对应 `h->cigar` 与 `h->cs`。

### 3.6 hapdiv（单倍型多样性）

把 end-to-end 模式应用于滑动的 101-mer，报告：命中的不同等位基因数、最大编辑距离、
完美匹配单倍型数、按编辑距离分桶的计数（`RB2_SW_MAX_ED=6`，桶 ed=0..6，≥6 的编辑
距离并入末桶；3.6 版为距离 5）。

## 4. 源码模块结构

| 模块 | 职责 |
|------|------|
| `build.c` | BWT 构建主流程（ropebwt2 动态 / libsais 批量 + merge，FMR/FMD/TREE/BRE 输出）|
| `mrope.c` / `mrope.h` | multi-rope：把 BWT 按**后缀首字符**分成 6 段（ACGTN$），每段一条 rope（`r[a]`）；`mr_insert1`/`mr_insert_multi` 序列插入（正反链）、`mr_rank2a` 跨 rope 定位计数 |
| `rope.c` / `rope.h` | B+-tree 风格动态 BWT 块结构（marginal counts、插入、rank）|
| `rld0.c` / `rld0.h` | RLE-BWT：run-length + delta 编码、rank、FMD 读写（`rld_t`）|
| `fm-index.c` / `fm-index.h` | FM-index：extend/backward search、SMEM、SSA、格式互转 |
| `search.c` | `mem`/`sw`/`hapdiv` 命令：SMEM（Gagie）、PAF 输出、gap/cov 统计 |
| `bwa-sw.c` | 修订版 BWA-SW：DAWG DP、候选集、F 状态、回溯 |
| `dawg.c` / `dawg.h` | query 的 DAWG 构建（BWA-SW 用）|
| `ssa.c` | 采样后缀数组构建（`ssa` 命令）|
| `bre.c` / `rle.c` | BRE（block-run encoding）与 RLE 工具 |
| `io.c` / `kseq.h` / `kalloc.c` / `kthread.c` | 序列 I/O、khash/ksort/kthread 等 lh3 公共库 |
| `sais-ss.c` | libsais 封装：SA→BWT 就地转换 + 32/64 位选择（`rb3_build_sais*`）|
| `libsais.c` / `libsais64.c` | 外部广义后缀数组库（IlyaGrebnov）|
| `main.c` | 命令分发 |
| `rb3tools.js` | mappability 过滤 / 简单 SNP 调用辅助脚本 |

## 5. 命令与关键参数

| 命令 | 功能 | 关键参数 |
|------|------|----------|
| `build` | 构建 BWT（FMR/FMD）| `-2/-s/-r` ropebwt2 算法（`-s` RLO、`-r` RCLO）、`-m` batch（默认 7G）、`-t` 线程、`-p` SA 并行线程、`-l/-n` 叶子块长/节点扇出、`-F/-R` 免正/反链、`-d/-b/-e/-T` 输出格式、`-i` 续建、`-S` 逐文件保存 |
| `mem` | SMEM 查找 | `-l` 最小长度（默认 19）、`-c` 最小出现、`-p` 输出位置、`--gap/--cov`、`--old-mem` 切回原始算法（默认 Gagie）|
| `sw` | BWA-SW 局部比对 | `-N` 每 DAWG 节点候选数（`n_best`）、`-k` 末端 k-mer（`end_len` 默认 11）、`-j` 启动 MEM 长度（`min_mem_len`）、`-m` min score、`-A/-B/-O/-E` 评分、`-C` rank 缓存、`-y` e2e 丢尾、`-e` 端到端（强制 `-k 1`）、`-b` 双链、`-u` 输出未比对、`-p` 多位置、`--seq` 输出 rs、`--all-e2e`/`-g`（也强制 `-e` + `-k 1`）|
| `suffix` | 找最长匹配后缀 | `-L` 输入单行一条序列 |
| `hapdiv` | 101-mer 单倍型多样性 | 内部调用 sw -e（`hapdiv_k=101, hapdiv_w=50`）|
| `ssa` | 采样后缀数组 | `-s` 采样率（每 2^INT 碱基一个 SA，默认 8）、`-t` 线程；输出文件 ≈ `64·(n/2^s + m)` 字节（n=符号数、m=序列数）|
| `get` | 按索引取序列 | `get <idx.fmr> <int> [...]` |
| `merge` | 合并多个 BWT（FMR）| `-t` 线程、`-o` 输出、`-S` 中途保存 |
| `kount` | 统计高出现 k-mer | `-k` k-mer 长、`-m` 最小出现（`kount -k 51 -m 100` 类）|
| `fa2line` | FASTX 转行（正反链）| `-R` 不含反向链 |
| `fa2kmer` | FASTX 抽 k-mer | `-k` 长度（默认 151）、`-w` 步长（默认 50；对末尾 k-mer 做截断处理）|
| `plain2fmd` | 纯文本 BWT → FMD | `-o` 输出 |
| `stat` | 报告序列数/符号数/run 数与 A/C/G/T/N 计数 | `-M` mmap（3.6 起 FMR 也支持）|

> **注意**：`mem`/`sw`/`hapdiv` 实际共享 `main_search` 入口，靠 `argv[0]` 子命令分发
> （`main.c`），`search` 是其通用别名。`main()` 共分发 15 个入口（含 `search` 别名与
> `version`）；`usage()` 帮助文本列出 14 个常规命令。

> **kount 与 pgr 的关联**：`kount` 在 FM-index 上做深度优先遍历统计（出现次数 ≥ `min_occ`）
> 的 k-mer，与 pgr 的 `kmer` 命令功能对应；若未来 pgr 引入 FM-index，kount 的"rank2a 定长
> 区间 DFS + 阈值剪枝"是可直接移植的计数骨架（`main_kount`，main.c:346-423，约 80 行 C）。
> 它接受**多个索引文件**，同步维护各索引的区间栈，输出"该 k-mer 在每个索引中的出现数"
> （`main_kount` 里对每个 `aux[i]` 独立 `rank2a` 并取交集判断），适合做多个样本间的
> k-mer 共享/差异统计。

## 6. 对 pgr 的启示

### 6.1 与 pgr 泛基因组路线的对照

| 维度 | pgr（PAF 隐式图）| ropebwt3（FM-index）|
|------|------------------|--------------------|
| 输入 | 外部比对器（FastGA/lastz）产出的 pairwise PAF | 原始序列集合 |
| 索引 | PAF 区间树 + CIGAR（`pgr paf index`）| RLE-BWT + 采样 SA |
| 查询 | 区间投影 / 传递 BFS | SMEM / BWA-SW 局部比对 |
| 查询输出 | PAF/MAF/VCF/GFA 图 | BED/PAF（单 hit/单倍型）|
| 规模 | 4 万大肠杆菌（pairwise 稀疏 + 传递闭包）| 7.3 Tb 细菌 / 472 人类组装（全量索引）|
| 优点 | 复用现有比对资产、图语义丰富 | 无 pairwise、全量可查、极致压缩 |
| 代价 | 依赖比对器、稀疏图需 BFS 推断 | 索引构建复杂、查询类型受限（无图）|

### 6.2 可借鉴的算法

1. **RLE-BWT + rank**：若 pgr 未来需要"序列集直接查询"（无 pairwise），rld0 的
   run-length + delta 编码 + 块级 frame 索引（近似 O(1) rank）是可移植的存储层
   设计；7.3 Tb → 30 GB 的压缩比是 PAF 索引无法比拟的。
2. **SMEM 的 Gagie 算法**：前向/后向交替扩展只输出长 MEM——pgr 若做原生种子
   检测（如简化 FastGA），可借鉴"只找长种子"的思路控制输出规模。
3. **BWA-SW 的 DAWG + 候选集 + F 状态**：FM-index 上的局部比对模板；pgr 的
   `ScalarAlignmentEngine` 是 O(nm) 内存全矩阵，BWA-SW 的"候选集 + 堆 + `fpar`(u128) F
   回溯"是低内存替代（但依赖 FM-index）。其 `-j` 的**长 MEM 门控预检**
   （`rb3_fmd_smem_present` 先判 query 有无足够长种子，再决定是否跑 SW）是通用思路，
   pgr 在昂贵的图遍历/比对前也可先做廉价种子门控。
4. **采样 SA 的序列 ID 编码**：SA 值低位存序列 ID、高位存偏移 + `r2i` 哨兵映射，
   是"稀疏坐标 + 序列归属"的紧凑方案，pgr 的 PAF 索引坐标定位可参考。
5. **两阶段构建流水线**：build 的"libsais 广义 SA → 就地 SA→BWT → 增量 merge 到 rope"
   配合 `-p` 的 SA/merge 双线程流水线，是"建索引与合并解耦、并行"的工程范式；
   `rb3_build_sais` 按 `len + sais_extra_len` 是否超过 `INT32_MAX` 在 32/64 位
   libsais 间切换，是省内存的取舍，pgr 大规模索引构建可参考。
6. **`ssa` 定位的 O(s/m) 算法**：`rb3_ssa_multi` 用堆维护待展开区间（优先小区间）、
   命中哨兵即得序列归属，把"区间内位置集合"定位加速到期望 O(s/m)（s 采样率、
   m 相似基因组数），对 pgr 在高度冗余集合上的坐标定位有借鉴价值。
7. **索引构建期工程细节**：rope 的"固定大小池 + 无释放"bump 分配（`mempool_t`）、
   每 2-bit 打包 + 查表 rank 的 query BWT（`rb3_bwtl_rank1a`）、以及"magic 字节 +
   统一 `rb3_fmi_restore` 自动识别 FMR/FMD"的多格式加载，都是 pgr 索引模块
   （`fa index`/PAF 索引）可复用的低成本模式；grlBWT 外挂 + `plain2fmd` 接入的
   "构建后端可插拔"设计则提示：pgr 的索引层也应把"核心数据结构"与"具体构建算法"
   解耦，便于将来替换或组合不同构建器。

### 6.3 结论

ropebwt3 代表与 pgr **互补而非替代**的路线：它不需要 pairwise 比对，但也不产出
比对图/变异语义。对 pgr 的 4 万大肠杆菌场景，两线可以共存——PAF 隐式图提供
图查询与 VCF/GFA，ropebwt3 式 FM-index 提供"全量序列集上的快速存在性/局部比对
查询"（如 `sw` 找相似单倍型）。但完整移植 FM-index 构建是大工程（rope + rld +
libsais 集成），当前优先级低于 pgr 已有路线的收尾；值得先落地的只是 RLE 存储与
Gagie SMEM 两个算法模式（各 ~200-300 行）。

## 7. 版本与许可

- 当前 3.10-r281（2025-11-25）：sw -p 多位置输出、libsais 更新。
- 关键历史：3.2 Gagie SMEM + BWA-SW + SSA；3.5 end-to-end + hapdiv + BLASTN 评分；
  3.8 位置定位 O(s/m) + rb3tools.js；3.9 部分 SA 与 merge 并行。
- LICENSE：ropebwt3 本体 MIT（LICENSE.txt，Dana-Farber Cancer Institute）；libsais 为
  **Apache License 2.0**（Ilya Grebnov，见 LICENSE.txt 末尾说明）——移植 libsais 集成时需
  留意许可差异。
- 参考：https://github.com/lh3/ropebwt3 ；相关：fermi/ropebwt2（FMR 格式源）、
  minimap2（BWA-SW 的现代替代）。
