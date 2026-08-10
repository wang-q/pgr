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
- 两种 BWT 格式：
  - **FMR**（ropebwt2/fermi 动态格式）：可增量追加/合并（`-i` 续建），BWT 构建
    用；同一 BWT 不保证同一 FMR。
  - **FMD**（静态格式）：结构简单、加载快、内存小、可 mmap，查询用。
  - 可互转（`build -i in.fmd -bo out.fmr` / 反向）。

## 3. 核心算法

### 3.1 BWT 构建（build.c，两种算法）

1. **ropebwt2 动态算法**（`-2`/`-s`/`-r` 选项）：
   - 把序列**反转**后逐条插入 multi-rope（`mr_insert_multi`）；
   - rope（rope.c）是 B+-tree 风格的块结构（`max_nodes=64`、`block_len=512`），
     每个节点保存 6 个 marginal counts（ACGTN$）；插入 = 在 BWT 中做 backward 插入，
     动态维护 BWT。
   - `-r` 用 RCLO（reverse-complement + 计数排序优化），适合短读。
2. **libsais 批量算法**（默认）：
   - 按 batch 读入序列（含正反链），用 **libsais 构建广义后缀数组**（libsais.c，
     外部库；`sais-ss.c` 的 `rb3_build_sais` 按 batch 长度选择：`len + sais_extra_len >=
     INT32_MAX` 时用 64 位 `libsais64_gsa`，否则用 32 位 `libsais_gsa` 省内存；
     "libsais16x64" 并不存在）→ 由 SA 生成部分 BWT；
   - **SA→BWT 直接构造**（`rb3_build_sais32/64` 的亮点）：不必对 SA 再排序，而是对每个
     后缀取其**前驱字符**填回原位——`SA[i] = T[SA[i]==0 ? len-1 : SA[i]-1]`（位置 0 的
     后缀取最后一个字符做循环哨兵），一趟把广义 SA 就地转成 BWT，省一次重排；
   - 多批次时把部分 BWT **增量合并**到已有 rope（`rb3_fmi_merge_plain`）；
   - `-p` 可让部分 SA 构建与 BWT 合并并行（3.9+，`worker_pipeline` 两步流水线，SA 用
     `sais_threads`、合并用 `n_threads - sais_threads` 线程）。

输出格式：`-b` FMR（rope dump）、`-d` FMD（rld 编码）、`-T` TREE、`-e` BRE
（block-run encoding）、默认 PLAIN 文本 BWT。

### 3.2 RLE-BWT 存储（rld0.c）

FMD 的核心数据结构 `rld_t`（rld0.h）：

- **run-length 编码**：BWT 的连续相同字符压缩为 (symbol, run-length)；
- **delta/变长编码**：run 长度按位块存储（23-bit 块粒度 + 16/32/64-bit 三种块宽
  `offset0`），`rld_dec` 迭代解压；
- **rank 支持**：`frame` 索引（块级前缀计数）+ `cnt/mcnt`（累计/边际计数），
  使 `rld_rank` 近似 O(1)——这是 FM-index 查询（occ）的基础；
- alphabet：DNA 模式 `asize=6`（ACGTN$）。

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
- **原始算法** `rb3_fmd_smem1`：对每个查询位置向后扩展区间（`rb3_fmd_extend`），
  区间为空时输出当前 MEM；每个位置尝试一次，合并得 SMEM。
- **Gagie 算法** `rb3_fmd_smem1_TG`（3.2 起默认，更快）：从 `x + min_len - 1`
  位置起向后扩展到 `min_len`，再**前向/后向交替扩展**，只输出长 MEM，省去大量
  短 MEM 的重复扫描。
- `mem` 输出 BED（query name, start, end, #hits），默认不输出位置（`-p` 可选，
  3.8 起输出半随机子集）；`--gap`/`--cov` 报告未覆盖/覆盖长度。

### 3.5 BWA-SW 局部比对（bwa-sw.c，修订版）

- **DAWG**（dawg.c）：对查询序列构建有向无环词图（`rb3_dawg_gen`），BWA-SW
  在 query 的 DAWG 上跑 DP，而不是 query 后缀数组；
- **候选集**：`sw_candset_t` 维护当前可扩展的 FM-index 区间集合，`sw_update_candset`
  去重/裁剪，`sw_heap_insert1` 用堆保留 top 候选；
- **F 状态**：`sw_track_F` 用 `fpar`（u128 位数组）跟踪 F 数组（前向延伸），
  `sw_cell_dedup` 删除 out-of-band cells（3.4 提速 20%）；
- **回溯**：`sw_backtrack1` / `sw_backtrack` 从最优位置回溯路径，输出 CIGAR
  （`cs` tag，3.7+）；
- **评分**：默认 BLASTN 风格参数（`rb3_swopt_init`：`match=1, mis=3, gap_open=5,
  gap_ext=2, min_sc=30`，`-A/-B/-G/-E/-m` 可覆盖）；`-N` 设每个 DAWG 节点保留的候选
  hit 数（`n_best`，默认 25），`-k` 要求比对末端有 k-mer 精确匹配（`end_len`，默认 11），
  `-j` 设启动比对所需的最小 MEM 长度（`min_mem_len`，默认 0）；`-e` end-to-end 模式
  输出相似单倍型（`--all-e2e` 紧凑输出）。
- **性能提示**：局部比对比 SMEM 慢数十倍，不用于高通量 reads。

### 3.6 hapdiv（单倍型多样性）

把 end-to-end 模式应用于滑动的 101-mer，报告：命中的不同等位基因数、最大编辑距离、
完美匹配单倍型数、按编辑距离分桶的计数（`RB2_SW_MAX_ED=6`，桶 ed=0..6，≥6 的编辑
距离并入末桶；3.6 版为距离 5）。

## 4. 源码模块结构

| 模块 | 职责 |
|------|------|
| `build.c` | BWT 构建主流程（ropebwt2 动态 / libsais 批量 + merge，FMR/FMD/TREE/BRE 输出）|
| `mrope.c` / `mrope.h` | multi-rope：6 条 rope（ACGTN$），序列插入（正反链）|
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
| `build` | 构建 BWT（FMR/FMD）| `-2/-s/-r` ropebwt2 算法、`-m` batch、`-p` SA 并行、`-d/-b` 格式、`-i` 续建 |
| `mem` | SMEM 查找 | `-l` 最小长度、`-c` 最小出现、`-p` 输出位置、`--gap/--cov`、`--old-mem` 切回原始算法（默认 Gagie）|
| `sw` | BWA-SW 局部比对 | `-N` 每 DAWG 节点候选数、`-k` 末端 k-mer（`end_len`）、`-j` 启动 MEM 长度、`-m` min score、`-e` 端到端、`-p` 多位置、`--all-e2e` |
| `suffix` | 找最长匹配后缀 | `-L` 输入单行一条序列 |
| `hapdiv` | 101-mer 单倍型多样性 | 内部调用 sw -e（`hapdiv_k=101, hapdiv_w=50`）|
| `ssa` | 采样后缀数组 | `-s` 采样率（每 2^INT 碱基一个 SA，默认 8）、`-t` 线程 |
| `get` | 按索引取序列 | `get <idx.fmr> <int> [...]` |
| `merge` | 合并多个 BWT（FMR）| `-t` 线程、`-o` 输出、`-S` 中途保存 |
| `kount` | 统计高出现 k-mer | `-k` k-mer 长、`-m` 最小出现（`kount -k 51 -m 100` 类）|
| `fa2line` | FASTX 转行（正反链）| `-R` 不含反向链 |
| `fa2kmer` | FASTX 抽 k-mer | `-k` 长度、`-w` 步长 |
| `plain2fmd` | 纯文本 BWT → FMD | `-o` 输出 |
| `stat` | 报告序列数/符号数/run 数（FMD）| `-M` mmap |

> **注意**：`mem`/`sw`/`hapdiv` 实际共享 `main_search` 入口，靠 `argv[0]` 子命令分发
> （`main.c`），`search` 是其通用别名。`main.c` 共注册 14 个命令。

> **kount 与 pgr 的关联**：`kount` 在 FM-index 上做深度优先遍历统计（出现次数 ≥ `min_occ`）
> 的 k-mer，与 pgr 的 `kmer` 命令功能对应；若未来 pgr 引入 FM-index，kount 的"rank2a 定长
> 区间 DFS + 阈值剪枝"是可直接移植的计数骨架（`main_kount`，约 80 行 C）。

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
   `ScalarAlignmentEngine` 是 O(nm) 内存全矩阵，BWA-SW 的"候选集 + 堆 + 位数组 F"
   是低内存替代（但依赖 FM-index）。
4. **采样 SA 的序列 ID 编码**：SA 值低位存序列 ID、高位存偏移 + `r2i` 哨兵映射，
   是"稀疏坐标 + 序列归属"的紧凑方案，pgr 的 PAF 索引坐标定位可参考。

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
- LICENSE：MIT（licenses.txt）；libsais 为 MIT（IlyaGrebnov）。
- 参考：https://github.com/lh3/ropebwt3 ；相关：fermi/ropebwt2（FMR 格式源）、
  minimap2（BWA-SW 的现代替代）。
