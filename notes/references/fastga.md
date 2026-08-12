# FastGA 源码与论文分析

> 整理于 2026-08，源自对 `FASTGA-main/` 目录源码（全部 `.c` 合计约 4.1 万行，
> `wc -l` 实测 40653 行）及 README 的通读。
> 2026-08-05 复核：`-S`（对称 adaptamer）已由 README 文档化，本机安装版二进制
> 更新为 `[-vkMS]`（与源码 V1.5 一致），§5/§8 相应修正。
> 2026-08-06 通读论文 PDF（Bioinformatics Advances 5(1):vbaf238）全文，
> 补 §12（灵敏度评估 + FastGA-gapfill）、§12.4（DToL 实验口径）及
> §2-§3 若干工程细节（GIX 体量口径、MSD 排序、链编码、merge 算法、
> trace points 参数、gap refinement 论文描述）。
> 目的：理解 FastGA 的快速全基因组比对算法（adaptive seeds + wave aligner + trace points），
> 为 pgr 的 pangenome 上游比对（verify-pangenome.sh 已用 `FastGA -psl/-pafx`）与对齐算法
> 提供参考。

## 1. FastGA 概览

- **工具定位**: 在两个高质量基因组之间（或一个基因组自比对）寻找全部局部 DNA 比对，
  默认输出 PAF，也可输出 PSL 或 ONEcode `.1aln` 格式。
- **作者/版本**: Gene Myers（daligner 作者）与 Chenxi Zhou；2023-05 首次发布，
  V1.5（2025-12-30）为当前版本。
- **核心假设**: 输入为近完整组装（至多几千 contig），序列质量 Q40+。
- **性能**: 2 Gbp 蝙蝠基因组 vs（8 核）约 5 分钟找到几乎所有 >100 bp、≥70% 相似区域；
  63.5 万个比对压缩到 44.5 MB `.1aln`（README 口径）。**论文口径（§摘要）**：
  2.1 分钟 / 8 线程 / 5.7 GB 内存 / 1.05M 条比对 / 覆盖各基因组 60%，
  ALN 66 MB → 1.03 GB PAF 仅 6 s——两处数字统计口径不同，都来自蝙蝠实验。
- **算法来源**: adaptive seed（Martin Frith 的 adaptamer 思想）+ 首个 wave-based
  local aligner（源自 daligner 2012）；数据编码用 Gene Myers 的 ONEcode 框架。

**与 pgr 的关系**：pgr 曾以 `FastGA -v -psl/-pafx A B` 作为 pangenome 上游比对器
（→ `pgr pl chainnet --syn`），FastGA 自己的 chaining 在 pgr 路线中不使用。2026-08
起 pgr 已原生实现同一管线（`pgr align pgi`），FastGA 变为外部对照；移植与实测
对照见 [[../design/pgi-align.md]]。

**实测内存（2026-08-02，MG1655 vs Sakai，-T8）**：FAtoGDB ~7 MB、GIXmake -T8
~160 MB、`FastGA -psl` 比对主进程 **332 MB**。全程无 mmap：GIX 用 `read()` 流式
读、序列默认 EXTERNAL 文件态。

## 2. 整体架构与数据流

```
FASTA/ONEcode
    │  FAtoGDB
    ▼
GDB (.1gdb 元数据 + 隐藏 .bps 2-bit 序列)     ← 随机访问、4 倍省 IO
    │  GIXmake
    ▼
GIX (.gix 稀疏 k-mer 索引, k=40 + (12,8) syncmer) ← 两个索引直接互查
    │  FastGA 主流程（种子扫描 → 链 → wave 对齐）
    ▼
.1aln（ONEcode trace point 编码，按 contig1→contig2→start 排序）
    │  ALNtoPAF / ALNtoPSL（线性时间）
    ▼
PAF / PSL
```

- 所有步骤可由 `FastGA` 一次触发；`FAtoGDB` / `GIXmake` / `ALNtoPAF` 等子进程允许
  分步控制，GDB/GIX 可持久化复用（-k 保留，多基因组重复比对时显著省时）。
- GIX 很大（每 Gbp 约 14 GB），建议批量比对前构建、之后用 `GIXrm` 清理，保留 GDB。
- 自比对模式（`FastGA A`）可检测基因组内部重复/单倍型间同源。
- **方向不对称**：adaptamer 依赖 source1 的种子，`FastGA A B` ≠ `FastGA B A`；
  `-S` 用两个基因组的 adaptamer 做对称（更慢，重复结构分析用；synteny 场景不建议；
  README:199-207 已文档化，见 §5）。
- **Soft mask**（V1.3+）：FASTA 小写=掩码，存入 GDB 的 `.1ano` 文件；默认忽略，
  `-M` 或 `#mask.1ano` 参数启用。

## 3. 核心算法

### 3.1 GDB：genome database（GDB.c）

- **两级结构**：scaffold → contig。`GDB_SCAFFOLD` 记录 scaffold 长度、首/末 contig、
  header 偏移；`GDB_CONTIG` 记录 contig 长度、scaffold 内起始、`.bps` 文件字节偏移。
- **序列存储**：2-bit 压缩（每个碱基 2 bit），存于隐藏文件 `.foo.bps`；元数据与序列分离，
  不需要序列的应用只读轻量 `.1gdb`。
- **N 处理**：FASTA 中 N 默认为 contig 间 gap；`-n` 指定阈值，短于阈值的 N 视为未知碱基
  （按 'a' 处理）。
- 派生自 daligner 的 GDB 代码。

### 3.2 GIX：syncmer 稀疏 k-mer 索引（GIXmake.c）

- 对每个 GDB 构建 k-mer 索引（`-k` 默认 40），但**不是全后缀数组**：只索引"以
  (12,8) syncmer 起始"的 40-mer（GIXmake.c: `TMER=12, SMER=8, SOFF=4`）。
  **实际检测是分布线程内的内联扫描**（GIXmake.c:220-311 / 480-601）：对 12-mer
  窗口内 **5 个重叠 8-mer**（TMER−SMER+1=5；8-mer 由两个相邻 4-mer 半字拼成，
  值 = `(nq[jq]<<8)|nh`，16-bit）取 canonical（正反链 `Comp`/`TMap`）最小值，
  当最小 8-mer 落在窗口**左端点**（`pos4 == i-SOFF`）时触发，40-mer 锚在
  `j = i-SOFF`（GIXmake.c:569）。
  **独立的 `is_syncmer` 函数（GIXmake.c:131）是 `#ifdef DEBUG_SYNCMERS` 下的
  死代码**，正常构建不参与（早前笔记误引为活动检测逻辑，已更正）。
- 每个索引条目 = 40-mer + 位置 + 掩码前缀信息；排序表存为 `-T` 个隐藏
  `.ktab.<int>` 分片（`.gix` 只是代理文件）。
- **新格式条目布局**（FastGA.c:5012-5016，`-DLCPs` 构建）：`kmer 后缀(KBYTES)`
  + `mask/len 字节(CBYTE)` + `lcp 字节(LBYTE = CBYTE+1)` + `contig/post 载荷
  (PAYOFF = LBYTE+1)`；GIXmake.c:1350 条目宽 `swide = MBYTES+PostBytes+ContBytes`，
  `MBYTES = KBYTES+1`。排序在 GIXmake 内调 `msd_sort`（MSDsort.c），GIXmake
  经 Makefile 以 `-DLCPs` 编译（记录相邻条目 lcp 到每条记录首字节）。
- 排序表按首字节 1024 桶（`NUM_BUCK`）预分 + `Ksplit` 均衡分片后多线程桶排序。
- **体量**：README 实测约 14 GB / Gbp；**论文口径 ~11 GB / Gbp、构建约
  15 s/Gbp**（M4 Max 16 核、8 线程）——两个数字统计口径不同，都远大于
  序列本身；人类 ~3 Gbp ≈ 33-42 GB。FastGA **默认在退出时自动删除自己
  创建的 GIX/GDB**（`Clean_Exit` 调 `GIXrm`，`-k` 才保留），运行时索引落在
  `TMPDIR` 或 `-P` 指定目录。
- 关键设计：FastGA **直接比较两个 GIX**（两个排序的 k-mer 位置流线性归并找相同
  40-mer），而不是把一方的序列在另一方的索引中逐条查询——这是速度来源之一。
- 索引用 2-bit 编码 + canonical 方向，正反链统一。
- **论文的索引模型**：排序 40-mer 表 + lcp 数组 = "截断到深度 K 的后缀数组"
  （truncated suffix array）；`(12,8)` syncmer 把 40-mer 数量削减一半以上，
  adaptamer 短于 s=12 的种子会丢失（论文认为此类种子与真实匹配相关性低）。
- **论文的 MSD 排序工程细节**（pgr 版 `radix_sort_u128_par` 的取舍对照）：
  1. 顶层按**前 4 bp 预分到 T 个文件**（各文件区间有序，排序后直接拼接），
     文件内只存相对上一个 40-mer 的**位置差**（<2 B/条）——GIX 磁盘格式；
  2. 排序时记录**相邻条目的 lcp 到每条记录首字节**（GIX 记录 = kmer 后缀 7 B
     + mask 1 B + lcp 1 B，见 [[../design/pgi-align.md]] §7.6）；
  3. 空分区跳过（随深度增加空桶渐多）、已就位条目不搬（省 ~10% 移动）、
     小分区退化为更简单排序。
  pgr 的 radix_sort 已含：公共前缀跳过、空分区跳过、cycle 置换、insertion
  fallback、并行桶排序（`src/libs/ds/radix_sort.rs`）；**未做**前 4 bp 预分片
  与 lcp 入字节（GIX 私有格式细节，pgr 的 `.pgi` 不需要，见 pgi-align.md §7.6）。

### 3.3 Adaptive seeds（adaptamer，libfastk.c / FastGA.c）

- **定义**：位置 p 的 adaptive seed = 从 p 开始、在另一基因组中也出现的最长字符串。
- **频率过滤**：若该字符串在 source2 出现次数 > `-f`（默认 10），视为重复、不作为种子。
- **最小化**：`is_minimal`（libfastk.c:590）把种子与其反向互补做字典序比较，保留更小者
  （canonical 方向），正反链统一——与 pgr 的 canonical minimizer/syncmer 思路一致。
- **种子命中**：adaptamer 在 source2 的每个出现位置 (p, q) 都是一个 seed hit。
- **论文的归并算法（§3.1, Algorithm 1）**：对 G 的每个 K-mer 维护
  `(fst, lst, L, cur, wall)` 五个 O(K) 辅助量——`[fst, lst)` 是 H 中与该
  K-mer 共享前缀 ≥ L 的条目区间、`cur` 是区间内第一个 ≥ α 的条目、`wall[l]`
  是共享前缀 ≥ l 的最小下标；相邻 G 条目靠 lcp（λ）增量推进，不用逐条目
  二分。Theorem 1：单趟线性归并找出 G 中全部非重复 adaptamer。实测（蝙蝠，
  2.23 Gbp vs 2.56 Gbp）：1.603B 条 syncmer 40-mer，交叠 2.646B 边、
  4.356B 次字符比较（**每 adaptamer < 3 次**），区间 `[fst,lst)` 通常大小为 1。
  pgr 的简化版（§3.3 lcp 起步窗口）语义等价、E. coli 上已无差距（见
  [[../design/pgi-align.md]] §3.3/§7.6）。
- `-S` 对称模式取两个基因组的 adaptamer 并集（`SYMMETRIC = flags['S']`，merge
  互换 T1/T2 跑双向再合并；README:199-207 已文档化，见 §5）。

### 3.4 种子链（chaining，FastGA.c align_contigs）

**术语：Tube**（出处 FastGA.c:3160，`align_contigs` 注释）——一个覆盖
≥ `CHAIN_MIN` 的种子链连同它的比对盒：anti 区间 `alow..ahgh` × 对角线带
`dgmin..dgmax`。Tube 不是正式 API/类型（源码无 `Tube` 结构、无
`tube_*` 函数），只是 Myers 在注释与调试宏（`DEBUG_TUBE`，FastGA.c:35/
3314，输出 "Did not reach top/bottom"）里对"链 + 盒"的称呼；它是给
wave 划定的搜索区域。pgr 移植时把 tube 提升为正式概念（`Tube` 结构体、
`chain_tubes`/`extend_tube`），语义同源。

种子 hit 按 **anti-diagonal（反对角线）空间** 排序扫描（`print_seeds` /
`align_contigs` 中维护按 `(ipost, apost)` 归并的种子流；注意这两个变量名
**装的是 anti**——reimport 时 `memcpy(_anti, ...)`，排序/归并键是
`anti = i+j`，不是坐标 `i`）。合法链需满足：

1. 所有种子落在宽度 128 的对角线带内；
2. 相邻种子间距 < `CHAIN_BREAK`（源码 2000 = 2×`-s` 1000，anti-diagonal 空间 2 倍）；
3. 链在两侧覆盖 ≥ `CHAIN_MIN`（源码 170 = 2×`-c` 85）个 anti-diagonal。

满足条件的链在"tube"（`alow..ahgh` × `dgmin..dgmax`）内触发 wave aligner
（`Local_Alignment`）。self 比对（`SELF && ctg1==ctg2 && !comp`，FastGA.c:3030
判定 `self`，对角线限制调用在 3245-3257；早前笔记误记 3220-3240，已更正）
对 tube 扩展做对角线限制：tube 带全正（`dgmin > 0`）时
`Local_Alignment(..., dgmin-1, -1)`、全负（`dgmax < 0`）时
`Local_Alignment(..., -1, -(dgmax+1))`，带跨 0 的 tube 整管跳过——
`Local_Alignment` 内 `minp = low-lbord` / `maxp = hgh+hbord` 成为 wave
`forward/reverse` 的对角线硬边界（`low >= minp` / `hgh <= maxp` 分支），
保证路径不跨越精确自同线（diag 0）。

**论文的链编码（§2.4）**：每个种子 hit 编码为六元组
`(ci, cj, ⌊(pi−pj)/D⌋, pi+pj, (pi−pj) mod D, t)`——对角线桶 b、anti 值、
桶内余数可重建两侧坐标；MSD 排序后同 contig 对、同对角线桶的种子按 anti
连续，相邻桶对（b, b+1）归并扫描即可覆盖宽 ≤ 2D 的链（D=64 编译期常量、
A=1000 为 anti 间距参数）。**冗余比对去除**：等价/包含（bounding box 内含）
的比对只留一个，输出按 source1 起始坐标排序——pgr 侧对应
`dedupe_contained`（0.95 阈值）与 PSL 排序输出（pgi-align.md §1.3.3）。

### 3.5 Wave-based local alignment（align.c）

**术语：Wave**——Myers wavefront 算法的波前：`V[k]` 是"编辑距离恰好为
d 时对角线 k 上的最远到达点"，随 d 逐波扩展（每波对角线 ±1）。"wave"
是算法本身的形象叫法（调试宏 `DEBUG_WAVE`/`SHOW_MATCH_WAVE`，
align.c:29-30），不是数据结构的名字。FastGA 用它做局部比对的两段式：
`forward_wave`（从 mid-line 锚点向高 anti 延伸）+ `reverse_wave`（镜像
序列上向低 anti 延伸，等价于原坐标反向），端点框定比对区间，再用
`dandc_nd` 回溯精确路径。

源自 daligner 的 wave-front 对齐：

- **forward_wave / reverse_wave**（align.c:352 / align.c:878）：沿对角线扩展 wave，维护
  - `V[k]`：对角线 k 的最远到达点（furthest reaching point）；
  - `M`：最近 TRIM_LEN 列的匹配数（位向量 1-bit 计数）；
  - `T`：隐含对齐最后列的位向量（用于轨迹）；
  - `Pebble` cells：wave cell 记录（ptr/diag/diff/mark）。
- **Local_Alignment**（align.c:1423）：在给定对角线带与 anti-diagonal 区间内做局部对齐，
  长度 ≥ `-l`（默认 100）、相似度 ≥ `-i`（默认 70%）。
- **Compute_Alignment**（align.c:5426）：divide-and-conquer trace（`dandc_nd` /
  `trace_nd` / `middle_np` / `iter_np`），用 sparse DP 在 wave 之间回溯完整比对路径，
  按 `tspace`（trace spacing）压缩轨迹。
- **Gap_Improver**（align.c:6714）：对 gap 区域做二次精修。论文 §3.2.2 给出
  算法描述：把 trace 转成 indel 数组，找同号且相邻值差 < R=50 的极大段
  （即一个 gap 簇）→ 在 D×L 梯形区域内用**压缩 wave**（gap 起始代价 1、
  延续 0，Hadlock 路径压缩）重算最少 gap 路径，O((G+S)D+L) 期望时间。
  实测蝙蝠：初始 64.6M 个 gap，检查 11.93M 个梯形，移除 13.7M 个 gap，
  平均 (G+S)D=88，精修约占 ALN→PAF 转换时间的 15%。**pgr 不接入**：
  该精修只对 tspace 采样产生的次优 trace 有效，pgr wave 用 `dandc_nd`
  全量精确回溯、无采样缺口，源码语义下必然无操作（完整论证见
  [[../design/pgi-align.md]] §3.1）。

> **术语对照（Wave / Tube / Cube）**：源码里只有 Wave（算法，`forward_/
> reverse_wave`、`Local_Alignment`）与 Tube（种子链 + 比对盒，注释语）。
> **没有 "Cube"**——FastGA 源码与文档均无此词（全仓 `rg -i cube` 无
> 结果）；若在别处看到 Cube，应为 Tube 的误记（t/c 同音或拼写混淆）。
> pgr 代码与笔记沿用 Wave/Tube 两个词，与上游一一对应。

### 3.6 Trace points 与 .1aln 编码（alncode.c / ONEaln.c）

- 每条比对记录为**轨迹点（trace points）**：按 tspace 采样的比对路径编码，配合
  diff/长度信息，在 ONEcode 二进制中极紧凑（63.5 万比对 → 44.5 MB）。
- **论文的编码模型（§3.2.1）**：A 侧按 δ=100 bp 面板划分，trace-point 数组 =
  每个面板对应 B 侧面板长度（bi 通常 ∈ [δ(1−ɛ), δ(1+ɛ)]，单字节）；配合四
  端点与 δ 可重建整条比对，**重建 O(n+δd)**（逐面板用 Wu et al. 1990 的
  skewed-wave）。1 万 bp 比对 ≈ 100 字节（与 ɛ 无关，CIGAR 要到 ɛ<1/100
  才同规模）；ALN→PAF 展开约 1 亿 aligned bases/s。
- **trace spacing 默认 `TSPACE = 100`**（FastGA.c:46 `#define TSPACE 100`；
  `New_Align_Spec(1.-ALIGN_RATE,100,...)`，FastGA.c:3760）——即每约 100 个对齐列采一个
  trace point。
- `.1aln` 头：schema 行 **`1 3 def 2 1`**（alncode.c:20 `alnSchemaText`，早前笔记误记为
  `1 3 aln 2 1`，已更正）、`!` 记录 FastGA 版本与参数（`oneAddProvenance`）、`<` 引用
  两个 GDB（`oneAddReference`，第三个引用是 cpath）、`t` 行记录 tspace（alncode.c:265-266）。
- **记录编码**（alncode.c `Write_Aln_Overlap`/`Write_Aln_Trace`）：每条比对一个 `A` 行
  （6 个 int：aread/abpos/aepos/bread/bbpos/bepos）+ 互补时 `R` 行 + `D` 行（diffs）；
  trace 用 `T` 行（逐 trace point 的坐标增量，取 trace 奇数下标）+ `X` 行（对应区间
  diff 数，偶数下标）交替编码。
- **排序**：按 source1 contig # → source2 contig # → source1 start 排序，便于线性扫描。
- ALNtoPAF / ALNtoPSL（多线程）在**线性时间**把轨迹展开为 PAF/PSL（含 CIGAR：
  `-pafx` = `=`/`X`，`-pafm` = `M`；`-pafs/S` = CS 字符串）。
- ONEaln.c 提供 C 库读取 .1aln（依赖 GDB/ONElib/alncode/align/gene_core 一起编译）。

### 3.7 ALNchain：基于 KD-tree 的局部链过滤（ALNchain.c，Chenxi Zhou）

`ALNchain` 是独立于主流程的**后处理工具**：读入一个 `.1aln`，只保留"最佳局部链"
里的比对（丢弃被更高分链覆盖/包含的冗余比对），输出过滤后的 `.1aln`
（默认 `<root>.chain.1aln`）。它是"从全比对集合提炼共线骨架"的参考实现。

- **数据结构**：比对盒（TNODE：a/b 侧 beg/end + score/clen）组织成 **KD-tree**，
  交替按 X=`aepos`、Y=`bepos` 分轴（`buildKDTree`，ALNchain.c:199），中位数用
  `quickSelect` 的 **median-of-medians** 保证 O(n) 最坏情况（`selectPivot`，
  ALNchain.c:162，`USE_MEDIAN_MEDIAN` 默认开）。
- **链式 DP（`localChain`→`KDRangeChain`，ALNchain.c:336）**：对每个节点做 KD-tree
  区间查询找前驱；节点作为前驱需满足 active 且两侧 ext>0、gap ≤ maxGap（默认 10000）、
  ovl ≤ maxOvl（默认 10000）、重叠小于查询盒自身跨度。得分
  `score = extX+extY − gap·penGap − ovl·penOvl`（penGap=penOvl=.10），保留
  `root->score + score` 最大的前驱。gap 即"相邻比对间的间距"、ovl 即"重叠量"，
  二者都按 X、Y 两侧分别计（gap<0 转成 ovl）。
- **贪心选链（`popLocalChain`→`backtrackLocal`，ALNchain.c:447）**：按 score 降序
  逐个回溯；链断裂条件为后继 `score > maxDrop + minScore`（`-z` maxDrop 默认 1000）
  或后继已被占用。链保留条件：`score ≥ 2·minScore`（`-s` 默认 10000，链得分按 X+Y
  双倍计故 ×2）且成员数 `clen ≥ minFrag`（`-n` 默认 1）。
- **跨链去冗余（`filterChain`，ALNchain.c:518）**：同一 scaffold 对内按 score 取主链
  （负链记录先 `reverseRangeStrand` 把 Y 侧翻正），其余链与其累计覆盖做 fuzzy merge
  （`-f` fzMerge 默认 1000 容差），若**两侧**覆盖率都 > `maxCov`（`-c` 默认 .50），
  或非重叠延伸 < `minExt`（`-e` 默认 0.0，占序列长度比例），则该链被过滤（置
  INTERNAL）；主链的覆盖随保留链累积，保证贪心但单调。
- **CLI**（ALNchain.c:48-52 Usage）：`-g` maxGap=10000、`-l` maxOvl=10000、
  `-p`/`-q` penGap/penOvl=.10、`-c` maxCov=.50、`-e` minExt=0.0、`-s` minScore=10000、
  `-n` minFrag=1、`-z` maxDrop=1000、`-f` fzMerge=1000、`-o` 输出、`-v`；输入为
  `.1aln` 路径。
- **与 pgr 的对照**：pgr 的 pgi/chainnet 里"保留主链 + 按覆盖去冗余"（如
  `dedupe_contained`、synteny 骨架挑选）思路同源，但 pgr 用的是区间树/贪心排序，
  ALNchain 用 KD-tree 区间查询做前驱 DP——KD-tree 在多维（X/Y 同时约束 gap 与 ovl）
  前驱查询上是更完整的形态，pgr 若需"双轴 gap+overlap 约束的链式 DP"可参考此实现。

## 4. 源码模块结构

| 模块 | 职责 |
|------|------|
| `FastGA.c`（~5.3k 行）| 主流程、参数解析、种子扫描、anti-diagonal 链、调用 aligner |
| `align.c`（~7.1k 行）| wave aligner（forward/reverse_wave、Local_Alignment）、Compute_Alignment、trace、Gap_Improver |
| `MSDsort.c`（~0.5k 行）| MSD radix 排序 `msd_sort`（GIXmake 的 k-mer 桶内排序用）|
| `RSDsort.c`（~0.4k 行）| radix 排序 `rmsd_sort`（FastGA 主流程种子流排序用，FastGA.c:149 声明）|
| `libfastk.c` / `FastKS.c` | FastK 生态 k-mer 计数库：读写 GIX 的 `.ktab` 表（Histogram / Kmer_Table / Kmer_Stream / Profile_Index）|
| `GDB.c` / `GDB.h` | genome database：scaffold/contig 两级结构、2-bit 序列随机访问 |
| `GIXmake.c` | syncmer 稀疏 k-mer 索引构建（k=40 + (12,8) syncmer，含 mask 支持）|
| `ONElib.c` / `ONEaln.c` | ONEcode 数据编码框架、.1aln 读取 C 库 |
| `alncode.c` | trace point 编解码 |
| `ALNtoPAF.c` / `ALNtoPSL.c` | 轨迹 → PAF/PSL（多线程、含 CIGAR 生成）|
| `ALNchain.c`（Chenxi Zhou）| 按 KD-tree 局部链过滤 .1aln 比对（见 §3.7）|
| `ALNreset.c` | 重设 .1aln 对 GDB 的内部引用 |
| `select.c` | 基因组选择表达式解析（**仅供展示/绘图工具** GDBshow/ALNshow/ALNplot/ANOshow 用，不参与 FastGA 主流程比对选择；语法：`@`scaffold、`.`contig、`:`position、`#`last、`-`range、`,`分隔，位置可带 G/M/k 后缀）|
| `PAFtoALN.c` / `PAFtoPSL.c` | 反向转换（PAF 带 X-CIGAR → .1aln/.psl）|
| `FAtoGDB.c` / `GDBtoFA.c` | FASTA/ONEcode ↔ GDB 互转 |
| `GDBshow`/`GDBstat`/`ALNshow`/`ALNplot`/`ANOshow`/`ANOstat` 等 | 查看/统计/绘图工具 |

> **libfastk 移植评估**：libfastk 是 FastK/GIX 私有格式（`.ktab.<int>` 分片）的访问库，
> 含 Histogram（k-mer 频率直方图）、Kmer_Table（加载排序表 + Fetch/Find）、
> Kmer_Stream（流式遍历 + GoTo 定位）、Profile_Index（raw reads profile）。pgr
> **不需要移植任何实现**：格式绑定 FastK 生态、pgr 不做 k-mer 计数/raw reads；
> 其中 `is_minimal`/`compress_norm`/`compress_comp`（2-bit 编码 + canonical 最小化）
> 与 pgr `nt.rs`/`syncmer.rs` 等价。仅 Kmer_Stream 的"排序流迭代 + 定位"接口形态
> 值得未来原生 `sd search --mode kmer` 的两流归并借鉴（Rust 自研实现）。

## 5. 关键参数（FastGA main 默认值，源码与 README 的对应）

| 参数 | 默认 | 含义 | 源码位置 |
|------|------|------|----------|
| `-f` | 10 | 最大种子频率（超过视为重复，不作为 adaptamer）| `FREQ = 10` |
| `-c` | 85 | 最小链覆盖 bp（源码 `CHAIN_MIN` 存 2×=170，anti-diagonal 空间）| `CHAIN_MIN = 170; <<= 1` |
| `-s` | 1000 | 相邻种子最大间距（源码 `CHAIN_BREAK` 存 2×=2000）| `CHAIN_BREAK = 2000` |
| `-l` | 100 | 最小局部比对长度 | `ALIGN_MIN = 100` |
| `-i` | 0.7 | 最小比对相似度（源码 `ALIGN_RATE = 1.-sim`，默认 .3；合法 [0.55,1)）| `ALIGN_RATE = .3` |
| `-k` | 40 | GIX k-mer 大小（GIXmake）| — |
| `-T` | 8 | 线程数 | `NTHREADS = 8` |
| `-S` | off | 对称 adaptamer（两个基因组的种子并集）；**README:199-207 已文档化**——用两方 adaptamer，稍慢但结果与 A/B 顺序基本无关；通常发现 B 中更多重复比对；**synteny 场景不建议，仅重复结构分析时用** | flags（`ARG_FLAGS("vkMS")`） |
| `-M` | off | 使用 GIX 中的 soft mask | flags |
| `-v` / `-L` | — | 详细模式 / 日志文件 | flags |

> **源码注释过时**：FastGA.c 中 `ALIGN_MIN` / `ALIGN_RATE` 的声明注释分别标 `// -a`、
> `// 1.-e`，与实际参数 `-l`（min alignment length）、`-i`（min similarity，内部折成
> `ALIGN_RATE = 1.−sim`）不符——是历史遗留注释，读源码勿据注释推断参数名，应以
> `main()` 的 `switch`（case `'c'/'f'/'i'/'l'/'s'`）为准。

## 6. 输出格式

- **PAF**（默认）：12 列标准 PAF + **恒附加两个 SAM tag**（ALNtoPAF.c:466/474）：
  `dv:f:`（query 相对 target 的分歧度 fraction）与 `df:i:`（最优比对差异数）——这是
  FastGA 自带的轻量分歧度口径，见 §7 第 6 条；`-pafx` 再追加 `cg:Z:`（`=`/`X`/`I`/`D`），
  `-pafm` 用 `M`，`-pafs/S` 追加 CS 字符串。
  注意 README 明确 `-m/-x/-s/-S` 会令 ALNtoPAF 时间 ×10、输出体积 ×~100（CIGAR/CS 展开
  开销）——pgr 的 PAF CIGAR 懒加载（BGZF vpos）正可规避此成本。
- **PSL**（`-psl`）：UCSC PSL 格式，可直接喂给 `pgr pl chainnet` / `pgr psl chain`。
  负链块约定：qStart/qEnd 用正链帧、内部 qStarts 用 RC 帧（`ALNtoPSL.c` 对 COMP
  记录按 `blen − pos` 反算）。
- **.1aln**（`-1:path`）：ONEcode 二进制（须指定输出文件），可用 ALNtoPAF/ALNtoPSL
  按需转换。

## 7. 对 pgr 的启示

1. **Adaptive seeds vs pgr 的 syncmer/minimizer**：FastGA 的 adaptamer 是"最长共享字符串"
  （长度自适应），pgr 的 closed syncmer（`src/libs/syncmer.rs`）是固定 k 的有界间隔采样。
  FastGA 的 `is_minimal` canonical 判断与 pgr 的 canonical rolling hash 同思路
  （`is_minimal` 是 canonical 方向判断，不是噪声抑制）。pgr 的 pgi 比对管线已移植
  其 merge 种子语义（见 [[../design/pgi-align.md]] §1.3.2）。
2. **Wave aligner vs pgr 的 ScalarAlignmentEngine**：pgr 的 POA 对齐是标量 O(nm) 矩阵 DP；
  FastGA 的 wave-front（V/M/T 位向量 + Pebble cells）在线性空间内扩展，与 WFA 同族。
  pgr 已移植该 wavefront 用于 `pgr align pgi` 的 tube 扩展（wave 依赖 tube 锚定
  上下文，不能单独替换 banded，见 [[../design/pgi-align.md]] §3.5.1）。pgr 的 POA
  对齐仍是标量 O(nm)，wave-front 对 pbit CIGAR 精修 / SD refine 仍是候选优化方向。
3. **Trace point 编码 vs pgr 的 PAF/MAF 存储**：FastGA 用轨迹点 + ONEcode 压缩比对集合
  （63.5 万比对 → 44.5 MB），支持线性时间重放为任意格式。pgr 的 paf index 已做
  CIGAR 懒加载（BGZF vpos），但完整比对的紧凑存储可借鉴 trace point 思路。
4. **GDB 2-bit 序列库 vs pgr 的 twobit**：FastGA 把 2-bit 序列存隐藏文件 + 元数据分离，
  与 pgr `TwoBitFile` 类似；pgr 的 `.loc`/BGZF 随机访问已覆盖同等需求。
5. **对称性语义**：`FastGA A B` ≠ `FastGA B A`（adaptamer 不对称）对 pgr 有直接影响——
  pangenome 管线里 FastGA 的 query/target 顺序会影响找到的比对集合；
  verify-pangenome.sh 固定 `FastGA(b,a)` 方向后 chainnet 统一精修，顺序影响被下游
  chain/net 部分吸收（但重复区域仍可能不对称）。`-S` 可消除顺序依赖，但 README
  明确"**synteny 场景不建议**，仅理解两基因组重复结构时用"——pgr 的 `align pgi`
  （synteny 用途）因此不实现 `-S`（见 [[../design/pgi-align.md]] §7.4）。
6. **轻量分歧度 tag**：ALNtoPAF 恒输出的 `dv:f:`（分歧度 fraction）/`df:i:`（差异数）
  两个 SAM tag（ALNtoPAF.c:466/474）是很轻量的分歧度口径，无需展开 CIGAR 即可随
  PAF 携带。pgr 的 PAF 输出（`pgr align` / paf 工具族）若需在比对时顺带标注分歧度，
  可参考该编码——比每次都解析 `cg:Z:` 便宜得多。

## 8. 版本与许可

- 当前 FASTGA-main 对应 V1.5（2025-12-30），含 ONEcode ANO 文件支持。
- 注意 FastGA.c:14 内部 `#define VERSION "0.1"` 是过时的占位字符串，不代表发布版本
  （README/发布版本为 V1.5），读源码勿据此判断版本。
- **版本核对（2026-08-05）**：本机安装版二进制（`~/.cbp/bin/FastGA`）帮助为
  `[-vkMS] [-L:<log:path>] [-T<int(8)>] [-P<dir($TMPDIR)>] [<format(-paf)>]`——
  支持 `-S`/`-M`/`-L`，与仓库源码 V1.5（`ARG_FLAGS("vkMS")`）一致；`-S` 语义见
  README:199-207（§5）。早前记录"安装版为 `[-vk]`、拒绝 `-S`/`-M`"已随二进制
  更新而过时，作废。
- LICENSE：MIT（ALNchain 单独标注 Chenxi Zhou，MIT）。
- 参考：https://github.com/thegenemyers/FASTGA ；ONEcode:
  https://github.com/thegenemyers/ONEcode ；daligner:
  https://github.com/thegenemyers/DALIGNER

## 9. GDB 与 pgr 存储格式对比

### 9.1 对比对象

| 格式 | pgr 侧实现 | 定位 |
|------|-----------|------|
| FastGA GDB（`.1gdb` + `.bps`）| — | 组装基因组数据库：元数据与 2-bit 序列分离 |
| pgr 2bit（`pgr fa to-2bit`）| `src/libs/fmt/twobit.rs` | 标准 UCSC 2bit：单级 contig + 内嵌 mask/N block |
| pgr loc 索引 FASTA | `src/libs/loc.rs` | 原样 FASTA/BGZF + `.loc` 偏移索引（`fa range`、paf `FastaStore`）|
| pbit 参考层 | `src/libs/pbit/` | 复用 2bit 记录格式（`read_2bit_record` / `write_2bit_record`）|

### 9.2 序列编码与空间效率

- **GDB**：2-bit 压缩（`COMPRESSED_LEN = ceil(len/4)`），编码表 **A=0, C=1, G=2, T=3**
  （libfastk.c 的 `code[128]`），N 在 FASTA 阶段即按 gap 拆分（`-n` 阈值）。
- **pgr 2bit**：同为 `ceil(len/4)` 字节，但编码表是 **UCSC 标准 T=00, C=01, A=10, G=11**
  （twobit.rs:118）。N 保留为 n_blocks（长度不变）。
- **空间结论**：两者序列密度等价（0.25 B/bp）。pgr loc+FASTA 是 1 B/bp（4 倍），
  但保留原文、无转换成本。
- **编码表差异是硬约束**：GDB 的 packed bytes 与 pgr 2bit 的 packed bytes 不能互读，
  且负链互补映射不同（GDB `comp` 表 vs pgr 的 T↔A、C↔G）。若要互操作必须走文本/ASCII
  再重新编码，不能直接搬 packed 数据。

### 9.3 结构与元数据模型

| 维度 | GDB | pgr 2bit | pgr loc |
|------|-----|----------|---------|
| 层级 | scaffold → contig 两级（N 即 gap）| 单级 contig 平铺 | 单级 FASTA record |
| 元数据 | 与序列**分离**（轻量 `.1gdb` 不需 `.bps`，只读骨架免载序列）| 一体（记录头含 dna_size + block 表）| 一体（原文）|
| N/gap | N 拆分 contig（组装语义，scaffold 保留 N 的 gap 长度）| N 保留为 n_block（序列语义，长度不变）| 原样保留 |
| mask | **外置** `.1ano` 区间文件（可多个 mask union，改 mask 需重建 GIX）| 内嵌 mask_blocks（保留 soft-mask 语义）| 原文大小写 |

- **GDB 的 scaffold 语义**对组装输入（contig + gap 长度估计）更友好；pgr 2bit 面向
  "序列就是序列"的通用场景，N 当未知碱基。
- **GDB 元数据/序列分离**是实际优势：统计 scaffold 数、长度、名称时不用碰 `.bps`。
  pgr 2bit 单文件一体化，读元数据要扫整条记录。
- **mask 模型**：GDB 外置 `.1ano` 可 union 多个 mask 且不改序列本体（但改 mask 必须
  重建 GIX）；pgr 2bit 内嵌（一次打包、不可变）。pgr 的 `pgr fa mask`（runlist 硬/软
  mask，长度保留）接近"外置改字符"，但语义不同。

### 9.4 随机访问

- **GDB**：`Get_Contig_Piece` 按 `boff + beg/4` 做 `fseeko`，读 `ceil((end-beg)/4)` 字节后
  `Uncompress_Read`；`seqstate != EXTERNAL`（`Load_Sequences` 整库 `Malloc`+`read` 载入）
  时 `Get_Contig` 直接 memcpy。区间读取是 **O(区间长)**。
  FastGA 全项目无 mmap 调用，序列访问只有 fseeko 与整库读入两种。
- **pgr 2bit**：`read_2bit_record` 同样 seek 到 `packed_dna_start + first_byte_idx`，
  只读区间字节后解压，**O(区间长)**——与 GDB 等价。
- **pgr loc**：`fetch_range_seq` 先 `fetch_record` **读整条 record** 再切片
  （除非 BGZF 虚拟位置做细分）。对长 contig 的短区间访问，loc 明显低效；
  这是 pgr 内部"区间提取优先用 2bit"（chain ScoreContext、pbit）的原因。
- **构建成本**：三者都是单遍扫描 FASTA；loc 最便宜（只记偏移），GDB/2bit 需编码。

### 9.5 输出形态与工具生态

- **GDB**：`Get_Contig` 支持 COMPRESSED / NUMERIC(0-4) / LOWER_CASE / UPPER_CASE 四种
  形态按需转换；ONEcode 生态（.1seq 输入、.1aln 输出）；`GDBtoFA` 可逆回 FASTA。
- **pgr 2bit**：`read_sequence(no_mask)` 控制大小写；标准 UCSC 2bit 可被 kent-tools /
  其他工具直接读取（pgr 与 UCSC 字节级一致，见 [[ucsc.md]]）；pbit 参考层复用同一记录格式。
- **pgr loc**：底层是标准 FASTA/BGZF，任何工具可读原文。
- **标准性**：pgr 2bit 互操作最强（UCSC 标准）；GDB 是 FastGA 生态私有格式，外部工具
  无法直接消费，属于生态锁定。

### 9.6 结论与建议

- **核心等价**：GDB 与 pgr 2bit 在"2-bit 压缩 + 字节偏移随机访问"上是同构设计，
  空间与区间读取性能相当。真正的差异在语义层（scaffold 两级 vs 平铺、mask 外置 vs
  内嵌、元数据分离 vs 一体）和互操作（私有 vs UCSC 标准）。
- **GDB 值得借鉴的两点**：
  1. **元数据/序列分离**：pgr 若需要"只看骨架不读序列"的场景（如大量基因组的名/长
     统计、`fa size` 的快速版），可参考 `.1gdb` 轻量件设计；
  2. **mask 外置（.1ano）**：pgr 的 mask 内嵌 2bit 后不可变，若要支持"同一参考不同
     mask 重复比对"，外置区间文件 + 读取时过滤（类似 GDB 的 mask 参数）更灵活。
- **不建议**：为互操作而模仿 GDB 的编码表/ONEcode——pgr 2bit 已是 UCSC 标准，
  与 kent-tools 字节级兼容是现有资产（ucsc pipeline 验证依赖它），不应为 FastGA
  私有格式放弃。
- **pgr loc 的定位**：适合"保留原文 + 便宜索引"（fa range、FastaStore），区间密集
  随机访问场景应继续用 2bit。

## 10. GIX 分析：好处与 pgr 借鉴评估

### 10.1 GIX 是什么

GIX（GIXmake.c）是每个基因组的**syncmer 稀疏 k-mer 索引**：只取"以 (12,8) canonical
syncmer 起始"的 k=40 的 k-mer（2-bit 压缩 10 字节），按字典序桶排序（首字节 1024 桶 +
`Ksplit` 均衡分割 + 多线程），排序后每个条目附位置信息（contig #、contig 内偏移、
方向位）与 lcp。README 实测体量约 **14 GB / Gbp**（`.gix` 代理 + `-T` 个 `.ktab.<int>` 隐藏
分片），但 **FastGA 默认退出时自动删除**（`Clean_Exit` → `GIXrm`，`-k` 才保留）；
构建/运行的临时分片在 `TMPDIR` 或 `-P` 目录。

> **与 pgr 的直接关联**：GIX 的 (12,8) syncmer 就是 [[syng.md]]/pgr 已实现的 closed
> syncmer 同族采样——FastGA 用 syncmer 稀疏化 40-mer 索引（密度约 2/(w+1)，大幅低于
> 全后缀数组），这正是 §10.4 建议 pgr"用 syncmer 稀疏化借鉴"的依据：**FastGA 自己
> 就是这么做的**，pgr 不必重新发明。

> **锚定差异（2026-08-07 核对确认，勿混淆）**：尽管同属"检测窗口最小 s-mer 在端点"
> 的 syncmer 族采样，GIX 的 match-mer 与 pgr 的 closed syncmer **锚点不同**：
> - GIX（`GIXmake.c`）把 seed k-mer 锚在**窗口起点**（`j = i-SOFF`，`GIXmake.c:569`）；
> - pgr（`syncmer.rs` / `pgi build`）锚在**最小 s-mer 端点**（为链对称性 / Mash-Jaccard，
>   见 `syncmer.rs` 注释）。
>
> 实证（smer=8/window=5）：两者采样密度相同（GIX=235258 vs pgr=235258），但 **~61%
> 位置锚点不同**，最大偏移 8 bp。端到端验证（self 比对 + 含 ~2% 替换/indel 的两基因组
> 比对）表明**最终对齐结果一致**（同 block 数、同覆盖度），锚点差异只影响 seed 层噪声，
> 不影响 `pgi align` 输出。`pgi build` 的 CLI 帮助与 `docs/pgi.md` / `docs/align-pgi.md`
> 均已注明此差异。

> **为何用户"排人类基因组没感觉生成多大的数据"**：FastGA 默认（不加 `-k`）在
> 结束时 `GIXrm` 清理自己创建的 GDB/GIX，只留下输出 PAF/PSL；GIX 分片只存在于
> 运行期（`TMPDIR`/`-P`，人类 ~42 GB），结束后即被删除。只有显式 `-k` 或预建
> GIX（`GIXmake`）才会留下这几十 GB。

### 10.2 GIX 的核心好处

1. **两个索引线性归并找同源（最重要）**：FastGA 不把 A 的序列在 B 的索引里逐条查询，
   而是把两个**已排序的 k-mer 位置流做一次归并**（FastGA.c 的 PAIR 文件流），相同
   40-mer 的两侧位置对 (i,j) 在一次 O(|A|+|B|) 线性扫描中全部发现。逐查询是
   O(|A|·log|B|)，归并消除了常数级和随机 IO 开销——这是 2 Gbp vs 5 分钟的主要来源。
2. **lcp 连续传播 → 种子长度自适应**：排序流中相邻相同 k-mer 的 lcp（最长公共前缀）
   直接给出共享长度，40-mer 命中自然扩展为任意长度的最长共享字符串（adaptamer），
   无需对不同 k 反复查询，且天然支持"频率过滤后最长种子"的语义。
3. **anti-diagonal 坐标（在种子流中计算，不在 GIX 里）**：`.ktab` 条目只存
   contig + 位置；两索引归并产种子命中 (i,j) 时**即时**算 `diag = i−j`、
   `anti = i+j` 写入种子流条目，链扫描（`align_contigs` 的 tube 逻辑）直接使用。
   "预编码"是归并阶段的 O(1) 实现选择，不是 GIX 的数据内容。
4. **流式 + 定宽条目**：`Post_List` 按块流式读入（POST_BLOCK），固定宽度条目
   （swide）顺序扫描，内存驻留可控。
5. **桶排序近线性构建**：MSD radix 风格（首字节 1024 桶 + Ksplit 负载均衡），
   多线程并行，构建接近线性；GIX 持久化后多轮比对复用（-k）。

### 10.3 代价与限制

- **空间巨大**：14 GB/Gbp（40-mer 表 + 位置 + 桶），细菌规模（5 Mb）约 70 MB 可接受，
  但大规模集合不可行。
- k=40 固定（可调但需 ≥12 且被 4 整除）；依赖完全匹配种子，靠 lcp 扩展与频率过滤
  （-f 10）控制重复区域。
- 私有格式，构建/读取绑定 FastGA 生态。

### 10.4 GIX 的设计价值与 pgr 的借鉴

GIX 的三个核心设计——(1) 两排序 k-mer 流线性归并找同源（O(|A|+|B|)，免逐查询）；
(2) 相邻条目 lcp 连续传播 → 种子长度自适应；(3) 归并时即时算 diag/anti 坐标 +
MSD 桶排序 + 定宽流式扫描——构成其高效种子检测的全部要点。代价是 14 GB/Gbp 的
空间与私有格式。

pgr 的 `.pgi` 索引与 `pgr align pgi` 比对管线以 **syncmer 稀疏 + 两流归并**的形式
落地了这三点（细菌级几十 MB，远小于 GIX），详见 [[../design/pgi-align.md]] 与
[[../benchmarks/bench-pgi-vs-gixmake.md]]；完整 GIX（.ktab 分片 + Kmer_Stream
流式）未移植。

## 11. 从 GIX 到 Wave align 的完整算法管线

> 本节把 FastGA 从索引到比对的完整数据流串起来，标注源码位置与关键参数，
> 便于整体理解或移植其中的算法模式。

```
FASTA
  │ 1. GDB 构建（FAtoGDB / GDB.c）
  ▼
GDB（.1gdb 元数据 + .bps 2-bit 序列，scaffold→contig 两级）
  │ 2. GIX 构建（GIXmake.c）
  ▼
GIX（.gix 代理 + N×.ktab.<int> 分片 = (12,8) syncmer 起始的 40-mer 排序表）
  │ 3. 归并找种子（FastGA.c new_merge_thread）
  ▼
PAIR 流（种子位置对：ipost/icont/jpost/jcont/lcp，按前缀面板归并）
  │ 4. 链扫描（FastGA.c align_contigs）
  ▼
"tube" 命中（对角线带 × anti-diagonal 区间，含链覆盖判定）
  │ 5. Wave 局部对齐（align.c Local_Alignment / forward_wave / reverse_wave）
  ▼
比对路径（start/end + diff + trace cells）
  │ 6. Trace 回溯（align.c Compute_Alignment / dandc_nd / trace_nd）
  ▼
trace points → .1aln（ONEcode，按 contig1→contig2→start 排序）
  │ 7. 格式输出（ALNtoPAF.c / ALNtoPSL.c，线性展开）
  ▼
PAF / PSL（含 CIGAR）
```

### 11.1 步骤 1-2：GDB + GIX（离线，可持久化复用）

- **GDB**（GDB.c）：FASTA → 2-bit 压缩（`COMPRESSED_LEN = ceil(len/4)`），
  scaffold/contig 两级，N 按 `-n` 阈值拆 gap；元数据（.1gdb）与序列（.bps）分离，
  序列默认 EXTERNAL 文件态（fseeko 随机访问，见 §9.4）。
- **GIX**（GIXmake.c）：只索引"以 (12,8) canonical syncmer 起始"的 k=40 k-mer
  （`TMER=12, SMER=8, SOFF=4`）；实际检测是分布线程内联扫描，对 12-mer 窗口
  内 5 个重叠 8-mer 取 canonical 最小值、锚在窗口起点（`j = i-SOFF`，
  GIXmake.c:569；独立的 `is_syncmer` 为 `DEBUG_SYNCMERS` 死代码，见 §3.2）。
  排序：首字节 1024 桶 → `Ksplit` 均衡分片 → 多线程桶排序，
  输出 `-T` 个 `.ktab.<int>` 分片（Kmer_Stream 流式读取，libfastk.c）。
- 参数：`-k 40`（k 大小）、`-f 10`（种子频率阈值，GIX 构建时同时算频率）。

### 11.2 步骤 3：归并找种子（两个排序 k-mer 流的 join）

`new_merge_thread`（FastGA.c:610，逐前缀面板并行）：

1. 两个 GIX 各自以 `Kmer_Stream` 流式迭代（字典序）。`adaptamer_merge`
   （FastGA.c:2280）把 16-bit k-mer 前缀空间切成 NTHREADS 段**连续面板区间**：
   按较长的 T1/T2 条目数均分取分割点 `parm[t].pbeg = (tp->cpre >> 8)`
   （FastGA.c:2308-2319），每线程处理一段 `[pbeg<<8, pend<<8)` 前缀面板，
   按 `Kmer_Stream.index` 定位起点。前缀面板分区 → 线程间零共享、天然无锁。
2. 对 T1 的每个前缀 `cpre`：T2 跳过前缀更小的条目，把前缀 == cpre 的 T2 条目
   载入小缓存（`cache`），然后做**前缀面板内的归并**。
3. 面板内相同的 40-mer（T1 条目 vs T2 缓存）产出种子位置对，写 PAIR 流；
   **频率过滤**：共享该 40-mer 的记录数 ≥ `FREQ`（-f 10，指针比较
   `hgh >= low + FREQ×kbyte`）则跳过。
4. 种子条目 = (ipost, icont, jpost, jcont, lcp)：两侧位置 + contig + 与下一个
   相同 k-mer 的**最长公共前缀**（lcp 连续传播 → adaptamer 可扩展到任意长度）。

复杂度：O(|A|+|B|) 线性归并（每个 k-mer 恰好处理一次），而非逐查询。

### 11.3 步骤 4：链扫描（anti-diagonal 空间）

`align_contigs`（FastGA.c:2973）把 PAIR 流转成 anti-diagonal 坐标
（`diag = i-j`、`anti = i+j`），按对角线桶（`diag >> BUCK_SHIFT`，桶宽 64）
组织后扫描：

1. 按对角线分三段（b/m/e）归并相邻对角线，保证链不跨对角线带。
2. 维护链的 `alow..ahgh`（anti 区间）与 `dgmin..dgmax`（对角线区间）"tube"；
   种子间距超过 `CHAIN_BREAK`（2000，=2×-s 1000）时结束当前链。
3. 链的 anti 覆盖 ≥ `CHAIN_MIN`（170，=2×-c 85）则触发"tube"处理
   （否则丢弃）。self 比对时跳过完全相同的对角线。
4. tube 处理：加载两个 contig 序列，按 `amid = alow + BUCK_ANTI` 分块调用
  `Local_Alignment`，每次对齐一个 anti-diagonal 子区间；tube 内阈值放宽为
  `alnMin = ALIGN_MIN − 50`、`alnRate = ALIGN_RATE + 0.05`（`align_contigs`），
  即默认 50 bp / 35% 差异容忍。

### 11.4 步骤 5：Wave 局部对齐（Myers wavefront）

`Local_Alignment`（align.c:1423）：

1. 分配 wave 数组（V/M/HA/NA/T，5×vlen）与 trace cells 空间。
2. `forward_wave` 从 mid-line 正向扩展；`reverse_wave` 从低端反向扩展；
   自比对时用 `minp=1/maxp=-1` 防止与自身完全重合。
3. `fshort`/`rshort`：若正向或反向扩展太短（< `DUB_TRIM`），调整边界后
   只重跑短的一侧。

`forward_wave`（align.c:336）核心：

- **0-wave 初始化**：对每个对角线 k，从 `x=(mida+k)/2` 做 snake（同向延伸
  匹配），`V[k]`=最远点、`T[k]`=轨迹位、`M[k]`=匹配数；每 `tspace` 个匹配
  生成一个 `Pebble` cell（ptr/diag/diff/mark）。
- **wave 推进**（`while more && lasta >= besta - TRIM_MLAG`）：
  - 每轮 `dif += 1`，对角线带 ±1 扩展；
  - 每个 k 按 Myers 三分支更新波前（`ac`/`am`/`ap` 三条候选取 max，
    mismatch +1、双 gap +2），再 snake 延伸并更新 M/T；
  - `TRIM_MLAG` 提前终止：最优波前推进超过一定滞后即停止。
- 阈值：长度 ≥ `ALIGN_MIN`（-l 100）、相似度 ≥ `1-ALIGN_RATE`（-i 70%）。

### 11.5 步骤 6：Trace 回溯与编码

`Compute_Alignment`（align.c:5426）根据任务类型组装：

- **DIFF_ONLY**：`split_nd` 只算差异数（用于种子阶段评估）。
- **PLUS_ALIGN**：`dandc_nd`（Hirschberg 风格分治）——`split_nd` 找中间点，
  递归左右，D==1 时输出单个 I/D/S 操作，得到完整路径。
- **DIFF_TRACE / PLUS_TRACE**：`trace_nd` 在中间点按 `tspace` 采样，把比对
  压缩为 trace points（每点记录到下一 trace point 的 diff 与坐标增量）。
- `Gap_Improver`（align.c:6714）对 gap 区域二次精修。

结果按 contig1 → contig2 → start 排序写入 `.1aln`（ONEcode 编码，
alncode.c）；ALNtoPAF/ALNtoPSL 多线程线性展开 trace → CIGAR（`-pafx` 的
`=`/`X` 或 `-pafm` 的 `M`）。

### 11.6 关键设计点总结

| 阶段 | 核心技巧 | 复杂度 |
|------|----------|--------|
| GIX | (12,8) syncmer 稀疏 + 首字节桶排序 | 近线性构建 |
| 归并 | 两排序流前缀面板 join + lcp 传播 | O(|A|+|B|) |
| 链 | anti-diagonal 坐标 + tube 扫描 | 线性（链稀疏）|
| Wave | Myers wavefront（V/M/T + Pebble cells）| 与差异数成正比，优于 O(nm) |
| Trace | Hirschberg 分治 + tspace 采样 | 线性于路径长 |

## 12. 论文实验：灵敏度评估与 FastGA-gapfill（2026-08-06 通读 PDF 补充）

> 来源：《FastGA: fast genome alignment》，Bioinformatics Advances 5(1):vbaf238
> （2025），PDF 全文通读。§5 Experimental results 是唯一讲"灵敏度"的实验章节；
> 其中 §5.2 的 FastGA-gapfill 就是 pgr 混合方案与 ALNfill 的论文原型。

### 12.1 模拟数据灵敏度（§5.1，真值已知）

- **数据构造**：模拟 A、B 两个基因组（各 84 Mb），由 10 kb 块组成；每个块 =
  前端"相似区"（长度 100/200/500/1000/2000/5000 bp × 分歧度 1%–65%）+ 随机序列；
  块顺序随机打乱，保证不存在跨块的长程比对。每个（长度, 分歧度）组合 100 个重复。
- **分歧引入**：在 B 上随机引入 80% 单碱基替换 + 10% 插入 + 10% 缺失。
- **灵敏度定义（核心）**：对每个（长度, 分歧度）组合统计"完整恢复的目标区数量"——
  某目标区被比对覆盖 **≥95%**（A、B 两侧基因组都要）才计为恢复；Fig. 6 的 y 轴就是
  100 个目标区里恢复几个。**不是按碱基覆盖率算，是按"目标区是否完整找回"算。**
- **特异性定义**：false aligned bases = 落在模拟目标区之外的 aligned bases（在 A 上
  统计）；FP 比对 = 跨多个目标区、或任一基因组上 >95% 的比对碱基在目标区外。
  实测 wfmash 74.26% 的比对碱基是假阳性，其余工具 <0.06%。
- **结论**：灵敏度随目标区变长而升、随分歧度升高而降；FastGA/minimap2 的衰减起点约
  为 1%/10%/15%/20%/25%/30% 分歧度（对应 100/200/500/1000/2000/5000 bp）；FastGA
  小目标区略逊 minimap2、大目标区更优；NUCmer 超过 200 bp 后落后；LastZ 全面最高，
  且是唯一在 40% 分歧度（2000/5000 bp 区）仍有合理结果的工具。

### 12.2 真实数据（§5.2，无真值 → 覆盖率代理指标）

- 五个哺乳动物基因组（人类 GRCh38、黑猩猩、长臂猿、猪、小鼠）对齐 CHM13；没有真值，
  用"被比对覆盖的基因组碱基数"作为间接灵敏度指标（只按比对 start/end 计覆盖，
  比对内部的 gap 也算覆盖——对 minimap2/wfmash 这类大 gap chaining 有利）。
- **覆盖率必须结合特异性看**：CHM13 chr16 52.3–96.3 Mb vs mouse chr8 的对照实验——
  直接比对时 wfmash 覆盖 28.4 Mb、LastZ 14.30 Mb；打乱中间 mouse 序列后 wfmash
  仍 bridging 出 22 Mb（多为假阳性），LastZ 掉到真实可比对区 3.69 Mb；反向互补中间
  序列后 LastZ 恢复 11.42 Mb。说明 wfmash 的覆盖虚高来自假阳性 bridging。

### 12.3 FastGA-gapfill：论文原版混合方案（§5.2）

- **流程**：FastGA 比对为锚点 → 对每对"**顺序一致、方向一致、不重叠**"、间隔
  ≤1 Mb（默认）的锚点，用两侧锚点在两个基因组上的 end/start 定义 bounding box →
  box 与锚点允许 **1 kb（默认）重叠**以利 LastZ 播种 → 只保留最小 box（没有更小的
  包含于其中）→ 每个 box 跑 LastZ → 合并 FastGA + LastZ 输出为最终结果。
- **结论**：FastGA-gapfill 灵敏度接近 LastZ，速度比 LastZ 快 19.3×–137.5×。
- **与 ALNfill 参数吻合**：c-zhou/alnfill 的 `alngap` 默认 `-l 100 / -m 1M / -e 1K`、
  reciprocal best（`-f 0.5`）就是论文 gapfill 的工程化实现；ALNfill README 也承认
  `-e` 造成的 FastGA/LastZ 重叠是已知问题（`-e 0` 又会漏掉从锚点延伸出去的比对）。
- **"顺序一致、方向一致"即共线性前提**——证实 [[../design/pgi-lastz-hybrid.md]]
  §3.4 的结论：hybrid 模式只适合 syntenic 搜索。
- **pgr 的原生实现**：`pgr align fill`（`src/cmd_pgr/align/fill.rs`）把论文 gapfill
  工程化为：按 (target, query, strand) 分组、按 target start 排序，对每对**相邻、不重叠、
  同链同向**的锚点，若 **target 与 query 两侧 gap 都在 [min_gap, max_gap]**（fill.rs:277-
  282）则生成 box，box 按 `--overlap` 在两侧外扩作 LASTZ 播种缓冲（fill.rs:285-289，即
  论文/ALNfill 的 1 kb seeding buffer 语义），再 rayon 并行逐 box 跑 LASTZ。与论文的
  一个差别：pgr 不做"只保留最小 box（无更小 box 包含）"这一步（论文 §5.2 用它削减重叠
  box），而是把重叠/冗余交给下游 `pgr pl chainnet` 统一去重。
- **`pgr align rest`**（rest.rs）是互补侧填充：对锚点未覆盖的"整基因组洞"做 target/query
  两侧 1D runlist 补集，再用 **syncmer/minimizer 预过滤**（`sample_hole`，rest.rs:406）配对
  可能的 hole 对——与 FastGA 用 (12,8) syncmer 稀疏采样的思路同族，只是这里用于
  "洞配对预筛选"而非种子检测。

### 12.4 其他物种与实验口径（§5.3 + §4 工程细节）

- **DToL 12 物种**（昆虫/鱼/鸟/爬行/哺乳/两栖各一对，基因组几百 Mb 到 24 Gb）：
  每个属做 (i) 种内单倍型 vs (ii) 种间比对两种。FastGA 种内覆盖
  85.6%–99.0%，种间覆盖 25%–99%（蛾子 A. psi vs A. aceris 最远缘，
  覆盖与 human vs mouse 相当）；两栖 L. vulgaris 24.2 Gb 是唯一全部
  aligner 都完成的（4611 CPU 分钟 / 29 GB），LastZ 有 4 个染色体对
  48 h 内未完成。
- **真实数据的实验口径**（做对照实验时值得照抄）：
  1. **soft-mask 只喂 LastZ**（重复区处理），其余 aligner 用未掩码序列；
  2. LastZ/NUCmer **按染色体对逐个跑**（输入长度限制），CIGAR 输出
     （FastGA `-pafx` / minimap2 `-c` / LastZ `--format=PAF:wfmash`）；
  3. 覆盖统计**只按比对 start/end 计**（比对内 gap 也算覆盖）；
  4. FastGA 索引构建（GIX 排序）峰值内存 ~29 GB，**比对主进程仅 ~1 GB**——
     内存大头在索引排序，不在对齐（pgr 的 `.pgi` 构建同理，见
     [[../benchmarks/bench-pgi-vs-gixmake.md]]）。

## 13. 并行化架构与 DALIGNER 遗留机制（2026-08-12 源码复核补充）

> 前文各节散见多线程的描述；本节把 FastGA 从种子归并到 .1aln 的**并行数据流**
> 完整串起来。核心结论：**并行度全部来自"前缀/contig 空间切分 + 每线程私有临时
> 文件"，线程间零共享、无锁**。这一点对 pgr 用 rayon 重写 pgi 管线的并行化取舍
> 很有参考价值。

```
GIX1 / GIX2（排序 k-mer 流）
  │ ① adaptamer_merge（FastGA.c:2280）：前缀面板空间切 NTHREADS 段 → 每线程一段
  ▼
NTHREADS² 个临时 PAIR 文件（正链 N_unit + 负链 C_unit，按源 contig Select[] 分发）
  │ ② reimport_thread（FastGA.c:2639）：每线程读自己的 PAIR 文件，算 diag/anti/band
  ▼
sarray（按源 contig 桶排布的大数组）→ RSDsort（rmsd_sort）每桶排序
  │ ③ search_seeds（FastGA.c:3715）：NTHREADS 线程各处理一段源 contig 区间
  ▼
align_contigs（每对 contig × 对角线段）→ 写本线程 Overlap 到 .las gather 文件
  │ ④ 每线程把 .las 排序 → 块文件（SORT_MAP：aread,abpos,bread,...）
  ▼
NTHREADS 个排序 .las 块
  │ ⑤ la_merge（FastGA.c:3991）：k-way 堆归并
  ▼
.1aln（写 header novl,tspace + 各 Overlap 的 trace）
```

1. **① 归并**：`adaptamer_merge` 把 16-bit k-mer 前缀空间按较长的 T1/T2 条目数均分
   成 NTHREADS 段连续面板区间（`parm[t].pbeg = tp->cpre >> 8`），每线程独立处理一段，
   把种子对写入**自己的一批临时文件**——按源 contig 分发到 NTHREADS² 个
   `PAIR` 文件（正链 `N_unit`、负链 `C_unit`，`Select[acont]` 定位目标），
   缓冲满才 `write()`（IOBuffer，~1 MB/块）。**写私有文件而非共享队列**避免了
   跨线程加锁，代价是后续要按文件再聚合。
2. **② 重导入 + 排序**：`reimport_thread` 每线程读回自己的 PAIR 文件，算
   `diag = i−j`、`anti = i+j`、`band = diag>>BUCK_SHIFT`（FastGA.c:2704-2716），
   按源 contig 桶填入 sarray；随后对每桶做 `rmsd_sort`（RSDsort.c）。这一步
   把"种子位置对"物化为 `align_contigs` 要用的 anti-diagonal 坐标流。
3. **③ 对齐**：`search_seeds` 每线程领一段**源 contig 区间**（`range->beg..end`），
   逐对 contig、逐对角线段调 `align_contigs`；每个 Overlap 记录（OVL_SIZE + trace）
   `fwrite` 到本线程的 gather 文件 `SORT_PATH/<root>_algn.<pid>.<tid>.las`
   （FastGA.c:3269）——**这里沿用了 DALIGNER 的 `.las` Overlap 格式**。
4. **④ 块排序**：每线程把 gather `.las` 按 `SORT_MAP`（aread → abpos → bread →
   COMP → ...，FastGA.c:3799）排序为"块文件" `ALGN_UNIQ.<tid>.las`。
5. **⑤ LAmerge**：`la_merge`（FastGA.c:3991）对 NTHREADS 个排序块做 **k-way 堆归并**
   （堆顶 Overlap + 每源按块 `ovl_reload` 流式补读），边归并边写 `.1aln`（先
   `open_Aln_Write` 写 header：版本/provenance/tspace/GDB 引用，再逐条写 Overlap 的
   trace）。**这一步是 DALIGNER LAmerge 的移植**——FastGA 内部把"并行对齐 → 汇聚"
   完整复用了 DALIGNER 的 `.las`+LAmerge 机件，最终才转成 ONEcode `.1aln`。

**对 pgr 的启示**：

- pgr 的 `pgr align pgi`（[[../design/pgi-align.md]]）与 `pgr sd` 的并行化可以对照
  这套"前缀/contig 空间切分 + 每线程私有输出 + 最后归并"模式：**不要在共享结构上加锁，
  而是按 key 空间预切分，让每线程写自己的分片，末尾再归并**。pgr 的 rayon
  `par_iter` + 每线程 `Vec` 分片收集正是等价形态。
- FastGA 复用 DALIGNER `.las`/LAmerge 说明：成熟的并行对齐工具往往内置一套
  "gather → 排序 → k-way 归并"的磁盘中间态，宁可多一次 IO 也要避免共享内存竞争。
  pgr 若做大规模并行比对输出，值得参考"每线程临时分片 + 排序 + 归并"而非共享写入。
- `.las` 中间态是 FastGA 的**内部格式**（非用户可见 API），pgr 无需复刻；其可借鉴点
  是"分片-归并"的数据流结构，不是 `.las` 编码本身。
