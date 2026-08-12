# ALNfill 源码分析（FastGA + LastZ 混合补 gap）

> 整理于 2026-08-06，通读 `alnfill-main/`（Chenxi Zhou，MIT，2024-12-19）全部源码。
> 定位：论文 FastGA-gapfill（[[fastga.md]] §12.3）的工程化实现，与
> [[../design/pgi-lastz-hybrid.md]] 的 pgi+LASTZ 混合方案直接对标。
> 目录自带 `.gitignore`（`*`），按 AGENTS.md 为纯参考、不入库。

## 1. 概览

- 两个程序：`alngap`（从 PAF 找 gap 区间）+ `alnfill`（对每个区间跑 LastZ 并
  坐标回移），配套 `k*.c/h`（lh3 klib）与 `rtree.c/h`（二维包围盒去嵌套）。
- 用法（README）：`FastGA qry.fa ref.fa > fga.paf` →
  `alngap -t8 fga.paf > intervals.txt` → `alnfill -t8 ref.fa qry.fa intervals.txt
  > laz.paf` → `cat fga.paf laz.paf > all.paf`。
- 输入 PAF 的 qname=qry、tname=ref；输出 PAF 保持同一约定，可直接 cat 合并。
- 只依赖 zlib + 系统 `lastz` 可执行（`-z` 可指路径）。

## 2. alngap：找 gap（alngap.c，634 行）

### 2.1 读取与去冗余（默认开，`-a` 关闭，`-f 0.5`）

- `read_pafs`：读一个或多个 PAF(.gz)，只取 12 列里的
  qn/ql/qs/qe/tn/tl/ts/te/ml（**不读 strand 列**，见 paf.c:35-65）。
- `reciprocal_best_aligns`：按 mlen **升序**排序，贪心保留——每条比对若与
  "已保留区间集"在 query 侧总重叠 ≤ mlen×0.5 且 target 侧 ≤ mlen×0.5 则保留，
  并把自身区间并入两侧的区间集合（rangeset_overlap/add，排序数组 + 二分）。
  效果是去冗余锚点（短比对优先，长比对若被已保留比对覆盖 >50% 任一侧则丢弃），
  并非严格意义的 reciprocal best hit。

### 2.2 找 gap（align_gaps）

- 按 RORDER（`alngap.c:263-297`，键序为 (qid, tid, qbeg, tbeg, qend, tend)，
  注意是 qbeg/tbeg 先、qend/tend 后）排序后按 (q,t) 序列对分组（`align_gaps` 的
  `ranges` 数组，`alngap.c:442-451` 记录每组起始+长度）；并行以**序列对分组**为
  单元（`kt_for` 调度 `gap_core`，见 §6.6）。
- 每组**首尾各加一个哨兵**（0 位与全长位，`alngap.c:349-354`）→ 染色体首端/尾端
  的 gap 也会被 box 化；端部哨兵的 box 上界/下界被钳到 0 / 全长，故不会越界。
- 对每个锚点 aln1，向后找 query 起点满足
  `qbeg2 ∈ [qend1+min_gap, qend1+max_gap)` 的锚点 aln2（query 侧 gap 在
  [min_gap, max_gap)）；再算 target 侧两区间间距
  `dist = max(tbeg1,tbeg2) − min(tend1,tend2)`，要求 `dist ∈ [min_gap, max_gap]`。
  **不检查链向**：混合方向的锚点对（target 顺序相反）只要间距达标同样成 box。
  - 扫描方式（`alngap.c:364-369`）：对每个 aln1 用 `aln2s`/`aln2e` 两个指针
    （while 前进）圈出 `[qend1+min_gap, qend1+max_gap)` 窗口内的 aln2，再逐个配对；
    因 aln2s 每次从 `aln1+1` 重启而非真正单调滑窗，最坏 O(n²)，但受 max_gap 窗口
    约束、实际开销集中在有 gap 的局部——pgr 若用排序数组 + 二分找窗口边界可做到
    真正线性。
- box 区间（`-e`=max_ovl，默认 1000）：
  - query：上游锚点 query 末端 − max_ovl 到下游锚点 query 起点 + max_ovl
    （即 `[max(abpos1, aepos1−max_ovl), min(aepos2, abpos2+max_ovl)]`，钳到锚点自身边界）；
  - target：同向时从上游锚点内端(bepos1) − max_ovl 到下游锚点内起点(bbpos2) + max_ovl
    （同样钳位）；用 `bepos1<bepos2` / `bbpos1>bbpos2` 条件分支兼容两锚点 target
    相对顺序，**方向相反（反链）时 box 覆盖两锚点之间的反链区域**——即代码并不
    要求两锚点 target 顺序与 query 一致。
  - 同时输出四侧**实际重叠长度**（qbol/qeol/tbol/teol，锚点短于 max_ovl 时钳位），
    供 alnfill 校验区间越界。
- **最小 box 去嵌套**：AORDER 按面积升序排序 → rtree 只插入"内部不含已插入节点"
  的 box（`rtree_exist_node_inside`），面积小的先插，最终只保留不包含任何更小
  box 的 box——与论文"retain only the minimal ones"一致。
- 输出 10 列：`#Q_NAME Q_BEG Q_END T_NAME T_BEG T_END Q_BEG_OVL Q_END_OVL
  T_BEG_OVL T_END_OVL`（头一行 `#Q_NAME...` 也是标准输出的一部分）。

### 2.3 默认参数

| 参数 | 默认 | 含义 |
|------|------|------|
| `-l` | 100 | 最小 gap（小于此不补，含 target 侧 dist）|
| `-m` | 1M | 最大 gap（大于此不补；`-l`/`-m`/`-f` 均走 `parse_num`，十进制 K/M/G）|
| `-e` | 1K | 两侧锚点内收重叠量（播种缓冲）|
| `-f` | 0.5 | 去冗余时单侧允许的最大被覆盖比例（按 mlen 计）|
| `-a` | off | 关闭 reciprocal best 去冗余 |
| `-t` | 1 | 线程数 |

> **CLI quirk**：`-f` 虽默认是浮点 0.5，但命令行传入时经 `parse_num`（`alngap.c:544`，
> `max_cov = (int) parse_num(opt.arg)`，help 也写成 `-f INT`）**转成整数**——传小数
> `-f 0.5` 会被截成 0（= 完全不约束重叠），实际只能传 `-f 0` 或 `-f 1` 这类整数值。
> 另 `-o`（`alngap.c:548-555`）可把输出写到文件（`freopen` 重定向 stdout，`-` 表示 stdout）；`-v` 为 verbose。
>
> **参数解析细节**：`parse_num2`（`alngap.c:494-509`）支持 K/M/G 后缀，但按**十进制**
> 乘（×1e3/1e6/1e9，非 1024），并 `+0.499` 取整；仅 `-l/-m/-f` 走 `parse_num`。
> `-e`（`alngap.c:545`）与 `-t`（`alngap.c:547`）用 `atoi`，**不支持 K/M/G 后缀**——
> 例如 `-e 1K` 会被解析成 1 而非 1000。pgr 若做同类 CLI 需统一后缀语义。

## 3. alnfill：跑 LastZ（alnfill.c，530 行）

- **整库入内存**：`make_sdict_from_fa`（sdict.c:131）用 kopen（bgzip 感知）读入
  全部序列并 `strdup` 保存——2 Gbp 级基因组 RAM 需求很高；且每条序列长度
  **不能 >4 Gb**（`>UINT32_MAX` 直接报错退出，sdict.c:150）。且 **ref 与 qry 各自
  调用一次** `make_sdict_from_fa`（alnfill.c:375-376），两份全基因组同时常驻内存。
  pgr 侧用 2bit/loc 区间提取可避免。
- **启动时自检**：`check_executable`（alnfill.c:83-93）用 `command -v <lastz>` 先校验
  lastz 可执行，缺失即报错退出；`run_system_cmd`（alnfill.c:74-81）调 `system()` 且带
  一次重试。pgr 若外包 lastz 进程可参考"先探测可执行、失败重试"的稳健性。
- 读 interval 文件：跳过 `#` 头行；解析 10 列（≥6 列可用）；用 ovl 列做越界校验
  （`qbeg ≥ qbol`、`qend+qeol ≤ qlen`、target 同理），非法区间跳过；
  **若 qname/tname 在已加载的 FASTA 中找不到则直接报错退出**（非跳过）。
  interval 文件经 `gzopen`（alnfill.c:387）读取，**支持 .gz 压缩**。
- 每线程一套临时文件（mkstemp 生成一个模板名，派生出 `_A.fna`/`_B.fna`/`_O.paf`
  三个具名临时文件，另有一个被 unlink 的匿名输出缓冲文件）；对每个 interval：
  临时文件目录由 **`-w`** 指定（`alnfill.c:327`，默认 `./`）；`-o` 写输出到文件、
  `-z` 指定 lastz 路径（默认 `lastz`）。
  1. 提取 `ref[tbeg,tend)`、`qry[qbeg,qend)` 写临时 FASTA（`>name` + 原序列名）；
  2. 跑 `lastz --format=PAF:wfmash --ambiguous=iupac --output=<pfile> <tfile> <qfile>`
     —— **target 在前**；除格式与 IUPAC 外全是 lastz 默认（无打分预设）；
  3. 读回 PAF 逐行**坐标回移**（paf_parse1）：ql/tl 换全长、qs/qe += qbeg、
     ts/te += tbeg，**第 9 列起（n_match/blk/mapq/CIGAR 等）原样拷贝**（PAF 第 9 列
     是 n_match 而非 CIGAR，CIGAR 在 `cg:` 标签里）；
  4. 删三个具名临时文件。
- 结果收集：各线程的匿名临时文件最后按线程序 cat 到 stdout；README 建议与
  FastGA PAF 直接 cat，**无去重**。

## 4. 与 pgr 混合方案的对照（详见 [[../design/pgi-lastz-hybrid.md]] §3.6）

1. **gap 大小过滤**：alngap 只补 [100, 1M] 的 gap（query、target 双侧都要达标）。
   pgr 的 holes 方案应加同样的区间长度过滤——太小的 gap 是普通 indel/碎片、
   太大的是真分歧/新序列（我们笔记 §4 已提"超长无锚点区间应跳过"，正好对应
   max_gap）。
2. **首尾 gap**：哨兵让染色体首端→第一个锚点、最后一个锚点→尾端的 gap 也会被补
   （同样限 [100, 1M]）；pgr holes 天然含首尾，但要按长度过滤而非全补。
3. **重叠缓冲**：`-e 1K` 是双侧（query+target）、两侧基因组各 1 kb；pgr §3.3 的
   trim 25–50 bp 只收 target 单侧且小两个数量级——实测时值得对比 50/500/1000。
4. **去冗余**：RBA 过滤对重复区的多映射锚点很有用；pgr 若用 pgi PSL 直接算 holes，
   重复区可能产生多余 box，syntenic 场景可先按链上锚点/最长块过滤。
5. **方向**：alngap 不读 strand、混合方向也成 box；论文描述是"一致顺序方向"。
   pgr 的 hybrid 定稿为仅共线性，可加方向过滤（或交给 `chainnet --syn`）。
6. **LastZ 选项**：alnfill 只有 `--format=PAF:wfmash --ambiguous=iupac`（默认打分）；
   pgr 复用 `pgr align lastz --preset`（set01..set07），**预设由用户选择**——泛基因组
   主场景差异小，默认贴近的 set01/set02，远缘比较才选 set06/set07，不默认远缘
   （2026-08-06 与用户讨论后调整，见 [[../design/pgi-lastz-hybrid.md]] §3.2/§3.6）。
7. **坐标回移**：alnfill 在进程内完成区间坐标 → 全长坐标（含 ql/tl 换全长）；
   pgr 已实现 `pgr align fill` / `pgr align rest`（2026-08-06 命名定稿，不做 `pgr pl align`），
   `fa range` 提取时须记录 offset，
   lastz 输出回移后再并 PSL（现有 `pgr align lastz` 只出 LAV、单序列限制）。
8. **合并**：ALNfill 直接 cat 不查重（README 承认 `-e` 造成的重叠是已知问题）；
   pgr 保留 >50% 重叠保长去重兜底。

## 5. 已知问题

- README：`-e` 使 FastGA/LastZ 结果重叠；`-e 0` 会漏掉从锚点延伸出的比对
  （与我们 §3.3"不缩 pgi PSL、只扩补集"的取舍一致）。
- `reciprocal_best_aligns` 的 mlen 升序 + 0.5 阈值是启发式，短比对优先保留；
  对长锚点链场景效果需实测。
- lastz 默认打分无预设，远缘/低相似 gap 的灵敏度受限。
- 全基因组读入内存，超大基因组（>24 Gb，论文 newt 场景）不现实。

## 6. 可迁移到 pgr 的通用算法（不限于补 gap）

以下技术散落在 alnfill 各模块，均已对照源码核实，pgr 的 chain/net/align 生态可
按需吸收（具体取舍见 §4）：

1. **rangeset 区间合并 + 覆盖统计**（alngap.c:104-182）：排序数组 + 二分
   （`find_last_small`/`find_first_large`）求"与已保留区间集的总重叠量"，再原位
   合并新区间（`rangeset_add`）。这是高效的区间覆盖查询，可用于 pgr 的锚点
   去冗余/重复区覆盖统计，语义等价于"单侧被覆盖 ≤ mlen×f"。
2. **rtree 最小 box 选择**（alngap.c:299-406 + rtree.c）：按面积升序插入二维
   包围盒，`rtree_exist_node_inside` 剔除"包含更小 box"的 box，只保留最小者。
   这正是 chain 里"只补不可被更小锚点覆盖的单 gap"的几何版，pgr chain 的 gap
   筛选可复用同一思路（用现有 coitrees 区间树即可，无需引入 rtree）。
3. **哨兵 + 相邻锚点对**（alngap.c:349-354）：每组 (q,t) 序列对首尾各加一个
   零长哨兵，使染色体首端/尾端的 gap 也被枚举；相邻配对用双指针滑窗而非二重
   全遍历。pgr 的 holes/chain 枚举 gap 时可照此处理端部。
4. **进程内逐字段坐标回移**（alnfill.c:97-134）：不重写整行，只对 PAF 前 8 列
   定点替换（ql/tl 换全长、qs/qe/ts/te 加偏移），第 9 列起原样透传。pgr 在
   `pgr align fill` 中回移 lastz 输出时可复用该"列级重写"而非逐列重建 PSL。
5. **多线程临时文件编排**（alnfill.c:431-496）：每线程一套独立临时文件，结果
   写入各自匿名缓冲，最后按线程序合并到 stdout——避免跨线程排序/加锁。
   与 pgr 惯用 rayon 并行、单写 outfile 不同，但"按线程隔离中间文件"的思路
   可迁移到并行子进程调用的场景（如并行跑 lastz 批次）。

