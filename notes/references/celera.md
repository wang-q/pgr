# Celera Assembler（wgs-8.3rc2）：原版 OLC 组装器源码分析

> 2026-08-12 整理，纯源码分析（`wgs-8.3rc2/`，revision 4627，2015-05-24 发布，
> 即 CABOG / wgs-assembler）。与 `references/canu.md` 配套：Canu fork 自 Celera
> **r4587**（更早），8.3rc2 是主干上更晚的 r4627——两条线随后分化。本文档记录
> 原版 OLC（overlap → unitig → consensus → scaffold）结构，并逐组件与 Canu
> 对照；pgr 视角沿用 canu.md §8 的设计意图（**多 k unitig 的 OLC 拼接**）。
> **实现状态（2026-08-12）**：已落地为 `pgr asm ovlp`/`layout`/`cns`/`olc`
>  四命令（`design/olc.md`），§10 的"pgr 借鉴"列按实现结果更新。

## 1. 概况

- **历史**：Myers 2000（Drosophila 论文，标准引用）→ CABOG（7.0 改架构）→
  Miller 2008（BOGART unitigger）→ Goldberg 2006 / Berlin 2014（454/PacBio
  支持）。8.3rc2 已内置 PacBio 生态组件（README：pbutgcns、pbdagcon、BLASR、
  FALCON 一部分）。
- **与 Canu 的关系**：Canu 从 Celera r4587 fork 并重写了 overlap（MHAP）与
  consensus（utgcns）；8.3rc2 主干保留 OlapFromSeeds + AS_CNS，但 runCA 已把
  pbdagcon/pbutgcns 作为可选 consensus——"模板 + DAG 共识"路线是主干演进方向，
  Canu 只是把它走完（详见 §7）。
- **语言/构建**：C + C++，内嵌 kmer 包 r1994（`kmer/`）、samtools、jellyfish 2、
  pbutgcns/pbdagcon/BLASR/FALCON 片段。
- **流水线**：merTrim（k-mer 修剪）→ overlap（mer overlapper）→ overlap error
  correction → finalTrim → unitig（utg/bog/bogart 三选一）→ consensus
  （cns/pbdagcon/pbutgcns）→ scaffold（CGW）→ 输出。

## 2. 仓库结构（`src/`）

```
wgs-8.3rc2/
├── kmer/                  # kmer 包 r1994（meryl 计数）
└── src/
    ├── AS_OVL/            # ★ overlap：OlapFromSeedsOVL.C（Delcher 2007，seed→扩展）
    │   ├── OlapFromSeedsOVL.C   # 5968 行，overlap 主程序（含 454 homopolymer 纠错）
    │   ├── SharedOVL.C / CorrectOlapsOVL.C / FragCorrectOVL.C  # 纠错/共享工具
    │   └── overlap_partition.C
    ├── AS_MER/            # meryl + merTrim（k-mer 判定修剪）、mercy（k-mer 纠错）
    ├── AS_OBT/            # ★ overlap-based trimming：finalTrim（bestEdge/evidenceBased/largestCovered）
    ├── AS_BOG/            # ★ unitigger（BOG，`unitigger` 程序，BuildUnitigs.C）
    ├── AS_BAT/            # ★ bogart（BOGART，Miller 2008，`bogart` 程序）
    │   ├── bogart.C                       # 主流程（阶段标记见 §6）
    │   ├── AS_BAT_BestOverlapGraph.C      # best edge 图
    │   ├── AS_BAT_PopulateUnitig.C        # greedy 建 unitig
    │   ├── AS_BAT_MergeSplitJoin.C        # ★ bubble pop / repeat split / join
    │   ├── AS_BAT_IntersectBubble.C / AS_BAT_ReconstructRepeats.C / AS_BAT_ExtendByMates.C
    │   ├── AS_BAT_PlaceContains.C / AS_BAT_PlaceZombies.C / AS_BAT_MoveContains.C
    │   ├── classifyMates*.C / AS_BAT_EvaluateMates.C / AS_BAT_InsertSizes.C
    │   └── splitUnitigs.C / markRepeatUnique.C / computeCoverageStat.C
    ├── AS_CGB/            # scaffold 阶段的 bubble popper（celagram 输出）
    ├── AS_CNS/            # ★ consensus（MultiAlign beads/columns + BaseCall 列投票）
    │   ├── utgcns.C                     # ★ 老版 consensus 驱动（与 Canu 的 utgcns 同名不同物！）
    │   ├── MultiAlignment_CNS.C         # MA 构建（Column/Bead/Fragment/MANode）
    │   ├── BaseCall.C                   # ★ BaseCallQuality/BaseCallMajority 列投票
    │   ├── MultiAlign.C / MultiAlignUnitig.C / MultiAlignStore.C / tigStore.C
    │   ├── AbacusRefine.C / MergeMultiAligns.C / ApplyAlignment.C / GetAlignmentTrace.C
    │   └── addReadsToUnitigs.C / RefreshMANode.C
    ├── AS_CGW/            # ★ scaffolder（Contig Graph，mate 链接，42 个 .C）
    │   ├── AS_CGW_main.C / ScaffoldGraph_CGW.C
    │   └── CIScaffoldT_*（Biconnected/Merge/Merge_AlignScaffold/Cleanup）
    ├── AS_ALN/            # 比对器（bruteforcedp / dpaligner / forcns）
    ├── AS_REZ/            # repeat/gap 相关（GapFillREZ、ConsistencyChecksREZ）
    ├── AS_GKP/            # gkpStore（read 数据库）
    ├── AS_OVS/            # ovlStore（overlap 数据库）
    ├── AS_RUN/            # runCA.pl（总流水线）+ runCA-overlapStoreBuild.pl
    ├── AS_TER / AS_LIN / AS_REF / AS_ARD / AS_UID / AS_VWR / AS_MSG / AS_ENV / AS_PER
    │                     # 终端/线性代数/参考评估/数据库/UID/查看器/消息/环境/性能
    └── AS_global.C / AS_global.H
```

## 3. 流水线（`AS_RUN/runCA.pl`）

关键参数族（`runCA.pl`）：

- `overlapper`（`:253`）：obt/ovl 共用，默认 **mer overlapper**——两阶段：
  seed 查找（`merOverlapperSeedBatchSize`）+ 扩展（`merOverlapperExtendBatchSize`），
  直接用 meryl 计数过滤（mer overlapper 读 meryl，`:2567`）。
- `unitigger`（`:688`）：**utg（默认，非 SFF）/ bog（SFF，Best Overlap Graph）/
  bogart（AS_BAT）** 三选一；`getUnitigger()`（`:1387`）按 gkpStore 特征自动选
  utg/bog。
- 误差率参数族：`utgErrorRate`（BOG/UTG）、`utgGraphErrorRate`/`utgMergeErrorRate`
  （bogart）、`cgwErrorRate`（scaffold 合并）、`obtErrorRate`（finalTrim）。
- `consensus`（`:820`）：默认 `cns`（AS_CNS），可选 pbdagcon/pbutgcns。
- `batOptions`（`:707`）：bogart 透传参数。

## 4. overlap：OlapFromSeedsOVL（Delcher，2007）

`AS_OVL/OlapFromSeedsOVL.C`（5968 行）是 8.3rc2 唯一的 overlap 检测器：

1. **seed 查找**（`Get_Seeds_From_Store` / `Read_Seeds`）：k-mer seed，默认
   `DEFAULT_KMER_LEN=9`（`OlapFromSeedsOVL.H`），候选按位置存储。
2. **扩展验证**（`Process_Seed`）：对每个 seed 候选做 **banded edit distance**
   （`Edit_Array`/`Edit_Space`/`Edit_Match_Limit`，`EDIT_DIST_PROB_BOUND=1e-4`），
   `Error_Bound[i] = i * MAXERROR_RATE` 判错；分类 dovetail / contain，
   支持 partial overlaps（G 选项，`Doing_Partial_Overlaps`）。
3. **454 homopolymer 纠错集成**（`Doing_Corrections`）：扩展同时输出
   `Set_New_*_Votes`（homopolymer 投票系列），由后续 `CorrectOlapsOVL` /
   `FragCorrectOVL` 消费。

**注意**：8.3rc2 **没有** `overlapInCore`（Canu 里那个 k-mer 哈希实现是 Canu
自己加的，不在 Celera 主干）。

**与 pgr 的关系**：seed → 扩展 骨架与 `pgr asm map`（seed → 精确验证）同源，
pgr 是精确版（无错误模型、无 454 纠错）。

## 5. trimming：merTrim + finalTrim

- `AS_MER/merTrim.C`：k-mer 计数判定修剪（mer total/threshold 参数族），
  `AS_MER/mercy.C` 做 k-mer 纠错。
- `AS_OBT/finalTrim.C`（+ bestEdge / evidenceBased / largestCovered）：基于
  overlap 证据修剪 read 末端（clear range），8.3 的 `obtErrorRate` 族控制。

## 6. unitig：三种 unitigger

### 6.1 `utg` / `bog`（AS_BOG，`unitigger` 程序）

`AS_BOG/BuildUnitigs.C`：**Best Overlap Graph** unitigger（Sanger/short-read
时代），mate（paired-end）驱动：

- `AS_BOG_BestOverlapGraph` / `AS_BOG_ChunkGraph` / `AS_BOG_PopulateUnitig`：
  与 AS_BAT 同构的 BOG + greedy 建 unitig。
- mate 逻辑：`AS_BOG_EvaluateMates`、`AS_BOG_InsertSizes`、`AS_BOG_MateChecker`、
  `AS_BOG_MateBubble`、`mate-based-splitting/`——用配对关系检测/拆 bubble 与
  嵌合。

### 6.2 `bogart`（AS_BAT，`bogart` 程序，Miller 2008）

`AS_BAT/bogart.C` 主流程（与 Canu 版对照，阶段行号）：

```
BUILDING UNITIGS(416)       按 chunk 长度序 populateUnitig（greedy，best edge 延伸）
BUILDING UNITIGS catching missed(425)  全量补漏
placeContains(434)          放 contained reads（best overlap 或 all overlaps）
placeZombies(460)           放 zombie（循环 contained / 问题 reads）
mergeSplitJoin(469)         ★ merge=bubble pop、split=repeat/unique junction 打断、
                              join=promiscuous unitig 连接
extendByMates(477, 可选)     用 mate 扩展（-E）
reconstructRepeats(486, 可选) 重建重复（-R）
cleanup(499)                splitDiscontinuous(499) + placeContains(501) + promoteToSingleton(512)
setParentAndHang(518) / 输出
```

- **bubble/repeat 判定用覆盖度证据**（`AS_BAT_MergeSplitJoin.C:46-47`）：
  `SPURIOUS_COVERAGE_THRESHOLD=6`（非 unitig reads 覆盖 >6 才算 repeat 区）、
  `ISECT_NEEDED_TO_BREAK=15`（确认 repeat junction 的最少 reads 数）——
  比 Canu 的相似度阈值更"证据驱动"。
- 辅助程序：`splitUnitigs.C`、`markRepeatUnique.C`、`computeCoverageStat.C`、
  `classifyMates*.C`（mate 分类）、`petey.C`。

### 6.3 AS_CGB（scaffold 阶段 bubble popper）

`AS_CGB/AS_CGB_Bubble_Popper.C`：对 scaffold 图做 bubble 检测/弹出
（celagram 输出），与 AS_BAT 的 bubble 处理互补。

### 6.4 与 Canu bogart 对照

Canu 继承 **AS_BAT 血统**（头注 r4587），但：

- **删掉** mate 相关（evaluateMates/extendByMates/reconstructRepeats/
  classifyMates/InsertSizes）——Canu 是单分子长读，无配对。
- **新增**：optimizePositions（全 overlap 精化坐标）、detectSpurs、
  coverage gap 检测、AssemblyGraph + markRepeatReads（repeat breaking 重构）、
  classifyTigsAsUnassembled（覆盖度/span 分类）、findCircularContigs。
- bubble 处理：Canu 用相似度/偏差阈值（`similarityBubble`/`deviationBubble`），
  原版用覆盖度证据——原版更接近"证据驱动"，与用户"不处理气泡"的裁定冲突更小。

### 6.5 mate-pair 处理（`classifyMates` / `InsertSizes` / `EvaluateMates`）

mate（paired-end）是原版 unitig/scaffold 的核心证据（Canu 长读版已全部删除）：

- `classifyMates.C:32-71`（`cmWorker`）：对每个 mate 依次跑 **spur → chimera → BFS
  → DFS → RFS → suspicious** 六层搜索，按 mate 间是否存在满意路径把 read 分为
  可接受/可疑，并对 innie/normal/anti/outtie 四向分别计数（`doSearchBFS/DFS/RFS/
  Suspicious`）。并行用 kmer 包的 `sweatShop` 线程池。
- `InsertSizes`（`AS_BAT_InsertSizes.H:30-52`）：每库用已可靠定位的 mate 距离累计
  `_mean/_stddev`（`valid` 要求样本 ≥10），供下游判定"该库 mate 距离是否合理"。
- 判定阈值（`AS_BAT_Datatypes.H:51-52`）：`BADMATE_INTRA_STDDEV=3`（同一 unitig 内
  mate 距均值超 3σ 判坏）、`BADMATE_INTER_STDDEV=5`（超 unitig 末端 5σ 判坏）。
- `evaluateMates`（`AS_BAT_EvaluateMates.C`）在各阶段后复核 mate 一致性；
  `extendByMates`（`bogart.C:477`）用未配对的单端向外扩 unitig。

## 7. consensus：AS_CNS（MultiAlign + 列投票）

老版 consensus 与 Canu 的 utgcns **同名不同物**：

1. **MA 构建**（`MultiAlignment_CNS.C`）：unitig 的 reads 按 overlap 的方向/
   hangs 把每条 read 的碱基（**beads**）对齐成 **columns**（`Column`/`Bead`/
   `Fragment`/`MANode` 结构，`SeedMAWithFragment:357`）。
   - 底层数据（`MultiAlignment_CNS_private.H`）：`Bead`（`soffset`/`foffset` +
     `prev/next/up/down` 双向十字链，`frag_index`/`column_index`，`:208-217`）、
     `Column`（`call` bead 指针 + `BaseCount count[]` 预聚合每列 4 碱基计数 +
     `depth`，`:254-262`）、`Fragment`（含 `is_contained/container_iid/manode`
     placement 元数据，`:228-244`）、`MANode`（若干列的集合，即一个 unitig 的
     MA 原子，`:269-275`）。`BaseCount` 预聚合计数让逐列投票 O(1) 读取。
2. **逐列投票**（`BaseCall.C`）：`BaseCallQuality:84`（质量加权投票，
   `cw[5]` consensus weight + QV）+ `BaseCallMajority:41`（简单多数），由
   `GetMANodeConsensus`（`MultiAlignment_CNS.C:393`）调用。
   - `BaseCallQuality` 不只是多数：把某列 beads 按 **best allele / 其他等位基因 /
     guides（非 read 的 unitig 序列）** 三组分类（`BaseCall.C:118-144`），支持
     等位基因拆分（`split_alleles`）与 **SNP phasing**（`CNS_OPTIONS_DO_PHASING_
     DEFAULT=1`，`MultiAlignment_CNS.H:38-40`）；majority 版才退化回"计数+平局
     用 QV 和打破"（`BaseCall.C:65-70`）。pgr 若做 polish，多数版已是够用的最简
     模型。
3. **精炼**：`AbacusRefine.C`（abacus）、`MergeMultiAligns.C`、`ApplyAlignment.C`、
   `RefreshMANode.C`。
4. 可选 **pbdagcon / pbutgcns**（runCA `consensus` 参数）——PacBio 的
   模板+DAG 共识作为外部选项已存在。

**与 Canu utgcns 对照**：老版 = 直接列投票（MA 从 overlap 方向构建，无重比对）；
Canu = template stitch → edlib 重比对到模板 → AlnGraphBoost POA-DAG →
bestPath + min-coverage 修剪（为高噪声长读设计）。**pgr 完美匹配场景下，
老版"列投票"反而是更直接的模型**（map 回来的 reads 天然无错位，直接多数投票
即可，不需要 DAG）。

## 8. scaffold：AS_CGW

- `AS_CGW_main.C` + `ScaffoldGraph_CGW.C`：基于 **mate 链接**的 scaffold 图，
  `CIScaffoldT_*` 系列（biconnected 组件、merge、align scaffold、cleanup）。
- 宏基因组参数先例：`cgwMergeMissingThreshold`（`:789`）——合并 scaffold 时
  允许一定比例 missing mates，注释明确提到 metagenomics 菌株保守区场景。

## 9. 数据存储与核心数据结构 + 工程技巧

### 9.1 术语映射：screed / overlap store / gatekeeper / bank

经典 OLC / AMOS 生态的存储层术语，与 8.3rc2 实际源码的对应（AMOS 本身不在本仓库，见 9.6）：

| 经典/AMOS 术语 | 8.3rc2 实现 | 说明 |
|---|---|---|
| screed（read 库旧称） | **gkpStore**（`AS_GKP` + `AS_PER_gkpStore.H`） | 7.0 之后 read DB 由 gatekeeper store 取代；**"screed" 一词在源码中已不存在**（grep 无结果） |
| overlap store | **ovlStore**（`AS_OVS`） | 位打包 overlap 数据库 |
| gatekeeper | **gkpStore / `gatekeeper` 程序**（`AS_GKP_main.C`） | "看门人"：校验 FRG/Library/Link/Placement 消息并写库 |
| bank（AMOS 的 unitig/assembly DB） | **tigStore**（`AS_CNS/tigStore.C`） | unitig/contig/consensus 的 MultiAlign 库 |

### 9.2 gkpStore：read 数据库（`AS_PER/gkFragment.H`）

- 三种 fragment 类型（`gkFragment.H:80-91`）：**Packed**（<256bp，
  `AS_READ_MAX_PACKED_LEN_BITS=8`）、**Normal**（≤2^18-1，
  `AS_READ_MAX_NORMAL_LEN_BITS=18`，`AS_global.H:220-227`）、**Strobe**（预留）。
  记录体用编译期位域打包，packed 型恰为 32bit，配 `#error ... size wrong`
  静态断言（`gkFragment.H:125-127`）。
- 每条 read 固定字段：UID / IID / **mateIID** / **libraryIID** / **orientation** /
  deleted / nonrandom / seqLen / **clearBeg / clearEnd**；Normal 型额外存
  sequence/quality 的字节偏移（`seqOffset/qltOffset`）。
- **clear range 版本化**（`gkFragment.H:45-65`）：`LATEST / CLR / VEC / MAX /
  TAINT / OBTINITIAL / OBTMERGE / OBTCHIMERA / ECR_0..8` 共 17 档。trimming / 纠错
  阶段逐层产出新区间（OBT 三档 + ECR 纠错 9 档），`LATEST` 指向当前生效区间——
  "trim 只改区间、不动序列"正是 OLC 精髓，pgr `fq trim` 可仿此维护多档区间。
- **mate 方向**（`gkFragment.H:29-37`）：`UNKNOWN / INNIE / OUTTIE / NORMAL /
  ANTINORMAL`（I/O/N/A；ANTI 非合法 mate 方向），unitig/scaffold 的 mate 判定均基于它。
- gatekeeper 入库校验阈值（`AS_GKP_include.H:35-37`）：`GATEKEEPER_MAX_ERROR_RATE
  =0.025`、QV 窗口 50bp、阈值 0.03——入库即校验，下游无需反复防御脏数据。

### 9.3 ovlStore：overlap 数据库（`AS_OVS/AS_OVS_overlap.H`）

- overlap 位打包进 64bit 或 3×32bit word（`AS_OVS_NWORDS` 依 read 长度编译期切换）：
  `a_hang / b_hang`（各 `AS_OVS_HNGBITS = AS_READ_MAX_NORMAL_LEN_BITS + 1`）、
  `flipped`、`orig_erate / corr_erate`（各 12bit）、`seed_value`（8bit）、`type`
  （2bit）。注释明示 "DO NOT rearrange"，`type` 字段靠 pad 位强制对齐
  （`AS_OVS_overlap.H:60-74`）。
- **error rate 编码**（`AS_OVS_overlap.H:49-52`）：`encodeQuality = 10000 × Q`，
  12bit 存 4 位有效数字、上限 40%（`AS_OVS_MAX_ERATE`）。pgr 若压缩 erate 可参考
  "定比 + 封顶"。
- 四种 `type`：`OVL / OBT / MER / UNS`（`AS_OVS_overlap.H:54-58`）——同一库按阶段
  存不同形态：overlap / 修剪区间 / k-mer seed / unassembled。
- 库结构（`AS_OVS_overlapStore.H:34-78`）：`OverlapStoreInfo` 头 + 按 iid 的稀疏索引
  `OverlapStoreOffsetRecord`（`a_iid / fileno / offset / numOlaps`），按
  `numOverlapsPerFile` 切成多文件；`AS_OVS_readOverlapsFromStore` 先定位 offset 再顺序读。
- 端判定辅助（`AS_OVS_overlap.H:271-319`）：由 `a_hang / b_hang / flipped` 组合推导
  5'/3' 端与 contain/container——不存冗余标志，纯位运算。

### 9.4 tigStore：unitig/consensus 数据库（`AS_CNS/tigStore.C`）

- 每条 unitig/contig 存为一个 **MultiAlign**（`MultiAlign.C`），dump 支持
  `properties / frags / unitigs / consensus / layout / multialign / matepair /
  sizes / coverage / thinoverlap / fmap` 12 种视图（`tigStore.C:34-45`）。其中
  **layout** 视图即 read 在 contig 上的 placement（方向 + 坐标），是 pgr
  `asm layout` 输出的对应物。
- BAT / CGW / consensus 各阶段读写同一 tigStore，靠 `unitig_status`
  （`AS_UNIQUE / AS_NOTREZ / AS_SEP / AS_UNASSIGNED`，`tigStore.C:107-110`）与
  `suggest_repeat` 标记推进。

### 9.5 BAT 内存 bank：FragmentInfo + OverlapCache

unitigger 把磁盘三库载入内存以提速（pgr `asm layout` 同款思路）：

- `FragmentInfo`（`AS_BAT_Datatypes.H:255-346`）：紧凑数组——每条 read 的
  `_fragLength / _mateIID / _libIID` + 每库 `_mean / _stddev`，
  `memoryUsage()` 自报 3×uint32/read。
- `OverlapCache`（`AS_BAT_OverlapCache.H:119-184`）：overlap 内存堆，`BAToverlapInt`
  位打包 **8B/条**（`AS_OVS_HNGBITS` hang + `AS_BAT_ERRBITS=7~12bit` erate +
  flipped + b_iid，`:51-57`），工作态展开为 32B 的 `BAToverlap`。用 **memory-mapped
  cache 文件**（`AS_BAT/memoryMappedFile.H`）避免重读 ovlStore；用 `_OVSerate→_BATerate`
  与 `_BATerate→error` 两张查找表做精度转换（`:179-180`）；带线程私有缓冲
  `OverlapCacheThreadData`。**位打包 + 查找表 + mmap 缓存**是 pgr 高性能内存表示的
  直接可借鉴样板。

### 9.6 AMOS

- **本仓库不含任何 AMOS 代码**（grep 无结果）。AMOS 是独立外部包，通过 Celera 输出的
  **FRG 消息格式 / bank 格式**与之互操作。Celera 侧的对称物是 gatekeeper 的 dump 家族：
  `dumpGateKeeperAsFRG / AsFasta / AsFastQ / AsNewbler`（`AS_GKP_include.H:99-147`）——
  即"读库 → 常见格式"的统一出口，pgr 可对齐这套命名。
- **mate-pair 处理完全在 Celera 内部**（`AS_BAT`，见 §6.5），不经 AMOS。

### 9.7 pgr 工程借鉴小结

| 源 | 技术 | pgr 落点 |
|---|---|---|
| `gkFragment` | 编译期位域 + `#error` 静态断言 | `asm ovlp/layout` 内存紧凑表示 |
| clear range 版本化 | 多档区间 + LATEST 指针 | `fq trim` 区间管理 |
| `OverlapCache` | 位打包 + erate 查找表 + mmap | `asm layout` 内存 overlap 表示 |
| `ovlStore` | erate 定比编码 + 按 iid 稀疏索引 | overlap 文件/索引布局 |
| `BaseCallMajority` | 列投票计数 | `asm cns` polish 雏形（§7） |

## 10. Canu vs Celera 差异总表 + pgr 视角

| 环节 | Celera 8.3rc2（r4627） | Canu 2.3（r4587 fork） | pgr 借鉴 |
|---|---|---|---|
| overlap | OlapFromSeedsOVL（seed k=9 + banded DP） | MHAP 默认 + overlapInCore 降级 | seed→verify 同源，pgr 是精确版（`asm ovlp`：canonical k-mer 种子 + ± 双向扩展，已实现） |
| trimming | merTrim + finalTrim | overlapBasedTrimming（同源） | 已由 fq 系列覆盖 |
| unitig | utg/bog（AS_BOG）+ bogart（AS_BAT） | bogart（AS_BAT 血统，去 mate） | greedy best-edge + 互惠检查已实现（`asm layout`）；覆盖度证据 repeat split（6/15）留 v1；气泡不处理 |
| consensus | AS_CNS 列投票（可选 pbdagcon） | utgcns（template + POA-DAG） | 精确缝合已实现（`asm cns`，v0 无投票）；列投票 = pgr polish 的雏形，留 v1 |
| scaffold | AS_CGW（mate 驱动） | 无（长读无 scaffold） | pgr 无 paired map 前不适用 |

**pgr 视角**（承接 canu.md §8 的多 k unitig OLC）：

- AS_BAT 是原版里最接近"unitig 层 OLC 拼接"的组件：greedy best-edge 建
  unitig + 覆盖度证据驱动的 repeat split——bubble 部分维持"不处理"裁定，
  repeat split 的"非 unitig reads 覆盖 >6 才算重复区"可作将来 contig 断裂
  的参考阈值。
- AS_CNS 的 BaseCallQuality 列投票是 pgr polish 的最简模型：完美匹配 map
  回来逐列多数投票即可，无需 Canu 的模板重比对。
- 8.3rc2 已把 pbdagcon 列为可选 consensus，佐证"模板+DAG"是 Celera 主干演进
  方向——对 pgr 的意义仍是"只取 consensus 后半段"，不整套搬 OLC。
- **实现映射（2026-08-12）**：`asm ovlp` = seed→verify 精确版；`asm layout`
  = greedy best-edge + 互惠检查（连接端语义见 canu.md §8.5）；`asm cns`
  = 精确缝合（overlap 已对齐坐标，投票无增量，留 v1）。合成数据 30× 下
  contigs 全部为基因组精确子串；6× 低覆盖出现重复区环形错装（top2 repeat
  近似不触发），v1 优先补 AS_BAT 的覆盖度证据（6/15）。

## 11. 关键文件清单（速查）

| 组件 | 文件:行 | 内容 |
|---|---|---|
| overlap | `AS_OVL/OlapFromSeedsOVL.C`（5968 行） | seed→banded DP 扩展 + 454 纠错 |
| overlap | `AS_OVL/OlapFromSeedsOVL.H:74/82` | `DEFAULT_KMER_LEN=9`、`EDIT_DIST_PROB_BOUND=1e-4` |
| trimming | `AS_MER/merTrim.C` / `AS_OBT/finalTrim.C` | k-mer / overlap 修剪 |
| unitig | `AS_BOG/BuildUnitigs.C` | BOG unitigger（utg/bog） |
| unitig | `AS_BAT/bogart.C:416-519` | BOGART 主流程（placeContains:434/placeZombies:462/mergeSplitJoin:471/extendByMates:477/reconstructRepeats:488） |
| unitig | `AS_BAT/AS_BAT_MergeSplitJoin.C:46-47` | bubble/repeat 覆盖度阈值 6/15 |
| unitig | `AS_BAT/classifyMates.C:32-71` | mate 分类（spur/chimera/BFS/DFS/RFS/suspicious） |
| unitig | `AS_BAT/AS_BAT_InsertSizes.H:30-52` | 每库 insert size mean/stddev |
| read 库 | `AS_PER/gkFragment.H` | gkpStore（Packed/Normal/Strobe + clear range 17 档 + mate 方向） |
| overlap 库 | `AS_OVS/AS_OVS_overlap.H` | ovlStore 位打包（OVL/OBT/MER/UNS + erate 12bit） |
| overlap 库 | `AS_OVS/AS_OVS_overlapStore.H:34-78` | ovlStore 头 + 按 iid 稀疏索引 |
| unitig 库 | `AS_CNS/tigStore.C:34-45` | tigStore / MultiAlign 12 种 dump 视图 |
| 内存 bank | `AS_BAT/AS_BAT_OverlapCache.H:119-184` | BAToverlapInt 8B + mmap + erate 查找表 |
| consensus | `AS_CNS/BaseCall.C:41/84` | 多数/质量加权列投票 |
| consensus | `AS_CNS/MultiAlignment_CNS.C:357/393` | MA 播种 / consensus 调用 |
| consensus | `AS_CNS/MultiAlignment_CNS_private.H:208-275` | Bead/Column/Fragment/MANode 结构 |
| scaffold | `AS_CGW/AS_CGW_main.C` + `ScaffoldGraph_CGW.C` | mate 驱动 scaffold |
| 流水线 | `AS_RUN/runCA.pl:688/820/1387` | unitigger/consensus 选择 |
