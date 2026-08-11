# Celera Assembler（wgs-8.3rc2）：原版 OLC 组装器源码分析

> 2026-08-12 整理，纯源码分析（`wgs-8.3rc2/`，revision 4627，2015-05-24 发布，
> 即 CABOG / wgs-assembler）。与 `references/canu.md` 配套：Canu fork 自 Celera
> **r4587**（更早），8.3rc2 是主干上更晚的 r4627——两条线随后分化。本文档记录
> 原版 OLC（overlap → unitig → consensus → scaffold）结构，并逐组件与 Canu
> 对照；pgr 视角沿用 canu.md §8 的设计意图（**多 k unitig 的 OLC 拼接**）。

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
placeContains(440)          放 contained reads（best overlap 或 all overlaps）
placeZombies(465)           放 zombie（循环 contained / 问题 reads）
mergeSplitJoin(471)         ★ merge=bubble pop、split=repeat/unique junction 打断、
                              join=promiscuous unitig 连接
extendByMates(473, 可选)     用 mate 扩展（-E）
reconstructRepeats(484, 可选) 重建重复（-R）
cleanup(~500)               splitDiscontinuous + placeContains + promoteToSingleton
setParentAndHang / 输出
```

- **bubble/repeat 判定用覆盖度证据**（`AS_BAT_MergeSplitJoin.C:46-48`）：
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

## 7. consensus：AS_CNS（MultiAlign + 列投票）

老版 consensus 与 Canu 的 utgcns **同名不同物**：

1. **MA 构建**（`MultiAlignment_CNS.C`）：unitig 的 reads 按 overlap 的方向/
   hangs 把每条 read 的碱基（**beads**）对齐成 **columns**（`Column`/`Bead`/
   `Fragment`/`MANode` 结构，`SeedMAWithFragment:357`）。
2. **逐列投票**（`BaseCall.C`）：`BaseCallQuality:84`（质量加权投票，
   `cw[5]` consensus weight + QV）+ `BaseCallMajority:41`（简单多数），由
   `GetMANodeConsensus`（`MultiAlignment_CNS.C:393`）调用。
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

## 9. Canu vs Celera 差异总表 + pgr 视角

| 环节 | Celera 8.3rc2（r4627） | Canu 2.3（r4587 fork） | pgr 借鉴 |
|---|---|---|---|
| overlap | OlapFromSeedsOVL（seed k=9 + banded DP） | MHAP 默认 + overlapInCore 降级 | seed→verify 同源，pgr 是精确版 |
| trimming | merTrim + finalTrim | overlapBasedTrimming（同源） | 已由 fq 系列覆盖 |
| unitig | utg/bog（AS_BOG）+ bogart（AS_BAT） | bogart（AS_BAT 血统，去 mate） | repeat split 的覆盖度证据可记；气泡不处理 |
| consensus | AS_CNS 列投票（可选 pbdagcon） | utgcns（template + POA-DAG） | 列投票 = pgr polish 的雏形 |
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

## 10. 关键文件清单（速查）

| 组件 | 文件:行 | 内容 |
|---|---|---|
| overlap | `AS_OVL/OlapFromSeedsOVL.C`（5968 行） | seed→banded DP 扩展 + 454 纠错 |
| overlap | `AS_OVL/OlapFromSeedsOVL.H` | `DEFAULT_KMER_LEN=9`、Error_Bound 语义 |
| trimming | `AS_MER/merTrim.C` / `AS_OBT/finalTrim.C` | k-mer / overlap 修剪 |
| unitig | `AS_BOG/BuildUnitigs.C` | BOG unitigger（utg/bog） |
| unitig | `AS_BAT/bogart.C:416-500` | BOGART 主流程 |
| unitig | `AS_BAT/AS_BAT_MergeSplitJoin.C:46-48` | bubble/repeat 覆盖度阈值 6/15 |
| consensus | `AS_CNS/BaseCall.C:41/84` | 多数/质量加权列投票 |
| consensus | `AS_CNS/MultiAlignment_CNS.C:357/393` | MA 播种 / consensus 调用 |
| scaffold | `AS_CGW/AS_CGW_main.C` + `ScaffoldGraph_CGW.C` | mate 驱动 scaffold |
| 流水线 | `AS_RUN/runCA.pl:688/820/1387` | unitigger/consensus 选择 |
