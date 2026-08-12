# Canu（2.3）：OLC 组装器源码分析（overlap / layout / consensus）

> 2026-08-12 整理，纯源码分析（`canu-2.3/`，版本 2.3，GitHub `marbl/canu`）。
> Canu 是 Celera Assembler r4587（`wgs-assembler`）的 fork，面向高噪声单分子
> 长读（PacBio CLR/ONT）。**与 pgr 的关系**：不是要对 reads 做 OLC，而是用户
> 的设计意图——**把不同 k 各自生成的 unitigs 当"伪 reads"，在 unitig 层做
> OLC 拼接**（见 §8）。本文档记录 Canu 的 OLC 三组件源码结构，以及该设计意图
> 对应的借鉴评估。
> **实现状态（2026-08-12）**：设计意图已落地为 `pgr asm ovlp`/`layout`/`cns`/
> `olc` 四命令（`design/olc.md`），§8.5 回写实现后的理解修正。

## 1. 概况

- **定位**：read 纠错 → 修剪 → OLC 组装（overlap-layout-consensus）的完整
  流水线，输出 contigs.fasta / GFA / BAM。
- **EOL**：Canu 已停止开发（README 原话 "has reached END OF LIFE... use Flye,
  Hifiasm or Verkko"）。它是最好的 OLC 参考实现，但不是未来方向——借鉴只取
  算法思想，不追随其做长读组装。
- **与 Celera 的关系**：源码头注 `This software is based on 'Celera Assembler'
  r4587`。overlap（MHAP）与 consensus（utgcns）是重写的；unitigger（bogart）
  从 Celera 继承骨架并持续维护（详见 §7）。
- **语言/构建**：C++（boilermake），内嵌 `meryl`（k-mer 计数）、`mhap-2.1.3.jar`
  （Java，MinHash overlap）、`utgcns/libpbutgcns`（PacBio pbdagcon 的 POA-DAG）、
  `htslib`（BAM 输出）、`libsnappy`、`libboost`。
- **流水线四步**（`pipelines/canu.pl`，1209 行）：①correction ②trimming
  ③unitigging（= assembly）④outputs；每步各自 meryl 计数 + overlap。

## 2. 仓库结构（OLC 相关）

```
canu-2.3/src/
├── pipelines/canu.pl        # 总流水线（cor/obt/utg 三阶段调度）
├── meryl/                   # k-mer 计数（每阶段先统计 abundance，供 overlap/纠错用）
├── overlapInCore/           # ★ Celera 经典 k-mer seed overlap（已降级）
│   ├── overlapInCore.C              # 入口 + 全局哈希表
│   ├── overlapInCore-Build_Hash_Index.C   # 全 reads k-mer 建哈希索引
│   ├── overlapInCore-Find_Overlaps.C      # ★ 滑窗查哈希 → 候选命中（Find_Overlaps:235）
│   ├── overlapInCore-Process_String_Overlaps.C  # ★ 对角合并 → Myers 扩展（:581/:442）
│   ├── edalign.C                   # Myers 位向量扩展（Extend_Alignment）
│   └── overlapPair.C               # overlap 记录/质量
├── mhap/                    # MHAP 2.1.3 jar + mhapConvert（默认 overlapper）
├── bogart/                  # ★ unitigger（BOG，Celera 血统）
│   ├── bogart.C                     # ★ 主流程（阶段标记见 §5）
│   ├── AS_BAT_BestOverlapGraph.C    # ★ best edge 图（findInitialEdges:67）
│   ├── AS_BAT_PopulateUnitig.C      # ★ greedy 双向延伸（互惠 best edge 种子）
│   ├── AS_BAT_ChunkGraph.C / AS_BAT_Unitig.C / AS_BAT_OptimizePositions.C
│   ├── AS_BAT_PlaceContains.C / AS_BAT_MergeOrphans.C  # contained/orphan/bubble
│   ├── AS_BAT_AssemblyGraph.C + AS_BAT_MarkRepeatReads.C  # ★ repeat breaking
│   ├── AS_BAT_DetectSpurs.C / AS_BAT_SplitDiscontinuous.C / AS_BAT_FindCircular.C
│   └── AS_BAT_PromoteToSingleton.C / AS_BAT_SetParentAndHang.C / AS_BAT_Outputs.C
├── utgcns/                  # ★ consensus（重写）
│   ├── unitigConsensus.C           # ★ template stitch + 重比对 + POA-DAG 投票
│   ├── utgcns.C                    # 分区/参数（默认 pbdagcon + edlib）
│   ├── layoutToPackage.C / unitigPartition.C / stashContains.C
│   └── libpbutgcns/AlnGraphBoost.C # ★ PacBio pbdagcon 的 POA-DAG（bestPath:490）
├── correction/              # read 纠错（Celera 8 没有）
│   ├── generateCorrectionLayouts.C / filterCorrectionLayouts.C
│   ├── falconConsensus.C           # ★ FALCON 移植：列 DP + link 回溯（getConsensus:69）
│   ├── falconsense.C               # 每 read 一个 consensus 的 CLI
│   └── computeGlobalScore.C / errorEstimate.C
├── overlapBasedTrimming/    # 修剪（trimReads-bestEdge / splitReads-*）
├── readErrorDetection/      # 读端错误谱估计（供 overlap 误差率调整）
├── overlapErrorAdjustment/  # 按每 read 的错误谱调整 overlap 误差率（unitig 段）
├── overlapCheck/            # overlap 质量检查
├── overlapAlign/            # overlap 的碱基比对（计算 alignments）
├── stores/                  # seqStore / ovStore / tigStore（二进制存储）
└── utility/                 # 基础设施（bits、strings、system、libbacktrace）
```

## 3. 流水线（`pipelines/canu.pl`）

三段各跑一遍 `meryl 计数 → overlap → 下游`（行号指 `main` 主流程的
调度点，即 `doCorrection/doTrimming/doUnitigging` 三个子程序被调用的位置；
子程序定义分别在 `canu.pl:926/968/1010`）：

```
doCorrection(1144)  meryl → overlap(cor, mhap) → buildCorrectionLayouts →
                    filterCorrectionLayouts → generateCorrectedReads(falconConsensus)
doTrimming(1165)    meryl → overlap(obt, mhap) → trimReads / splitReads
doUnitigging(1179)  meryl → overlap(utg, mhap) → readErrorDetection →
                    overlapErrorAdjustment → unitig(bogart) → consensus(utgcns)
                    → generateOutputs
```

- overlapper 默认全部 `mhap`（`canu.pl:125-127`）；`overlapInCore` 走
  `ovlOverlapper=ovl`，`canu.pl:710-720` 有醒目的 `DO-NOT-USE` /
  `LUDICROUSLY SLOW` 警告框（原笔记写 `:718`，实际是 710–720 的整块警告）。
- correction 是 Canu 相对 Celera 8 的新增阶段：Celera 8 直接用原始 reads
  做 overlap；Canu 先把每条 read 按 overlap 布局纠错成 consensus 再组装。
- 每阶段按 `-canuIterationMax`（默认 1）可迭代重跑：各 `*Check` 子程序以
  `foreach (1..canuIterationMax+1)` 循环（`canu.pl:1146-1147` 等），
  迭代会改写脚本并重新提交。

## 4. overlap 检测

### 4.1 overlapInCore（Celera 经典，k-mer seed → 扩展）

1. **建索引**（`Build_Hash_Index.C`）：所有 reads 的每个 k-mer 以 2-bit 编码滑窗
   哈希，存入哈希表 + check 位向量（`Hash_Check_Array` 快速排除）。
2. **候选命中**（`Find_Overlaps.C:235`）：对每条 read 滑窗，命中即记
   `(read, offset, ref)`；`HOPELESS_MATCH` 过滤两端无意义的命中。
3. **对角合并**（`Process_String_Olaps:581`）：连续 k-mer 命中合并成对角候选，
   用共享 k-mer 数做期望值过滤（`computeExpected(kmerSize, ovlLen, erate):22`）。
4. **扩展验证**（`:442` `Extend_Alignment`，edalign = Myers 位向量）：候选对角
   向两端扩展，错误数超过 `Error_Bound[Olap_Len]` 则丢弃（`:463`），分类
   dovetail / contain，算 quality / evalue。

**与 pgr 的关系**：这是 `pgr asm map` 的"祖先"——seed → verify 骨架相同，
pgr 把扩展步收成"整条精确相等"（无错误模型）。overlapInCore 的 k-mer 计数
过滤、evalue 评分等错误率机制在 pgr 完美匹配路线下不需要。

### 4.2 MHAP（默认）

- 二进制只有 `mhap-2.1.3.jar` + `mhapConvert.C`（MHAP 源码是独立仓库）。
- 算法（Canu 论文 + MHAP 文档）：MinHash sketch + **自适应 k-mer 加权**
  （tf-idf 式，高丰度 k-mer 降权）→ LSH 找候选对 → banded 比对验证错误率。
- 为长读高噪声设计；pgr 有精确 seed 索引（`asm map`、`align pgi`），不适用。

## 5. layout / unitig（bogart）

`bogart.C` 主流程（阶段标记行号）：

```
LOADING AND FILTERING OVERLAPS(491)  OverlapCache 按误差率/长度过滤
BUILDING GREEDY TIGS(540)            BestOverlapGraph → ChunkGraph →
                                     populateUnitig 按 chunk 长度序 greedy 建 tig
                                     optimizePositions(555) → splitDiscontinuous(566)
                                     → detectSpurs(573)
PLACE CONTAINED READS(589)           placeUnplacedUsingAllOverlaps(600) + 二次优化
MERGE ORPHANS(624)                   孤儿并入（无相似度阈值）
MARK SIMPLE BUBBLES(636)             mergeOrphans(deviationBubble, similarityBubble)
GENERATING ASSEMBLY GRAPH(692)       AssemblyGraph：只用与现有 tig 兼容的边
BREAK REPEATS(713)                   markRepeatReads(721)：双定位 reads → 打断
CLEANUP MISTAKES(736)                splitDiscontinuous + promoteToSingleton
GENERATE OUTPUTS(754)                findCircularContigs → setParentAndHang → 写 store
```

关键点：

- **BestOverlapGraph**（`AS_BAT_BestOverlapGraph.C:67` `findInitialEdges`）：
  每条 read 的 5'/3' 端各选一个 best edge（按质量/长度/误差率排序）；若 best
  edge 覆盖的 reads 比例不足 `minReadsBest`（默认 0.8），自动放宽误差率到
  `erateMax` 重算（自适应阈值）。建图前的过滤链（都在 BestOverlapGraph 内）：
  先 `findContains()` 剔除 contained，再按标签过滤 coverage gap / lopsided /
  spur / high-error 四类"问题 reads"，最后才 `findEdges()`。
- **OverlapCache**（`AS_BAT_OverlapCache.C`）：`loadOverlaps:485` 从 ovlStore
  读入后按 `_maxEvalue`（＝误差率）与 `_minOverlap` 过滤（`filterOverlaps:416`），
  再去重（`filterDuplicates:320`）并把 overlap 对称化（`symmetrizeOverlaps:635`）；
  内存由 `-M` 预算约束（`computeOverlapLimit:202`）。→ 对 pgr 的意义：
  OLC 前先按长度/质量过滤 overlap 集、去重、保证对称，能大幅减小 layout 图。
- **greedy 建 tig**（`AS_BAT_PopulateUnitig.C`）：种子 read 要求两端 best edge
  互惠（`edgeTo5 && edgeTo3`，`:154-155`），沿 best edge 单向延伸、已放置即停；
  非互惠的 read 不种子（避免错装）。种子按 **chunk 长度降序**取
  （`bogart.C:545` `nextReadByChunkLength`）——最长 read 先建 tig，与 pgr
  `layout` 的"按 unitig 长度降序取种子"一致。
- **repeat breaking**（`AS_BAT_MarkRepeatReads.C:973` `markRepeatReads`）：
  内部是证据驱动的一条流水线（每 tig）：
  `annotateRepeatsOnRead:83`（从 AssemblyGraph 收集"外部 reads 对 tig 的 overlap"
  → `mergeAnnotations:117`（按证据 read 折叠成 tig 坐标区间）
  → `discardSpannedRepeats`（被 tig 内 read 完整跨越的区间不算重复）
  → `mergeAdjacentRegions`（先外扩 `MIN_ANCHOR_HANG` 再合并相邻区间）
  → `findConfusedEdges`（找"confused edge"：某 read 的 best edge 落在重复区，
  但存在强度相近的 near-best edge 指向 tig 外）
  → `buildBreakPoints` + `splitTigAtReadEnds`（打断）。
  关键常量在文件顶：`MIN_ANCHOR_HANG=500`（重复区边界至少要锚定这么多碱基）、
  `REPEAT_OVERLAP_MIN=50`、`REPEAT_FRACTION=0.5`。**这是"外部 reads 证据 +
  无跨越则打断"的图级实现**——注意 Canu 2.3 并无固定的
  `SPURIOUS_COVERAGE_THRESHOLD`/`ISECT_NEEDED_TO_BREAK` 两个常量（那是 Celera
  8.3 bogart 血统的，见 §8.5），Canu 用 confused-edge + deviationRepeat 判定。
- **气泡**（`mergeOrphans` 第二遍，`bogart.C:639`）：按相似度阈值
  `similarityBubble` 合并平行路径——正是用户裁定"不处理"的那种启发式。

### 5.1 关键参数默认值（`bogart.C:87-128`）

> 注意：`bogart.C` 的 usage help 文本（`:376-423`）与代码实际默认值有**多处不符**
> （help 文本过时），下表以代码默认值为准。

| 参数 | 代码默认 | 含义 |
|---|---|---|
| `erateGraph` / `erateMax` / `erateForced` | 0.075 / 0.100 / 1.0 | best-edge 建图误差率 / 放宽上限 / 强制值（<1.0 才用） |
| `minReadsBest` | 0.8 | 有 best edge 的 reads 比例下限，低于则放宽到 erateMax（help 写 0.9，代码是 0.8） |
| `deviationGraph` / `deviationBubble` / `deviationRepeat` | 6.0 / 6.0 / 3.0 | 三场景下与均值相差的标准差倍数 |
| `confusedAbsolute` / `confusedPercent` | 2500 / 15.0 | repeat 检测里 confused-edge 的绝对碱基差 / 百分比差（help 写 2100/200，代码是 2500/15.0） |
| `minOverlapLen` / `minIntersectLen` / `maxPlacements` | 500 / 500 / 2 | 最小 overlap / 最小交集（建 unitig）/ 最大放置数 |
| `spurDepth` | 3 | spur 检测回溯深度 |
| `lopsidedDiff` | 25.0 | 判定 lopsided read 的 5'/3' 百分比差 |
| `fewReadsNumber` / `tooShortLength` / `spanFraction` / `lowcovFraction` / `lowcovDepth` | 2 / 0 / 1.0 / 0.5 / 3 | 未组装（unassembled）分类的覆盖/长度/跨度阈值（`classifyTigsAsUnassembled:681`） |

→ 对 pgr 的意义：`deviationGraph=6`（best-edge 用均值±6σ 的宽松窗口）、
`confusedAbsolute/Percent`（repeat 区判定）与 `minReadsBest`（自适应放宽）
都是**经验参数**，pgr 的 `asm olc` 目前用固定阈值（top2 边长度比 ≥0.9），
这些是 v1 调参的直接参考来源。

## 6. consensus（utgcns）

`unitigConsensus::generate`（`unitigConsensus.C:1466`）先做
`switchToUncompressedCoordinates:164`（把 homopolymer 压缩坐标解回），再按
consensus 算法字符分发：

- **默认 `P` = pbdagcon**（`generatePBDAG:854`）：
  1. `generateTemplateStitch:259`：按 layout 顺序用 edlib（Myers O(ND)）把 reads
     逐步"缝"成一条模板（带 homopolymer 压缩坐标）。关键细节：对每个待加
     read，只用**期望 overlap 长度**（layout 坐标差）取模板的 80%（`templateSize=0.90`
     实际是 90%）与 read 的 10%（`extensionSize`）做带 band 的 edlib 对齐，
     失败则降 min overlap、升误差率重试（`alignAgain` 标签，`:427` 起）。
  2. 所有 reads 用 edlib 重比对到模板（`alignEdLib:675`：band 递增、错误率
     递增重试；`ERROR_RATE_FACTOR=4`、`NUM_BANDS=2` ⇒ `MAX_RETRIES=8`，
     见 `unitigConsensus.C:35-38`）。
  3. 比对建成 PacBio 的 AlnGraphBoost POA-DAG（`libpbutgcns/AlnGraphBoost.C`），
     `mergeNodes:188` 合并等价节点。
  4. `bestPath:490`（heaviest path）→ `consensusNoSplit:386`：沿 best path 取
     base，低于 `minCoverage` 的段截掉（Canu 的 min-coverage 修剪）。
  → 对 pgr 的意义：Canu 的 template stitch 是"逐步缝合"，每一步只用一个
  局部 banded 对齐来决定接缝位置；pgr 的 `asm cns` 因 overlap 全精确，坐标
  已由 layout 对齐，缝合不需要 edlib——但 Canu 的"取期望 overlap 的局部窗口
  对齐来确定接缝"思想，在 pgr 未来引入错配 overlap 时可直接复用。
- **`Q` = quick**（`generateQuick:976`）：只做 template stitch，不做 consensus——
  拼贴序列直接当 contig（文档说适合做中间检查/抛光输入）。
- **`S` = singleton**（`generateSingleton:1012`）：单 read 原样输出。
- **correction 的 consensus**（`correction/falconConsensus.C:69`）：FALCON 移植，
  证据 reads 对齐到模板后按列 DP + link 回溯打分（`getConsensus`），不是
  AlnGraphBoost。
- **对齐器**：`utgcns.C:316-317` 明说 edlib 是唯一（也是默认）的 aligner——
  `-pbdagcon` 选项只切换 POA 图的构建/遍历算法，碱基对齐一律走 edlib。

## 7. Canu vs Celera（改进清单）

| 环节 | Celera 8 经典 | Canu 2.3 | 结论 |
|---|---|---|---|
| overlap | `overlapInCore`（k-mer seed + Myers） | 默认 MHAP（MinHash + 自适应加权）；overlapInCore 保留但警告不用 | 重写，明显更好 |
| correction | 无独立纠错 | correction 阶段（falconConsensus） | 新增 |
| layout/unitig | BOGART unitigger | bogart（同源，头注 r4587） | 继承 + 维护，骨架未变 |
| consensus | AS_CNS（multi-align + early exit） | utgcns（template stitch + edlib + POA-DAG + bestPath） | 重写，明显更好 |

## 8. pgr 视角：多 k unitig 的 OLC（用户设计意图，2026-08-12）

### 8.1 设计

**不是对 reads 做 OLC**：reads 已经用 DBG（`asm contig`/`unitig`，tadpole/bcalm
血统）压缩成 unitigs。**不同 k（如 21/51/81）各生成一套 unitigs，把
unitigs 当"伪 reads"，直接在 unitig 层做 OLC 拼接**。

### 8.2 为什么合理

- **数据量**：unitigs 数量远小于 reads（宏基因组下可能少 1~2 个数量级），
  OLC 的 all-pairs overlap 成本从不可行变为可行。
- **规避气泡**：unitig 语义 = 最大无分支路径（bcalm graph3 移植，无气泡），
  OLC 拼接只处理 unitig 间 overlap，不引入"平行路径选哪条"的启发式——
  与用户"气泡不如不处理"的裁定一致。
- **多 k 互补**：小 k 连通性好（低覆盖区、重复边界），大 k 特异性强
  （区分重复/菌株）；不同 k 的 unitigs 共享精确子串，天然有 overlap 证据，
  不需要 SPAdes 式多 k 图合并。
- **pgr 基础设施现成**：`asm unitig`（生成，含 `--links`/`--gfa` 方向规则）、
  `asm map`（unitig → 参考的精确匹配，可作 overlap 验证与 consensus 回放）、
  `sam to-rg` + `rg coverage`（覆盖度/投票）。

### 8.3 可借鉴的 Canu 组件

- **consensus 后半段（最值得）**：layout 已知 → 把"伪 reads"（各 k 的 unitigs，
  或原始 reads）放回投票。Canu 在 POA-DAG 上走 bestPath + min-coverage 修剪；
  pgr 是线性完美比对，更简单：`asm map` 多 hit 天然标记重复，单列多数投票
  即 consensus，`--min-coverage` 截断语义与 Canu `consensusNoSplit` 呼应。
- **repeat breaking 思路**（`MarkRepeatReads`）：unitig 双定位（同时匹配两个
  contig 位置）→ 投影 → 无跨越则打断；可作将来 contig 验证的参考。
- **覆盖度/span 分类**（`fewReadsNumber/tooShortLength/spanFraction/lowcovDepth`）
  与 metaMDBG 渐进丰度过滤是同一思想（覆盖度驱动）的两种实现，可并进宏基因
  组低覆盖单元处理方案。

### 8.4 不适用 / 待决

- **MHAP / overlapInCore 的 Myers 扩展**：pgr 走精确匹配；unitig 间 overlap
  是否允许少量错配（unitig 序列可能有 DBG 错误）待数据验证后再定。
- **气泡/孤儿合并**（`mergeOrphans` bubble 遍）：用户裁定不做。
- **难点**：unitig 边界/包含关系（Canu `PlaceContains` 场景）、重复混淆、
  unitig 无方向（正反链都要查）、以及不同 k 的 unitig 冗余（同一区域多套
  表示，overlap 图会稠密，需去重/合并包含）。其中"unitig 无方向"已解决
  （canonical 索引 + 双向扩展验证，见 §8.5）；重复混淆与冗余去重留 v1。
- **Canu EOL**：只取算法思想，不引入其代码/依赖。

### 8.5 实现后的理解回写（2026-08-12，`design/olc.md`）

实现四命令后，对 bogart / utgcns 的几个语义有了更具体的理解：

- **互惠检查在连接端而非自由端**：Canu 互惠种子（`edgeTo5 && edgeTo3`）的
  语义在延伸时落实为——target 的**连接端**（junction end）的 best edge 必须
  指回当前 unitig，而自由端留给下一步扩展。实现中若误查自由端，线性链会
  在每一步断掉（延伸边和回程边天然指向不同 unitig 的两端）。
- **repeat 双定位的单元化近似**：v0 用"unitig 某端 top2 best 边长度近等
  （≥ 0.9×）且指向不同 unitig → 该端标记 repeat、停止延伸"，是
  `markRepeatReads` 图级双定位在 unitig 层的简化，不携带覆盖度证据。
  合成数据观察：6× 低覆盖时随机基因组仍出现重复区**环形错装**（contig 在
  重复处把基因组两段接反；当时该 unitig 只有唯一 best 边且连的就是错误副本，
  top2 近似不触发）——印证覆盖度证据的必要性，v1 优先补。
  **勘误**：`SPURIOUS_COVERAGE_THRESHOLD=6` / `ISECT_NEEDED_TO_BREAK=15`
  这两个常量**不在 Canu 2.3 源码里**（本笔记初版误标为"Canu"），它们属于
  **Celera 8.3 bogart** 血统（`wgs-8.3rc2/src/AS_BAT/AS_BAT_MergeSplitJoin.C:46-47`：
  "Need to have more than this coverage … to call it a repeat area" /
  "Need at least this number of reads confirming a repeat junction"）。
  Canu 2.3 重写的 repeat breaker（`AS_BAT_MarkRepeatReads.C`，见 §5）改用
  confused-edge + `deviationRepeat=3` + `MIN_ANCHOR_HANG=500` 锚定，无固定
  覆盖度阈值。两者共性（**外部 reads 证据区分重复/唯一**）才是 pgr v1 要补
  的方向；`design/olc.md` §10 沿引的这两常量应改指 celera.md 而非 canu.md。
- **consensus：精确缝合即共识**：unitig 无错 + overlap 全精确时，overlap 已
  把坐标完全对齐，缝合即共识（`asm cns`，重叠区不一致会友好报错）；Canu 的
  template stitch + edlib 重比对 + POA-DAG bestPath 是为高噪声长读设计的。
  §8.3 建议的"asm map 回放 + 多数投票"在 v0 简化为直接缝合（坐标已对齐，
  投票无增量），列投票留 v1（引入错配 overlap 或 junction 不一致时启用）。
- **unitig 无方向的解决**：canonical k-mer 索引天然双向（`canonical_keys`），
  边界 k-mer 查回候选后对 ± 两个方向各做一次扩展验证（`comp` 互补表），
  一条 overlap 记录同时推出两个方向的连接边（flip 位），layout 沿 flip 翻转
  链上后续 unitig 的 strand——对应 Celera/Canu 的 overlap orientation 处理。
- **同 k unitig 通常无端到端精确重叠**：唯一区满足 DBG 最大路径性质（有重叠
  就会被合并）；重复区例外（两个副本共享序列，可产生同 k 重叠，正是 repeat
  检测的用武之地）。跨 k 的 contain/延伸重叠天然存在，使 overlap 图稀疏；
  合成数据 30× 下 3 条 contigs 全部为基因组精确子串。

### 8.6 工程实现细节（对 pgr 的借鉴）

- **内存有界的 all-pairs overlap（`overlapInCore.C:204-219`）**：不是一次把
  所有 reads 的 k-mer 装进哈希表，而是**按 read id 分块**——循环
  `Build_Hash_Index(bgnHashID, endHashID)` 建一块的哈希，再对该块内的 reads
  查 `[bgnRefID, endRefID]` 的参考范围，处理完即释放下一块。→ pgr 的
  `asm ovlp` 目前把全部 unitigs 一次建索引（unitig 数远小于 reads），
  宏基因组数据量大时可参照此"分块建索引 + 滑窗查询"以控内存峰值。
- **并行按块划分（`bogart.C:977-978` / `MarkRepeatReads` blockSize）**：
  Canu 用 `tiLimit` 与线程数把 tig 分成近似均匀的块（`blockSize`），
  repeat 检测等阶段按 tig 并行。pgr 用 rayon 按 unitig 并行，思路一致。
- **overlap 对称化**：`OverlapCache::symmetrizeOverlaps` 保证图边双向一致，
  避免 layout 方向歧义。pgr 的 `ovlp` 已通过 canonical 索引 + 双向验证
  天然对称（§8.5 第 4 条）。
- **二进制存储中间件（`stores/` seqStore/ovlStore/tigStore）**：Canu 把
  reads/overlaps/tigs 存为二进制 store 以便断点续跑与跨进程传递。pgr 的
  `asm olc` 阶段间走内存（`--keep-dir` 才落地文本中间件），当前规模足够，
  暂不需要二进制 store。

## 9. 关键文件清单（速查）

| 组件 | 文件:行 | 内容 |
|---|---|---|
| 流水线 | `pipelines/canu.pl:1144/1165/1179` | cor/obt/utg 三段（main 调度点；子程序定义 `:926/968/1010`） |
| overlap | `overlapInCore/overlapInCore.C:204-219` | 分块建哈希索引（内存有界） |
| overlap | `overlapInCore/overlapInCore-Find_Overlaps.C:235` | k-mer 滑窗查哈希 |
| overlap | `overlapInCore/overlapInCore-Process_String_Overlaps.C:581/442` | 对角合并 / Myers 扩展 |
| layout | `bogart/bogart.C:491-754` | BOG → greedy → repeat breaking |
| layout | `bogart/AS_BAT_BestOverlapGraph.C:67` | best edge + 自适应误差率 |
| layout | `bogart/AS_BAT_PopulateUnitig.C:154-155` | 互惠种子 + 单向延伸 |
| layout | `bogart/AS_BAT_MarkRepeatReads.C:973` | repeat breaking（confused-edge 证据） |
| layout | `bogart/AS_BAT_OverlapCache.C:485` | overlap 过滤/去重/对称化 |
| consensus | `utgcns/unitigConsensus.C:854` | template stitch + edlib 重比对 + DAG |
| consensus | `utgcns/libpbutgcns/AlnGraphBoost.C:386/490` | bestPath + min-coverage 修剪 |
| correction | `correction/falconConsensus.C:69` | FALCON 列 DP consensus |
