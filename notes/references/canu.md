# Canu（2.3）：OLC 组装器源码分析（overlap / layout / consensus）

> 2026-08-12 整理，纯源码分析（`canu-2.3/`，版本 2.3，GitHub `marbl/canu`）。
> Canu 是 Celera Assembler r4587（`wgs-assembler`）的 fork，面向高噪声单分子
> 长读（PacBio CLR/ONT）。**与 pgr 的关系**：不是要对 reads 做 OLC，而是用户
> 的设计意图——**把不同 k 各自生成的 unitigs 当"伪 reads"，在 unitig 层做
> OLC 拼接**（见 §8）。本文档记录 Canu 的 OLC 三组件源码结构，以及该设计意图
> 对应的借鉴评估。

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
├── overlapAlign/            # overlap 的碱基比对（计算 alignments）
├── stores/                  # seqStore / ovStore / tigStore（二进制存储）
└── utility/                 # 基础设施（bits、strings、system、libbacktrace）
```

## 3. 流水线（`pipelines/canu.pl`）

三段各跑一遍 `meryl 计数 → overlap → 下游`：

```
doCorrection(1144)  meryl → overlap(cor, mhap) → buildCorrectionLayouts →
                    filterCorrectionLayouts → generateCorrectedReads(falconConsensus)
doTrimming(1165)    meryl → overlap(obt, mhap) → trimReads / splitReads
doUnitigging(1179)  meryl → overlap(utg, mhap) → readErrorDetection →
                    overlapErrorAdjustment → unitig(bogart) → consensus(utgcns)
                    → generateOutputs
```

- overlapper 默认全部 `mhap`（`canu.pl:125-127`）；`overlapInCore` 走
  `ovlOverlapper=ovl`，`canu.pl:718` 有醒目的 `DO-NOT-USE` / `LUDICROUSLY SLOW`
  警告。
- correction 是 Canu 相对 Celera 8 的新增阶段：Celera 8 直接用原始 reads
  做 overlap；Canu 先把每条 read 按 overlap 布局纠错成 consensus 再组装。

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
  `erateMax` 重算（自适应阈值）。
- **greedy 建 tig**（`AS_BAT_PopulateUnitig.C`）：种子 read 要求两端 best edge
  互惠（`edgeTo5 && edgeTo3`），沿 best edge 单向延伸、已放置即停；非互惠
  的 read 不种子（避免错装）。
- **repeat breaking**（`AS_BAT_MarkRepeatReads.C:973`）：用 AssemblyGraph 找
  "与 tig 外 read 也有 overlap"的 read（双定位），投影回 tig，若无 read 跨越
  该区域则打断。这是 OLC 处理重复的核心手段。
- **气泡**（`mergeOrphans` 第二遍，`bogart.C:639`）：按相似度阈值
  `similarityBubble` 合并平行路径——正是用户裁定"不处理"的那种启发式。

## 6. consensus（utgcns）

`unitigConsensus::generate`（`unitigConsensus.C:1466`）按算法分发：

- **默认 `P` = pbdagcon**（`generatePBDAG:854`）：
  1. `generateTemplateStitch:259`：按 layout 顺序用 edlib（Myers O(ND)）把 reads
     逐步"缝"成一条模板（带 homopolymer 压缩坐标，`switchToUncompressedCoordinates:164`）。
  2. 所有 reads 用 edlib 重比对到模板（`alignEdLib`：band 递增、错误率递增重试，
     `MAX_RETRIES`）。
  3. 比对建成 PacBio 的 AlnGraphBoost POA-DAG（`libpbutgcns/AlnGraphBoost.C`），
     `mergeNodes:188` 合并等价节点。
  4. `bestPath:490`（heaviest path）→ `consensusNoSplit:386`：沿 best path 取
     base，低于 `minCoverage` 的段截掉（Canu 的 min-coverage 修剪）。
- **`Q` = quick**（`generateQuick:976`）：只做 template stitch，不做 consensus——
   拼贴序列直接当 contig（文档说适合做中间检查/抛光输入）。
- **`S` = singleton**（`generateSingleton:1012`）：单 read 原样输出。
- **correction 的 consensus**（`correction/falconConsensus.C:69`）：FALCON 移植，
  证据 reads 对齐到模板后按列 DP + link 回溯打分（`getConsensus`），不是
  AlnGraphBoost。

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
  表示，overlap 图会稠密，需去重/合并包含）。
- **Canu EOL**：只取算法思想，不引入其代码/依赖。

## 9. 关键文件清单（速查）

| 组件 | 文件:行 | 内容 |
|---|---|---|
| 流水线 | `pipelines/canu.pl:1144/1165/1179` | cor/obt/utg 三段 |
| overlap | `overlapInCore/overlapInCore-Find_Overlaps.C:235` | k-mer 滑窗查哈希 |
| overlap | `overlapInCore/overlapInCore-Process_String_Overlaps.C:442` | Myers 扩展 |
| layout | `bogart/bogart.C:491-754` | BOG → greedy → repeat breaking |
| layout | `bogart/AS_BAT_BestOverlapGraph.C:67` | best edge + 自适应误差率 |
| layout | `bogart/AS_BAT_PopulateUnitig.C` | 互惠种子 + 单向延伸 |
| layout | `bogart/AS_BAT_MarkRepeatReads.C:973` | repeat breaking |
| consensus | `utgcns/unitigConsensus.C:854` | template stitch + edlib 重比对 + DAG |
| consensus | `utgcns/libpbutgcns/AlnGraphBoost.C:386/490` | bestPath + min-coverage 修剪 |
| correction | `correction/falconConsensus.C:69` | FALCON 列 DP consensus |
