# pgr 遮蔽版重复标记实现计划（Dfam 全库）

> 设计笔记，日期：2026-08-03。场景笔记见 [[../repeat-masking.md]]（ir/rept/trf + fa mask
> 的现状与 SD 关系）；FastK 分析见 [[../references/fastk.md]]；RepeatMasker 源码梳理见
> 附录 A（依据仓库内 `RepeatMasker/` 目录，open-4.2.4）。

## 1. 决策

**只实现遮蔽版重复标记**：把 Dfam consensus 全库与目标基因组做比对，输出
runlist 区间交给 `pgr fa mask` 做遮蔽。不做 family/class 注释、不做
`.out/.tbl` 报告、不构建自有重复库。

理由：

*   RepeatMasker 的价值一半在库（Dfam 几十年积累），一半在注释后处理（二十年
    打磨）——两者都不追：**短期无法达到 Dfam 的库积累，注释后处理是泥潭**。
*   直接用 Dfam consensus FASTA 当候选集，用 pgr **原生比对命令**
    `pgr align pgi`（或 `pgr align lastz`）找区间——`pgr sd search --engine pgi`
    已验证 pgi 用于重复区间检测的模式；遮蔽管道 `fa mask` 已存在。
*   遮蔽版只回答一个问题："基因组上哪些区间是重复的"，质量取决于比对敏感度
    （验证见 §5）。

## 2. 方案

### 2.1 核心思路：Dfam 全库 + 一套通用比对

RepeatMasker 完整流程 = 搜索（库 vs 基因组）+ 注释后处理（碎片整合、边界精修、
family/class、K2P、报告，详见附录 A.5）。我们只做**搜索侧的遮蔽**。

**简化点：不按物种分类。** RepeatMasker 的物种逻辑只影响两件事——库筛选
（FamDB 按 `-species` 从 Dfam 取子集）和搜索配方选择（按 primates/rodentia
跑专项 stage），两者都是为注释精度服务的。对纯遮蔽目标，直接**把 Dfam
全库 consensus 当候选集**，一套通用搜索参数跑完：

*   不做 Taxonomy / FamDB / 物种库解析；
*   不按谱系选配方，RepeatMasker 的 17 个配方（附录 A.3）收敛为 1–2 套通用
    参数（可参考 `general_search_parameters` 的量级，放宽 min-len/identity）；
*   这正好对应 RepeatMasker 自己的 `-lib` 模式——官方支持"用户给全库 +
    `general_search_parameters` 一套流程"的简化路径（附录 A.2）。

**代价**：全库比对会让近缘家族的 consensus 互相 hit（cross-family），造成
一定过度遮蔽；物种特异配方（如灵长类年轻 Alu 的"先切除再补搜"）的灵敏度
也会丢失。对"宁可多遮不漏"的遮蔽目标可接受，但要在验证里量化（§5）。

### 2.2 明确不做（RepeatMasker 的注释后处理）

*   碎片整合（`cycleReJoin`）：遮蔽不在乎拼回完整元件，只在乎别漏区域；
*   边界精修：区间级覆盖即可；
*   family/class 注释、K2P %div、`.out/.tbl` 报告：全部不做。

### 2.3 现状保留为对照与兜底

`ir + trf + fa mask` 继续可用：方案落地后与之对比覆盖区间；低复杂度缺口
（polyA 等，§4）由 `pl trf` 兜底。

## 3. 实现步骤

不做注释时实现很轻，基础设施全在：

1.  **库**：直接取 Dfam consensus FASTA 全库（不做物种筛选）；可选加简单重复
    条目，或交给 `pl trf` 兜底低复杂度。
2.  **比对**：使用 `pgr align pgi`（pgr 原生归并比对，输入 FASTA/2bit/.pgi）
    跑"全库 vs 基因组"；或 `pgr align lastz`。`pgr align pgi` 是独立命令，
    不是从 `sd search` 借用的引擎——`sd search --engine pgi` 只是它的一个
    调用场景。注意 `sd search` 的过滤器（>1 kb、>90% identity）是给 SD 调的，
    转座子拷贝分歧大（70–90%），直接跑 `pgr align pgi` 时需自行放宽
    min-len / identity 过滤。
3.  **区间合并**：`spanr cover / merge / fill` 一行管道。
4.  **输出**：runlist → `pgr fa mask`。

工作量比完整 RepeatMasker 小一个数量级。遮蔽版需要的能力映射见附录 A.7。

## 4. 关键风险

*   **比对敏感度**：k-mer（k=17）对高分歧拷贝会漏；pgi 的 syncmer 种子对
    70% identity 的拷贝同样不轻松。这是决定遮蔽质量的核心，必须实测。
*   **全库 cross-family 假阳性**：不做物种筛选后，保守的转座子区域可能让
    基因/其他序列被误遮蔽（over-masking）。遮蔽场景可接受，但需在验证中
    对比"全库 vs 物种库"的遮蔽量差异。
*   **低复杂度缺口**：RepeatMasker 默认屏蔽 low complexity（polyA、卫星、
    homopolymer）。现有 `ir` 只管库内散在重复，`trf` 覆盖串联重复，polyA 这类
    不一定被覆盖。这是遮蔽质量上更实际的差距，与用 k-mer 还是比对无关。
*   **验证基准**：E. coli 几乎无转座子，无参考价值。需用拟南芥/玉米等
    转座子丰富基因组，与 RepeatMasker 的 masked 输出对比 recall。

## 5. 验证实验（实施前调参）

方向已定（§1），验证的目的是评估比对敏感度、确定过滤参数：

1.  取转座子丰富基因组（拟南芥或玉米）；
2.  用 Dfam consensus FASTA **全库**经 lastz/pgi 对一遍，放宽 hit 过滤；
3.  对比数据：
    *   时间与 hits 数量；
    *   覆盖区间 vs 现有 `ir` 的差异；
    *   与 RepeatMasker masked 输出的 recall；
    *   （可选）全库 vs 按物种取库的遮蔽量差异，评估 over-masking 代价。
4.  依据 recall / over-masking / 耗时确定最终过滤参数（min-len、identity、
    gap 模型），落成新命令（如 `pgr fa mask` 的比对模式）的默认值。

## 附录 A：RepeatMasker 源码梳理

> 2026-08-03 依据仓库内 `RepeatMasker/` 目录源码（open-4.2.4，Arian Smit &
> Robert Hubley）梳理，作为 §2 方案的源码证据。

### A.1 概览

RepeatMasker 筛查 DNA 序列中的**散在重复**（interspersed repeats）与**低复杂度
序列**，输出详细注释（family/class、%div/%del/%ins）以及 masked 序列。

*   **实现语言**：Perl（主驱动 + 模块化库），比对本体交给外部搜索引擎
    （RMBlast / crossmatch / NHMMER / ABBLAST），自身只做流程编排与后处理。
*   **主要文件**：
    *   `RepeatMasker`（267 KB）：主驱动——分片、多阶段搜索、汇总 `.cat`。
    *   `ProcessRepeats`（364 KB）：后处理——碎片整合、注释、距离估计、输出。
    *   `DupMasker`：片段重复（segmental duplication）检测。
    *   `RepeatProteinMask`：转座子编码蛋白的蛋白级比对。
    *   `DateRepeats`：按 %div 估计重复插入时间。
    *   `TRF.pm` / `TRFResult.pm`：串联重复（tandem repeat）封装。
    *   `Libraries/RepeatAnnotationData.pm`：family/class 注释数据。
*   **搜索引擎抽象**：`SearchEngineI.pm` 定义接口，
    `CrossmatchSearchEngine` / `WUBlastSearchEngine` / `HMMERSearchEngine` /
    `NCBIBlastSearchEngine`（RMBlast = 特制 blastn）四个实现。
*   **库来源**：Dfam（`DFAM.pm`）/ RepBase（`RepbaseEMBL.pm`）经 FamDB 按物种
    生成 consensus/HMM 库；也可用 `-lib` 直接给 FASTA。库 header 解析
    family/class（如 `>XXX#LTR/ERV`）。

### A.2 主驱动 RepeatMasker

**顶层流程**：

1. **参数解析**：Getopt::Long；配置优先级 CLI > 环境变量 > 配置文件
   （`RepeatMaskerConfig.pm`）。
2. **分片**：`fragmentSize = 60000`、`overlapLen = 2000`（`-frag` 可覆盖，
   且必须 ≥ 2×overlap）。分片边界会产生 edge effect，后处理专门处理。
3. **库解析**：FamDB 按物种取库，或 `-lib` 自定义；`refineableHash.dat`
   记录可精修的元素 id。
4. **并行**：`-parallel N` fork 子进程，每个 batch 独立跑搜索。
5. **每 batch 搜索**：`runSearchStages`（blast 类引擎）或
   `runHMMERSearchStages`（NHMMER）。
6. **汇总**：各 batch 的 `.cat` 拼接成总 `.cat`（含版本/引擎/库 header +
   `## RAW Annotations:` 段）。
7. **后处理**：调用 `ProcessRepeats`（`-nopost` 跳过，供手动分步跑）。

**多阶段搜索（runSearchStages）**：搜索不是一次 blastn，而是**一串配方化搜索 +
迭代切除**：

1. **TRF 阶段**：`runTRFStage` 找 perfect simple repeats（`-nolow` 跳过）。
2. **用户库/物种库**：用 `general_search_parameters` 配方（`-cutoff` 改
   minscore），输出 `*.tmp.custom`。**注意：`-lib` 自定义库时根本不走物种
   分类，直接用这一套通用配方**——这就是正文 §2.1"全库 + 统一参数"简化
   所对应的官方路径。
3. **物种特异阶段**：按分类学（primates/rodentia…）跑对应配方——如灵长类
   先 `cut_young_sines_in_primates` 切 Alu，再 `mask_*` 遮蔽其余 SINE；
   每阶段后 `postProcessSearch` 把已找到的 hit **从当前序列中切除**
   （`excise`），下一阶段只搜剩余序列。

`runHMMERSearchStages` 的注释给出了完整阶段顺序（哺乳类）：sinecutlib →
shortcutlib → cutlib →（人类再扫一次 sinecutlib）→ shortlib → longlib →
mirs.lib → mir.lib → retro.lib → l1.lib → simple.lib → at.lib。每个阶段都
是"搜索 → masklevel 过滤 → 结果筛选 → 切除 → 写 cat"。

**"两轮搜索"的真实含义**：早期版本是两遍 blastn（高严格度找 anchor + 低严格度
补漏），现版本是"多阶段配方 + 切除迭代"——每次切除后剩余序列变少，后续阶段
只处理未覆盖区域。

### A.3 搜索配方（getSearchRecipes）

共 17 个配方，每个定义一套完整搜索参数：

| 配方 | 用途要点 | minscore | minmatch | matrix | bandwidth | masklevel | excise |
| :--- | :--- | ---: | :--- | :--- | ---: | ---: | :--- |
| `perfect_simple_repeats` | 完美简单重复（TRF 替代） | 180 | [8,9,10,11] | simple1 | 1 | 1 | 1 |
| `general_search_parameters` | 通用库搜索（默认路径） | 225 | [8,9,11,13] | 20p##g | 14 | 90 | 0 |
| `cut_young_sines_in_primates` | 灵长类年轻 SINE 切除 | 1200 | [7,8,10,12] | 14p##g | 20 | 1 | 1 |
| `mask_young_sines_in_primates` | 灵长类年轻 SINE 遮蔽 | 1500 | [7,8,10,12] | 14p##g | 20 | 80 | 0 |
| `mask_sines_in_primates` | 灵长类其余 SINE | 225 | [7,8,8,9] | 20p##g | 14 | 80 | 0 |
| `mask_sines_in_non_primate_mammals` | 非灵长哺乳类 SINE | 225 | [6,7,8,10] | 18p##g | 14 | 1 | 1 |
| `general_full_length_repeats` | 全长元件（大带宽跨大 gap） | 300 | [9,10,11,13] | 18p##g | 40 | 1 | 1 |
| `complete_3end_of_young_line1s` | 年轻 LINE1 的 3' 端 | 300 | [9,10,11,13] | 18p##g | 40 | 90 | 1 |
| `older_ALUs_in_primates` | 古老 Alu | 800 | [7,8,10,12] | 14p##g | 20 | 1 | 0 |
| `more_ALUs_in_primates` | 其余 Alu | 400 | [7,8,9,11] | 18p##g | 14 | 10 | 0 |
| `short_repeats_and_satellites_rodents` | 啮齿类短重复/卫星 | 210 | [7,8,9,10] | 25p##g | 14 | 90 | 0 |
| `short_repeats_and_satellites` | 通用短重复/卫星 | 225 | [7,8,10,12] | 20p##g | 14 | 90 | 0 |
| `long_interspersed_repeats` / `ancient_repeats` / `tough_ancient_repeats` / `retroviruses` / `tough_line1s_in_eutheria` | 长散在/古老/反转录病毒等分级搜索 | 各异 | 各异 | 18p##g | 各异 | 各异 | 各异 |
| `simple_repeats_again` / `simple_repeats_flanking` | 低复杂度补搜与侧翼 | 各异 | 各异 | simple 系矩阵 | 小 | 高 | 0 |

`minmatch` 是 4 个值的数组，按严格度档位选（`selectParameter`，受 `-s`/
`-q`/`-qq` 速度档影响）；`matrix` 形如 `20p##g.matrix`（`#` 数量按 GC 含量
调整，见 `runStage` 的矩阵选择逻辑）；`masklevel` 配合 RMBlast 的
complexity-adjusted scoring；`raw=1` 时用 basic scoring。

### A.4 搜索引擎与 RMBlast 参数

`search()` 把配方翻译成引擎参数：

*   `NCBIBlastSearchEngine` → `rmblastn`：
    *   `-word_size 14`（默认；`minmatch` 经 `-word_size` 传入）；
    *   `-dust no`（低复杂度过滤交给 masklevel，而不是 blast 的 dust）；
    *   `-min_raw_gapped_score` = minscore；
    *   `-xdrop_ungap / -xdrop_gap_final / -xdrop_gap` 由 minscore 推导
        （refine 阶段带宽传 `-1` 放宽）；
    *   `-matrix`（配方矩阵）、`-gapopen / -gapextend`（由 gap init/ext 换算）；
    *   默认 complexity-adjusted score mode（`raw=0`），`raw=1` 时 basic mode。
*   **失败重试**：bandwidth > 14 时缩到 14，==4 时缩到 1，minmatch < 10 时
    增大，仍失败则放弃该 batch。
*   crossmatch/WUBlast 用各自的矩阵目录（`Matrices/crossmatch/`、
    `Matrices/wublast/aa/`、`Matrices/ncbi/nt/`）。

### A.5 后处理 ProcessRepeats

输入 `.cat`（raw annotations），输出 `.out` / `.tbl` / `.masked` /
`.align` / GFF / HTML。核心是 `processSequence` 的多 cycle 流水线：

*   **cycle 1**：`cycleReJoin` —— 把被打断的转座子碎片按打分连成链
    （left/right linked hit 结构）。
*   **cycle 2**：移除 edge effect 注释（分片边界产物）、移除 masklevel 违规、
    卫星重复重命名、构建 DNA 转座子等价结构（`%chainBeg/%chainEnd`）。
*   **cycle 3–5**：进一步去碎片化（CYCLE5 维护最近 21 个 join 打分做决策）、
    LTR/LINE 配对（`scoreLINEPair`：按 query/consensus 重叠与名字兼容性打分）、
    等价类间传播。
*   **边界精修**：对 `refineableHash` 中的元素，用更宽松带宽的搜索
    （`-xdrop_gap_final` 放宽）重对齐边界，`replaceRMFragmentChainWithRefinement`
    替换原链。
*   **统计**：`calcKimuraDivergence`（K2P，CpG 位点校正）给出
    %div/%del/%ins；simple/low_complexity 不记 divergence。
*   **输出**：`generateOutput`（.out 标准表）、`generateTableOutput`（.tbl
    汇总）、masked 序列（默认 N，`-xsmall` 小写）、`-a` 时 source alignments。

### A.6 辅助程序

*   `DupMasker`：片段重复检测（独立于库搜索的另一条路径）。
*   `RepeatProteinMask`：用 RepeatPeps 蛋白库跑蛋白比对，找转座子编码蛋白。
*   `DateRepeats`：根据 %div 估计插入时间（转座子年代学）。
*   `TRF.pm`：封装 TRF 4.09+，`runTRFStage` 用它替代 simple consensus 搜索。

### A.7 与 pgr 的对应关系

**遮蔽版需要的能力**（正文 §3）：

| 步骤 | pgr 现状 | 判断 |
| :--- | :--- | :--- |
| 库-基因组比对 | `pgr align pgi`（原生）或 `pgr align lastz` | 高可行，预计比 RMBlast 快一个数量级 |
| 区间合并/覆盖 | `spanr cover / merge / fill` | ✅ 已有 |
| 输出遮蔽 | `pgr fa mask --runlist` | ✅ 已有 |
| 低复杂度兜底 | `pgr pl trf` | 已有（缺口见正文 §4） |

**遮蔽版明确不做**（对应 RepeatMasker 的注释后处理，A.5）：

| 步骤 | 说明 |
| :--- | :--- |
| 碎片整合（cycleReJoin） | 遮蔽不在乎拼回完整元件，只在乎别漏区域 |
| 边界精修 | 区间级覆盖即可 |
| family/class 注释 | 不做 |
| K2P %div | 不做 |
| `.out/.tbl` 报告 | 不做 |

**其余启示**：

*   搜索侧 = 配方化 blastn + 迭代切除，用 pgi/lastz 替换理论上可行且更快；
*   真正的工程量和价值在 `ProcessRepeats` 的 6+ cycle 后处理——这正是遮蔽版
    不做的部分；
*   低复杂度（Simple_repeat / Low_complexity）在 RepeatMasker 里由
    TRF + simple.lib 覆盖，对应 pgr 的 `pl trf`，是遮蔽质量的实际差距。
