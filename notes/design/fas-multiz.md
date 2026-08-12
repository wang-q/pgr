# fas-multiz 设计与实现

> **状态：已实现**（2026-08-03）——`libs::fas_multiz` 与 `pgr fas multiz` CLI 均已落地； DP 引擎为
> multiz yama 的直接移植（§3.4、§4.2），合并次序与打分已确定化（§4.1）。

本文档是 `libs::fas_multiz` 的设计与实现说明，基于对 multiz 源码的分析 （见
[multiz.md](../references/multiz.md)）。

## 1. 定位：与 multiz 的关系

`pgr` 里与 multiz 相关的工作有两条路径，定位不同，但**不存在策略层面的对立**：

| 路径                      | 做法                     | DP  | 输出             | 适用场景                     |
|---------------------------|--------------------------|-----|------------------|------------------------------|
| `pgr pl p2m` + `fas join` | 区间交集 + 机械堆叠      | 无  | 严格交集（core） | 参考骨架一致，只要共同覆盖区 |
| `pgr fas multiz`          | yama 直译的 profile 合并 | 有  | union/mesh       | 参考 gap 不一致，要完整比对  |

`p2m` 之所以"快"，不是因为它采用了另一种比对策略，而是主动放弃了 multiz 会做的两件事——
**gap 冲突求解**与**边缘覆盖保留**——只回答"所有输入共同覆盖的 区域长什么样"。`fas-multiz` 则与
multiz 是同一套策略：以参考坐标为主轴做一体化 DP，在合并时实时解决参考 gap 冲突，保留边缘覆盖。

选择建议：只要严格交集 → `p2m`；要 multiz 同款完整比对 → `pgr fas multiz`。
两者是"快而粗"与"慢而全"的取舍。

### 1.1 p2m 快速交集捷径（对照）

`pgr pl p2m` 内部依次：`fas cover` 提取各输入参考区间（`--trim 10` 切边） → `spanr intersect`
取共有区域 → `fas slice` 按交集切片 → `fas join` 按参考坐标 Key 机械堆叠 → `fas name` / `fas subset`
统一物种列顺序。

局限：不处理 gap——它假设各输入在同一坐标下的参考片段完全一致，否则直接堆叠 会让非参考序列错位。
`fas refine` 是独立的可选后处理，不在 `p2m` 流程内。

## 2. 设计

### 2.1 核心概念：profile 合并

把 multiz（以及 fas-multiz）理解为 **profile 对 profile 的渐进比对** 是最直观
的方式：

- 一个 **profile** 就是一段多序列比对：参考序列 + 若干物种，按列排列（例如
  一个 `.fas` block 或一段 MAF alignment）。
- 合并两个输入时，对齐的最小单位是**整列**，而不是单个碱基对：C/D/I 三状态
  分别表示"A 列配 B 列 / A 列插入、B 全 gap / B 列插入、A 全 gap"。
- 打分是 **sum-of-pairs**：对列内所有物种对（K×L 笛卡尔积）求和，列内多个
  物种共同决定哪两列对齐——这是 profile–profile 比对与序列两两比对的根本区别。
- 多物种不是一次全对齐，而是**渐进合并**：按合并次序（multiz 用 guide tree，
  pgr 用内容驱动的贪心，§4.1）反复两两合并，每次把一个 profile 并进累计的
  profile，最终得到全物种比对。

一句话：multiz = 以参考坐标为锚的 profile–profile 渐进比对；`libs::fas_multiz`
复刻同一语义，只是载体从 MAF 换成 block FA。

### 2.2 目标与输入/输出

- **目标**：把 `k`（≥2）个 block FA 文件（如 pairwise 派生 `.fas`），在共享 参考坐标系下合并为
  union/mesh 风格的 `.fas`。
- **输入**：`k` 个 `.fas` 文件，block 中均含同名参考序列；可选一个交集区域 限制计算范围。
  所有输入的参考序列应来自同一基因组版本（建议相同的 masking/裁剪流程）。
- **输出**：合并后的 `.fas`——交集区域内与 `p2m + join` 兼容；边缘/非完全 交集区域尽量保留各输入的对齐
  （union 行为）。

### 2.3 合并语义（union）

行为固定为 union（与 multiz 一致，不提供 core 模式）：

- 窗口：任一输入覆盖即保留（窗口推导见 §3.2）。
- 物种：取各输入的并集，缺失输入以 gap 填充。
- 单步合并失败时：起始两块失败 → 回退保守合并；后续某步失败 → 跳过该输入
  继续（见 §3.3）。

### 2.4 与 multiz 的异同

**继承的部分**：

- 以参考坐标为主轴，在参考坐标上定义窗口/段落，窗口内做 profile 合并。
- 带状 DP 限制搜索空间，而非全空间 MSA。
- union/mesh 语义：尽量保留真实比对关系，只在必要时取舍冲突列。

**pgr 的调整**：

- **工作层级**：multiz 直接对齐 MAF profile；fas-multiz 对齐 block FA，上游由 `pgr axt/maf to-fas`
  规整，语法与元数据更简单。
- **合并次序**：multiz 按 guide tree 渐进合并；fas-multiz 用内容驱动的贪心次序（§4.1）。
- **职责边界**：multiz 是独立 C 项目；fas-multiz 是 pgr 的 libs 模块，与
  `fas cover/slice/join/refine` 共享 `.fas` 生态。

DP 引擎本身未做简化——`banded_align.rs` 是 yama 的直接移植（§3.4、§4.2）。

### 2.5 命令行接口

`pgr fas multiz -r <ref> <in.fas>... [选项]`

- `-r, --ref-name <NAME>`（必需）：参考序列名，须存在于所有输入。
- `<infiles>...`：至少 2 个 block FA 文件。
- `--radius <INT>`（默认 30）：带状 DP 半径；同时用于窗口推导时的区间扩展
  （§3.2）。
- `--min-width <INT>`（默认 1）：最小窗口宽度。
- `-o, --outfile`：输出文件（默认 stdout）。

打分不提供任何 CLI 参数，硬编码 multiz 原版值：HOX70（= `hoxd55`）矩阵 +
gap open 400 / extend 30（§3.5）。

### 2.6 与原版 CLI 的行为差异（待决策）

当前实现与原版 multiz 的行为差异，是否对齐原版尚未决定：

- **未使用块与单行块输出**：原版支持 `out1 out2`（收集未参与合并的块）与
   `all`（默认不输出单行 block，指定 `all` 才输出）；pgr 当前没有这两个机制——
   合并失败的窗口直接丢弃，单输入窗口的 block 固定保留。
   候选：实现 `out1/out2` 与 `all`，对齐原版。

## 3. 实现

### 3.1 模块结构

| 文件              | 职责                                                                                                              |
|-------------------|-------------------------------------------------------------------------------------------------------------------|
| `mod.rs`          | 类型（`FasMultizConfig`/`Window`）与文件级入口 `merge_fas_files` / `merge_fas_files_auto_windows` |
| `windows.rs`      | 窗口推导                                                                                                          |
| `merge.rs`        | `merge_window`：progressive 合并与保守回退                                                                        |
| `banded_align.rs` | 单步 profile–profile 带状 DP（yama 直译）                                                                         |
| `tests.rs`        | 13 个单元测试                                                                                                     |

依赖：`fmt::fas`（block 数据结构）、`chain::SubMatrix`（打分，硬编码 `hoxd55`）、
`ds`（区间合并/覆盖计数）。不依赖 `libs::alignment`。

### 3.2 窗口推导（merge_fas_files_auto_windows）

1. 提取所有输入参考序列的 `Range`，按 `radius` 扩展后按染色体合并重叠区间。
2. 过滤宽度小于 `min_width` 的窗口。
3. 无需显式覆盖度过滤：每个窗口都源自某个输入参考区间按 `radius` 扩展，天然被该输入
   覆盖（零宽单碱基块也能保留，`merge_intervals` 不会丢弃）。

### 3.3 单窗口合并流程（merge_window）

1. 每个输入取窗口内第一个与参考重叠的 block；无 block 的输入直接跳过。
2. block ≥ 2 时尝试 progressive DP 合并（§3.4）：
    - 按内容驱动的确定性次序（§4.1）两两合并；
    - 单步合并先要求参考去 gap 后相同（ungapped equal）→ 走带状 DP；参考
      去 gap 不同 → 尝试 crossover 拼接（需共享非参考物种打分，见 §4.3）；
    - 起始两块合并失败 → 整体回退保守合并；后续某步失败 → 跳过该输入继续，
      已合并的部分保留。
3. 保守合并：所有 block 的参考 entry 完全相同（含 gap）时直接堆叠，取物种并集；
   否则放弃该窗口。

### 3.4 核心算法：yama 直译（banded_align_refs）

单步把两个 block 当作 profile，在参考×参考的带状网格上做 DP：

- **C/D/I 三状态**：C=替换，D=删除（A 列配 B 全 dash），I=插入（B 列配 A 全 dash）；每格记录
  `flag_c | flag_d<<2 | flag_i<<4` 供回溯。
- **准自然 GAP 查表**：16 种"最后两条边"构型中 6 种收 gap_open（`mz_scores.c` 直译）。
- **全体物种对打分（K×L 笛卡尔积）**：base–base 用替换矩阵原始值（`hoxd55`），
  base–gap 收一次 gap_extend（-30），gap–gap 为 0；I/D 按"插入列非 dash 数 ×
  对方行数 × gap_extend"收费。
- **端部 gap 免费**：I 在末行、C 在起点列、D 在起点/终点列不收 gap-open （extend 照收）。
- **参考锚定 LB/RB**：参考去 gap 逐位配对，平滑"香肠"扩展后每行只在 `[lb[i], rb[i]]` 内计算
  （变宽带）。
- **边缘裁剪**：只删两端"全物种都是 gap"的列（旧"单侧列"规则会误删真实内容）。

### 3.5 打分参数

打分全部硬编码 multiz 原版值（不提供 CLI 参数）：

- 替换矩阵：`chain::SubMatrix::hoxd55()`（= multiz HOX70，A-A 91 / A-C -114
  等原始值，不做缩放）。
- gap：open 400 / extend 30；`SS('-',x) = -extend`，GAP 表按构型收 open。

## 4. 实现演进记录

> 以下按时间记录关键决策；当前行为以 §3 为准。

### 4.1 阶段 1：全体物种 SP 打分 + 确定性合并次序（2026-08-03）

- **打分升级**：从"只对两个 block 的共同物种打分"改为"共享物种自对 + 参考交叉对" （星型拓扑）。
  两个实证发现（S288c 三输入，参考去 gap 3826 bp）：
    - 朴素全笛卡尔积会拖偏参考锚点：真实数据上产生 155 列错位（3981 bp vs 3826），因 base–gap
      逐对罚分远大于错配罚分，DP 倾向移位避开 gap 列。
    - 对角线上 gap 贡献改为 0，gap 成本全部由转移承担，避免与替换分数失衡。
    （该打分方案在阶段 1.5 被"全体物种对（K×L 笛卡尔积）"替换，见 §4.2；当前行为以 §3.4 为准。）
- **确定性合并次序**：不按输入顺序，而按内容贪心——先选物种最多的块，之后反复
  选与累计物种集重叠最大的块，并列时按内容键（参考区间 + 物种名）打破；输入 顺序无关已有回归测试。
- **配套修正**：合并参考只保留第一块的参考序列（不再用第二块碱基填 gap）， 保证合并参考恒等于输入参考
  （3826 bp）。

### 4.2 阶段 1.5：yama 引擎直译（2026-08-03）

对照 `multiz-multiz/mz_yama.c` / `mz_preyama.c` / `mz_scores.c` 升级单步 DP：

- C/D/I 三状态 + 准自然 GAP 查表 + 全体物种对打分（§3.4）。
- **LB/RB 必要性实证**：无锚定时自由端 gap 会把列数差整段堆到块端，随后被 边缘裁剪删掉真实内容
  （Spar 4057→3724）；锚定后 overhang 被限制在 radius 内。
- **验证**：S288c 三输入合并 4193 列，参考去 gap 恒 3826， RM11_1a 3834 / Spar 4057 / YJM789
  3822，与输入逐碱基一致、零丢失；`merge_window_preserves_species_content` 等测试全绿。

### 4.3 P2 状态（2026-08-03）

- **打分参数硬编码 ✅**：`--score-scheme` / `--gap-model` / `--align-gap-open` /
  `--align-gap-extend` / `--match-score` / `--mismatch-score` / `--gap-score`
  已全部删除，矩阵与 gap 罚分硬编码 multiz 原版值（§3.5）。
- **hox70 别名（不做）**：multiz 的 HOX70 与 pgr 的 `hoxd55` 数值完全相同 （91/-114/-31/-123，gap
  open 400/extend 30），直接用 `hoxd55` 即可，不引入 `hox70` 别名。
- **v=0 模式（未做）**：需要第二次 yama 对齐参考行；pgr 渐进合并以累计块 参考为锚（等价 multiz
  v=1），价值不明，不做推测性实现。
- **multiz 回归对比（部分受阻）**：MAF vs fas 格式、分块语义不同，字节级对比 不可行。
  实测暴露一个预先存在的输入模型限制：部分重叠块若只共享参考且参考 去 gap 不等，合并被拒绝
  （crossover 需要共享非参考物种打分）。主流程（切片 输入）不受影响；修复需重新设计"仅参考共享"的
  crossover，留待真实需求。

## 5. 已知局限与后续方向

- **渐进合并视野**：每步 DP 只看到"当前累计块 + 一个新输入"，尚未实现一次 决策看到全体输入的 SP-DP。
- **DP 网格**：参考×参考二维带状网格（与 multiz yama 的成对 profile 合并 一致）；非参考物种通过 SP
  打分参与决策，但没有独立坐标轴。
- **多块部分重叠输入**：见 §4.3。
- **演进方向**（仅当出现真实需求）：
    1. 局部小 K 多输入 DP（K≤3）：困难区域用 exact 多维 DP 精修。
    2. 多轨迹 2D 网格近似：以"参考坐标 × 合法状态"综合多输入打分。
    3. DP 失败时更智能的降级策略（退回 `fas join`、标记未合并等）。

## 6. 真实数据验证（2026-08-12，multi6 4 路酵母 multiz）

承接 §4.3 的"留待真实需求"与 `chain-algorithms.md` §12.3：用本地真实多基因组
比对（`~/data/egaz/multi6/`，S288c + RM11_1a/YJM789/Spar 的 4 路 multiz）
验证 `pgr fas multiz` 的合并与 `best_crossover` 拼接。

**数据与转换**：三个 pairwise synNet MAF（`Pairwise/*/mafSynNet/<chr>.synNet.maf.gz`，
与真实 4 路 multiz 的输入完全一致，见 multi4_mz MAF 头部的 multiz 命令行）经
`pgr maf to-fas` 转 block FA；真实 4 路输出 `multi4_mz/*.maf.gz` 同样转换作对照。

### 6.1 发现 1：原始 pairwise 输入合并覆盖丢失 87–100%

16 条染色体（MITO 无 pairwise 数据）直接喂 `pgr fas multiz -r S288c`，输出参考
覆盖 vs 三输入参考覆盖并集：

| chr | 输入并集 bp | 合并保留 bp | 保留率 |
| :--- | ---: | ---: | ---: |
| I | 223,412 | 1,593 | 0.7% |
| II | 792,899 | 771 | 0.1% |
| III | 309,858 | 268 | 0.1% |
| IV | 1,495,312 | 1,878 | 0.1% |
| V | 561,433 | 0 | 0.0% |
| VI | 260,888 | 3,457 | 1.3% |
| VII | 1,048,553 | 131,423 | 12.5% |
| VIII | 543,091 | 0 | 0.0% |
| IX | 431,058 | 245 | 0.1% |
| X | 718,535 | 58,023 | 8.1% |
| XI | 665,092 | 269 | 0.0% |
| XII | 1,026,841 | 33,246 | 3.2% |
| XIII | 905,346 | 0 | 0.0% |
| XIV | 764,150 | 6,128 | 0.8% |
| XV | 1,055,602 | 1,357 | 0.1% |
| XVI | 916,193 | 0 | 0.0% |

保留率仅 0–12.5%，多数 <1%（VII/X 例外：这些染色体存在两输入 block 参考坐标
恰好一致的区段，可走正常 DP 堆叠）。与设计 §2.3 承诺的 union 语义（"任一输入
覆盖即保留"）严重不符。

### 6.2 根因（两个叠加的输入模型限制）

1. **每窗口每输入只取一个 block**：`merge_window` 用 `group.iter().find(...)`
   取第一个与窗口重叠的 block；真实 pairwise MAF 的 block 结构在窗口内是
   多块、彼此参考坐标不连续的（如 YJM789 在 chr I 12907–17133 内碎成
   12907–13005 / 13018–13062 / 13063–17133 多块），只有第一块参与合并，
   其余整段丢失。
2. **不同物种对输入既无相同参考坐标也无共享非参考物种**：正常 DP 路径要求两
   block 参考去 gap 相等（`ungapped_equal`），`best_crossover` 路径要求共享
   非参考物种打分；真实 pairwise 输入（S288c+A、S288c+B、S288c+C）两者都不
   满足，只有参考坐标恰好一致的块能堆叠（§6.1 的 VII/X 例外即来自此）。

切片到三输入公共交集（`fas slice`，chr I 184,085 bp）后仍大量丢弃——切片保留
各输入自己的 block 结构，不改变根因 1/2。对照：真实 multiz 用同一组 pairwise
MAF 产出完整 4 路（`multi4_mz`），证明输入数据本身没问题。

### 6.3 发现 2：真实数据上无自然参考冲突

chr I 三输入间 163 对重叠参考区间，去 gap 后全部一致（0 冲突）——所有输入对齐
同一 S288c 参考，`best_crossover` 在干净真实数据上不会自然触发（符合预期）。
真实场景要触发它需要不同参考版本/构建的输入。

### 6.4 发现 3：`best_crossover` 机制本身在真实尺度下正确（受控测试）

取 chr VII 真实 130 kb 连续 block（S288c.VII 405390–535219 + RM11_1a），构造
受控参考冲突（block0 右半 10% SNP、block1 左半 10% SNP；注意语义：crossover
取 block0 左半 + block1 右半，故 block0 左半须可信、block1 右半须可信）：

* 合并成功不丢弃；输出参考与原始序列 **100% 一致**（129,995/129,995 bp），
  切点落在中部；
* 共享物种（RM11_1a）与仅单侧存在的物种（YJM789，只出现在 block0）**100%
  保留**（去 gap 后逐碱基一致）——§4.3 前的"单侧物种跨切点塌成 gap"数据丢失
  bug 类在真实尺度下确认未复发。

### 6.5 结论

`best_crossover` 拼接机制：**真实尺度验证通过**（§12.3 该项闭环）。
`pgr fas multiz` 的窗口合并编排：**真实多两两输入下不满足设计 union 承诺**，
覆盖丢失 87–100%（§6.1）；修复需按 §5 演进方向重做窗口合并（每输入取全部
重叠块、或"仅参考共享"的 crossover 打分），已登记待办。

### 6.6 修复（2026-08-12，对照 multiz.c 块流合并重写）

按 multiz 原版（`multiz-multiz/multiz.c` + `notes/references/multiz.md`
§2.1/§2.2）重写窗口合并为**逐重叠区的块流合并**，16 条染色体覆盖全部恢复：

**实现（`src/libs/fas_multiz/`）**：

* `merge_window` 收集每输入的**全部**重叠块并按参考坐标排序，形成单覆盖
  block 流（canonical 输入模型）；流按首块 content key 排序后逐对合并，
  与输入文件顺序无关；
* `merge_two_streams` 镜像 `multiz()` 主循环：非重叠前端块直接输出 →
  前端部分（`slice_part`，含两侧插入列）→ 重叠区 `[beg,end]` 切片
  （`slice_overlap`，pre_yama 语义）后交给既有 DP 合并 → `keep_from`
  （keep_ali）保留未消费尾部 → 尾随插入列单独输出；
* 参考坐标→列映射 `ref_pos_to_col`（mafPos2Col）、列切片 `slice_block_cols`
  （make_part_ali_col + mafColDashRm）、`leading_insertions`/`trailing_insertions`
  （gap-front / tail）；
* 无参考行的插入块在后续流合并中直接透传（canonical 里 unused 走独立文件，
  pgr 输出流等价处理）；
* 流推进后补 canonical 的 `if (a1->end < a2->start) continue` 重检——这是
  覆盖丢失的**最后根因**：a2 前端循环推进后 a1 位于新 a2 之前时，直接进重叠
  段会把 a1 剩余块整个丢掉。

**验证（multi6 4 路酵母，对照真实 multiz）**：

* **参考覆盖**：16/16 染色体 = 输入并集 **100%**（修复前 0–12.5%；chr
  I 223,412 bp 逐 bp 保留，耗时 1.2s vs 修复前 >240s）；
* **参考碱基**：与真实 4 路 multiz **100% 一致**（10,996,637/10,996,637）；
* **物种内容**：参考对齐位逐碱基保真 RM11_1a 99.999%、YJM789 99.999%、
  Spar 99.995%（残余为插入边界坐标漂移）；
* **列级对照**：RM11_1a 100%、YJM789 99.7–100%、Spar 97.9–98.5% 与真实
  multiz 一致（差异 = 渐进合并次序不同导致，pgr 用 content 序、multiz 用
  guide tree；参考碱基与物种内容均保真，差异限于插入列摆放与少量列对齐
  选择）。

**行为变化（需知）**：输出从"每窗口一块"变为 canonical 多块结构（前端/
  重叠/尾部各自成块），与真实 multiz 输出结构一致；合成测试数据中 ref
  range 与序列长度不一致的两处已按 MAF 语义修正（range 大小 = ref 非 gap
  碱基数）。

### 6.7 残余列差异的本质（2026-08-12 分析，结论：非误差、不解决）

修复后 pgr 与真实 multiz 的列级对照：RM11_1a 100%、YJM789 99.7–100%、
Spar 97.9–98.5%。对 Spar 的 ~2% 差异做归因分析：

* **差异构成**（全染色体 188,399 个差异位）：65% 插入列摆放（gap-vs-base）、
  35% 真实错配（base-vs-base，65,992 个）；
* **错配成段**（chr IV 8,620 个聚成 2,103 段，83% 为连续段、最长 38 位，
  pgr/multiz 交替"对/错"）——典型**对齐滑移**：合并时 indel 列选择不同，
  导致一段列整体错位，属等价对齐而非错误；
* **决定性判据——输入保真度**：在 65,992 个真实错配处，**99.997%
  （65,977 个）pgr 的碱基与输入 pairwise（S288cvsSpar）一致，multiz 侧仅
  7 个**。全染色体参考对齐位统计：pgr 偏离输入 RM 0.000% / YJM 0.000% /
  Spar 0.001%；multiz 偏离输入 Spar **0.600%**（66,111 位）。

**结论**：残余差异不是 pgr 的误差——pgr 是输入 pairwise 的**忠实 union**
（逐碱基保留每个输入的对齐），差异全部来自 multiz 渐进合并时的**列重排**
（其 yama DP 重新选择了列摆放）。参考碱基 100% 一致、物种内容 99.995%+
保真，说明 pgr 的"参考坐标上的序列语义"完全正确。要消除差异只能复刻
multiz 的 guide tree 合并次序（pgr 无树输入，不可行）或复刻其列重排（会
破坏 pgr 的输入保真优势且无客观收益）。**判定：不解决，记录差异本质**
（下游 to-vcf/to-fas 以参考坐标为准，列滑移不影响 SNP 语义）。另注：
`pgr fas multiz` 的输出照管线惯例还要过 `pgr fas refine`（builtin POA /
mafft 对每个 block 重新 MSA），真实 multiz 输出同样如此——列摆放差异在
refine 阶段会被重新优化，进一步确认其对最终产物无影响。

## 附录：上游链路 LASTZ 与链化

UCSC 典型 WGA 流程中，multiz 位于 "pairwise 比对 + 链化 + net + mafFromNet" 之后，只消费整理好的
MAF。`pgr` 的前置封装：

- `pgr align lastz`：LASTZ 前端，内置 UCSC preset（`set01`..`set07`）、自动加
  `--format=lav --markend --ambiguous=iupac` 等选项，支持目录递归与 rayon 并行。
- `pgr psl chain`：axtChain 风格链化，`SubMatrix` 打分 + `GapCalc` gap 曲线
  （`--gap-model loose|medium` 或显式仿射），KD-tree 加速前驱搜索。

两者覆盖 "blastz/lastz 比对 + axtChain 链化" 两步，为 fas-multiz 提供 pairwise 对齐基础。
