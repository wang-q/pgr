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

### 2.1 目标与输入/输出

- **目标**：把 `k`（≥2）个 block FA 文件（如 pairwise 派生 `.fas`），在共享 参考坐标系下合并为
  union/mesh 风格的 `.fas`。
- **输入**：`k` 个 `.fas` 文件，block 中均含同名参考序列；可选一个交集区域 限制计算范围。
  所有输入的参考序列应来自同一基因组版本（建议相同的 masking/裁剪流程）。
- **输出**：合并后的 `.fas`——交集区域内与 `p2m + join` 兼容；边缘/非完全 交集区域尽量保留各输入的对齐
  （union 行为）。

### 2.2 合并语义（union）

行为固定为 union（与 multiz 一致，不提供 core 模式）：

- 窗口：任一输入覆盖即保留（窗口推导见 §3.2）。
- 物种：取各输入的并集，缺失输入以 gap 填充。
- DP 失败时跳过该输入继续合并，而不是整体放弃。

### 2.3 与 multiz 的异同

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

### 2.4 命令行接口

`pgr fas multiz -r <ref> <in.fas>... [选项]`

- `-r, --ref-name <NAME>`（必需）：参考序列名，须存在于所有输入。
- `<infiles>...`：至少 2 个 block FA 文件。
- `--radius <INT>`（默认 30）：带状 DP 半径。
- `--min-width <INT>`（默认 1）：最小窗口宽度。
- `--score-scheme <文件|预设>`：替换矩阵（LASTZ 格式文件，或 `hoxd55` 预设； 默认 `hoxd55`）。
- `--gap-model constant|medium|loose`（默认 `medium`）；`--align-gap-open` / `--align-gap-extend`
  显式覆盖。
- `--match-score`（默认 2）/ `--mismatch-score`（默认 -1）/ `--gap-score`（默认 -2）。
- `-o, --outfile`：输出文件（默认 stdout）。

### 2.5 与原版 CLI 的行为差异（待决策）

以下两点是当前实现与原版 multiz 的行为差异，是否对齐原版尚未决定：

1. **打分参数**：原版硬编码 HOX70（= `hoxd55`）+ gap open 400 / extend 30；
   pgr 当前暴露 `--score-scheme` / `--gap-model` / `--align-gap-open` /
   `--align-gap-extend` / `--match-score` / `--mismatch-score` / `--gap-score`
   供配置（其中 `--mismatch-score` 目前只存在于配置中、不参与打分）。
   候选：删除这些参数，硬编码原版值。
2. **未使用块与单行块输出**：原版支持 `out1 out2`（收集未参与合并的块）与
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
| `tests.rs`        | 10 个单元测试                                                                                                     |

依赖：`fmt::fas`（block 数据结构）、`chain::SubMatrix`/`GapCalc`（打分）、 `ds`（区间合并/覆盖计数）。
不依赖 `libs::alignment`。

### 3.2 窗口推导（merge_fas_files_auto_windows）

1. 提取所有输入参考序列的 `Range`，按 `radius` 扩展后按染色体合并重叠区间。
2. 过滤宽度小于 `min_width` 的窗口。
3. 过滤覆盖度：只保留被至少一个输入覆盖的窗口。

### 3.3 单窗口合并流程（merge_window）

1. 每个输入取窗口内与参考重叠的 block；无 block 的输入直接跳过。
2. block ≥ 2 时尝试 progressive DP 合并（§3.4）：
    - 按内容驱动的确定性次序（§4.1）两两合并；
    - 任一步失败 → 跳过该输入继续；全部失败则回退保守合并。
3. 保守合并：所有 block 的参考 entry 完全相同（含 gap）时直接堆叠，取物种并集；
   否则放弃该窗口。

### 3.4 核心算法：yama 直译（banded_align_refs）

单步把两个 block 当作 profile，在参考×参考的带状网格上做 DP：

- **C/D/I 三状态**：C=替换，D=删除（A 列配 B 全 dash），I=插入（B 列配 A 全 dash）；每格记录
  `flag_c | flag_d<<2 | flag_i<<4` 供回溯。
- **准自然 GAP 查表**：16 种"最后两条边"构型中 6 种收 gap_open（`mz_scores.c` 直译）。
- **全体物种对打分（K×L 笛卡尔积）**：base–base 用替换矩阵（/50 缩放）， base–gap 收一次 gap_extend，
  gap–gap 为 0；I/D 按"插入列非 dash 数 × 对方行数 × gap_extend"收费。
- **端部 gap 免费**：I 在末行、C 在起点列、D 在起点/终点列不收 gap-open （extend 照收）。
- **参考锚定 LB/RB**：参考去 gap 逐位配对，平滑"香肠"扩展后每行只在 `[lb[i], rb[i]]` 内计算
  （变宽带）。
- **边缘裁剪**：只删两端"全物种都是 gap"的列（旧"单侧列"规则会误删真实内容）。

### 3.5 打分参数

- 替换矩阵：`chain::SubMatrix`，DP 中除以 50 缩放（默认 `hoxd55`，可 `--score-scheme` 覆盖）。
- gap：`--gap-model` 取 `GapCalc` 的 quasi-natural 曲线（`medium`/`loose`） 反推仿射参数，
  `constant` 直接用 `gap_score`；也可 `--align-gap-open` /`--align-gap-extend` 显式指定。所有值按
  `match_score/100` 缩放。
- `--mismatch-score` 目前只存在于配置中，实际 base–base 打分由替换矩阵决定。

## 4. 实现演进记录

> 以下按时间记录关键决策；当前行为以 §3 为准。

### 4.1 阶段 1：全体物种 SP 打分 + 确定性合并次序（2026-08-03）

- **打分升级**：从"只对两个 block 的共同物种打分"改为"共享物种自对 + 参考交叉对" （星型拓扑）。
  两个实证发现（S288c 三输入，参考去 gap 3826 bp）：
    - 朴素全笛卡尔积会拖偏参考锚点：真实数据上产生 155 列错位（3981 bp vs 3826），因 base–gap
      逐对罚分远大于错配罚分，DP 倾向移位避开 gap 列。
    - 对角线上 gap 贡献改为 0，gap 成本全部由转移承担，避免与替换分数失衡。
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

## 附录：上游链路 LASTZ 与链化

UCSC 典型 WGA 流程中，multiz 位于 "pairwise 比对 + 链化 + net + mafFromNet" 之后，只消费整理好的
MAF。`pgr` 的前置封装：

- `pgr align lastz`：LASTZ 前端，内置 UCSC preset（`set01`..`set07`）、自动加
  `--format=lav --markend --ambiguous=iupac` 等选项，支持目录递归与 rayon 并行。
- `pgr psl chain`：axtChain 风格链化，`SubMatrix` 打分 + `GapCalc` gap 曲线
  （`--gap-model loose|medium` 或显式仿射），KD-tree 加速前驱搜索。

两者覆盖 "blastz/lastz 比对 + axtChain 链化" 两步，为 fas-multiz 提供 pairwise 对齐基础。
