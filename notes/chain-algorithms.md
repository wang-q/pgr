# pgr chain 模块运行逻辑

> 本文件记录 `src/libs/chain` 与 `src/cmd_pgr/chain` 中各算法的实际运行流程，便于后续复用/修改时快速理解。

## 1. 总体流程

```
PSL / chain 文件
    │
    ├─ psl chain ──→ 分组 ──→ KD-tree DP ──→ Chain 文件
    │
    ├─ chain sort ──→ 按 score 降序
    │
    ├─ chain split ──→ 按 target/query 序列名拆文件
    │
    ├─ chain stitch ──→ 同 ID 片段合并
    │
    ├─ chain anti-repeat ──→ 过滤低复杂度/重复区
    │
    ├─ chain pre-net ──→ 预过滤（标记已覆盖区域）
    │
    └─ chain net ──→ 构建 Gap/Fill 树（Net）
```

## 2. Chain 数据结构

见 `src/libs/chain/record.rs`。

- `ChainHeader`：score, t_name/t_size/t_strand/t_start/t_end, q_name/q_size/q_strand/q_start/q_end, id。
  - `t_strand` 始终为 `'+'`。
  - `q_strand` 可为 `'+'` 或 `'-'`。
  - 坐标为 0-based、半开区间 `[start, end)`。
- `ChainData`：每个 block 的 `size`，以及到下一个 block 的 `dt`（target 侧 gap）和 `dq`（query 侧 gap）。
- `Chain` 可通过 `to_blocks()` 把相对坐标展开为绝对 `Block`（t_start/t_end/q_start/q_end）。

## 3. PSL → Chain（`src/libs/chain/psl_chain.rs` + `connect.rs`）

入口：`chain_psl(reader, writer, gap_calc, min_score, score_context)`。

### 3.1 分组

- 按 `(target_name, query_name, query_strand)` 分组。
- 对每条 PSL 的每个 block，生成一个 `ChainableBlock`。
- 若提供了 `ScoreContext`（含 2bit 序列 + 替换矩阵），则调用 `calc_block_score` 按实际序列重新计算 block score；否则用 `size * 100.0` 的近似值。
- 忽略 target strand 为负的 PSL 记录。

### 3.2 KD-tree 动态规划（`chain_blocks`）

目标：把同组 block 连成一条或多条链，使总分最大。

#### 3.2.1 DP 定义

对每个 block `i`：

- `dp_entries[i].total_score`：以该 block 结尾的最佳链总分。
- `dp_entries[i].best_pred`：最佳前驱 block 的索引（`None` 表示无）。

初始化：`total_score = block.score`。

#### 3.2.2 KD-tree 前驱搜索（通用算法）

见 `src/libs/ds/kdtree.rs`。这是一个**二维区间上的带权最长路径前驱搜索**，不依赖于生物序列语义，可复用到任何需要"找左下侧最优前驱"的场景。

**问题抽象**：给定平面上带权矩形（`[t_start,t_end) × [q_start,q_end)`，权重为 `score`），对每个矩形找另一个完全位于其左下方（`t_end <= target.t_start` 且 `q_end <= target.q_start`）的矩形，使 `pred_total_score + target.score - cost` 最大。

**数据结构**：
- 2D KD-tree，维度交替为 `q_start` 和 `t_start`。
- 每个节点维护子树内最大的 `q_end`、`t_end` 和 `max_score`（当前已知的最佳链尾总分）。

**核心操作**：
1. `build`：按中位数递归分割，O(n log n)。
2. `update_scores`：DP 更新某个 item 的总分后，沿途更新节点的 `max_score`，O(log n)。
3. `best_predecessor`：查询当前 item 的最佳前驱，带两类剪枝：
   - 分数上界剪枝：`node_max_score + target_score <= best_score`。
   - 几何距离下界剪枝：`node_max_score + target_score - lower_bound_cost(dq, dt) <= best_score`。
   - 搜索顺序：先搜包含目标坐标的近子树，再判断远子树是否可能更优。

**通用复用点**：
- 把 `ChainItem` trait 实现到任意带 `(x_start, x_end, y_start, y_end, score)` 的结构即可复用。
- `cost_func` 和 `lower_bound_func` 用闭包传入，可适配不同代价模型。

#### 3.2.3 转移函数

对候选前驱 `cand` 和当前 block `target`：

1. 必须满足 `cand.t_start <= target.t_start` 且 `cand.q_start <= target.q_start`（单调性）。
2. 计算 `dt = target.t_start - cand.t_end`，`dq = target.q_start - cand.q_end`。
3. 若 `dt < 0` 或 `dq < 0`，说明存在重叠：
   - `overlap_len = max(|dt|, |dq|)`。
   - `density = max(cand.score/len_cand, target.score/len_target)`。
   - `overlap_penalty = overlap_len * density`。
4. `cost = gap_calc.calc(dq, dt)`。
5. 新总分 = `cand.total_score + target.score - cost - overlap_penalty`。

注意：即使候选与当前 block 坐标有重叠，也会尝试连接并用 penalty 惩罚，不是直接丢弃。

#### 3.2.4 提取链（peeling）

- 按 `total_score` 降序排列 block。
- 从最高分未命中 block 开始回溯前驱，形成一条链，并标记这些 block 为 `hit`。
- 对链内 block 做后处理：
  - `remove_exact_overlaps`：删除坐标完全重复的 block，保留第一个，分数不累加。
  - `merge_abutting_blocks`：合并首尾相接的 block，分数相加。
  - 若有 `ScoreContext`，调用 `trim_overlaps` 用实际序列找最佳 cut 点，把相邻 block 的重叠部分切开。
- 用 `score_chain` 重新计算链总分（精确序列打分或近似打分）。
- 若 `score <= 0` 则丢弃。
- 生成 `Chain` 头信息并分配递增 ID。

### 3.3 分段插值 Gap 代价（通用算法）

见 `src/libs/ds/gap_calc.rs`。这是一个**一维标量函数的分段查表 + 线性插值 + 外推**实现，可复用到任何需要"根据距离查非线性代价"的场景。

**问题抽象**：给定一组控制点 `(pos, value)`，实现函数 `f(x)`：
- `x < small_size`：查预计算整数表。
- `small_size <= x <= max_pos`：在控制点间线性插值。
- `x > max_pos`：用最后一段斜率外推。

**当前用法**：
- 维护 query gap、target gap、同时 gap 三类代价表。
- 内建 `medium()`（mouse/human）和 `loose()`（distant species）。
- 也支持 `affine(open, extend)`，此时代价 = `open + extend * len`。
- `calc(dq, dt)`：
  - `dt == 0`：查 query gap 表。
  - `dq == 0`：查 target gap 表（通常与 query 相同）。
  - 否则：用 `max(dq, dt)` 查 both gap 表。
- 负距离会被截断为 0。

**通用复用点**：
- 把 `pos` 和 `vals` 换成其他物理量的控制点，即可得到任意分段线性函数。
- 预计算小表加快热路径访问；插值/外推保证大值也有定义。

### 3.4 带替换矩阵的最优重叠剪切（通用算法）

见 `src/libs/chain/connect.rs` 中的 `trim_overlaps` / `find_crossover`。这是一个**在两个重叠区间上找最佳分界点**的线性扫描算法。

**问题抽象**：两个对象 A、B 在一段长度为 L 的区域上重叠。对 cut 位置 `i ∈ [0, L]`，A 保留前 `i` 部分、B 保留后 `L-i` 部分。定义评分函数 `score_A(i) + score_B(L-i)`，求使总分最大的 `i`。

**核心步骤**：
1. 读取 A 的右端重叠区与 B 的左端重叠区序列。
2. 预处理 B 区总分 `r_score = sum(score_B[k])`。
3. 从左到右扫描 `i = 0..=L`：
   - `current_l = sum(score_A[0..i])`。
   - `current_r = r_score - sum(score_B[0..i])`。
   - 记录 `current_l + current_r` 的最大值及对应 `i`。
4. 按最优 `i` 调整 A/B 边界。

**当前用法**：用 DNA 替换矩阵给每对碱基打分；query 负链需先做反向互补。

**通用复用点**：
- 任意有两个重叠片段、需要按某种打分函数选择最佳 cut 的场景都可复用。
- 时间 O(L)，空间 O(L)（用于 char 数组）。

## 4. chain sort

- 读入所有 chain，按 `header.score` 降序排序。
- 默认重新编号 ID（从 1 开始）；`--save-id` 保留原 ID。
- 注意：不接受 stdin，必须通过文件或 `--input-list` 提供。

## 5. chain split

- 默认按 `t_name` 拆分，可用 `--by-query` 改为 `q_name`。
- 输出目录中的文件名为 `<seq>.chain`。
- `--lump N`：
  - 扫描序列名中第一个连续数字段，取 `val % N` 作为桶名（3 位零填充）。
  - 无数字时退化为 `fxhash64(name) % N`。
  - 实际输出文件数可能小于 N，取决于输入中不同桶的数量。
- 会检查序列名不含路径分隔符且不以 `.` 开头，防止路径穿越。

## 6. chain stitch

- 按 `chain.header.id` 分组。
- 同组内要求 `t_name`、`q_name`、`q_strand` 一致，不一致则跳过并 warning。
- 把所有片段的 block 合并后按 `(t_start, q_start)` 排序。
- 调用 `Chain::from_blocks` 重建 data，header 范围自动更新。
- 分数累加所有片段分数。
- 输出按 score 降序排列。
- 不检查片段间是否重叠。

## 7. chain anti-repeat

- 输入为 chain 文件 + target/query 2bit 文件。
- 对每条 chain：
  - score >= `--no-check-score`（默认 200000）直接保留。
  - score < `--min-score`（默认 5000）直接丢弃。
  - 否则执行两个检查：
    1. **Degeneracy 检查**：统计匹配碱基中 ACGT 的分布。若前两种碱基占比 > 80%，按超出比例降低 score；低于 `min_score` 则丢弃。
    2. **Repeat 检查**：统计 lowercase（软屏蔽）碱基比例。`adjusted_score = score * 2 * (total - rep) / total`，低于 `min_score` 则丢弃。

### 7.1 Top-K 成分纯度检测（通用算法）

Degeneracy 检查本质是一个**Top-K 类别占比惩罚**算法。

见 `src/libs/ds/top_k_purity.rs`。

**问题抽象**：给定一个多类别计数向量 `counts[0..C-1]`，总次数 `total`。若最大的 K 项之和占比超过阈值 `ok_ratio`，则按超出比例对分数进行惩罚。

**核心步骤**：
1. 统计每个类别的出现次数。
2. 取最大的 K 项之和 `best_k`。
3. `observed = best_k / total`。
4. 若 `observed <= ok_ratio`，不惩罚。
5. 否则 `adjust_factor = 1.01 - (observed - ok_ratio) / (1 - ok_ratio)`，用 `score * adjust_factor` 作为调整后分数。

**当前用法**：`C = 4`（T/C/A/G），`K = 2`，`ok_ratio = 0.80`，只统计匹配位置（target == query）。

**通用复用点**：
- 任意需要检测序列/样本是否被少数类别主导的场景。
- 可调整 K、类别数、阈值和惩罚曲线。

**anti-repeat 其他实现细节**：
- 读取序列时保留大小写（`include soft masks`）。
- 负链 query 需要反向互补后比较；使用独立的 `nt_val` 映射（T=0, C=1, A=2, G=3）以便用 `(v + 2) % 4` 求互补。

## 8. chain pre-net 与位图范围覆盖（通用算法）

见 `src/libs/ds/bitmap.rs` 与 `src/libs/chain/pre_net.rs`。`BitMap` 是一个**紧凑的位向量范围集合**，支持 O(⌈L/64⌉) 的范围设置与全置检查。

**问题抽象**：在 `[0, size)` 上维护一个布尔集合，支持：
- `set_range(start, len)`：把区间全部置 1。
- `is_fully_set(start, len)`：判断区间是否全为 1。

**实现要点**：
- 用 `Vec<u64>` 存储位，每个 word 64 位。
- 设置时按首/尾 word 的边界构造掩码；中间 word 直接写 `!0u64`。
- 查询时同理检查掩码是否全部被置位。
- 所有边界计算使用 `saturating_add` 并 clamp 到 `size`，防止溢出。

**当前用法**：
- pre-net 为每条序列建一个 `BitMap`。
- 按 score 降序处理 chain；若 chain 的某个 block 还有未覆盖区域则保留，并把 block 扩展 `pad` 后标记为已覆盖；否则丢弃。
- 缺失 sizes 的序列会报错。

**通用复用点**：
- 任何需要"按优先级首次覆盖并去重"的场景都可用此模式（例如区间去重、已访问区域过滤）。
- 若需要计数而不仅仅是布尔，可扩展为分块计数器。

## 9. chain net（最复杂）

入口：`src/cmd_pgr/chain/net.rs` 调用 `src/libs/chain/net/` 下模块。

### 9.1 数据结构（`types.rs`）

- `Chrom`：代表一条 target/query 染色体，含 `root`（根 Gap）、`spaces`（可搜索空间索引）。
- `Gap`：未对齐区间，可包含多个 `Fill`。
- `Fill`：对齐区间，引用来源 chain，可包含多个子 `Gap`。
- `Space`：`Chrom` 上用于定位 Gap 的索引项， key 为 `start`，value 指向所属 Gap。

初始化时 `Chrom.root` 覆盖 `[0, size)`，对应一个 `Space(0, size, root)`。

### 9.2 构建 ChainNet（`builder.rs`）

`ChainNet::new(sizes)` 为每条染色体创建空的 net 树。

#### 9.2.1 插入 chain（target 视角）

`add_chain(chain, min_space, min_fill, min_score)`：

1. 过滤 score < min_score 的 chain。
2. 取出 chain 所在 target 染色体的 `Chrom`。
3. 调用 `add_chain_core(chrom, chain, blocks, is_q=false, ...)`。

#### 9.2.2 插入 chain（query 视角）

`add_chain_as_q(...)`：

1. 把 chain 的 block 按 query 坐标翻转（负链需先 reverse）。
2. 以 `q_name` 为染色体键插入，相当于把 query 当作 target 来建 net。

#### 9.2.3 add_chain_core 流程

对 chain 覆盖的每个连续 space：

1. 在 chain 的 block 列表中找到与该 space 相交的所有 block，得到 `[first_idx, last_idx]`。
2. 计算 fill 的 `[fill_start, fill_end]` 为这些 block 与 space 交集的并集范围。
3. 若 fill 长度 < min_fill，跳过。
4. 调用 `fill_space(...)` 在 space 中创建 Fill。

### 9.3 层次化区间填充 fill_space（通用算法）

见 `src/libs/chain/net/builder.rs`。这是一个**按优先级将区间插入到剩余空间并生成嵌套结构**的算法。

**问题抽象**：染色体上有一组互不重叠的初始空间（初始为 `[0, size)`）。按优先级顺序处理若干"请求区间"，每个请求会：
1. 找到当前所有与其重叠的空间。
2. 对每个空间，用请求的交集填充一个"Fill"，并把原空间分裂为左右剩余空间（若长度足够）。
3. 若请求内部还有子间隙，则在这些间隙处创建子空间，供后续更低优先级的请求填充。

**核心数据结构**：
- `Chrom.spaces`：`BTreeMap<u64, Space>`，key 为 space 起点，用于快速查找与某区间重叠的所有 space。
- `find_spaces(start, end)`：扫描 `range(..end)` 中 `space.end > start` 的项。

**fill_space 步骤**：
1. 从 `chrom.spaces` 删除原 space。
2. 计算 Fill 在另一侧（query/target）的坐标 `o_start/o_end`，处理负链 reverse。
3. 创建 `Fill` 并挂到原 space 的父 Gap 上。
4. 若 fill 左侧/右侧还有足够长度 >= min_space，则生成新的 space 并插入索引。
5. 对相邻 block 之间的内部 gap：
   - 若 gap 长度 >= min_space，创建新的 `Gap`。
   - 把该 gap 作为 space 插入 `chrom.spaces`，并挂到当前 Fill 的 `gaps` 下。

**当前用法**：chain 按 score 降序插入，高分 chain 先占据顶层空间，低分 chain 只能在剩余 gap 中形成嵌套 Fill。

**通用复用点**：
- 任何"按优先级在染色体/序列上放置区间并分裂剩余空间"的问题（例如基因组注释叠加、区间布局、嵌套区间树）。
- 把 `Fill`/`Gap` 换成领域对象即可。

### 9.4 最终化（`finalize.rs`）

所有 chain 插入后：

1. `sort_net`：每层 Fill/Gap 按 `start` 排序，保证输出稳定。
2. `calc_other_fill`：根据实际 chain data 重新计算每个 Fill 的 `o_start/o_end`。
   - target 视角：遍历 chain data，找与 `[fill_start, fill_end)` 重叠的 block，取 query 方向的最小/最大坐标，负链 reverse。
   - query 视角：类似，但按 query 坐标 clip 后反推 target 坐标。

### 9.5 同线性分类（`syntenic.rs`）

目标：为每个 Fill 打标签 `top / syn / inv / nonSyn`，并计算 `qDup / qOver / qFar`。

#### 9.5.1 区间深度统计 DupeTree（通用算法）

见 `src/libs/ds/dupe_tree.rs`。这是一个**带符号区间覆盖深度统计**工具。

**问题抽象**：给定若干带符号区间 `[start, end)`，权重为 `+1`（增加覆盖）或 `-1`（减少覆盖），回答查询：区间 `[q_start, q_end)` 内覆盖深度 >= `threshold` 的总长度是多少？

**核心步骤**：
1. `add(start, end)` / `subtract(start, end)`：收集事件 `(start, +d)` 和 `(end, -d)`。
2. `build()`：
   - 按坐标排序事件；同坐标时，结束事件（负 delta）排在开始事件（正 delta）之前，保证相邻区间不互相渗透。
   - 扫描线生成 constant-depth `Segment` 列表。
3. `count_over(start, end, threshold)`：二分定位起始 segment，然后线性遍历到 `start >= end`，累加深度达标的重叠长度。

**当前用法**：
- Net 中每个 Fill 对其 query 范围 `+1`。
- 每个 Fill 内部的 Gap（子节点会占据的区域）`-1`，从而得到"真正重复覆盖"的深度。
- 用 `count_over(..., 2)` 计算 `qDup`。

**通用复用点**：
- 任意需要统计区间叠加深度（例如重叠区域、嵌套区间覆盖）的场景。
- 可扩展为按深度分层输出区间。

#### 9.5.2 classify_syntenic

- 先调用 `r_calc_dupes` 给所有 query 染色体建立 DupeTree。
- 再调用 `r_net_syn` 自顶向下遍历 net：
  - 无 parent 的 Fill 标记为 `top`。
  - 否则比较 parent Fill 的 `o_chrom` 和 `o_strand`：
    - 染色体不同 → `nonSyn`。
    - 染色体相同、链相同 → `syn`。
    - 染色体相同、链不同 → `inv`。
  - `qDup` = DupeTree 中深度 >= 2 的重叠碱基数。
  - `qOver` = 当前 Fill 与父 Gap 在 query 上的重叠长度。
  - `qFar` = 不重叠时到父 Gap 的距离，重叠时为 0。

### 9.6 输出（`writer.rs`）

`write_net(chrom, writer, is_q, min_score, min_fill)`：

- 对每个 Fill：
  - 若 `Fill.chain` 存在，调用 `subchain_info` 按当前 Fill 的 clip 范围重新计算 `sub_size`（ali）和 `sub_score`。
  - 否则用存储的 `ali/score`。
  - 仅当 `sub_score >= min_score` 且 `sub_size >= min_fill` 时输出。
- 输出顺序：先 `net <chrom> <size>`，然后递归输出 Fill/Gap，缩进递增。
- 字段顺序遵循 UCSC Net 格式：
  - fill：`fill tStart tLength qName qStrand qStart qLength id chainId score ali [qOver] [qFar] [qDup] type <class> [tN/qN/tR/qR/tTrf/qTrf]`。
  - gap：`gap tStart tLength qName qStrand qStart qLength [tN/qN/tR/qR/tTrf/qTrf]`。

### 9.7 Net 过滤（`filter.rs`）

`FilterCriteria` 支持按 score、ali、size、染色体名、 synteny 类型等过滤。

- `prune_gap`：递归遍历 Gap，对不通过 `filter_one` 的 Fill 直接丢弃；对不通过 `min_gap` 的 Gap 直接丢弃。
- `syn_filter`：判断是否属于 syntenic，依据 score、size、ali、`qFar`、class 等。
- `do_syn` 只保留 syntenic，`do_nonsyn` 只保留 non-syntenic。

## 10. 其他工具函数

- `range_intersection`：`[a,b)` 与 `[c,d)` 的交集长度。
- `lump_name`：从序列名提取第一个整数段做桶号，无整数时用 fxhash64。
- `is_haplotype`：序列名含 `_hap` 或 `_alt`。

## 11. 关键不变量/注意事项

1. Chain 坐标半开，`t_strand` 恒为 `'+'`。
2. `chain sort` / `chain net` / `chain pre-net` 都要求输入按 score 降序，否则会报错。
3. `GapCalc` 对负 gap 会截断为 0，DP 中重叠惩罚单独计算。
4. Net 构建时高 score chain 先插入并占据 space，低 score chain 只能在剩余 gap 中填充，因此输入顺序直接影响层次结构。
5. `anti-repeat` 中的 `nt_val` 映射是独立的，不要与 `crate::libs::nt::NT_VAL` 混用。
6. `stitch` 不检测片段重叠，依赖调用方保证。
7. `remove_exact_overlaps` 不会累加重复 block 的分数。
