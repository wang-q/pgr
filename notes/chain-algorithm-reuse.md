# pgr chain 模块算法的可复用场景分析

> 本文档梳理 `src/libs/chain` 中各通用算法，分析它们在 pgr 其他模块中的潜在复用场景，为后续能力下沉/重构提供参考。
> 本文只输出分析结论，不立即修改代码。

## 1. 引言

`src/libs/chain` 实现了 UCSC Chain 格式处理的一整套算法：从 PSL/Chain 读入、KD-tree 动态规划链式化、Gap 代价计算、重复/退化过滤，到 Net 构建与分类。这些算法中相当一部分与“生物序列”语义解耦，可以泛化为通用计算模式。

本文的目标是把它们抽象成“问题 + 当前实现 + 可复用模块/文件 + 优先级”，方便后续在合适时机把能力下沉到 `libs/` 的通用位置，减少重复实现。

## 2. chain 模块通用算法清单

下面 8 个算法被认为具有跨模块复用价值。具体实现细节见 `notes/chain-algorithms.md`，本文只描述其通用抽象。

| 算法 | 位置 | 通用问题抽象 | 当前依赖 |
|------|------|--------------|----------|
| KD-tree 前驱搜索 | `src/libs/ds/kdtree.rs` | 二维区间带权最长路径：为每个矩形找左下方的最优前驱 | `KdTreeItem` trait |
| 分段插值 Gap 代价 | `src/libs/ds/gap_calc.rs` | 一维标量函数：查表 + 线性插值 + 外推 | 无（纯数值） |
| 最优重叠剪切 | `src/libs/chain/connect.rs::find_crossover` | 两个重叠区间上按评分函数选最佳 cut 点 | `SubMatrix` + `SequenceReader` |
| Top-K 成分纯度检测 | `src/libs/chain/anti_repeat.rs` | 多类别计数中前 K 类占比过高时惩罚 | 无（纯计数） |
| 位图范围覆盖 | `src/libs/ds/bitmap.rs` | `[0, size)` 上布尔范围集合：O(⌈L/64⌉) 置位/全置检查 | 无 |
| 层次化区间填充 | `src/libs/chain/net/builder.rs::fill_space` | 按优先级在剩余空间中插入区间并分裂为嵌套结构 | `BTreeMap` + `Rc<RefCell>` |
| 区间深度统计 | `src/libs/ds/dupe_tree.rs::DupeTree` | 带符号区间覆盖深度统计 | 无 |
| 排序/拆分/合并模式 | `src/libs/chain/sort.rs` / `stitch.rs` / `src/cmd_pgr/chain/split.rs` | 按 score 排序、按序列名分桶、同 ID 片段合并 | `Chain` 结构 |

## 3. 模块级复用机会

### 3.1 PAF 模块 (`src/libs/paf/`)

#### 3.1.1 `src/libs/paf/query.rs`：syntenic filter 可用 `BitMap`/`DupeTree` 替换线性扫描

当前 `query.rs` 加载 UCSC chain 文件后，用 `HashMap<(String, String), Vec<(u64, u64)>>` 保存每对 `(t_name, q_name)` 的 query 区间，过滤时做线性 `any`：

```rust
spans.iter().any(|&(cs, ce)| qs < ce && qe > cs)
```

当 chain 文件很大（全基因组比对）时，每个 region 都要重复扫描大量区间。可改为：

- 若已知 q_chrom 长度，为每对 `(t, q)` 建一个 `BitMap`，`set_range(q_start, q_end - q_start)`，查询 `is_fully_set` 或判断是否有覆盖。
- 若需要“覆盖比例”而不仅是布尔，可用 `DupeTree` 的 `count_over(..., 1)` 统计重叠长度。

收益：把 O(N_regions × N_spans) 降到 O(N_regions × ⌈L/64⌉)。

#### 3.1.2 `src/libs/paf/graph/segment.rs::split_alignment`：切分后可链式化

PAF 按 CIGAR 在 indel 处切分后得到一组 `Segment`（每个 segment 有 t_start/t_end/q_start/q_end）。这些 segment 天然满足 `ChainableBlock` 的语义。未来若 PAF 需要“链式化”低质量/碎片化比对，可直接复用 `connect.rs::chain_blocks` 的 KD-tree DP。

收益：避免重新实现一套 DP + 剪枝。

#### 3.1.3 `src/libs/paf/graph/builder.rs`：区间管理可借鉴 `fill_space`

GFA 图构建时需要管理 segment 在染色体上的布局、合并相邻 gap、用 DSU 做传递闭包。当前实现用排序 + `HashSet` 去重。若未来需要处理嵌套区间或按优先级放置 segment，可直接复用 `ChainNet` 的 `fill_space` 模式（`BTreeMap` 空间分裂）。

收益：统一嵌套区间插入 API，减少手工维护空间索引。

### 3.2 FAS / multiz 模块 (`src/libs/fas_multiz/` + `src/cmd_pgr/fas/`)

#### 3.2.1 `src/libs/fas_multiz/banded_align.rs`：已复用 `GapCalc` + `SubMatrix`（最佳案例）

该文件已经把 `chain::GapCalc` 和 `chain::SubMatrix` 拿来将 medium/loose gap 模型转换为 affine open/extend：

```rust
use crate::libs::chain::GapCalc;
use crate::libs::chain::sub_matrix::SubMatrix;
```

这是 chain 算法复用的正面案例。后续若要把 `GapCalc` 进一步下沉为通用 1D 插值模块，这里不需要改动。

#### 3.2.2 `src/libs/fas_multiz/merge.rs`：冲突合并可借鉴 `find_crossover`

成对 `FasBlock` 合并时，若参考序列在重叠区不完全相同，当前实现直接放弃合并。可改为在重叠区调用 `connect.rs::find_crossover` 的思想：用替换矩阵对重叠碱基打分，选择总分最大的 cut 点，把参考序列切开后再合并。

收益：提高合并成功率，减少因局部冲突而丢弃整个窗口的情况。

#### 3.2.3 `src/libs/fas_multiz/windows.rs`：覆盖深度可用 `DupeTree`/`BitMap`

当前 `derive_windows_from_blocks` 手动排序、合并区间，然后逐个判断每个输入是否覆盖窗口。`core` 模式要求所有输入都覆盖，`union` 模式要求至少一个输入覆盖。这些判断本质上是区间覆盖深度统计：

- 用 `DupeTree` 收集所有输入的参考区间，`build()` 后 `count_over(window.start, window.end, required_inputs)` 即可判断是否满足 core/union 条件。
- 若只关心布尔覆盖且染色体很大，可用 `BitMap`。

收益：统一覆盖深度 API，减少手写扫描逻辑。

#### 3.2.4 `src/cmd_pgr/fas/cover.rs`：大规模覆盖可用 `BitMap` 替代 `IntSpan`

该命令用 `intspan::IntSpan` 聚合多个 block 在参考序列上的覆盖。对于单条染色体很大的场景（哺乳动物基因组），`BitMap` 的位向量布局比 `IntSpan` 的区间列表更紧凑，且 `is_fully_set`/`set_range` 是 O(⌈L/64⌉)。`benches/bitmap_intspan_benchmark.rs` 的实测结果：100M 染色体、约 10% 覆盖时，`BitMap` 内存约为 `IntSpan` 的 1/3（12.5 MB vs 38 MB），完全覆盖查询速度提升约 1500 倍，混合构建+查询速度提升约 530 倍；构建速度在小染色体上 `IntSpan` 略快，但在 100M 规模两者接近。

收益：显著降低内存占用，并在大型染色体上大幅提升覆盖查询与聚合速度。

### 3.3 Alignment 模块 (`src/libs/alignment/`)

#### 3.3.1 `src/libs/alignment/trim.rs`：indel 区间运算可用 `DupeTree` 统一

`trim.rs` 用 `IntSpan` 处理 indel 区间，做 union/intersect。`DupeTree` 的“区间加减”能力与其等价，且额外提供深度统计。若后续需要按深度过滤 indel，可直接复用 `DupeTree`。

#### 3.3.2 `src/libs/alignment/stat.rs`：MSA 列质量可借鉴 Top-K 纯度检测

`anti_repeat.rs` 的退化度检测统计匹配碱基中 ACGT 分布，前两类占比 > 80% 则降权。这个思想可直接用于 MSA 列过滤：若某列被少数碱基主导，说明该列保守性差或存在系统误差，可降权或标记。

收益：为 MSA 质量评估增加一个通用指标。

#### 3.3.3 `src/libs/alignment/msa.rs`：列过滤可借鉴退化度检测

与 `stat.rs` 类似，MSA 输出前的列过滤可以引入 Top-K 纯度惩罚，提高输出序列质量。

### 3.4 POA 模块 (`src/libs/poa/`)

#### 3.4.1 `src/libs/poa/poa.rs`：序列加入顺序可借鉴 `chain_blocks`  peeling

当前 POA 按输入顺序依次 `add_sequence`。若输入顺序有偏（例如先加入低质量序列），图结构会偏向这些序列。可借鉴 `chain_blocks` 的 peeling 策略：先对所有序列做一次全局打分/排序，按 score 降序选择“骨架序列”优先加入，再补充其他序列。

收益：提高 MSA/共识序列的稳定性。

注意：这属于理念借鉴，需要较多改造，优先级中低。

### 3.5 Clust 模块 (`src/libs/clust/`)

#### 3.5.1 `src/libs/clust/hier.rs`：代表序列选择可借鉴 DP peeling

层次聚类完成后，需要为每个簇挑选代表序列。`chain_blocks` 的 peeling 策略（按 score 降序提取不重叠/最优链）与“从簇中提取代表性序列”目标一致。

收益：统一“按分数提取代表”的算法模式。

#### 3.5.2 `src/libs/clust/tree_cut/dynamic.rs`：区间扫描可参考 chain 的区间操作

该模块处理高度序列的分割/合并，本质上是一维区间扫描。虽然 `BitMap`/`DupeTree` 不完全适用（因为高度是连续值而非布尔），但“带约束的区间合并”思想可参考。

### 3.6 Net / PSL 模块

#### 3.6.1 `src/cmd_pgr/net/class.rs` 与 `filter.rs`：已消费 chain 库

这两个命令直接调用 `libs/chain::net` 的 `read_nets`、`collect_stats_gap`、`filter_chrom`、`prune_gap` 等函数。说明 `libs/chain::net` 本身已经被视为通用 Net 处理库。

未来若要支持其他层级区间结构（例如 GFF 的 exon-intron 树、VCF 的 SV 嵌套区间），可直接套用 `FilterCriteria` + 树遍历模式。

#### 3.6.2 `src/libs/chain/psl_chain.rs`：分组逻辑可泛化为通用转换器

`group_psl_blocks` 按 `(target_name, query_name, query_strand)` 分桶，并把每个 PSL block 转换为 `ChainableBlock`。这个“成对比对格式 → 链式块”的转换逻辑可抽象为通用函数，供 PAF 等格式复用。

收益：PAF 若要做链式化，可直接调用该抽象，避免重复分桶/块构建代码。

### 3.7 Pipeline 模块 (`src/libs/pl/` + `src/cmd_pgr/pl/ucsc.rs`)

#### 3.7.1 `src/cmd_pgr/pl/ucsc.rs`：子命令调用可下沉为库函数调用

UCSC pipeline 已大量调用 chain 子命令（`chain sort`、`chain pre-net`、`chain net`、`chain stitch` 等）。当前通过外部进程调用，存在序列化/反序列化开销。未来可评估把部分调用下沉为直接调用 `libs/chain` 函数：

- `libs/chain::sort::sort_chains`
- `libs/chain::pre_net::pre_net`
- `libs/chain::net::builder::ChainNet`
- `libs/chain::stitch::stitch_chains`

收益：减少外部进程开销，提高 pipeline 吞吐量。

注意：工作量大，属于长期重构，优先级低。

## 4. 推荐优先级

优先级以“能否减少重复代码”和“是否已有明显痛点”为主要标准。

| 优先级 | 算法/能力 | 目标模块/文件 | 理由 |
|--------|-----------|---------------|------|
| 高 | `GapCalc` 分段插值 | 已复用于 `fas_multiz/banded_align.rs`；建议下沉为通用 1D 插值 | 复用价值已被验证，抽象成本低 |
| 高 | `BitMap` 范围覆盖 | `libs/paf/query.rs` syntenic filter | 线性扫描是明显痛点，位图可显著提速 |
| 高 | `DupeTree` 区间深度 | `libs/fas_multiz/windows.rs` | 手写覆盖判断重复且易错，统一 API 收益大 |
| 中 | `find_crossover` 最优剪切 | `libs/fas_multiz/merge.rs` | 能提高合并率，但需引入序列读取和替换矩阵 |
| 中 | KD-tree 前驱搜索 | `libs/paf/graph/segment.rs`（未来链式化） | 通用性强，但 PAF 当前未明确需要链式化 |
| 中 | Top-K 纯度检测 | `libs/alignment/stat.rs` / `msa.rs` | 思想通用，但需与现有质量指标整合 |
| 低 | `ChainNet` 层次填充 | `libs/paf/graph/builder.rs` | 理念匹配，但改造范围大 |
| 低 | sort/stitch/split 模式 | `libs/paf/to_bed.rs` / `to_maf.rs` | 收益有限，当前已有排序实现 |
| 低 | Pipeline 子命令下沉 | `cmd_pgr/pl/ucsc.rs` | 长期工程，需保持 CLI 与库行为一致 |

## 5. 后续可执行动作

若后续决定真正复用，可按以下顺序推进：

1. **创建通用 1D 插值模块**
   - 将 `src/libs/chain/gap_calc.rs` 中的插值逻辑提取为 `src/libs/ds/gap_calc.rs`。
   - `GapCalc` 保持链式 gap 代价语义，同时可被其他模块直接复用。
   - 验证 `fas_multiz/banded_align.rs` 仍正常工作。

2. **创建通用位图模块**
   - 将 `src/libs/chain/bitmap.rs` 提升为 `src/libs/ds/bitmap.rs`。
   - 在 `src/libs/paf/query.rs` syntenic filter 中试用，替换 `HashMap<(t,q), Vec<(u64,u64)>>`。
   - 补充基准测试，对比 `IntSpan` 在密集场景下的内存/速度。

3. **创建通用区间深度模块**
   - 将 `src/libs/chain/net/syntenic.rs::DupeTree` 提升为 `src/libs/ds/dupe_tree.rs`。
   - 替换 `src/libs/fas_multiz/windows.rs` 中的手动覆盖判断。
   - 考虑扩展为按深度分层输出区间的 API。

4. **最优剪切泛化**
   - 将 `src/libs/chain/connect.rs::find_crossover` 的核心“重叠区间最佳 cut”逻辑提取为通用函数。
   - 在 `src/libs/fas_multiz/merge.rs` 中试点，处理参考序列冲突合并。

5. **KD-tree 泛化**
   - 将 `src/libs/chain/kdtree.rs` 的 `ChainItem` trait 重命名为更通用的 `RectItem` 或 `IntervalItem2D`。
   - 评估是否在 PAF 链式化或 POA 序列排序中使用。

6. **文档与测试**
   - 为下沉后的通用模块补充单元测试。
   - 更新 `notes/chain-algorithms.md` 与 `notes/chain-algorithm-reuse.md`，说明能力已下沉到新位置。

## 6. 注意事项

- `BitMap`、`GapCalc`、`KdTree`、`DupeTree` 已下沉到 `src/libs/ds`，可由其他模块通过 `crate::libs::ds::*` 使用。`chain` 模块仍通过 `pub use` 保持向后兼容。
- 下沉时应避免破坏 chain 模块的现有 API，优先采用“新通用模块 + chain 端薄封装”的方式。
- 任何下沉都应伴随单元测试和至少一个实际使用场景验证，避免为抽象而抽象。
