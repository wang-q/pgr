# Spoa 源码分析（SIMD POA）

> 整理于 2026-08，源自对 `spoa-4.1.5/` 目录源码（约 5000 行 C++）及 README 的通读。
> 目的：理解 Spoa 的偏序比对（POA）算法、图结构与 SIMD 加速机制，为 pgr 的
> Rust 原生 POA（`libs/poa/`，`pgr fas consensus` / `pgr fas refine`）提供
> 参考与差异对照；移植与实现状态见 §8（原 `design/spoa_port.md` 已合并）。

## 1. Spoa 概览

- **工具定位**: POA（Partial Order Alignment）算法的 C++ 实现，用于从一组序列
  生成一致性序列（consensus）与多序列比对（MSA）。算法源自 Lee 2002
  （Bioinformatics 18:452，POA）与 Lee 2003（Bioinformatics 19:1236，consensus）。
- **版本**: 4.1.5（2024 年发布，4.x 系列）。
- **比对模式**（3 种）: 局部 Smith-Waterman（kSW）、全局 Needleman-Wunsch
  （kNW）、半全局 Overlap（kOV）。
- **gap 罚分模式**（3 种）: linear、affine、convex（分段 affine，双罚分函数取 min）。
- **SIMD 支持**: SSE4.1 与 AVX2（README 自评 "marginally faster due to high
  latency shifts"）；[SIMDe](https://github.com/simd-everywhere/simde) 提供
  非 x86 可移植；运行时分派按 AVX2 > SSE4.1 > SSE2。
- **交付形态**: 库（`libspoa`）+ 独立二进制 `spoa`（FASTA/FASTQ → consensus /
  MSA / GFA / DOT，stdin/stdout 友好）。
- **依赖**: bioparser/biosoup（序列 I/O）、zlib；可选 cereal（图序列化）、
  simde（可移植 SIMD）、cpu_features（x86 分派）、gtest（测试）。
- **线程模型**: 库与 CLI 均为**单线程**——并行化由使用者完成（spoa 本身只提供
  `Graph::Subgraph` 这类分治原语）。pgr 移植时的"单线程版本"即对应此现状。

## 2. 架构与数据流

```
FASTA/FASTQ（bioparser 流式读取）
  → 逐序列: AlignmentEngine::Align(seq, graph) → Alignment
  → Graph::AddAlignment(alignment, seq, [weight|quality])
  → 循环直至全部序列入图
  → Graph::GenerateConsensus / GenerateMultipleSequenceAlignment / GFA / DOT
```

- **Alignment** = `vector<pair<int32, int32>>`：`(图节点 id, 序列位置)` 对，
  -1 表示 gap（`first == -1` 为插入到图的新节点，`second == -1` 为序列被删）。
- **双层结构**：
  - `Graph`：偏序图，增量插入序列并维护拓扑序；
  - `AlignmentEngine`：抽象基类 + 具体引擎，把一条序列比对到图，返回
    Alignment（纯函数，不修改图）。

## 3. Graph：偏序图结构与增量插入

### 3.1 数据结构（`include/spoa/graph.hpp`，324 行）

- **Node**: `id`（图内序号）、`code`（碱基整数编码，0..num_codes-1，经
  `coder_`/`decoder_` 双向映射表与原始字符互转）、`inedges`/`outedges`、
  `aligned_nodes`。
- **Edge**: `tail`/`head`、`labels`（覆盖此边的序列 id 列表）、`weight`
  （int64，边覆盖权重累加）。
- **aligned_nodes（失配分支）**: POA 的核心机制——两条序列比对上同一位置但
  碱基不同时，两个节点通过 `aligned_nodes` 双向互连。失配节点与主节点共享
  MSA 列，同时保持 DAG（不引入边，无环）。

### 3.2 AddAlignment 插入逻辑（`src/graph.cpp`，693 行）

1. 空 alignment：整条序列作为新链线性插入（`AddSequence`）。
2. 有 alignment：
   - 未比对前缀/后缀（alignment 首尾之外的序列段）→ `AddSequence` 加链；
   - 对齐部分逐对处理：图节点 `code == 碱基` → 复用节点；不同 → 在
     `aligned_nodes` 中找同 code 节点，无则新建节点并与整个 clique（原节点 +
     其全部 aligned_nodes）双向互连；
   - 相邻对齐节点间 `AddEdge`：同 head 边合并（追加 label、累加 weight）。
3. **权重语义**: 全部落在边上——`AddSequence` 建边时权重为
   `weights[i-1] + weights[i]`（相邻两节点各贡献一次）；节点本身无 weight。
   支持统一 weight、逐位权重向量、Phred quality（`quality[i] - 33`）三种入参。
4. 每次插入后全图 `TopologicalSort()`（O(V+E)）。

### 3.3 拓扑排序

- 手写 DFS 栈式算法（非 Kahn），三态标记 + `ignored` 标志。
- **关键细节**: 节点入 rank 时，其 `aligned_nodes` 被"忽略"不单独入 rank，
  而是随主节点一起排在相邻位置——这保证 MSA 中失配节点共享同一列。
- `assert` 失败即图非 DAG（正常流程下不可能）。

### 3.4 Consensus：TraverseHeaviestBundle

- 逐 rank 节点取**最大权重入边**（tie 时比较前驱分数），累加得到到该节点的
  heaviest path 分数；全程取全局最大终点。
- 若终点仍有出边（路径可继续分支延伸），`BranchCompletion` 从该 rank 之后
  重算，继续追最大分数路径。
- 回溯得到 `consensus_` 节点序列；`GenerateConsensus(min_coverage)` 用
  `Node::Coverage()`（in+out 边 labels 的并集大小）过滤低覆盖节点。

### 3.5 MSA

- `InitializeMultipleSequenceAlignment`：rank 顺序 + aligned_nodes 共享列映射
  （每个节点一个列号，失配节点同列）。
- 每序列沿 `Successor(label)`（按边 labels 找到该序列对应的下一节点）走图，
  填碱基，其余位置补 `-`。

### 3.6 子图与坐标映射

- `ExtractSubgraph` / `Subgraph` / `UpdateAlignment`：提取两个节点之间的子图
  （含 aligned_nodes），并建立子图 ↔ 原图节点 id 映射，供局部重比对场景使用
  （主 CLI 未启用）。

## 4. AlignmentEngine：工厂、参数与标量 DP（SISD）

### 4.1 工厂与参数校验（`src/alignment_engine.cpp`，112 行）

- `Create(type, m, n, g[, e[, q, c]])`：gap open/extend 必须非正，否则抛异常。
- **subtype 判定**: `g >= e` → linear；`g <= q || e >= c` → affine；否则
  convex。linear 时折叠 `e = g`；affine 时折叠 `q = g, c = e`。
- `WorstCaseAlignmentScore` 预检潜在溢出（`kNegativeInfinity = i32::MIN + 1024`）。
- 工厂先尝试 `CreateSimdAlignmentEngine`，SIMD 不可用（无指令集编译）时回退
  `SisdAlignmentEngine`。

### 4.2 标量 DP（`src/sisd_alignment_engine.cpp`，928 行）

- **矩阵布局**: `(图节点数+1) × (序列长+1)`，按 `rank_to_node`（图拓扑序）逐
  行推进——DP 行 = 图节点，列 = 序列位置。
- **状态矩阵数**随 gap 模式增长：linear 1 个（H）、affine 3 个（H/F/E）、
  convex 5 个（H/F/E/O/Q，O/Q 为第二 affine 罚分函数）。
- **sequence_profile**: `alphabet × 序列长` 的 match/mismatch 查分表（预计算，
  逐碱基查表代替比对）。
- **递推（affine 为例，逐节点处理）**：
  - 图 gap（序列碱基插入图）: `F[j] = max(H_pred[j] + g, F_pred[j] + e)`
  - match/mismatch: `H[j] = H_pred[j-1] + profile[j]`
  - 序列 gap（图碱基被删）: `E[j] = max(H[j] + g, E[j-1] + e)`
  - 多入边（图分支）: 对每个前驱重复取 max
  - SW 每格 clamp 0；NW/OV 的行/列边界初始化不同（首列 gap 罚分 vs 0）
- **回溯**: SW 从最高分回溯到 0；NW 从 (末行, 末列) 到 (0,0)；OV 到行/列边界。
  输出 Alignment（图节点 id / 序列位置，gap = -1），最后逆序。

## 5. SIMD 引擎：垂直并行与分派

`src/simd_alignment_engine_implementation.hpp`（2065 行）+ 分派/调度
（dispatcher/dispatch，共 136 行）。

### 5.1 垂直并行（lane = 序列位置）

- `InstructionSet<A, T>` 模板抽象指令集，`T ∈ {i16, i32}`：
  - AVX2（256-bit）: i16 → 16 lane、i32 → 8 lane；
  - SSE4.1（128-bit）: i16 → 8 lane、i32 → 4 lane。
- **矩阵宽度压缩**为 `ceil(seq_len / kNumVar)` 个向量列：一个寄存器向量 =
  连续 kNumVar 个序列位置的 DP 值；图节点（行）循环在外层，每行处理若干
  向量列——即经典 striped SW（Farrar 2007）思想在 POA 上的变体。

### 5.2 向量内对角线搬移

- 对角线依赖 `H_pred[j-1]` 通过字节移位实现：`_mmxxx_slli_si(向量, kLSS)` 把
  每个 lane 左移一位（低位补 0），再与上一向量列提取出的尾元素（
  `_mmxxx_srli_si` 保存的 `x`）OR 拼接，随列推进滚动。

### 5.3 前缀最大值（gap 链）

- 行内 gap 延伸 `max_k(H[j-k] + k·g)` 是前缀最大值：`_mmxxx_prefix_max` 用
  `log kNumVar` 步的 mask + shift + add + max 完成（penalties 为 2 的幂倍 g），
  避免每列 O(kNumVar) 的串行扫描。README 所称"high latency shifts"即指这些
  字节移位/变量移位指令的延迟。

### 5.4 类型与架构选型

- **动态 lane 类型**: `WorstCaseAlignmentScore` 低于 i32 下界 → 全程 i32 lane；
  否则用 i16 lane（更宽、更快）。
- **编译期实例化**: `dispatch.cpp` 按 `__AVX2__` / `__SSE4_1__` / 默认 选
  `Architecture::kAVX2/kSSE4_1/kSSE2` 实例化模板。
- **运行时分派**: `dispatcher.cpp` 用 cpu_features（或手写 cpuid）检测
  AVX2 → SSE4.1 → SSE2；非 x86 目标用 SIMDe 翻译。
- 回溯与标量一致：向量列内定位最大分用 `_mmxxx_index_of` / `_mmxxx_value_at`
  （store 回标量数组后线性找）。

**与 pgr HV 分派形式的对照**：spoa 是"编译期三档模板（AVX2/SSE4.1/SSE2）×
运行期 cpuid 检测 + SIMDe 兜底非 x86"；pgr HV（`libs/hv.rs`，见
`benchmarks/bench-simd-hv-jaccard.md` §2）是更简的"单一 AVX2 手写 intrinsic +
`is_x86_feature_detected!` 运行时检测 + portable `wide` 回退"——不保留
SSE4.1 中间档、不用 SIMDe、AVX-512 仅留在基准对照。差异来源：spoa 的
SIMD 是 SSE 时代的遗产（当时无 wide 类可移植库），多档实例化是历史必需；
pgr 的 `wide` 回退自动映射到 NEON/SSE/标量，SSE4.1 中间档的额外实例化与
cpuid 复杂度收益小（README 自评 SIMD 增益 "marginal"）。**pgr 若未来为
POA 做 SIMD，分派应沿用 HV 式**（见 §7）。

## 6. 与 pgr Rust 移植的对照（`libs/poa/`，1602 行）

| 维度 | spoa-4.1.5 | pgr `libs/poa/` |
|---|---|---|
| 比对模式 | kSW / kNW / kOV | Local / Global / SemiGlobal ✓ |
| gap 模式 | linear / affine / convex | 仅 affine（gap_open + gap_extend） |
| DP 实现 | SISD + SIMD（SSE4.1/AVX2/SIMDe） | 标量 + SIMD（AVX2 手写 / `wide` 回退，2026-08-09） |
| 线程模型 | 库/CLI 单线程 | 单线程（pgr 侧 `--parallel` 按 block 并行） |
| 权重存储 | 全在边上（labels + weight） | 节点 weight + 边 weight |
| 权重入参 | 统一 / 逐位 / Phred quality | 统一（默认 1） |
| min_coverage 过滤 | ✓ | ✗ |
| strand-ambiguous | ✓ | ✗ |
| Subgraph / UpdateAlignment | ✓ | ✗ |
| GFA / DOT 输出 | ✓ | ✗（pgr 无此需求） |
| 图实现 | 自研 vector + 裸指针 | petgraph DiGraph + NodeData/EdgeData |
| 拓扑排序 | 手写 DFS 栈（aligned 处理） | 逐行移植同算法 |
| consensus | heaviest bundle（边权重） | heaviest bundle（边 + 节点权重） |

差异说明：

- **权重语义**: C++ 版权重全在边上（`weights[i-1] + weights[i]` 隐含节点贡献）；
  Rust 版额外把 coverage 存为节点 `weight` 并在 heaviest path 得分中计入
  （`consensus.rs`），使多数碱基胜过首序列骨架——输出一致性已由移植验证
  （§8 双引擎对照）。
- **SIMD 落地（2026-08-09）**: 移植初期只做标量（§8 实施阶段）；
  后补 `libs/poa/simd.rs`（垂直并行，分派沿用 HV 式——AVX2 手写 +
  `is_x86_feature_detected!` + `wide` 回退，§5.4 对照），`Poa` 默认引擎
  已切换。基准：120 bp ~8.7×、600 bp ~12.3×（`benches/poa_benchmark.rs`）。

## 7. 参考价值与可借鉴点

1. **SIMD 垂直并行方案**（**已落地 2026-08-09**，`libs/poa/simd.rs`）：
   lane = 序列位置、向量内移位取对角线、前缀最大值扫描（affine 的 E 状态
   一阶依赖，无需 i16/i32 动态选型，统一 i32 lane 与标量 `neg_inf` 语义
   对齐）。分派沿用 HV 式（§5.4 对照）——单一 AVX2 手写路径 +
   `is_x86_feature_detected!` + `wide` 可移植回退，不复制 spoa 的 SSE4.1
   中间档与 SIMDe。三路基准（Ryzen 9 7945HX，`benches/poa_benchmark.rs`）：
   120 bp 对齐 scalar 981 µs / wide 205 µs / avx2 124 µs；600 bp scalar
   23.5 ms / wide 3.85 ms / avx2 1.97 ms——wide ~4.8–6.1×、avx2 ~7.9–12×
   （相对标量），avx2 相对 wide 再快 ~1.7–2×。与 pgr HV 的 AVX2 经验
   （`bench-simd-hv-jaccard.md`）相互印证：移位指令延迟高是共同限速因素
   （spoa README 自评收益 "marginal"；pgr HV 的 i16 变量移位实验同样否决），
   故移植初期先做标量是合理决策。CLI 实测（2026-08-09，单线程，100 条 ×
   1 kb block FA）：`pgr fas consensus --engine builtin` ~1.9 s vs 外部
   `spoa` 4.1.4 ~7.0 s——Rust SIMD 反超外部 C++ SIMD 约 3.5×（输出逐字节
   一致）。
2. **Subgraph / UpdateAlignment 分治原语**：若未来需要长序列/大 block 并行
   或局部重比对，这是 spoa 提供的现成分解思路（pgr 移植未含）。
3. **quality 权重与 min_coverage**：Rust 版 `AddAlignment` 只支持统一权重；
   若未来处理 FASTQ 输入或需要低覆盖节点过滤，可参考 C++ 的 Phred 转换与
   `Node::Coverage` 语义补齐。
4. **参数语义对照**: subtype 判定规则（linear/affine/convex 折叠）与
   `WorstCaseAlignmentScore` 溢出预检是 pgr `ScalarAlignmentEngine` 尚未
   完整对齐的部分（Rust 版仅 affine 路径）。

## 8. 移植与实现状态（原 `design/spoa_port.md`，2026-08-09 合并）

- **目标**：将 Spoa 移植为 Rust 原生 POA（`libs/poa/`），集成到
  `pgr fas consensus` / `pgr fas refine`。保留外部 `spoa` 二进制作为引擎
  选项——移植初期外部 C++ SIMD 更快，后被内置反超（§7 CLI 实测 ~3.5×）。
- **架构决策**：图复用 `petgraph::DiGraph<NodeData, EdgeData>`（`NodeData` =
  碱基 + 失配分支 `aligned_nodes`，`EdgeData` = 边权重）；`AlignmentType`
  Global/Local/SemiGlobal + 仿射罚分；模块布局 = `mod/poa/graph/align/
  consensus/msa`（后加 `simd`）。
- **实施阶段**：① 图 + 标量 DP + 拓扑排序（对照 Spoa 验证）；② consensus/MSA
  与双引擎（`--engine builtin|spoa`、`--msa builtin|spoa`）集成，删除临时
  `cmd_pgr/poa/`；③ SIMD 垂直并行（§7，2026-08-09）。
- **现状**：功能齐备 + SIMD 加速，输出与外部 `spoa` 一致（§6 对照表）；
  命令用法见 `docs/fas.md`。
