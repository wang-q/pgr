# minigraph 分析笔记

本文档记录 minigraph 项目的架构、核心算法与数据结构，并分析其对 `pgr` 项目的启示。 minigraph
（Heng Li 开发，约 1.3 万行 C 代码）是一个**构建参考型泛基因组图**的工具，核心是用 rGFA
格式建模泛基因组图，通过增量映射-增强的方式把多个 assembly 合并进图。

与 `pgr` 的对照参考见 [[paf-route.md]]（路线决策）、[[graph-design.md]]（图构建层设计）、
[[impg.md]]（隐式图路线）、[[cactus.md]]（Caf 退火-熔化）。

---
## 1. 项目概览

### 1.1 设计哲学：参考锚定的增量图

minigraph 走的是"参考锚定 + 增量增强"路线，与 impg（隐式图）和 pggb（POA 物化图）都不同：

- **参考锚定**：第一个输入 assembly 作为参考骨架，后续 assembly 通过映射定位到现有图
- **增量增强**：每次映射后，把 query 中"映射不良"的区段作为新 segment 插入图
- **rGFA 坐标系**：每个 segment 携带 `SN:Z:`/`SO:i:`/`SR:i:` 三 tag，标记其在某个 stable sequence
  （参考路径）上的偏移，提供稳定坐标

这种设计让 minigraph **极快**（k-mer 免比对，不需要 lastz/wfmash）且**线性可扩展** （每个新
assembly 只需一次映射），但代价是**依赖参考选择**——参考选不好会导致图偏斜。

### 1.2 与 pgr 起点的差异

| 维度     | minigraph                                  | pgr                                |
|----------|--------------------------------------------|------------------------------------|
| 输入     | FASTA（需自己跑比对）                      | 已有 MAF/PAF（复用 pairwise 资产） |
| 比对方式 | minimizer 种子 + chaining                  | 已有 lastz→chain→net→axt→maf       |
| 图模型   | 显式 GFA（`gfa_t`）                        | 隐式图（PAF 边集 + 区间树）        |
| MSA 方式 | reference-guided 线性插入（`mg_path2seq`） | POA（`libs/poa/`）                 |
| 坐标系   | rGFA 三 tag（SN/SO/SR）                    | PAF target 坐标（0-based forward） |
| 外部依赖 | 零（纯 C，自包含）                         | 零（纯 Rust POA）                  |

pgr 已在 [[paf-route.md §1]] 论证：pgr 不需要 minigraph 的比对能力，因为 pgr 已有更成熟的 pairwise
基础设施。本笔记关注 minigraph 的**图算法层**对 pgr 的启示。

---
## 2. 模块分层与代码地图

```
入口层    main.c / options.c            命令分发、选项
IO 层     bseq.c / gfa-io.c / format.c  FASTA/FASTQ/GFA/GAF 读写
索引层    index.c / sketch.c            minimizer 索引
映射层    map-algo.c / lchain.c / gchain1.c / galign.c / miniwfa.c
          种子→线性链→图链→精细对齐
图构建层  ggen.c / ggsimple.c / gfa-aug.c   增量增强
图算法层  gfa-base.c / gfa-ed.c / gfa-bbl.c / shortk.c
          基础操作/GWFA/Bubble/K最短路径
后处理层  gcmisc.c / cal_cov.c / asm-call.c 排序/过滤/覆盖度/变异调用
```

各文件行数（粗略）：

- `ggsimple.c`（700 行）— 增量图构建核心
- `map-algo.c`（500 行）— 序列到图映射
- `gfa-ed.c`（600 行）— GWFA 图编辑距离
- `gchain1.c`（600 行）— 图 chaining DP
- `gfa-io.c`（400 行）— GFA 读写
- `gfa-base.c`（450 行）— GFA 基础操作
- `miniwfa.c`（700 行）— mini WFA 实现
- `lchain.c`（450 行）— 线性 chaining
- `gfa-aug.c`（300 行）— 图增强
- `gfa-bbl.c`（370 行）— Bubble calling
- `shortk.c`（250 行）— K 最短路径
- 其余文件均 < 300 行

---
## 3. 核心数据结构

### 3.1 `gfa_t`（GFA 图）

定义在 `gfa-priv.h`，操作在 `gfa-base.c`：

- `gfa_seg_t` 数组：节点（segment），字段 `seq`/`len`/`rank`/`snid`/`soff`
- `gfa_arc_t` 数组：有向边，字段 `v_lv`（长度）/`ov`/`ow`（overlap）
- `gfa_sseq_t` 数组：stable sequences（参考路径），字段 `name`/`rank`/`min`/`max`

**顶点编码**：`v = seg_id << 1 | strand`，每个 segment 有正反向两个顶点。 `v^1` 取反向顶点，`v>>1`
取 segment id。这个编码贯穿整个 minigraph，简化了 strand 处理。

**关键约束**：minigraph 不支持 overlap segments（`mg_gfa_overlap` 检查），所有 arc 的 `ov`/`ow`
必须为 0。这与 GFA spec 允许 overlap 的设计不同，是 minigraph 的简化。

### 3.2 rGFA 三 tag

rGFA 是 GFA 1.0 的扩展，给 segment 加三个 tag：

- `SN:Z:` — stable sequence name（参考路径名）
- `SO:i:` — stable offset（在参考路径上的偏移）
- `SR:i:` — rank（0=参考路径，> 0=非参考）

这三个 tag 提供**稳定坐标系**：即使图后续被增强（插入新 segment、分割旧 segment），
参考路径上的坐标仍然可追溯。这是 minigraph 区别于普通 GFA 工具的核心。

### 3.3 `mg_idx_t`（minimizer 索引）

定义在 `mgpriv.h`，操作在 `index.c`：

- 分桶哈希表（`mg_idx_bucket_t`），桶数 `1<<b`
- **特殊编码**：出现 1 次的 minimizer 直接存 value；出现多次的存位置数组 `p[]`
- `gfa_edseq_t`：每个顶点的正反向序列缓存，供 GWFA 使用

这种"出现 1 次特殊编码"的设计在 minimap2/minigraph 中一脉相承，省内存。

### 3.4 `mg_gchains_t`（图链集合）

图映射结果的顶层容器，三层结构：

- `gc[]`（gchain）：图链，跨多个 segment
- `lc[]`（llchain）：线性链，单 segment 内的种子链
- `a[]`（anchor）：单个种子（minimizer 命中）

每层通过 `off`/`cnt` 字段引用下一层。这种"三级索引"结构在 `gcmisc.c` 的 `mg_gchain_restore_order`/
`mg_gchain_restore_offset` 中维护。

---
## 4. 核心算法流程

### 4.1 增量图构建（`mg_ggen_aug`，ggen.c）

```
for each input assembly:
  1. mg_index(g)              对当前图建索引
  2. ggen_map()               将 assembly 映射到图
  3. mg_ggsimple()            识别映射不良区段
  4. gfa_augment()            将新序列作为 segment 插入图
```

**关键点**：每加入一个 assembly 都要**重建索引**（因为图变了）。这是 minigraph 线性但
非增量的代价——索引不能复用。pgr 的 PAF 索引是静态的（构建一次查询多次），无此问题。

**粗框架哲学（≥100bp SV 过滤）**：minigraph 在第 3 步 `mg_ggsimple` 只把长度差 ≥ `min_var_len`
的变异插入图（[ggsimple.c#L213](file:///Volumes/ExtHome/Scripts/pgr/minigraph-master/ggsimple.c#L213)，
源码默认 50，论文 L153/L384 称 100bp）。论文 L601-609 给了四条理由：

1. 图会爆炸——"composed of millions of short segments"
2. minigraph 会失败——"Not indexing minimizers across segments, minigraph will fail to seed"
3. 小变体用标准方法更易分析——"small variants are easier to analyze with the standard methods"
4. 无算法能为数百人类基因组构建这种复杂图

**关键区分**：这是**图构建层**的过滤，不是查询层。minigraph **保留**完整比对（base-level CIGAR），
只是不把小变体变成图节点。pgr 的隐式图天然避开这个问题——V1-V3 查询层全量返回同源区段
（[[paf-route.md]] §2.3），V4 GFA 物化时才需引入同等过滤（[[graph-design.md]] §4.3）。

### 4.2 序列到图映射（`mg_map_frag`，map-algo.c）

```
1. collect_minimizers    对 query 做 minimizer sketch（sketch.c）
2. collect_seed_hits     在图索引中查每个 minimizer 的位置
3. mg_lchain             线性 chaining（同 segment 内，lchain.c）
4. mg_gchain1_dp         图 chaining（跨 segment，经 arc 连接，gchain1.c）
5. miniwfa / GWFA        精细对齐（segment 内 / 跨边界，miniwfa.c / gfa-ed.c）
```

**线性 chaining**（`lchain.c`）：DP 和 RMQ 两种实现，把同 segment 内的种子连成链。 **图 chaining**
（`gchain1.c`）：把线性链当节点，用最短路径算链间图距离，DP 找最优组合。**精细对齐**：segment 内用
miniwfa（WFA），跨边界用 GWFA（图扩展 WFA）。

### 4.3 GWFA：图扩展 WFA（gfa-ed.c）

GWFA 是 WFA 的图上扩展，wavefront 推进时跨 segment 边界。核心数据结构：

- `gwf_diag_t`：对角线（vertex, diagonal, k, traceback）
- `gwf_intv_t`：对角线区间（处理 reach-end-of-vertex 的情况）
- `gwf_trace_t`：traceback 栈

**关键算法**：

- `gwf_ed_extend`：wavefront 推进，处理四种情况（中间/vertex 末/query 末/双末）
- `gwf_dedup`：wavefront 去重（interval merge + diagonal dedup）
- `gwf_prune`：剪枝（去除远落后于最远 wavefront 的对角线）

GWFA 的"forbidden bands"机制（`gwf_mixed_dedup`）用区间合并处理 vertex 边界， 这是图扩展 WFA
区别于线性 WFA 的核心。

### 4.4 Bubble calling（gfa-bbl.c）

用 Tarjan SCC 算法识别图中的 bubble 结构：

- `gfa_scc1`：单源 SCC，返回 `gfa_sub_t`（子图）
- `gfa_bubble`：遍历所有 stable sequence 的起点，找 bubble
- `bb_n_paths`：数 bubble 内的路径数（DP）

每个 bubble 记录：

- `vs`/`ve`：起止 vertex
- `ss`/`se`：起止 stable offset
- `len_min`/`len_max`：最短/最长路径长度
- `n_paths`：路径数
- `seq_min`/`seq_max`：最短/最长路径序列
- `is_bidir`：是否涉及双链（inversion）

### 4.5 K 最短路径（`mg_shortest_k`，shortk.c）

AVL 树 + Dijkstra 的 K 最短路径实现：

- 每个顶点维护大小为 `max_k` 的 max-heap，存到达该顶点的 K 条最短路径
- 用 `target_dist` + `target_hash` 支持目标导向搜索
- 返回 `mg_pathv_t[]` 回溯数组

用于图 chaining 中计算线性链之间的图距离。`MG_MAX_SHORT_K` 是上限。

### 4.6 图增强（gfa-aug.c）

`gfa_augment` 把插入（insertion）应用到图：

1. 分割现有 segment（如果插入点在中间）
2. 创建新 segment（插入序列）
3. 更新 arc（删除旧 arc，添加新 arc）

`gfa_ins_adj` 调整插入坐标，处理相邻插入的边界情况。

---
## 5. 与 pgr 路线的对照

[[paf-route.md §2]] 已明确 pgr 的核心决策，下面分析 minigraph 各部分对 pgr 的适用性。

### 5.1 pgr 不复用 minigraph 的 `gfa_t` 数据结构

[[graph-design.md §2]] 已分析：在 Rust 中重建 `gfa_t` 需要：

- 节点定义（CIGAR block vs reference 坐标切分）
- 边定义（strand → 四种边方向）
- 路径定义（P-line + SN/SO/SR tag 管理）
- 节点去重（`gfa_aux_update_cv`/`gfa_sort_ref_arc`）
- 自环/重复处理（`mg_gchain_set_parent`）

工作量远超 100 行，且 V1 只需坐标输出不需 GFA。**结论：V1 不物化 GFA，推迟到 V4**。

### 5.2 pgr 不复用 minigraph 的映射算法

minigraph 的映射（minimizer → linear chain → gchain → GWFA）是为"在已有图上定位 query" 设计的。pgr
的场景是**已有 pairwise 比对（MAF/PAF）**，不需要重新做比对。

[[paf-route.md §1.2]] 已论证：pgr 不需要 `--sparsify`，不需要 wfmash，不需要 minimizer chaining，
因为 MAF 里的每对已经跑过 pairwise 了。

### 5.3 minigraph 对 pgr 仍有价值的部分

#### (1) Bubble 模型作为查询层后处理过滤

minigraph 的 bubble（[gfa-bbl.c](file:///Volumes/ExtHome/Scripts/pgr/minigraph-master/gfa-bbl.c)）
用 Tarjan SCC 识别，[asm-call.c](file:///Volumes/ExtHome/Scripts/pgr/minigraph-master/asm-call.c)
基于 bubble 做变异调用。pgr 虽然不构建 GFA，但**PAF 传递闭包的连通分量**在概念上等价于 bubble。

[[paf-route.md §4.4]] 已记录 Caf 的 melting 过滤维度可作为传递闭包后处理。 minigraph 的 bubble 提供
**正交视角**：

- Caf 是离线全局过滤（图构建时）
- minigraph bubble 是在线局部结构（图查询时）
- pgr 的 BFS 传递闭包结果天然就是"隐式 bubble"

可借鉴 minigraph bubble 的指标作为 pgr 传递闭包的过滤维度：

- `n_paths`（路径数）→ pgr 的 `--min-degree N`
- `len_min`/`len_max`（长度区间）→ pgr 的 `--min-chain-length N`
- `is_bidir`（双链）→ pgr 的 inversion 标注

**注意**：这些是传递闭包的**后处理过滤**，不是 BFS 本身的中断条件（查询时无法做全图 SCC）。

#### (2) `mg_path2seq` 的 reference-guided 思路

`mg_path2seq`（[ggen.c](file:///Volumes/ExtHome/Scripts/pgr/minigraph-master/ggen.c)）本质是
"在参考序列上按位置依次插入 query 序列段"，不是 all-vs-all MSA。算法循环：

```
while (1) {
  1. 找 rs ≤ r ≤ re 的得分最高 chain
  2. 有 → 写 ref 片段 + query 序列，前进 v
  3. 无 → 写剩余 ref 片段，结束
}
```

[[graph-design.md §1.2]] 已指出这启示 pgr：当 cohort 有明确 reference 时，
**reference-guided 线性 MSA 比 POA 更快且无分支膨胀**。pgr 可在 `pgr paf query -o maf` 中根据是否有
`--reference` 参数选择后端：

- 有 reference → `fas multiz`（banded DP，reference-guided）
- 无 reference → `fas consensus`（SPOA，无参考）

#### (3) K 最短路径（`mg_shortest_k`）对图距离的启发

[shortk.c](file:///Volumes/ExtHome/Scripts/pgr/minigraph-master/shortk.c) 的 K 最短路径用于 图
chaining 中计算线性链之间的图距离。pgr 的 BFS 传递闭包目前只做"可达性"，不计算"图距离"。

如果未来需要**按同源紧密度排序**传递闭包结果，可借鉴 minigraph 的：

- `target_dist` + `target_hash` 机制（目标导向搜索）
- AVL 树 + max-heap 的 K 最短路径实现

但这是**远期需求**，V1 不需要。

#### (4) 覆盖度计算模型

[cal_cov.c](file:///Volumes/ExtHome/Scripts/pgr/minigraph-master/cal_cov.c) 的 `mg_coverage_asm` 计算：

- segment 覆盖度（区间合并）
- link 覆盖度（arc 计数）

pgr 的 PAF 区间树已支持区间查询，可类似地计算"每个 query 区间被多少 pairwise 比对覆盖"， 作为
**传递闭包置信度**的量化指标。这与 [[paf-route.md §4.4]] 的 Degree 过滤对应。

#### (4b) 粗框架过滤的两种正交维度

minigraph 的 `--min-var-len`（默认 100）是按**变异长度**过滤的粗框架（§4.1）。seqwish 提供另一种
正交维度：`--repeat-max` / `--min-repeat-dist` 按**重复拷贝数**过滤——限制同一序列在图同一位置的
拷贝数，避免高拷贝重复把图吹爆（详见 [[seqwish.md]] §3.5、§6.2）。

pgr V4a 物化粗 GFA 时可同时启用两种过滤：

- `--min-var-len 100`（minigraph 风格）——过滤 < 100bp 的小变体
- `--repeat-max N`（seqwish 风格）——限制重复拷贝数

两者维度不同，互不替代。详见 [[graph-design.md]] §4.3.5。

#### (5) GAF 紧凑路径编码

[format.c](file:///Volumes/ExtHome/Scripts/pgr/minigraph-master/format.c) 的 `mg_write_gaf` 实现
**紧凑路径编码**：

- 当连续 segment 属于同一 stable sequence 且连续时，合并为 `chr:start-end`
- 否则展开为 `>seg1<seg2...`

pgr 的 `pgr paf query` 输出同源区间列表时，可借鉴这种"能合并就合并，不能就展开"的双模式输出。

---
## 6. pgr 相对 minigraph 的独有优势

### 6.1 复用 pairwise 资产

minigraph 必须自己跑比对（minimizer chaining），pgr 复用已有 MAF/PAF。 见 [[paf-route.md §1.2]]。

### 6.2 Chain/Net syntenic 验证

minigraph 没有 UCSC Chain/Net 体系，pgr 可用 Chain/Net 做同源置信度标注。 见 [[paf-route.md §2.1]]。
这是"复用已有 pairwise 基础设施"的深层含义：不仅复用比对数据，还复用比对数据的**质量注释**。

### 6.3 查询层挑选

minigraph 的图是预先构建的（构建时就要决定参数），pgr 的隐式图支持查询时按 `--min-identity`
等参数动态过滤。见 [[paf-route.md §2.3]]。

### 6.4 MSA 质量可能更优

pgr 的 `fas_multiz.rs`（banded DP）对 core 区段比 minigraph 的 reference-guided 线性插入 更精确。
见 [[paf-route.md §2.4]]。

---
## 7. 结论与行动建议

### 7.1 结论

minigraph 的核心价值在于**证明了一条完整的"图构建→图映射→图增强"管道可行**， 但其 `gfa_t` 数据结构和
minimizer chaining 算法对 pgr **不直接适用**——pgr 已有更成熟的 pairwise 基础设施和 PAF 隐式图。

minigraph 值得借鉴的是**算法思想**而非具体实现：

- Bubble 作为传递闭包后处理过滤的结构化指标
- reference-guided vs POA 的 MSA 后端选择策略
- K 最短路径用于按紧密度排序同源结果（远期）
- GAF 的紧凑路径输出双模式
- 覆盖度作为置信度量化指标

**物化图的两条路径**：minigraph（增量增强，输入 FASTA 自跑 minimizer chaining）与 seqwish（PAF 诱导，
输入已有 pairwise 比对）是物化 GFA 的两条不同路径。pgr V4a 输入是 PAF，与 seqwish 同源，因此
算法骨架（spanning tree → BFS → DSU → compact → links → GFA）直接复用 seqwish（详见
[[seqwish.md]] §6.2）；minigraph 的 `--min-var-len` 粗框架过滤哲学则作为正交补充（§4.1、§5.3(4b)）。

### 7.2 行动建议

对 pgr V3（`pgr paf query -o maf`）的影响：

- **不影响**：V3 继续走 PAF → POA → MSA 路线（[[graph-design.md §3]]），约 150 行新代码
- **可借鉴**：`pgr paf query` 输出格式可参考 GAF 紧凑路径编码
- **可借鉴**：`pgr paf query --transitive` 的后处理过滤可参考 minigraph bubble 指标

对 pgr V2/V3 的影响：

- V2 评估是否引入 bubble 指标作为传递闭包过滤维度
- V3 评估是否引入 rGFA 标准（届时再考虑 Rust 版 `gfa_t`）

---
## 8. 附录：与其他 notes 文档的引用关系

```
minigraph.md (本文档) ─ 架构参考 ──┐
    │ §5.1 gfa_t 不复用 → graph-design.md §2  │
    │ §5.3(1) bubble → paf-route.md §4.4       │
    │ §5.3(2) reference-guided → graph-design.md §1.2
    │ §6 pgr 优势 → paf-route.md §1-2          │
    │                                           │
cactus.md ────────────── 架构参考 ─────────────┤
    │ §8 Caf 退火-熔化 → paf-route.md §4.4     │
    │ §3 Minigraph-Cactus → paf-route.md §4.5  │
    │                                           │
impg.md ──────────────── 路线参考 ─────────────┤
    │ §4 传递闭包 → paf-route.md §2.4          │
    │ §9 启示 → paf-route.md                   │
    │                                           │
paf-route.md (路线决策) ──────────────────────┤
    │ §1 起点差异                               │
    │ §2 核心决策                               │
    │ §4 存量资产优势                           │
    │                                           │
graph-design.md (图构建层设计) ───────────────┘
    │ §1 三种图构建路线
    │ §2 V1 不做 GFA
    │ §3 V1 最小可行实现
```

