# PAF 隐式图路线决策

本文档记录 pgr 走向泛基因组时的**路线选择与理由**。它回答"为什么走这条路"，不涉及具体实现细节。
实现参考见 [[paf-implementation.md]]，第一步行动计划见 [[pairwise-selection.md]]。

参考文档：[[impg.md]]（隐式图路线与传递闭包）、[[cactus.md]]（Caf 退火-熔化与 Minigraph-Cactus
分治）、[[cactus_lastz.md]]（pgr 已有 lastz 比对链的能力证明）、[[project-understanding.md]]（pgr
现状基线）。

---
## 0. 项目核心目标

pgr 走向泛基因组的核心目标是：
**复用已有的 pairwise 比对基础设施，构建 PAF 隐式图， 按需回答"哪些序列的哪些区段同源"，而非物化一张泛基因组图。**

### 0.1 四条原则

1. **隐式图优先，粗 GFA 作为可选投影**——默认用 PAF 索引 + 区间树 + BFS 传递闭包做"按需图遍历"，
   不物化 GFA。隐式图的优势是"不构建整张图"（[[impg.md]] §1.1.2），同一份索引可服务不同严格度的查询。
   粗全局 GFA 是索引的**显式投影**（V4a），数据源仍是 PAF 索引，不引入新的真实源——这与 minigraph
   必须物化图（它没有隐式图）是本质区别。**粒度差异**：seqwish 的传递闭包是全局一次性（一次性算出
   全部等价类再写图），pgr 是局部按需（每次查询从一个区间出发 BFS，只算相关等价类）——
   详见 [[seqwish.md]] §5。
2. **复用 pairwise 资产，大 cohort 用 Mash KNN sparsify**——pgr 已有成熟的 pairwise 比对链
   （`pgr lav lastz` → chain → net → axt → maf），这是区别于 impg（从 FASTA 跑 wfmash）和
   pggb（从 FASTA 跑 wfmash）的独特起点。**小 cohort + 已有 MAF**时直接复用，不需要 sparsify
   （§1.2）；**大 cohort + 无先验**（如 4 万 E. coli）时用 Mash KNN sparsify 把 N² 降到 N×K，
   传递闭包推断未比对的对（§1.2）。
3. **查询层全量，图构建层粗框架**——两个层次分离（§2.3）：
    - **查询层**（V1-V3）：全量返回同源区段，由用户用 `--merge-distance` 等参数控制粗细
    - **图构建层**（V4a+）：物化 GFA 时才引入 `--min-var-len`（默认 100）粗框架过滤，对齐 minigraph
      （[[minigraph.md]] §4.1）
4. **pipe 友好，两段式 GFA**——遵循 Unix 哲学。V1 默认 PAF 输出、`-o bed` 可选（坐标投影的轻量产物，
   最 pipe 友好），`bed → fa range → fas consensus` 的 MSA 路径已通，不需要在 query 层重复实现 POA
   （§5.2）。GFA 采用两段式（§0.3 V4）：粗全局 GFA 提供"地图"（哪里有大 SV），区域精细 GFA 按需做
   碱基级分析——这比 minigraph（只有粗）和 impg（只有精）更完整。

### 0.2 能力栈与当前进度

```
索引层 ✅  |  查询层 ✅  |  图构建层 V4a ✅（[[graph-design.md]] §4.3.1）|  应用层 ← 远期
```

impg 的四层能力栈（[[impg.md]] §1.1.3）在 pgr 的现状：pairwise 比对层资产比 impg 更成熟，
索引层和查询层已补齐，下一步聚焦图构建层。

### 0.3 V1-V5 路线一览

| 阶段    | 内容                                                                   |  状态  |
|---------|------------------------------------------------------------------------|:------:|
| **V1**  | `pgr paf query -o paf`（默认）+ `-o bed` + `-b regions.bed` 批查 + `--min-degree`/`--min-chain-length` 后处理过滤 | ✅ 已完成 |
| **V2**  | `pgr paf to-maf`（pairwise MAF，按 CIGAR 直接还原，需 `-f TSV`）       | ✅ 已完成 |
| **V3**  | `pgr paf to-maf --msa`（POA MSA，多序列合并，需 `--transitive` + POA） | ✅ 已完成 |
| **V4a** | 粗全局 GFA（`pgr paf graph -f refs.fa --min-var-len 100`，seqwish DSU 风格） | ✅ 已完成 |
| **V4b** | 区域精细 GFA（`pgr paf query -o gfa -r region`，impg 风格）            | 待评估 |
| **V5**  | 区域 GFA → MAF/VCF（精细分析输出）                                     | 待评估 |

注：早期规划中的 `-o fasta`（未比对序列提取）经评估不做——`pgr fa range` 已提供独立提取能力，
用户场景不需要在 query 路径耦合裸序列输出（见 [[graph-design.md]] §4.2）。

详细路线与代码量估计见 [[graph-design.md]] §4。V4 采用两段式（粗全局 + 区域精细）， V4a 与 V4b
可独立实现，互不依赖。

---
## 1. pgr 与 impg 的起点差异

pgr 走向泛基因组时，面对的问题与 impg **完全不同**。impg 的起点是"只有 FASTA，没有 pairwise 比对"，
需要先选对、再比对、再索引。pgr 的起点是"已有成熟的 pairwise 比对基础设施"，需要的是 "复用已有资产，
补上缺失的图遍历层"。

### 1.1 差异对照

| 维度     | impg                       | pgr                              |
|----------|----------------------------|----------------------------------|
| 比对来源 | 从 FASTA 跑 wfmash/sweepga | 已有两序列 MAF（可转 PAF）       |
| 挑选时机 | align 阶段（无先验）       | 可借已有 MAF 先验                |
| 网络结构 | 直接 all-vs-all            | 已有 pairwise 基础设施           |
| 核心问题 | 选哪些对比对               | 复用已有 pairwise，做 PAF 隐式图 |
| 比对工具 | wfmash/FastGA              | pgr 已有 `pgr lav lastz` 全套    |

### 1.2 pgr 与 `--sparsify` 的关系：分场景

impg 的 `--sparsify auto`（[[impg.md]] §6.4）用 Mash KNN 从 N 个基因组中选 K 个近邻做比对， 把 N²
降到 N×K。**这是隐式图架构避免 N² 爆炸的核心机制**——稀疏比对 + 查询时 BFS 传递闭包 推断未比对的对
（A↔B、B↔C 已比对，查询 C 时经 B 到达 A，即使 A↔C 从未直接比对）。

pgr 是否需要 sparsify **取决于 cohort 规模和已有资产**，分两种场景：

| 场景                     | 规模          | 已有资产          | 是否需要 sparsify                       |
|--------------------------|---------------|-------------------|-----------------------------------------|
| **小 cohort + 已有 MAF** | 几十基因组    | 现成 pairwise MAF | **不需要** — MAF 里的对已跑过比对       |
| **大 cohort + 无先验**   | 27000 E. coli | 只有 FASTA        | **必需** — 否则 27000² ≈ 3.6 亿对不可行 |

**小 cohort 场景**（本文档 §1 原本针对的）：

- **不需要选对** — MAF 里的每对已经跑过 pairwise 了
- **不需要 wfmash** — 即使要补充新比对，pgr 有完整的 `pgr lav lastz` 链（见 §3.2）
- **挑选发生在查询层** — pgr 的"挑选"是"查询时用 `--min-identity` 等参数过滤 PAF 记录"，
  不是"选哪些对跑比对"

**大 cohort 场景**（如 [[ecoli-cohort.md]] 的 4 万 E. coli）：

- **必需 sparsify** — Mash KNN 把 27000² 降到 27000×K（K≈50，约 135 万对，~270 倍缩减）
- **用 FastGA 比对** — 对 K×N 条边跑 FastGA（external，细菌基因组优化，比 lastz 快）， 产 MAF → 转
  PAF。FastGA 贯穿去冗余精比（步骤 1b）和 sparsify 比对（步骤 2c），保持工具一致。`pgr lav lastz`
  仍作为 pgr 已有的 fallback，但大 cohort 首选 FastGA。
- **传递闭包补全** — 稀疏比对的缺口由查询时 BFS 推断，不需要补全 A↔C 直接比对

> **与 PGGB 的对照**：PGGB 走 all-vs-all 无偏路线（`wfmash` + `seqwish` + `smoothxg`），O(N²)
> 限制规模到 ~百单倍型。pgr 用 Mash KNN sparsify + 传递闭包是同一无偏思想的可扩展替代路径。
> 详见 [[pangenome-tools.md]] §3.2。

> **spanning tree 优化（远期）**：seqwish 在传递闭包前用最大权生成树剪枝，把 N(N-1)/2 边压缩
> 到 N-1 边（[[seqwish.md]] §3.1）。pgr 查询层 BFS 在 135k 边的稠密对齐网络上若性能瓶颈显现，
> 可在加载 PAF 阶段预计算生成树，查询时优先沿树边走，显著减少 BFS 边遍历数。当前 V1 不做，
> 待性能数据出来再评估。

### 1.3 能力栈映射

impg 的四层能力栈（[[impg.md]] §1.1.3）在 pgr 的现状：

```
索引层 — ✅ pgr paf index + build_multi（多文件合并）
查询层 — ✅ pgr paf query + --transitive BFS
图构建层 V4a — 已实现（[[graph-design.md]] §4.3.1，seqwish DSU 风格）
应用层 — 远期
```

pgr 在 pairwise 比对层的资产比 impg 更成熟，索引层和查询层现已补齐。 下一步聚焦图构建层。
当前实现记录见 [[pairwise-selection.md]]。

---
## 2. 核心决策

以下决策是后续所有行动的**不变前提**。

### 2.1 用 PAF 作隐式图边集，不用 Chain

**决策**：PAF 是图的边，Chain/Net 是查询层的 syntenic 过滤器。

**理由**（详见 [[impg.md]] §9.1）：

- Chain 是 star topology（ref↔query_i），做传递闭包时 ref 成为必经枢纽，ref 缺失区段会断开间接同源路径
- Chain 已被 UCSC 流程过滤（score 阈值、syntenic 净化），不是原始比对，丢失了 paralog/低质量区间
- Chain 的 gap-less tBlock 分段结构不适合做图边——转换为 PAF 会丢信息且无收益

**Chain/Net 的正确角色**：

1. PAF 边集提供"所有可能同源"（全量装入，不过滤）
2. Chain/Net 提供 syntenic 验证——如果两条 PAF 声称 A↔B 同源但该区段不在 Chain/Net 中，标记为低置信度
3. 查询时用户可选择过滤级别：
    - `--syntenic-filter strict` — 只保留 Chain/Net 验证过的
    - `--syntenic-filter lenient` — 保留全部但标注置信度
    - `--syntenic-filter none` — 不看 Chain/Net

这个三角关系是 pgr 独有的优势——impg 没有 UCSC Chain/Net 体系，无法做这种 syntenic 验证。
这是"复用已有 pairwise 基础设施"的深层含义：不仅复用比对数据，还复用比对数据的**质量注释**。

### 2.2 PAF 来源：MAF → PAF 转换，不跑新比对

**决策**：pgr 已有的两序列 MAF 直接转换为 PAF。不引入 wfmash。

**理由**：两序列 MAF 的每个 block 等价于一条 pairwise alignment——`s` 行给出坐标和链向， 可直接映射到
PAF 的 12 列。这是 pgr 复用已有 pairwise 基础设施通往 PAF-based 隐式图的 天然桥梁（[[impg.md]]
§9.2）。

### 2.3 索引全量装入，挑选发生在查询层

**决策**：PAF 索引时不做过滤，所有记录全量装入区间树。过滤参数只在查询时生效。

**理由**：impg 的核心哲学——"比对即图"，索引保留所有边，挑选推迟到查询。同一份索引
可服务不同严格度的查询（[[impg.md]] §1.1.3）。对应 impg 的 `Index` 命令只有文件路径和 index-mode
参数，而 `QueryOpts` 才有 `-d`/`--min-result-identity`/`-l` 等过滤开关。

**与 minigraph 粗框架的关系**：minigraph 的 ≥100bp SV 过滤是**图构建层**的事（[[minigraph.md]]
§4.1），不是查询层。pgr V1-V3 查询层全量返回同源区段，由用户用 `--merge-distance` 等参数控制
粗细；V4 GFA 物化时才在 graph engine 内部做 `min_var_len` 过滤（[[graph-design.md]] §4.3）。
两个层次不能混淆——查询层全量是"让用户决定粗细"，图构建层粗框架是"避免图爆炸"。

**索引层工程优化（4 万大肠杆菌规模可考虑）**：seqwish 的 `SeqIndex` 用 FM-index 索引序列名 +
`SparseBitVec`（只存 1-bit 位置数组，select O(1)）记录序列边界，比 HashMap 省内存；`PosT` 把
offset+方向打包进单 u64，单棵区间树同时表达正反链对齐。详见 [[seqwish.md]] §2.1、§2.2、§6.1。
pgr 当前用 HashMap + coitrees 在 4 万大肠杆菌规模尚可，若扩到 HPRC 规模（数百单倍型、Gb 级），
可借鉴此方案。

### 2.4 传递闭包是图遍历，不是多序列比对

**决策**：传递闭包做"图遍历可达性"，不产出多重比对。找到所有同源片段后，如需 MSA， 再调用
`fas consensus`（SPOA）或 `fas multiz`（banded DP）。

**理由**（[[impg.md]] §4.3）：图遍历和 MSA 是正交的两个步骤——

- 传递闭包告诉你"哪些序列的哪些区段同源"
- MSA 告诉你"这些同源区段具体如何对齐"
- pgr 的 MSA 基础设施（`libs/poa/` + `libs/fas_multiz.rs`）已经就绪，不需要把两者耦合

pgr 的 MSA 质量可能优于 impg 的 per-bubble POA——`fas_multiz.rs` 实现了 banded DP 合并
（`FasMultizMode::Core`），对 core 区段比纯 POA 更精确。

### 2.5 第一步不物化 GFA

**决策**：第一步扩展 `pgr paf query -o fasta`/`-o maf`，不物化 GFA。GFA 物化推迟到 V3。
详见 [[graph-design.md]] §3-§4。

**理由**："先物化再分析"对 pgr 是过载的。隐式图的优势正在于"按需计算图遍历，不构建整张图"（[[impg.md]
] §1.1.2）。minigraph 的 chain → GFA 路线需要 `gfa_t` 数据结构（`gfa-priv.h`）的重建，impg 的 POA
→GFA 路线需要 spoa_rs + gfaffix + gfasort 外部依赖链。pgr 的 `libs/poa/` 纯 Rust 引擎可直接输出
MSA，零外部依赖。

---
## 3. 三层挑选问题

pgr 的"挑选"不是 impg 的"选哪些对跑比对"，而是分三层（按实现优先级排序）：

### 3.1 第一层：从已有 MAF/PAF 挑选（查询层，无需新比对）

**问题**：cohort 有 N 个基因组，pgr 已有两序列 MAF（ref↔query_i 或 query_i↔query_j）。
要回答"query_A 的某区段在 cohort 中有哪些同源"，需要在 PAF 网络上跨记录传递。

**机制**：impg 的传递闭包（`-x` BFS，[[impg.md]] §4.2）——若 A↔B、B↔C 在同一区段有比对， 则 A↔C
间接同源。所有 pairwise 比对当作图的边集，从目标区间出发做 BFS，自动发现所有直接和间接同源片段。

**这一层不需要新比对**，只需把已有 MAF 转成 PAF 装入区间树，做查询层挑选。

### 3.2 第二层：补充 pairwise 比对的挑选（align 层）

**问题**：已有 MAF 只覆盖已跑过的对。以下场景 MAF 缺失或不足：

- cohort 加入新基因组，与已有基因组未跑过 pairwise
- 已有 MAF 在某区段断开（low identity 区段被过滤），但可能有间接同源
- 某些 sample 对需要更精细的 region-level 重比对

**是否需要补充**取决于应用：

- 单 locus 查询（HLA/C4）：第一层传递闭包通常够用
- 全 cohort 泛基因组图构建：需要补充 A↔B 直接比对

**补充比对的五种策略**：

| 策略                | 来源                   | 适用                   | pgr 实现门槛                     |
|---------------------|------------------------|------------------------|----------------------------------|
| 已有 PAF 覆盖度先验 | pgr 独有               | 已有部分 PAF 的 cohort | **推荐**，复用已有 PAF           |
| `pgr lav lastz`     | pgr + Cactus 风格      | 特定 pair 需要新比对   | `pgr lav lastz`（不含 `--self`） |
| 系统发育树引导      | Cactus 风格            | 有 phylogeny           | 复用 `pgr nwk` 模块              |
| Mash KNN            | impg `--sparsify auto` | 无先验全选             | 需引入 mash crate                |

**已有 PAF 覆盖度先验策略**（pgr 推荐）：

1. 对每个 query_i，统计其在已有 PAF 上的覆盖区间集合 C_i
2. 对 query_i、query_j，计算 |C_i ∩ C_j| / |C_i ∪ C_j|（Jaccard）
3. 选 Jaccard 高于阈值且**尚未跑过 pairwise**的对补充比对
4. 补充比对用 `pgr lav lastz`，结果转 PAF 并入区间树

这样把 N² 降到"PAF 覆盖度共享的子集"。

**lastz 策略**：

pgr 可以通过 `pgr lav lastz`（不含 `--self`）为特定 pair 生成 pairwise 比对。
`--self` 模式（`src/cmd_pgr/lav/lastz.rs:328-358`）是 Cactus 风格的**重复屏蔽**
管道的一部分——它要求 target 和 query 是同一文件（碎片化的单序列 chunks），
目的是检测基因组内重复区段，而非生成泛基因组的 pairwise 比对网络。`--self` 的正确用途是
`pgr fa window → pgr lav lastz --self → pgr lav to-psl → pgr psl lift → spanr coverage → pgr fa mask`，
详见 `notes/cactus_lastz.md` §5.6。

当 cohort 完全没有任何 pairwise 比对时，需要逐个 pair 跑 `pgr lav lastz` （不带 `--self`），然后
LAV → PSL → PAF 进入泛基因组管道。

### 3.3 第三层：region 级精细比对挑选

**问题**：已有 MAF 是粗粒度的（受原始比对参数限制）。某些 region（HLA、KIR、C4）需要更精细
pairwise，但全基因组精细比对代价高。

**挑选机制**：

- 从已有 PAF 的 gap/low-identity 区段筛选候选 region
- 对候选 region 跑 `pgr lav lastz`（Cactus 风格，已有）
- 合并回 PAF 网络

**这一层是第一层的补充**，不是泛基因组的核心路径，按需开启。

---
## 4. pgr 的存量资产优势（比最初估计更强）

通读 notes/ 下全部文档并分析 pgr 源码后，对 pgr 已有资产的认识持续深化。
以下三项发现显著降低了第一步的实现门槛。

### 4.1 `loc.rs` — pgr 的 IO 层比 impg 更成熟

分析了 `src/libs/loc.rs`（202 行）与 impg `paf.rs`（417 行）的对应关系（详见 [[paf-implementation.md]
] §10）。核心发现：

- **`Input` enum 比 impg 的 `PafHandle` 更强**：多了 `Buf` 变体（支持 stdin）， 且 `Bgzf` 变体使用
  `IndexedReader`（自带索引，seek 无需外部 `.gzi` 文件）
- **`read_offset()` 可直接替代 impg 的 `read_cigar_data()`**：同样是 seek+read+返回字节， pgr
  的实现更简洁（11 行 match + 2 行 I/O vs impg 的 46 行分支）
- **pgr 已有 BGZF 行读取能力**（`create_loc` 中对 `Input::Bgzf` 调用 `read_line`）， 只是需要抽象为
  `Input::read_line` 方法供 PAF 解析使用

**结论**：PAF 模块中最棘手的 IO 部分（多格式输入、BGZF 随机访问、CIGAR 懒加载） pgr 已经解决了 80%。
真正需要从零写的只有三样：区间树索引、PAF 行解析、CIGAR 编解码。这比 [[paf-implementation.md]]
最初估计的实现量减少了约 30%。

### 4.2 `IndexedReader` 自带索引能力，不需要 impg 的 GZI 机制

impg 的 `parse_paf_bgzf_with_gzi` 需要外部 `.gzi` 索引文件来做多线程解压，且需要
显式 `bgzf::VirtualPosition::from(offset)` 转换。pgr 的 `bgzf::io::IndexedReader`
在内部处理了这一切——调用者只需传字节偏移量。

这意味着 pgr 的 BGZF PAF 支持可以**跳过 impg 的模式 3**（GZI 两遍扫描）， 直接用 `IndexedReader`
做到同等性能。唯一需要注意的：CIGAR 懒加载需要记录 vpos，当前 `Input::Bgzf` 未暴露
`virtual_position()`，需要加一个方法。

### 4.3 pgr 已有的比对生成能力

写 §3.2 第二层时，最初假设 pgr **只能复用已有 MAF**。读完 `cactus_lastz.md` 和
`src/cmd_pgr/lav/lastz.rs` 后认识到 pgr 有完整的 lastz 封装（7 套预设参数、并行执行），可以为特定
pair 生成 pairwise 比对。`--self` 模式是重复屏蔽管道的 一部分（碎片自比对），不是泛基因组比对工具。

### 4.4 Cactus Caf 的"退火-熔化"循环对 pgr 挑选机制的直接参考

`cactus.md` §8 详细分析了 Caf 模块（`caf.c`、`annealing.c`、`melting.c`）的迭代循环：

- **Annealing（加法）**：把两两比对捏合成 Pinch Graph 中的 Block。关键约束是
  `stCaf_annealBetweenAdjacencyComponents`——"只连接不同连通分量的序列对"，避免在
  同一连通区域形成复杂环
- **Melting（减法）**：按 Degree、Tree Coverage、Chain Length 进行多维过滤，
  `stCaf_getBlocksInChainsLessThanGivenLength` 丢弃破碎短链

对应的过滤维度可以搬到 pgr 的查询层：

| Caf 过滤维度   | pgr 对应参数           | 语义                       | 状态          |
|----------------|------------------------|----------------------------|---------------|
| Degree         | `--min-degree N`       | 过滤支持序列数 < N 的区间  | ✅ V1 已实现  |
| Tree Coverage  | `--min-tree-coverage`  | 过滤进化树上分布稀疏的区间 | 待实现        |
| Chain Length   | `--min-chain-length N` | 过滤总长 < N bp 的传递链   | ✅ V1 已实现  |
| Block End Trim | `--end-trim N`         | 切除比对边缘不可靠的 N bp  | ⏭️ 推迟       |

**但要警惕**：Caf 的 melting 在**图构建时**做（离线、全局视角），而 pgr 的挑选在**查询时**做
（在线、局部视角）。查询时无法做全图 Tree Coverage 计算。因此这些 Caf 过滤维度更适合作为传递闭包的
**后处理过滤**，而非 BFS 本身的中断条件。V1 已实现 `--min-degree` 与 `--min-chain-length`（见
[[paf-implementation.md]] §8）；`--end-trim` 推迟——它需要 per-interval 修剪 CIGAR 两端，
与当前"区间整体投影"的输出模型不兼容，待 V2 引入序列输出时一并处理。`--min-tree-coverage` 需要
进化树上下文，留待后续阶段。

### 4.5 Minigraph-Cactus 分治策略对 pgr partition 的启示

`cactus.md` §3.1 详述了 Minigraph-Cactus 的五阶段流程：

```
Minigraph 骨架构建 → 图映射定位 → rgfa-split 切分 → 批量 Cactus 比对 → 合并
```

与 impg 的 Partition + Lace 模式（[[impg.md]] §6.3）对比：

| 维度     | Minigraph-Cactus               | impg Partition + Lace         |
|----------|--------------------------------|-------------------------------|
| 拆分依据 | Minigraph SV 图连通分量        | 传递闭包 + masking 去重       |
| 拆分粒度 | 染色体级（MB）                 | locus 级（KB-MB，窗口可控）   |
| 局部处理 | Cactus full pipeline (Caf+Bar) | 独立 GFA 构建 (seqwish/crush) |
| 合并方式 | HAL/VG join                    | lace（坐标驱动重新拼装）      |

**对 pgr 的启示**：如果 pgr 未来需要 partition（处理 > 100 基因组的 cohort）， 建议走
**Minigraph-Cactus 的"先粗后细"路线**：

- **粗拆分**：用已有的 Chain/Net syntenic 信息做染色体级拆分（类似 `cactus_graphmap_split.py` 的
  heuristic contig selection：regex + size + dropoff，见 `cactus.md` §2.4.2）
- **细拆分**：在每个大区块内，用传递闭包 BFS + masking 去重切分成 per-locus 批次
- **比对**：per-locus 跑 `pgr lav lastz`（不含 `--self`）生成局部 pairwise
- **合并**：per-locus PAF 汇总回全 cohort 区间树

这个路线复用了 pgr 的三项已有资产：Chain/Net syntenic 信息、`pgr lav lastz`、PAF 区间树。 但这是
**第二步或第三步的任务**，第一步不需要 partition。

### 4.6 文档间的引用关系

```
cactus.md ────────────── 架构参考 ──────┐
    │ §1.11 Chain/Net↔Flower 对应       │
    │ §8 Caf 退火-熔化 → §4.4 过滤维度  ├── paf-route.md
    │ §3 Minigraph-Cactus → §4.5 partition│   (本文档：路线决策)
    │                                    │
cactus_lastz.md ─────── 能力证明 ───────┤
    │ §5.6 完整 lastz 流程 → §3.2 策略  │
    │                                    │
impg.md ──────────────── 路线参考 ──────┘
    │ §4 传递闭包 → query --transitive
    │ §6.3 Partition → 未来第二步
    │ §9 启示 → 整体方向
    │
project-understanding.md ─ 现状基线
    │ §4 pgr 核心库 → 已有资产清单
    │ §6.3 待补全 → 第一步填补的目标
```

---
## 5. 暂不实现

以下功能明确排除在第一步之外。每条给了触发条件以防止 scope creep。

| 暂不实现                        | 理由                                       | 重新评估的触发条件                                                    |
|---------------------------------|--------------------------------------------|-----------------------------------------------------------------------|
| 补充 pairwise 比对（第二层）    | 第一步复用已有 MAF 已足够                  | 传递闭包覆盖率不足                                                    |
| Mash KNN pair-selection         | 小 cohort 有 MAF 先验时不需要              | 大 cohort 无先验（如 4 万 E. coli，见 [[ecoli-cohort.md]]）           |
| `pgr lav lastz --self` 全自比对 | 此 flag 用于重复屏蔽管道，非泛基因组比对   | 需要全新 cohort 的 pairwise 比对时评估 `pgr lav lastz`（不含 --self） |
| syng 免比对后端                 | [[impg.md]] §1.1 已明确不参考              | 永不                                                                  |
| GFA 物化（seqwish/crush）       | 隐式图 + 按需 MSA 已覆盖第一步场景         | 用户需要全图统计或 variant calling                                    |
| crush bubble 压缩               | pgr 还没有 GFA 图（[[impg.md]] §9.3）      | GFA 构建管道就绪后                                                    |
| partition / lace / refine       | 处理 >100 基因组的 cohort 时需要           | N > 50                                                                |
| stage DSL                       | 单命令不需要管道化                         | 出现三个以上 stage 串联                                               |
| 基因分型（genotype/infer）      | 能力栈顶端，依赖图构建层（[[impg.md]] §7） | 图构建层就绪后                                                        |

### 5.2 paf query 的输出格式策略

**PAF 是默认输出，BED 为 `-o bed` 可选**：impg `query` 默认 `-o bed`
（[main.rs#L4873](file:///Volumes/ExtHome/Scripts/pgr/impg-0.4.1/src/main.rs#L4873)），
README L20-23 明确："It outputs BED / BEDPE / PAF — ready to feed FASTA extraction, multiple
sequence alignment... can also emit GFA directly"。"can also" 表明 GFA/MAF 是**可选附加**，
不是核心。pgr V1 选择 **PAF 为默认**（含 CIGAR/gi/bi 完整比对记录），BED 通过 `-o bed` 可选——
这与 impg 的"BED 默认"不同，理由是 pgr 既有测试已断言 PAF 输出，且 PAF 对需要完整比对记录的场景更
直接。BED 三列（`name start end`）是坐标投影的轻量产物，最 pipe 友好，用 `-o bed` 显式切换。
详见 [[graph-design.md]] §3.1 默认输出决策。

**FAS（block FASTA with shared coords）不输出**：FAS 格式的核心假设是所有序列共享一个统一 坐标系
（通常以 reference 为锚），这在泛基因组场景不成立——PAF query 结果是各基因组**独立坐标系**
下的同源区段列表。

**fasta/maf 是可选附加，按依赖链后移**：对照 impg 源码后修订（见 [[graph-design.md]] §3.2）。 impg
的 11 种输出按"是否需要序列文件"分两类——坐标类（`bed`/`bedpe`/`paf`，不需 `-f`）是核心，序列类/MSA
类/图类（`fasta`/`maf`/`gfa` 等，需 `-f`）是可选。pgr 按此分阶段：

```bash
# V1：坐标输出（PAF 默认，BED 可选，pipe 友好，不需 -f）
pgr paf query ... --transitive -o paf   # 默认，完整比对记录
pgr paf query ... --transitive -o bed   # 轻量坐标，喂给 fa range

# V2：未比对序列（可选，需 -f）
pgr paf query ... --transitive -o fasta -f genomes.fa

# V3：POA MSA（可选，需 -f）
pgr paf query ... --transitive -o maf -f genomes.fa
```

**为何不在 query 层重复实现 MSA 后端**：`-o maf` 内部调用 `libs/poa`（pgr 已有的纯 Rust
POA），不重复 `pgr fas consensus` 的成熟设施。但 V1/V2 阶段，用户要 MSA 时走 pipe 路径已通：
`pgr paf query -o bed` → `pgr fa range` → `pgr fas consensus`。`pgr fas consensus` 支持 builtin
POA + 外部 spoa、可配分矩阵、并行处理、outgroup 支持——作为 MSA 后端已足够成熟。

**与 §2.4 的一致性**：传递闭包是图遍历（找全同源片段），`-o maf` 是图遍历后的 MSA 物化。 两者在
query 命令内串联（V3+），但逻辑上是正交步骤——BFS 结果可输出为 BED/PAF（无 MSA），也可输出为 MAF
（有 MSA），由用户选择。V1 先做核心的坐标输出。

### 5.3 impg 的做法：查询输出即最终输出

impg 有 8 种输出格式（BED/BEDPE/PAF/GFA/VCF/FASTA/MAF/FASTA_ALN）， 但它们都是
**直接从 `AdjustedInterval` 格式化输出**（`main.rs:11849-12444`），没有任何中间桥接层。
查询返回什么就输出什么——PAF 是最常用的格式，GFA 用于图构建，其他按需。

pgr 当前的 `-o paf`（默认）/`-o bed` 输出遵循相同模式。不需引入 FAS 作为中间格式。

### 5.4 pgr 已有的 MSA 资产（供后续阶段按需使用）

以下组件不在查询层使用，但在方向 D（图构建）或下游分析中可以直接调用：

| 组件               | 源码                          | 后续用途                                           |
|--------------------|-------------------------------|----------------------------------------------------|
| POA 引擎           | `libs/poa/poa.rs`             | 图构建阶段的 per-bubble 共识/比对                  |
| Banded DP          | `libs/fas_multiz.rs`          | partition 内多 pairwise 合并（比 impg POA 更精确） |
| `get_subs`         | `libs/alignment.rs:214`       | MSA 上的变体检测                                   |
| 裁剪函数           | `libs/alignment.rs:1351-1687` | BFS 结果边界清理                                   |
| crossbeam 并行管道 | `consensus.rs:250`            | `build_multi` 并行化                               |

但这些都是**独立的 CLI 命令或库函数**，通过 Unix pipe 组合，不与 `paf query` 耦合。

---
## 6. 附录：与其他文档的对照

### impg.md

| impg.md 章节                              | 本文档对应                           |
|-------------------------------------------|--------------------------------------|
| §1.1.3 能力栈四层（索引→查询→图→应用）    | §1.3 能力栈映射                      |
| §1.1.4 名词解释 pair-selection            | §3 三层挑选问题                      |
| §4.2 传递闭包 BFS（`-x`/`-m`/`-d`）       | §2.4（图遍历≠MSA）                   |
| §4.3 传递闭包 ≠ 多序列比对                | §2.4                                 |
| §6.4 避免 all-vs-all 机制（`--sparsify`） | §1.2（pgr 不需要）                   |
| §9.1 隐式图 vs 物化图                     | §2.1（PAF 边集）、§2.5（不物化 GFA） |
| §9.2 PAF/MAF 作为隐式图边集来源           | §2.2（MAF→PAF 转换）                 |
| §9.4 第一步最小原型                       | [[pairwise-selection.md]]            |
| §9.4 "为何 pgr 不需要 --sparsify"         | §1.2                                 |

### cactus.md

| cactus.md 章节                  | 本文档对应                           |
|---------------------------------|--------------------------------------|
| §1.11 Chain/Net↔Flower 对应分析 | §2.1（Chain/Net 是 syntenic 过滤器） |
| §3 Minigraph-Cactus 五阶段分治  | §4.3（pgr partition 参考）           |
| §8 Caf 退火-熔化循环            | §4.2（过滤维度映射到查询层）         |

### cactus_lastz.md

| cactus_lastz.md 章节         | 本文档对应                                            |
|------------------------------|-------------------------------------------------------|
| §5.2 `pgr lav lastz` 能力    | §3.2（lastz 策略，`--self` 仅用于重复屏蔽）           |
| §5.6 完整 lastz 重复屏蔽流程 | §3.2（`--self` 的正确用途：重复屏蔽，非泛基因组比对） |

