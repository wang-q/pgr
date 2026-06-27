# pgr 泛基因组第一步：挑选哪些要做 pairwise 比对

本文档承接 [`impg.md`](impg.md) 的分析，聚焦 pgr 走向泛基因组时面临的**第一个决策**： 在 cohort
场景下，挑选哪些基因组（或区段）对做 pairwise 比对。对应 impg 的 align 阶段 +`--sparsify` 机制
（impg.md §6.4）。

## 1. 问题定义

### 1.1 为什么需要"挑选"

cohort 级泛基因组分析的底层是 pairwise 比对网络。N 个基因组的 all-vs-all 比对复杂度是 O(N²)：

- N=10 → 45 对（可接受）
- N=50 → 1225 对（吃力）
- N=100+ → 4950+ 对（不可行）

impg 的 `--sparsify auto` 用 Mash KNN 把 N² 降到 N×K（impg.md §6.4）。pgr 需要等价机制。

### 1.2 pgr 与 impg 的根本差异

| 维度     | impg                       | pgr                              |
|----------|----------------------------|----------------------------------|
| 比对来源 | 从 FASTA 跑 wfmash/sweepga | 已有两序列 MAF（可转 PAF）       |
| 挑选时机 | align 阶段（无先验）       | 可借已有 MAF 先验                |
| 网络结构 | 直接 all-vs-all            | 已有 pairwise 基础设施           |
| 核心问题 | 选哪些对比对               | 复用已有 pairwise，做 PAF 隐式图 |

**关键认识**：pgr 已有两序列 MAF（可转 PAF），天然避开了 all-vs-all 比对（impg.md §9.4 已论证）。
pgr 的"挑选"问题因此**不是**"选哪些对跑 wfmash"，而是分三层：

## 2. 三层挑选问题

### 2.1 第一层：从已有 MAF/PAF 挑选（查询层，无需新比对）

**问题**：cohort 有 N 个基因组，pgr 已有两序列 MAF（ref↔query_i 或 query_i↔query_j）。
要回答"query_A 的某区段在 cohort 中有哪些同源"，需要在 PAF 网络上跨记录传递。

**机制**：impg 的传递闭包（`-x` BFS，impg.md §4.2）正是为此设计—— 若 A↔B、B↔C 在同一区段有比对，则
A↔C 间接同源。

**PAF 来源**：pgr 已有的两序列 MAF 可通过 `pgr maf to-paf` 转成 PAF（每个 MAF block = 一条 pairwise
alignment，序列化为 PAF 行 with `=`/`X` CIGAR）。这是 pgr 复用已有 pairwise 基础设施 通往 PAF-based
隐式图的天然桥梁（impg.md §9.2）。

**pgr 的挑选参数**（对应 impg `QueryOpts`，impg.md §1.1.3）：

- `--merge-distance`：单跳能吸收的最大 gap/SV 长度
- `--min-identity`：最低 gap-compressed identity
- `--min-output-length`：最短输出区间
- `--max-depth`：BFS 深度（impg 默认 2，即 A↔B↔C 一跳到 C）

**这一层不需要新比对**，只需要把已有 MAF 转成 PAF 装入区间树，做查询层挑选。 对应 impg.md §9.4
的第一步最小原型（`pgr paf query --transitive`）。

### 2.2 第二层：补充 pairwise 比对的挑选（align 层）

**问题**：已有 MAF 只覆盖已跑过的对。以下场景 MAF 缺失或不足：

- cohort 加入新基因组，与已有基因组未跑过 pairwise
- 已有 MAF 在某区段断开（low identity 区段被过滤），但可能有间接同源
- 某些 sample 对需要更精细的 region-level 重比对

**是否需要补充**取决于应用：

- 单 locus 查询（HLA/C4）：第一层传递闭包通常够用
- 全 cohort 泛基因组图构建：需要补充 A↔B 直接比对

**若需补充，挑选策略**（参考 impg §6.4）：

| 策略                | 来源                   | 适用                   | pgr 实现门槛               |
|---------------------|------------------------|------------------------|----------------------------|
| Mash KNN            | impg `--sparsify auto` | 无先验全选             | 需引入 mash crate          |
| 已有 PAF 覆盖度先验 | pgr 独有               | 已有部分 PAF 的 cohort | **推荐**，复用已有 PAF     |
| 系统发育树引导      | Cactus 风格            | 有 phylogeny           | 复用 pgr `nwk` 模块        |
| syncmer anchor      | impg syng 后端         | 免比对                 | pgr 不参考（impg.md §1.1） |

**已有 PAF 覆盖度先验策略**（pgr 推荐）：

1. 对每个 query_i，统计其在已有 PAF 上的覆盖区间集合 C_i
2. 对 query_i、query_j，计算 |C_i ∩ C_j| / |C_i ∪ C_j|（Jaccard）
3. 选 Jaccard 高于阈值且**尚未跑过 pairwise**的对补充比对（`pgr lav lastz`）
4. 已有 `fas multiz`、`psl`、`maf` 工具链可直接处理结果

这样把 N² 降到"PAF 覆盖度共享的子集"，复用 pgr 已有基础设施。

### 2.3 第三层：region 级精细比对挑选

**问题**：已有 MAF 是粗粒度的（受原始比对参数限制）。某些 region（HLA、KIR、C4）需要更精细
pairwise，但全基因组精细比对代价高。

**挑选机制**：

- 从已有 PAF 的 gap/low-identity 区段筛选候选 region
- 对候选 region 跑 `pgr lav lastz`（Cactus 风格，已有）
- 合并回 PAF 网络

**这一层是第一层的补充**，不是泛基因组的核心路径，按需开启。

## 3. 第一步推荐方案

### 3.1 聚焦第一层：PAF 传递闭包

**理由**：

- pgr 已有两序列 MAF（可转 PAF），无需新比对
- 对应 impg.md §9.4 最小原型，验证"隐式图 + 按需物化"闭环
- 单 locus 查询场景（HLA/C4）第一层已足够
- 即使后续需要第二层，第一层也是基础（传递闭包结果是补充比对的优先级依据）

### 3.2 第一步的具体目标

对应 impg.md §9.4，复述要点：

1. **`pgr maf to-paf`** — 新增子命令，两序列 MAF → PAF（`=`/`X` CIGAR）
2. **PAF 索引格式**（`.pgr.paf.idx`）— 按 target seq 建区间树，全量装入不过滤
3. **`pgr paf query <region>`** — 区间投影，暴露过滤参数（`--merge-distance`/`--min-identity`/
   `--min-output-length`/`--subset-sequence-list`）
4. **`pgr paf query <region> --transitive`** — 传递闭包 BFS，`--max-depth`/`--max-gap` 控制
5. **`pgr paf query <region> --transitive -o maf`** — 传递闭包 + 局部 MSA（复用 `fas_multiz.rs`）

### 3.3 暂不做

- **第二层补充 pairwise**：待第一层验证后，根据"传递闭包覆盖率是否够用"决定
- **Mash KNN**：pgr 有 MAF 先验，不必引入 mash 依赖
- **syng 后端**：impg.md §1.1 已明确不参考
- **stage DSL**：第一层是单命令，不需要管道化

## 4. 待澄清的问题

1. **MAF → PAF 的 CIGAR 精度**：两序列 MAF 的每个 block 对应一段连续 match/mismatch/indel。 转
   PAF 时需把 MAF block 的逐 base 对齐展开成 `=`/`X`/`I`/`D` CIGAR。需确认 pgr 现有 `maf`模块
   （目前只有 `to_fas`）能否支持，或需新增 `maf to-paf` 子命令。
2. **多 PAF 文件的统一索引**：cohort 可能有多个 PAF 文件（每对一个），如何装入单一区间树？ impg 的
   `MultiImpg`（impg.md §3.4）协调 per-file 子索引，pgr 可参考。
3. **MAF 的链向信息**：MAF block 内的 strand 字段（+/-）对应 PAF 的 strand。需确认 pgr MAF
   解析正确处理反向互补 block（特别是 inversion 场景）。
4. **PAF 质量过滤**：PAF 无 score 字段（不像 Chain）。impg 的做法是"索引层不过滤，查询层用
   `--min-result-identity` 过滤"（impg.md §1.1.3）。pgr 的 PAF 由 MAF 转来，identity 可从 CIGAR
   计算，查询层过滤即可。

## 5. 与 impg.md 的对照

| impg.md 章节                      | 本文档对应                    |
|-----------------------------------|-------------------------------|
| §1.1.3 能力栈四层                 | §1.2 pgr 与 impg 的根本差异   |
| §1.1.4 名词解释 pair-selection    | §1.1 为什么需要"挑选"         |
| §4.2 传递闭包                     | §2.1 第一层 PAF 传递闭包      |
| §6.4 避免 all-vs-all 机制         | §2.2 第二层补充 pairwise 挑选 |
| §9.2 PAF/MAF 作为隐式图边集       | §2.1 PAF 来源                 |
| §9.4 第一步最小原型               | §3 第一步推荐方案             |
| §9.4 "为何 pgr 不需要 --sparsify" | §2.2 已有 PAF 覆盖度先验策略  |

