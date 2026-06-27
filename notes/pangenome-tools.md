# 泛基因组工具生态与 pgr 的定位

> 整理于 2026-06-28，源自对当前泛基因组图构建工具的网调。 目的：把 pgr 放进工具生态里对照，
> 明确差异化与护城河。

## 0. 为什么要梳理

pgr 的隐式图路线借鉴自 impg，但 impg 并非唯一参考。当前泛基因组图工具按**构建思想**可分为 5+ 类，
每类背后是一种关于"如何从基因组集合得到图"的回答。pgr 的选择（PAF 索引 + 隐式 interval tree +
BFS 传递闭包）在这个谱系里有明确位置，也有明确的竞争对手。本文梳理全景，并在 §3 给出对 pgr
路线的具体启示。

## 1. 按构建思想分类

### 1.1 参考 + VCF 变异增强（早期范式）

把已知的 VCF 变异 threading 到线性参考上，形成增强图。

- **vg (Variation Graph toolkit)** — 第一个综合开源框架。双向循环图，支持 SV/倒位。依赖外部变异
  caller（minimap2+Sniffles、MUMmer+SVMU 等），常用 SURVIVOR 合并 VCF。
- **GrapTyper** — 嵌入式图，对短读长迭代重比对，专注小变异 calling。
- **Graph Genome pipeline (Rakocevic 2019)** — 商业，仅人类，不开放。

**特点**: 成熟、生态完整，但**依赖参考 + VCF**，受参考偏倚和 VCF 表达力限制（嵌套 SV、
复杂重排表达弱）。

### 1.2 参考 + 直接构建（无需 VCF）

输入参考 + 全部基因组，直接输出图，跳过 VCF 中间步。

- **Minigraph** — 扩展 minimap2 的 chaining，渐进式把 > 50bp SV 加进图。3 小时跑 20 个人类基因组，
  **不识别 SNP**，只保留大 SV 的"粗框架"。
- **Minigraph-Cactus** — Minigraph 粗框架 + Cactus 碱基级 MSA。HPRC 三大工具之一。

**特点**: 比 vg 简单，但**仍需参考 + 输入顺序**，有顺序偏倚。Minigraph 的"≥100bp 粗框架过滤"哲学是
pgr V4a 的直接借鉴来源（见 [[minigraph.md]]）。

### 1.3 无参考 all-vs-all 比对（无偏范式）

不假设参考、不假设系统发育，从 all-vs-all 比对诱导图。

- **PGGB** — `wfmash` (all-vs-all) → `seqwish` (诱导图) → `smoothxg`+`GFAffix` (归一化)。
    - **wfmash**: MashMap 找同源块 + wflign (WFA) 碱基级比对
    - **seqwish**: **用 implicit interval tree**把 PAF 物化成变异图
    - **smoothxg**: 消除 all-vs-all 产生的复杂 looping motif，归一化
    - **无偏**: 从 SNP 到 SV 全尺度，不依赖输入顺序
- **Progressive Cactus** — 渐进式，用 NJ 树引导，无单一参考但依赖树拓扑（见 [[cactus.md]]）。

**特点**: 最无偏、碱基级最完整，但**O(N²) all-vs-all**限制规模。HPRC 90 单倍型已是上限。

### 1.4 de Bruijn 图（k-mer 范式，alignment-free）

固定 k 长度，把基因组拆成 k-mer 建 cdBG，避免比对。

- **TwoPaCo** — Bloom filter + hash table 直接构建 compacted de Bruijn 图，100 个人类基因组/天。
- **Bifrost** — 高度并行构建 colored cdBG，识别 junction k-mers，无需中间 uncompacted 图。
- **mdbg** — minimizer-space de Bruijn 图，**可扩展到 60 万细菌基因组**。
- **PanTools** — generalized DBG 存入 Neo4j 图数据库，层次化（基因组 + 注释 + 同源群），
  从病毒到人类都能跑。

**特点**: 极致可扩展，但**固定 k 长度**丢失结构信息，难以做 SV 分析。适合变异检测/基因分型，
不适合结构变异研究。

### 1.5 基因内容图（gene-level，非碱基级）

节点是基因而非碱基段，关注基因有无与顺序。

- **Pangene (李恒 2024)** — 蛋白序列 → miniprot 比对 → 基因图（节点 = 基因，边 =
  邻接），每个基因组是图中一条 walk。`bibubbles` 捕获基因拷贝数/顺序/方向变化。
  填补了人类泛基因组基因级分析的空白。
- **PANSEQ / PGAP / Roary / Panaroo** — 细菌泛基因组传统工具，BLAST 聚类核心/附属基因。

**特点**: 关注**基因有无**而非碱基变异，适合细菌但丢失位置细节。Pangene 是 pgr 的上层互补（不同尺度）。


### 1.6 索引/可视化/工具链（非构建工具）

- **odgi** — 动态图接口，排序/剪枝/可视化大图（Gbp 级）
- **xg** — succinct 静态索引（节点/边/路径），实现 libhandlegraph API
- **GBWT** — 基于 PBWT 的路径索引，1 bit/kbp 压缩 1000G
- **PGVF** — GFAv1 硬 fork，支持图对图比对
- **pgge** — 泛基因组图评估 pipeline

## 2. 五类范式横向对比

| 范式       | 代表工具               | 参考依赖 | 规模上限   | 变异尺度  | 主要用途               |
|------------|------------------------|----------|------------|-----------|------------------------|
| 参考+VCF   | vg, GrapTyper          | 强       | 数十       | SNP~SV    | 变异 calling、基因分型 |
| 参考+直接  | Minigraph, MC          | 中       | 数十~百    | ≥50bp SV  | 粗框架图、可视化       |
| all-vs-all | PGGB, ProgCactus       | 无       | ~百        | SNP~SV    | 无偏泛基因组参考       |
| de Bruijn  | TwoPaCo, mdbg, Bifrost | 无       | **数十万** | SNP/indel | 大规模索引、比对       |
| 基因内容   | Pangene, Roary         | 无       | 数千~万    | 基因 PAV  | 基因级群体分析         |

## 3. 与 pgr 的关系

### 3.1 思想同源：seqwish

**seqwish 是 pgr 最直接的思想同源工具**。两者都用**implicit interval tree**从 PAF 处理比对：

- **seqwish**: 把 PAF **物化**成 GFA 变异图，目标是产出可持久化的图文件
- **pgr**: 保持**隐式**，按需查询/物化，目标是避免全局图爆炸

差异化的边界：pgr 不产出全局 GFA，只在查询时走 BFS 传递闭包，在 V4b 按需产出局部 GFA。
这一点 [[graph-design.md]] §4.3.4 已明确"局部 GFA 不合并回全局"。

**粒度差异**是核心：seqwish 的传递闭包是**全局、一次性**的——一次性算出全部等价类再写图
（[transclosure.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/transclosure.rs) 的
spanning tree → BFS discovery → DSU union-find 流程）；pgr 是**局部、按需**的——每次查询从一个
区间出发 BFS，只算相关等价类。这两种粒度对应不同应用场景：全图统计 vs 单点查询。详见
[[seqwish.md]] §5 的对照表与 §6 的版本启示。

### 3.2 反例：PGGB 的 O(N²) 瓶颈

PGGB 的 all-vs-all 是 pgr 的**反面教材**：

- PGGB: 90 单倍型已是 HPRC 上限，all-vs-all 比对 + smoothxg 归一化成本极高
- pgr: 用**Mash KNN sparsify**把 N² 降到 N×K，靠**传递闭包**补全稀疏比对的缺口

这条路径在 [[paf-route.md]] §1.2 已确立，[[ecoli-cohort.md]] 给出了 4 万 E. coli 的具体落地（27000²
→ 27000×50）。

### 3.3 护城河：相对 mdbg 的碱基级优势

mdbg 可扩展到 60 万细菌基因组，是 pgr 在大规模场景的潜在竞争对手。pgr 的护城河：

- **保 SV**: mdbg 固定 k 长度丢失结构信息，pgr 基于 PAF 保留碱基级比对
- **保坐标投影**: pgr 支持任意基因组间的坐标互查，mdbg 无此能力
- **保传递性**: pgr 的 BFS 传递闭包可推断间接同源，mDBG 的 k-mer 合并无法区分直接/间接

代价：pgr 的扩展性不如 mdbg。若未来场景扩到百万级细菌基因组，可能需要 k-mer 范式作为预处理层（先用
mdbg 聚类，再对子集跑 PAF）。

### 3.4 互补：Pangene 的基因级视图

Pangene 与 pgr 在不同尺度：

- **pgr**: 碱基级，关注 SV/indel/SNP
- **Pangene**: 基因级，关注基因 PAV/拷贝数/顺序

两者互补而非竞争。pgr V5+ 可考虑把 Pangene 作为上层接口——查询某基因时，先用 Pangene 定位子图，再用
pgr 查碱基级细节。当前路线不纳入，但记录为远期可能。

## 4. 对 pgr 路线的启示

1. **seqwish 对照已完成** — [[seqwish.md]] 详细拆解了 seqwish 的 6 阶段流程（序列索引 → 比对索引 →
   传递闭包 → 节点压缩 → 边派生 → GFA 输出），明确 pgr V4a 可直接复用其算法骨架（spanning tree →
   BFS → DSU → compact → links → GFA），V1 可借鉴 `PosT` 编码与 `SparseBitVec`，V4b 可借鉴
   orphan recovery，查询层可考虑预计算生成树优化 BFS。
2. **PGGB 是 sparsify 动机的最佳反例** — [[paf-route.md]] §1.2 已写大 cohort 场景，可把 PGGB 的
   O(N²) 作为显式对照。
3. **mdbg 定义了扩展性上限** — pgr 在 4 万 E. coli 场景护城河清晰，但百万级需警惕，可考虑 k-mer
   预处理作为远期 fallback。
4. **Pangene 是上层互补** — 不纳入当前 V1-V5，但作为 V6+ 候选记录。

## 5. 参考文献

- Garrison et al. *Pangenome graphs*. Annu Rev Genomics Hum Genet 2020.
- Garrison & Guarracino. *Unbiased pangenome graphs (seqwish)*. Bioinformatics 2023.
- Andreace et al. *Comparing methods for constructing and representing human pangenome graphs*.
  Genome Biology 2023.
- Li. *Minigraph*. Nat Commun 2020.
- Hickey et al. *Pangenome graph construction from genome assemblies with Minigraph-Cactus*. Nat
  Biotechnol 2023.
- Guarracino et al. *PGGB*. bioRxiv 2023.
- Minkin et al. *TwoPaCo*. 2016.
- Holley & Melsted. *Bifrost*. Genome Biology 2020.
- Ekim et al. *mdbg*. Nat Biotechnol 2021.
- Li, Marin, Farhat. *Exploring gene content with pangene graphs*. Bioinformatics 2024.
- Jonkheer et al. *PanTools*. Bioinformatics 2016.

