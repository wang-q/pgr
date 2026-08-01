# 4 万大肠杆菌基因组泛基因组路线

本文档记录 pgr 处理 4 万个大肠杆菌（*E. coli*）基因组的**完整路线**。它从原始 cohort 出发，
最终到达隐式图查询与可选的 GFA 物化。前段（去冗余）是 pgr 的**上游预处理**，后段（PAF 索引、
查询、GFA）已在 [[paf-pangenome.md]] 规划，本文档把它们串成一条端到端 pipeline。

参考文档：[[paf-pangenome.md]]（PAF 隐式图核心目标、路线、query / to-maf / graph / to-gfa / to-vcf / stat 实现）、[[impg.md]]（impg 的
sparsify 与去冗余差异）、[[minigraph.md]]（粗框架哲学）。

---
## 1. 端到端流程概览

```
4 万 E. coli 基因组
    │
    │  [步骤 1: 去冗余] —— Mash 近邻 + FastGA 精比 + 1e-5 阈值
    ↓
~27000 非冗余基因组
    │
    │  [步骤 2: Mash KNN sparsify] —— 27000² 降到 27000×K（见 §2.1）
    ↓
pairwise 比对 (MAF/PAF)
    │
    │  [步骤 3: PAF 索引] —— pgr paf index（已实现）
    ↓
PAF 索引 (隐式粗图)
    │
    │  [步骤 4: 查询] —— pgr paf query --transitive（query / to-bed / to-maf）
    ↓
同源区段列表 (BED/PAF/FASTA/MAF)
    │
    │  [步骤 5: GFA 物化] —— graph 粗全局 + to-gfa 区域精细（可选）
    ↓
rGFA（粗全局）/ 局部 GFA（区域）→ MAF/VCF
```

**关键边界**：步骤 1-2 是 pgr 的**上游**，步骤 3-5 是 pgr 的**核心**。本文档重点记录步骤 1
（用户已明确），步骤 2-5 引用已有规划。

> MAF 是上游工具（FastGA/lastz）的**输出格式**，在 pgr 入口处转成内部格式：
> 隐式图路径 `maf to-paf` → `paf index`，多序列比对/核心基因组路径 `maf to-fas`
> → `pgr pl p2m`（`fas` 层 cover/slice/join/refine）。内部链路不使用 MAF。

---
## 2. 步骤 1：去冗余（4 万 → ~27000）

### 2.1 目标

从 4 万个原始基因组中找出 **non-redundant** 集合。判据：**全基因组核苷酸差异 < 1e-5**
（即 0.001%，约每 10 万 bp 1 个差异）的基因组视为冗余，只保留一个代表。

实测结果：4 万 → **~27000** 非冗余基因组（去冗余率约 32.5%）。

### 2.2 两阶段方法

去冗余分两个阶段，对应不同的精度/速度权衡：

#### 阶段 1a：Mash 近邻（alignment-free，快速筛）

- **工具**：Mash（external，[mash sketch](https://mash.readthedocs.io/) + `mash dist`）
- **原理**：k-mer sketch 算基因组间距离，无需比对
- **作用**：在 4 万² ≈ 8 亿对中，为每个基因组找出 K 个最近邻，形成 KNN graph
- **判据**：Mash distance 低于阈值（待定）的对进入阶段 1b 精比
- **输出**：候选冗余对列表（远小于 8 亿对）

**为什么需要 Mash**：直接对 4 万² 对跑 FastGA 不可行。Mash 是 alignment-free 的，4 万基因组
sketch + 全对距离矩阵可在数小时内完成，把 N² 降到可管理的规模。

#### 阶段 1b：FastGA 精比（alignment-based，精确判冗余）

- **工具**：FastGA（external，专为细菌/小基因组优化的比对工具）
- **原理**：对 Mash 找出的近邻对跑完整全基因组比对
- **作用**：算出近邻对的**真实核苷酸差异率**
- **判据**：全基因组差异 < 1e-5 → 判为冗余
- **输出**：冗余对列表 + 每对的差异率

### 2.3 代表基因组选取

判定冗余对后，需要为每个冗余 cluster 选一个**代表基因组**：

- 用 `pgr clust` 的 k_medoids（已实现）可选 cluster 中心作为代表
- 或按用户指定的优先级（如 NCBI RefSeq 优先、组装质量优先）

### 2.4 与 impg sparsify 的区别

**这是两个不同的层次**，不能混淆（见 [[impg.md]] §6.4）：

| 层次 | 工具 | 作用 | 输入 | 输出 |
|------|------|------|------|------|
| **去冗余**（本文 §1）| Mash + FastGA | 哪些基因组可以代表 cluster | 4 万基因组 | ~27000 非冗余 |
| **pair-selection**（impg sparsify）| Mash KNN | 哪些 pair 跑比对 | N 个非冗余基因组 | 候选比对对列表 |

impg 的 `--sparsify` 假设输入已是非冗余的，它在 N² 对里选 K 近邻边跑 wfmash。pgr 的步骤 1
是 impg 假设的**前置条件**——先把 4 万去冗余到 27000，再考虑 27000² 里选对。

### 2.5 pgr 的角色

步骤 1 采用**方案 A（外部预处理）**：Mash + FastGA 在 pgr 外跑，pgr 只接收非冗余基因组列表。
这与 [[paf-pangenome.md]] §1 的"不重新做比对"原则一致——pgr 的核心价值在步骤 3-5
（PAF 索引 + 隐式图），步骤 1 是上游数据准备。

步骤 2 同理：Mash KNN sparsify + FastGA 比对在 pgr 外跑，pgr 只接收产出的 PAF 做索引。
FastGA 贯穿步骤 1b（去冗余精比）和步骤 2c（sparsify 比对），保持工具一致性。

未来若要端到端封装，可演进到方案 B（`pgr pl dedup` + `pgr pl sparsify` 调用 Mash + FastGA），
但当前不优先。

| 方案 | 含义 | 状态 |
|------|------|------|
| **A. 外部预处理** | Mash + FastGA 在 pgr 外跑；pgr 接收非冗余基因组（步骤 1）和 pairwise PAF（步骤 2c） | ✅ 当前采用 |
| **B. pgr 封装** | pgr 新增 `pgr pl dedup` / `pgr pl sparsify` 命令，内部调用 Mash + FastGA | 远期 |
| **C. pgr 原生** | 用 Rust 重写 Mash sketch + 用 lastz 替代 FastGA | 不推荐（lastz 比 FastGA 慢）|

---
## 3. 后续步骤

步骤 2-5 的细节引用已有文档，本文档不重复。

### 3.1 步骤 2：Mash KNN sparsify（27000² 降到 27000×K）

**这是隐式图架构避免 N² 爆炸的核心机制**。隐式图的价值就是"稀疏比对 + 传递闭包推断全量"：
A↔B、B↔C 已比对，查询 C 的同源区段时 BFS 经 B 到达 A——即使 A↔C 从未直接比对。
因此步骤 2 **不需要**全量 all-vs-all，只需稀疏覆盖。

| 子步骤 | 做什么 | 工具 |
|--------|--------|------|
| 2a | Mash sketch 全部 27000 基因组，算 N² 距离矩阵 | Mash（alignment-free，快）|
| 2b | 每个基因组取 K 个最近邻，形成 KNN graph | Mash KNN（impg `--sparsify auto` 同款）|
| 2c | 对 K×N 条边跑 FastGA，产 MAF → 转 PAF | FastGA（external，细菌基因组优化，比 lastz 快）|
| 2d | PAF 索引（步骤 3，已实现）| `pgr paf index` |

**规模估算**：K=50 时约 135 万对 FastGA，比全量 3.6 亿对缩减 ~270 倍。

**为什么是 sparsify 而非全量**：见 [[paf-pangenome.md]] §1.2——大 cohort + 无 MAF 先验时，
sparsify 是必需的（不是可选的）。传递闭包负责推断 sparsify 遗漏的对。

**为什么用 FastGA 而非 lastz**：FastGA 专为细菌/小基因组优化，比 lastz 快；135 万对的规模下
速度差异显著。lastz 是 pgr 已有的 fallback，但步骤 2c 首选 FastGA。这与步骤 1b（去冗余精比）
用 FastGA 一致——同一工具贯穿去冗余和 sparsify 比对两个阶段。

**为什么不借鉴 seqwish 的 `--sparse-factor` 哈希稀疏化**：seqwish 在
[alignments.rs](../../seqwish-master/src/alignments.rs) 的 `keep_sparse`
用哈希函数随机丢弃对齐，是廉价的随机稀疏化。pgr 走 Mash KNN 稀疏化，按生物学相似性选近邻边，
质量更高，能保证每个基因组的 K 个最近邻都被覆盖。详见 [[seqwish.md]] §8。

**待调参数**（不是方向问题，是调优问题）：

- K 值：K 越大覆盖越好但计算量线性增长，K=50 是 impg 默认，需按 E. coli 实际多样性调整
- 传递闭包覆盖率：跑完后用 `pgr paf query --transitive` 抽样验证 BFS 能否到达预期对
- FastGA 耗时：135 万对 FastGA 的实际耗时需测，可能需并行化（FastGA 本身支持并行）

### 3.2 步骤 3：PAF 索引（已实现）

- `pgr paf index`（见 [[paf-pangenome.md]]）
- 输入：步骤 2 产出的 pairwise MAF/PAF
- 输出：PAF 索引（隐式粗图）

### 3.3 步骤 4：查询（query / to-bed / to-maf 已实现）

- `pgr paf query` / `to-bed` / `to-maf`（见 [[paf-pangenome.md]] §3）
- query ✅: PAF 默认输出 + `-b` 批查 + `to-bed` 独立子命令
- to-maf ✅: `to-maf` pairwise MAF（按 CIGAR 还原，需 `-f TSV`）
- to-maf --msa ✅: `to-maf --msa` POA 多序列 MSA（需 `--transitive`）

### 3.4 步骤 5：GFA 物化与 VCF（graph / to-gfa / to-vcf 已实现）

- graph ✅: `pgr paf graph [-f refs.fa] --min-var-len 100` 粗全局 GFA（seqwish DSU 风格）；`-f` 可选，拓扑模式零序列依赖
- to-gfa ✅: `pgr paf to-gfa` 区域精细 GFA（impg 风格，unchop 默认开，`--crush` 可选）
- to-vcf ✅: `pgr paf to-vcf` POA MSA → VCF（SNP + INS/DEL，1bp anchor）

**graph 算法骨架**：seqwish 风格段级 DSU（CIGAR 切分 → 段对 → DSU 传递闭包 → 节点序列 → 路径
+ novel 段补全 → 边派生 → GFA 输出），详见 [[seqwish.md]] §7.2。pgr 输入 PAF 与 seqwish 一致，
相对 minigraph（需自跑 minimizer chaining）更天然适配。

**实现简化项**（相对 seqwish）：无 disk-backed interval tree、无 SparseBitVec、无 lock-free DSU
（单线程足够），路径方向恒 `+`（反向已翻转坐标到正链），rGFA SN/SO/SR tag 已补全（见
[[paf-pangenome.md]] §3.3）。

**4 万大肠杆菌规模超 RAM 的兜底**：seqwish 的 `AdaptiveTree` 提供磁盘/内存双后端区间树
（[intervaltree.rs](../../seqwish-master/src/intervaltree.rs)），数据
落盘 mmap，内存吃紧时切到 disk-backed 模式，牺牲速度换可跑性。pgr paf graph 处理 4 万大肠杆菌全图物化
（Gbp 级）时可考虑引入此兜底机制。详见 [[seqwish.md]] §3.3、§7.2。

---
## 4. 小 cohort 先行验证（3 → 10 基因组，2026-08-02）

在 4 万数据可用前，用 MG1655 × Sakai × SE11（11 replicon，数据见 [[ecoli-genome.md]]）
跑通端到端管线，验证 `pgr paf` 在真实数据上的行为并固化回归。**流程（用户指定）**：

```
FastGA -psl（3 对）→ pgr pl chainnet --syn（每对）→ maf to-paf → paf index
→ query --transitive → graph / stat
```

要点：FastGA 产出 PSL 后**必须过 `pgr pl chainnet --syn`**（syntenic 过滤），
不直接用原始 PAF。完整命令见 `scripts/verify-pangenome.sh`。

**结果**：chainnet --syn 把 10 条序列过滤到 5 条（只留有 syntenic 链的 replicon，
如 sakai 的 NC_002127 被滤掉）；`query --transitive` 从 `mg1655.NC_000913:100000-110000`
出发，BFS 深度 2 找到 Sakai（`sakai.NC_002695`）与 SE11（`se11.NC_011415`）染色体及两者
间中转连接——三方同源成立；`graph`/`stat`/`to-vcf`/`to-maf --msa`/`to-gfa`/`to-fas`
全部跑通。回归脚本一键 PASS。

**踩坑记录**（4 万扩展同样适用）：

- chainnet 的 MAF 序列名带文件名前缀（`mg1655.NC_000913`），查询 region 必须用带
  前缀的名字。
- `paf query` 只认 `.paf.idx` 后缀的索引（`ends_with(".paf.idx")`），其它后缀会被当
  PAF 文件解析报 UTF-8 错。
- 不同 pair 的 MAF 可能同名（都是 `NC_000913.maf`），转 PAF 时输出名必须按 pair
  区分，否则互相覆盖。
- 小 replicon（质粒）的弱链会被 chainnet --syn 过滤，小 cohort 的图以染色体为主。
- **PAF 方向要保证 seed 在 target 侧**：`paf query` 直接查询只查 `trees`
  （seed 作为 target 的记录），`--transitive` 的 BFS 额外走 `reverse_trees`
  （镜像），但镜像**只为 '+' 链记录创建**。若 seed 落在 '-' 链记录的 query 侧，
  BFS 无法跨过它（实测：mg1655 作 query 时 ec2011c_3493 的覆盖记录是 '-' 链，
  传递查询到达不了）。因此脚本统一 `FastGA(b,a) → chainnet(a,b)`，让 mg1655
  （种子所在序列）恒为每对的 target，深度 1 即可直达所有株。这是 pgr 的一个
  已知缺口：`paf index` 应对 '-' 记录也建镜像（需处理反向 CIGAR 语义）。
- 扩展脚本曾把 PAF 输出名写成裸 basename（`NC_004431.paf`），同一 replicon 在
  多个 pair 里作 target 时被互相覆盖，导致 mg1655 侧的记录丢失、BFS 够不到对应
  菌株——必须带 pair 前缀（`$pa-$pb-*.paf`）。

**对 4 万扩展**：脚本遍历 chainnet 输出目录内所有 maf（不硬编码染色体），把 PAF 来源
换成步骤 1-2 的去冗余 + sparsify 产出即可直接复用；规模相关待调参数见 §3.1。

**cohort 扩展（3 → 10）**：数据侧已预选 7 株覆盖主要 pathotype 的典型菌株
（CFT073 / E2348/69 / 042 / 2011C-3493 / E24377A / EC958 / Nissle 1917，
下载命令见 [[ecoli-genome.md]] "Typical strains"），10 基因组验证**已完成并 PASS**
（2026-08-02，45 对全跑，约 7 分钟）：`mg1655.NC_000913:100000-110000` 的传递查询
深度 1 直达全部 10 株染色体（含多质粒菌株 E24377A / 2011C-3493 / EC958 /
E2348/69 / 042，质粒弱链被 --syn 滤掉、图以染色体为主）；graph / stat 一并通过。
`scripts/verify-pangenome.sh` 已按 seed 在 target 侧 + pair 前缀命名固化。

**10 基因组下游验证（to-maf --msa / to-vcf / to-fas）**：区域查询在 10 基因组规模
暴露了 MSA 输入爆炸——传递 BFS 对 10kb 区域返回 73 条结果（同一条序列经不同路径
反复到达），`build_msa_entries` 不去重就把 ~90 条序列喂给 POA，10kb 区域 release
跑 5 分钟不结束。修复：同名同链结果先精确去重、再合并重叠区间（排序键
`(name, strand, qs)`，否则 '-' 链区间会插在两条 '+' 区间中间挡住合并），输入降到
每株 1-2 条。修复后 10kb MSA / VCF / block FA 均在 ~54s 完成（release），MSA 行
序列与源文件逐字一致；VCF 4875 行。剩余瓶颈是 POA 引擎本身（O(节点数×序列长)，
progressive 逐条加入），属 [[fas-multiz.md]] 已记录的 profile-profile DP 优化方向。
回归测试：`command_paf_to_maf_msa_dedup_transitive_duplicates`（菱形拓扑 → 3 行）。

**EC958 数据注意**：必须用 GCF_000285655.3 完整版（EC958.v1，chr + 2 质粒）；
`.2` 是 WGS scaffold（240 条 contig），不要用。

---
## 5. 与 pgr 核心目标的关系

回到 [[paf-pangenome.md]] 开篇的核心目标——"复用 pairwise 资产，构建 PAF 隐式图"——4 万 E. coli
是这条路线的**应用场景**，不是目标本身。步骤 1（去冗余）是场景的前置条件，步骤 2（pair-selection）
是场景的规模约束，步骤 3-5 是 pgr 的核心能力。

**优先级**：pgr 的开发资源应集中在步骤 3-5（query / to-maf / graph / to-gfa / to-vcf 路线）。步骤 1 作为外部预处理，待 pgr 核心
能力稳定后再考虑是否封装（方案 B/C）。

---
## 6. 变更日志

- 2026-08-02：§4 完成 10 基因组验证（45 对全跑 PASS）——cohort 从 3 扩到 10，
  脚本固化 FastGA(b,a) → chainnet(a,b) 方向与 pair 前缀 PAF 命名。
- 2026-08-02：§4 补充下游验证与 MSA 输入修复——传递 BFS 重复结果导致 POA 爆炸，
  `build_msa_entries` 增加同名同链去重 + 重叠区间合并（10kb MSA/VCF/FAS 从 >5min
  降到 ~54s），新增菱形拓扑回归测试。
- 2026-08-02：§4 补充 cohort 扩展说明（3 → 10 基因组）——预选 7 株典型 pathotype
  菌株（CFT073 / E2348/69 / 042 / 2011C-3493 / E24377A / EC958 / Nissle 1917），
  下载命令见 [[ecoli-genome.md]] "Typical strains"。
- 2026-08-02：新增 §4 小 cohort 先行验证（3 基因组）——吸收原独立文档
  `ecoli-pangenome-3way.md`（已删除），记录 FastGA → chainnet --syn 流程的执行结果
  与踩坑。

- 2026-06-28：初稿。记录步骤 1（去冗余）的 Mash + FastGA + 1e-5 方法（已落地，4 万→~27000），
  步骤 2-5 引用已有文档。
- 2026-06-28：步骤 1 已由用户实际跑通，删除原 §2 开放问题章节，章节顺次前移。
- 2026-06-28：步骤 2c 明确用 FastGA（贯穿步骤 1b 和 2c，保持工具一致），不用 lastz。
- 2026-06-28：基于 [[seqwish.md]] 分析，补充 §2.1 不借鉴哈希稀疏化的理由，§2.4 graph 算法骨架
  与磁盘后端兜底。
