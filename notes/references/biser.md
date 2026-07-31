# BISER 源码与论文分析

> 整理于 2026-07，源自对 `biser-master/` 目录源码及 published paper 的通读。目的：理解 BISER 在
> segmental duplication (SD) 检测与分解中的算法设计，并为 pgr 中重复/同源区域分析提供参考。

## 1. BISER 概览

### 1.1 工具定位

- **工具名称**: BISER (Brisk Inference of Segmental duplication Evolutionary stRucture)，版本 v1.4。
- **功能**: 在单个或多个基因组中快速检测 segmental duplications (SDs)，并将检测到的 SDs 分解为
  elementary SDs 与 core duplicons。
- **核心假设**: 输入基因组应预先 soft-masked（低复杂度/重复序列用小写表示）；BISER
  会在内部将其转换为 hard-masked 序列进行分析。
- **与 SEDEF 的关系**: BISER 是 SEDEF 的继任者，继承了其 SD error model，但在算法上改用线性
  plane-sweep 替代 MinHash，从而获得数倍加速。

### 1.2 输入输出

**命令行接口（由 Python wrapper 提供）**:

```bash
biser -o <output> -t <threads> <genome1.fa> [genome2.fa ...]
```

主要阶段（`biser/__main__.py`）:

1. `mask`: 若输入非 hard-masked，将 lowercase bases 过滤，生成 hard-masked 基因组。
2. `search`: 对每条染色体做 putative SD detection。
3. `align`: 对 putative SD pairs 做局部比对与边界精修。
4. `cross_search` / `cross_align`: 多基因组时，将每个基因组的 SD 映射到其他基因组。
5. `cluster` / `decompose`: 聚类重叠 SD 并分解为 elementary SDs。
6. `translate`: 将 hard-masked 坐标映射回原基因组坐标。

**输出格式**:

- 主输出为 BEDPE 格式，描述 SD mate 的坐标、链向、CIGAR、error rate 等。
- `.elem` 文件记录 elementary SD decomposition 结果，包括 core duplicon 标记。

## 2. BISER 论文算法

BISER 发表于*Algorithms for Molecular Biology* (2022) 与*WABI 2021*。论文将 SD
分析问题拆成三个子问题：

1. **SD detection**: 在单个基因组内找到所有合法 SD pairs。
2. **Cross-species conservation detection**: 将一个基因组的 SD copies 映射到其他相关基因组。
3. **SD decomposition**: 将 SDs 分解为 elementary SDs，并识别 core duplicons。

### 2.1 SD 定义与错误模型

**定义 1** (segmental duplication): 给定错误阈值 ε，SD 是一对 paralog 序列 (G_i, G_j)，满足：

1. `err(G_i, G_j) ≤ ε`
2. 最优比对长度 `ℓ ≥ 1000`
3. 两个 paralog 之间的重叠不超过 `ε · n`

其中 `err(s, s') = E(s, s') / ℓ` 为编辑错误率。

**SD error model**: BISER 假设 SD 的突变由两个独立过程叠加而成：

- **PSV (paralogous sequence variants)**: 背景点突变，随机分布。在 75% 同源性时贡献约 ≤15% 错误。
- **Block edits**: 大规模块插入/删除/重排，非随机分布。在 75% 同源性时贡献约 10% 错误。

因此总错误率 `ε = ε_P + ε_B`。只有 PSV 部分适合用 Poisson 模型建模。

### 2.2 Putative SD Detection：有序 Jaccard + Plane Sweep

由于全基因组局部比对的二次复杂度不可行，BISER 先用 k-mer 相似性做过滤，快速得到*putative SD*区域对，
然后再做精确比对。

**Jaccard 下界**: 在 SD error model 下，两个 paralog 序列的期望 Jaccard index 满足：

```
τ ≥ (1 - ε_B) / (1 + ε_B) · 1 / (2e^(k·ε_P) - 1)
```

**有序 Jaccard (ordered Jaccard)**: 经典 Jaccard 允许任意 k-mer 交叉匹配，可能引入 edit distance
下不可能的交叉比对。BISER 引入 `s ⊛ s'`，表示 `s` 与 `s'` 之间最大*colinear* k-mer matching 的大小；
并定义有序 Jaccard：

```
Ĵ(s, s') = (s ⊛ s') / |K(s) ∪ K(s')|
```

**Lemma 1**: 在 SD error model 假设下（共享 k-mer 在 copy 事件前已共享），有序 Jaccard 等于经典
Jaccard。因此可用 `Ĵ ≥ τ` 作为过滤条件。

**Plane-sweep 算法**（图 1）:

1. 构建全基因组 k-mer 索引 `I_G`。
2. 从左到右扫描基因组，维护一个已发现 putative SD 的有序列表 `L`。
3. 在位置 `x` 处，查询 `I_G` 得到当前 k-mer 的所有出现位置 `K`。
4. 对 `K` 中每个 `y`，判断它：
    - (1) 开启一个新的 putative SD；
    - (2) 延伸 `L` 中已有的 putative SD；
    - (3) 被已有 putative SD 覆盖。
5. 若 `L` 中某个 SD 与当前 `y` 距离过远，则将其提升为最终 putative SD（若满足 `Ĵ ≥ τ`）。

该算法在实践中近似线性，因为 `|L|` 受距离阈值限制，`K` 中高频 k-mer 也会被过滤。

**Winnowing**: 为加速索引构建，BISER 不索引所有 k-mer，而是对每个大小为 `w` 的滑动窗口取字典序最小的
k-mer（tie 时取最右）。期望指纹大小为 `2|G|/(w+1)`。论文/代码默认 `k=14`, `w=16`。

### 2.3 局部比对：Seed-and-Extend + Chaining + Refinement

对每对 putative SD，BISER 执行：

1. **Generate anchors**: 在 putative SD 的两个 mate 之间找 10-mer 精确匹配锚点，并做简单延伸。
2. **Chaining**: 使用 Priority Search Tree (PST) 在 `O(n log n)` 时间内找到得分最高的锚点链。
   该链给出 SD 的粗略边界。
3. **Refinement**: 在锚点链之间的 gap 上执行稀疏动态规划（sparse DP），精修边界与 CIGAR。代码中通过
   `align.refine` 里的 DP 合并 anchors，并用 `bio.seq.align`（SIMD/分块 DP）处理具体 gap。

### 2.4 SD Decomposition：k-mer Chaining + Set Cover

**Elementary SD**: SD 常由更古老的 SD 片段拷贝拼接而成。每个 SD 可分解为若干*elementary SDs*的拼接，
其中每个 elementary SD 在不同 SD 中有多个相似拷贝。

**分解算法**:

1. 将所有 SD 覆盖区域 `R` 取出，按重叠关系用 union-find 聚类（`cluster.codon`）。
2. 对每个 cluster，构建 `R` 上所有 k-mer 的索引 `I_k`（不使用 winnowing，默认 `k=10`）。
3. 用与 search 类似的 plane-sweep，在 `R` 上扫描相同 k-mer 的多个位置，通过距离阈值 `d_g = 50`
   将位置链式合并，得到 putative elementary SDs。
4. 当一个 elementary SD 不再能延伸时，输出长度超过 `μ`（默认 100bp）的拷贝集合。

**Core duplicon**: 定义为能覆盖所有 SD 的 elementary SD 最小集合。BISER 使用贪心 set-cover 近似算法
（`cover.py` 中的 `greedy_set_cover`）识别 core duplicons，并在 `.elem` 文件中标记为 `CORE`。

### 2.5 多基因组扩展

对每个基因组独立执行 detection 与 alignment 后，BISER 将每个基因组的 SD 映射到其他基因组
（cross_search）。这避免了把两个基因组之间的保守区域误判为 SD。cross_search 本质上仍是 plane-sweep：
以基因组 A 为参考建索引，用基因组 B 中已提取的 SD 区域作为 query 进行扫描。

## 3. BISER 代码实现

### 3.1 构建与架构

- **实现语言**: 核心算法用 [Codon](https://github.com/exaloop/codon/)（带 Seq plugin）编写，编译为
  `biser/exe/biser.exe`；流程调度与多进程用 Python (`biser/__main__.py`) 完成。
- **构建**: `setup.py` 中的 `CustomBuild` 调用
  `codon build -plugin seq biser/codon/__init__.codon -release -o biser/exe/biser.exe`，并拷贝
  Codon runtime 动态库。
- **依赖**: Python 侧需要 `tqdm`、`ncls`、`multiprocess`；运行时需要 `samtools faidx`。

### 3.2 源码模块与论文算法对应

- `biser/codon/__init__.codon`（入口）：子命令分发、参数解析、全局常量（`MAX_ERROR`、`KMER_SIZE`、
  `WINNOW_SIZE` 等）。
- `biser/codon/search.codon`（Putative SD detection）：k-mer 索引、winnowing、plane-sweep、putative
  SD 输出；含 `cross_search`。
- `biser/codon/hit.codon`（数据结构）：`Hit` / `Chromosome` / `Locus` 定义、CIGAR 操作、错误率计算、
  `merge`。
- `biser/codon/chain.codon`（Chaining）：`PrioritySearchTree` 实现、`chain()` 函数完成 `O(n log n)`
  锚点链构建。
- `biser/codon/align.codon`（Alignment refinement）：`generate_anchors()` 10-mer 锚点生成、
  `refine()` 基于 DP 的边界精修。
- `biser/codon/cluster.codon`（SD clustering）：重叠 SD 的 union-find 风格聚类、提取每个 cluster 的
  FASTA。
- `biser/codon/decompose.codon`（SD decomposition）：k-mer chaining 分解 elementary SDs。
- `biser/codon/mask.codon`（预处理/后处理）：hard-mask 生成、hard-masked 与原基因组坐标互转。
- `biser/cover.py`（Core duplicon）：用 `ncls` 做区间重叠查询，再用贪心 set cover 找 core duplicons。
- `biser/__main__.py`（Pipeline）：多进程任务拆分、临时目录管理、各阶段串接。

### 3.3 关键实现细节

#### 3.3.1 Plane-sweep in `search.codon`

`build_index()` 同时承担两种角色：

- `build_index=True`: 扫描参考基因组，建立 k-mer → 位置列表的索引。
- `find_sds=True`: 用已建索引对 query 序列做 plane-sweep，更新 `ListNode` 链表并输出 hits。

**2-bit 编码**: 代码中对碱基做 `(int(si) & 3)` 得到 2-bit 值，即 `A=0, C=1, G=2, T=3`。
滚动哈希 `h = ((h << 2) | base) & ((1 << 28) - 1)`（`KMER_SIZE=14`）生成 28-bit k-mer hash。

**Winnowing 实现**: `build_index()` 中用单调栈/队列维护 `(hash, pos)`：
- 新 k-mer 入队前，从队尾弹出 hash 不小于新 hash 的元素，保证队首为窗口最小值。
- 从队首弹出位置超出当前窗口 (`pos < i - KMER_SIZE + 1 - WINNOW_SIZE`) 的元素。
- 窗口填满 (`i - KMER_SIZE + 1 >= WINNOW_SIZE`) 后，队首 hash 即为一个 fingerprint；
  只有 fingerprint 变化时才进入 plane-sweep/索引插入，避免同一 hash 连续处理。

**Plane-sweep 链表 `update_list()`**: 维护按 `(chr, first)` 排序的 `ListNode` 链表，
每个节点保存：
- query 区间 `(first, last)` 与 reference 对应区间 `(ref, ref_last)`；
- 已扫描步数 `age`、命中次数 `count`；
- 是否曾经满足阈值 `potentional`。

对当前 query 位置 `current.loc` 对应的 reference 位置列表 `loci`（已按 `(chr, loc)` 排序），
`update_list()` 依次处理三种情况：
1. **延伸**: 若 `loci[lidx]` 与某 walker 同属一条 chromosome，且距离在 `MAX_DISTANCE` 内
   (`loci.loc - MAX_DISTANCE ≤ walker.last < loci.loc`)，则延伸 walker 的 `last`/`ref_last`
   并增加 `count`。
2. **插入**: 若 `loci[lidx]` 落在当前 walker 与下一个 walker 之间，则插入新节点，
   年龄初始为 0（本轮结束时再 `age += 1`）。
3. **老化**: 若当前 walker 未被延伸或插入，则 `age += 1`，并检查 `count ≥ ceil(age · τ)`。
   - 若满足且长度 `< MAX_SD_LEN`，标记 `potentional`。
   - 若不满足或长度超限，且此前为 `potentional`，并满足 `last - first > QUERY_THRESHOLD`
     与 `current.loc - walker.ref ≥ REF_THRESHOLD`，则调用 `save_sd()` 输出该 hit。

**Tau 计算**: `tau()` 中先算 `gap_error = min(1.0, (MAX_ERROR - MAX_EDIT_ERROR) / MAX_EDIT_ERROR * MAX_EDIT_ERROR)`，
再返回
`((1 - gap_error) / (1 + gap_error)) * (1 / (2 * exp(KMER_SIZE * MAX_EDIT_ERROR) - 1))`。
此即论文中的 Jaccard 下界。

**输出 `save_sd()`**: 输出前对边界做 `MAX_EXTEND` 填充（同 chromosome 时避免两个 mate 重叠），
并过滤掉完全自重叠（same species/name/strand 且坐标相同）的情况。最终生成 `Hit` 存入按
`(species1, chr1, species2, chr2, complement)` 分组的 `result` 字典。

#### 3.3.2 Chaining in `chain.codon`

`PrioritySearchTree` 是 chaining 的核心数据结构。它基于 y 坐标建静态二叉搜索树，每个叶子对应一个 anchor 的
y 位置 `(y + l - 1, anchor_idx)`。支持：

- `activate(x, score)`: 激活一个点并赋值。
- `deactivate(x)`: 将点置为 `-INF`。
- `rmq(lo, hi)`: 在区间 `[lo, hi]` 内找得分最高的激活点。

`chain()` 流程：

1. 对 reference 上的每个 anchor 左端点 `x` 和右端点 `x + l` 生成事件，按 x 排序扫描。
2. 在左端点事件时：
   - 先 deactivate 距离过远（`x - (anchor_j.x + anchor_j.l) > MAX_CHAIN_GAP`）的旧 anchor。
   - 在 PST 中查询 y 区间 `[ay - MAX_CHAIN_GAP, ay - 1]` 内得分最高的已激活 anchor `j`。
   - 若 `j` 存在，计算 gap cost `dx + dy`（`dx = ax - (jx + jl)`，`dy = ay - (jy + jl)`），
     更新 `dp[i] = max(MATCH_SCORE + (MATCH_SCORE//2)*(al-1) + dp[j] - gap, 自身得分)`，
     并设置 `prev[i] = j`。
3. 在右端点事件时，以 `dp[i] - gap_to_end` 激活当前 anchor，其中 `gap_to_end` 是从当前 anchor 到
   序列末端的最大可能距离惩罚。
4. 按 DP 得分排序并重构链，过滤掉 span 小于 `MIN_UPPERCASE_MATCH` 或 `MIN_READ_SIZE * (1 - MAX_ERROR)` 的链。

anchor 得分模型：第一个碱基得 `MATCH_SCORE`（默认 4），每多一个匹配碱基加 `MATCH_SCORE // 2`（即 2）。
gap 惩罚为 reference gap 与 query gap 之和（`dx + dy`），无 open 项。

#### 3.3.3 Alignment refinement in `align.codon`

`generate_anchors(h, KMER_SIZE=10)` 在两个 mate 序列间找 10-mer 精确匹配锚点：
- 先在 reference mate (`h.y`) 上建 hash 表 `ref_hashes: kmer -> [pos]`。
- 再扫描 query mate (`h.x`)，对每个 k-mer 查找匹配位置。
- 使用 `slide[d]` 数组（`d = len(x) + yi - xi`）避免同一对角线上被包含的短匹配：
  只保留每个对角线上能向右延伸最远的 anchor。
- 对保留的匹配向右延伸，直到遇到 `N` 或碱基不一致。
- 同染色体时跳过与对角线距离小于等于 `KMER_SIZE` 的 trivial 匹配。
- 只使用出现次数 `0 < count < 1000` 的 k-mer。

`refine(orig_h, hits, max_iter=500)` 对 anchors 做第二次 DP 精修：
- anchors 先按 `(x.chr, y.chr, x.end, x.start, y)` 排序。
- 每个 anchor 的自身得分 `score = MATCH * matches - MISMATCH * mismatches - GAP * indel_bp`。
- DP 仅向后查找最多 `max_iter=500` 个候选前驱；对每对 `(aj, ai)` 计算：
  - `mi = min(c_xs - p.x.end, c_ys - p.y.end)`
  - `ma = max(c_xs - p.x.end, c_ys - p.y.end)`
  - 若 `ma ≥ MAX_GAP` 或同染色体且中间无 gap，则跳过。
  - 更新 `dp[ai] = max(dp[ai], dp[aj] + score[ai] - MISMATCH * mi - GAPOPEN - GAP * (ma - mi))`。
- 按 `dp` 得分降序回溯，得到不相交的 anchor 链；合并被已有结果完全覆盖的链。
- 对最终链，用 `Hit(orig_h, [anchors], SIDE_ALIGN)` 拼接 anchors：
  - 相邻 anchors 之间的 gap 调用 `align_gap()`，小 gap（≤1000）用 `bio.seq.align()` 做精确比对；
    大 gap 用两端各比对 1000bp 并取分高者，中间用 `I`/`D` 表示。
  - 在链两端各取 `SIDE_ALIGN=500` 做 `ltrim()` / `rtrim()`：从端点向内侧扫描，找到累积比对得分
    最大的位置作为新边界。

`hit.codon` 中的 `Hit.align()` 会调用 `bio.seq.align()`，并将返回的 `M` 操作细分为 `=`（匹配）与
`X`（错配），CIGAR 操作仅使用 `=`, `X`, `I`, `D`。

#### 3.3.4 Decomposition in `decompose.codon`

分解阶段读取 `cluster` 输出的每个 cluster FASTA，对每个序列扫描 `k=10` 的 k-mer。

**索引构建**: 对 cluster 中所有序列（来自 `cluster.codon` 输出的 `.fa`）建立完整 k-mer 索引，
`kmer.as_int()` 作为 key，记录每个 k-mer 出现的 `(chr_id, loc)` 列表。同样用频率阈值过滤最高的
0.1% k-mer。

**`update_list()` 与 search 的差异**: search 阶段每个节点只跟踪一段 putative SD 的两个 mate；
decompose 阶段需要跟踪同一 elementary SD 在多个序列上的多个拷贝，因此每个 `ListNode` 额外维护：
- `mappings: Dict[int, int]`：记录该节点对应 elementary SD 在每个 chromosome 上的当前最右边界；
- `gap`：自上次命中以来的未命中步数；
- `score`：命中计数。

处理当前 k-mer 的位置列表 `index`（按 `(chr, loc)` 降序排列）时：
1. 对尚未进入链表的新位置，在链表头部插入新节点，并继承当前 `mappings`。
2. 对可延伸的节点，更新 `end`、`score`、`gap=0`、`count += 1`。
3. 对老化节点（`gap >= diff=50`），若此前为 `potentional` 且长度超过 `MIN_MATE_LEN`，
   调用 `process()` 输出：
   - 若 `mappings` 为空，输出 `(chr, begin, end, score)`；
   - 否则对每个 chromosome 输出从 `begin` 到 `mappings[chr]` 的区间，并更新 `begin` 为 `mappings[chr]+1`。
4. 清空已输出区域对应的 `visited` 标记，避免同一碱基被多次分解。

**合并 `merge()`**: 对相邻的 elementary SD 集合，若它们在所有 chromosome 上连续且间隔不超过 500bp，
则合并为一个集合。

**输出**: 每个 elementary SD 集合输出为 BED 行，格式为
`species\tchrom\tbegin\tend\tset_id\tlength\tscore\tstrand`；其中 strand 来自序列名中的 `+`/`-`。

## 4. 相关研究：T2T-CHM13 SD 注释来源与下游分析

本节记录 T2T-CHM13 基因组中 SD 注释的产生方式，以及基于这些注释开展的 SNV / IGC 研究。

### 4.1 T2T-CHM13 SD 注释的算法来源（Vollger et al. 2022）

> Vollger M. R. et al. *Segmental duplications and their variation in a complete human genome*
> . Science. 2022;376:eabj6965.
> [https://doi.org/10.1126/science.abj6965](https://doi.org/10.1126/science.abj6965)

T2T-CHM13 v1.1 的 SD 注释直接继承自 Vollger et al. 2022 对完整 T2T 基因组所做的 SD 注释。
该注释的生成流程如下：

1. **输入**: T2T-CHM13 v1.0 组装，并拼接上 GRCh38 的 chrY（用于包含 Y 染色体 SD 信息）。
2. **屏蔽重复序列**: 使用**Tandem Repeats Finder (TRF)**与**RepeatMasker**对组装进行 masking，
   从而只保留“常见重复序列之外”的区域进行同源搜索。这一步避免了卫星 DNA、简单重复等序列被误判为 SD。
3. **SD 检测**: 使用**SEDEF v1.1-31g68de243**（BISER 的前身）对屏蔽后的组装进行全基因组自比对。
    - SEDEF 的核心算法与 BISER 类似：基于**Jaccard similarity**做快速过滤，再用**local chaining**
      精修比对边界。
    - 与早期 WGAC 只能处理约 10% pairwise error 不同，SEDEF 可捕获**高达 25% pairwise error**的 SD，
      从而能检测更古老的重复事件。
4. **输出**: 将 SEDEF 检测到的、经 TRF/RepeatMasker 屏蔽后剩余的同源片段作为
   **T2T-CHM13 v1.0 的最终 SD 注释**。
5. **数据发布**: 完整注释流程以 Snakemake workflow 形式发布在 Zenodo：
   `10.5281/zenodo.5498988/workflows/sedef`

**结果规模**: 在 T2T-CHM13 中鉴定出约**208 Mbp 非冗余 SD 序列**（包含来自 GRCh38 chrY 的 15.6 Mbp），
使人类基因组中 SD 占比从 GRCh38 的约 5.4% 提升到约**7%**（218 Mbp / 3.1 Gbp）。其中约三分之二的
acrocentric 短臂序列由 SD 构成。

### 4.2 基于单倍型组装的 SD/IGC 分析（Vollger et al. 2023）

> Vollger M. R. et al. *Increased mutation and gene conversion within human segmental duplications*
> . Nature. 2023;617:335–344.
> [https://doi.org/10.1038/s41586-023-05895-y](https://doi.org/10.1038/s41586-023-05895-y)

这篇 Nature 文献的核心目标不是“发现新的 SD”，而是在已鉴定的 SD 区域内系统比较 SNV 模式并检测
**interlocus gene conversion (IGC)**。其方法亮点在于利用高质量单倍型组装建立 1:1 orthologous
alignment，从而绕过短读比对在 SD 区域的定位难题。

#### 4.2.1 数据与 SD 区域定义

- **样本**: 102 个人类 haplotype-resolved 基因组，主要来自 HPRC（94 个）和其他已发表组装（8 个）。
- **参考**: T2T-CHM13 v1.1（含完整 SD 区域）。
- **SD 定义**: 直接采用 T2T-CHM13 v1.1 的 SD 注释：所有非等位基因间或染色体内成对比对，长度
  `>1 kbp`、序列同一性 `>90%`，且不完全由常见重复或卫星序列组成。
- **Unique 区域定义**: T2T-CHM13 中不属于 SD、ancient SD（<90% 同一性）、着丝粒或卫星阵列的区域。
- **过滤**: 排除 Tandem Repeats Finder 鉴定的串联重复区域，并用 RepeatMasker 做额外重复类别注释。

#### 4.2.2 同线性 1:1 alignment 策略

为避免拷贝数变异区域带来的多对多比对歧义，研究只保留 1:1 同线性块：

1. 用**minimap2 v2.24**将每个 query 单倍型比对到 T2T-CHM13 v1.1：
  ```bash
  minimap2 -a -x asm20 --secondary=no -s 25000 -K 8G ref.fa query.fa
  ```
2. 用**rustybam v0.1.29**处理 PAF：
    - `trim-paf`: 去除 query 上的冗余比对。
    - `break-paf`: 在 `>10 kbp` 的结构变异处断开比对。
3. 保留连续比对长度 `>1 Mbp` 的区块作为**syntenic 1:1 alignment**。

该策略使研究能聚焦于拷贝数基本不变的 SD 区域（约 120 Mbp）及其侧翼 unique 序列，排除大尺度 SV 干扰。

#### 4.2.3 IGC 检测：双重比对法

IGC 会导致某个 haplotype 上的 SD 序列与其参考位置上的 ortholog 不一致，而与另一个 paralogue 更相似。
检测思路是比较“基于侧翼信息的同线性比对”和“不依赖侧翼的独立重比对”。

**算法流程**:

1. **同线性比对（Alignment 1）**: 用 minimap2 得到 query 单倍型与参考的 1:1 syntenic alignment。

2. **窗口化独立重比对（Alignment 2）**: 将 syntenic block 切成**1-kbp 窗口**，以**100-bp 步长**
   滑动，独立地重新比对回参考，找到每个窗口的 single best alignment position。

3. **候选 IGC 窗口**: 若某窗口在 Alignment 2 中的最佳位置与 Alignment 1 中的同线性位置**不重叠**，
   则标记为候选 IGC 窗口。

4. **合并**: 当相邻候选窗口在 donor 和 acceptor 序列上都连续重叠时，合并为更大的 IGC interval。

5. **SNV 支持计数**: 利用 CIGAR 字符串分别计算 donor 位点（Alignment 2）和 acceptor 位点（Alignment
   1）的匹配/错配碱基数，统计支持该转换事件的 SNV 数量。

6. **置信度评估**: 用累积二项分布计算每个候选 IGC 的*P*值：

  ```
  P(X ≤ k) = B(k, n, p)
  ```
其中*n*是两个 paralogue 之间的 informative site 数，*k*是支持 acceptor（未转换）序列的 site 数，
   *p = 0.5*（假设支持碱基可来自 donor 或 acceptor）。

**结果**: 在 102 个 haplotype 中平均每个个体检测到约 1,193 个 putative IGC 事件，高置信度（*P* <
0.05）callset 中约 4.3 Mbp 的 SD 序列受影响。

#### 4.2.4 质量控制策略

论文用多层独立数据验证组装的 SD 区域可靠性：

- **拷贝数**
    - **方法**: FastCN + 31-mer Illumina WGS
    - **关键结论**: 893 个测试中 756 个完全一致，差异多为 1 个拷贝
- **碱基准确性**
    - **方法**: Merqury + Illumina
    - **关键结论**: SD 区域平均 QV = 53（约 <1 SNV / 200 kbp）
- **组装连续性**
    - **方法**: GAVISUNK + ONT
    - **关键结论**: SD 区域错误率 0.11%，与 unique 区域 0.14% 相当
- **相位错误**
    - **方法**: Verkko 重新组装 CHM1
    - **关键结论**: 结果与 HPRC trio hifiasm 一致，相位误差可忽略

#### 4.2.5 与 BISER 的对比

- **核心目标**
    - **BISER**: *de novo* SD 检测与分解
    - **Vollger et al. 2023**: 已鉴定 SD 中的 SNV/IGC 分析
- **输入**
    - **BISER**: 单个或多个基因组 FASTA
    - **Vollger et al. 2023**: 102 个人类单倍型组装 + T2T-CHM13 注释
- **SD 定义**
    - **BISER**: 自定义错误模型（可低至 75% 同一性）
    - **Vollger et al. 2023**: 基于 T2T-CHM13 注释（> 1 kbp, > 90% 同一性）
- **关键算法**
    - **BISER**: ordered Jaccard + plane sweep + chaining
    - **Vollger et al. 2023**: minimap2 1:1 alignment + 窗口化重比对 + 二项检验
- **输出**
    - **BISER**: SD pairs / elementary SDs / core duplicons
    - **Vollger et al. 2023**: SNV 密度、IGC 事件、hotspots、突变谱

## 5. 对 pgr 的启示

1. **低同源性重复检测的过滤策略**: BISER 的 ordered Jaccard + plane-sweep 提供了一种不依赖 MinHash
   的线性过滤思路。若 pgr 未来需要检测古老重复或低同源性 segdups，可考虑类似 colinear k-mer
   matching + 扫描线的设计。
2. **Winnowing 作为采样手段**: BISER 用 winnowing 将 k-mer 采样率降到约
   `2/(w+1)`，同时保证同一 windows 内必然命中。这与 pgr sketching 的有损采样不同，
   但为“在保敏感度的前提下降低索引规模”提供了可借鉴的参数化方法。
3. **Priority Search Tree 用于 chaining**: `chain.codon` 中的 PST 实现是一个清晰的 `O(n log n)`
   chaining 模板，可迁移到 pgr 的 PAF/anchor chaining 场景中。
4. **多阶段 pipeline 的薄壳设计**: BISER 将复杂算法放在 Codon 核心（`biser/codon/`），Python
   仅负责任务分发与文件 I/O，这与 pgr 的 `libs/` + `cmd_pgr/` 分层理念一致。
5. **长读组装 + 1:1 syntenic alignment 用于 SD 下游分析**: Vollger et al. 2023 的方法说明，
   对于高同一性 SD，可以通过 minimap2 + rustybam 的 `trim-paf` / `break-paf` 获得可靠的 1:
   1 orthologous alignment。若 pgr 未来需要比较多个组装在重复区域的变异，可借鉴其 `>1 Mbp`
   连续同线性块与 `>10 kbp` SV 断开的策略。
6. **窗口化重比对检测基因转换**: 该文献的 IGC 检测思路（1-kbp 窗口独立重比对 + donor/acceptor
   CIGAR 比较 + 二项检验）为 pgr 提供了一个可复用的“寻找 ectopic 最佳匹配”模板，
   可用于检测组装间的非等位基因转换或重复片段的迁移事件。

## 6. 迁移到 PGR 的可执行方案

> 本节基于对 BISER 源码的逐行阅读与 PGR 现有代码的深入对比，给出将 BISER 体系迁移到 PGR 的
> 具体算法映射、可复用组件清单、缺失模块清单以及分阶段实施计划。所有可复用组件均与代码库中的
> 实际 API 对应，目标是让迁移工作可以直接按图索骥。

### 6.1 迁移目标与边界

- **目标**: 在 PGR 中新增 `pgr sd` 命令族，实现与 BISER 等价的功能：
  putative SD 检测、局部比对精修、跨基因组映射、SD 聚类、elementary SD 分解、
  core duplicon 识别、hard-mask 坐标与原基因组坐标互转。
- **边界**: 本次迁移聚焦算法实现与命令接口；多进程调度、临时目录管理、 resume 等工程特性
  可在核心算法稳定后按 PGR 已有模式（如 `pgr pl` pipeline）补充。
- **原则**: 复杂算法放 `src/libs/`，`cmd_pgr/` 仅做参数解析、I/O 转换与调用。单命令专用的
  复杂逻辑也放 `libs/`。

### 6.2 坐标系统约定

PGR 内部不同模块混用 0-based half-open 与 1-based inclusive 两种约定。迁移前必须逐
个核对调用点，否则 `mask`、`translate`、`loc` 三个环节极易出现 `off-by-one` 或区间
开闭错误。

#### 6.2.1 0-based half-open `[start, end)`

- **BISER 内部**: putative SD 区间、anchor 坐标、elementary SD 边界、hard-masked 坐标
  均使用 0-based half-open。
- **`src/libs/chain/record.rs`**: `ChainHeader` 与 `Block` 的 `t_start/t_end/q_start/q_end`
  均为 0-based half-open；`ChainData` 的 `size/dt/dq` 也是相对增量。
- **`src/libs/chain/connect.rs`**: `ChainableBlock` 的 `t_start/t_end/q_start/q_end` 为
  0-based half-open；`chain_blocks` 内部 gap 计算 `dt = target.t_start - cand.t_end`、
  `dq = target.q_start - cand.q_end`，负值表示重叠。
- **`src/libs/ds/kdtree.rs`**: `KdTreeItem` 要求 `x_start/y_start` 为 0-based inclusive，
  `x_end/y_end` 为 0-based exclusive。自定义 anchor 类型实现该 trait 时务必注意这一
  开闭差异。
- **`src/libs/paf/cigar.rs`**: `slice_cigar_by_target(cigar, target_start, ts, te)` 的
  `ts/te` 为 0-based half-open；`project` 内部投影也按 half-open 处理。
- **`src/libs/paf/record.rs`**: `PafRecord` 的 `query_start/target_start` 为 0-based
  inclusive，`query_end/target_end` 为 0-based exclusive，符合 PAF 规范。
- **`src/libs/paf/fasta.rs::FastaStore::fetch_range(name, start, end)`**: `start/end` 为
  0-based half-open 的 `i32`，函数内部转成 `noodles` 的 1-based inclusive position。
- **`src/libs/io.rs::SequenceReader::read_sequence(name, start, end)`**: `start/end` 为
  `Option<usize>`，0-based half-open；`None` 表示从头/到尾。
- **`src/libs/fmt/twobit.rs::TwoBitFile::read_sequence(name, start, end, no_mask)`**: 参数
  为 0-based half-open；`no_mask=true` 时返回 uppercase，否则保留 soft-mask。
- **`src/libs/ds/bitmap.rs::BitMap`**: `set_range(start, len)` 与 `is_fully_set(start, len)`
  都按 0-based half-open `[start, start+len)` 解释。
- **`src/libs/ds/dupe_tree.rs::DupeTree`**: `add(start, end)`、`subtract(start, end)`、
  `count_over(start, end, threshold)` 均为 0-based half-open。
- **`src/libs/alignment/coords.rs`**: `reverse_range_pair(start, end, size)` 与
  `reverse_range_1based_pair(start, end, size)` 分别处理 0-based half-open 与
  1-based inclusive 链向反转。

#### 6.2.2 1-based inclusive `[start, end]`

- **`src/libs/fmt/fa.rs::mask_sequence(seq, spans, hard)`**: `spans` 是
  `intspan::IntSpan`，1-based inclusive；函数内部通过 `offset = lower - 1` 转成切片
  索引。该函数**保留序列长度**，与 BISER 的 hard-mask（删除 lowercase）不同。
  现有命令 `pgr fa mask` 也是基于 runlist 做长度保留的 hard/soft mask，不能替代
  BISER 的 `mask`。
- **`src/libs/fmt/fa.rs::find_masked_regions(seq, gap_only)`**: 返回 0-based inclusive
  的 `(begin, end)` 对，与 `mask_sequence` 的输入约定不同，不要直接混用。现有命令
  `pgr fa masked` 将该输出转换为 1-based inclusive 显示。
- **`src/libs/fmt/fa.rs::windows`**: 内部切片是 0-based half-open，但输出名称格式为
  `name:start-end`，其中 `start = 原 start + 1`（1-based inclusive），`end` 保持为
  切片末尾位置（同样按 1-based inclusive 显示）。
- **`src/libs/loc.rs`**: `intspan::Range`（支持 `chr:start-end` 与 `chr(-):start-end`）
  为 1-based inclusive；`slice_record` 与 `fetch_range_seq` 接收该 `Range`。
  当 `start == 0` 时，`fetch_range_seq` 返回整条序列。`slice_record` 在负链时会对切片做
  reverse complement。
- **`src/libs/alignment/coords.rs`**: `chr_to_align` / `align_to_chr` 的输入位置是
  1-based inclusive，且配合 `IntSpan`（由 `seq_intspan` 从对齐序列的 gap 列生成）使用。
  这两个函数只适用于**带 gap 的对齐坐标**与基因组坐标之间的转换，不适用于 BISER
  hard-mask 后产生的简单偏移映射。

#### 6.2.3 SD 模块内部建议

- 在 `src/libs/sd/` 内部统一使用 **0-based half-open**，仅在以下边界做显式转换：
  - 调用 `loc::fetch_range_seq` / `slice_record` 时，把 0-based half-open 区间转为
    `intspan::Range` 的 1-based inclusive；
  - 调用 `fa::mask_sequence` 时，把 0-based half-open 区间转为 `IntSpan` 的
    1-based inclusive；
  - 输出 BED/文件名时，若需与 BISER 保持一致，再决定使用 0-based 还是 1-based。
- hard-masked 坐标 ↔ original 坐标的映射表建议保存为 0-based half-open：
  `Vec<(orig_start, orig_end, masked_start, masked_end)>`。

### 6.3 BISER 算法阶段与 PGR 组件映射

#### 6.3.1 Hard-masking（`mask.codon`）

- **BISER 实现**
  - 文件: `biser/codon/mask.codon:7-15`
  - 行为: 读取 FASTA，只保留 uppercase A/C/G/T，其余字符（包括 lowercase a/c/g/t、N、IUPAC ambiguity、gap 等）全部删除，按 `width=80` 输出为新的 hard-masked FASTA。
    因此 hard-masked 序列长度会变短，需要额外记录 uppercase run 在原基因组中的坐标，供后续 translate 使用。
- **PGR 可复用组件**
  - `src/libs/fmt/fa.rs`:
    - `reader(infile) -> Result<fasta::io::Reader<...>>`: 顺序读取 FASTA，支持 stdin 与 gzip。
    - `writer_with_wrap(outfile, 80) -> Result<fasta::io::Writer<...>>`: 按 80 bp 换行输出（与 BISER 的 `width=80` 一致）。
    - `new_record(name, seq)` / `new_record_preserving_desc(name, source, seq)`: 构造 FASTA record。
    - `find_fasta_files(path) -> Vec<PathBuf>`: 递归收集 `.fa` 与 `.fa.gz` 文件，输入为文件时返回单元素 vec。
    - `build_gzi_index(path)`: 为 BGZF FASTA 构建 `.gzi` 索引，`FastaStore::new` 需要该索引做随机访问。
    - `find_masked_regions(seq, gap_only) -> Vec<(usize, usize)>`: 返回 0-based inclusive 的 masked 区间。`gap_only=false` 时返回 lowercase 或 `nt::is_n` 为 true 的字符（即 N/n 与 IUPAC ambiguity）所在区间；`gap_only=true` 时只返回 N/n 与 IUPAC ambiguity 所在区间。**注意** lowercase 的 A/C/G/T 属于 `gap_only=false` 的返回范围，但它们不是 `nt::is_n`。
    - `mask_sequence(seq, spans, hard) -> Result<String>`: `seq` 为 `&str`，`spans` 为 1-based inclusive 的 `IntSpan`；函数将区间内字符替换为 `N`（hard）或小写（soft）。它**保留序列长度**，与 BISER 的 hard-mask（删除 lowercase）行为不同，因此不能直接用于生成 hard-masked FASTA。
  - `src/libs/nt.rs`:
    - `NT_VAL: &[usize; 256]`: 将 ASCII 字节映射到碱基编码。**A/a→0, C/c→1, G/g→2, T/t→3, U/u→3**（U/u 与 T/t 共用编码 3）；M/R/W/S/Y/K/V/H/D/B 及其小写，以及 N/n，映射到 4；其余字符（包括 gap `-`、`*` 等）映射到 255 (Invalid)。
      **关键注意**：lowercase a/c/g/t/u 在 `NT_VAL` 中同样映射到 0/1/2/3，因此在做 BISER 风格的 2-bit 滚动哈希前，必须先完成 hard-mask（删除 lowercase），不能直接用 `NT_VAL[b] & 3` 处理原始序列。
    - `is_n(b) -> bool`: 当 `NT_VAL[b] == 4` 时返回 true，即 N/n 与所有 IUPAC ambiguity codes。lowercase a/c/g/t/u 返回 false。
    - `is_lower(b) -> bool`: 判断字符是否为小写 ASCII。
    - `to_nt(nt) -> Nt`: 将字节映射到 `Nt` 枚举（A/C/G/T/N/Invalid）。
    - `count_n(seq) -> usize`: 统计 `is_n` 为 true 的字符数量。
    - `complement(seq) -> impl DoubleEndedIterator<Item = u8>`: 正序互补迭代器。
    - `rev_comp(seq) -> impl Iterator<Item = u8>`: 反向互补迭代器，在 cluster 阶段构造反向链序列时可直接使用。该迭代器保留原字符大小写，因此 lowercase 碱基在反向互补后仍为 lowercase。
  - `src/libs/fmt/twobit.rs::Blocks::from_dna`: 在打包 DNA 为 2bit 时，非 A/C/G/T 字符被记为 N-block，lowercase A/C/G/T 被记为 soft-mask block。这是 2bit 写入时的 mask 语义，与 BISER hard-mask（删除字符）不同，但可作为参考。
    - **重要**：2bit 内部位编码为 `T=00, C=01, A=10, G=11`，而 BISER 的 2-bit 滚动哈希使用 `A=0, C=1, G=2, T=3`。如果直接读取 2bit 的 packed bytes 做哈希，必须重新映射位；若通过 `TwoBitFile::read_sequence` 读取字符串后再用 `NT_VAL` 编码，则天然得到 BISER 编码。
  - `src/libs/fmt/twobit.rs::TwoBitFile`:
    - 实现 `SequenceReader` trait 时固定调用 `read_sequence(..., no_mask=false)`，
      即通过 trait 接口读取会保留 soft-mask（返回 lowercase）和 N-block。SD 流程若要从 2bit 获得 hard-masked 后的 uppercase 序列，
      应调用其 inherent 方法 `TwoBitFile::read_sequence(name, start, end, no_mask=true)`，而不是 trait 方法。
    - `TwoBitFile::read_sequence` 的 `start/end` 为 `Option<usize>` 的 0-based half-open 区间；`no_mask=true` 时 soft-mask block
      也会被转为 uppercase，`no_mask=false` 时保留 lowercase；N-block 始终返回 `N`。
- **需要新增的实现**
  - BISER 的 `mask` 是**删除 lowercase bases**并输出 hard-masked FASTA（序列长度变短），没有现成函数。
  - 实现方式：读取 record 时只保留 uppercase A/C/G/T（即 `NT_VAL[b] <= 3 && !b.is_ascii_lowercase()`）。小写 a/c/g/t、N、IUPAC ambiguity 以及 gap 字符均删除。
  - 同步记录每个 uppercase run 在原序列中的 `[orig_start, orig_end)` 边界以及对应 hard-masked 坐标 `[masked_start, masked_end)`，存入 `Vec<(orig_start, orig_end, masked_start, masked_end)>` 供 `translate` 使用。
  - 输出用 `fa::writer_with_wrap(outfile, 80)`。
  - 注意：`fa::reader` 会一次性将整条 record 读入内存（`noodles_fasta` 的 `Record.sequence()` 返回完整序列）。人类尺度染色体（~250 Mbp）尚可接受，但若要在更大基因组或内存受限场景下处理，建索引阶段可改用 `src/libs/fmt/twobit.rs::TwoBitFile` 顺序扫描，区间提取再用 `twobit` 或 `loc`。
  - 与现有 `pgr fa mask` / `pgr fa masked` 的关系：`pgr fa mask` 基于 runlist 做长度保留的 hard/soft mask；`pgr fa masked` 只负责找出 masked 区域。BISER 的 `mask` 是独立功能，建议放在 `pgr sd mask` 中实现。

#### 6.3.2 Putative SD detection（`search.codon`）

- **BISER 实现**
  - 文件: `biser/codon/search.codon:189-343`
  - 核心: 2-bit 滚动哈希 + winnowing + plane-sweep 链表 + tau 阈值 + 输出候选 hit。
- **PGR 可复用组件**
  - `src/libs/nt.rs`:
    - `NT_VAL: &[usize; 256]` 将 A/a→0, C/c→1, G/g→2, T/t→3, U/u→3（U/u 与 T/t 共用 3）。该表可直接用于 BISER 风格的 2-bit 滚动哈希，但**前提是序列已经过 hard-mask**：只保留 uppercase A/C/G/T，其余字符（lowercase a/c/g/t/u、N、IUPAC ambiguity、gap 等）均已删除。
    - 编码方式：对 hard-masked 后的字节 `b`，若 `NT_VAL[b] <= 3 && !b.is_ascii_lowercase()`，则 2-bit 值为 `NT_VAL[b] & 3`（即 A=0, C=1, G=2, T=3，与 BISER 一致）；否则应作为无效字符跳过。注意 2bit 文件内部位编码与 BISER 不同，不要直接对 2bit packed bytes 使用 `NT_VAL`。
  - `src/libs/fmt/fa.rs`:
    - `reader()`: 顺序读取 FASTA，适合建 k-mer 索引。注意它会将整条 record 载入内存。
  - `src/libs/ds/bitmap.rs`:
    - `BitMap::new(size)` + `set_range(start, len)` + `is_fully_set(start, len)`: 0-based 位图，可用于标记 plane-sweep 或 decomposition 中已访问/已输出的基因组位置，避免同一碱基被重复命中。
- **不可直接复用**
  - `src/libs/hash.rs`: 提供 canonical minimizer 采样（`seq_sketch`、`JumpingMinimizer`）与 Jaccard/Mash 距离计算。
    `seq_sketch` 返回 `Vec<MinimizerInfo>`，包含 hash、seq_id、pos、strand，形式上类似 BISER 的
    winnowing 输出，但本质不同：
    - `hash.rs` 使用 fxhash/rapidhash/murmurhash 等哈希函数，且支持 canonical k-mer；
    - BISER search/decompose 依赖 exact 2-bit k-mer + winnowing（非 canonical、非 hash-based）。
    因此 `hash.rs` 的 minimizer 流程不能直接复用，仅可作为 sketch 验证或后续扩展使用。
- **需要新增的实现**
  - `src/libs/sd/kmer_index.rs`: exact 2-bit k-mer 滚动哈希、winnowing 采样、
    `kmer -> Vec<(chr_id, pos)>` 索引、频率阈值过滤（0.1%）。
  - `src/libs/sd/plane_sweep.rs`: `ListNode` 链表、`update_list()` 三种分支逻辑、
    `save_sd()` 输出、tau 计算。
  - `src/libs/sd/hit.rs`: SD hit 数据结构（坐标、species、chromosome、strand、CIGAR、
    error rate），可参考 `src/libs/chain/record.rs` 的 `Chain` / `Block` 设计。

#### 6.3.3 Alignment refinement（`align.codon` + `hit.codon`）

- **BISER 实现**
  - 文件: `biser/codon/align.codon:5-112`、`biser/codon/hit.codon:325-348`
  - 核心: 10-mer anchor 生成、PST chaining、sparse DP refine、CIGAR 精修。
- **PGR 可复用组件**
  - `src/libs/ds/kdtree.rs`:
    - `KdTree::build(indices, items)` + `update_scores(leaf_idx, score, items)` +
      `best_predecessor(target_idx, current_score, items, cost_func, lower_bound_func)` 是底层 chaining 引擎。
    - `KdTreeItem` trait 要求 `x_start/y_start` 为 0-based inclusive，`x_end/y_end` 为 0-based exclusive；
      `score` 用于叶子初始得分。在 `src/libs/chain/connect.rs` 的实现中 `x` 对应 query、`y` 对应 target；
      BISER 的 anchor 可映射为 `x_start=q_start, x_end=q_end, y_start=t_start, y_end=t_end`。
    - **关键限制：`KdTree` 不支持 deactivate**。BISER 的 PST 在扫描锚点时需要按事件激活/ deactivate 锚点（当锚点与当前扫描位置距离超过 `MAX_CHAIN_GAP` 时置为 `-INF`）。`KdTree::update_scores` 只会把叶子和祖先节点的 `max_score` 向上提升，不会向下衰减；一旦某个内部节点的 `max_score` 被设为高分，后续即使该子树下的所有叶子都应 deactivate，该节点仍保留旧高分，导致剪枝边界和 `best_predecessor` 结果错误。因此**不能直接用 `KdTree` 实现 BISER 的 event-driven PST chaining**。
    - 可行的替代方案：
      1. **按扫描线顺序维护一个“当前窗口内活跃锚点”集合**，对该集合重建（或增量维护）一棵 KD-tree/Fenwick tree/线段树。由于 BISER 的 `MAX_CHAIN_GAP` 有限，窗口内锚点数量通常可控，但每次扫描线推进都重建 KD-tree 的复杂度是 `O(w log w)`，总复杂度会上升到 `O(n·w·log w)`，不适合大规模数据。
      2. **在 y 坐标离散化后使用线段树或树状数组（Fenwick tree）维护每个 y 位置的最大 DP 得分**，扫描线从左到右推进时，在 y 区间 `[ay - MAX_CHAIN_GAP, ay - 1]` 内查询最大值，并用单点更新写入当前锚点的 DP 值。这是实现 BISER PST 的最简洁路径，时间复杂度 `O(n log n)`，且天然支持 deactivate（当锚点滑出窗口时将其对应 y 位置重置为 `-INF`）。
      3. 如果坚持使用 `KdTree`，只能用于“所有锚点同时激活、无 deactivate”的 chaining 场景，此时 gap 惩罚仍需在 `cost_func` 中按 `dx + dy` 计算，并通过返回 `None` 过滤超出 `MAX_CHAIN_GAP` 的前驱；但这与 BISER 的 sweep + PST 逻辑不等价。
  - `src/libs/chain/connect.rs`:
    - `chain_blocks(blocks, gap_calc, score_ctx, ...) -> Result<Vec<Chain>>` 是已经实现的完整 chaining DP，
      包含去重、merge、trim、score recalc。但它的打分模型面向 UCSC `axtChain`（`GapCalc` 取 `max(dq, dt)`、
      有 overlap trim 等），与 BISER 的 PST chaining 不完全等价。
    - `ScoreContext { t_2bit, q_2bit, matrix }` 提供序列读取与替换矩阵，用于 overlap trim 和最终 score recalc。
      当 `score_ctx` 为 `Some` 时，链构建完成并去重/merge 后，会调用 `trim_overlaps`：对相邻 block 的 target 重叠区，
      用 `SubMatrix` 分别计算把重叠区全归左 block 或全归右 block 的得分，取得分更高的切分点，然后重新计算链总得分。
      负链 query 的序列会通过 `reverse_range_pair` 取反向互补后再参与评分。
    - **不建议用 `chain_blocks` 做 BISER 原型**：由于 gap 模型差异（`max(dq, dt)` vs `dx + dy`），其输出链与
      BISER 会有系统性偏差，无法验证 PST chaining 逻辑是否正确。原型阶段应使用 y 离散化 + 线段树的 PST 实现。
  - `src/libs/ds/gap_calc.rs`:
    - `GapCalc::medium()` / `GapCalc::loose()` / `GapCalc::affine(open, extend)`: 预计算 gap cost 表。
    - **重要差异**: `GapCalc::calc(dq, dt)` 在 `dq > 0 && dt > 0` 时使用 `max(dq, dt)` 查表，而 BISER chaining
      要求 `dx + dy`。因此 BISER chaining 不能通过 `GapCalc` 表达，需要在 segment-tree/Fenwick PST 的
      查询更新逻辑中直接按 `dx + dy` 计算。
    - BISER alignment refinement 的 sparse DP 使用 `GAPOPEN` / `GAP`，单轴 gap 可用 `GapCalc::affine(gap_open, gap_extend)`
      近似；但 BISER refine DP 对双 gap 的惩罚公式特殊（`MISMATCH * mi + GAPOPEN + GAP * (ma - mi)`），
      需要在新模块中重新实现，不能简单套用 `GapCalc`。
  - `src/libs/poa/align.rs` + `src/libs/poa/graph.rs` + `src/libs/poa/poa.rs`:
    - `src/libs/poa/mod.rs` 只导出 `AlignmentParams`、`AlignmentType`、`Poa`。`ScalarAlignmentEngine` 和 `PoaGraph`
      需分别通过 `poa::align::ScalarAlignmentEngine` 和 `poa::graph::PoaGraph` 使用。
    - `ScalarAlignmentEngine::new(AlignmentParams { match_score, mismatch_score, gap_open, gap_extend }, AlignmentType::Local)
      提供 Smith-Waterman 局部比对；也支持 `SemiGlobal` 和 `Global`。
    - `Alignment { score, path }`：`score` 为最佳路径总得分；`path` 为最佳路径上的步骤序列。
    - `AlignmentParams::default()` 的默认值为 `match_score=5, mismatch_score=-4, gap_open=-8, gap_extend=-6`。
      这些默认值来自 SPOA，适合做小片段 pairwise alignment；若要与 BISER 的 `MATCH_SCORE=4` 等参数对齐，
      应显式构造 `AlignmentParams`。
    - **三种模式的精确语义（对线性 POA graph 做 pairwise 比对时）**：
      - `Local`：query 与 graph 均允许自由起点/终点，返回得分最高的局部子比对。适合在较大 gap 区域内寻找最优局部对齐块。
      - `SemiGlobal`：query 必须完整对齐，graph 的起点/终点自由。适合把一段 query 锚定到 reference 的任意子区间。
      - `Global`：query 必须完整对齐，graph 起点固定（必须从第一个节点开始），但 graph 终点自由（可以在任意节点结束）。因此它**不是**传统 Needleman-Wunsch 的双端固定全局比对；若要做两条序列完全对齐，需要在得到 alignment 后手动检查 path 是否覆盖 graph 首尾，或改用 `SemiGlobal` 并在后续截断。
    - 做 pairwise alignment 时，把其中一条 mate 构建成线性 POA graph：
      ```rust
      use pgr::libs::poa::graph::PoaGraph;
      use pgr::libs::poa::align::{ScalarAlignmentEngine, AlignmentParams, AlignmentType};

      let mut graph = PoaGraph::new();
      let mut prev = None;
      for &base in ref_seq {
          let node = graph.add_node(base);
          if let Some(p) = prev { graph.add_edge(p, node, 1); }
          prev = Some(node);
      }
      let engine = ScalarAlignmentEngine::new(params, AlignmentType::Local);
      let aln = engine.align(qry_seq, &graph);
      ```
    - `Alignment.path: Vec<(Option<usize>, Option<NodeIndex>)>` 的精确语义（path 已经按正向顺序排列）：
      - `(Some(seq_idx), Some(node_idx))`：序列碱基 `seq_idx` 与 graph 节点 `node_idx` 匹配/错配；
      - `(Some(seq_idx), None)`：序列碱基 `seq_idx` 在 graph 中对应位置为插入（CIGAR `I`）；
      - `(None, Some(node_idx))`：graph 节点碱基在序列中对应位置为删除（CIGAR `D`）。
      按 path 顺序遍历，同时推进 `ref_seq`（graph 节点碱基）和 `qry_seq`（序列索引），即可输出
      `=`/`X`/`I`/`D` CIGAR。
    - 从 `Alignment.path` 生成 CIGAR 的模板（使用 `CigarOp::try_new`，`CigarOp::new` 为 `pub(crate)` 不可外部调用）：
      ```rust
      use pgr::libs::paf::cigar::{CigarOp, format_cigar};
      use petgraph::graph::NodeIndex;

      fn path_to_cigar(path: &[(Option<usize>, Option<NodeIndex>)], ref_seq: &[u8], qry_seq: &[u8], graph: &PoaGraph) -> anyhow::Result<String> {
          let mut ops: Vec<CigarOp> = Vec::new();
          let mut push = |op: char| {
              if let Some(last) = ops.last_mut() {
                  if last.op() == op {
                      *last = CigarOp::try_new(last.len() + 1, op).unwrap();
                      return;
                  }
              }
              ops.push(CigarOp::try_new(1, op).unwrap());
          };
          for &(q_idx, node) in path {
              match (q_idx, node) {
                  (Some(i), Some(n)) => {
                      let rb = graph.graph[n].base;
                      if qry_seq[i] == rb {
                          push('=');
                      } else {
                          push('X');
                      }
                  }
                  (Some(_), None) => push('I'),
                  (None, Some(_)) => push('D'),
                  (None, None) => {}
              }
          }
          Ok(format_cigar(&ops))
      }
      ```
    - 也可以先把 `ref_seq` 加入 `Poa`，再对 `qry_seq` 调用 `engine.align`：
      ```rust
      let mut poa = Poa::new(params, AlignmentType::Local);
      poa.add_sequence(ref_seq);
      let aln = engine.align(qry_seq, poa.graph());
      ```
      注意 `poa.add_sequence()` 会通过 `PoaGraph::add_alignment` 修改 graph，因此若只想做**不改变 graph 的 pairwise alignment**，应手动构建线性 `PoaGraph`。
    - `PoaGraph::add_alignment` 的具体行为（理解 consensus/MSA 的基础）：
      - 对 `Alignment.path` 中每个 `(seq_idx, node_idx)` 步骤，若 `seq_idx` 存在则消费 query 碱基；
        `(None, Some)` 或 `(None, None)` 表示 deletion，不消费 query 碱基。
      - 在消费 query 碱基前，先把 path 中未对齐的 query 前缀（`(Some, None)` 之前的独立碱基）逐个加入 graph，
        作为新节点并用 weight=1 的边连接。
      - `(Some(seq_idx), Some(node_idx))` 且 graph 节点碱基与 query 碱基一致时，该节点 `weight += 1`，
        predecessor 边 weight += 1。
      - `(Some(seq_idx), Some(node_idx))` 但碱基不一致时，先检查 `node_idx.aligned_nodes` 中是否已有相同碱基的节点；
        若有则复用该节点，否则**新建一个节点**保存 query 碱基，并将其加入 `node_idx.aligned_nodes` 与原有 clique 合并；
        随后该节点 weight += 1。
      - 连续的 query 碱基通过 `add_edge` 连接，边权重累加；`add_edge` 对重复边做 weight 累加而非新建多重边。
      - 因此 `PoaGraph` 中节点 weight 表示某碱基在已加入序列中出现的次数，边 weight 表示相邻关系出现次数；
        `generate_consensus` 与 `generate_msa` 基于此做多数表决/回溯。
    - **性能注意**: `ScalarAlignmentEngine` 是标量 O(nm) 实现，无 SIMD/banded。小 gap（≤1000 bp）可直接使用；
      人类全基因组尺度批量调用可能成为瓶颈，届时再评估 `parasail-rs` 或 banded 优化。
    - 若需要多序列 consensus/MSA（例如 cluster 内多个拷贝），可用 `Poa::new(params, align_type)` +
      `add_sequence()` + `consensus()` / `msa()`，这比直接调 `ScalarAlignmentEngine` 更方便。
  - `src/libs/paf/cigar.rs`:
    - `parse_cigar(s) -> Result<Vec<CigarOp>>` / `format_cigar(ops) -> String`: CIGAR 解析与格式化。
    - `extract_cigar(tags)`: 从 PAF tags 提取 CIGAR。
    - `reverse_cigar(ops)`: 反向并交换 I/D。
    - `cigar_from_alignment(ref, qry)`: 从两条**等长对齐序列**（含 gap 字符 `'-'`）生成 `=` / `X` / `I` / `D` CIGAR，
      比对大小写不敏感。不能直接从 `ScalarAlignmentEngine` 的 `Alignment.path` 调用，需要先展开成等长对齐字符串。
    - `cigar_stats(ops) -> CigarStats` / `gap_compressed_identity(ops) -> f64` / `block_identity(ops) -> f64`:
      统计与 identity 计算。BISER 的 error rate 是编辑错误率（`E/ℓ`），应基于 CIGAR 计算
      `(mismatches + ins_bp + del_bp) / (matches + mismatches + ins_bp + del_bp)`，等价于
      `1 - block_identity(ops)`。`gap_compressed_identity` 把每个 indel 只计为一个事件，会低估
      错误率，不适合直接作为 BISER error rate。
    - `slice_cigar_by_target(cigar, target_start, ts, te)`: 按 target 子区间切 CIGAR。
    - `CigarOp` 使用 bit-packed `u32`：高 3 bits 存 op code（`=`/`X`/`I`/`D`/`M`），低 29 bits 存长度，最大单 op 长度约 512 Mbp。
  - `src/libs/alignment/stat.rs`:
    - `pair_d(seq1, seq2) -> Result<f32>`: 计算两条**等长对齐序列**的 divergence（忽略 gap 列，IUPAC 按 N 处理）。
      适合从展开后的对齐字符串计算，而不是直接用于 CIGAR。
    - `alignment_stat(seqs) -> Result<(i32, i32, i32, i32, i32, f32)>`: 多序列对齐列统计（长度、可比列、差异列、gap 列、模糊列、平均差异）。
  - `src/libs/alignment/msa.rs`:
    - `align_seqs(seqs, "builtin")` 调用内置 POA 做多序列对齐，返回 MSA 字符串（含 `'-'`），**不是 pairwise alignment**。
      也支持 `"spoa"`、`"clustalw"`、`"muscle"`、`"mafft"` 等外部 aligner。
    - `align_seqs_quick(seqs, aligner, pad, fill)` 在已有粗对齐（所有序列长度相同）的基础上，仅对
      head/tail 和 gap 邻近区域调用外部 aligner 重新对齐，再拼回原位。适合先快速得到整体对齐框架，
      再局部精修的场景。
    - `get_consensus_poa_builtin(seqs, match_score, mismatch_score, gap_open, gap_extend, algo_code)` /
      `get_consensus_poa_external(seqs, ...)` 直接用 POA 或外部 spoa 生成 consensus 字符串。
    - 对于 pairwise 小 gap，优先直接用 `ScalarAlignmentEngine`；对于多拷贝 consensus 或 cluster MSA，
      可复用 `align_seqs(..., "builtin")`、`align_seqs_quick` 或 `poa::Poa`。
  - `src/libs/chain/sub_matrix.rs`:
    - `SubMatrix` 提供 256×256 字节替换矩阵（含大小写）与 gap open/extend。
    - `SubMatrix::default()` 是一个简化的 identity-like 矩阵：A/C/G/T 匹配得 100，错配 -100，N 相关 -100，
      `gap_open=400`, `gap_extend=30`。它**不是** lastz 默认矩阵，仅作为通用 fallback。
    - `SubMatrix::hoxd55()` 才是 lastz 默认的 HoxD55 矩阵（A-A=91, C-C=100, G-G=100, T-T=91，非对角线负值），
      同样 `gap_open=400`, `gap_extend=30`。
    - `SubMatrix::from_name(name)` 支持 `"hoxd55"` 预设，其他名字按 BLAST 格式从文件解析。
    - `SubMatrix::get_score(c1, c2)` 按字符 ASCII 值查表，大小写均可。
    - `chain_blocks` 的 `ScoreContext` 使用 `SubMatrix` 在 overlap trim 时重新计算匹配得分。注意 `ScalarAlignmentEngine`
      **不使用** `SubMatrix`，它只接受简单的 `match_score`/`mismatch_score`；若 BISER 精修需要 HoxD55 等复杂矩阵，
      需改用外部 aligner（如 lastz）或自行实现支持替换矩阵的 DP。
  - `src/libs/lastz.rs`:
    - 提供 lastz 的预设评分矩阵与参数（`PRESETS`、`find_preset`、`run_lastz`）。若 `ScalarAlignmentEngine`
      性能不足或需要更复杂的评分矩阵，可将小 gap 区域提取后调用 lastz 作为外部 fallback。
  - `src/libs/fas_multiz/banded_align.rs`:
    - `banded_align_refs(...)` 实现 banded DP + affine gap，但当前紧密绑定 `FasBlock` 输入，不能直接复用。
      若后续需要 banded pairwise align，可参考其索引函数 `idx(i, j)` 和 band 半径计算逻辑，提取为通用函数。
  - `src/libs/alignment/trim.rs`:
    - 提供 `trim_head_tail`、`trim_complex_indel`、`trim_outgroup` 等函数，用于 MSA 后处理。
    - BISER 的 `ltrim`/`rtrim` 是从链两端向内扫描、找累积比对得分最大的边界，与这些通用 trim 不同，需要自行实现。
- **需要新增的实现**
  - `src/libs/sd/anchor.rs`: 在 putative SD 的两个 mate 间生成 10-mer exact-match anchors，
    包含 `slide[d]` 去重、向右延伸、过滤高频 k-mer、过滤 trivial self-overlap。
  - `src/libs/sd/refine.rs`: 基于 y 坐标离散化 + segment tree / Fenwick tree 的 event-driven
    PST chaining + sparse DP 精修、大 gap 处理、两端 score-based `ltrim`/`rtrim`、生成最终 CIGAR。
    小 gap（≤1000 bp）调用 `ScalarAlignmentEngine` 并把 path 转为 CIGAR；大 gap 按 BISER 策略
    （两端各比对 1000 bp，中间用 `I`/`D`）处理。

#### 6.3.4 SD clustering（`cluster.codon`）

- **BISER 实现**
  - 文件: `biser/codon/cluster.codon:53-165`
  - 核心: 对 hit 的四个端点做区间 coloring，等价于 union-find 找重叠 hit 的连通分量，
    然后提取每个 cluster 的序列 FASTA。
- **PGR 可复用组件**
  - `src/libs/paf/graph/dsu.rs`:
    - 已实现 union-by-rank + path compression 的 `Dsu`，但它是 `pub(super)`，仅在
      `paf::graph` 模块内部可见，**不能直接作为公共 API 使用**。
    - 根据项目约束，纯数据结构应放在 `src/libs/ds/`。建议将 `Dsu` 迁移到 `src/libs/ds/dsu.rs` 并公开为 pub，
      原 `paf/graph/dsu.rs` 通过 `pub use` 保持 API 兼容。SD 聚类直接复用 `src/libs/ds/dsu.rs::Dsu`。
  - `src/libs/ds/dupe_tree.rs`:
    - `DupeTree::add(start, end)` + `build()` + `count_over(start, end, threshold)`: 0-based 区间深度树，
      可用于统计 hit 端点或 cluster 区域的覆盖深度，识别高重复 hotspot。
  - `src/libs/ds/bitmap.rs`:
    - `BitMap::set_range` / `is_fully_set`: 0-based 位图，可用于标记已被 cluster 覆盖的碱基，避免重复提取。
  - `src/libs/fmt/fa.rs`:
    - `reader()` / `new_record()` / `writer()` / `writer_with_wrap()`: 读取基因组并构造 cluster FASTA。
  - `src/libs/fmt/fas.rs`:
    - block FA 读写与 `FasBlock` 数据结构。若 SD cluster 阶段需要输出多序列比对块（类似 MAF/block FA），
      可参考该模块，但它紧密围绕 block FA 格式设计，不是通用 MSA 容器。
  - `src/libs/nt.rs`:
    - `rev_comp(seq)`: 生成反向互补序列，返回迭代器。
  - `src/libs/loc.rs`:
    - `create_loc(infile, locfile, is_bgzf) -> Result<()>`: 为 plain 或 BGZF FASTA 创建 `.loc` 索引。
      普通代码更常用 `open_indexed` 自动创建。
    - `open_indexed(infile, force_update) -> Result<(Input, IndexMap<name, (offset, size)>)>`: 打开带 `.loc` 索引的 FASTA
      （plain 或 BGZF 均可），不存在时自动创建索引。内部通过 `is_bgzf` 判断压缩类型。
    - `open_input(infile, is_bgzf) -> Result<Input>`: 打开 FASTA 为 `Input::File` 或 `Input::Bgzf`。
    - `fetch_record(reader, loc_of, name) -> Result<fasta::Record>`: 按名字读取完整 record。
    - `fetch_range_seq(reader, loc_of, rg) -> Result<String>`: 按 `intspan::Range`（1-based inclusive，支持 `chr(-):start-end` 链向）
      提取子序列。
    - `slice_record(record, rg) -> Result<fasta::record::Sequence>`: 从已加载 record 中按 1-based Range 切片，负链会返回 reverse complement。
    - `get_seq_loc(file, range) -> Result<String>`: 便捷函数，对无效 range 或找不到的 chromosome 返回空字符串而非报错；
      测试/脚本中可用，但生产代码建议用 `open_indexed` + `fetch_range_seq` 以明确错误处理。
  - `src/libs/io.rs`:
    - `SequenceReader` trait: `read_sequence(name, start, end) -> Result<String>`，定义 0-based half-open 的随机访问接口。
      `TwoBitFile` 实现该 trait，`chain_blocks` 的 `ScoreContext` 也依赖它，因此 SD 流程中可用统一接口切换 2bit / indexed FASTA。
  - `src/libs/fmt/twobit.rs`:
    - `TwoBitFile::read_sequence(name, start, end, no_mask)`: 0-based half-open 随机访问 2bit 序列，
      适合区间提取；`no_mask=true` 时返回 uppercase，否则保留 mask（N-blocks 变 N，soft-mask 变 lowercase）。
    - `TwoBitFile` 实现 `SequenceReader` trait，可直接传给 `ScoreContext` 等需要序列读取的接口。
  - `src/libs/paf/fasta.rs`:
    - `load_fasta_tsv(path) -> Result<IndexMap<name, path>>`: 读取 `name\tbgzf_fasta_path` 格式的 TSV，
      用于将 PAF/SD 中的序列名映射到 BGZF FASTA 文件路径。
    - `prepare_store(tsv_path, idx) -> Result<FastaStore>`: 加载 TSV、校验覆盖所有 `idx.names`、
      构造 `FastaStore` 的一站式函数；SD 流程若用 TSV 管理多基因组输入，可直接复用。
    - `load_all_seqs(tsv_path) -> Result<HashMap<name, seq>>`: 一次性加载 TSV 中所有序列到内存，
      适合 cluster/decompose 阶段加载单个 cluster 的小规模 FASTA。
    - `FastaStore::new(seq_to_file)` + `fetch_range(name, start, end)` + `fetch_full(name)`: 管理多个 **BGZF FASTA**
      文件，带 `.loc` 索引与 LRU 缓存，适合多基因组 cross_search/cross_align 时批量提取 mate 序列。
    - **限制**: `FastaStore::new` 内部使用 `noodles_bgzf::io::indexed_reader`，因此输入必须是 BGZF 压缩的 FASTA。
      普通 gzip 或未压缩 FASTA 应先用 `loc`（plain/BGZF 通用）或 `twobit` 处理。
    - **注意**: `FastaStore` **没有**实现 `SequenceReader` trait，它提供自己的 `fetch_range` / `fetch_full` API。
      若函数签名要求 `&mut dyn SequenceReader`（如 `chain_blocks` 的 `ScoreContext`），应使用 `TwoBitFile` 而非 `FastaStore`。
- **需要新增的实现**
  - `src/libs/sd/cluster.rs`: 将 hit 端点排序，用 `Dsu` 合并重叠端点，输出每个 cluster 的
    FASTA（序列名采用 `species#chrom+/-#start#end` 格式）。

#### 6.3.5 SD decomposition（`decompose.codon`）

- **BISER 实现**
  - 文件: `biser/codon/decompose.codon:52-277`
  - 核心: 对 cluster FASTA 建 10-mer 完整索引，再用 plane-sweep + mappings 输出
    elementary SD 集合。
- **PGR 可复用组件**
  - `src/libs/sd/kmer_index.rs`: 复用 exact 10-mer 索引（调整 `k` 参数即可）。
  - `src/libs/sd/plane_sweep.rs`: 复用链表扫描框架，但 decompose 需要支持多拷贝 mappings，
    因此 `PlaneSweepState` 需要泛化或单独实现一个 `MultiCopyPlaneSweep`。
  - `src/libs/ds/bitmap.rs`: 用于标记已输出的 `visited` 区域。
- **需要新增的实现**
  - `src/libs/sd/decompose.rs`: 在 cluster FASTA 上调用 10-mer 索引 + 多拷贝 plane-sweep + merge，
    输出 `.elem` 格式 BED。

#### 6.3.6 Core duplicon identification（`cover.py`）

- **BISER 实现**
  - 文件: `biser/cover.py:1-102`
  - 核心: 用 `ncls` 建立 elementary SD → 覆盖的 SD 列表映射，再用贪心 set cover 找出
    能覆盖所有 SD 的最小 elementary SD 集合，标记为 `CORE`。
- **PGR 可复用组件**
  - `coitrees`（已通过 `src/libs/paf/index/builder.rs` 使用）:
    - 直接建立 `Interval<ElemMetadata>` 区间树，查询重叠 interval，替代 `ncls`。
    - 比 `PafIndex` 更轻量；`PafIndex` 的 `query()` 与 `query_transitive_bfs()` 展示了
      “interval tree + CIGAR 投影”模式，可作为复杂场景参考。
  - `src/libs/paf/index/query.rs`:
    - `project(ts, te, metadata, cigar)` 实现了 target 子区间到 query 坐标的投影，
      输入 `ts/te` 为 0-based half-open，返回 query 区间也是 0-based half-open。
      在将 elementary SD 区间映射到 SD 覆盖时思路类似，但注意它处理的是带 CIGAR 的对齐投影，
      不是简单的坐标偏移。
- **需要新增的实现**
  - `src/libs/sd/set_cover.rs`: 用二叉堆实现 `greedy_set_cover`，输入为
    `elementary_id -> Vec<sd_id>`，输出被选中的 core elementary IDs。
  - `src/libs/sd/cover.rs`: 组装 interval overlap + set cover 流程，在 `.elem` 文件中追加
    `CORE` 标记。

#### 6.3.7 Coordinate translation（`mask.codon:32-154`）

- **BISER 实现**
  - 文件: `biser/codon/mask.codon:32-154`
  - 核心: 记录原始序列中大写区段的上边界 `uppers` 与原始坐标 `lowers`，在 hard-masked 坐标
    与原始坐标之间双向映射，并同步更新 CIGAR 中的 `M`/`I`/`D` 为 `M`/`S`/`N`。
- **PGR 可复用组件**
  - `src/libs/paf/cigar.rs`:
    - `parse_cigar()` / `format_cigar()` / `slice_cigar_by_target()`: 解析、切片、重构 CIGAR。
  - `src/libs/io.rs`:
    - `reader()` / `writer()`: 读取原始 FASTA 与输出转换后的结果。
  - `src/libs/fmt/fa.rs`:
    - `reader()` / `writer_with_wrap()`: 读取/输出带 coordinate translation 的 FASTA（若需要）。
- **不可直接复用**
  - `src/libs/alignment/coords.rs`: 这里的 `align_to_chr()` / `chr_to_align()` 是通过**带 gap 的对齐列**
    （用 `IntSpan` 记录非 gap 位置）在对齐坐标与基因组坐标之间转换，适用于 MSA/POA 输出。
    它**不适用于** BISER 的 hard-masked ↔ original 坐标映射，因为后者只是简单的“删除 lowercase 碱基”
    后产生的坐标偏移，不存在对齐 gap。
  - `src/libs/fmt/fa.rs::mask_sequence`: 该函数保留序列长度，只替换字符，不能用于 translate
    阶段的坐标映射。
- **需要新增的实现**
  - `src/libs/sd/translate.rs`: 在 `mask` 阶段记录 uppercase run 列表
    `[(orig_start, orig_end, masked_start, masked_end)]`（0-based half-open），通过二分查找实现
    hard-masked ↔ original 坐标双向映射，并同步改写 CIGAR。
  - 具体映射规则（需对照 BISER 源码精确实现）：
    - hard-masked 坐标下的 match/mismatch 区间，映射回 original 坐标时对应 uppercase 区段，CIGAR 保持 `=`/`X`/`M`。
    - 当 hard-masked 的 gap（`I`/`D`）跨越 original 中的 lowercase/masked 区域时，需要把该段替换为
      `S`（soft-mask）或 `N`（hard-mask），以反映原始基因组中这些碱基并非真正参与比对。
    - BISER 内部将 CIGAR 操作重新归类为 `M`/`S`/`N`；translate 阶段应产生与 BISER 输出语义一致的 CIGAR。

### 6.4 其他值得复用的工具

除了按算法阶段映射的组件外，以下通用工具在 SD 流程中也 likely 有用：

- `src/libs/par.rs`
  - 并行 pipeline 原语：`spawn_writer_and_pool` 创建 writer thread + rayon pool；
    `resolve_paths` / `load_entries` / `load_two_sets` / `par_run_pairs` 支持列表解析、批量加载与
    成对并行迭代。SD pipeline 的 search/align 阶段若需要“生产者-消费者”式并行输出，可直接复用。
- `src/libs/io.rs`
    - `reader(input)` / `writer(output)`: 通用缓冲读写，支持 `stdin`、普通文件、`.gz`。
    - `read_names(path)` / `read_sizes(path)`: 名单读取、大小文件读取。
    - `is_bgzf(path)`: 通过读取文件头判断是否为 BGZF 格式，`FastaStore` 与 `loc::open_indexed` 内部均用此决定打开方式。
    - `read_runlist(path)`: 安全读取 runlist JSON 并转为 `BTreeMap<String, IntSpan>`，避免原 `intspan` API 在错误输入上 panic。
    - `get_basename(path) -> Option<String>`: 提取文件基本名（去掉路径与扩展名），`lastz.rs` 与多个命令用它生成输出文件名。
    - `SequenceReader` trait: 统一 0-based half-open 随机访问接口，`TwoBitFile` 已实现。
- `src/libs/ds/bitmap.rs::BitMap`
  - 固定大小的 0-based 位图，支持 `set_range(start, len)` 和 `is_fully_set(start, len)`。
  - 用途: 标记 plane-sweep 或 decomposition 中已访问/已输出的基因组位置；避免重复命中。
- `src/libs/ds/dupe_tree.rs::DupeTree`
  - 一维 0-based 区间深度树，支持 `add/subtract` 后 `build()`，再 `count_over(start, end, threshold)`。
  - 用途: 统计 SD hit 端点或 cluster 区域在基因组上的覆盖深度，识别高重复 hotspot。
- `src/libs/ds/top_k_purity.rs::TopKPurity`
  - 跟踪离散类别计数，计算 top-K 类别占总观测的比例，并在比例过高时返回 penalty factor。
  - 用途: 在扫描 k-mer 或序列窗口时检测类别分布过于集中的低复杂度区域（不限于 AT 富集），作为额外过滤条件。
- `src/libs/fasta/stat.rs::count_bases`
  - 统计序列中 A/C/G/T/N 数量（IUPAC ambiguous codes 计为 N，其他非标准字符不计入长度）。
  - 用途: 快速评估 hard-masked 后有效碱基比例，或过滤 N 含量过高的 putative SD。
- `src/libs/fmt/twobit.rs::Block`
  - 2bit 内部使用的 0-based half-open mask block 类型。若 SD 流程需要记录 hard-masked 区间或
    N-block，可参考其区间重叠查询实现 `Blocks::overlaps`。
- `src/libs/paf/fasta.rs::FastaStore`
  - 多 **BGZF FASTA** 管理器，支持 `fetch_range(name, start, end)`（0-based half-open，参数为 `i32`）
    与 `fetch_full(name)`，带 `.loc` 索引与 LRU 缓存。
  - 用途: 多基因组 cross_search/cross_align 时批量、高效地提取 mate 序列。注意输入必须是 BGZF 压缩 FASTA。
- `src/libs/paf/persist.rs`
  - `PafIndex::save(path)` / `PafIndex::load(path)`: 将 interval tree + CIGAR 索引持久化为 `.paf.idx`。
  - 用途: 若 SD 流程需要将 k-mer/anchor 索引或 hit 索引缓存到磁盘，可参考其 bincode + version + magic 的序列化模式。
- `src/libs/lastz.rs`
  - 提供 UCSC lastz 预设（`PRESETS`、`find_preset`、`run_lastz`）与评分矩阵。可作为 `ScalarAlignmentEngine` 性能不足时的外部局部比对 fallback。
  - 内置矩阵：`MATRIX_DEFAULT`（HoxD55， Human/Mouse/Macaque/Cow）、`MATRIX_DISTANT`（Human/Zebrafish/Opossum）、
    `MATRIX_SIMILAR`（Human/Chimp）、`MATRIX_SIMILAR2`（Human/Primate，更敏感）。
  - 内置 preset `set01`–`set07` 的参数（如 `O=400 E=30 K=3000 L=2200` 等）直接来自 UCSC pipeline；
    `run_lastz` 会处理 target/query 笛卡尔积、文件名去重与 `--self` 模式。
- `src/libs/chain/stitch.rs`
  - 按 chain ID 合并 fragments（`chainStitchId` 语义），要求 fragments 之间不重叠。SD 流程若输出带 ID 的 chain fragments，可参考；但它不是通用的“合并相邻 chain”工具。
- `src/libs/alignment/trim.rs`
    - `trim_pure_dash(seqs)`：删除所有序列在该列均为 gap（`-`）的列。
    - `trim_outgroup(seqs)`：删除 outgroup-only 的插入列（当前实现取 gap 列的并集与交集关系判断）。
    - `trim_head_tail(seqs)` / `trim_complex_indel(seqs)`：其他 MSA 后处理辅助函数。
    - SD refine 阶段若需要对齐后修剪，可参考这些函数，但 BISER 的 score-based `ltrim`/`rtrim` 仍需自行实现。
- `src/libs/fas_multiz/banded_align.rs`
  - 实现 banded DP + affine gap，但当前紧密绑定 `FasBlock` 输入，不能直接复用。若后续需要 banded
    pairwise align，可参考其索引函数 `idx(i, j)` 和 band 半径计算逻辑，提取为通用函数。

### 6.5 建议的模块与命令结构

#### 6.5.1 `src/libs/sd/` 目录（新增）

- `kmer_index.rs`: exact 2-bit k-mer 索引 + winnowing + 频率过滤。
- `plane_sweep.rs`: plane-sweep 链表与 hit 输出。
- `hit.rs`: SD hit 数据结构（坐标、species、strand、CIGAR、error rate）。
- `anchor.rs`: 10-mer exact-match anchor 生成。
- `refine.rs`: 基于 y 坐标离散化 + segment tree / Fenwick tree 的 event-driven PST
  chaining + sparse DP 的比对精修；包含 `path_to_cigar` 辅助函数。
- `cluster.rs`: 重叠 hit 聚类并输出 cluster FASTA。
- `decompose.rs`: elementary SD 分解。
- `set_cover.rs`: 贪心 set cover。
- `cover.rs`: core duplicon 标记流程。
- `translate.rs`: hard-masked 与原基因组坐标互转。注意 `src/libs/translate.rs` 已存在（蛋白质翻译），
  新增 SD 的 `translate.rs` 位于 `src/libs/sd/translate.rs`，不会冲突。

#### 6.5.2 `src/cmd_pgr/sd/` 目录（新增）

- `mod.rs`: 子命令注册与分发。
- `mask.rs`: `pgr sd mask <genome.fa> -o <masked.fa>`。
- `search.rs`: `pgr sd search <genome.fa> -o <hits.bed>`。
- `align.rs`: `pgr sd align <genome.fa> <hits.bed> -o <hits.align.bed>`。
- `cluster.rs`: `pgr sd cluster <genomes...> <hits.align.bed> -o <clusters.dir>`。
- `decompose.rs`: `pgr sd decompose <cluster.fa> -o <cluster.elem.bed>`。
- `cover.rs`: `pgr sd cover <hits.align.bed> <elems.txt> -o <elems.covered.txt>`。
- `translate.rs`: `pgr sd translate <hits.align.bed> <genomes...> -o <out.bed>`。
- `run.rs`: `pgr sd run <genomes...> -o <out.bed>`，按 BISER 顺序串接上述步骤。

### 6.6 分阶段实施计划

#### 第一阶段：索引与 plane-sweep（验证：human chr21 自比对命中数与 BISER 一致）

1. 实现 `libs/sd/kmer_index.rs`：2-bit k-mer + winnowing + 频率过滤。
2. 实现 `libs/sd/plane_sweep.rs`：ListNode、update_list、tau、save_sd。
3. 实现 `libs/sd/hit.rs`：hit 数据结构。
4. 实现 `cmd_pgr/sd/search.rs`：单基因组 putative SD 检测。
5. 验证: 在相同参数下，PGR 输出的 hit 数量与坐标与 BISER `search` 子命令高度一致。

#### 第二阶段：比对精修（验证：hit 的 CIGAR 与 error rate 与 BISER align 一致）

1. 实现 `libs/sd/anchor.rs`：10-mer anchor 生成。
2. 实现 `libs/sd/refine.rs`：
   - 用 y 坐标离散化 + segment tree / Fenwick tree 实现 BISER 风格的 event-driven PST chaining：
     扫描线按 x 坐标推进，右端点事件时以 `dp[i] - gap_to_end` 单点更新 y 位置；左端点事件时先
     在 y 区间 `[ay - MAX_CHAIN_GAP, ay - 1]` 查询最大前驱得分，再计算 gap 惩罚 `dx + dy` 更新
     `dp[i]`，最后把滑出 `MAX_CHAIN_GAP` 窗口的旧锚点对应 y 位置重置为 `-INF`。
     `chain_blocks` 的 `GapCalc` 取 `max(dq, dt)`，不能 bit-exact 匹配 BISER。
   - 小 gap（≤1000 bp）用 `ScalarAlignmentEngine`（把一条 mate 构建成线性 POA graph），
     并自己实现 `path_to_cigar` 从 `Alignment.path` 生成 `=`/`X`/`I`/`D` CIGAR。
   - 大 gap 按 BISER 策略处理（两端各比对 1000 bp，中间用 `I`/`D`）。
   - 用 `paf/cigar.rs::format_cigar` 输出 CIGAR，用 `block_identity` 计算 BISER 风格的 error rate
     （`1 - block_identity`）；不要使用 `gap_compressed_identity`，因为它会低估 indel 错误。
3. 实现 `cmd_pgr/sd/align.rs`。
4. 验证: 对同一组 putative hits，PGR 与 BISER 输出的 alignment span、CIGAR、error rate
   差异 < 1%。

#### 第三阶段：聚类与分解（验证：elementary SD 集合与 BISER 一致）

1. 前置：将 `src/libs/paf/graph/dsu.rs` 的 `Dsu` 迁移到 `src/libs/ds/dsu.rs` 并公开为 pub，
   原文件通过 `pub use` 保持兼容。
2. 实现 `libs/sd/cluster.rs`：用 `ds::Dsu` 聚类并输出 cluster FASTA（复用 `fa::new_record`、
   `nt::rev_comp`、`loc::fetch_range_seq`、`twobit::TwoBitFile` 或 `paf/fasta.rs::FastaStore` 提取子序列）。
3. 实现 `libs/sd/decompose.rs`：多拷贝 plane-sweep。
4. 实现 `cmd_pgr/sd/cluster.rs` 与 `cmd_pgr/sd/decompose.rs`。
5. 验证: cluster 数量与覆盖范围、elementary SD 数量与 `.elem` 内容一致。

#### 第四阶段：core duplicon 与坐标转换（验证：CORE 标记与 translate 后坐标一致）

1. 实现 `libs/sd/set_cover.rs` 与 `libs/sd/cover.rs`（用 `coitrees` 做区间重叠查询）。
2. 实现 `libs/sd/translate.rs`（基于 `mask` 阶段记录的 uppercase run 列表做二分查找映射，
   不要误用 `alignment/coords.rs` 的 gapped-alignment 坐标函数）。
3. 实现对应子命令。
4. 验证: CORE 标记集合与 BISER 一致；translate 后的坐标与 CIGAR 可通过一致性检查。

#### 第五阶段：pipeline 与跨基因组（验证：多基因组输入与 BISER 最终输出一致）

1. 实现 `cmd_pgr/sd/run.rs` 串联所有步骤。
2. 实现 `cross_search` / `cross_align` 等价逻辑（复用 search/align，只是 query 为另一基因组）。
3. 验证: 多基因组场景下最终 `out.bed` 与 `.elem.txt` 与 BISER 等价。

### 6.7 风险与注意事项

- **比对器依赖（基本可解决，但性能有差异）**: BISER 的 alignment 精修依赖 `bio.seq.align()` 这个内置 SIMD aligner。
  PGR 的 `src/libs/poa/align.rs::ScalarAlignmentEngine` 提供 Global、Local、SemiGlobal 三种模式以及自定义
  match/mismatch/gap_open/gap_extend 参数。注意 `src/libs/poa/mod.rs` 只导出 `AlignmentParams`、`AlignmentType`、`Poa`；
  `ScalarAlignmentEngine` 和 `PoaGraph` 需分别通过 `poa::align::ScalarAlignmentEngine` 和 `poa::graph::PoaGraph` 使用。
  把其中一条 mate 构建成线性 POA graph 后即可做 pairwise alignment；小 gap（≤1000 bp）可直接复用 `libs/poa`。
  但它是标量 O(nm) 实现，无 SIMD/banded；在人类全基因组尺度批量调用可能显著慢于 BISER。后续若性能成为瓶颈，
  再评估 SIMD aligner（如 `parasail-rs`）或提取 `fas_multiz/banded_align.rs` 的 banded DP 为通用函数。
- **`chain_blocks` 不能直接复用于 BISER chaining**: `src/libs/chain/connect.rs::chain_blocks` 的 `cost_func` 是内部硬编码的，
  使用 `GapCalc::calc(dq, dt)`。`GapCalc` 在 `dq > 0 && dt > 0` 时取 `max(dq, dt)` 查表，而 BISER 要求 `dx + dy`。
  因此不能通过传入参数让 `chain_blocks` bit-exact 匹配 BISER；必须用 y 坐标离散化 + segment tree / Fenwick tree
  自行实现 BISER 的 event-driven PST chaining。现有 `src/libs/ds/kdtree.rs::KdTree` 不支持 deactivate，不能直接用。
- **BISER refine DP 的 gap 模型需自行实现**: BISER 的 sparse DP 使用特殊公式
  `MISMATCH * mi + GAPOPEN + GAP * (ma - mi)`，不能简单套用 `GapCalc::affine`。
- **`Dsu` 需要迁移到 `src/libs/ds/`**: `src/libs/paf/graph/dsu.rs::Dsu` 是 `pub(super)`，不能作为公共 API。
  根据项目约束，纯数据结构应放 `src/libs/ds/`。建议迁移到 `src/libs/ds/dsu.rs` 并公开，原文件通过 `pub use` 保持兼容。
- **`hash.rs` 不是 exact k-mer 索引**: `src/libs/hash.rs` 提供基于哈希的 canonical minimizer 采样（`seq_sketch`、`JumpingMinimizer`），
  而 BISER 的 search/decompose 依赖 exact 2-bit k-mer + winnowing，因此 k-mer 索引必须重新实现，
  `hash.rs` 仅可作为 sketch 验证或后续扩展使用。
- **坐标系统不一致**: PGR 内部不同模块使用不同坐标约定：
  - `chain`、`paf`、`twobit`、`FastaStore`、`BitMap`、`DupeTree`、`io::SequenceReader` 使用 0-based half-open。
  - `loc.rs` 的 `intspan::Range` 和 `slice_record` 使用 1-based inclusive，支持链向。
  - `fa::mask_sequence` 的 `IntSpan` 是 1-based inclusive；`fa::windows` 输出的坐标也是 1-based inclusive。
  - BISER 内部大量 0-based half-open。迁移时应在 SD 模块内部统一使用 0-based half-open，仅在调用 `loc` 或输出时转换。
- **hard-masked ↔ original 坐标映射不要误用 `alignment/coords.rs`**: `alignment/coords.rs` 是处理
  带 gap 对齐列的坐标转换，而 BISER 的 translate 只是“删除 lowercase 碱基”后的简单偏移映射。
  应在 `mask` 阶段记录 uppercase run 边界，在 `translate` 中做二分查找。
- **大染色体内存问题**: `fa::reader` 通过 `noodles_fasta` 顺序读取时，会把整条 record（完整染色体序列）
  读入内存（`Record.sequence()` 返回完整序列）。人类尺度染色体（~250 Mbp）尚可，但更大基因组或内存受限时，
  建索引阶段应改用 `src/libs/fmt/twobit.rs::TwoBitFile` 顺序扫描，或按 `MAX_CHROMOSOME_SIZE` 切片读取。
  区间提取 mate 序列时，优先使用 `loc` 或 `twobit`，避免加载完整染色体。
- **反向互补**: BISER 在 chromosome 级别同时索引 forward 与 reverse complement，
  PGR 中可用 `nt::rev_comp` 生成反向链序列，或按 BISER 方式在索引阶段同时扫描两条链。
- **性能**: BISER 的 Codon 实现经编译为原生代码，plane-sweep 是其性能关键。
  Rust 实现应可达到相近性能，但需对 k-mer 索引做内存优化（如 `u64` key + `Vec<(u32, u32)>`）。

## 7. 参考文献

- Išerić H, Alkan C, Hach F, Numanagić I.
  *Fast characterization of segmental duplication structure in multiple genome assemblies*
  . Algorithms Mol Biol. 2022;17:4.
  [https://doi.org/10.1186/s13015-022-00210-2](https://doi.org/10.1186/s13015-022-00210-2)
- Išerić H, Alkan C, Hach F, Numanagić I.
  *BISER: Fast Characterization of Segmental Duplication Structure in Multiple Genome Assemblies*
  . WABI 2021. LIPIcs, Vol. 201, 15:1–15:18.
  [https://drops.dagstuhl.de/opus/volltexte/2021/14368/pdf/LIPIcs-WABI-2021-15.pdf](https://drops.dagstuhl.de/opus/volltexte/2021/14368/pdf/LIPIcs-WABI-2021-15.pdf)
- Vollger M. R. et al. *Segmental duplications and their variation in a complete human genome*
  . Science. 2022;376:eabj6965.
  [https://doi.org/10.1126/science.abj6965](https://doi.org/10.1126/science.abj6965)
- Vollger M. R. et al. *Increased mutation and gene conversion within human segmental duplications*
  . Nature. 2023;617:335–344.
  [https://doi.org/10.1038/s41586-023-05895-y](https://doi.org/10.1038/s41586-023-05895-y)
- BISER GitHub Repository: [https://github.com/0xTCG/biser](https://github.com/0xTCG/biser)
