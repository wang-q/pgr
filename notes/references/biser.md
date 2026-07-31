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

PGR 内部各模块使用的坐标约定不统一，迁移时必须明确每个调用点的坐标系统，避免
`off-by-one` 或区间开闭错误。

- **0-based half-open `[start, end)`**
  - BISER 内部大量坐标（putative SD 区间、anchor 坐标、elementary SD 边界）使用 0-based half-open。
  - `src/libs/chain/record.rs` 的 `ChainHeader` / `Block`。
  - `src/libs/chain/connect.rs` 的 `ChainableBlock`。
  - `src/libs/paf/cigar.rs` 的 `slice_cigar_by_target` 与 `project` 函数。
  - `src/libs/paf/fasta.rs::FastaStore::fetch_range(start, end)`。
  - `src/libs/io.rs::SequenceReader::read_sequence(start, end)`。
  - `src/libs/fmt/twobit.rs::TwoBitFile::read_sequence(start, end, no_mask)`。
  - `src/libs/ds/bitmap.rs::BitMap` 与 `src/libs/ds/dupe_tree.rs::DupeTree`。
  - `src/libs/fmt/fa.rs::windows` 内部切片使用 0-based half-open，只是输出名称中把 `start` 转换为 1-based、把 `end` 作为 1-based inclusive 显示。

- **1-based inclusive `[start, end]`**
  - `src/libs/fmt/fa.rs::mask_sequence(seq, spans, hard)` 的 `spans` 是 `intspan::IntSpan`，1-based inclusive。
  - `src/libs/loc.rs` 的 `intspan::Range`（支持 `chr:start-end` 与 `chr(-):start-end`）。
  - `src/libs/loc.rs::slice_record` 与 `fetch_range_seq` 接收 1-based inclusive 的 `Range`。
  - `src/libs/alignment/coords.rs` 的 `chr_to_align` / `align_to_chr` 配合 `IntSpan` 使用 1-based inclusive。

- **SD 模块内部建议**
  - 在 `src/libs/sd/` 内部统一使用 0-based half-open，仅在调用 `loc`、`fa::mask_sequence`、输出 BED/文件名时做显式转换。
  - hard-masked 坐标 ↔ original 坐标的映射表也用 0-based half-open：`[(orig_start, orig_end, masked_start, masked_end)]`。

### 6.3 BISER 算法阶段与 PGR 组件映射

#### 6.2.1 Hard-masking（`mask.codon`）

- **BISER 实现**
  - 文件: `biser/codon/mask.codon:7-15`
  - 行为: 读取 FASTA，只保留 uppercase A/C/G/T，其余字符（包括 lowercase a/c/g/t、N、IUPAC ambiguity、gap 等）全部删除，按 `width=80` 输出为新的 hard-masked FASTA。
    因此 hard-masked 序列长度会变短，需要额外记录 uppercase run 在原基因组中的坐标，供后续 translate 使用。
- **PGR 可复用组件**
  - `src/libs/fmt/fa.rs`:
    - `reader(infile) -> Result<fasta::io::Reader<...>>`: 顺序读取 FASTA，支持 stdin 与 gzip。
    - `writer_with_wrap(outfile, 80) -> Result<fasta::io::Writer<...>>`: 按 80 bp 换行输出（与 BISER 的 `width=80` 一致）。
    - `find_masked_regions(seq, gap_only=false) -> Vec<(usize, usize)>`: 返回 lowercase 或 N/n 区域的 0-based inclusive 区间（`gap_only=true` 时只统计 N/n）。
    - `mask_sequence(seq, spans, hard=true) -> Result<String>`: 将 `spans`（1-based inclusive `IntSpan`）指定区间替换为 `N`（hard）或小写（soft）。注意它**保留序列长度**，与 BISER 的“删除 lowercase”行为不同，因此不能直接用于生成 hard-masked FASTA。
  - `src/libs/nt.rs`:
    - `NT_VAL: &[usize; 256]`: A/a/C/c/G/g/T/t/U/u 映射到 0/1/2/3，IUPAC ambiguity codes 与 N/n 映射到 4，其余为 255 (Invalid)。
    - `is_lower(b) -> bool`: 判断小写碱基。
    - `is_n(b) -> bool`: 判断 N 或 IUPAC ambiguous（含 M/R/W/S/Y/K/V/H/D/B）。
    - `rev_comp(seq) -> impl Iterator<Item = u8>`: 反向互补迭代器，在 cluster 阶段构造反向链序列时可直接使用。
- **需要新增的实现**
  - BISER 的 `mask` 是**删除 lowercase bases**并输出 hard-masked FASTA（序列长度变短），没有现成函数。
  - 实现方式：读取 record 时只保留 uppercase A/C/G/T（即 `NT_VAL[b] <= 3 && !b.is_ascii_lowercase()`）。小写 a/c/g/t、N、IUPAC ambiguity 以及 gap 字符均删除。
  - 同步记录每个 uppercase run 在原序列中的 `[orig_start, orig_end)` 边界以及对应 hard-masked 坐标 `[masked_start, masked_end)`，存入 `Vec<(orig_start, orig_end, masked_start, masked_end)>` 供 `translate` 使用。
  - 输出用 `fa::writer_with_wrap(outfile, 80)`。
  - 注意：`fa::reader` 会一次性将整条 record 读入内存（`noodles_fasta` 的 `Record.sequence()` 返回完整序列）。人类尺度染色体（~250 Mbp）尚可接受，但若要在更大基因组或内存受限场景下处理，建索引阶段可改用 `src/libs/fmt/twobit.rs::TwoBitFile` 顺序扫描，区间提取再用 `twobit` 或 `loc`。

#### 6.2.2 Putative SD detection（`search.codon`）

- **BISER 实现**
  - 文件: `biser/codon/search.codon:189-343`
  - 核心: 2-bit 滚动哈希 + winnowing + plane-sweep 链表 + tau 阈值 + 输出候选 hit。
- **PGR 可复用组件**
  - `src/libs/nt.rs`:
    - `NT_VAL: &[usize; 256]` 将 A/a/C/c/G/g/T/t/U/u 映射到 0/1/2/3，可直接用于 BISER 风格的 2-bit 滚动哈希；
      遇到 `NT_VAL[b] > 3` 时跳过（与 BISER hard-mask 后只保留 A/C/G/T 的行为一致）。
  - `src/libs/fmt/fa.rs`:
    - `reader()`: 顺序读取 FASTA，适合建 k-mer 索引。注意它会将整条 record 载入内存。
  - `src/libs/ds/bitmap.rs`:
    - `BitMap::new(size)` + `set_range(start, len)` + `is_fully_set(start, len)`: 0-based 位图，可用于标记 plane-sweep 或 decomposition 中已访问/已输出的基因组位置，避免同一碱基被重复命中。
- **不可直接复用**
  - `src/libs/hash.rs`: 提供 canonical minimizer 采样（`seq_sketch`、`JumpingMinimizer`）与 Jaccard/Mash 距离计算。
    但 BISER search/decompose 依赖 exact 2-bit k-mer + winnowing（非 canonical、非 hash-based），
    因此 `hash.rs` 的 minimizer 流程不能直接复用，仅可作为 sketch 验证或后续扩展使用。
- **需要新增的实现**
  - `src/libs/sd/kmer_index.rs`: exact 2-bit k-mer 滚动哈希、winnowing 采样、
    `kmer -> Vec<(chr_id, pos)>` 索引、频率阈值过滤（0.1%）。
  - `src/libs/sd/plane_sweep.rs`: `ListNode` 链表、`update_list()` 三种分支逻辑、
    `save_sd()` 输出、tau 计算。
  - `src/libs/sd/hit.rs`: SD hit 数据结构（坐标、species、chromosome、strand、CIGAR、
    error rate），可参考 `src/libs/chain/record.rs` 的 `Chain` / `Block` 设计。

#### 6.2.3 Alignment refinement（`align.codon` + `hit.codon`）

- **BISER 实现**
  - 文件: `biser/codon/align.codon:5-112`、`biser/codon/hit.codon:325-348`
  - 核心: 10-mer anchor 生成、PST chaining、sparse DP refine、CIGAR 精修。
- **PGR 可复用组件**
  - `src/libs/ds/kdtree.rs`:
    - `KdTree::build(indices, items)` + `update_scores(leaf_idx, score, items)` +
      `best_predecessor(target_idx, current_score, items, cost_func, lower_bound_func)` 是底层 chaining 引擎。
    - `KdTreeItem` trait 要求 `x_start/x_end/y_start/y_end/score`。在 `src/libs/chain/connect.rs` 的实现中
      `x` 对应 query、`y` 对应 target；BISER 的 anchor 可映射为
      `x_start=q_start, x_end=q_end, y_start=t_start, y_end=t_end`。
    - 由于 BISER chaining 的 gap 惩罚是 `dx + dy`（无 open 项），而 `chain_blocks` 的 `cost_func` 是内部硬编码、
      使用 `GapCalc::calc(dq, dt)` 的，因此**不能直接复用 `chain_blocks` 来严格匹配 BISER**。应直接使用 `KdTree`，
      在自定义的 `cost_func` 中返回 `None` 过滤不满足 `MAX_CHAIN_GAP` 的前驱，并按 `dx + dy` 计算惩罚。
  - `src/libs/chain/connect.rs`:
    - `chain_blocks(blocks, gap_calc, score_ctx, ...) -> Result<Vec<Chain>>` 是已经实现的完整 chaining DP，
      包含去重、merge、trim、score recalc。但它的打分模型面向 UCSC `axtChain`（`GapCalc` 取 `max(dq, dt)`、
      有 overlap trim 等），与 BISER 的 PST chaining 不完全等价。
    - 如果只想快速验证“KD-tree chaining”在 PGR 中是否可行，可用 `chain_blocks` 做原型；
      若要 bit-exact 匹配 BISER，必须基于 `KdTree` 自己实现。
  - `src/libs/ds/gap_calc.rs`:
    - `GapCalc::medium()` / `GapCalc::loose()` / `GapCalc::affine(open, extend)`: 预计算 gap cost 表。
    - **重要差异**: `GapCalc::calc(dq, dt)` 在 `dq > 0 && dt > 0` 时使用 `max(dq, dt)` 查表，而 BISER chaining
      要求 `dx + dy`。因此 BISER chaining 不能通过 `GapCalc` 表达，需要在 `KdTree::best_predecessor` 的
      `cost_func` 中直接按 `dx + dy` 计算。
    - BISER alignment refinement 的 sparse DP 使用 `GAPOPEN` / `GAP`，单轴 gap 可用 `GapCalc::affine(gap_open, gap_extend)`
      近似；但 BISER refine DP 对双 gap 的惩罚公式特殊（`MISMATCH * mi + GAPOPEN + GAP * (ma - mi)`），
      需要在新模块中重新实现，不能简单套用 `GapCalc`。
  - `src/libs/poa/align.rs` + `src/libs/poa/graph.rs` + `src/libs/poa/poa.rs`:
    - `src/libs/poa/mod.rs` 只导出 `AlignmentParams`、`AlignmentType`、`Poa`。`ScalarAlignmentEngine` 和 `PoaGraph`
      需分别通过 `poa::align::ScalarAlignmentEngine` 和 `poa::graph::PoaGraph` 使用。
    - `ScalarAlignmentEngine::new(AlignmentParams { match_score, mismatch_score, gap_open, gap_extend }, AlignmentType::Local)`
      提供 Smith-Waterman 局部比对；也支持 `SemiGlobal` 和 `Global`。
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
    - `Alignment.path: Vec<(Option<usize>, Option<NodeIndex>)>` 描述序列位置与 graph 节点的对齐关系。
      需要自己编写 `path_to_cigar(ref_seq, qry_seq, &aln.path) -> Vec<CigarOp>`：按 path 顺序遍历，
      同时推进 `ref_seq`（graph 节点碱基）和 `qry_seq`（序列索引），输出 `=`/`X`/`I`/`D`。
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
      统计与 identity 计算，可直接用于 error rate 估算。
    - `slice_cigar_by_target(cigar, target_start, ts, te)`: 按 target 子区间切 CIGAR。
    - `CigarOp` 使用 bit-packed `u32`：高 3 bits 存 op code（`=`/`X`/`I`/`D`/`M`），低 29 bits 存长度，最大单 op 长度约 512 Mbp。
  - `src/libs/alignment/stat.rs`:
    - `pair_d(seq1, seq2) -> Result<f32>`: 计算两条**等长对齐序列**的 divergence（忽略 gap 列，IUPAC 按 N 处理）。
      适合从展开后的对齐字符串计算，而不是直接用于 CIGAR。
    - `alignment_stat(seqs) -> Result<(i32, i32, i32, i32, i32, f32)>`: 多序列对齐列统计（长度、可比列、差异列、gap 列、模糊列、平均差异）。
  - `src/libs/alignment/msa.rs`:
    - `align_seqs(seqs, "builtin")` 调用内置 POA 做多序列对齐，返回 MSA 字符串（含 `'-'`），**不是 pairwise alignment**。
      也支持 `"spoa"`、`"clustalw"`、`"muscle"`、`"mafft"` 等外部 aligner。
    - 对于 pairwise 小 gap，优先直接用 `ScalarAlignmentEngine`；对于多拷贝 consensus 或 cluster MSA，
      可复用 `align_seqs(..., "builtin")` 或 `poa::Poa`。
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
  - `src/libs/sd/refine.rs`: 基于 `KdTree` 的 anchor chaining + sparse DP 精修、大 gap 处理、
    两端 score-based `ltrim`/`rtrim`、生成最终 CIGAR。小 gap（≤1000 bp）调用 `ScalarAlignmentEngine`
    并把 path 转为 CIGAR；大 gap 按 BISER 策略（两端各比对 1000 bp，中间用 `I`/`D`）处理。

#### 6.2.4 SD clustering（`cluster.codon`）

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
  - `src/libs/nt.rs`:
    - `rev_comp(seq)`: 生成反向互补序列，返回迭代器。
  - `src/libs/loc.rs`:
    - `open_indexed(infile, force_update) -> Result<(Input, IndexMap<name, (offset, size)>)>`: 打开带 `.loc` 索引的 FASTA
      （plain 或 BGZF 均可），不存在时自动创建索引。
    - `fetch_record(reader, loc_of, name) -> Result<fasta::Record>`: 按名字读取完整 record。
    - `fetch_range_seq(reader, loc_of, rg) -> Result<String>`: 按 `intspan::Range`（1-based inclusive，支持 `chr(-):start-end` 链向）
      提取子序列。
    - `slice_record(record, rg) -> Result<fasta::record::Sequence>`: 从已加载 record 中按 1-based Range 切片。
  - `src/libs/io.rs`:
    - `SequenceReader` trait: `read_sequence(name, start, end) -> Result<String>`，定义 0-based half-open 的随机访问接口。
      `TwoBitFile` 实现该 trait，`chain_blocks` 的 `ScoreContext` 也依赖它，因此 SD 流程中可用统一接口切换 2bit / indexed FASTA。
  - `src/libs/fmt/twobit.rs`:
    - `TwoBitFile::read_sequence(name, start, end, no_mask)`: 0-based half-open 随机访问 2bit 序列，
      适合区间提取；`no_mask=true` 时返回 uppercase，否则保留 mask（N-blocks 变 N，soft-mask 变 lowercase）。
    - `TwoBitFile` 实现 `SequenceReader` trait，可直接传给 `ScoreContext` 等需要序列读取的接口。
  - `src/libs/paf/fasta.rs`:
    - `FastaStore::new(seq_to_file)` + `fetch_range(name, start, end)` + `fetch_full(name)`: 管理多个 **BGZF FASTA**
      文件，带 `.loc` 索引与 LRU 缓存，适合多基因组 cross_search/cross_align 时批量提取 mate 序列。
    - **限制**: `FastaStore::new` 内部使用 `noodles_bgzf::io::indexed_reader`，因此输入必须是 BGZF 压缩的 FASTA。
      普通 gzip 或未压缩 FASTA 应先用 `loc`（plain/BGZF 通用）或 `twobit` 处理。
- **需要新增的实现**
  - `src/libs/sd/cluster.rs`: 将 hit 端点排序，用 `Dsu` 合并重叠端点，输出每个 cluster 的
    FASTA（序列名采用 `species#chrom+/-#start#end` 格式）。

#### 6.2.5 SD decomposition（`decompose.codon`）

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

#### 6.2.6 Core duplicon identification（`cover.py`）

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
      在将 elementary SD 区间映射到 SD 覆盖时思路类似。
- **需要新增的实现**
  - `src/libs/sd/set_cover.rs`: 用二叉堆实现 `greedy_set_cover`，输入为
    `elementary_id -> Vec<sd_id>`，输出被选中的 core elementary IDs。
  - `src/libs/sd/cover.rs`: 组装 interval overlap + set cover 流程，在 `.elem` 文件中追加
    `CORE` 标记。

#### 6.2.7 Coordinate translation（`mask.codon:32-154`）

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
- **需要新增的实现**
  - `src/libs/sd/translate.rs`: 在 `mask` 阶段记录 uppercase run 列表
    `[(orig_start, orig_end, masked_start, masked_end)]`（0-based half-open），通过二分查找实现
    hard-masked ↔ original 坐标双向映射，并同步改写 CIGAR。
  - 具体映射规则（需对照 BISER 源码精确实现）：
    - hard-masked 坐标下的 match/mismatch 区间，映射回 original 坐标时对应 uppercase 区段，CIGAR 保持 `=`/`X`/`M`。
    - 当 hard-masked 的 gap（`I`/`D`）跨越 original 中的 lowercase/masked 区域时，需要把该段替换为
      `S`（soft-mask）或 `N`（hard-mask），以反映原始基因组中这些碱基并非真正参与比对。
    - BISER 内部将 CIGAR 操作重新归类为 `M`/`S`/`N`；translate 阶段应产生与 BISER 输出语义一致的 CIGAR。

### 6.3 其他值得复用的工具

除了按算法阶段映射的组件外，以下通用工具在 SD 流程中也 likely 有用：

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
- `src/libs/paf/fasta.rs::FastaStore`
  - 多 **BGZF FASTA** 管理器，支持 `fetch_range(name, start, end)`（0-based half-open）与 `fetch_full(name)`，带 `.loc` 索引与 LRU 缓存。
  - 用途: 多基因组 cross_search/cross_align 时批量、高效地提取 mate 序列。注意输入必须是 BGZF 压缩 FASTA。
- `src/libs/lastz.rs`
  - 提供 UCSC lastz 预设（`PRESETS`、`find_preset`、`run_lastz`）与评分矩阵。可作为 `ScalarAlignmentEngine` 性能不足时的外部局部比对 fallback。
- `src/libs/chain/stitch.rs`
  - 按 chain ID 合并 fragments（`chainStitchId` 语义），要求 fragments 之间不重叠。SD 流程若输出带 ID 的 chain fragments，可参考；但它不是通用的“合并相邻 chain”工具。
- `src/libs/alignment/trim.rs`
  - 提供 MSA 后处理的 trim 函数。SD refine 阶段若需要对齐后修剪，可参考，但 BISER 的 score-based `ltrim`/`rtrim` 仍需自行实现。

### 6.4 建议的模块与命令结构

#### 6.4.1 `src/libs/sd/` 目录（新增）

- `kmer_index.rs`: exact 2-bit k-mer 索引 + winnowing + 频率过滤。
- `plane_sweep.rs`: plane-sweep 链表与 hit 输出。
- `hit.rs`: SD hit 数据结构（坐标、species、strand、CIGAR、error rate）。
- `anchor.rs`: 10-mer exact-match anchor 生成。
- `refine.rs`: 基于 `KdTree` chaining + sparse DP 的比对精修；包含 `path_to_cigar` 辅助函数。
- `cluster.rs`: 重叠 hit 聚类并输出 cluster FASTA。
- `decompose.rs`: elementary SD 分解。
- `set_cover.rs`: 贪心 set cover。
- `cover.rs`: core duplicon 标记流程。
- `translate.rs`: hard-masked 与原基因组坐标互转。注意 `src/libs/translate.rs` 已存在（蛋白质翻译），
  新增 SD 的 `translate.rs` 位于 `src/libs/sd/translate.rs`，不会冲突。

#### 6.4.2 `src/cmd_pgr/sd/` 目录（新增）

- `mod.rs`: 子命令注册与分发。
- `mask.rs`: `pgr sd mask <genome.fa> -o <masked.fa>`。
- `search.rs`: `pgr sd search <genome.fa> -o <hits.bed>`。
- `align.rs`: `pgr sd align <genome.fa> <hits.bed> -o <hits.align.bed>`。
- `cluster.rs`: `pgr sd cluster <genomes...> <hits.align.bed> -o <clusters.dir>`。
- `decompose.rs`: `pgr sd decompose <cluster.fa> -o <cluster.elem.bed>`。
- `cover.rs`: `pgr sd cover <hits.align.bed> <elems.txt> -o <elems.covered.txt>`。
- `translate.rs`: `pgr sd translate <hits.align.bed> <genomes...> -o <out.bed>`。
- `run.rs`: `pgr sd run <genomes...> -o <out.bed>`，按 BISER 顺序串接上述步骤。

### 6.5 分阶段实施计划

#### 第一阶段：索引与 plane-sweep（验证：human chr21 自比对命中数与 BISER 一致）

1. 实现 `libs/sd/kmer_index.rs`：2-bit k-mer + winnowing + 频率过滤。
2. 实现 `libs/sd/plane_sweep.rs`：ListNode、update_list、tau、save_sd。
3. 实现 `libs/sd/hit.rs`：hit 数据结构。
4. 实现 `cmd_pgr/sd/search.rs`：单基因组 putative SD 检测。
5. 验证: 在相同参数下，PGR 输出的 hit 数量与坐标与 BISER `search` 子命令高度一致。

#### 第二阶段：比对精修（验证：hit 的 CIGAR 与 error rate 与 BISER align 一致）

1. 实现 `libs/sd/anchor.rs`：10-mer anchor 生成。
2. 实现 `libs/sd/refine.rs`：
   - 用 `KdTree` 直接实现 BISER 风格的 anchor chaining（`cost_func` 中按 `dx + dy` 计算 gap 惩罚，
     并通过返回 `None` 过滤超出 `MAX_CHAIN_GAP` 的前驱）。`chain_blocks` 的 `GapCalc` 取 `max(dq, dt)`，不能 bit-exact 匹配 BISER。
   - 小 gap（≤1000 bp）用 `ScalarAlignmentEngine`（把一条 mate 构建成线性 POA graph），
     并自己实现 `path_to_cigar` 从 `Alignment.path` 生成 `=`/`X`/`I`/`D` CIGAR。
   - 大 gap 按 BISER 策略处理（两端各比对 1000 bp，中间用 `I`/`D`）。
   - 用 `paf/cigar.rs::format_cigar` 输出 CIGAR，用 `gap_compressed_identity` / `block_identity` 计算 error rate。
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

### 6.6 风险与注意事项

- **比对器依赖（基本可解决，但性能有差异）**: BISER 的 alignment 精修依赖 `bio.seq.align()` 这个内置 SIMD aligner。
  PGR 的 `src/libs/poa/align.rs::ScalarAlignmentEngine` 提供 Global、Local、SemiGlobal 三种模式以及自定义
  match/mismatch/gap_open/gap_extend 参数。注意 `src/libs/poa/mod.rs` 只导出 `AlignmentParams`、`AlignmentType`、`Poa`；
  `ScalarAlignmentEngine` 和 `PoaGraph` 需分别通过 `poa::align::ScalarAlignmentEngine` 和 `poa::graph::PoaGraph` 使用。
  把其中一条 mate 构建成线性 POA graph 后即可做 pairwise alignment；小 gap（≤1000 bp）可直接复用 `libs/poa`。
  但它是标量 O(nm) 实现，无 SIMD/banded；在人类全基因组尺度批量调用可能显著慢于 BISER。后续若性能成为瓶颈，
  再评估 SIMD aligner（如 `parasail-rs`）或提取 `fas_multiz/banded_align.rs` 的 banded DP 为通用函数。
- **`chain_blocks` 不能直接复用于 BISER chaining**: `src/libs/chain/connect.rs::chain_blocks` 的 `cost_func` 是内部硬编码的，
  使用 `GapCalc::calc(dq, dt)`。`GapCalc` 在 `dq > 0 && dt > 0` 时取 `max(dq, dt)` 查表，而 BISER 要求 `dx + dy`。
  因此不能通过传入参数让 `chain_blocks` bit-exact 匹配 BISER；必须用 `src/libs/ds/kdtree.rs::KdTree` 自行实现 chaining。
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
