# BISER 源码与论文分析

> 整理于 2026-07，源自对 `biser-master/` 目录源码及 published paper 的通读。目的：理解 BISER 在
> segmental duplication (SD) 检测与分解中的算法设计，并为 pgr 中重复/同源区域分析提供参考。

> **实施状态（2026-08-02）**：§6.6 第一阶段已落地——`pgr sd search`（LASTZ-based putative SD
> 检测）已实现并验证，见 §6.6 第一阶段；§6.8 的 chain/net 链路改用 `pgr pl chainnet`（原生实现，
> 与 UCSC 字节级一致，替代有 Linux 崩溃风险的 `pgr pl ucsc`）。

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

> **注**：BISER 的 error model 面向低至 75% 同一性的古老 SD。PGR 采用 T2T-CHM13 SD 标准（> 1 kbp、
> > 90% 同一性，见 4.2.1），不追求 75% 同一性场景；上述模型仅作为 BISER 算法背景理解，不作为 PGR
> 的检测阈值。

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

**2-bit 编码**: 代码中对碱基做 `(int(si) & 3)` 得到 2-bit 值，即 `A=0, C=1, G=2, T=3`。 滚动哈希
`h = ((h << 2) | base) & ((1 << 28) - 1)`（`KMER_SIZE=14`）生成 28-bit k-mer hash。

**Winnowing 实现**: `build_index()` 中用单调栈/队列维护 `(hash, pos)`：

- 新 k-mer 入队前，从队尾弹出 hash 不小于新 hash 的元素，保证队首为窗口最小值。
- 从队首弹出位置超出当前窗口 (`pos < i - KMER_SIZE + 1 - WINNOW_SIZE`) 的元素。
- 窗口填满 (`i - KMER_SIZE + 1 >= WINNOW_SIZE`) 后，队首 hash 即为一个 fingerprint； 只有
  fingerprint 变化时才进入 plane-sweep/索引插入，避免同一 hash 连续处理。

**Plane-sweep 链表 `update_list()`**: 维护按 `(chr, first)` 排序的 `ListNode` 链表， 每个节点保存：

- query 区间 `(first, last)` 与 reference 对应区间 `(ref, ref_last)`；
- 已扫描步数 `age`、命中次数 `count`；
- 是否曾经满足阈值 `potentional`。

对当前 query 位置 `current.loc` 对应的 reference 位置列表 `loci`（已按 `(chr, loc)` 排序），
`update_list()` 依次处理三种情况：

1. **延伸**: 若 `loci[lidx]` 与某 walker 同属一条 chromosome，且距离在 `MAX_DISTANCE`
   内 (`loci.loc - MAX_DISTANCE ≤ walker.last < loci.loc`)，则延伸 walker 的 `last`/`ref_last`
   并增加 `count`。
2. **插入**: 若 `loci[lidx]` 落在当前 walker 与下一个 walker 之间，则插入新节点， 年龄初始为 0
   （本轮结束时再 `age += 1`）。
3. **老化**: 若当前 walker 未被延伸或插入，则 `age += 1`，并检查 `count ≥ ceil(age · τ)`。
    - 若满足且长度 `< MAX_SD_LEN`，标记 `potentional`。
    - 若不满足或长度超限，且此前为 `potentional`，并满足 `last - first > QUERY_THRESHOLD` 与
      `current.loc - walker.ref ≥ REF_THRESHOLD`，则调用 `save_sd()` 输出该 hit。

**Tau 计算**: `tau()` 中先算
`gap_error = min(1.0, (MAX_ERROR - MAX_EDIT_ERROR) / MAX_EDIT_ERROR * MAX_EDIT_ERROR)`，再返回
`((1 - gap_error) / (1 + gap_error)) * (1 / (2 * exp(KMER_SIZE * MAX_EDIT_ERROR) - 1))`。
此即论文中的 Jaccard 下界。

**输出 `save_sd()`**: 输出前对边界做 `MAX_EXTEND` 填充（同 chromosome 时避免两个 mate 重叠），
并过滤掉完全自重叠（same species/name/strand 且坐标相同）的情况。最终生成 `Hit` 存入按
`(species1, chr1, species2, chr2, complement)` 分组的 `result` 字典。

#### 3.3.2 Chaining in `chain.codon`

`PrioritySearchTree` 是 chaining 的核心数据结构。它基于 y 坐标建静态二叉搜索树，每个叶子对应一个
anchor 的 y 位置 `(y + l - 1, anchor_idx)`。支持：

- `activate(x, score)`: 激活一个点并赋值。
- `deactivate(x)`: 将点置为 `-INF`。
- `rmq(lo, hi)`: 在区间 `[lo, hi]` 内找得分最高的激活点。

`chain()` 流程：

1. 对 reference 上的每个 anchor 左端点 `x` 和右端点 `x + l` 生成事件，按 x 排序扫描。
2. 在左端点事件时：
    - 先 deactivate 距离过远（`x - (anchor_j.x + anchor_j.l) > MAX_CHAIN_GAP`）的旧 anchor。
    - 在 PST 中查询 y 区间 `[ay - MAX_CHAIN_GAP, ay - 1]` 内得分最高的已激活 anchor `j`。
    - 若 `j` 存在，计算 gap cost `dx + dy`（`dx = ax - (jx + jl)`，`dy = ay - (jy + jl)`），
      更新 `dp[i] = max(MATCH_SCORE + (MATCH_SCORE//2)*(al-1) + dp[j] - gap, 自身得分)`，并设置
      `prev[i] = j`。
3. 在右端点事件时，以 `dp[i] - gap_to_end` 激活当前 anchor，其中 `gap_to_end` 是从当前 anchor 到
   序列末端的最大可能距离惩罚。
4. 按 DP 得分排序并重构链，过滤掉 span 小于 `MIN_UPPERCASE_MATCH` 或 `MIN_READ_SIZE * (1 - MAX_ERROR)`
   的链。

anchor 得分模型：第一个碱基得 `MATCH_SCORE`（默认 4），每多一个匹配碱基加 `MATCH_SCORE // 2`（即
2）。gap 惩罚为 reference gap 与 query gap 之和（`dx + dy`），无 open 项。

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
`kmer.as_int()` 作为 key，记录每个 k-mer 出现的 `(chr_id, loc)` 列表。同样用频率阈值过滤最高的 0.1%
k-mer。

**`update_list()` 与 search 的差异**: search 阶段每个节点只跟踪一段 putative SD 的两个 mate；
decompose 阶段需要跟踪同一 elementary SD 在多个序列上的多个拷贝，因此每个 `ListNode` 额外维护：

- `mappings: Dict[int, int]`：记录该节点对应 elementary SD 在每个 chromosome 上的当前最右边界；
- `gap`：自上次命中以来的未命中步数；
- `score`：命中计数。

处理当前 k-mer 的位置列表 `index`（按 `(chr, loc)` 降序排列）时：

1. 对尚未进入链表的新位置，在链表头部插入新节点，并继承当前 `mappings`。
2. 对可延伸的节点，更新 `end`、`score`、`gap=0`、`count += 1`。
3. 对老化节点（`gap >= diff=50`），若此前为 `potentional` 且长度超过 `MIN_MATE_LEN`， 调用
   `process()` 输出：
    - 若 `mappings` 为空，输出 `(chr, begin, end, score)`；
    - 否则对每个 chromosome 输出从 `begin` 到 `mappings[chr]` 的区间，并更新 `begin` 为
      `mappings[chr]+1`。
4. 清空已输出区域对应的 `visited` 标记，避免同一碱基被多次分解。

**合并 `merge()`**: 对相邻的 elementary SD 集合，若它们在所有 chromosome 上连续且间隔不超过 500bp，
则合并为一个集合。

**输出**: 每个 elementary SD 集合输出为 BED 行，格式为
`species\tchrom\tbegin\tend\tset_id\tlength\tscore\tstrand`；其中 strand 来自序列名中的 `+`/`-`。

## 4. 相关研究：T2T-CHM13 SD 注释来源与下游分析

本节记录 T2T-CHM13 基因组中 SD 注释的产生方式，以及基于这些注释开展的 SNV / IGC 研究。

### 4.1 T2T-CHM13 SD 注释的算法来源（Vollger et al. 2022）

> Vollger M. R. et al. *Segmental duplications and their variation in a complete human genome* . Science.
> 2022;376:eabj6965.[https://doi.org/10.1126/science.abj6965](https://doi.org/10.1126/science.abj6965)

T2T-CHM13 v1.1 的 SD 注释直接继承自 Vollger et al. 2022 对完整 T2T 基因组所做的 SD 注释。
该注释的生成流程如下：

1. **输入**: T2T-CHM13 v1.0 组装，并拼接上 GRCh38 的 chrY（用于包含 Y 染色体 SD 信息）。
2. **屏蔽重复序列**: 使用****与**RepeatMasker**对组装进行 masking，
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

1. 用**rustybam v0.1.29**处理 PAF：
    - `trim-paf`: 去除 query 上的冗余比对。
    - `break-paf`: 在 `>10 kbp` 的结构变异处断开比对。
2. 保留连续比对长度 `>1 Mbp` 的区块作为**syntenic 1:1 alignment**。

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
    - **Vollger et al. 2023**: 基于 T2T-CHM13 注释（> 1 kbp, > 90% 同一性）✅ **PGR 采用此标准**
- **关键算法**
    - **BISER**: ordered Jaccard + plane sweep + chaining
    - **Vollger et al. 2023**: minimap2 1:1 alignment + 窗口化重比对 + 二项检验
- **输出**
    - **BISER**: SD pairs / elementary SDs / core duplicons
    - **Vollger et al. 2023**: SNV 密度、IGC 事件、hotspots、突变谱

## 5. 对 pgr 的启示

1. **SD 同一性标准采用 T2T-CHM13 定义**: PGR 采用 > 1 kbp、> 90% 同一性的 SD 定义（见 4.2.1），
   不追求 BISER 默认的 75% 同一性（30% 错误率）场景。BISER 的 ordered Jaccard + plane-sweep
   虽然支持低同源性检测，但该能力不在 PGR 当前需求范围内；若未来确需检测古老重复，可再参考其
   colinear k-mer matching + 扫描线设计。
2. **Winnowing 作为采样手段**: BISER 用 winnowing 将 k-mer 采样率降到约 `2/(w+1)`，
   同时保证同一 windows 内必然命中。这与 pgr sketching 的有损采样不同，
   但为“在保敏感度的前提下降低索引规模”提供了可借鉴的参数化方法。pgr 现已落地 closed syncmer
   （`src/libs/syncmer.rs`）：同属有界间隔采样（连续点距 ≤ `2(w-1)`、密度约 `2/(w+1)`），且
   canonical 哈希使其链向对称，Jaccard/Mash 距离偏差小于 minimizer （Edgar 2021）；未来原生 BISER
   search 可优先以其替代 winnowing。
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
> 具体算法映射、可复用组件清单、缺失模块清单以及分阶段实施计划。所有可复用组件均与代码库中的 实际
> API 对应，目标是让迁移工作可以直接按图索骥。

### 6.1 迁移目标与边界

- **目标**: 在 PGR 中新增 `pgr sd` 命令族，实现与 BISER 等价的功能： putative SD 检测、局部比对精修、
  跨基因组映射、SD 聚类、elementary SD 分解、core duplicon 识别、hard-mask 坐标与原基因组坐标互转。
- **策略调整（简化）**: 不从头实现 BISER 的 k-mer 索引 + plane-sweep 搜索，而是先用 LASTZ
  （或基于 LASTZ 的覆盖度方法）生成候选 SD 区间，再经 UCSC chain/net 流程（`pgr pl ucsc`）
  做 chaining/refine，最后接入 BISER 后续的 cluster / decompose / translate 算法。
  **chaining 统一使用 UCSC chain/net，不用 lastz/FastGA 内置的 chaining，也不实现 BISER 的 PST refine**；
  SD 场景下 `pgr pl ucsc` 不加 `--syn`（不做共线性筛选），以保留伴随重排的 SD。chain/net 输出的
  MAF 经 `pgr maf to-paf` 转为 PAF 后进入下游。这样可以用 PGR 已成熟的 lastz/UCSC 封装快速验证端到端
  pipeline，同时保留未来替换为原生 k-mer search 的接口空间。
- **边界**: 本次迁移聚焦算法实现与命令接口；多进程调度、临时目录管理、 resume 等工程特性
  可在核心算法稳定后按 PGR 已有模式（如 `pgr pl` pipeline）补充。
- **原则**: 复杂算法放 `src/libs/`，`cmd_pgr/` 仅做参数解析、I/O 转换与调用。单命令专用的
  复杂逻辑也放 `libs/`。

### 6.2 坐标系统约定

PGR 内部不同模块混用 0-based half-open 与 1-based inclusive 两种约定。迁移前必须逐 个核对调用点，
否则 `mask`、`translate`、`loc` 三个环节极易出现 `off-by-one` 或区间 开闭错误。

#### 6.2.1 0-based half-open `[start, end)`

- **BISER 内部**: putative SD 区间、anchor 坐标、elementary SD 边界、hard-masked 坐标 均使用
  0-based half-open。
- **`src/libs/chain/record.rs`**: `ChainHeader` 与 `Block` 的 `t_start/t_end/q_start/q_end` 均为
  0-based half-open；`ChainData` 的 `size/dt/dq` 也是相对增量。
- **`src/libs/chain/connect.rs`**: `ChainableBlock` 的 `t_start/t_end/q_start/q_end` 为
  0-based half-open；`chain_blocks` 内部 gap 计算 `dt = target.t_start - cand.t_end`、
  `dq = target.q_start - cand.q_end`，负值表示重叠。
- **`src/libs/ds/kdtree.rs`**: `KdTreeItem` 要求 `x_start/y_start` 为 0-based inclusive，
  `x_end/y_end` 为 0-based exclusive。自定义 anchor 类型实现该 trait 时务必注意这一 开闭差异。
- **`src/libs/paf/cigar.rs`**: `slice_cigar_by_target(cigar, target_start, ts, te)` 的 `ts/te` 为
  0-based half-open；`project` 内部投影也按 half-open 处理。
- **`src/libs/paf/record.rs`**: `PafRecord` 的 `query_start/target_start` 为 0-based inclusive，
  `query_end/target_end` 为 0-based exclusive，符合 PAF 规范。
- **`src/libs/paf/fasta.rs::FastaStore::fetch_range(name, start, end)`**: `start/end` 为 0-based
  half-open 的 `i32`，函数内部转成 `noodles` 的 1-based inclusive position。
- **`src/libs/io.rs::SequenceReader::read_sequence(name, start, end)`**: `start/end` 为
  `Option<usize>`，0-based half-open；`None` 表示从头/到尾。
- **`src/libs/fmt/twobit.rs::TwoBitFile::read_sequence(name, start, end, no_mask)`**: 参数 为
  0-based half-open；`no_mask=true` 时返回 uppercase，否则保留 soft-mask。
- **`src/libs/ds/bitmap.rs::BitMap`**: `set_range(start, len)` 与 `is_fully_set(start, len)` 都按
  0-based half-open `[start, start+len)` 解释。
- **`src/libs/ds/dupe_tree.rs::DupeTree`**: `add(start, end)`、`subtract(start, end)`、
  `count_over(start, end, threshold)` 均为 0-based half-open。
- **`src/libs/alignment/coords.rs`**:
    - `reverse_range<T: Copy + Sub>(start: &mut T, end: &mut T, size: T)`: 原地反转 0-based
      half-open 区间，`new_start = size - old_end`, `new_end = size - old_start`。
    - `reverse_range_pair(start, end, size) -> (T, T)`: 非原地版本，返回反转后的 0-based half-open
      区间。
    - `reverse_range_1based(start: &mut usize, end: &mut usize, size: usize)`: 原地反转 1-based
      inclusive 区间，`new_start = size - old_end + 1`, `new_end = size - old_start + 1`。
    - `reverse_range_1based_pair(start, end, size) -> (usize, usize)`: 非原地版本。
    - `indel_intspan(seq: &[u8]) -> IntSpan`: 从带 gap 的对齐序列生成 1-based inclusive 的 gap
      位置集合。
    - `seq_intspan(seq: &[u8]) -> IntSpan`: 从带 gap 的对齐序列生成 1-based inclusive 的非 gap
      位置集合（即 `indel_intspan` 的补集）。
    - `chr_to_align(ints: &IntSpan, pos: i32, chr_start: i32, strand: &str) -> Result<i32>`:
      把基因组坐标（1-based inclusive）映射到对齐列坐标；`ints` 为 `seq_intspan` 结果，`chr_start`
      为该序列在染色体上的起始（1-based），`strand` 为 `"+"` 或 `"-"`。
    - `align_to_chr(ints: &IntSpan, pos: i32, chr_start: i32, strand: &str) -> Result<i32>`:
      把对齐列坐标映射回基因组坐标；当 `pos` 落在 gap 列时，会 pin 到左侧最近非 gap 碱基。
    - 以上函数均按 `size - ...` 计算，因此把 forward 区间映射到 reverse strand 坐标系，负链 block
      转换后 `start < end` 仍然成立。
    - **注意**: `chr_to_align` / `align_to_chr` 只适用于带 gap 的对齐序列与基因组坐标之间的转换；
      BISER hard-mask 后的坐标映射是简单偏移，需要 `sd/translate.rs` 自行实现。

#### 6.2.2 1-based inclusive `[start, end]`

- **`src/libs/fmt/fa.rs::mask_sequence(seq, spans, hard)`**: `spans` 是 `intspan::IntSpan`，1-based
  inclusive；函数内部通过 `offset = lower - 1` 转成切片 索引。该函数**保留序列长度**，与 BISER 的
  hard-mask（删除 lowercase）不同。现有命令 `pgr fa mask` 也是基于 runlist 做长度保留的 hard/soft
  mask，不能替代 BISER 的 `mask`。
- **`src/libs/fmt/fa.rs::find_masked_regions(seq, gap_only)`**: 返回 0-based inclusive 的
  `(begin, end)` 对，与 `mask_sequence` 的输入约定不同，不要直接混用。现有命令 `pgr fa masked`
  将该输出转换为 1-based inclusive 显示。
- **`src/libs/fmt/fa.rs::windows`**: 内部切片是 0-based half-open，但输出名称格式为
  `name:start-end`，其中 `start = 原 start + 1`（1-based inclusive），`end` 保持为 切片末尾位置
  （同样按 1-based inclusive 显示）。
- **`src/libs/loc.rs`**: `intspan::Range`（支持 `chr:start-end` 与 `chr(-):start-end`） 为
  1-based inclusive；`slice_record` 与 `fetch_range_seq` 接收该 `Range`。当 `start == 0` 时，
  `fetch_range_seq` 返回整条序列。`slice_record` 在负链时会对切片做 reverse complement。
- **`src/libs/alignment/coords.rs`**: `chr_to_align` / `align_to_chr` 的输入位置是 1-based
  inclusive，且配合 `IntSpan`（由 `seq_intspan` 从对齐序列的 gap 列生成）使用。这两个函数只适用于
  **带 gap 的对齐坐标**与基因组坐标之间的转换，不适用于 BISER hard-mask 后产生的简单偏移映射。

#### 6.2.3 SD 模块内部建议

- 在 `src/libs/sd/` 内部统一使用**0-based half-open**，仅在以下边界做显式转换：
    - 调用 `loc::fetch_range_seq` / `slice_record` 时，把 0-based half-open 区间转为
      `intspan::Range` 的 1-based inclusive；
    - 调用 `fa::mask_sequence` 时，把 0-based half-open 区间转为 `IntSpan` 的 1-based inclusive；
    - 输出 BED/文件名时，若需与 BISER 保持一致，再决定使用 0-based 还是 1-based。
- hard-masked 坐标 ↔ original 坐标的映射表建议保存为 0-based half-open：
  `Vec<(orig_start, orig_end, masked_start, masked_end)>`。

### 6.3 BISER 算法阶段与 PGR 组件映射

#### 6.3.1 Hard-masking（`mask.codon`）

- **BISER 实现**
    - 文件: `biser/codon/mask.codon:7-15`
    - 行为: 读取 FASTA，只保留 uppercase A/C/G/T，其余字符（包括 lowercase a/c/g/t、N、IUPAC
      ambiguity、gap 等）全部删除，按 `width=80` 输出为新的 hard-masked FASTA。因此 hard-masked
      序列长度会变短，需要额外记录 uppercase run 在原基因组中的坐标，供后续 translate 使用。
- **PGR 可复用组件**
    - `src/libs/fmt/fa.rs`:
        - `reader(infile) -> Result<fasta::io::Reader<Box<dyn BufRead>>>`: 顺序读取 FASTA，支持
          `stdin`、普通文件与 `.gz`。
        - `writer(outfile) -> Result<fasta::io::Writer<Box<dyn Write>>>`: 单行序列输出（不换行）。
        - `writer_with_wrap(outfile, line_base_count) -> Result<fasta::io::Writer<Box<dyn Write>>>`:
          按指定 bp 换行输出（与 BISER 的 `width=80` 一致）。
        - `writer_from_writer(writer) -> fasta::io::Writer<W>`: 把已有 `Write` 包装为单行 FASTA
          writer。
        - `new_record(name, seq) -> fasta::Record`: 从名称与序列字节构造 FASTA record。
        - `new_record_preserving_desc(name, source, seq) -> fasta::Record`: 构造 record 并保留
          `source` 的 description。
        - `find_fasta_files(path) -> Vec<PathBuf>`: 递归收集 `.fa` 与 `.fa.gz` 文件，
          输入为文件时返回单元素 vec。
        - `build_gzi_index(path) -> Result<()>`: 为 BGZF FASTA 构建 `.gzi` 索引，`FastaStore::new`
          需要该索引做随机访问。
        - `find_masked_regions(seq, gap_only) -> Vec<(usize, usize)>`: 返回 0-based inclusive 的
          masked 区间。`gap_only=false` 时返回 lowercase 或 `nt::is_n` 为 true 的字符（即 N/n 与
          IUPAC ambiguity）所在区间；`gap_only=true` 时只返回 N/n 与 IUPAC ambiguity 所在区间。
          **注意** lowercase 的 A/C/G/T 属于 `gap_only=false` 的返回范围，但它们不是 `nt::is_n`。
        - `mask_sequence(seq, spans, hard) -> Result<String>`: `seq` 为 `&str`，`spans` 为
          1-based inclusive 的 `IntSpan`；函数将区间内字符替换为 `N`（hard）或小写（soft）。它
          **保留序列长度**，与 BISER 的 hard-mask（删除 lowercase）行为不同，因此不能直接用于生成
          hard-masked FASTA。
    - `src/libs/nt.rs`:
        - `NT_VAL: &[usize; 256]`: 将 ASCII 字节映射到碱基编码。**A/a→0, C/c→1, G/g→2, T/t→3, U/u→3**
          （U/u 与 T/t 共用编码 3）；M/R/W/S/Y/K/V/H/D/B 及其小写，以及 N/n，映射到 4；其余字符
          （包括 gap `-`、`*` 等）映射到 255 (Invalid)。**关键注意**：lowercase a/c/g/t/u 在
          `NT_VAL` 中同样映射到 0/1/2/3，因此在做 BISER 风格的 2-bit 滚动哈希前，必须先完成
          hard-mask（删除 lowercase），不能直接用 `NT_VAL[b] & 3` 处理原始序列。
        - `is_n(b) -> bool`: 当 `NT_VAL[b] == 4` 时返回 true，即 N/n 与所有 IUPAC ambiguity codes。
          lowercase a/c/g/t/u 返回 false。
        - `is_lower(b) -> bool`: 判断字符是否为小写 ASCII。
        - `to_nt(nt) -> Nt`: 将字节映射到 `Nt` 枚举（A/C/G/T/N/Invalid）。
        - `count_n(seq) -> usize`: 统计 `is_n` 为 true 的字符数量。
        - `complement(seq) -> impl DoubleEndedIterator<Item = u8>`: 正序互补迭代器。
        - `rev_comp(seq) -> impl Iterator<Item = u8>`: 反向互补迭代器，在 cluster
          阶段构造反向链序列时可直接使用。该迭代器保留原字符大小写，因此 lowercase
          碱基在反向互补后仍为 lowercase。
    - `src/libs/fmt/twobit.rs::Blocks::from_dna`: 在打包 DNA 为 2bit 时，非 A/C/G/T 字符被记为
      N-block，lowercase A/C/G/T 被记为 soft-mask block。这是 2bit 写入时的 mask 语义，与 BISER
      hard-mask（删除字符）不同，但可作为参考。
        - **重要**：2bit 内部位编码为 `T=00, C=01, A=10, G=11`，而 BISER 的 2-bit 滚动哈希使用
          `A=0, C=1, G=2, T=3`。如果直接读取 2bit 的 packed bytes 做哈希，必须重新映射位；若通过
          `TwoBitFile::read_sequence` 读取字符串后再用 `NT_VAL` 编码，则天然得到 BISER 编码。
    - `src/libs/fmt/twobit.rs::TwoBitFile`:
        - 实现 `SequenceReader` trait 时固定调用 `read_sequence(..., no_mask=false)`，
          即通过 trait 接口读取会保留 soft-mask（返回 lowercase）和 N-block。SD
          流程若要从 2bit 获得 hard-masked 后的 uppercase 序列，应调用其 inherent 方法
          `TwoBitFile::read_sequence(name, start, end, no_mask=true)`，而不是 trait 方法。
        - `TwoBitFile::read_sequence` 的 `start/end` 为 `Option<usize>` 的 0-based half-open 区间；
          `no_mask=true` 时 soft-mask block 也会被转为 uppercase，`no_mask=false` 时保留 lowercase；
          N-block 始终返回 `N`。
- **需要新增的实现**
    - BISER 的 `mask` 是**删除 lowercase bases**并输出 hard-masked FASTA（序列长度变短），
      没有现成函数。
    - 实现方式：读取 record 时只保留 uppercase A/C/G/T（即 `NT_VAL[b] <= 3 && !b.is_ascii_lowercase()`）。
      小写 a/c/g/t、N、IUPAC ambiguity 以及 gap 字符均删除。
    - 同步记录每个 uppercase run 在原序列中的 `[orig_start, orig_end)` 边界以及对应 hard-masked 坐标
      `[masked_start, masked_end)`，存入 `Vec<(orig_start, orig_end, masked_start, masked_end)>` 供
      `translate` 使用。
    - 输出用 `fa::writer_with_wrap(outfile, 80)`。
    - 注意：`fa::reader` 会一次性将整条 record 读入内存（`noodles_fasta` 的 `Record.sequence()`
      返回完整序列）。人类尺度染色体（~250 Mbp）尚可接受，但若要在更大基因组或内存受限场景下处理，
      建索引阶段可改用 `src/libs/fmt/twobit.rs::TwoBitFile` 顺序扫描，区间提取再用 `twobit` 或
      `loc`。
    - 与现有 `pgr fa mask` / `pgr fa masked` 的关系：`pgr fa mask` 基于 runlist 做长度保留的
      hard/soft mask；`pgr fa masked` 只负责找出 masked 区域。BISER 的 `mask` 是独立功能，建议放在
      `pgr sd mask` 中实现。
- **与 LASTZ-based search 的关系**
    - 原生 BISER search 必须在 hard-masked 序列上建 k-mer 索引，因此 `mask` 是前置步骤。
    - 改用 LASTZ 后，`mask` 不再是 search 的前置条件：可以直接对原始（或 soft-masked）基因组跑
      lastz，再通过 RepeatMasker 差集或 PAF 过滤去除 TE 命中。
    - 但如果希望下游 `align` / `refine` / `cluster` 保持与 BISER 完全一致的 hard-masked 坐标系，
      仍建议实现 `pgr sd mask`，并基于 hard-masked 基因组做 lastz 自比对；此时 `translate`
      步骤负责把 hard-masked 坐标映射回原基因组。

#### 6.3.2 Putative SD detection（`search.codon`）

- **BISER 实现**
    - 文件: `biser/codon/search.codon:189-343`
    - 核心: 2-bit 滚动哈希 + winnowing + plane-sweep 链表 + tau 阈值 + 输出候选 hit。
- **简化策略：用 LASTZ 替代原生 search**
    - 基于当前实施计划简化要求，第一阶段不实现 BISER 的 k-mer 索引与 plane-sweep，而是复用 PGR
      已有的 lastz 基础设施生成候选 SD 区间。两种具体形态可选：
        - **形态 A：全基因组自比对（`lastz --self`）**
            - 直接调用 `pgr lav lastz --self <genome.fa> <genome.fa>`，输出 LAV；
            - 经 `pgr lav to-psl` 转为 PSL，得到 pairwise alignments（chaining/refine 交由 UCSC
              chain/net，见 6.3.3 与 6.8）；
            - 过滤短命中（< 1 kbp）和同一性 < 90% 的命中（T2T-CHM13 SD 标准，见 4.2.1），得到
              putative SD pairs。
        - **形态 B：滑动窗口覆盖度（`scripts/pgr-repeat.sh`）**
            - 把基因组切成 200 bp 重叠窗口，用 lastz 回贴到基因组；
            - 计算每个碱基的覆盖深度，取深度 ≥ 4 的区域作为重复区；
            - 用 RepeatMasker 区间做差集，过滤 TE；
            - 剩余区间作为候选 SD 区，再相互比对生成 pairs。
    - 形态 A 更直接，输出本身就是 pairwise hits，可跳过"候选区 → pairs"的二次匹配；形态 B 更稳健，
      能显式控制 TE 污染，但需要额外一步候选区自匹配。
- **PGR 可复用组件**
    - `src/libs/lastz.rs`:
        - `PRESETS` 提供 lastz 参数集；对 SD 自比对建议降低 `K/L`（如 `K=1500 L=1500`）以提高灵敏度，
          或直接使用 `--self` 模式。
        - `run_lastz` 负责并行调用 lastz 并输出 LAV。
    - `src/libs/fmt/lav.rs` / `src/cmd_pgr/lav/to_psl.rs`:
        - 解析 LAV 并转为 PSL，供 `pgr pl ucsc` 消费；PSL 不直接转 PAF，最终 PAF 由 UCSC chain/net
          输出的 MAF 经 `pgr maf to-paf` 产出（见 6.8、6.9）。
    - `src/libs/fmt/psl.rs`:
        - `PslRecord` 数据结构、`to_chain` 等方法可用于坐标转换与链向处理。
    - `src/libs/paf/cigar.rs` / `src/libs/paf/parser.rs`:
        - 若最终输出 PAF，可用 `parse_cigar`、`block_identity` 计算 error rate。
    - `src/libs/ds/bitmap.rs`:
        - `BitMap::new(size)` + `set_range(start, len)` + `is_fully_set(start, len)`: 0-based 位图，
          可用于标记已输出的基因组位置，避免同一碱基被重复命中；在 TE 差集与候选区合并时亦可复用。
- **需要新增的实现**
    - `src/libs/sd/search_lastz.rs`（或复用 `src/libs/sd/from_lastz.rs`）:
        - 封装"lastz 自比对 → 格式转换 → 过滤 → 输出 putative hits"的完整流程。
        - 输入：基因组 FASTA、lastz 参数、最小 hit 长度、最大 error rate。
        - 输出：原始 pairwise alignments（统一为 PSL：lastz LAV 经 `pgr lav to-psl`，FastGA/minimap2
          PAF 需 PAF→PSL），未经 chaining；最终供下游消费的 PAF 由 UCSC chain/net refine 阶段产出
          （见 6.3.3、6.8、6.9）。
    - `src/libs/sd/coverage.rs`（形态 B 需要）:
        - 封装 `pgr-repeat.sh` 中的窗口化、自比对、lift、覆盖度计算逻辑，输出候选重复区 BED。
    - `src/libs/sd/subtract_repeatmasker.rs`（形态 B 需要）:
        - 读取 RepeatMasker 输出，与候选区做差集。
    - `src/libs/sd/hit.rs`: SD hit 数据结构（坐标、species、chromosome、strand、CIGAR、error
      rate），可参考 `src/libs/chain/record.rs` 的 `Chain` / `Block` 设计。
- **原生 BISER search 的保留信息（未来可替换）**
    - 若后续要替换为原生 k-mer plane-sweep，需要实现 `src/libs/sd/kmer_index.rs` 与
      `src/libs/sd/plane_sweep.rs`。`src/libs/nt.rs::NT_VAL` 与 `src/libs/fmt/fa.rs::reader`
      可复用；`src/libs/hash.rs` 的 minimizer 流程与 BISER exact 2-bit k-mer 不同，不能直接复用。
      `src/libs/syncmer.rs` 的 DNA 路径已基于 2-bit canonical rolling hash（复用 `nt::NT_VAL`），
      编码与 BISER 一致，未来原生 search 的 k-mer 索引可优先基于它扩展。
    - 为保持接口统一，建议 `pgr sd search` 设计为 `--mode lastz|coverage|kmer`，默认 `lastz`， 未来
      `kmer` 模式输出格式与 `lastz` 模式完全一致。

#### 6.3.3 Alignment refinement（`align.codon` + `hit.codon`）

- **BISER 实现**
    - 文件: `biser/codon/align.codon:5-112`、`biser/codon/hit.codon:325-348`
    - 核心: 10-mer anchor 生成、PST chaining、sparse DP refine、CIGAR 精修。
    - **外部比对路线不实现本阶段**：lastz/FastGA 路线的 chaining/refine 统一交由 UCSC chain/net 流程
      （`pgr pl ucsc`，SD 不加 `--syn`），其 MAF 输出经 `pgr maf to-paf` 转为 PAF（见 6.8、6.9）。
      本节 PST refine 仅用于原生 BISER 路线，属延后实现项。
- **PGR 可复用组件**
    - `src/libs/ds/kdtree.rs`:
        - `KdTree::build(indices, items)` + `update_scores(leaf_idx, score, items)` +
          `best_predecessor(target_idx, current_score, items, cost_func, lower_bound_func)` 是底层
          chaining 引擎。
        - `KdTreeItem` trait 要求 `x_start/y_start` 为 0-based inclusive，`x_end/y_end`
          为 0-based exclusive；`score` 用于叶子初始得分。在 `src/libs/chain/connect.rs`
          的实现中 `x` 对应 query、`y` 对应 target；BISER 的 anchor 可映射为
          `x_start=q_start, x_end=q_end, y_start=t_start, y_end=t_end`。
        - **关键限制：`KdTree` 不支持 deactivate**。BISER 的 PST 在扫描锚点时需要按事件激活
          / deactivate 锚点（当锚点与当前扫描位置距离超过 `MAX_CHAIN_GAP` 时置为 `-INF`）。
          `KdTree::update_scores` 只会把叶子和祖先节点的 `max_score` 向上提升，不会向下衰减；
          一旦某个内部节点的 `max_score` 被设为高分，后续即使该子树下的所有叶子都应
          deactivate，该节点仍保留旧高分，导致剪枝边界和 `best_predecessor` 结果错误。因此
          **不能直接用 `KdTree` 实现 BISER 的 event-driven PST chaining**。
        - 可行的替代方案：
            1. **按扫描线顺序维护一个“当前窗口内活跃锚点”集合**，对该集合重建（或增量维护）
               一棵 KD-tree/Fenwick tree/线段树。由于 BISER 的 `MAX_CHAIN_GAP` 有限，
               窗口内锚点数量通常可控，但每次扫描线推进都重建 KD-tree 的复杂度是 `O(w log w)`，
               总复杂度会上升到 `O(n·w·log w)`，不适合大规模数据。
            2. **在 y 坐标离散化后使用线段树或树状数组（Fenwick tree）维护每个 y 位置的最大 DP 得分**，
               扫描线从左到右推进时，在 y 区间 `[ay - MAX_CHAIN_GAP, ay - 1]` 内查询最大值，
               并用单点更新写入当前锚点的 DP 值。这是实现 BISER PST 的最简洁路径，时间复杂度
               `O(n log n)`，且天然支持 deactivate（当锚点滑出窗口时将其对应 y 位置重置为 `-INF`）。
            3. 如果坚持使用 `KdTree`，只能用于“所有锚点同时激活、无 deactivate”的 chaining 场景，
               此时 gap 惩罚仍需在 `cost_func` 中按 `dx + dy` 计算，并通过返回 `None` 过滤超出
               `MAX_CHAIN_GAP` 的前驱；但这与 BISER 的 sweep + PST 逻辑不等价。
    - `src/libs/chain/connect.rs`:
        - `chain_blocks(blocks, gap_calc, score_ctx, ...) -> Result<Vec<Chain>>` 是已经实现的完整
          chaining DP，包含去重、merge、trim、score recalc。但它的打分模型面向 UCSC
          `axtChain`（`GapCalc` 取 `max(dq, dt)`、有 overlap trim 等），与 BISER 的 PST chaining
          不完全等价。
        - `ScoreContext { t_2bit, q_2bit, matrix }` 提供序列读取与替换矩阵，用于 overlap trim
          和最终 score recalc。当 `score_ctx` 为 `Some` 时，链构建完成并去重/merge 后，会调用
          `trim_overlaps`：对相邻 block 的 target 重叠区，用 `SubMatrix` 分别计算把重叠区全归左
          block 或全归右 block 的得分，取得分更高的切分点，然后重新计算链总得分。负链 query
          的序列会通过 `reverse_range_pair` 取反向互补后再参与评分。
        - **不建议用 `chain_blocks` 做 BISER 原型**：由于 gap 模型差异（`max(dq, dt)` vs
          `dx + dy`），其输出链与 BISER 会有系统性偏差，无法验证 PST chaining 逻辑是否正确。
          原型阶段应使用 y 离散化 + 线段树的 PST 实现。
    - `src/libs/chain/record.rs`:
        - `Chain { header: ChainHeader, data: Vec<ChainData> }`: UCSC chain 格式的内存表示。
        - `ChainHeader` 字段：`score`, `t_name`, `t_size`, `t_strand`, `t_start`, `t_end`, `q_name`,
          `q_size`, `q_strand`, `q_start`, `q_end`, `id`。`t_start/t_end/q_start/q_end` 为
          0-based half-open；`t_strand` 恒为 `'+'`，`q_strand` 可为 `'+'`/`'-'`。
        - `ChainData { size, dt, dq }`: 相对增量，最后一个 block 的 `dt=dq=0`。
        - `Chain::to_blocks() -> Vec<Block>`: 把相对 `ChainData` 转为绝对坐标
          `Block { t_start, t_end, q_start, q_end }`（0-based half-open）。
        - `Chain::from_blocks(header, blocks) -> Vec<ChainData>`: 从绝对 block 重建相对 data，
          并更新 header 的 `t_start/t_end/q_start/q_end`。
        - `Chain::subset(t_start, t_end) -> Option<Chain>`: 按 target 子区间切 chain，返回新的
          `Chain`。
        - `Chain::write(writer)`: 按 UCSC 格式写出。
        - `ChainReader<R>` / `read_chains(reader)`: 顺序读取 UCSC chain 文件；非 chain 行被忽略，
          `#` 行存入 `header_comments`。
        - BISER refine 阶段输出的 anchors/CIGAR 若需导出为 chain 格式（例如与 `pgr chain`
          生态互操作），可直接使用这些结构。
    - `src/libs/chain/sort.rs`:
        - `sort_chains(chains, renumber)`: 按 score 降序排序 chain；`renumber=true` 时从 1
          开始重新赋值 id。
    - `src/libs/chain/stitch.rs`:
        - `stitch_chains(reader, writer)`: 按 chain ID 合并 fragments（UCSC `chainStitchId` 语义）。
          要求同一 ID 的 fragments 在 target/query/strand 上一致且 block 不重叠；输出按 score 降序。
        - SD 流程若输出带 ID 的 chain fragments 并需要合并，可参考；但它不是通用的“合并相邻
          chain”工具。
    - `src/libs/ds/gap_calc.rs`:
        - `GapCalc::medium()` / `GapCalc::loose()` / `GapCalc::affine(open, extend)`: 预计算 gap
          cost 表。
        - 该类型也在 `src/libs/chain/mod.rs` 中通过 `pub use crate::libs::ds::GapCalc;` 重新导出为
          `chain::GapCalc`，`src/libs/fas_multiz/banded_align.rs` 等模块即通过此路径使用。
        - **重要差异**: `GapCalc::calc(dq, dt)` 在 `dq > 0 && dt > 0` 时使用 `max(dq, dt)` 查表，
          而 BISER chaining 要求 `dx + dy`。因此 BISER chaining 不能通过 `GapCalc` 表达，需要在
          segment-tree/Fenwick PST 的 查询更新逻辑中直接按 `dx + dy` 计算。
        - BISER alignment refinement 的 sparse DP 使用 `GAPOPEN` / `GAP`，单轴 gap 可用
          `GapCalc::affine(gap_open, gap_extend)`近似；但 BISER refine DP 对双 gap 的惩罚公式特殊
          （`MISMATCH * mi + GAPOPEN + GAP * (ma - mi)`），需要在新模块中重新实现，不能简单套用
          `GapCalc`。
    - `src/libs/poa/align.rs` + `src/libs/poa/graph.rs` + `src/libs/poa/poa.rs`:
        - `src/libs/poa/mod.rs` 只导出 `AlignmentParams`、`AlignmentType`、`Poa`。
          `ScalarAlignmentEngine` 和 `PoaGraph`需分别通过 `poa::align::ScalarAlignmentEngine` 和
          `poa::graph::PoaGraph` 使用。

        - `ScalarAlignmentEngine::new(AlignmentParams { match_score, mismatch_score, gap_open, gap_extend }, AlignmentType::Local) 提供 Smith-Waterman 局部比对；也支持  `SemiGlobal
          `和`Global`。

        - `Alignment { score, path }`：`score` 为最佳路径总得分；`path` 为最佳路径上的步骤序列。

        - `AlignmentParams::default()` 的默认值为
          `match_score=5, mismatch_score=-4, gap_open=-8, gap_extend=-6`。这些默认值来自 SPOA，
          适合做小片段 pairwise alignment；若要与 BISER 的 `MATCH_SCORE=4` 等参数对齐，应显式构造
          `AlignmentParams`。

        - **三种模式的精确语义（对线性 POA graph 做 pairwise 比对时）**：

            - `Local`：query 与 graph 均允许自由起点/终点，返回得分最高的局部子比对。适合在较大 gap
              区域内寻找最优局部对齐块。
            - `SemiGlobal`：query 必须完整对齐，graph 的起点/终点自由。适合把一段 query 锚定到
              reference 的任意子区间。
            - `Global`：代码注释为 "Needleman-Wunsch"，但实际实现是
              **query 完整对齐、graph 起点固定（从第一个节点开始）、graph 终点自由**。因此它
              **不是**传统 Needleman-Wunsch 的双端固定全局比对；若要做两条序列完全对齐，需要在得到
              alignment 后手动检查 path 是否覆盖 graph 首尾，或改用 `SemiGlobal` 并在后续截断。
              `align_seqs(..., "builtin")` 即使用此 `Global` 模式依次加入序列并输出 MSA。
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
        - `Alignment.path: Vec<(Option<usize>, Option<NodeIndex>)>` 的精确语义（path
          已经按正向顺序排列）：
            - `(Some(seq_idx), Some(node_idx))`：序列碱基 `seq_idx` 与 graph 节点 `node_idx`
              匹配/错配；
            - `(Some(seq_idx), None)`：序列碱基 `seq_idx` 在 graph 中对应位置为插入（CIGAR `I`）；
            - `(None, Some(node_idx))`：graph 节点碱基在序列中对应位置为删除（CIGAR `D`）。 按 path
              顺序遍历，同时推进 `ref_seq`（graph 节点碱基）和 `qry_seq`（序列索引），即可输出 `=`/
              `X`/`I`/`D` CIGAR。
        - 从 `Alignment.path` 生成 CIGAR 的模板（使用 `CigarOp::try_new`，`CigarOp::new` 为
          `pub(crate)` 不可外部调用）：
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
注意 `poa.add_sequence()` 会通过 `PoaGraph::add_alignment` 修改 graph，因此若只想做
          **不改变 graph 的 pairwise alignment**，应手动构建线性 `PoaGraph`。
        - `PoaGraph::add_alignment` 的具体行为（理解 consensus/MSA 的基础）：
            - 对 `Alignment.path` 中每个 `(seq_idx, node_idx)` 步骤，若 `seq_idx` 存在则消费 query
              碱基；`(None, Some)` 或 `(None, None)` 表示 deletion，不消费 query 碱基。
            - 在消费 query 碱基前，先把 path 中未对齐的 query 前缀（`(Some, None)` 之前的独立碱基）
              逐个加入 graph，作为新节点并用 weight=1 的边连接。
            - `(Some(seq_idx), Some(node_idx))` 且 graph 节点碱基与 query 碱基一致时，该节点
              `weight += 1`，predecessor 边 weight += 1。
            - `(Some(seq_idx), Some(node_idx))` 但碱基不一致时，先检查 `node_idx.aligned_nodes`
              中是否已有相同碱基的节点；若有则复用该节点，否则**新建一个节点**保存 query 碱基，
              并将其加入 `node_idx.aligned_nodes` 与原有 clique 合并；随后该节点 weight += 1。
            - 连续的 query 碱基通过 `add_edge` 连接，边权重累加；`add_edge` 对重复边做 weight
              累加而非新建多重边。
            - 因此 `PoaGraph` 中节点 weight 表示某碱基在已加入序列中出现的次数，边 weight
              表示相邻关系出现次数；`generate_consensus` 与 `generate_msa` 基于此做多数表决/回溯。
        - **性能注意**: `ScalarAlignmentEngine` 是标量 O(nm) 实现，无 SIMD/banded。小 gap（≤1000
          bp）可直接使用；人类全基因组尺度批量调用可能成为瓶颈，届时再评估 `parasail-rs` 或 banded
          优化。
        - 若需要多序列 consensus/MSA（例如 cluster 内多个拷贝），可用 `Poa::new(params, align_type)` +
          `add_sequence()` + `consensus()` / `msa()`，这比直接调 `ScalarAlignmentEngine` 更方便。
    - `src/libs/paf/cigar.rs`:
        - `parse_cigar(s) -> Result<Vec<CigarOp>>` / `format_cigar(ops) -> String`: CIGAR
          解析与格式化。
        - `extract_cigar(tags)`: 从 PAF tags 提取 CIGAR。
        - `reverse_cigar(ops)`: 反向并交换 I/D。
        - `cigar_from_alignment(ref, qry)`: 从两条**等长对齐序列**（含 gap 字符 `'-'`）生成 `=`
          / `X` / `I` / `D` CIGAR，比对大小写不敏感。不能直接从 `ScalarAlignmentEngine` 的
          `Alignment.path` 调用，需要先展开成等长对齐字符串。
        - `cigar_stats(ops) -> CigarStats` / `gap_compressed_identity(ops) -> f64`
          / `block_identity(ops) -> f64`:统计与 identity 计算。BISER
          的 error rate 是编辑错误率（`E/ℓ`），应基于 CIGAR 计算
          `(mismatches + ins_bp + del_bp) / (matches + mismatches + ins_bp + del_bp)`，等价于
          `1 - block_identity(ops)`。`gap_compressed_identity` 把每个 indel 只计为一个事件，会低估
          错误率，不适合直接作为 BISER error rate。
        - `slice_cigar_by_target(cigar, target_start, ts, te)`: 按 target 子区间切 CIGAR。
        - `CigarOp` 使用 bit-packed `u32`：高 3 bits 存 op code（`=`/`X`/`I`/`D`/`M`），低 29 bits
          存长度，最大单 op 长度约 512 Mbp。
    - `src/libs/alignment/stat.rs`:
        - `pair_d(seq1: &[u8], seq2: &[u8]) -> Result<f32>`: 计算两条**等长对齐序列**的 divergence。
          只统计 `NT_VAL` 均 ≤ 3 的位置（即 A/C/G/T，IUPAC ambiguity 按 N 排除）；忽略 gap 列；
          比较时忽略大小写。`comparable == 0` 时返回错误。
        - `alignment_stat(seqs: &[&[u8]]) -> Result<(i32, i32, i32, i32, i32, f32)>`:
          多序列对齐列统计，返回 `(length, comparable, difference, gap, ambiguous, mean_d)`。
            - `comparable`: 该列所有序列都是 A/C/G/T 的列数；
            - `difference`: comparable 列中至少有一个序列与第一个序列不同的列数；
            - `gap`: 该列包含 `'-'` 的列数；
            - `ambiguous`: 其余列（含 IUPAC ambiguity 或 N）；
            - `mean_d`: 所有序列对之间 `pair_d` 的平均值（序列数 < 2 时为 0.0）。
    - `src/libs/alignment/msa.rs`:
        - `align_seqs(seqs, "builtin")` 调用内置 POA 做多序列对齐，返回 MSA 字符串（含 `'-'`），
          **不是 pairwise alignment**。内部使用 `AlignmentType::Global` 将所有序列依次加入 POA
          graph，再调用 `poa.msa()` 输出。也支持 `"spoa"`、`"clustalw"`、`"muscle"`、`"mafft"`
          等外部 aligner。
        - `align_seqs_quick(seqs, aligner, pad, fill)` 在已有粗对齐（所有序列长度相同）
          的基础上，仅对 head/tail 和 gap 邻近区域调用外部 aligner 重新对齐，再拼回原位。
          适合先快速得到整体对齐框架，再局部精修的场景。
        - `get_consensus_poa_builtin(seqs, match_score, mismatch_score, gap_open, gap_extend, algo_code)`
          /`get_consensus_poa_external(seqs, ...)` 直接用 POA 或外部 spoa 生成 consensus 字符串。
        - 对于 pairwise 小 gap，优先直接用 `ScalarAlignmentEngine`；对于多拷贝 consensus 或 cluster
          MSA，可复用 `align_seqs(..., "builtin")`、`align_seqs_quick` 或 `poa::Poa`。
    - `src/libs/chain/sub_matrix.rs`:
        - `SubMatrix` 提供 256×256 字节替换矩阵（含大小写）与 gap open/extend。
        - `SubMatrix::default()` 是一个简化的 identity-like 矩阵：A/C/G/T 匹配得 100，错配 -100，
          N 相关 -100，`gap_open=400`, `gap_extend=30`。它**不是** lastz 默认矩阵，仅作为通用
          fallback。
        - `SubMatrix::hoxd55()` 才是 lastz 默认的 HoxD55 矩阵（A-A=91, C-C=100, G-G=100, T-T=91，
          非对角线负值），同样 `gap_open=400`, `gap_extend=30`。
        - `SubMatrix::from_name(name)` 支持 `"hoxd55"` 预设，其他名字按 BLAST 格式从文件解析。
        - `SubMatrix::get_score(c1, c2)` 按字符 ASCII 值查表，大小写均可。
        - `chain_blocks` 的 `ScoreContext` 使用 `SubMatrix` 在 overlap trim 时重新计算匹配得分。
          注意 `ScalarAlignmentEngine`**不使用**`SubMatrix`，它只接受简单的 `match_score`/
          `mismatch_score`；若 BISER 精修需要 HoxD55 等复杂矩阵，需改用外部 aligner（如 lastz）
          或自行实现支持替换矩阵的 DP。
    - `src/libs/lastz.rs`:
        - 提供 lastz 的预设评分矩阵与参数（`PRESETS`、`find_preset`、`run_lastz`）。若
          `ScalarAlignmentEngine`性能不足或需要更复杂的评分矩阵，可将小 gap 区域提取后调用 lastz
          作为外部 fallback。
    - `src/libs/fas_multiz/banded_align.rs`:
        - `banded_align_refs(...)` 对两个 `FasBlock` 的 reference entry 做 banded DP + affine gap，
          返回两 reference 序列列与列之间的对齐映射 `(Vec<Option<usize>>, Vec<Option<usize>>)`。
        - DP 会同时考虑两个 block 中所有共有的 species，对每个 `match/mismatch` 位置累加所有物种对的
          `SubMatrix::get_score` 得分（gap 只算一次），因此它本质上是多序列 banded 对齐，不是单纯的
          pairwise。
        - 当前实现紧密绑定 `FasBlock` 输入，不能直接复用于 BISER 的 pairwise gap 精修。若后续需要
          banded pairwise align，可参考其索引函数 `idx(i, j)` 和 band 半径计算逻辑，提取为通用函数。
        - 该模块使用 `crate::libs::chain::GapCalc`，它是 `src/libs/ds/gap_calc.rs::GapCalc` 的
          re-export。
    - `src/libs/alignment/trim.rs`:
        - `trim_pure_dash(seqs: &mut [String])`: 删除所有序列共同为 gap 的列（交集）。
        - `trim_head_tail(seqs: &mut [String])`: 从两端删除纯 gap 列，直到遇到任一序列的非 gap
          字符。
        - `trim_outgroup(seqs: &mut [String]) -> Result<()>`: 要求至少 3 条序列，最后一条为
          outgroup；删除 outgroup 有插入而 ingroup 共同为 gap 的区域（ingroup gap 的并集是 ingroup
          gap 交集的超集时才删）。
        - `trim_complex_indel(seqs: &mut [String]) -> Result<IntSpan>`: 在 `trim_outgroup` 后使用，
          识别并删除 ingroup 内部复杂 indel 区域，返回被删除的区域。
        - 这些函数面向多序列对齐后处理；BISER 的 `ltrim`/`rtrim` 是从链两端向内扫描、
          找累积比对得分最大的边界，与这些通用 trim 不同，需要自行实现。
    - `src/libs/alignment/slice.rs`:
        - `slice_block(block, name, set, writer) -> Result<()>`: 按参考物种的 chromosome runlist
          对一个 `FasBlock` 做切片，每个子区间按 `>range\nseq\n` 输出每个物种。内部使用
          `chr_to_align` / `align_to_chr` 做带 gap 的坐标转换。
        - 若 BISER 流程需要按区域从 MSA 中提取子序列，可参考其坐标转换逻辑。
- **需要新增的实现**
    - `src/libs/sd/anchor.rs`: 在 putative SD 的两个 mate 间生成 10-mer exact-match anchors， 包含
      `slide[d]` 去重、向右延伸、过滤高频 k-mer、过滤 trivial self-overlap。
    - `src/libs/sd/refine.rs`: 基于 y 坐标离散化 + segment tree / Fenwick tree 的 event-driven PST
      chaining + sparse DP 精修、大 gap 处理、两端 score-based `ltrim`/`rtrim`、生成最终 CIGAR。
      小 gap（≤1000 bp）调用 `ScalarAlignmentEngine` 并把 path 转为 CIGAR；大 gap 按 BISER 策略
      （两端各比对 1000 bp，中间用 `I`/`D`）处理。

#### 6.3.4 SD clustering（`cluster.codon`）

- **BISER 实现**
    - 文件: `biser/codon/cluster.codon:53-165`
    - 核心: 对 hit 的四个端点做区间 coloring，等价于 union-find 找重叠 hit 的连通分量， 然后提取每个
      cluster 的序列 FASTA。
- **PGR 可复用组件**
    - `src/libs/paf/graph/dsu.rs`:
        - 已实现 union-by-rank + path compression 的 `Dsu`，但它是 `pub(super)`，仅在 `paf::graph`
          模块内部可见，**不能直接作为公共 API 使用**。
        - 根据项目约束，纯数据结构应放在 `src/libs/ds/`。建议将 `Dsu` 迁移到 `src/libs/ds/dsu.rs`
          并公开为 pub，原 `paf/graph/dsu.rs` 通过 `pub use` 保持 API 兼容。SD 聚类直接复用
          `src/libs/ds/dsu.rs::Dsu`。
    - `src/libs/ds/dupe_tree.rs`:
        - `DupeTree::add(start, end)` + `build()` + `count_over(start, end, threshold)`: 0-based
          区间深度树，可用于统计 hit 端点或 cluster 区域的覆盖深度，识别高重复 hotspot。
    - `src/libs/ds/bitmap.rs`:
        - `BitMap::set_range` / `is_fully_set`: 0-based 位图，可用于标记已被 cluster 覆盖的碱基，
          避免重复提取。
    - `src/libs/fmt/fa.rs`:
        - `reader()` / `new_record()` / `writer()` / `writer_with_wrap()`: 读取基因组并构造 cluster
          FASTA。
    - `src/libs/fmt/fas.rs`:
        - block FA 读写与 `FasBlock` 数据结构。若 SD cluster 阶段需要输出多序列比对块（类似
          MAF/block FA），可参考该模块，但它紧密围绕 block FA 格式设计，不是通用 MSA 容器。
    - `src/libs/nt.rs`:
        - `rev_comp(seq)`: 生成反向互补序列，返回迭代器。
    - `src/libs/loc.rs`:
        - `create_loc(infile, locfile, is_bgzf) -> Result<()>`: 为 plain 或 BGZF FASTA 创建 `.loc`
          索引。普通代码更常用 `open_indexed` 自动创建。
        - `open_indexed(infile, force_update) -> Result<(Input, IndexMap<name, (offset, size)>)>`:
          打开带 `.loc` 索引的 FASTA （plain 或 BGZF 均可），不存在时自动创建索引。内部通过
          `is_bgzf` 判断压缩类型。
        - `open_input(infile, is_bgzf) -> Result<Input>`: 打开 FASTA 为 `Input::File` 或
          `Input::Bgzf`。
        - `fetch_record(reader, loc_of, name) -> Result<fasta::Record>`: 按名字读取完整 record。
        - `fetch_range_seq(reader, loc_of, rg) -> Result<String>`: 按 `intspan::Range`（1-based
          inclusive，支持 `chr(-):start-end` 链向）提取子序列。
        - `slice_record(record, rg) -> Result<fasta::record::Sequence>`: 从已加载 record 中按
          1-based Range 切片，负链会返回 reverse complement。
        - `get_seq_loc(file, range) -> Result<String>`: 便捷函数，对无效 range 或找不到的
          chromosome 返回空字符串而非报错；测试/脚本中可用，但生产代码建议用 `open_indexed` +
          `fetch_range_seq` 以明确错误处理。
    - `src/libs/io.rs`:
        - `SequenceReader` trait: `read_sequence(name, start, end) -> Result<String>`，定义 0-based
          half-open 的随机访问接口。`TwoBitFile` 实现该 trait，`chain_blocks` 的 `ScoreContext`
          也依赖它，因此 SD 流程中可用统一接口切换 2bit / indexed FASTA。
    - `src/libs/fmt/twobit.rs`:
        - `TwoBitFile::read_sequence(name, start, end, no_mask)`: 0-based half-open 随机访问 2bit
          序列，适合区间提取；`no_mask=true` 时返回 uppercase，否则保留 mask（N-blocks 变 N，
          soft-mask 变 lowercase）。
        - `TwoBitFile` 实现 `SequenceReader` trait，可直接传给 `ScoreContext` 等需要序列读取的接口。
    - `src/libs/paf/fasta.rs`:
        - `load_fasta_tsv(path) -> Result<IndexMap<name, path>>`: 读取 `name\tbgzf_fasta_path`
          格式的 TSV，用于将 PAF/SD 中的序列名映射到 BGZF FASTA 文件路径。
        - `prepare_store(tsv_path, idx) -> Result<FastaStore>`: 加载 TSV、校验覆盖所有 `idx.names`、
          构造 `FastaStore` 的一站式函数；SD 流程若用 TSV 管理多基因组输入，可直接复用。
        - `load_all_seqs(tsv_path) -> Result<HashMap<name, seq>>`: 一次性加载 TSV 中所有序列到内存，
          适合 cluster/decompose 阶段加载单个 cluster 的小规模 FASTA。
        - `FastaStore::new(seq_to_file)` + `fetch_range(name, start, end)` + `fetch_full(name)`:
          管理多个**BGZF FASTA**文件，带 `.loc` 索引与 LRU 缓存，适合多基因组
          cross_search/cross_align 时批量提取 mate 序列。
        - **限制**: `FastaStore::new` 内部使用 `noodles_bgzf::io::indexed_reader`，因此输入必须是
          BGZF 压缩的 FASTA。普通 gzip 或未压缩 FASTA 应先用 `loc`（plain/BGZF 通用）或 `twobit`
          处理。
        - **注意**: `FastaStore` **没有**实现 `SequenceReader` trait，它提供自己的 `fetch_range`
          / `fetch_full` API。若函数签名要求 `&mut dyn SequenceReader`（如 `chain_blocks` 的
          `ScoreContext`），应使用 `TwoBitFile` 而非 `FastaStore`。
- **需要新增的实现**
    - `src/libs/sd/cluster.rs`: 将 hit 端点排序，用 `Dsu` 合并重叠端点，输出每个 cluster 的 FASTA
      （序列名采用 `species#chrom+/-#start#end` 格式）。

#### 6.3.5 SD decomposition（`decompose.codon`）

- **BISER 实现**
    - 文件: `biser/codon/decompose.codon:52-277`
    - 核心: 对 cluster FASTA 建 10-mer 完整索引，再用 plane-sweep + mappings 输出 elementary SD
      集合。
- **PGR 可复用组件**
    - `src/libs/sd/kmer_index.rs`: 复用 exact 10-mer 索引（调整 `k` 参数即可）。
    - `src/libs/sd/plane_sweep.rs`: 复用链表扫描框架，但 decompose 需要支持多拷贝 mappings， 因此
      `PlaneSweepState` 需要泛化或单独实现一个 `MultiCopyPlaneSweep`。
    - `src/libs/ds/bitmap.rs`: 用于标记已输出的 `visited` 区域。
- **需要新增的实现**
    - `src/libs/sd/decompose.rs`: 在 cluster FASTA 上调用 10-mer 索引 + 多拷贝 plane-sweep + merge，
      输出 `.elem` 格式 BED。

#### 6.3.6 Core duplicon identification（`cover.py`）

- **BISER 实现**
    - 文件: `biser/cover.py:1-102`
    - 核心: 用 `ncls` 建立 elementary SD → 覆盖的 SD 列表映射，再用贪心 set cover 找出 能覆盖所有
      SD 的最小 elementary SD 集合，标记为 `CORE`。
- **PGR 可复用组件**
    - `coitrees`（已通过 `src/libs/paf/index/builder.rs` 使用）:
        - 直接建立 `Interval<ElemMetadata>` 区间树，查询重叠 interval，替代 `ncls`。
        - 比 `PafIndex` 更轻量；`PafIndex` 的 `query()` 与 `query_transitive_bfs()` 展示了
          “interval tree + CIGAR 投影”模式，可作为复杂场景参考。
    - `src/libs/paf/index/query.rs`:
        - `project(ts, te, metadata, cigar) -> Option<(qs, qe, ts, te)>` 实现了 target 子区间到
          query 坐标的投影，输入 `ts/te` 为 0-based half-open，返回 query 区间也是 0-based
          half-open。在将 elementary SD 区间映射到 SD 覆盖时思路类似，但注意它处理的是带 CIGAR
          的对齐投影，不是简单的坐标偏移。
    - `src/libs/paf/index/bfs.rs`:
        - `PafIndex::query_transitive_bfs(...)`: 从 seed target 区间出发做双向 BFS，遍历 alignment
          graph。
        - `PafIndex::merge_results(results, max_gap, fasta_store)`: 按
          `(query_id, target_id, strand)` 分组，合并 query 区间间隔不超过 `max_gap` 的相邻结果；
          若提供 `FastaStore`，合并不同 record 时会通过 FASTA 重新计算 CIGAR。
        - SD 流程中若需把多个重叠/相邻 hit 合并成连续区间（例如 cluster 后合并同一 chain 上的
          fragments），可参考 `merge_results` 的分组-排序-合并逻辑，但它依赖 PAF/CIGAR 语义，
          需要适配为 SD 的 interval merge。
- **需要新增的实现**
    - `src/libs/sd/set_cover.rs`: 用二叉堆实现 `greedy_set_cover`，输入为
      `elementary_id -> Vec<sd_id>`，输出被选中的 core elementary IDs。
    - `src/libs/sd/cover.rs`: 组装 interval overlap + set cover 流程，在 `.elem` 文件中追加 `CORE`
      标记。

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
    - `src/libs/alignment/coords.rs`: 这里的 `align_to_chr()` / `chr_to_align()` 是通过
      **带 gap 的对齐列**（用 `IntSpan` 记录非 gap 位置）在对齐坐标与基因组坐标之间转换，
      适用于 MSA/POA 输出。它**不适用于** BISER 的 hard-masked ↔ original 坐标映射，
      因为后者只是简单的“删除 lowercase 碱基”后产生的坐标偏移，不存在对齐 gap。
    - `src/libs/fmt/fa.rs::mask_sequence`: 该函数保留序列长度，只替换字符，不能用于 translate
      阶段的坐标映射。
- **需要新增的实现**
    - `src/libs/sd/translate.rs`: 在 `mask` 阶段记录 uppercase run 列表
      `[(orig_start, orig_end, masked_start, masked_end)]`（0-based half-open），通过二分查找实现
      hard-masked ↔ original 坐标双向映射，并同步改写 CIGAR。
    - 具体映射规则（需对照 BISER 源码精确实现）：
        - hard-masked 坐标下的 match/mismatch 区间，映射回 original 坐标时对应 uppercase 区段，
          CIGAR 保持 `=`/`X`/`M`。
        - 当 hard-masked 的 gap（`I`/`D`）跨越 original 中的 lowercase/masked 区域时，需要把该段替换为
          `S`（soft-mask）或 `N`（hard-mask），以反映原始基因组中这些碱基并非真正参与比对。
        - BISER 内部将 CIGAR 操作重新归类为 `M`/`S`/`N`；translate 阶段应产生与 BISER 输出语义一致的
          CIGAR。

### 6.4 其他值得复用的工具

除了按算法阶段映射的组件外，以下通用工具在 SD 流程中也 likely 有用：

- `src/libs/par.rs`
    - 并行 pipeline 原语，无 clap 依赖：
        - `spawn_writer_and_pool(outfile, num_threads) -> Result<(Sender<String>, JoinHandle<()>)>`：
          创建 writer thread + 按 `num_threads` 配置全局 rayon pool。
        - `resolve_paths(infile, is_list) -> Result<Vec<String>>`：`is_list=true`
          时按行读取路径列表，否则返回单元素 vec。
        - `load_entries<E, F>(paths, load_fn) -> Result<Vec<E>>`：顺序加载每个路径的条目。
        - `load_two_sets<E, F>(infiles, is_list, load_fn) -> Result<(Vec<E>, Vec<E>)>`：单输入时返回
          `(clone, original)` 用于自比较；双输入时分别加载。
        - `par_run_pairs(entries1, entries2, sender, pair_fn)`：rayon 并行遍历笛卡尔积，每 1000
          条结果通过 `sender` 批量发送。
    - SD pipeline 的 search/align 阶段若需要“生产者-消费者”式并行输出，可直接复用。
- `src/libs/io.rs`
    - `reader(input) -> Result<Box<dyn BufRead>>`: 通用缓冲读，支持 `stdin`、普通文件、`.gz`。
    - `writer(output) -> Result<PgrWriter>`: 通用缓冲写，支持 `stdout` 与普通文件；`PgrWriter` 在
      drop 时会尝试 flush，失败则向 stderr 输出警告。
    - `read_lines(input) -> Result<Vec<String>>`: 一次性读取所有行。
    - `read_names<T: FromIterator<String>>(path) -> Result<T>`: 读取第一列（按空白分割），可收集为
      `Vec<String>` 或 `HashSet<String>`。
    - `read_sizes<T: FromStr>(path) -> Result<BTreeMap<String, T>>`: 读取 `name\tsize`，
      支持任意数值类型。
    - `is_bgzf(path) -> bool`: 通过读取文件头判断是否为 BGZF 格式，`FastaStore` 与
      `loc::open_indexed` 内部均用此决定打开方式。
    - `read_runlist(path)`: 安全读取 runlist JSON 并转为 `BTreeMap<String, IntSpan>`，避免原
      `intspan` API 在错误输入上 panic。
    - `get_basename(path) -> Option<String>`: 提取文件基本名（去掉路径与扩展名），`lastz.rs`
      与多个命令用它生成输出文件名。
    - `SequenceReader` trait: 统一 0-based half-open 随机访问接口，`TwoBitFile` 已实现。
- `src/libs/ds/bitmap.rs::BitMap`
    - 固定大小的 0-based 位图，支持 `set_range(start, len)` 和 `is_fully_set(start, len)`。
    - 用途: 标记 plane-sweep 或 decomposition 中已访问/已输出的基因组位置；避免重复命中。
- `src/libs/ds/dupe_tree.rs::DupeTree`
    - 一维 0-based 区间深度树，支持 `add/subtract` 后 `build()`，再
      `count_over(start, end, threshold)`。
    - 用途: 统计 SD hit 端点或 cluster 区域在基因组上的覆盖深度，识别高重复 hotspot。
- `src/libs/ds/top_k_purity.rs::TopKPurity`
    - 跟踪离散类别计数，计算 top-K 类别占总观测的比例，并在比例过高时返回 penalty factor。
    - 用途: 在扫描 k-mer 或序列窗口时检测类别分布过于集中的低复杂度区域（不限于 AT 富集），
      作为额外过滤条件。
- `src/libs/fasta/stat.rs::count_bases`
    - 统计序列中 A/C/G/T/N 数量（IUPAC ambiguous codes 计为 N，其他非标准字符不计入长度）。
    - 用途: 快速评估 hard-masked 后有效碱基比例，或过滤 N 含量过高的 putative SD。
- `src/libs/fmt/twobit.rs::Block`
    - 2bit 内部使用的 0-based half-open mask block 类型。若 SD 流程需要记录 hard-masked 区间或
      N-block，可参考其区间重叠查询实现 `Blocks::overlaps`。
- `src/libs/paf/fasta.rs::FastaStore`
    - 多**BGZF FASTA**管理器，支持 `fetch_range(name, start, end)`（0-based half-open，参数为
      `i32`）与 `fetch_full(name)`，带 `.loc` 索引与 LRU 缓存。
    - 用途: 多基因组 cross_search/cross_align 时批量、高效地提取 mate 序列。注意输入必须是 BGZF
      压缩 FASTA。
- `src/libs/paf/persist.rs`
    - `PafIndex::save(path)` / `PafIndex::load(path)`: 将 interval tree + CIGAR 索引持久化为
      `.paf.idx`。
    - 用途: 若 SD 流程需要将 k-mer/anchor 索引或 hit 索引缓存到磁盘，可参考其 bincode + version +
      magic 的序列化模式。
- `src/libs/paf/record.rs` + `src/libs/paf/parser.rs`
    - `PafRecord { query_name, query_length, query_start, query_end, strand, target_name, target_length, target_start, target_end, matches, block_length, mapq, tags }`:
      `query_start/target_start` 为 0-based inclusive，`query_end/target_end` 为 0-based
      exclusive，符合 PAF 规范；`tags` 为可选 SAM 风格标签字符串列表。
    - `parse_paf(reader) -> Result<Vec<PafRecord>>`: 读取完整 PAF 文件，跳过空行与 `#` 注释行。
    - `parse_paf_line(line) -> Result<PafRecord>`: 解析单行，字段数不足 12 或 strand 非法时返回错误。
    - `write_paf_record(writer, rec)`: 写出完整 PAF 记录（含 tags）。
    - 用途: 若 SD 流程需要以 PAF 作为中间格式（例如把 hit 的 CIGAR 存为 PAF 后再投影/合并），
      可直接复用。
- `src/libs/fasta/stat.rs`
    - `count_bases(seq: &[u8]) -> (usize, [usize; 5])`: 统计序列中 A/C/G/T/N 的数量，IUPAC ambiguity
      按 N 计数，gap 等非标准字符不计入 `len`。可用于 SD 流程中快速计算 GC 含量或 N 含量。
    - `calc_n50_stats(lens, opt_nx, opt_genome) -> N50Stats`: 从序列长度列表计算 N50/Nx/E-size
      等组装统计。
    - `transpose<T>(v: Vec<Vec<T>>) -> Vec<Vec<T>>>`: 矩阵转置，输出工具函数。
- `src/libs/fasta/filter.rs`
    - `pass_filters(seq, minsize, maxsize, maxn, is_uniq, seen, name) -> bool`: 按长度、N 数、
      名称唯一性过滤 FASTA record。`NO_LIMIT`（`usize::MAX`）表示不限制。
    - `format_sequence(seq, is_dash, is_iupac, is_upper) -> String`: 可选删除 `'-'`、把 IUPAC
      折叠为 `N`、转大写。
    - 可用于 SD 流程的输入预处理（例如过滤太短/太 N 多的序列）。
- `src/libs/fasta/chunk.rs`
    - `SizeChunker`: 按累计序列大小或每 2 条记录切换输出文件的状态机；`max_files_exceeded()`
      判断是否达到最大文件数。
    - 主要用于 `pgr fa split`，与 BISER 关系不大。
- `src/libs/fasta/dedup.rs`
    - `record_signature(name, desc, seq, opts) -> Result<u64>`: 基于 rapidhash 计算 FASTA record
      的签名，支持按 name/description/sequence 去重，`is_both` 模式同时考虑正链与反向互补。
    - 若 SD 流程需要去除重复 contig/record，可直接复用。
- `src/libs/ds/top_k_purity.rs`
    - `TopKPurity::new(num_classes, k, ok_ratio)`: 跟踪类别计数，检测前 K 类是否超过可接受比例。
    - `penalty_factor()`: 若分布过于集中，返回 `Some(1.01 - (observed - ok_ratio) / (1 - ok_ratio))`
      作为 score 惩罚因子。
    - 被 `chain/anti_repeat.rs` 用于低复杂度过滤；SD 流程若需过滤高重复 k-mer 或低复杂度 hit
      可参考。
- `src/libs/fmt/psl.rs`
    - `Psl` 结构体：UCSC PSL 格式的内存表示，字段包括 `match_count`、`mismatch_count`、
      `rep_match`、`n_count`、`q_num_insert`、`q_base_insert`、`t_num_insert`、`t_base_insert`、
      `strand`、`q_name`、`q_size`、`q_start`、`q_end`、`t_name`、`t_size`、`t_start`、`t_end`、
      `block_count`、`block_sizes`、`q_starts`、`t_starts`。
    - `Psl::from_align(...) -> Option<Psl>`: 从两条等长对齐字符串构建 PSL 记录，会跳过两端纯 indel
      列。
    - `Psl::from_str(line)`: 从 PSL 行解析；支持 `"+"`、`"-"`、`"++"`、`"+-"` 等 strand。
    - 若 SD 流程需要与 PSL 格式互操作（例如把 POA alignment 导出为 PSL 再转 chain），可参考。
- `src/libs/lastz.rs`
    - 提供 UCSC lastz 预设（`PRESETS`、`find_preset`、`run_lastz`）与评分矩阵。可作为
      `ScalarAlignmentEngine` 性能不足时的外部局部比对 fallback。
    - 内置矩阵：`MATRIX_DEFAULT`（HoxD55， Human/Mouse/Macaque/Cow）、
      `MATRIX_DISTANT`（Human/Zebrafish/Opossum）、`MATRIX_SIMILAR`（Human/Chimp）、
      `MATRIX_SIMILAR2`（Human/Primate，更敏感）。
    - 内置 preset `set01`–`set07` 的参数（如 `O=400 E=30 K=3000 L=2200` 等）直接来自 UCSC pipeline；
      `run_lastz` 会处理 target/query 笛卡尔积、文件名去重与 `--self` 模式。
- `src/libs/chain/stitch.rs`
    - 按 chain ID 合并 fragments（`chainStitchId` 语义），要求 fragments 之间不重叠。SD 流程若输出带
      ID 的 chain fragments，可参考；但它不是通用的“合并相邻 chain”工具。
- `src/libs/chain/anti_repeat.rs`
    - `check_chain(chain, t_2bit, q_2bit, min_score) -> bool`: 对 UCSC chain 做 degeneracy 与
      repeat 过滤。
    - `check_degeneracy`: 用 `TopKPurity(4, 2, 0.80)` 检测低复杂度（1–2 个碱基占比过高），
      若占比过高则按 penalty factor 降低 chain score。
    - `check_repeat`: 检测 chain 中 soft-masked（lowercase）碱基比例，若过高则降低 score。
    - 内部使用与 `NT_VAL` 不同的 2-bit 编码（T=0, C=1, A=2, G=3），因负链 complement 用
      `(v + 2) % 4`；不能直接用 `crate::libs::nt::NT_VAL`。
    - SD 流程若需过滤低复杂度或重复区域 hits，可参考其思路，但 `check_chain` 强依赖 `TwoBitFile` 与
      chain 结构，需要适配。
- `src/libs/chain/pre_net.rs`
    - `pre_net(reader, writer, t_hash, q_hash, opts)`: 实现 UCSC `chainPreNet` 语义，按 score
      降序遍历 chain，保留 target/query 尚未被完全覆盖的 chain，并用 `BitMap` 标记已用区间（可带
      pad）。
    - 输入 chain 必须已按 score 降序排序，否则会报错。
    - `is_haplotype(name)`: 判断名称是否含 `_hap` 或 `_alt`。
    - 该逻辑与 BISER 的 coverage/decomposition 有相似之处，但针对 chain 格式；SD
      的“去重复覆盖”步骤可参考其实现。
- `src/libs/chain/psl_chain.rs`
    - `group_psl_blocks(reader, score_ctx)`: 读取 PSL 记录，按 `(target, query, strand)` 分组为
      `ChainableBlock`。
    - `chain_psl(reader, writer, gap_calc, min_score, score_ctx)`: 对每组调用 `chain_blocks` 生成
      chain，过滤低分链后输出。
    - 展示了 PSL → chain 的完整 pipeline，但使用 `chain_blocks` 的 `max(dq, dt)` gap 模型，不适用于
      BISER refine。
- `src/libs/alignment/trim.rs`
    - `trim_pure_dash(seqs)`：删除所有序列在该列均为 gap（`-`）的列（交集）。
    - `trim_head_tail(seqs)`：从两端删除纯 gap 列，直到遇到任一序列的非 gap 字符。
    - `trim_outgroup(seqs)`：要求至少 3 条序列且最后一条为 outgroup，删除 outgroup-only 的插入列。
    - `trim_complex_indel(seqs)`：在 `trim_outgroup` 后使用，识别并删除 ingroup 内部复杂 indel
      区域。
    - SD refine 阶段若需要对齐后修剪，可参考这些函数，但 BISER 的 score-based `ltrim`/`rtrim`
      仍需自行实现。
- `src/libs/alignment/variation.rs`
    - `get_subs` / `get_indels`: 从 MSA 中检测替换与 indel；`Substitution` / `Indel` 结构体记录位置、
      碱基、频率、模式。
    - 主要用于 MSA 变异检测，不是 BISER SD 检测核心，但可作为 MSA 后处理参考。
- `src/libs/fas_multiz/banded_align.rs`
    - `banded_align_refs(...)` 对两个 `FasBlock` 的 reference entry 做 banded DP + affine gap，
      返回两 reference 序列列与列之间的对齐映射 `(Vec<Option<usize>>, Vec<Option<usize>>)`。
    - DP 会同时考虑两个 block 中所有共有的 species，对每个 `match/mismatch` 位置累加所有物种对的
      `SubMatrix::get_score` 得分（gap 只算一次），因此它本质上是多序列 banded 对齐，不是单纯的
      pairwise。
    - 当前实现紧密绑定 `FasBlock` 输入，不能直接复用于 BISER 的 pairwise gap 精修。若后续需要
      banded pairwise align，可参考其索引函数 `idx(i, j)` 和 band 半径计算逻辑，提取为通用函数。
    - 该模块使用 `crate::libs::chain::GapCalc`，它是 `src/libs/ds/gap_calc.rs::GapCalc` 的
      re-export。
- `src/libs/fas_multiz/merge.rs` + `src/libs/fas_multiz/windows.rs`
    - `merge_blocks_with_dp` / `merge_two_blocks_with_dp`: 基于 `banded_align_refs` 对多个
      `FasBlock` 做 progressive 合并，`Core` 模式要求所有输入都包含该物种，`Union` 模式允许缺失。
    - `derive_windows_from_blocks`: 从 reference ranges 推导合并窗口，先按染色体合并重叠区间（扩展
      `radius`），再按 `min_width` 与覆盖输入数过滤。
    - 面向多物种 block FA 合并，与 BISER 单基因组 SD 检测无直接复用关系，但可作为“区间合并 +
      progressive 对齐”模式参考。

### 6.5 建议的模块与命令结构

#### 6.5.1 `src/libs/sd/` 目录（新增）

**第一阶段（LASTZ-based search，优先实现）**

- `search_lastz.rs`: 封装 lastz 自比对 → 格式转换 → 过滤 → putative hits 输出。
- `coverage.rs`（可选）: 滑动窗口覆盖度重复区检测，复用 `pgr-repeat.sh` 逻辑。
- `subtract_repeatmasker.rs`（可选）: 用 RepeatMasker 结果过滤候选重复区。
- `hit.rs`: SD hit 数据结构（坐标、species、strand、CIGAR、error rate）。

**第二阶段及以后（BISER 后续算法）**

- `anchor.rs`: 10-mer exact-match anchor 生成。
- `refine.rs`: 基于 y 坐标离散化 + segment tree / Fenwick tree 的 event-driven PST chaining +
  sparse DP 的比对精修；包含 `path_to_cigar` 辅助函数。
- `cluster.rs`: 重叠 hit 聚类并输出 cluster FASTA。
- `decompose.rs`: elementary SD 分解。
- `set_cover.rs`: 贪心 set cover。
- `cover.rs`: core duplicon 标记流程。
- `translate.rs`: hard-masked 与原基因组坐标互转。注意 `src/libs/translate.rs` 已存在（蛋白质翻译），
  新增 SD 的 `translate.rs` 位于 `src/libs/sd/translate.rs`，不会冲突。

**未来可选（原生 BISER search）**

- `kmer_index.rs`: exact 2-bit k-mer 索引 + winnowing + 频率过滤。
- `plane_sweep.rs`: plane-sweep 链表与 hit 输出。

#### 6.5.2 `src/cmd_pgr/sd/` 目录（新增）

- `mod.rs`: 子命令注册与分发。
- `mask.rs`: `pgr sd mask <genome.fa> -o <masked.fa>`。
- `search.rs`: `pgr sd search <genome.fa> -o <hits.psl> [--mode lastz|coverage|kmer]`，输出原始
  pairwise alignments（lastz LAV 经 `pgr lav to-psl`；FastGA PAF 需 PAF→PSL），供 `pgr sd align` 的
  chain/net 消费。默认 `lastz`，未来可扩展 `kmer` 原生模式。
- `align.rs`: `pgr sd align <genome.fa> <hits.psl> -o <hits.align.paf>`，封装 `pgr pl ucsc`（不加
  `--syn`）+ `pgr maf to-paf`，用 UCSC chain/net 做 chaining/refine。
- `cluster.rs`: `pgr sd cluster <genomes...> <hits.align.paf> -o <clusters.dir>`。
- `decompose.rs`: `pgr sd decompose <cluster.fa> -o <cluster.elem.bed>`。
- `cover.rs`: `pgr sd cover <hits.align.paf> <elems.txt> -o <elems.covered.txt>`。
- `translate.rs`: `pgr sd translate <hits.align.paf> <genomes...> -o <out.paf>`。
- `run.rs`: `pgr sd run <genomes...> -o <out.bed>`，按 BISER 顺序串接上述步骤。

### 6.6 分阶段实施计划

#### 第一阶段：LASTZ-based putative SD 检测（验证：输出格式正确、下游可消费）

1. 实现 `libs/sd/hit.rs`：统一的 SD hit 数据结构，支持从 PAF / chain 初始化。
2. 实现 `libs/sd/search_lastz.rs`：
    - 调用 `pgr::libs::lastz::run_lastz` 做 `--self` 自比对，或按用户指定做 target/query 比对；
    - 用 `pgr::libs::fmt::lav` 解析 LAV，或先 `pgr lav to-psl` 再读 PSL；
    - 输出 PSL（lastz LAV 经 `pgr lav to-psl`；PAF 输入则需 PAF→PSL），供 `pgr sd align` 的
      chain/net 消费，不做 chaining；
    - 过滤：长度 < 1 kbp、error rate > 0.1（即同一性 < 90%，采用 T2T-CHM13 SD 标准，见 4.2.1；用
      `block_identity` 计算）、低复杂度（`TopKPurity`）。
3. 实现 `cmd_pgr/sd/search.rs`：
    - CLI：`pgr sd search <genome.fa> -o <hits.psl> [--mode lastz]`；
    - 默认输出 PSL（lastz LAV 经 `pgr lav to-psl`；FastGA PAF 需 PAF→PSL），坐标 0-based half-open；
    - 序列名按 `species#chr` 编码（为 `cluster` 阶段输出 FASTA 头做准备）。
4. 验证:
    - `pgr sd search` 能在人类 chr21 上成功运行并输出合法 PSL；
    - 输出 PSL 可被 `pgr sd align`（封装 `pgr pl ucsc`）消费；
    - hit 数量级与 BISER `search` 同数量级（不必 bit-exact，但不应差一个数量级）。

> **已实现（2026-08-02）**：`pgr sd search <genome.fa> -o hits.psl` 落地
> （`src/libs/sd/search_lastz.rs` + `src/cmd_pgr/sd/search.rs`）。流程：lastz --self
> （.gz 输入自动解压）→ LAV → PSL → 按 block_len ≥ `--min-len`（默认 1000）且
> block_identity ≥ `--min-identity`（默认 0.90，`(matches+rep)/block_len`，含 insert 碱基）
> 过滤。MG1655 实测：81 秒，264 条 putative hits。下游链路已验证：
> `pgr pl chainnet`（**非 --syn**，原生实现替代 `pgr pl ucsc`）→ `pgr maf to-paf`
> 产出 90 条 PAF，可直接接 cluster/decompose。`pgr sd align`（封装上述 chainnet + to-paf）
> 属第二阶段，尚未实现。

**可选（若选择形态 B：pgr-repeat.sh 覆盖度路线）**

1. 实现 `libs/sd/coverage.rs`：封装窗口化、lastz 回贴、lift、覆盖度计算，输出候选重复区 BED。
2. 实现 `libs/sd/subtract_repeatmasker.rs`：读取 RepeatMasker `.out`/BED，做区间差集。
3. 在 `cmd_pgr/sd/search.rs` 中增加 `--mode coverage`：
    - 先调用 `coverage.rs` 得到候选区；
    - 减 TE；
    - 提取候选区序列，调用 lastz/minimap2 做 all-vs-all 自比对；
    - 输出与 `--mode lastz` 格式一致的 PSL，供 `pgr sd align` 的 chain/net 消费（不直接产出 PAF）。

#### 第二阶段：比对精修（验证：hit 的 CIGAR 与 error rate 与 BISER align 一致）

> **注**：外部比对路线（lastz/FastGA）的精修由 UCSC chain/net（`pgr pl ucsc` 不加 `--syn`）+
> `pgr maf to-paf` 完成，不实现下列 PST refine；下列 PST refine 仅用于原生 BISER 路线，可延后。
> 外部路线只需在 `pgr sd align` 中封装 `pgr pl ucsc`。

1. 实现 `libs/sd/anchor.rs`：10-mer anchor 生成。
2. 实现 `libs/sd/refine.rs`：
    - 用 y 坐标离散化 + segment tree / Fenwick tree 实现 BISER 风格的 event-driven PST chaining：
      扫描线按 x 坐标推进，右端点事件时以 `dp[i] - gap_to_end` 单点更新 y 位置；左端点事件时先 在 y
      区间 `[ay - MAX_CHAIN_GAP, ay - 1]` 查询最大前驱得分，再计算 gap 惩罚 `dx + dy` 更新 `dp[i]`，
      最后把滑出 `MAX_CHAIN_GAP` 窗口的旧锚点对应 y 位置重置为 `-INF`。`chain_blocks` 的 `GapCalc`
      取 `max(dq, dt)`，不能 bit-exact 匹配 BISER。
    - 小 gap（≤1000 bp）用 `ScalarAlignmentEngine`（把一条 mate 构建成线性 POA graph）， 并自己实现
      `path_to_cigar` 从 `Alignment.path` 生成 `=`/`X`/`I`/`D` CIGAR。
    - 大 gap 按 BISER 策略处理（两端各比对 1000 bp，中间用 `I`/`D`）。
    - 用 `paf/cigar.rs::format_cigar` 输出 CIGAR，用 `block_identity` 计算 BISER 风格的 error rate
      （`1 - block_identity`）；不要使用 `gap_compressed_identity`，因为它会低估 indel 错误。
3. 实现 `cmd_pgr/sd/align.rs`。
4. 验证: 对同一组 putative hits，PGR 与 BISER 输出的 alignment span、CIGAR、error rate 差异 < 1%。

#### 第三阶段：聚类与分解（验证：elementary SD 集合与 BISER 一致）

1. 前置：将 `src/libs/paf/graph/dsu.rs` 的 `Dsu` 迁移到 `src/libs/ds/dsu.rs` 并公开为 pub，
   原文件通过 `pub use` 保持兼容。
2. 实现 `libs/sd/cluster.rs`：用 `ds::Dsu` 聚类并输出 cluster FASTA（复用 `fa::new_record`、
   `nt::rev_comp`、`loc::fetch_range_seq`、`twobit::TwoBitFile` 或 `paf/fasta.rs::FastaStore`
   提取子序列）。
3. 实现 `libs/sd/decompose.rs`：多拷贝 plane-sweep。
4. 实现 `cmd_pgr/sd/cluster.rs` 与 `cmd_pgr/sd/decompose.rs`。
5. 验证: cluster 数量与覆盖范围、elementary SD 数量与 `.elem` 内容一致。

#### 第四阶段：core duplicon 与坐标转换（验证：CORE 标记与 translate 后坐标一致）

1. 实现 `libs/sd/set_cover.rs` 与 `libs/sd/cover.rs`（用 `coitrees` 做区间重叠查询）。
2. 实现 `libs/sd/translate.rs`（基于 `mask` 阶段记录的 uppercase run 列表做二分查找映射， 不要误用
   `alignment/coords.rs` 的 gapped-alignment 坐标函数）。
3. 实现对应子命令。
4. 验证: CORE 标记集合与 BISER 一致；translate 后的坐标与 CIGAR 可通过一致性检查。

#### 第五阶段：pipeline 与跨基因组（验证：多基因组输入与 BISER 最终输出一致）

1. 实现 `cmd_pgr/sd/run.rs` 串联所有步骤。
2. 实现 `cross_search` / `cross_align` 等价逻辑（复用 search/align，只是 query 为另一基因组）。
3. 验证: 多基因组场景下最终 `out.bed` 与 `.elem.txt` 与 BISER 等价。

### 6.7 风险与注意事项

- **比对器依赖（基本可解决，但性能有差异）**: BISER 的 alignment 精修依赖 `bio.seq.align()`
  这个内置 SIMD aligner。PGR 的 `src/libs/poa/align.rs::ScalarAlignmentEngine` 提供 Global、Local、
  SemiGlobal 三种模式以及自定义 match/mismatch/gap_open/gap_extend 参数。注意 `src/libs/poa/mod.rs`
  只导出 `AlignmentParams`、`AlignmentType`、`Poa`；`ScalarAlignmentEngine` 和 `PoaGraph`
  需分别通过 `poa::align::ScalarAlignmentEngine` 和 `poa::graph::PoaGraph` 使用。把其中一条 mate
  构建成线性 POA graph 后即可做 pairwise alignment；小 gap（≤1000 bp）可直接复用 `libs/poa`。
  但它是标量 O(nm) 实现，无 SIMD/banded；在人类全基因组尺度批量调用可能显著慢于 BISER。
  后续若性能成为瓶颈，再评估 SIMD aligner（如 `parasail-rs`）或提取 `fas_multiz/banded_align.rs`
  的 banded DP 为通用函数。
- **`chain_blocks` 不能直接复用于 BISER chaining**: `src/libs/chain/connect.rs::chain_blocks` 的
  `cost_func` 是内部硬编码的，使用 `GapCalc::calc(dq, dt)`。`GapCalc` 在 `dq > 0 && dt > 0` 时取
  `max(dq, dt)` 查表，而 BISER 要求 `dx + dy`。因此不能通过传入参数让 `chain_blocks` bit-exact
  匹配 BISER；必须用 y 坐标离散化 + segment tree / Fenwick tree 自行实现 BISER 的 event-driven PST
  chaining。现有 `src/libs/ds/kdtree.rs::KdTree` 不支持 deactivate，不能直接用。
- **BISER refine DP 的 gap 模型需自行实现**: BISER 的 sparse DP 使用特殊公式
  `MISMATCH * mi + GAPOPEN + GAP * (ma - mi)`，不能简单套用 `GapCalc::affine`。
- **`Dsu` 需要迁移到 `src/libs/ds/`**: `src/libs/paf/graph/dsu.rs::Dsu` 是 `pub(super)`，
  不能作为公共 API。根据项目约束，纯数据结构应放 `src/libs/ds/`。建议迁移到 `src/libs/ds/dsu.rs`
  并公开，原文件通过 `pub use` 保持兼容。
- **`hash.rs` 不是 exact k-mer 索引（当前阶段不实现）**: `src/libs/hash.rs` 提供基于哈希的
  canonical minimizer 采样（`seq_sketch`、`JumpingMinimizer`），而 BISER 的 search/decompose
  依赖 exact 2-bit k-mer + winnowing。由于第一阶段改用 LASTZ-based search，k-mer 索引与
  plane-sweep 暂时不需要实现；未来若要替换为原生 BISER search，再新增 `src/libs/sd/kmer_index.rs`
  与 `src/libs/sd/plane_sweep.rs`，`hash.rs` 仅可作为 sketch 验证或后续扩展使用。采样底座上，
  `src/libs/syncmer.rs` 的 DNA 路径已基于 2-bit canonical rolling hash（复用 `nt::NT_VAL`），比
  `hash.rs` 更贴近 BISER 的 2-bit k-mer。
- **坐标系统不一致**: PGR 内部不同模块使用不同坐标约定：
    - `chain`、`paf`、`twobit`、`FastaStore`、`BitMap`、`DupeTree`、`io::SequenceReader` 使用
      0-based half-open。
    - `loc.rs` 的 `intspan::Range` 和 `slice_record` 使用 1-based inclusive，支持链向。
    - `fa::mask_sequence` 的 `IntSpan` 是 1-based inclusive；`fa::windows` 输出的坐标也是 1-based
      inclusive。
    - BISER 内部大量 0-based half-open。迁移时应在 SD 模块内部统一使用 0-based half-open，仅在调用
      `loc` 或输出时转换。
- **hard-masked ↔ original 坐标映射不要误用 `alignment/coords.rs`**: `alignment/coords.rs` 是处理
  带 gap 对齐列的坐标转换，而 BISER 的 translate 只是“删除 lowercase 碱基”后的简单偏移映射。应在
  `mask` 阶段记录 uppercase run 边界，在 `translate` 中做二分查找。
- **大染色体内存问题**: `fa::reader` 通过 `noodles_fasta` 顺序读取时，会把整条 record
  （完整染色体序列）读入内存（`Record.sequence()` 返回完整序列）。人类尺度染色体（~250 Mbp）尚可，
  但更大基因组或内存受限时，建索引阶段应改用 `src/libs/fmt/twobit.rs::TwoBitFile` 顺序扫描，
  或按 `MAX_CHROMOSOME_SIZE` 切片读取。区间提取 mate 序列时，优先使用 `loc` 或 `twobit`，
  避免加载完整染色体。
- **反向互补**: BISER 在 chromosome 级别同时索引 forward 与 reverse complement， PGR 中可用
  `nt::rev_comp` 生成反向链序列，或按 BISER 方式在索引阶段同时扫描两条链。
- **LASTZ-based search 的灵敏度与计算代价**: 用 lastz 替代 BISER k-mer plane-sweep 可以快速拿到
  putative hits。PGR 采用 > 90% 同一性标准（见 4.2.1），lastz 默认 `set01`（`K=3000 L=2200`）
  正好匹配，无需为低同源性额外调参；如需略提灵敏度可适度降低 `K/L` 或使用 `--self` 模式。另外，
  全基因组 lastz 自比对 + 后续 all-vs-all 候选区比对的计算量通常大于 BISER 的 plane-sweep，
  人类尺度基因组需充分并行化。
- **性能**: BISER 的 Codon 实现经编译为原生代码，plane-sweep 是其性能关键。 当前阶段不实现
  plane-sweep，因此该风险暂时规避；未来若替换为原生 BISER search，Rust 实现需对 k-mer 索引做内存优化
  （如 `u64` key + `Vec<(u32, u32)>`）。

## 6.8 外部全基因组自比对替代路线：lastz --self 与 FastGA

> 本节从“复用 PGR 已有的外部比对基础设施”出发，讨论是否可以用 lastz --self 或 FastGA 替代
> BISER 内部的 `search` + `align` 阶段，直接把外部比对结果经 UCSC chain/net（`pgr pl ucsc`，
> SD 不加 `--syn`）精炼后转为 PAF（见 6.9），再接入 `cluster` / `decompose` / `translate`。
> **lastz 与 FastGA 内置的 chaining 均不采用，chaining 统一由 UCSC chain/net 承担**。所有分析均结合
> BISER 源码中的实际输入输出格式与 PGR 现有命令能力，并给出可落地的模块与命令设计。

> **注意**：在 6.1/6.3.2/6.6 的简化实施策略中，`lastz --self` 已被选为第一阶段 putative SD 检测的
> 主要实现路径。本节的技术细节（LAV/PSL/PAF 转换、字段映射、tag schema）因此成为当前迁移方案的核心
> 参考，而不再只是“替代路线”。

### 6.8.1 BISER search/align 的内部数据契约

要判断外部比对能否替代 BISER 的 `search` + `align`，首先必须明确这两个阶段的数据契约。

**`biser search` 的输出**

- 文件：`biser/codon/search.codon:61-86` 的 `save_sd()`。
- 输出时机：plane-sweep 链表中的 walker 满足 `count >= ceil(age * tau)` 且长度超过
  `QUERY_THRESHOLD=100` 与 `REF_THRESHOLD=500` 时被提升为最终 hit。
- 输出格式：7 列 BEDPE（`x_name\tx_start\tx_end\ty_name\ty_start\ty_end\tspecies1:species2`），
  **无 CIGAR、无链向列**。
- 坐标系：hard-masked 基因组上的 0-based half-open；同染色体时会通过 `pad_sd()` 做 `MAX_EXTEND=5000`
  填充，并保证两个 mate 不重叠（`biser/codon/search.codon:62-69`）。

**`biser align` 的输入与输出**

- 输入：即 `search` 输出的 7 列粗 BEDPE。
- 处理：
    1. 用 `align.generate_anchors()`（`biser/codon/align.codon:5-28`）在 mate 之间找 10-mer
       精确匹配锚点；
    2. 用 `chain.chain()`（`biser/codon/chain.codon:171-244`）做 PST chaining，得到粗略锚点链；
    3. 用 `align.refine()`（`biser/codon/align.codon:31-112`）做 sparse DP 精修边界并生成 CIGAR。
- 输出：14 列 `.align` BEDPE，格式由 `hit.Hit.__str__()`（`biser/codon/hit.codon:155-171`）定义。

**下游阶段只依赖 `.align` 文件**

- `biser cluster`（`biser/codon/cluster.codon:53-165`）读取 `.align` 后：
    - 用第 1–3 列与第 4–6 列构建 `SD` 的 `mate1` / `mate2` 区间；
    - 用第 7 列拆分 species（`l[6].split(':')`）；
    - 用第 10 列判断反向（`l[9] == '-'`）；
    - **不解析第 13 列 CIGAR**。
- `biser translate`（`biser/codon/mask.codon:32-154`）读取 `.align` 后：
    - 用第 1–6 列做 hard-masked → original 坐标映射；
    - **逐字符解析第 13 列 CIGAR**，把 `M`/`I`/`D` 映射为 original 坐标系下的 `M`/`S`/`N`；
    - 因此 CIGAR 的语义必须与 BISER 完全一致。

**结论**：外部比对只要按 6.9 的映射生成标准 PAF，就可以无缝替代 `search` + `align`，`cluster` /
`decompose` / `translate` 直接以 PAF 为输入。

### 6.8.2 为什么外部比对可以替代 BISER 的 search/align

BISER 内部 pipeline 在单基因组场景下的数据流如下（`biser/__main__.py:480-499`）：

1. `biser search`（`biser/codon/search.codon:265`）：在 hard-masked 基因组上建立 2-bit k-mer 索引，
   用 plane-sweep 找出 putative SD pairs，输出无 CIGAR 的粗 BEDPE。
2. `biser align`（`biser/codon/__init__.codon:75-104`）：读取粗 BEDPE，对每对 mate 用
   `align.generate_anchors` + `chain.chain` + `align.refine` 精修边界并生成 CIGAR，输出带完整 CIGAR
   的 14 列 BEDPE（`.align`）。
3. 后续 `cluster` / `decompose` / `translate` 只依赖 `.align` 文件的坐标、链向与 CIGAR。

因此，
**只要能把“全基因组自我比对”结果按 6.9 转换为标准 PAF，就可以跳过 BISER 自研的 search/align，直接复用下游阶段**。
PGR 已经具备成熟的外部比对调用与格式转换能力，这为该替代路线提供了工程基础。

### 6.8.3 BISER `.align` 文件格式（参考：6.9 的 PAF 映射依据）

`biser/codon/hit.codon:155-171` 定义了 `Hit.__str__` 输出的 14 列：

- 第 1–3 列：`x_name`, `x_start`, `x_end`（query mate，0-based half-open）
- 第 4–6 列：`y_name`, `y_start`, `y_end`（target mate，0-based half-open）
- 第 7 列：`species1:species2`
- 第 8 列：总错误率 `err() * 100`（`gap_err + mis_err`，`hit.codon:237-238`）
- 第 9–10 列：`x_strand`, `y_strand`（`+` 或 `-`；BISER 保证 `x` 永远为正链，`hit.codon:100`）
- 第 11 列：`max(x_end - x_start, y_end - y_start)`
- 第 12 列：`span()`（CIGAR 总长度）
- 第 13 列：`simple_cigar`（把连续的 `=`/`X` 合并为 `M`，`I`/`D` 保留，如 `100M50I200M`）
- 第 14 列：`X=...;ID=...`（ mismatch 错误率与 gap 错误率）

下游阶段对列的使用：

- `cluster`（`biser/codon/cluster.codon:53-84`）只读取第 1–7 列与第 10 列（`l[9] == '-'` 判断反向），
  不解析 CIGAR。
- `translate`（`biser/codon/mask.codon:69-147`）会读取第 13 列 CIGAR，并在 hard-masked
  坐标与原基因组坐标之间转换。
- 因此 PAF 转换器必须保证坐标、链向、CIGAR 三要素与 BISER 语义一致，错误率按 6.9.6 的 `er`/`xm`/
  `id` tag 计算（沿用 `gap_err + mis_err` 与 `block_identity` 思路，见 6.3.3 对 `paf/cigar.rs`
  的讨论）。

### 6.8.4 lastz --self 路线

PGR 已有 `pgr lav lastz` 命令（`src/cmd_pgr/lav/lastz.rs:86-194`），其 `RunLastzOptions.is_self`
会在 target 与 query 为同一文件时调用 lastz 的 `--self` 标志（`src/libs/lastz.rs:138-247`），
避免冗余计算。

**完整数据流**：

1. **Hard-mask 输入基因组**（与 BISER 保持一致）：

  ```bash
  pgr sd mask genome.fa -o genome.hard.fa
  ```
这一步对应 BISER `mask.codon:7-15`，只保留 uppercase A/C/G/T，删除 lowercase 与 N；PGR
   实现需同步记录 uppercase run 映射表供后续 `translate` 使用（见 6.3.7）。

2. **lastz 自比对**：

  ```bash
  pgr lav lastz genome.hard.fa genome.hard.fa --self --preset set01 -o lav_out/
  ```
输出为 LAV 文件（`src/libs/fmt/lav.rs`）。lastz 的 `--self` 会自动跳过 query 与 target
   完全相同且同方向的 trivial 比对，只保留有意义的同源区段。

3. **LAV → PSL**：

  ```bash
  pgr lav to-psl lav_out/[genome]vs[genome].lav -o self.psl
  ```
对应 `src/cmd_pgr/lav/to_psl.rs`，调用 `fmt::lav::lav_to_psl`（`src/libs/fmt/lav.rs:403-487`）。
   该函数解析 LAV 的 `s {}/h {}/a {}` stanzas，把 `l` 行（`l t_start q_start t_end q_end percent_id`）
   转换为 0-based half-open 的 `Block`，再汇总为 `Psl`。

4. **PSL → UCSC chain/net（必经，`pgr pl ucsc`，SD 不加 `--syn`）**：

  ```shell
  # target=query=自身基因组，做 self chain/net；默认（不带 --syn）即不做共线性筛选
  pgr pl ucsc genome.hard.fa genome.hard.fa self.psl -o ucsc_out/
  # 输出 ucsc_out/*.maf（chain/net 精炼后的 pairwise alignments）
  ```
对应 `src/cmd_pgr/pl/ucsc.rs`，内部串接 `axtChain`（`-psl`、`-linearGap=loose|medium`、
   `-minScore`）→ `chainAntiRepeat` → `chainMergeSort` → `chainPreNet` → `chainNet` + `netSyntenic`
   → `netChainSubset` + `chainStitchId` → `netSplit` → `netToAxt` + `axtSort` → `axtToMaf`。
   **chaining 统一走此流程，不用 lastz 内置 chaining，也不用 `pgr psl chain` 单独建链**。
   SD 场景必须使用默认（非 `--syn`）路径：`--syn` 会经 `netFilter -syn` 只保留共线性比对，
   丢掉伴随重排的 SD。`axtChain` 的 gap 模型为 `max(dq, dt)`（与 BISER 的 `dx + dy` 不同），
   这是本路线接受的已知差异。

5. **MAF → PAF（`pgr maf to-paf`）**：

    - `pgr pl ucsc` 输出的 `ucsc_out/*.maf` 直接喂给
      `pgr maf to-paf`（`src/libs/paf/maf_import.rs::maf_block_to_paf`），按 6.9 的映射输出标准 PAF
      （12 列 + `cg:Z:` + 推荐 tag）。
    - 坐标：MAF block 的坐标天然来自 hard-masked 基因组上的 chain/net，直接填入 PAF 的
      `query_start/end` 与 `target_start/end`。
    - CIGAR：`pgr maf to-paf` 从 MAF block 的逐列比对直接生成 `=`/`X`/`I`/`D` 写入 `cg:Z:`，
      无需回查 FASTA。
    - 错误率：用 `block_identity` 计算 `1 - identity` 作为总错误率写入 `er:f:`，并拆分 mismatch
      与 gap 比例写入 `xm:f:` / `id:f:`（见 6.9.6）。不要误用 `gap_compressed_identity`（6.3.3
      已说明其会低估 indel）。
6. **接入下游**：转换后的 PAF 可直接作为 `pgr sd cluster` 的输入，后续 `decompose` / `translate` 按
   6.9.3 消费 PAF。

### 6.8.5 FastGA / PAF 路线

FastGA 在 PGR 文档中被描述为 impg 支持的备选 aligner（`docs/paf.md:39`），输出 PAF 格式，磁盘占用低，
适合小基因组/细菌基因组的 all-vs-all 场景（`notes/ecoli-cohort.md`）。

**完整数据流**：

1. 同样需要先做 `pgr sd mask` 得到 hard-masked FASTA。

2. **FastGA 自比对**：

  ```bash
  FastGA -pafx genome.hard.fa genome.hard.fa > self.paf
  ```
`-pafx` 表示输出 PAF 并包含 CIGAR（`cg:Z:` tag）。

3. **PAF → PSL → UCSC chain/net**：与 lastz 路线一样，FastGA 内置的 chaining 不采用，chaining
   统一交由 UCSC chain/net。

    - 先把 FastGA PAF 转为 PSL（PAF 与 PSL 列语义可互转，供 `axtChain -psl` 消费；该 PAF→PSL
      转换需新增）；
    - 再走 `pgr pl ucsc genome.hard.fa genome.hard.fa self.psl -o ucsc_out/`（默认非 `--syn`），
      输出 MAF；
    - 经 `pgr maf to-paf` 转为 PAF，并补齐 6.9.6 的推荐 tag（`sp:Z:`、`er`/`xm`/`id`、`ms`/`sp2`）。
    - PAF 的 `target_start/query_start` 为 0-based inclusive，`target_end/query_end` 为 0-based
      exclusive（`src/libs/paf/record.rs`），与 BISER 的 0-based half-open 一致，坐标无需转换。
4. **接入下游**：与 lastz 路线相同，PAF 进入 `pgr sd cluster`。

相比 lastz，FastGA 在小基因组上自比对更快、磁盘占用更低；但两条路线在 chaining 之后完全合流
（都走 UCSC chain/net + `pgr maf to-paf`），不再有"FastGA 步骤更短"的优势。PGR 目前对 FastGA
没有原生命令封装，需要用户在 PGR 外部安装并调用。

### 6.8.6 两种路线与 BISER 原生的技术对比

- **灵敏度**
    - **BISER**：ordered Jaccard + plane-sweep 针对 SD error model 设计，
      默认允许总错误率高达 `--max-error=30`（即低至 70% 同一性），并区分 `MAX_ERROR` 与
      `MAX_EDIT_ERROR`（`biser/__main__.py:407-413`）。
    - **lastz/FastGA**：默认参数面向种间/种内保守比对（通常 > 90% 同一性）。PGR 采用 T2T-CHM13 SD
      标准（> 1 kbp、> 90% 同一性，见 4.2.1），lastz 默认参数正好匹配，无需为低同源性额外调参。
- **速度**
    - **BISER**：原生 plane-sweep 近似线性。
    - **lastz --self**：在全基因组上寻找所有局部比对开销大，但 PGR 的 `pgr lav lastz` 已用 rayon
      并行化。
    - **FastGA**：针对细菌规模优化，在人类尺度染色体上的行为需实测。
- **坐标一致性**
    - **BISER**：内部全程在 hard-masked 坐标系下工作。
    - **lastz --self / FastGA**：若直接作用于 hard-masked FASTA，输出坐标自然对齐；若作用于原始
      FASTA，则必须在转换器或 `mask` 阶段做坐标映射，增加复杂度。
- **CIGAR 语义**
    - **BISER**：CIGAR 只含 `M`/`I`/`D`（`simple_cigar`），且 `M` 同时代表 match 与 mismatch。
    - **lastz/PSL/PAF**：通常含 `=`/`X`/`I`/`D`。PGR 下游（`cluster`/`translate`）直接消费
      `cg:Z:`（`=`/`X`/`I`/`D`），仅在导出 BISER `.align` 时才合并 `=`/`X` 为 `M`（见 6.9.4）。
- **下游兼容性**
    - `cluster` 对坐标/链向敏感，对 CIGAR 不敏感；`translate` 对 CIGAR 敏感。
      因此外部比对路线至少能保证 clustering 正确，但 translate 阶段必须验证 CIGAR 在 hard-masked ↔
      original 映射时不出错。

### 6.8.7 参数选择：匹配 T2T-CHM13 SD 标准

PGR 采用 T2T-CHM13 SD 标准（> 1 kbp、> 90% 同一性，见 4.2.1），不追求 BISER 默认的 30% 错误率
（低至 70% 同一性）场景。lastz 的 UCSC preset（如 `set01`）通常面向 > 90% 同一性，正好匹配该标准，
无需为低同源性额外调参：

- **lastz 方向**
    - `set01` 默认参数（`K=3000 L=2200`）即可满足 > 90% 同一性 SD 检测；如需略提灵敏度可适度降低
      `K/L`（如 `K=1500 L=1500`），但不必放宽到 75% 同一性级别。
    - 提高 `--querydepth` 上限，避免 lastz 在高重复区域提前截断。
    - 使用 `--self` 时，lastz 自动跳过完全相同的 query/target 同方向 trivial 比对，这与 BISER 过滤
      self-overlap 一致。
- **FastGA 方向**
    - FastGA 的 `-p` 选项控制映射参数；按默认 identity 阈值（面向高同一性比对）即可，
      无需为低同一性重复放宽。
    - 确保输出包含 CIGAR（`-pafx` 或等价选项），否则无法填入 `cg:Z:` 与 error rate tag。
- **共同策略**
    - 所有外部比对都应作用于 `pgr sd mask` 输出的 hard-masked FASTA，确保坐标系一致。
    - 转换后按 `block_identity` 计算 `er`/`xm`/`id` tag（见 6.9.6），而不是直接采用 lastz/FastGA
      的内部 identity，因为它们的 identity 计算方式可能与 `block_identity` 不同。

### 6.8.8 建议新增的实现

若要在 PGR 中正式支持该替代路线，建议新增以下模块/命令：

- `src/libs/sd/from_lastz.rs`
    - 输入：LAV 文件（经 `pgr lav to-psl` 转 PSL）与 hard-masked FASTA。
    - 流程：`PSL → pgr pl ucsc（默认非 --syn）→ MAF → pgr maf to-paf → PAF`。
    - 输出：标准 PAF（12 列 + `cg:Z:` + 6.9.6 推荐 tag）。
    - 核心：调用 `pgr pl ucsc` 完成 UCSC chain/net（chaining 不依赖 lastz 内置结果），
      再用 `pgr maf to-paf`（`src/libs/paf/maf_import.rs::maf_block_to_paf`）转 PAF，并用
      `block_identity` 风格补 `er`/`xm`/`id` tag。
- `src/libs/sd/from_paf.rs`
    - 输入：PAF 文件（来自 FastGA / minimap2 / wfmash）与 hard-masked FASTA。
    - 流程：`PAF → PSL → pgr pl ucsc（默认非 --syn）→ MAF → pgr maf to-paf → PAF`（chaining 不依赖
      aligner 内置结果）。
    - 输出：补齐 6.9.6 推荐 tag 的标准 PAF。
    - 核心：先把 PAF 转为 PSL（供 `axtChain -psl` 消费，需新增 PAF→PSL），其余与 `from_lastz.rs`
      一致；不产出 `.align`，仅在需要导出给外部 BISER 工具时按 6.9.7 的 `to-align` 另写。
- `src/libs/sd/align_validator.rs`（可选）
    - 对转换后的 PAF 做一致性校验：
        - 验证 `cg:Z:` 总长度与 query/target 坐标差是否一致（类似 `biser/codon/mask.codon:18-29` 的
          `validate_cigar`）；
        - 验证 CIGAR 展开后 query/target 长度分别等于 `query_end - query_start` 与
          `target_end - target_start`。
- `src/cmd_pgr/sd/import.rs`
    - CLI 入口，例如：
        - `pgr sd import-lastz genome.hard.fa lav_out/ -o hits.paf`（内部：LAV→PSL→`pgr pl ucsc`→MAF→
          `pgr maf to-paf`）
        - `pgr sd import-paf genome.hard.fa self.paf -o hits.paf`（内部：PAF→PSL→`pgr pl ucsc`→MAF→
          `pgr maf to-paf`）
    - 命令内部可校验 uppercase run 映射表是否存在。
- 与 `pgr sd run` 的集成
    - 在 `pgr sd run` 中增加 `--aligner biser|lastz|fastga` 选项，允许用户选择用原生 BISER 算法、
      lastz --self 还是 FastGA 生成 PAF。
    - 选择 `lastz`/`fastga` 时，跳过 BISER 原生 search/align，由 `from_lastz.rs` / `from_paf.rs`
      调用 `pgr pl ucsc`（非 `--syn`）完成 chaining/refine，再进入 `cluster` / `decompose` /
      `translate`（均以 PAF 为中间格式，见 6.9）。

### 6.8.9 与 6.3 自研路线的关系

本节并非推翻 6.3 的迁移方案，而是提供**同一目标下的第二条实现路径**：

- 若追求与 BISER bit-exact 一致，按 6.3 自研 `kmer_index.rs` + `plane_sweep.rs` + `refine.rs`。
- 若追求快速复用 PGR 已有外部比对生态，按本节路线用 lastz --self 或 FastGA 生成 self-alignment，
  经 UCSC chain/net（`pgr pl ucsc`，SD 不加 `--syn`）精炼后按 6.9 转为标准 PAF，`cluster` /
  `decompose` / `translate` 直接以 PAF 为输入。
- 两条路线的最终验证标准相同：转换后的 PAF 能被 `pgr sd cluster` 正确读取，且 `translate`
  后坐标/CIGAR 与 BISER 原生输出一致。

## 6.9 用标准 PAF 替代 BISER `.align` 格式

> 本节说明为何用 PGR 原生支持的标准 PAF（Pairwise mApping Format）替代 BISER 私有的 14 列 `.align`
> BEDPE，并给出完整的字段映射与迁移方案。PAF 已定为 SD 流程的中间格式；由于 PGR 已有成熟的 PAF
> 生态，这比保留 `.align` 更符合项目架构。下面结合 BISER 下游源码与 PGR 的 PAF 实现展开。

### 6.9.1 为什么 PAF 可以替代 `.align`

BISER 的 `.align` 文件本质上是带 CIGAR 的成对比对记录，下游阶段只使用其中一部分字段：

- `cluster`（`biser/codon/cluster.codon:53-84`）使用：
    - 第 1–6 列：两个 mate 的染色体与坐标；
    - 第 7 列：`species1:species2`；
    - 第 10 列：`y` 链向（`+`/`-`）。
    - **不使用** CIGAR、错误率、max span、span。
- `translate`（`biser/codon/mask.codon:69-147`）使用：
    - 第 1–6 列：坐标；
    - 第 7 列：物种对；
    - 第 10 列：`y` 链向；
    - 第 13 列：`simple_cigar`（`M`/`I`/`D`），用于把 hard-masked 坐标映射回 original 基因组。
- `decompose`（`biser/codon/decompose.codon:226-277`）不直接读取 `.align`，而是读取
  `cluster` 输出的 FASTA；其头部格式为 `{sp}#{ch}{sr}#{st}#{ed}`，因此它依赖的是 `cluster`
  对物种/染色体/链向的编码能力，而非 `.align` 的具体列结构。

PAF 的 12 列本身就包含：query 名称/长度/起止、target 名称/长度/起止、链向、匹配碱基数、block 长度、
mapq。BISER 额外需要的物种对、CIGAR、错误率、max span、span 均可通过 PAF 的**可选 SAM-like tag**
承载。因此，**PAF 在信息量上完全覆盖 `.align`**。

### 6.9.2 `.align` → PAF 字段映射

将 BISER `.align` 的 14 列映射到 PAF：

- **query 侧（x mate）**
    - `x_name` → PAF 第 1 列 `query_name`
    - `x_start` → PAF 第 3 列 `query_start`（0-based inclusive）
    - `x_end` → PAF 第 4 列 `query_end`（0-based exclusive）
- **target 侧（y mate）**
    - `y_name` → PAF 第 6 列 `target_name`
    - `y_start` → PAF 第 8 列 `target_start`（0-based inclusive）
    - `y_end` → PAF 第 9 列 `target_end`（0-based exclusive）
- **链向**
    - BISER 保证 `x_strand` 恒为 `+`（`hit.codon:100`），因此 PAF 的 strand 列直接对应 `y_strand`：
        - `y_strand == '+'` → PAF 第 5 列为 `+`
        - `y_strand == '-'` → PAF 第 5 列为 `-`
- **物种对**
    - `.align` 第 7 列 `species1:species2` 需要额外 tag，例如 `sp:Z:hg38:panTro6`。
    - 或者把物种编码进序列名：`query_name = "hg38#chr1"`，`target_name = "panTro6#chr2"`。这与
      `cluster` 输出 FASTA 头 `{sp}#{ch}...` 的约定一致。
- **max span 与 span**
    - `.align` 第 11 列 `max(x_end - x_start, y_end - y_start)` 可用 tag `ms:i:` 记录。
    - `.align` 第 12 列 `span()`（CIGAR 总长度）可从 `cg:Z:` 直接推导，也可用 tag `sp2:i:` 缓存。
- **CIGAR**
    - `.align` 第 13 列 `simple_cigar` 使用 `M`/`I`/`D`，其中 `M` 同时代表 match 与 mismatch。
    - 标准 PAF 使用 `cg:Z:` tag，操作符为 `=`/`X`/`I`/`D`（`minimap2 --eqx` / `wfmash` 输出）。
    - 两种方案：
        1. **标准方案**：输出 `cg:Z:`（`=`/`X`/`I`/`D`），在需要 BISER-style CIGAR 时把 `=`/`X`
           合并为 `M`。
        2. **兼容方案**：同时输出 `cg:Z:` 和 `bc:Z:`（BISER-compatible CIGAR，`M`/`I`/`D`），
           `translate` 直接读取 `bc:Z:`。
- **错误率**
    - `.align` 第 8 列总错误率 `err() * 100` 可用 tag `er:f:` 记录。
    - 第 14 列 `X=...;ID=...` 可拆为 `xm:f:`（mismatch rate）与 `id:f:`（gap rate）。
    - 也可用标准 tag `nm:i:`（edit distance）配合 `block_length` 反推。

### 6.9.3 下游阶段如何直接消费 PAF

若把中间格式改为 PAF，`cluster` 与 `translate` 需要相应调整：

**`cluster`**

- 读取 PAF 12 列即可工作：
    - `query_name/target_name` → 染色体名；
    - `query_start/query_end/target_start/target_end` → 坐标；
    - `strand` → 反向判断；
    - 物种对通过序列名前缀（`species#chr`）或 `sp:Z:` tag 获得。
- 不需要 CIGAR，因此 PAF 记录是否带 `cg:Z:` 不影响 clustering。
- 输出 FASTA 头保持 `{sp}#{ch}{sr}#{st}#{ed}`，与现有 `decompose` 输入兼容。

**`translate`**

- 需要 PAF 带 CIGAR tag：
    - 若使用 `bc:Z:`（BISER-style `M`/`I`/`D`），可直接复用 `biser/codon/mask.codon:79-147` 的逻辑。
    - 若只使用 `cg:Z:`（`=`/`X`/`I`/`D`），需先把 `=`/`X` 合并为 `M`，再进行 hard-masked → original
      的坐标映射。
- PGR 已有 `src/libs/paf/cigar.rs` 可以解析 `cg:Z:`，并把 `=`/`X`/`I`/`D` 折叠为 `M`/`I`/`D`。
- 坐标映射逻辑与 `.align` 相同：因为 PAF 的 `query_start/target_start` 是 0-based inclusive、
  `query_end/target_end` 是 0-based exclusive，与 BISER 坐标系一致。

**`decompose`**

- 无需改动，因为它只读取 `cluster` 输出的 FASTA。

### 6.9.4 坐标与 CIGAR 的精确处理

**坐标系统**

- BISER `.align` 使用 0-based half-open `[start, end)`。
- PAF 同样使用 0-based：start inclusive、end exclusive（`src/libs/paf/record.rs:13-26`）。
- 因此坐标可直接拷贝，无需任何转换。
- 唯一注意点：BISER 中 `x_strand` 恒为 `+`，而 PAF 的 strand 表示 query 相对 target 的方向。当
  `y_strand == '-'` 时，PAF strand 为 `-`；target 坐标仍需按 forward 坐标系给出（与 PAF 规范一致），
  `translate` 阶段再按需反向。

**CIGAR 转换**

- 从标准 `cg:Z:` 到 BISER `simple_cigar`：
    - `=` 和 `X` 都合并为 `M`；
    - `I` 和 `D` 保持不变；
    - 连续同操作符需要合并。
- PGR 实现可直接调用 `src/libs/paf/cigar.rs`：
    - `parse_cigar()` 解析 `cg:Z:`；
    - 遍历 `CigarOp`，把 `=`/`X` 统一当作 `M` 输出；
    - 对 `I`/`D` 直接输出。
- 反向转换（BISER → PAF）需要实际序列，因为把 `M` 拆成 `=`/`X` 必须知道碱基是否相同。
  因此若要保持信息无损，输出 PAF 时应直接记录 `cg:Z:`（`=`/`X`/`I`/`D`），而不是事后再从
  `simple_cigar` 反推。

### 6.9.5 物种编码策略

PAF 本身没有物种列，需要额外约定：

- **方案 A：序列名编码物种（推荐）**
    - `query_name = "hg38#chr1"`，`target_name = "panTro6#chr2"`。
    - `cluster` 解析 `name.split('#')` 得到 species 与 chromosome。
    - 优点：不依赖自定义 tag，与 `cluster` 输出 FASTA 头 `{sp}#{ch}...` 天然一致，也便于 `pgr paf`
      索引直接复用。
    - 缺点：要求输入 FASTA 的序列名必须含物种前缀；跨项目需要统一命名规范。
- **方案 B：自定义 `sp:Z:` tag**
    - 每条 PAF 记录带 `sp:Z:species1:species2`。
    - 优点：显式、不污染序列名。
    - 缺点：需要所有下游工具解析该 tag；`pgr paf` 现有索引不识别此 tag。
- **方案 C：单独物种映射表**
    - 一个 TSV 文件：`sequence_name<TAB>species`。
    - 优点：最灵活，不改 PAF。
    - 缺点：多一个文件，容易遗漏或不同步。

**建议**：在 PGR 的 SD 流程中采用**方案 A**（序列名含物种），并在文档与 CLI 中强制要求；`sp:Z:` tag
作为可选冗余信息，方便与外部 PAF 互操作。

### 6.9.6 建议的 PAF tag schema

PAF 已定为中间格式，建议统一以下 tag：

- **必需**
    - `cg:Z:`：标准 CIGAR（`=`/`X`/`I`/`D`），用于精确坐标投影与图构建。
- **推荐**
    - `sp:Z:species1:species2`：物种对（当序列名未编码物种时使用）。
    - `er:f:0.15`：总错误率（对应 `.align` 第 8 列 / 100）。
    - `xm:f:0.10`：mismatch 错误率（对应 `.align` 第 14 列 `X=` 部分 / 100）。
    - `id:f:0.05`：gap（insertion/deletion）错误率（对应 `.align` 第 14 列 `ID=` 部分 / 100）。
    - `ms:i:10000`：max span（对应 `.align` 第 11 列）。
    - `sp2:i:10500`：span（CIGAR 总长度，对应 `.align` 第 12 列）。
- **可选**
    - `bc:Z:100M50I200M`：BISER-style CIGAR（`M`/`I`/`D`），供 `translate` 直接消费，避免运行时合并
      `=`/`X`。
    - `nm:i:1500`：编辑距离，可与 `block_length` 互推 error rate。
    - `gi:f:0.95`：gap-compressed identity（PGR `paf` 已有惯例）。

所有 tag 类型遵循 SAM specification：`Z` 表示字符串，`f` 表示 float，`i` 表示 integer。

### 6.9.7 迁移路径与模块设计

既然 PAF 已定为中间格式，建议按以下模块实现：

- `src/libs/sd/to_paf.rs`
    - 输入：BISER `.align` 14 列 BEDPE + hard-masked / original FASTA（用于把 `M` 拆分为 `=`/`X`）。
    - 输出：标准 PAF（12 列 + `cg:Z:` + 推荐 tag）。
    - 核心：
        - 坐标直接映射；
        - 序列名按 `species#chr` 编码；
        - 对 `simple_cigar` 的每个 `M` 块，从 FASTA 提取 query/target 子序列，逐碱基比较拆分为 `=`/
          `X`；
        - 计算 `er`/`xm`/`id` tag；
        - 写出 `cg:Z:`，可选写出 `bc:Z:`。
- `src/libs/sd/from_paf.rs`（升级版）
    - 在 6.8.8 的 FastGA/minimap2 PAF 导入基础上，统一补齐 6.9.6 的推荐 tag：
        - 读取 `cg:Z:`，用 `cigar_stats` 统计 match/mismatch/ins/del；
        - 读取 `sp:Z:` 或序列名前缀得到物种对；
        - 计算 `er`/`xm`/`id`/`ms`/`sp2` tag；
        - 输出标准 PAF。仅当需要导出给外部 BISER 工具时，才额外产出 14 列 `.align`（见下文
          `to-align`）。
- `src/libs/sd/cluster_paf.rs`
    - `cluster` 的 PAF 版本：直接读取 PAF，按 `query/target` 坐标与 `strand` 做 coloring，输出
      `{sp}#{ch}{sr}#{st}#{ed}` FASTA。
    - 优势：可直接利用 PGR 的 PAF 索引（`src/libs/paf/index.rs`）做区间查询，
      未来可扩展为只处理指定区域。
- `src/libs/sd/translate_paf.rs`
    - `translate` 的 PAF 版本：读取 PAF，用 `cg:Z:` 或 `bc:Z:` 做 hard-masked → original 映射，
      输出更新后的 PAF（坐标与 CIGAR 均已转换）。
    - 输出 tag 中 `bc:Z:` 自动转换为 `M`/`S`/`N`（对应 `mask.codon` 的 `S`/`N` 语义），坐标恢复
      original 基因组坐标。
- `src/cmd_pgr/sd/` 子命令调整
    - `pgr sd align`：输出 PAF。
    - `pgr sd cluster`：接受 PAF 输入。
    - `pgr sd translate`：接受 PAF 输入，输出 PAF。
    - `pgr sd to-align` / `pgr sd from-align`：仅作为与外部 BISER `.align` 互操作的转换命令，不进入
      PGR 自身 pipeline。

### 6.9.8 与 6.8 外部比对路线的关系

将 `.align` 替换为 PAF、并固定 chaining 由 UCSC chain/net 承担后，6.8 节的 lastz --self / FastGA
路线会变得更自然：

- **统一 chaining**：lastz 与 FastGA 内置的 chaining 都不采用，两条路线在自比对之后都走
  `pgr pl ucsc`（默认非 `--syn`，SD 不做共线性筛选）。
- **lastz --self**：`LAV → pgr lav to-psl → PSL → pgr pl ucsc → MAF → pgr maf to-paf → PAF`。
- **FastGA / minimap2 / wfmash**：`PAF → PSL → pgr pl ucsc → MAF → pgr maf to-paf → PAF`（需
  PAF→PSL 转换；chaining 同样不依赖 aligner 内置结果）。
- **不能跳过 chain/net**：SD 检测需要的 chaining/边界精修由 UCSC chain/net 提供，不能在 `cluster`
  阶段靠 overlap 合并替代。

因此，**UCSC chain/net 统一了各 aligner 的 chaining，PAF 统一了下游中间表示**；`pgr pl ucsc`（非
`--syn`）+ `pgr maf to-paf` 成为所有外部比对路线共用的精炼段。

### 6.9.9 优势、风险与落地建议

**优势**

- 标准格式：可被 `pgr paf`、`pgr maf to-paf`、minimap2、wfmash、FastGA 等直接消费。
- 减少私有格式维护：无需为 `.align` 单独写解析器、验证器、文档。
- 坐标系统一致：PAF 与 BISER 均为 0-based half-open，避免转换错误。
- 易于扩展：新增 tag 即可携带更多元信息，不影响旧工具（旧工具会忽略未知 tag）。

**风险**

- 需要更新 `cluster` / `translate`：虽然改动范围有限，但需仔细验证 hard-masked → original 的 CIGAR
  映射。
- 物种编码必须统一：若用户输入 FASTA 未按 `species#chr` 命名，需要强制重命名或要求 `sp:Z:` tag。
- CIGAR 语义差异：标准 `cg:Z:` 与 BISER `simple_cigar` 的转换点必须清晰文档化，避免 `translate`
  输出错误。

**落地建议**

- PAF 为下游统一中间格式：`pgr sd align`（chain/net 后）/ `pgr sd cluster` / `pgr sd translate`
  使用 PAF；`pgr sd search` 输出原始 PSL 供 chain/net 消费。PGR 自身 pipeline 不产出 `.align`。
- `pgr sd to-paf`（`.align` → PAF，用于导入 BISER 原生输出）与 `pgr sd to-align`（PAF → `.align`，
  用于导出给外部 BISER 工具）仅作互操作，不进入主流程。
- 迁移验证：对照 BISER 原生 `.align` 输出，确认 PAF 的坐标、链向、CIGAR、错误率 tag 语义一致（见
  6.9.4）。

## 6.10 历史项目参考：App-Egaz 与 intspan/cmd_linkr

> 本节分析作者早期两个项目（`~/Scripts/App-Egaz`、`~/Scripts/intspan/src/cmd_linkr/`）中与 SD
> 检测相关的流程，找出可与 BISER 迁移方案相互参照或借鉴的部分。这两个项目采用“lastz/blastn 找种子 +
> link 图聚类 + MSA 精修”的路线，与 BISER 的“k-mer plane-sweep + PST chaining + 分解”路线差异较大，
> 但在 pipeline 组织、图聚类、区间操作上仍有参考价值。

### 6.10.1 App-Egaz 的 SD 检测流程

`App-Egaz`（Perl 项目，`lib/App/Egaz/Command/`）是一个围绕 UCSC chain/net 与 lastz/blastn
的基因组比对流程工具。其自我 SD 检测流程在 `share/3_self.tt2.sh` 与 `doc/Scer-self.md` 中有完整描述，
可概括为以下阶段：

**阶段 1：lastz 自比对与 chain/net 精炼**

- 调用 `egaz lastz --isself --set set01 -C 0` 做全基因组自比对（`share/1_self.tt2.sh:15-19`）。
- 调用 `egaz lpcnam` 跑通 `lav → psl → chain → net → axt` 流程
  （`lib/App/Egaz/Command/lpcnam.pm:101-393`）。
    - 使用 kent-tools 的 `axtChain`、`chainAntiRepeat`、`chainMergeSort`、`chainNet`、`netToAxt`
      等外部程序。
    - 参数 `--lineargap loose --minscore 1000` 用于人类尺度；酵母等小基因组可更宽松。
- 输出为 `axtNet/*.axt.gz`，即经 net 过滤后的 pairwise 局部比对。

**阶段 2：提取初始精确/近精确拷贝**

- `fasr axt2fas` 把 axtNet 转成 block FA（多序列比对块）。
- `fasr filter --ge 1000` 保留长度 ≥1000 bp 的块（`share/3_self.tt2.sh:47-48`）。
- `fasr link axt.fas` 从 block FA 中提取 bilateral links，每行两个 range：`chr(strand):start-end`。
- 输出 `links.lastz.tsv`。

**阶段 3：blastn 扩展寻找更多 paralogs**

- 从初始拷贝中去除重复、去除含过多 N 的序列，得到 `axt.gl.fasta`（`share/3_self.tt2.sh:58-62`）。
- `egaz blastn axt.gl.fasta genome.fa` 把候选 paralogs 比对回全基因组。
- `egaz blastmatch` 按 coverage ≥0.95 提取命中区域（`lib/App/Egaz/Command/blastmatch.pm:63-227`）。
- 提取命中序列后与原始候选合并为 `axt.all.fasta`。
- 再次 `egaz blastn axt.all.fasta axt.all.fasta` 做 all-vs-all 比对。
- `egaz blastlink -c 0.95` 把 blastn 结果转成 links（`lib/App/Egaz/Command/blastlink.pm:61-141`）。

**阶段 4：link 图清理与聚类**

- `linkr sort`：对 links 去重并排序。
- `linkr clean`：合并链向、去除嵌套 links、按 `--bundle 500` 合并重叠 links
  （`~/Scripts/intspan/src/cmd_linkr/clean.rs`）。
- `rgr merge -c 0.95`：按双向 coverage ≥0.95 合并重叠 range（`~/Scripts/intspan/src/cmd_rgr/merge.rs`）。

- `linkr clean -r links.merge.tsv --bundle 500`：用 merge 结果替换原始 ranges 后再次清理。
- `linkr connect -r 0.9`：把 bilateral links 连接成 multilateral connected components
  （`~/Scripts/intspan/src/cmd_linkr/connect.rs`）。
- `linkr filter -r 0.8`：按长度差异 ratio ≥0.8 过滤（`~/Scripts/intspan/src/cmd_linkr/filter.rs`）。

**阶段 5：MSA 精修**

- `fasr create genome.fa links.filter.tsv -o multi.temp.fas`：按 links 从基因组提取序列并生成 block
  FA。
- `fasr refine multi.temp.fas -o multi.refine.fas --msa mafft -p 8 --chop 10`：用 mafft
  做多序列比对精修。
- `fasr link multi.refine.fas`：从精修后的 block FA 重新提取 links。
- 对 pairwise best links 再次精修，得到最终 `pair.refine.fas`。

### 6.10.2 intspan/cmd_linkr 的 link 图操作

`intspan`（Rust 项目）中的 `cmd_linkr` 模块提供了一套区间 link 操作，可直接视为 SD 聚类的图工具：

- **Link 文件格式**
    - bilateral link：`range_0\trange_1` 或 `range_0\trange_1\thit_strand`。
    - multilateral link：`range_0\trange_1\t...\trange_n`。
    - range 格式：`chr(strand):start-end`，例如 `chr1(+):1000-2000`。坐标为 1-based inclusive。
- **`linkr sort`**（`src/cmd_linkr/sort.rs:28-61`）
    - 用 `BTreeSet` 去重；
    - 调用 `intspan::sort_links` 对每行内 ranges 及行之间排序。
- **`linkr clean`**（`src/cmd_linkr/clean.rs:59-391`）
    - 把带链向的 ranges 统一规范化为正链；
    - 去除完全嵌套的 links（两个 range 均被另一条 link 的两个 range 包含）；
    - `--replace`：用 `rgr merge` 的结果替换 ranges；
    - `--bundle N`：当两条 link 在两个端点上都有 ≥N bp 重叠时，用图连通分量合并为一条 link；
    - 去除 self-link（同染色体重叠 > 50%）。
- **`linkr connect`**（`src/cmd_linkr/connect.rs:53-297`）
    - 把 bilateral links 当作无向图边，用 `petgraph::algo::tarjan_scc` 找连通分量；
    - 根据边权重（`+`/`-`）为每个连通分量内的节点分配链向；
    - `--ratio 0.9`：若两节点长度差异过大则断开边；
    - 输出 multilateral links（每个连通分量一行）。
- **`linkr filter`**（`src/cmd_linkr/filter.rs:49-103`）
    - `--number`：按每行 range 数量过滤（如 `--number 2-10` 只保留 copy number 2–10 的 SD）；
    - `--ratio`：按行内最大/最小 range 长度比过滤。
- **`rgr merge`**（`src/cmd_rgr/merge.rs:61-212`）
    - 对单个染色体上的 ranges 构建重叠图；
    - 若两个 range 的双向 coverage 均 ≥ `--coverage 0.95`，则连边；
    - 输出 `original_range\tmerged_range` 替换表，供 `linkr clean --replace` 使用。

### 6.10.3 与 BISER 的对比

**整体策略差异**

- **BISER**：
    - 在 hard-masked 基因组上建 exact 2-bit k-mer 索引；
    - plane-sweep 快速找出 putative SD pairs；
    - PST chaining + DP 精修边界与 CIGAR；
    - interval coloring 聚类；
    - k-mer frequency 分解 elementary SDs。
- **App-Egaz/intspan**：
    - 用 lastz 做全基因组自比对，经 UCSC chain/net 过滤；
    - 从精炼比对中提取长 block（≥1000 bp）作为种子；
    - blastn 把种子扩展为更大 paralog 集合；
    - all-vs-all blastn 生成 pairwise links；
    - linkr 图操作（clean/connect/filter）完成聚类；
    - mafft MSA 精修对齐边界。

**灵敏度与准确性**

- **BISER 优势**：
    - k-mer + plane-sweep 近似线性，适合人类尺度大基因组；
    - 默认允许较高 error rate（30%，低至 70% 同一性），能检测古老 SD（PGR 当前采用 > 90% 标准，
      该能力不在需求范围内，见 4.2.1）；
    - 有系统化的 decomposition 步骤，输出 elementary SDs。
- **App-Egaz 劣势（也是用户观察到效果不如 BISER 的原因）**：
    - 依赖 lastz `--set set01` 默认参数（面向 > 90% 同一性），在 PGR 的 > 90% 标准下正好匹配，
      但无法覆盖 BISER 的低同源性场景；
    - 中间步骤过多（lav/psl/chain/net/axt/blastn），每个环节都可能丢失信息；
    - blastmatch/blastlink 仅按 coverage 过滤，没有 BISER 的 error model；
    - chain/net 若启用共线性筛选（`netFilter -syn`）会去除大量非共线性重复，而 SD 本身常伴随重排
      （SD 路线因此走非 `--syn`，见 6.8）；
    - 没有类似 BISER `decompose` 的 elementary SD 分解，只能得到 multilateral link clusters。

**可借鉴之处**

尽管效果不如 BISER，以下流程与工具对 PGR 的 SD 迁移仍有参考价值：

1. **lastz 自比对作为种子来源**
    - App-Egaz 的 `egaz lastz --isself` 与 6.8 节的 `pgr lav lastz --self` 本质相同；
    - 可作为 BISER `search` 阶段的一种低灵敏度、高特异度替代，尤其适合小基因组或细菌；
    - 其 `lpcnam` 流程（lav → psl → chain → net → axt）可被 PGR 的 `pgr lav to-psl` + `pgr pl ucsc`
      替代（后者内部串接 axtChain/chainNet/netToAxt 等步骤，见 6.8.4）。
2. **blastn 扩展策略**
    - 用 lastz/blastn 的可靠长 hits 作为种子，再用 blastn 回贴基因组寻找更多 paralogs，
      是一种经典的“seed-and-expand”策略；
    - 在 PGR 中，若 BISER 的 k-mer search 漏掉某些低复杂度或高度分化区域，可考虑用 lastz/blastn
      作为补充种子源。
3. **link 图聚类思想**
    - BISER `cluster` 用 interval coloring 把重叠 SDs 分到同一组；
    - App-Egaz/linkr 用图连通分量（clean/connect/filter）实现类似功能，且更直观；
    - PGR 的 SD 聚类阶段可参考 `linkr connect` 的图方法，把 PAF 记录当作边，构建 SD 图，再输出
      clusters。
4. **range 操作工具链**
    - `linkr clean` 的嵌套去除、bundle 合并、self-link 过滤；
    - `linkr filter` 的 copy number 与长度比过滤；
    - `rgr merge` 的双向 coverage 合并；
    - 这些操作本质上是对 SD hit 集合的后处理，可在 PGR 中复现为 `pgr sd clean` / `pgr sd filter`
      等子命令。
5. **分区处理大染色体**
    - App-Egaz 的 `partition` 命令（`lib/App/Egaz/Command/partition.pm:57-83`）把大染色体切成带
      overlap 的小段（默认 10 Mbp + 10 kbp overlap），以便 lastz 并行处理；
    - BISER 的 `MAX_CHROMOSOME_SIZE` 也有类似切片逻辑（`biser/__main__.py:114`），PGR 实现
      `kmer_index.rs` 时可参考两者，选择按固定长度或按 N/gap 边界分区。
6. **MSA 精修边界**
    - App-Egaz 用 `fasr refine --msa mafft` 精修 cluster 内序列边界；
    - PGR 已有 POA 引擎（`src/libs/poa`），可在 `pgr sd refine` 中替代 mafft，实现不依赖外部
      aligner 的边界精修；
    - 但 App-Egaz 的“提取候选序列 → 多序列比对 → 重新生成 links”循环值得借鉴，可用于改进 BISER
      `cluster` 后、`decompose` 前的边界质量。

### 6.10.4 对 PGR SD 迁移的启示

- **不要直接复刻 App-Egaz 流程**：其 lastz + chain/net + blastn 的多步组合、coverage-only 过滤、
  缺少 elementary SD 分解，使整体灵敏度不如 BISER；PGR SD 复用其中的 chain/net 段（非 `--syn`），
  但下游接 BISER 的 cluster/decompose/translate，不照搬 App-Egaz。
- **可取用的组件化思路**：
    - `src/libs/sd/seed.rs`：封装 lastz/blastn 种子生成；
    - `src/libs/sd/cluster_graph.rs`：用图连通分量替代 BISER interval coloring，输出 multilateral
      SD clusters；
    - `src/cmd_pgr/sd/clean.rs` / `filter.rs`：提供 `linkr clean` / `linkr filter` 风格的 SD hit
      后处理；
    - `src/libs/sd/refine_poa.rs`：在 cluster 阶段用 POA 精修边界（补充 BISER 的 `align.refine`）。
- **与 6.3 / 6.8 / 6.9 的关系**：
    - 6.3 的 BISER 原生路线（k-mer search + PST refine）为远期目标；
    - 6.8 的 lastz/FastGA 外部比对路线（UCSC chain/net + PAF）是近期主路径，可看作 App-Egaz 阶段 1
      的现代化、标准化版本；
    - 6.9 的 PAF 替代 `.align` 可让 App-Egaz/linkr 的图工具更直接地消费 BISER 输出；
    - 6.10 的历史项目则为后处理（聚类、过滤、精修）提供了额外思路。

### 6.10.5 需要避免的历史项目缺陷

- **共线性筛选与 SD**：UCSC chain/net 的 `netFilter -syn` 是为 syntenic alignment 设计的，
  会系统性地丢弃非共线性重复；SD 路线因此使用 `pgr pl ucsc` 默认（非 `--syn`）路径，只做 chain/net
  结构化与边界精修，不做共线性筛选。
- **coverage-only 过滤**：blastmatch/blastlink 仅按 coverage 过滤，没有区分 match/mismatch/gap，
  容易保留低质量 hits；应使用 BISER 的 `block_identity` 或 PAF 的 `er`/`gi` tag 做质量控制。
- **1-based inclusive 坐标**：App-Egaz 的 range 格式使用 1-based inclusive，与 BISER/PGR 的 0-based
  half-open 不同，迁移时需小心转换。
- **缺少 elementary SD 分解**：App-Egaz 只能输出 clusters，没有 BISER `decompose` 的 nested SD
  拆解能力；PGR 应保留 BISER 的 decomposition 阶段。

## 6.11 基于覆盖度的重复区检测：`pgr-repeat.sh` 与 BISER `search` 的关系

`scripts/pgr-repeat.sh` 是 PGR 自带的一个 Cactus-style 重复屏蔽示例脚本。用户提出：
该脚本找出的"重复区"包含转座子与 SD 两类，若用 RepeatMasker 区间做差集得到候选 SD 区，
再对这些候选区做相互匹配，是否就等价于 BISER 流程前面"找到潜在 SD 区间"的步骤？本节分析这一思路与
BISER `search` 阶段的异同，并给出集成建议。

### 6.11.1 `pgr-repeat.sh` 的工作流程

该脚本本质上是一个"滑动窗口自比对 + 覆盖度阈值"的重复检测器，具体流程如下：

1. **窗口化**：`pgr fa window -l 200 -s 100` 把基因组切成大量 200 bp、步长 100 bp 的重叠窗口，
   输出头为 `>seq_name:start-end`（1-based inclusive）。
2. **自比对**：`pgr lav lastz` 以基因组染色体为 target、窗口为 query 做全基因组比对；默认使用
   `--preset set01`（`C=0 E=30 K=3000 L=2200 O=400 Y=3400 Q=similar`，源自 UCSC Human vs Chimp）。
3. **格式转换**：`pgr lav to-psl` 把 LAV 转为 PSL；`pgr psl lift --q-sizes`
   将窗口坐标提升回原始基因组坐标。
4. **提取范围**：`pgr psl to-range` 把 query 端比对坐标转为 `.rg` 范围（1-based inclusive）。
5. **覆盖度过滤**：`spanr coverage -m 4` 计算每个碱基被多少条比对覆盖，输出深度 ≥ 4 的区域。

其逻辑是：2x 覆盖窗口自身会产生约 2 的深度基线；若某区域存在 paralog，
则来自同源拷贝的额外比对会使深度 ≥ 4，从而把重复/ duplicated 区域标记出来。

### 6.11.2 "重复区"与 SD 的关系

`pgr-repeat.sh` 检测的是广义的**高深度重复信号**，它至少包含两类：

- **可移动元件（TE）**：转座子、逆转录转座子等，通常已被 RepeatMasker 标注。
- **Segmental duplications（SD）**：较大的（≥1 kbp）、高同一性的 paralog 片段。

此外还可能混入：

- 低复杂度序列 / 卫星 DNA；
- 多拷贝基因家族的小外显子；
- 假基因、rDNA 簇等。

因此脚本输出的 `mask_regions.json` 并不是纯 SD 集合。用户提出的"用 RepeatMasker 区间做差集"是合理的
**第一步过滤**：减去已知的 TE 区间后，剩余的高深度区间更有可能是 SD 候选区。但这仍不足以区分 SD
与其他非 TE 重复，需要后续步骤做长度、同一性、成对关系验证。

### 6.11.3 与 BISER `search` 阶段的等价性分析

**概念层面：部分等价**

BISER `search` 阶段的输出是**putative SD pairs**（成对的同源区间），
其目标是快速缩小后续精确比对的搜索空间。`pgr-repeat.sh` + TE 差集后得到的是**候选 SD 区间集合**
（单端区间），它回答了"基因组哪些地方像是重复的"，但还没有回答"这些重复区之间如何成对匹配"。

因此：

- 若把 BISER `search` 理解为"找出基因组中可能属于 SD 的候选区域"，那么 `pgr-repeat.sh` 的输出（经
  TE 过滤后）可以作为其替代品。
- 若把 BISER `search` 严格理解为"输出成对的 putative SD mates 并带有近似坐标"，那么 `pgr-repeat.sh`
  还需要额外一步：把这些候选区间相互比对，才能产生可比拟的成对关系。

**算法层面：不等价**

- **BISER `search`**：基于 exact k-mer（`KMER_SIZE=14`）+ winnowing + plane-sweep，使用有序 Jaccard
  下界做过滤。对 error rate 容忍度高（默认 `MAX_ERROR=0.3`，即低至 70% 同一性），能检测较古老、
  分化较严重的 SD（PGR 采用 > 90% 标准，该能力不在需求范围内）。
- **`pgr-repeat.sh`**：基于 lastz 局部比对 + 覆盖度阈值。lastz 的 `set01` 参数（`K=3000 L=2200`）
  面向 > 90% 同一性，与 PGR 的 SD 标准匹配；200 bp 窗口分辨率限制了边界精度。

**输出层面：不等价**

- BISER `search` 输出的是 bedPE-like 的 putative hit pairs，包含坐标、链向、近似长度。
- `pgr-repeat.sh` 输出的是单端深度区间（JSON runlist），没有直接给出哪些区间互为 paralog。

### 6.11.4 "覆盖度 + TE 差集 + 自匹配"路线的可行性

用户提出的完整路线可以概括为：

1. 运行 `pgr-repeat.sh` 得到高深度重复区；
2. 用 RepeatMasker 区间做差集，得到候选 SD 区间；
3. 把这些候选区间从基因组提取出来，相互做 all-vs-all 自比对；
4. 将自比对结果作为 putative SD pairs，供后续 refine / cluster / decompose 使用。

**这条路是可行的，且与 BISER `search` 在功能上互补，但需要注意以下几点**：

1. **候选区间的长度与质量**
    - `pgr-repeat.sh` 的 200 bp 窗口会得到大量碎片化的候选区间，需要合并相邻或重叠区间（可用
      `spanr merge` 或 PGR 的 range 工具）。
    - 建议过滤掉过短（如 < 1 kbp）的区间，因为 SD 定义通常要求 ≥ 1 kbp。
    - 候选区间内部可能仍残留未标注的 TE 片段，建议结合 `TopKPurity`（`src/libs/ds/top_k_purity.rs`）
      做低复杂度过滤。
2. **自比对的计算代价**
    - 若候选 SD 区间数量为 N，all-vs-all 自比对是 O(N²)。人类基因组经 TE
      过滤后通常仍有数万个候选区间，直接全矩阵比对开销巨大。
    - 可以先对候选区间建 k-mer 索引（复用 `src/libs/sd/kmer_index.rs` 的设计），
      只让可能相似的区间对进入 lastz/minimap2 精确比对。
    - 也可以把候选区间作为 query，全基因组作为 target 做 one-vs-all 回贴（类似 App-Egaz 的 blastn
      扩展策略），降低组合爆炸。
3. **lastz 参数需要调整**
    - `set01` 的 `K=3000` 阈值较高，对 SD 检测偏严格。建议对自匹配步骤使用更敏感的参数：
        - 降低 `K/L`（如 `K=1500 L=1500` 或更低）；
        - 使用 `--self` 模式（`pgr lav lastz --self`）而非 target/query 模式；
        - 或改用 `minimap2 -DP -k19 -w19 -m200`（如 `doc/Scer-self.md` 中的示例）输出 PAF，再按 6.8
          经 `PAF → PSL → pgr pl ucsc（非 --syn）→ MAF → pgr maf to-paf → PAF` 完成精炼。
4. **坐标系统一致性**
    - `pgr fa window`、`pgr psl lift`、`pgr psl to-range` 均使用 1-based inclusive 坐标；
    - BISER 内部与 PAF 使用 0-based half-open；
    - 在把候选区间喂给自比对工具前，必须统一转换为 0-based half-open，避免 off-by-one 错误。
5. **无法直接替代 BISER `search` 的 colinearity 保证**
    - BISER 的 ordered Jaccard 通过 k-mer 顺序隐式保证 putative pairs 大致共线；
    - 覆盖度方法只告诉你"这里重复"，不保证两个候选区间之间的 match 是共线的；
    - 自比对后的结果经 UCSC chain/net（`pgr pl ucsc`，非 `--syn`）做 chaining/边界精修（见 6.8）；
      SD 不做共线性筛选，重排/inversion 形式的 SD 作为独立 chain/net 条目保留，交由 `cluster`/
      `decompose` 处理。

### 6.11.5 优势与局限

**优势**

- **复用现有 PGR 命令**：`pgr fa window`、`pgr lav lastz`、`pgr psl lift`、`pgr psl to-range`
  都已存在，无需新增核心算法。
- **TE 控制更直接**：通过 RepeatMasker 差集显式去除已知转座子，比 BISER 的 soft-mask 依赖更清晰、
  可调试。
- **结果直观**：覆盖度深度是生物学上易解释的指标，便于阈值调优。
- **与 Cactus RepeatMasking 流程一致**：该脚本本来就是为 Cactus 重复屏蔽设计的，可与现有基因组比对
  pipeline 无缝衔接。

**局限**

- **灵敏度受 lastz preset 限制**：`set01` 面向 > 90% 同一性，与 PGR SD 标准匹配，但无法覆盖 BISER
  的低同源性场景（PGR 当前不需要）。
- **分辨率受窗口大小限制**：200 bp 窗口会模糊 SD 边界，需要后续 refine。
- **计算成本高**：全基因组 lastz + 候选区间 all-vs-all 自比对的总开销通常高于 BISER 的 k-mer
  plane-sweep。
- **不能直接输出 pairs**：需要额外一步自匹配才能产生 putative SD mates。
- **假阳性来源多**：未标注 TE、低复杂度区、基因家族都会表现为高深度，需要多层过滤。

### 6.11.6 对 PGR SD 流程的集成建议

若要在 PGR 中把 `pgr-repeat.sh` 的思路固化为 BISER `search` 的替代或补充路线，建议新增以下模块/命令：

- `src/libs/sd/coverage.rs`
    - 封装"滑动窗口 → 自比对 → lift → 覆盖度"逻辑，输出候选 SD 区间（0-based half-open BED）。
    - 参数化窗口大小、步长、深度阈值、lastz preset，方便针对不同基因组调参。
- `src/libs/sd/subtract_repeatmasker.rs`
    - 读取 RepeatMasker `.out` / `.gff` / BED，与候选 SD 区间做差集，输出过滤后的候选区。
    - 可复用 `intspan` 的区间操作思路，但统一使用 0-based half-open 坐标。
- `src/cmd_pgr/sd/coverage.rs`
    - CLI：`pgr sd coverage <genome.fa> -o candidates.bed`，内部调用 `libs/sd/coverage.rs`。
- `src/cmd_pgr/sd/search.rs` 的 `--mode coverage` 选项
    - 在 BISER 原生 k-mer search 之外，提供 `coverage` 模式作为备选；
    - 输出格式与 `--mode lastz` 一致（PSL，见 6.9），供 `pgr sd align` 的 UCSC chain/net 消费；
      不直接产出 PAF。
- `src/cmd_pgr/sd/mask.rs`
    - 把 `pgr-repeat.sh` 当前输出 `mask_regions.json` 转换为标准的 BED/runlist，便于与 RepeatMasker
      结果做集合运算。

### 6.11.7 与 6.3 / 6.8 / 6.10 的关系

- **6.3 的 BISER 原生路线**：远期目标（k-mer search + PST refine），在外部路线稳定后再实现；
  理论上更适合大基因组场景。PGR 当前采用 > 90% 同一性标准（见 4.2.1），BISER 的高 error-rate
  能力不在需求范围内。
- **6.8 的外部比对路线**：近期主路径（UCSC chain/net + PAF）；`pgr-repeat.sh` 可视为 6.8 中
  `lastz --self` 思想的延伸，但它从“找重复区”出发，而不是直接输出 pairwise alignments。
- **6.10 的 App-Egaz 种子扩展**：App-Egaz 用 lastz/blastn 找可靠长 hits 作种子，再用 blastn
  回贴扩展；`pgr-repeat.sh` 的"候选区自匹配"步骤与此类似，只是候选区来源从 lastz hits
  改为覆盖度高深区。

**结论**：`pgr-repeat.sh` 经过 TE 差集和候选区自匹配后，
**可以在功能上近似 BISER `search` 阶段输出的 putative SD 区间集合**，但算法机制、灵敏度、
输出形式均有差异。最稳妥的做法不是直接替换，而是把它作为 BISER k-mer search 的**补充验证路线**或
**小基因组快速原型路线**，在 PGR 中封装为 `pgr sd coverage` / `pgr sd search --mode coverage`。

## 7. 参考文献

- Išerić H, Alkan C, Hach F, Numanagić I.
  *Fast characterization of segmental duplication structure in multiple genome assemblies*
  . Algorithms Mol Biol. 2022;17:4.
  [https://doi.org/10.1186/s13015-022-00210-2](https://doi.org/10.1186/s13015-022-00210-2)
- Išerić H, Alkan C, Hach F, Numanagić I.
  *BISER: Fast Characterization of Segmental Duplication Structure in Multiple Genome Assemblies*
  . WABI 2021. LIPIcs, Vol. 201, 15:1–15:18.
  [https://drops.dagstuhl.de/opus/volltexte/2021/14368/pdf/LIPIcs-WABI-2021-15.pdf](https://drops.dagstuhl.de/opus/volltexte/2021/14368/pdf/LIPIcs-WABI-2021-15.pdf)
- Vollger M. R. et al. *Segmental duplications and their variation in a complete human genome* . Science.
  2022;376:eabj6965.[https://doi.org/10.1126/science.abj6965](https://doi.org/10.1126/science.abj6965)
- Vollger M. R. et al. *Increased mutation and gene conversion within human segmental duplications*
  . Nature. 2023;617:335–344.
  [https://doi.org/10.1038/s41586-023-05895-y](https://doi.org/10.1038/s41586-023-05895-y)
- BISER GitHub Repository: [https://github.com/0xTCG/biser](https://github.com/0xTCG/biser)
