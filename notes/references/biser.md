# BISER 源码与论文分析

> 源自对 `biser-master/` 目录源码及 published paper 的通读。目的：理解 BISER 在
> segmental duplication (SD) 检测与分解中的算法设计，并为 pgr 中重复/同源区域分析提供参考。

> **与 pgr 的关系**：pgr 的 SD 管线（`pgr sd`）已按 BISER 的阶段结构落地为
> `search → align → cluster → decompose → cover`（外加 `cross` 与 `run`）。检测阶段用 pgr
> 自身引擎替代 BISER 原生 search/align：`--engine pgi` 为**原生**（pgr 自研种子链对齐，
> 默认 k=40 + syncmer 8/5），`--engine lastz` 为**外部** lastz；`cluster`/`decompose`/`cover`
> 则直接复刻 BISER 的算法与常量。CLI 与设计决策见 `docs/sd.md` 与 [[../design/sd.md]]；
> 本文档仅分析 BISER 本身。

## 1. BISER 概览

### 1.1 工具定位

- **工具名称**: BISER (Brisk Inference of Segmental duplication Evolutionary stRucture)，版本 v1.4。
- **功能**: 在单个或多个基因组中快速检测 segmental duplications (SDs)，并将检测到的 SDs 分解为
  elementary SDs 与 core duplicons。
- **核心假设**: 输入基因组应预先 soft-masked（低复杂度/重复序列用小写表示）；BISER
  会在内部将其转换为 hard-masked 序列进行分析。
- **与 SEDEF 的关系**: BISER 是 SEDEF 的继任者，继承了其 SD error model，但在算法上改用线性
  plane-sweep 替代 MinHash，从而获得数倍加速。
- **版本演进**（README Changelog，v1.4，2023-03）: 切换为 Codon 实现；并
  "change of alignment refinement heuristics (should be faster now)"——即 §3.3.3 的 `refine()`
  稀疏 DP + 分块 `bio.seq.align` 精修。README 提示 v1.4 产出的 SD 可能与旧版略有差异。

### 1.2 输入输出

**命令行接口（由 Python wrapper 提供）**:

```bash
biser -o <output> -t <threads> <genome1.fa> [genome2.fa ...]
```

主要阶段（`biser/__main__.py`）:

1. `mask`: 若输入非 hard-masked（`--hard/-H` 可跳过），将 lowercase bases 过滤，
   生成 hard-masked 基因组。
2. `search`: 对每条染色体做 putative SD detection。
3. `align`: 对 putative SD pairs 做局部比对与边界精修。
4. `cross_search` / `cross_align`: 多基因组时，将每个基因组的 SD 映射到其他基因组
   （先对每个基因组的已比对 hits 跑 `extract` 子命令，把 SD 覆盖区域抽成 `<sp>.bed.regions.txt`
   作为后续扫描的 query 区间；`cross_search` 只对无序对 `g1 < g2` 各执行一次）。
5. `cluster` / `decompose`: 聚类重叠 SD 并分解为 elementary SDs（`cluster` 在
   `decompose` 阶段内部先执行）。
6. `translate`: 将 hard-masked 坐标映射回原基因组坐标。

> `--no-decomposition` 可跳过阶段 5（此时不生成 `.elem` 文件）；`biser test`（内部 `hello`
> 子命令）做构建自检。`extract` 是 cross_search 内部的辅助子命令，不单独出现在主流程列表里。

**输出格式**:

- 主输出为 BEDPE 格式，描述 SD mate 的坐标、链向、CIGAR、error rate 等。
- `.elem` 文件记录 elementary SD decomposition 结果，包括 core duplicon 标记。

**BEDPE 字段**（`hit.codon` `Hit.__str__`，行 155-171）:

| 字段 | 含义 |
|------|------|
| 前 6 列 | 标准 BEDPE：`chr1, start1, end1, chr2, start2, end2`（半开区间、0-indexed） |
| `reference` | 两 mate 基因组名，以 `:` 分隔（如 `hs:hs`） |
| `score` | 总 alignment error（%）= `mis_err + gap_err` |
| `strand1/2` | `+` / `-` |
| `max_len` | 较长 mate 的长度 |
| `aln_len` | `span()` = CIGAR 各 op 长度之和 |
| `cigar` | 由 `simple_cigar()`（行 173-184）生成，把 `=`/`X` 连续段折叠为单个 `M` |
| `optional` | `X=<mismatch%>;ID=<gap%>` 两个字段 |

> 注意输出 CIGAR 用 `M`（而非 `=`/`X`）表示匹配/错配；`=`/`X` 的细分仅用于 error rate 统计
> （`mis_err`/`gap_err`/`err`，行 233-238）。

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

> 注：代码中的 winnowing 实现与上述“滑动窗口最小值”语义有偏差（见 §3.3.1），实测指纹密度约为
> `2|G|/(w+1)` 的一半，且不再保证每个窗口内必然命中。

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

1. 将所有 SD 覆盖区域 `R` 取出，按重叠关系做区间染色聚类（`cluster.codon`）。
2. 对每个 cluster，构建 `R` 上所有 k-mer 的索引 `I_k`（不使用 winnowing，默认 `k=10`）。
3. 用与 search 类似的 plane-sweep，在 `R` 上扫描相同 k-mer 的多个位置，通过距离阈值 `d_g = 50`
   将位置链式合并，得到 putative elementary SDs。
4. 当一个 elementary SD 不再能延伸时，输出长度超过 `μ`（默认 100bp）的拷贝集合。

**Core duplicon**: 定义为能覆盖所有 SD 的 elementary SD 最小集合。BISER 使用贪心 set-cover 近似算法
（`cover.py` 中的 `greedy_set_cover`）识别 core duplicons，并在 `.elem` 文件中标记为 `CORE`。
（每个 core 需覆盖至少 2 个 SD，`len(sds) < 2` 的集合被跳过。实现细节见 §3.3.5。）

### 2.5 多基因组扩展

对每个基因组独立执行 detection 与 alignment 后，BISER 将每个基因组的 SD 映射到其他基因组
（cross_search）。这避免了把两个基因组之间的保守区域误判为 SD。cross_search 本质上仍是 plane-sweep：
以基因组 A 为参考建索引，用基因组 B 中已提取的 SD 区域作为 query 进行扫描。

## 3. BISER 代码实现

### 3.1 构建与架构

- **实现语言**: 核心算法用 [Codon](https://github.com/exaloop/codon/)（带 Seq plugin）编写，编译为
  `biser/exe/biser.exe`；流程调度与多进程用 Python (`biser/__main__.py`) 完成。
- **构建**: `setup.py` 中的 `CustomBuild` 调用
  `codon build -plugin seq biser/codon/__init__.codon -release -o <build_lib>/biser/exe/biser.exe`
  （安装后即 `biser/exe/biser.exe`），并拷贝 Codon runtime 动态库（`libcodonrt`、`libomp`）。
- **依赖**: Python 侧需要 `tqdm`、`ncls`、`multiprocess`；运行时需要 `samtools faidx`。

### 3.2 源码模块与论文算法对应

- `biser/codon/__init__.codon`（入口）：子命令分发、参数解析；常量实际定义于 `search.codon`，此处
  `import` 并用 `global` 声明以便命令行覆盖（`MAX_ERROR`、`KMER_SIZE`、`WINNOW_SIZE` 等）；`align`
  子命令对 span > 300 kb 的 putative SD 直接跳过（其中 `max_iter = 50` 因紧跟 `continue` 而不会生效）。
- `biser/codon/search.codon`（Putative SD detection）：k-mer 索引、winnowing、plane-sweep、putative
  SD 输出；含 `cross_search`。
- `biser/codon/hit.codon`（数据结构）：`Hit` / `Chromosome` / `Locus` 定义、CIGAR 操作、错误率计算、
  `merge`。
- `biser/codon/chain.codon`（Chaining）：`PrioritySearchTree` 实现、`chain()` 函数完成 `O(n log n)`
  锚点链构建。
- `biser/codon/align.codon`（Alignment refinement）：`generate_anchors()` 10-mer 锚点生成、
  `refine()` 基于 DP 的边界精修。
- `biser/codon/cluster.codon`（SD clustering）：重叠 SD 的区间染色聚类（interval coloring，
  颜色合并表近似 union-find）、提取每个 cluster 的 FASTA。
- `biser/codon/decompose.codon`（SD decomposition）：k-mer chaining 分解 elementary SDs。
- `biser/codon/mask.codon`（预处理/后处理）：hard-mask 生成、hard-masked 与原基因组坐标互转；
  `translate()` 在映射回原坐标时会把 CIGAR 的 `M` 按 lowercase 区间重新切分为 `M`/`S`/`N`
  （`I`/`D` 同理切分为 `I`/`S`、`D`/`N`）。
- `biser/cover.py`（Core duplicon）：用 `ncls` 做区间重叠查询，再用贪心 set cover 找 core duplicons。
- `biser/__main__.py`（Pipeline）：多进程任务拆分、临时目录管理、各阶段串接。

### 3.3 关键实现细节

#### 3.3.1 Plane-sweep in `search.codon`

> 关键常量 `MAX_ERROR=0.3`、`MAX_EDIT_ERROR=0.15`、
> `KMER_SIZE=14`、`WINNOW_SIZE=16`、`MAX_DISTANCE=250`、`MAX_SD_LEN=2_000_000`、
> `QUERY_THRESHOLD=100`、`REF_THRESHOLD=500`、`MAX_EXTEND=5_000`；§3.3.2 的
> `MATCH_SCORE=4`、`MAX_CHAIN_GAP=210`、`MIN_UPPERCASE_MATCH=90`、
> `MIN_READ_SIZE=700`；§3.3.3 的 `MATCH=10/MISMATCH=1/GAP=0.5/GAPOPEN=100`、
> `SIDE_ALIGN=500`、`MAX_GAP=10*1024`、`MIN_READ=900`；此外还有
> `MERGE_DIST=500`、`MAX_CHROMOSOME_SIZE=300_000_000`、`INDEX_CUTOFF=0.001`
> （search/decompose 共用，按 0.1% 累积频率过滤高频 k-mer）；`hit.codon` 打分
> `MATCH=5/MISMATCH=-4/GAPO=-40/GAPE=-1`。

search 按（染色体, chunk）逐段处理：Python 侧按 `--max-chromosome-size` 将每条染色体切成若干段，
每段一个 search job（默认跳过含 `_` 的 contig 与 `chrM`，`--keep-contigs` 关闭该过滤）。对目标
染色体先建前向 + 反向互补两个拷贝的索引（写入同一个 k-mer 索引），长度截断到
`MAX_CHROMOSOME_SIZE + MAX_SD_LEN`；反向互补拷贝作为 query 扫描时只接受来自前向拷贝
（`chr == current.chr - 1`）的 loci，避免 rc↔rc 重复配对。同一染色体超出截断的部分以及排在
`-c chr` 之后的染色体只作为 query 侧扫描，不再扩索引。

`build_index()` 同时承担两种角色：

- `build_index=True`: 扫描参考基因组，建立 k-mer → 位置列表的索引。
- `find_sds=True`: 用已建索引对 query 序列做 plane-sweep，更新 `ListNode` 链表并输出 hits。

**2-bit 编码**: `build_index` 遍历 Codon `bio.seq` 的 `Seq`（元素为 2-bit 数值
`A=0, C=1, G=2, T=3`），`(int(si) & 3)` 是防御性掩码。滚动哈希
`h = ((h << 2) | base) & ((1 << 28) - 1)`（`KMER_SIZE=14`）生成 28-bit k-mer hash。

**Winnowing 实现**: `build_index()` 用单调队列维护 `(hash, pos)`，但实现与论文的“滑动窗口最小值”
语义存在偏差（疑似 `winnow[-1][1]` 应为 `winnow[0][1]`）：

- 新 k-mer 入队前，从队尾弹出 hash 不小于新 hash 的元素，使队列 hash 严格递增、队首为最小。
- “窗口过期”循环 `while winnow and winnow[-1][1] < (i - KMER_SIZE + 1) - WINNOW_SIZE:
  winnow.pop(0)` 检查的是**队尾**位置却弹出**队首**；条件只依赖不随弹首变化的队尾，因此一旦触发会
  一直弹到队列清空（模拟实测 100% 以空队列结束）。队首因此不是真正的窗口最小值，而是自上次清空以来
  的 running minimum（记录最小值）；触发条件等价于队尾（最右的右到左最小值）落后当前 k-mer 超过
  `WINNOW_SIZE` 个位置。
- 后果：指纹密度约为论文期望值 `2|G|/(w+1)` 的一半（随机/GC41 序列实测约 `0.056|G|` vs
  `0.118|G|`，`w=16`），且“每个窗口内必有指纹”的保证被破坏（相邻指纹间距实测最大 94，约 50% 的
  间距 > `w+1`）。
- 窗口填满 (`i - KMER_SIZE + 1 >= WINNOW_SIZE`) 后，队首 hash 即为一个 fingerprint；只有
  fingerprint 变化时才进入 plane-sweep/索引插入，避免同一 hash 连续处理（首个 k-mer 因
  `len(index) == 0` 特判总是被处理）。

**Plane-sweep 链表 `update_list()`**: 维护按 `(chr, first)` 排序的 `ListNode` 链表， 每个节点保存：

- reference（索引侧）区间 `(first, last)` 与 query（扫描侧）区间 `(ref, ref_last)`；
- 已扫描步数 `age`、命中次数 `count`；
- 是否曾经满足阈值 `potentional`。

对当前 query 位置 `current.loc` 对应的 reference 位置列表 `loci`（已按 `(chr, loc)` 排序），
`update_list()` 依次处理三种情况：

1. **延伸**: 若 `loci[lidx]` 与某 walker 同属一条 chromosome，且
  `walker.ref != current.loc`、距离在 `MAX_DISTANCE`
  内 (`loci.loc - MAX_DISTANCE ≤ walker.last < loci.loc`)，则延伸 walker 的 `last`/`ref_last`
   并增加 `count`。
2. **插入**: 若 `loci[lidx]` 落在当前 walker 与下一个 walker 之间，则插入新节点， 年龄初始为 0
   （本轮结束时再 `age += 1`；若新位置恰落在当前 walker 自身区间内，walker 的 `age` 还会额外 +1）。
3. **老化**: 若当前 walker 未被延伸或插入，则 `age += 1`，并检查 `count ≥ ceil(age · τ)`。
    - 若满足且长度 `< MAX_SD_LEN`，标记 `potentional`。
    - 若不满足或长度超限，且此前为 `potentional`，并满足 `last - first > QUERY_THRESHOLD` 与
      `current.loc - walker.ref ≥ REF_THRESHOLD`，则调用 `save_sd()` 输出该 hit。

另外，`update_list()` 遍历前会把 `loci` 中位于首节点之前的元素从链表头部逐个插入，遍历结束后把剩余
`loci` 追加到链表尾部；`build_index()` 在 `find_sds=True` 扫描结束时，对 `L` 中剩余满足
`potentional` 且 (`last - first > QUERY_THRESHOLD` 或 `count >= 4`) 的节点调用 `save_sd()` 兜底输出。

**Tau 计算**: `tau()` 先算 `ratio = (MAX_ERROR - MAX_EDIT_ERROR) / MAX_EDIT_ERROR`（=1），
`gap_error = min(1.0, ratio * MAX_EDIT_ERROR)`（=0.15），再返回
`((1 - gap_error) / (1 + gap_error)) * (1 / (2 * exp(KMER_SIZE * MAX_EDIT_ERROR) - 1))`。
此即论文中的 Jaccard 下界。

**输出 `save_sd()`**: 输出前对边界做 `MAX_EXTEND` 填充（同 chromosome 时避免两个 mate 重叠），
并过滤掉完全自重叠（same species/name/strand 且坐标相同）的情况。最终生成 `Hit` 存入按
`(species1, chr1, species2, chr2, complement)` 分组的 `result` 字典。

**`cross_search()`**: 以基因组 A 的单个目标染色体 `ch` 为参考（`find_sds=False`，只建索引），对其建
前向 + 反向互补两个拷贝的索引；过滤高频 k-mer 后，读取基因组 B 的 `<sp>.bed.regions.txt`（每行
`(contig, st, ed)` 区间），把 B 每条 contig 追加为一个 query 染色体，再对其中每个 `(st, ed)` 区间
`i.seq[st:ed]` 调 `build_index(find_sds=True, st=st)` 做 plane-sweep（`st` 作为 `Locus` 的坐标偏移）。
与单基因组 `search` 一致，反向互补 query 只接受 `chr == current.chr - 1` 的 loci。

**search 输出的 merge**: `__init__.codon` 对每组的 hits 调 `hit.Hit.merge(r, MERGE_DIST=500)`
再写盘（`hit.codon` 的 `merge`，与 §3.3.4 decompose 阶段的 `merge()` 是两处独立实现）。流程：先按
`h.x.start` 排序；当 `ph.x.end + dist < h.x.start`（x 轴间隔超过 `dist`）时开启新组并 flush 旧组，
否则把当前 hit 并入当前组——组内维护按 `y.end` 二分排序的 `y_sorted/y_hits`，从 `h.y.start - dist`
起找候选，对每个满足
`x.end+dist ≥ h.x.start && y.end+dist ≥ h.y.start && y.start ≤ h.y.end+dist`
（x/y 双轴都在 `dist` 内接近）的候选取并集合并边界并删除之，随后把当前 hit 按其 `y.end` 插回。
flush 时经 `_y()`：对 `x == y` 的自匹配 hit，若 `span() > 1000` 则复制一份并把 `y.start += 1000`
（把第二个 mate 起点推后，避免与自身完全重叠），其余 hit 原样输出。

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
4. 按 DP 得分排序并重构链，保留满足 `span >= MIN_UPPERCASE_MATCH or span >=
   int(MIN_READ_SIZE * (1 - MAX_ERROR))` 的链；`int(700 * 0.7) = 490 > 90`，OR 条件等价于
   `span >= 490`，实际起约束作用的是 490 bp。

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

- anchors 先排序：`Hit.__lt__` 按 `(x.chr, y.chr, x.start, x.end, y.start, y.end)`（Interval
  dataclass 的字典序，源码中的注释 lambda 是旧键，未生效）。
- 每个 anchor 的自身得分 `score = MATCH * matches - MISMATCH * mismatches - GAP * indel_bp`。
- DP 从 `ai-1` 向前回溯，且 `it` 只对通过上述 interval 过滤的候选计数，最多处理 `max_iter=500`
  个合法前驱；对每对 `(aj, ai)` 计算：
    - `mi = min(c_xs - p.x.end, c_ys - p.y.end)`
    - `ma = max(c_xs - p.x.end, c_ys - p.y.end)`
    - 若 `ma ≥ MAX_GAP` 则跳过。
    - 同染色体时另有双重过滤：(a) 单个 anchor 的 x/y 区间若几乎重叠
      （`max(span_x, span_y) - max(qo, 0) < SIDE_ALIGN`，即贴近对角线的自匹配）直接跳过该 anchor；
      (b) 若两个 anchor 在两个轴上都不相交（`max(0, min(c_xs, c_ys) - max(p.x.end, p.y.end)) ≥ 1`）
      则跳过该前驱。
    - 更新 `dp[ai] = max(dp[ai], dp[aj] + score[ai] - MISMATCH * mi - GAPOPEN - GAP * (ma - mi))`。
- 按 `dp` 得分降序回溯，得到不相交的 anchor 链；候选链先过滤 `est_size < MIN_READ - SIDE_ALIGN`
  （400），再丢弃相对已有结果在两个轴上都被覆盖到 `SIDE_ALIGN` 以内的链；最终链 span 须
  `≥ MIN_READ`（900）。
- 对最终链，用 `Hit(orig_h, [anchors], SIDE_ALIGN)` 拼接 anchors：
    - 重叠的 anchor 先经 `extend()` 合并，其余相邻 anchors 之间的 gap 调用 `align_gap()`：
      令 `mi = min(Δx, Δy)`（较短侧长度），若 `max(Δx, Δy) ≤ 1000` 直接用 `bio.seq.align()`
      精确比对；否则两端各比对 `mi` bp 并取分高者，中间 `max(Δx, Δy) - mi` 用 `I`/`D` 表示。
    - 在链两端各取 `SIDE_ALIGN=500` 做 `ltrim()` / `rtrim()`：从端点向内侧扫描，找到累积比对得分
      最大的位置作为新边界（实际向两侧取的宽度为 `off = min(SIDE_ALIGN, x.start, y.start)`，端点
      靠近 0 时取更小值）。

`hit.codon` 中的 `Hit.align()` 按 `MAX_ALIGN=60 kb` 分块调用 `bio.seq.align()`，并将返回的 `M`
操作细分为 `=`（匹配）与 `X`（错配），CIGAR 操作仅使用 `=`, `X`, `I`, `D`。

#### 3.3.4 Decomposition in `decompose.codon`

分解阶段读取 `cluster` 输出的每个 cluster FASTA，对每个序列扫描 `k=10` 的 k-mer。

**索引构建**: 对 cluster 中所有序列（来自 `cluster.codon` 输出的 `.fa`）建立完整 k-mer 索引，
`kmer.as_int()` 作为 key，记录每个 k-mer 出现的 `(chr_id, loc)` 列表。同样用频率阈值过滤最高的 0.1%
k-mer。`find_elems()` 扫描时跳过频率 ≥ threshold 的 k-mer，遇到已 `visited` 位置时先 flush 当前
链表 `L` 再继续。

**`update_list()` 与 search 的差异**: search 阶段每个节点只跟踪一段 putative SD 的两个 mate；
decompose 阶段需要跟踪同一 elementary SD 在多个序列上的多个拷贝，因此每个 `ListNode` 额外维护：

- `mappings: Dict[int, int]`：记录该节点对应 elementary SD 在每个 chromosome 上的当前最右边界；
- `gap`：自上次命中以来的未命中步数；
- `score`：命中计数。

处理当前 k-mer 的位置列表 `index`（按 `(chr, loc)` 升序构建，`update_list` 自列表末尾向前遍历）时，
对每个节点依次判分支：

1. **延伸**: 若 `n.chr == index[i][0]` 且 `index[i][1] - diff ≤ n.end < index[i][1]`（位置在当前
   节点右侧 `diff` 内），则延伸——更新 `end`、`score += 1`、`gap=0`、`count += 1`；若本次延伸使
   节点长度首次跨过 `MIN_MATE_LEN` 且 `mappings` 中已跟踪多个节点（该 elementary SD 有多个拷贝），
   先触发一次 `process(L, visited, [], n.mappings)` 输出各拷贝已覆盖的前缀区间。
2. **被覆盖**: 若 `n.chr == index[i][0]` 且 `n.end + diff > index[i][1] ≥ n.begin`（位置已落在
   当前节点区间内），只 `score += 1` 并跳过该位置（不延伸、不计数）。
3. **插入**: 若位置不属于当前节点、且应插在当前节点之后（无后继，或后继同 chr 且其
   `end + diff < index[i][1]`），则 `insert_after` 新节点并赋 `score=1`、继承 `mappings`；对尚未
   进入链表的新位置则在链表头部插入（开头的 while 与末尾的 while）并同样继承当前 `mappings`。
4. **老化**: 每轮 `gap += 1, age += 1`；若 `gap < diff` 置 `potentional=True`，否则若此前为
   `potentional` 且长度超过 `MIN_MATE_LEN`，调用 `process()` 输出并清空整个链表（`L = None`）：
    - 若 `mappings` 为空，输出 `(chr, begin, end, score)`；
    - 否则对每个 chromosome 输出从 `begin` 到 `mappings[chr]` 的区间，并更新 `begin` 为
      `mappings[chr]+1`。
5. 节点长度超过 `MIN_MATE_LEN` 时把 `mappings[id] = end`，记录各拷贝的当前最右边界。
6. 将已输出区域对应的 `visited` 位置置为已访问（`process()` 中 `visited[c][i] = True`），避免同一
   碱基被多次分解。注意 `process()` 仅在输出区域数 ≠ 1 时才标记 `visited` 并返回结果；若恰好只
   输出 1 个区域，则返回空列表，调用方不会记录该 elementary SD（该区域也不标记 visited）。

**合并 `merge()`**: 对相邻的两个 elementary SD 集合 `e1/e2`，若区域数相同、逐区域 chromosome 一致、
且 `e1[i]` 的 end 恰落在 `e2[i]` 的 start 之前 `500 bp` 内（`e1[i].end ≤ e2[i].start` 且
`e1[i].end + gap ≥ e2[i].start`，`gap=500`），则合并为一个集合：每个区域取 `(chr, e1.start, e2.end,
score 求和)`；否则各自独立输出。

**输出**: 每个 elementary SD 集合输出为 BED 行，格式为
`species\tchrom\tbegin\tend\tset_id\tlength\tscore\tstrand`；其中 strand 来自序列名中的 `+`/`-`。

> **strand 编码细节**: cluster.codon 写出的 FASTA 头为 `>{sp}#{chr}{strand}#{st}#{ed}`，其中
> `strand` 取 `'+_'[int(reversed)]`——正向为 `+`、反向为 `_`（`cluster.codon:139` 构造 `c2`
> 处）。decompose 按 `#` 拆分后取 `chr[-1]` 作为 strand，凡非 `+` 均按反向处理并翻转坐标
> （`decompose.codon:265-270`）。因此 decompose 的输入头、输出 strand 与 search 的 BEDPE 链向
> (`+`/`-`) 是两套独立的表示，仅在 cross 阶段由 cluster 头统一桥接。

#### 3.3.5 Core duplicon cover in `cover.py`

`cover(bed, elems)` 读取 `final.bed` 与 `.elem.txt`，为每个 `(species, chrom, strand)` 键用 `ncls`
建关于 elementary SD 坐标的区间树（`start, end, elem_id`）。对每条 SD 行，取其两个 mate 的
`(species, chrom, strand, start, end)`，用 `find_overlap` 找到覆盖该 mate 的所有 elementary SD，维护：
- `elem_sd[elem_id]`：该 elementary SD 覆盖到的 SD 行索引集合；
- `sd_elem[sd_li]`：覆盖该 SD 的所有 elementary SD id 集合。

随后 `greedy_set_cover(subsets=elem_sd, parent_set=sd_elem)` 把**所有 SD 行索引**作为被覆盖全集，用
最大堆（按 `max - len(s - result_set)` 打分）每次挑选能覆盖最多*尚未覆盖* SD 的 elementary SD，贪心
近似"覆盖所有 SD 的最小 elementary 集合"。凡覆盖 SD 数 `< 2` 的集合被跳过（`len(sds) < 2`），剩余
core 在 `.elem` 输出中给对应行（`len(l) <= 8`）追加 `CORE` 标记。

#### 3.3.6 Python 管线编排（`__main__.py`）

Codon 二进制只做单条序列对/单 cluster 的计算，跨染色体/chunk、跨物种的调度全在 Python 层：

- **align 分桶**：`align()`/`cross_biser()` 把每条染色体 search 出的 hits 按 `span`
  （= 两 mate 区间较大者）排序，再 **round-robin**（`i % nbuckets`）分到
  `nbuckets = max(threads*2, 50)` 个 bucket，每个 bucket 单独跑一次 `align`（`__main__.py:142-185`），
  避免单次 align 输入过大；cross_align 的 nbuckets 恒为 `max(threads*2, 50)`（`:261`）。
- **任务缓存 / 断点续跑**：`run_biser` 以 `(exe, args)` 的 **md5** 在 `{tmp}/status/<subcmd>_<md5>`
  写缓存文件，已成功的子命令直接跳过（`__main__.py:47-76`）。配合 `--keep-temp` 保留的临时目录，
  `--resume` 可复用该目录续跑未完成部分（`--resume` 会自动置 `keep_temp=True`）。
- **CLI 选项**（`__main__.py:354-438`）：除 `--hard/-H`、`--keep-contigs`、`--no-decomposition`、
  `--max-chromosome-size`（默认 300_000_000）外，还有 `--temp/-T`、`--threads/-t`（默认 1）、
  `--output/-o`（必填）、`--keep-temp/-k`、`--resume`、`--max-error`（默认 30，限 `0<v<100`）、
  `--max-edit-error`（默认 15，同上）、`--kmer-size`（14）、`--winnow-size`（16）、
  `--ld-path`、`--gc-heap`（调试用）、`--version/-v`。
- **任务划分**：基因组名 = `Path(path).stem`；按 `.fai` 把每条染色体切成 `--max-chromosome-size`
  一段，每段一个 `search -c chr <start>` job；默认跳过含 `_` 的 contig 与 `chrM`
  （`valid_chr`，`--keep-contigs` 关闭过滤）。
- **输出路径**：`!hard` 时对 `final.bed` 跑 `translate` 映射回原基因组坐标；`hard` 时直接
  `shutil.copy` final.bed（及 `.elem.txt`）。进程统一设 `env OMP_NUM_THREADS=1` 限制 Codon 的
  OpenMP；`--gc-heap` 设 `GC_INITIAL_HEAP_SIZE`。

#### 3.3.7 Clustering in `cluster.codon`

`cluster(genomes, bed, output)` 把重叠的 SD 按"区间染色"归入若干 group，每个 group 抽成一个
`<output>/<i>.fa`，供 decompose 逐 cluster 处理（`cluster.codon:53-165`）：

1. 把每条 SD 的两个 mate 端点投影到"点集合"：收集所有 `get_first()`/`get_last()` 的 `(chr, pos)`
   并排序去重得 `help_all_dots`（行 82-84），`colors` 数组记录每个点的当前颜色。
2. 对每条 SD，取其两 mate 的起始/结束在 `colors` 中的区间下标，取并集得 `colors_avail`（行 99）。
   - 若全是 0（未染色）：给两段区间赋一个新颜色 `color`，并在 `colors_dict` 记 `color→color`
     （行 100-105）。
   - 否则：取 `colors_avail` 中最小的非零色 `l[0]` 为代表色，把 `colors_dict` 中所有出现过的颜色
     映射到 `l[0]`（近似 union-find 的路径压缩，行 106-113），再给两段区间赋 `l[0]`。
3. 染色完成后，逐 SD 沿 `colors_dict` 做并查集 `parent` 上溯到根（行 120-128），把共享同一根色的
   SD 归入同一 group；按 group 写 FASTA（行 133-165），每行 `{sp}#{chr}{strand}#{st}#{ed}`，
   反向 mate 取反向互补序列。

> 该染色法用"离散点着色"代替区间树：把每个 mate 端点映射到排序后的唯一坐标下标，用数组区间赋值
> 实现 O(区间跨度) 的区间覆盖染色，`while colors_dict[parent] != parent` 上溯即并查集 find（无
> 路径压缩，但整体仍接近线性）。pgr 的 `sd cluster`（`src/libs/sd/cluster.rs:43` `cluster_paf()`）
> 直接复刻了这一染色聚类逻辑。

## 4. 相关研究：T2T-CHM13 SD 注释来源与下游分析

本节记录 T2T-CHM13 基因组中 SD 注释的产生方式，以及基于这些注释开展的 SNV / IGC 研究。

### 4.1 T2T-CHM13 SD 注释的算法来源（Vollger et al. 2022）

> Vollger M. R. et al. *Segmental duplications and their variation in a complete human genome* . Science.
> 2022;376:eabj6965.[https://doi.org/10.1126/science.abj6965](https://doi.org/10.1126/science.abj6965)

T2T-CHM13 v1.1 的 SD 注释直接继承自 Vollger et al. 2022 对完整 T2T 基因组所做的 SD 注释。
该注释的生成流程如下：

1. **输入**: T2T-CHM13 v1.0 组装，并拼接上 GRCh38 的 chrY（用于包含 Y 染色体 SD 信息）。
2. **屏蔽重复序列**: 使用 **TRF** 与 **RepeatMasker** 对组装进行 masking，
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

> **pgr 集成现状**：pgr 的 `sd` 子命令族（`src/cmd_pgr/sd/`，见 `docs/sd.md`）已按 BISER
> 的阶段结构落地为 `search → align → cluster → decompose → cover`（外加 `cross` 与 `run`）。
> 其中 `cluster`/`decompose`/`cover` 直接复刻了 BISER 的算法与关键常量——区间重叠 union-find
> 聚类、k=10 + 50 bp gap + 100 bp 下限的 k-mer 链式分解、贪心 set-cover 标 `CORE`；而 `search`
> 阶段用 pgr 自研 `pgi` 种子链引擎（原生，默认 k=40 + syncmer 8/5）或外部 `lastz`
> （`--engine lastz`）替代了 BISER 的 ordered Jaccard + plane-sweep。下文按“对 pgr 的借鉴点”组织。

1. **SD 同一性标准采用 T2T-CHM13 定义**: PGR 采用 > 1 kbp、> 90% 同一性的 SD 定义（见 4.2.1），
   不追求 BISER 默认的 75% 同一性（30% 错误率）场景。BISER 的 ordered Jaccard + plane-sweep
   虽然支持低同源性检测，但该能力不在 PGR 当前需求范围内；若未来确需检测古老重复，可再参考其
   colinear k-mer matching + 扫描线设计。
2. **Winnowing 作为采样手段**: BISER 用 winnowing 将 k-mer 采样率降到约 `2/(w+1)`，
   并声称同一 windows 内必然命中（注：代码实测指纹密度约 `2|G|/(w+1)` 的一半且该保证被破坏，见
   §3.3.1）。这与 pgr sketching 的有损采样不同，
   但为“在保敏感度的前提下降低索引规模”提供了可借鉴的参数化方法。pgr 现已落地 closed syncmer
   （`src/libs/syncmer.rs`）：同属有界间隔采样（连续点距 ≤ `2(w-1)`、密度约 `2/(w+1)`），且
   canonical 哈希使其链向对称，Jaccard/Mash 距离偏差小于 minimizer （Edgar 2021）。这一替代在
   pgr 的 `pgi` search 引擎中已实际生效（默认 k=40 + syncmer 8/5 作种子采样），从而规避了
   BISER winnowing 实现密度减半/覆盖不完整的缺陷。
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
7. **软屏蔽→硬屏蔽与坐标回翻译**: BISER 要求 soft-masked 输入，内部先 `mask()` 滤掉小写碱基生成
   hard-masked 基因组，检测/分解全程在硬屏蔽序列上做，最后 `translate()`（`mask.codon:32`）把坐标
   与 CIGAR 映射回原基因组：将 CIGAR 的 `M` 按 lowercase 区间重新切分为 `M`/`S`/`N`，`I`→`I`/`S`、
   `D`→`D`/`N`（用 `uppers`/`lowers` 两张表做坐标偏移 `conv`，行 33-49）。这与 pgr `rept`
   （`src/cmd_pgr/rept/masker.rs`）"屏蔽后分析、再回映射"是同类的工程模式；若 pgr 需在 soft-masked
   输入上直接跑 SD 检测，可复用其 lowercase 区间映射思路。
8. **高频 k-mer 频率阈值过滤**: search 与 decompose 共用 `INDEX_CUTOFF=0.001` 的"累积频率"过滤
   （`search.codon` `filter_index` 行 247-262、`decompose.codon` 行 243-254）：先统计每个 k-mer
   出现次数分布，按出现次数降序累加，直到累计占比 ≤ 0.1% 即把该频次作为阈值，更高频的 k-mer 在
   扫描时跳过（search 侧 `max_frequency=threshold`，decompose 侧 `len(index[k]) >= threshold`）。
   这一"按累积频率而非绝对频次"的自适应阈值，对 pgr 的重复/低复杂度 k-mer 过滤有直接借鉴价值。
9. **染色体分块 + round-robin 分桶约束单任务规模**: search 按 `--max-chromosome-size` 分块，
   align/cross_align 把 hits 按 span 排序后 round-robin 分到 `max(threads*2, 50)` 个 bucket
   （`__main__.py:142-185`、`:261`），每个 bucket 独立跑一次 align。这与 pgr `libs/sd/`
   （`src/libs/sd/cluster.rs`、`decompose.rs`）按染色体的任务划分一致，是控制内存峰值/单核负载的
   实用工程模式。

## 6. 参考文献

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
