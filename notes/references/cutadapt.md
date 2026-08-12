# cutadapt: Mott/BWA 式质量修剪与接头去除

> 整理于 2026-08，源自对 `cutadapt-main/`（v5.2, 2025-10-23）源码的分析。最初承接 [sickle.md](./sickle.md) 中"现代质量修剪算法对比"的结论——滑窗仍是默认，cutadapt 的 `-q` 提供了**唯一真正不同的算法方向**（Mott/BWA 累积质量法），故初稿聚焦**按质量分数修剪**，为 pgr 的 `fq trim-qual` 提供算法与 CLI 参考。后按需求补充了**接头（adapter）去除**的完整算法（indel-aware 半全局比对、多适配器匹配、修饰器与配对端处理），本文现同时覆盖质量修剪与接头去除两条主线。pgr 当前主要关心质量修剪，接头部分为算法参考。

## 1. 简介

`cutadapt`（Marcel Martin, 2011）是 FASTQ 预处理的事实标准工具，以**接头/引物去除**见长，同时提供质量修剪、长度过滤、poly-A 修剪、expected-errors 过滤等。

- **质量修剪入口**：CLI 的 `-q, --quality-cutoff [5'CUTOFF,]3'CUTOFF`（见 [cli.py](file:///home/wangq/Scripts/pgr/cutadapt-main/src/cutadapt/cli.py#L268)）。
- **核心实现**：`quality_trim_index`（[qualtrim.pyx](file:///home/wangq/Scripts/pgr/cutadapt-main/src/cutadapt/qualtrim.pyx#L22)），Cython 加速，单条读段 O(n)。
- **算法来源**：注释明确说明与 **BWA 的 `bwa_trim_read`** 相同（`qualtrim.pyx:29-33`）。这与经典 Mott 算法同源（BWA `-q` 即 Mott 算法），累计和取最小。
- **变体**：`nextseq_trim_index`（NextSeq polyG 暗循环）、`poly_a_trim_index`（poly-A/poly-T）、`expected_errors`（Edgar 2015 期望错误数）。
- **接头去除主线**：见 §3——Cython 半全局比对器 `Aligner`（indel-aware、混合 cost/score）+ k-mer 预过滤 + 多适配器索引，配以 `AdapterCutter` 修饰器与配对端变体。

> **范围说明 / 落地状态**：初稿时 pgr 的 `fq trim-qual` 还是设计提案；现已实现于 `src/cmd_pgr/fq/trim_qual.rs` + `src/libs/fq/trim.rs`，`quality_trim_index`（Mott/BWA 累积质量法）作为 `--method mott`，与继承自 sickle 的滑窗（`--method sliding`）并存。§5 已从"设计提案"改写为"已实现对照 + 剩余借鉴点"。

## 2. 核心算法：`quality_trim_index`（Mott/BWA 累积质量法）

这是 cutadapt 质量修剪的核心，约 50 行 Cython。**5' 端与 3' 端分别独立计算切点**。

### 2.1 算法步骤（`qualtrim.pyx:22-73`）

```
参数：qualities（ASCII 质量串）、cutoff_front、cutoff_back、base（默认 33）

5' 端：
  s = 0; max_qual = 0; start = 0
  for i in 0..n:
    s += cutoff_front - (qual[i] - base)     # cutoff 减实际质量
    if s < 0: break                          # 累积和为负，停止
    if s > max_qual: max_qual = s; start = i+1

3' 端（反向）：
  s = 0; max_qual = 0; stop = n
  for i in reversed(0..n):
    s += cutoff_back - (qual[i] - base)
    if s < 0: break
    if s > max_qual: max_qual = s; stop = i

if start >= stop: start, stop = 0, 0         # 无有效区间
return (start, stop)
```

**关键数学直觉**：定义 `score(i) = sum_{j<i} (cutoff - q_j)`（对 5' 端）。当累积和**首次由正变负**（`s < 0`）时停止；切点选在累积和**达到最大**（`max_qual`）之后的位置。等价于 BWA/Mott：在"前面质量足够好"的连续段内，找到累积质量相对 cutoff 的**局部最大点**作为修剪边界。

- **5' 端**：从 5' 到 3' 扫描，`s` 是 `cutoff - q` 的累积。`s < 0` 表示从此处开始质量持续低于 cutoff，应停止 5' 修剪；`s > max_qual` 时更新 `start = i+1`，即质量相对高点的下一个碱基。
- **3' 端**：从 3' 到 5' 反向扫描，对称逻辑，`stop = i`。
- **区间判定**：`start >= stop` 时整段无效，返回 `(0,0)`（即不修剪，或全部丢弃由上层决定——注意这里**不是**丢弃整条，而是置空区间，`QualityTrimmer` 会返回空 read）。

### 2.2 与滑窗（sickle/Trimmomatic）的本质区别

| 维度 | 滑窗（sickle） | Mott/BWA（cutadapt `-q`） |
| :--- | :--- | :--- |
| 判定单位 | 固定窗口内**平均**质量 | 逐碱基累积，**每碱基**都参与 |
| 切点 | 窗口首/末首个越阈值碱基 | 累积和局部最大点 |
| 中间低质量区 | 无法处理（只看两端窗口） | **可定位并修剪中间低质量区** |
| 复杂度 | O(n)（差分滑动） | O(n)（单遍累积） |
| 阈值语义 | 平均质量阈值 | 逐碱基 cutoff（`cutoff - q` 的累积） |

**优势**：Mott 能处理"读段中部有局部低质量区"的情况，滑窗只能修两端。**劣势**：语义更抽象（累积和），不如"窗口平均质量"直观。

### 2.3 三个变体（`qualtrim.pyx`）

1. **`nextseq_trim_index`**（`--nextseq-trim`）：NextSeq 双色编码中"暗循环"（无颜色）通常被读成高质量 G，出现在读段 3' 端。算法与 `quality_trim_index` 的 3' 端相同，但把 **G 碱基的质量强制设为 `cutoff - 1`**（`qualtrim.pyx:107-108`），使其不贡献正累积，从而把 polyG 尾巴当低质量去掉。
2. **`poly_a_trim_index`**（`--poly-a`）：poly-A/poly-T 尾巴检测，'A'(或'T') 得 +1，其他碱基 −2，累计 score 最大处为切点。错误率上限 0.2 的校验是**位置相关的**：5' 端（polyT head）为 `errors * 5 <= i+1`、3' 端（polyA tail）为 `errors * 5 <= n-i`（`qualtrim.pyx:147,161`），即错误按"当前已扫到的尾巴长度"而非整条读长计。长度 < 3 的尾巴忽略（`best_index < 3` → 0；`best_index > n-3` → n）。
3. **`expected_errors`**（`--max-ee`）：用 Edgar et al. (2015) 公式从 Phred 质量计算期望错误数 `sum(10^(-Q/10))`，用于按总错误数过滤（非修剪）。C 实现 `expected_errors_from_phreds`（`qualtrim.pyx:15`）。

## 3. adapter 修剪算法：indel-aware 半全局比对

cutadapt 的接头去除核心是一个 Cython 实现的**混合 cost/score 半全局比对器** `Aligner`（[_align.pyx](file:///home/wangq/Scripts/pgr/cutadapt-main/src/cutadapt/_align.pyx#L93)），配合 k-mer 预过滤（`KmerFinder`）与多适配器索引（`AdapterIndex`）。pgr 的接头去除走 BBDuk k-mer 路线（`clean --ref --k`），**未采用**本节算法，此处为算法参考。

### 3.1 类层次：`Adapter` / `Match`（adapters.py）

- **`Adapter`**（抽象）→ **`SingleAdapter`**（抽象）→ 具体类型：`FrontAdapter`(5')、`BackAdapter`(3')、`RightmostFront/BackAdapter`、`AnywhereAdapter`、`NonInternalFront/BackAdapter`、`PrefixAdapter`(锚定 5')、`SuffixAdapter`(锚定 3')，以及 **`LinkedAdapter`**（5'+3' 连接的复合接头）。
- **`Match`**（抽象）→ `SingleMatch` → **`RemoveBeforeMatch`**（剪掉 match 之前的序列，用于 5'/front，保留 `read[rstop:]`）/ **`RemoveAfterMatch`**（剪掉 match 之后的序列，用于 3'/back，保留 `read[:rstart]`）；`LinkedMatch` 组合前后两段。

`SingleAdapter` 关键参数（adapters.py:564）：`sequence`（自动转大写、U→T、I→N）、`max_errors=0.1`、`min_overlap=3`（与接头长度取 min）、`adapter_wildcards=True`、`read_wildcards=False`、`indels=True`。要点：
- **max_errors ≥ 1 时按非 N 碱基数折算成比率**（`max_errors /= len - N`，adapters.py:580-581）。
- **错误率 = 错误数 / 与接头比对的那段长度**（非整条 read），见 _align.pyx:144 的 `errors / (reference_stop - reference_start)`。
- indel 通过 `_make_aligner` 把 `indel_cost` 设为 1（允许）或 100000（`--no-indels` 时，adapters.py:605），后者等价于禁止 indel。
- 匹配类 `SingleMatch` 的 `score` 与 `errors` 分别记录比对得分与错误数（用于多接头竞争与报告）。

### 3.2 比对核心 `Aligner`（_align.pyx）

半全局：允许以零代价跳过 reference（接头）/query（read）的前缀和后缀，由 4 个 bit flag 控制（`EndSkip`：REFERENCE_START=1、QUERY_START=2、REFERENCE_END=4、QUERY_STOP=8；全设即 SEMIGLOBAL=15）。adapter 类型通过 `Where` 枚举映射到 flag 组合（见 §3.3）。

- **双矩阵**：每格同时维护 `cost`（编辑距离，mismatch=1、indel=indel_cost）和 `score`（match=+1、mismatch=−1、insertion/deletion=−2，_align.pyx:16-19）。`cost` 用于错误率约束，`score` 用于最优化选择。
- **错误率上限**：`k = int(max_error_rate * m)`（m=接头长）。候选需满足 `cost <= effective_length * max_error_rate` 且 `length >= min_overlap`，其中 `effective_length` 扣除了 N 通配碱基（_align.pyx:557-560）。
- **最优化判据**（_align.pyx:146-154）：① error_rate 不超上限；② 其中 score 最高；③ score 相同时错误数最少；④ 仍相同时取 read 内最左位置。
- **单列 DP + origin**：内存中仅保留一个 `_Entry` 列（`column`，含 cost/score/origin），`origin` 记录比对起点（负=起点在 reference 内，正=在 query 内），据此回溯得到 `(ref_start, ref_stop, query_start, query_stop)`（_align.pyx:579-587）。
- **Ukkonen 带状剪枝**：`last` 追踪"cost ≤ k 的最远行"，只计算错误带内的格；`start_in_query=0` 时列上限 `min(n, m+k)`，`stop_in_query=0` 时列下限 `max(0, n-m-k)`（_align.pyx:343-352）。
- **结束位置搜索**：若 `stop_in_query`，当整条 reference 比对完（`last==m`）时检查 `column[m]` 作为候选（match 结束于 read 内部）；若 `stop_in_reference`，列扫完后在最后一列倒序搜索（match 结束于 read 末尾）（_align.pyx:494-572）。
- **N 通配计数的前缀数组**：`_set_reference`（_align.pyx:250）预计算 `n_counts[i] = reference[:i] 中 N 的个数`，并缓存 `effective_length = m - N`；当仅匹配 reference 的子段时，用 `length - (n_counts[m] - n_counts[m-length])` 重算该子段的非 N 有效长度（_align.pyx:504-510, 543-549），以此对错误率约束"按非通配长度"折算。用前缀数组把任意子区间 N 计数降到 O(1)。
- **字符表翻译与快速比较**：`_reference` 在初始化时通过 `_match_tables.py` 的查找表（`_upper_table`/`_acgt_table`/`_iupac_table`）一次性转成紧凑字节表示（`translate`，_align.pyx:43）；无通配时走 `compare_ascii` 直接按字节相等比较，有 IUPAC 通配时用**位掩码交集** `(s1[i-1] & s2[j-1]) != 0` 判断兼容（_align.pyx:322-328, 442-445）。用"空间换比较速度"，避免每次逐字符查通配表。
- **精确匹配提前终止**：当候选满足 `cost == 0 且 origin >= 0`（零错误精确命中）时立即 `break` 跳过整条 read 的剩余扫描（_align.pyx:531-533）。
- 返回 `(ref_start, ref_stop, query_start, query_stop, score, errors)`；无满足条件者返回 `None`。

> 配套的纯 Python 模块 `align.py` 提供 `edit_distance`（经典 O(mn) DP）、`hamming_environment`/`naive_edit_environment`/`slow_edit_environment` 等供测试对照的实现，`edit_environment` 的 Cython 版则在 §3.6 的 `AdapterIndex` 索引构建中实际使用。

**`PrefixComparer` / `SuffixComparer`**（_align.pyx:594,696）：当锚定适配器禁用 indel 时使用，不做 DP 矩阵，只逐位统计前缀/后缀错误数，满足 `errors <= max_k` 且 `length >= min_overlap` 即返回固定区间元组（更快；`PrefixAdapter` 在 `indels=False` 时走此路径）。

### 3.3 各适配器类型对应的 Where 标志（adapters.py:39-53）

| 类型 | Where 标志 | 语义 |
| :-- | :-- | :-- |
| `FrontAdapter` (5') | FRONT = Q_START\|Q_STOP\|R_START | 接头对齐读段 5' 端，剪掉接头之前 |
| `BackAdapter` (3') | BACK = Q_START\|Q_STOP\|R_END | 接头对齐读段 3' 端，剪掉接头之后 |
| `AnywhereAdapter` | ANYWHERE = SEMIGLOBAL | 若 `rstart==0` 视为 5'，否则视为 3' |
| `NonInternalFrontAdapter` | FRONT_NOT_INTERNAL = R_START\|Q_STOP | 不允许内部匹配 |
| `NonInternalBackAdapter` | BACK_NOT_INTERNAL = Q_START\|R_END | 不允许内部匹配 |
| `PrefixAdapter` (锚定 5') | PREFIX = Q_STOP | 必须从读段起始对齐（`^`） |
| `SuffixAdapter` (锚定 3') | SUFFIX = Q_START | 必须对齐到读段末尾（`$`） |
| `RightmostFront/BackAdapter` | 反转序列后以 BACK/FRONT 比对 | 偏好最右匹配 |

`AnywhereAdapter.match_to` 里靠 `alignment[2] == 0`（rstart 为 0）判 5' 还是 3'，分别产出 `RemoveBeforeMatch` / `RemoveAfterMatch`（adapters.py:930-934）。

### 3.4 接头规格解析与 CLI 类型映射（parser.py）

`make_adapter` / `AdapterSpecification` 把 `-a/-g/-b` 的字符串解析成具体类。放置限制记号（parser.py:292-319）：
- `^ADAPTER` → `PrefixAdapter`（锚定 5'）
- `ADAPTER$` → `SuffixAdapter`（锚定 3'）
- `XADAPTER`（X 前缀）→ `NonInternalFrontAdapter`
- `ADAPTERX`（X 后缀）→ `NonInternalBackAdapter`
- `...` 区分前后接头与 anywhere（linked 由两个子接头组成）
- `;rightmost`、`name=...`、`min_overlap=`（`o=`）等参数经 `parse_search_parameters` 解析；`{...}` 花括号展开（`expand_braces`）。
- 同一条规格里**不能**同时出现多个放置限制。

### 3.5 k-mer 预过滤（_kmer_finder.pyx + kmer_heuristic.py）

真正的 DP 比对前先做**可命中性检查** `KmerFinder.kmers_present()`（adapters.py:715 等）：对每条 adapter 生成一组"位置 → k-mer 集合"，要求 read 中至少一个 k-mer 出现在指定位置附近，否则直接判不匹配、跳过昂贵的 DP。`create_positions_and_kmers`（kmer_heuristic.py:118）按错误率把 adapter 分成若干段，每段取 `max_errors+1` 个互不重叠 chunk 作为必需 k-mer（例：AAAAATTTTT 允许 1 错时，AAAAA 或 TTTTT 至少一个必须存在）；front/back 适配器另加部分重叠的短 k-mer（`create_back_overlap_searchsets`）。k-mer 过长无法建索引时退回 `MockKmerFinder`（恒真，退化为纯 DP）。

- **shift-and 位并行多模式匹配**（_kmer_finder.pyx）：`KmerFinder.kmers_present`（:170）把同一搜索位置的一组必需 k-mer **首尾相接拼进单个 64-bit 机器字**，每碱基占 1 位；`init_mask` 标记各 k-mer 的起点位、`found_mask` 标记终点位，再对 read 单遍跑 shift-and（`shift_and_multiple_is_present`，:241）：`R = (R<<1 | init_mask) & needle_mask[base]`，`R & found_mask` 非零即命中。多条 k-mer 同时检测，**O(read_len) 且常数极小**；总长超过 64 位时在 `__cinit__` 内层循环自动拆成多个 bitmask 分组（:131-149）。`needle_mask` 是 128 项查表（按 ASCII 值索引），IUPAC 兼容匹配通过 `_match_tables.matches_lookup` 预合并到位掩码里。对 pgr 以 k-mer 为中心的路线（如 `clean --k` 的 BBDuk 计数、`paf` k-mer 索引）是值得参考的位并行技巧——相比逐 k-mer `find()`，它把多模式匹配常数压到最低。

### 3.6 多适配器匹配：`MultipleAdapters` 与 `AdapterIndex`

- **`MultipleAdapters.match_to`**（adapters.py:1265）：依次对每个 adapter 调 `match_to`，选 **score 最高、score 相同取 error 最少** 者为最佳匹配——这是多接头竞争的基本规则。
- **`AdapterIndex`**（adapters.py:1289）：对**锚定**（Prefix/Suffix）适配器建索引加速。索引构建（`_make_index`，:1396）把每个 adapter 的错误数 ≤k 的**编辑环境**内所有字符串预先展开存入 dict：允许 indel 时用 `edit_environment(seq, k)`（_align.pyx:785，一个按 DP + Ukkonen 带状约束遍历 edit distance ≤k 的字符串、同时计数匹配数的生成器），不允许 indel 时用 `hamming_sphere`（_align.pyx:717，k=1/2 有专门展开、k>2 递归）。查询时直接 O(1) 查 read 的固定后缀/前缀。限制：k≤3、adapter 不允许通配、read 不允许通配、适配器须为锚定类型。多长度时按长→短依次查，用"匹配数不可能超过已找到的 `best_m`"提前终止（:1507）。出现歧义（两个 adapter 对同一字符串同分，`matches` 相等）的字符串会被**删除**——含歧义序列的 read 将**不修剪**（并有日志警告，:1444-1466）。k=3 且有 indel 时索引可能巨大，日志会提示 `--no-indels`/`--no-index`（:1407-1412）。
- **read 中 N 通配的兜底**：限制"不允许通配"指 adapter 序列；read 侧遇 N 时 `_lookup_with_n`（adapters.py:1535）先把 N 替换成 A 查索引，命中后**对该 affix 重跑一次真比对**（`adapter.match_to`）修正 errors/score——避免把 N 误算成匹配。
- `AdapterCutter._regroup_into_indexed_adapters`（modifiers.py:127）把用户适配器里"可被索引的锚定"拆出来建 `IndexedPrefixAdapters` / `IndexedSuffixAdapters`，其余进 `MultipleAdapters`。

### 3.7 修饰器 `AdapterCutter` 与 action（modifiers.py:82）

`AdapterCutter` 负责反复找接头并按 `action` 处理：
- `times`（`-n/--times`，默认 1；`-e` 是 `--error-rate`）：每轮只剪最佳匹配，剪完再搜下一轮，直到某轮无匹配为止（modifiers.py:225-231）。
- `action`（`--action`）：`trim`（默认，删除接头及上下游）、`mask`（将被删区替换为 N）、`lowercase`（转小写）、`retain`（保留接头本体，只删其上下游）、`crop`（只保留接头匹配区）、`none`（`--no-trim`，只记录不删除）。`retain`/`crop` 不能与 `times>1` 组合。
- 统计：`with_adapters`、每个 adapter 的 `EndStatistics`（按 removed 长度 × 错误数计数、`adjacent_bases`）供报告。

### 3.8 paired-end 处理（modifiers.py + cli.py）

- **默认**：R1、R2 各自的 `AdapterCutter` 独立工作（`PairedEndModifierWrapper` 包成 `(modifier1, modifier2)`），接头可只给一端（另一端为 `None`）。
- **`--revcomp`**：`ReverseComplementer`（单端）/ `PairedReverseComplementer`（双端）——同时跑正向与反向互补（双端即 R1/R2 交换）两种方案，取**总 score 更高**者；反向时给 read 名加 `" rc"` 后缀并置 `is_rc` 标记（modifiers.py:278-405）。
- **`--pair-adapters`**：`PairedAdapterCutter`（modifiers.py:412）——R1/R2 的接头必须成对给（两列表长度相同），只有"第 i 个接头在 R1、第 i 个在 R2 都匹配"才修剪，选总 score 最高的接头对；与 `--revcomp` 互斥（cli.py:1086-1087）。
- **过滤步（steps.py + predicates.py）**：修饰完成后按序执行 `SingleEndFilter`/`PairedEndFilter`。配对过滤有 `pair_filter_mode`：`any`（任一命中即丢）、`both`（都命中才丢）、`first`（仅看 R1）。谓词：`TooShort`(`-m`)、`TooLong`(`-M`)、`TooManyExpectedErrors`(`--max-ee`)、`TooHighAverageErrorRate`(`--max-er`)、`TooManyN`(`--max-n`)、`CasavaFiltered`、`IsUntrimmed`/`IsTrimmed`（`--discard-untrimmed`/`--discard-trimmed`）。过滤顺序：rest/info/wildcard 写出器先见所有 read，随后是上述过滤器，最后是 sink。

## 4. 质量修剪的修饰器封装与 CLI

### 4.1 `QualityTrimmer` 修饰器（`modifiers.py:840`）

cutadapt 采用"修饰器（modifier）"流水线架构——每个操作是一个 `SingleEndModifier`，对 `read` 返回 `read[slice]`：

```python
class QualityTrimmer(SingleEndModifier):
    def __call__(self, read, info):
        start, stop = quality_trim_index(read.qualities, self.cutoff_front, self.cutoff_back, self.base)
        self.trimmed_bases += len(read) - (stop - start)
        return read[start:stop]
```

- 统计修剪碱基数（`trimmed_bases`）供报告。
- 返回 `read[start:stop]`，即保留 `[start, stop)` 区间（**左闭右开**）。

### 4.2 CLI 参数（`cli.py:268`）

- `-q, --quality-cutoff [5'CUTOFF,]3'CUTOFF`：可指定单个（默认只修 3' 端）或 `5,3` 两个值。注：cutoff 为 `0` 时该端禁用修剪（`cli.py:1059` 的 `cutoff != "0"` 判断）。
- `--quality-base N`：Phred 偏移，默认 33（Sanger）。
- `-Q`：R2 的独立 cutoff（配对，默认继承 R1）。
- 双端时 `-q 5` 未给 R2 则 R2 复制 R1 的修剪器（`cli.py:1065`）。

**`parse_cutoffs`**（`cli.py:419`）解析 `"5"` → `(0,5)`（只修 3'），`"6,7"` → `(6,7)`。

### 4.3 执行顺序（pipeline）

cutadapt 的 read 修饰器按**固定顺序**依次对 read 生效，`make_pipeline_from_args` 依序 append（`cli.py:937-980`）。实际顺序是：

1. `--cut` 无条件切头尾（`UnconditionalCutter`）
2. `--nextseq-trim`（NextSeq polyG 修剪）
3. `-q/--quality-cutoff` 质量修剪（`QualityTrimmer`）
4. 去接头（`-a/-g/-b`，`AdapterCutter`）
5. `--poly-a` poly-A/poly-T 修剪
6. `-l/--length` 缩短（`Shortener`）
7. 两端通用修饰：`--trim-n`（`NEndTrimmer`）、`--length-tag`、`--strip-suffix`、`-x/-y` 前后缀、`--zero-cap`（`ZeroCapper`）
8. `--rename` 重命名（`Renamer`/`PairedEndRenamer`）

**注意：质量修剪在去接头之前**（`cli.py:947-967`，`make_quality_trimmers` 先于 `make_adapter_cutter`），与直觉相反——先按质量把低质量区剪掉，再做接头比对，这样接头搜索不受低质量 3' 尾干扰。pgr 移植纯质量修剪时无此依赖，但若日后加 `--method` 组合（如质量修剪 + 去接头）应保持"先质量后接头"的顺序。

## 5. 对 pgr 的启示：`fq trim-qual` 与 `fq clean`

> **落地状态**：pgr 的 `fq trim-qual` 已实现（`src/cmd_pgr/fq/trim_qual.rs` + `src/libs/fq/trim.rs`），本文 §2 的 `quality_trim_index` 以 `Method::Mott` 落地为 `--method mott`。本节从初稿的"设计提案"改写为"已实现对照 + 剩余借鉴点"。

### 5.1 已落地：`mott_cut`（`libs/fq/trim.rs`）

`Method::Mott` 在 `trim_interval`（`libs/fq/trim.rs:263`）中调用 `mott_cut`（`libs/fq/trim.rs:192`），与 Cython 版 `quality_trim_index` 逐行对应，约 30 行，仍是"单遍累积 + 局部最大"：

```rust
fn mott_cut(qual: &[u8], base: u8, cutoff_front: f64, cutoff_back: f64) -> (usize, usize) {
    // 5' end: s += cutoff_front - (q - base); s<0 break; s>max_s => start=i+1
    // 3' end(reverse): 对称，s<0 break; s>max_s => stop=i
    if start >= stop { (0, 0) } else { (start, stop) }
}
```

与 Cython 原版的差异：

- **int → f64**：cutadapt 的 `-q` 是整数 cutoff；pgr 的 `--qual-threshold` 是 `f64`，score 累积用 f64，阈值语义更宽（可给小数）。数值上等价于整数版。
- **单阈值两端共用**：`trim_interval` 把 `--qual-threshold` 同时作为 front 与 back cutoff；`--no-fiveprime` 时 front 置 `0.0` 禁用 5' 修剪。**未实现** cutadapt 的 `5,3` 独立 cutoff（见 §5.3.1）。
- **零 panic 校验**：`process_record` 先调 `validate_quality`（`trim.rs:245`），把质量字符限制在 `[base, base+93]`，越界即 `bail!` 报错——正是初稿 §5.4 的建议，已落实。
- **区间判定一致**：`start >= stop → (0,0)`；`trim_interval` 再按 `--length-threshold` 决定丢弃整条。
- **`--polyg-right`（polyg_end, `trim.rs:227`）**：3' 端数连续 G 跑，达标才剪（见 §5.3.2 与 cutadapt 的差异）。

### 5.2 与滑窗共存（已实现）

CLI 已提供 `--method sliding`（默认）与 `--method mott`（`trim_qual.rs:71-77`），对应 §2.2 的两类算法，验证了 [sickle.md](./sickle.md) §4.5 的结论。两者共用 `--qual-threshold`/`--length-threshold` 与长度过滤。

### 5.3 其余值得借鉴但尚未落地的点

1. **5'/3' 独立 cutoff**：cutadapt `-q 5,3` 允许两端不同阈值；pgr `mott_cut` 目前单阈值（加 `--no-fiveprime` 整体禁 5'）。若需"前端高严格、后端低严格"，可加 `--qual-front/--qual-back` 或 `--quality-cutoff FRONT,BACK`。
2. **cutadapt 的 NextSeq polyG（`nextseq_trim_index`）未移植**：pgr 的 `--polyg-right`（`polyg_end`）是 **BBDuk 式简单实现**——从 3' 端数连续 G 跑，长度达标才剪；`fq clean` 的 `--trim-poly-g-right` 同理。cutadapt 则是把 G 的质量强制设为 `cutoff-1` 后**嵌入累积算法**，能处理"夹着少量错配/低质量非 G 段"的 G 尾。对"纯连续 G 尾"两者等价，对"含噪声的 G 尾"cutadapt 式更鲁棒。若 pgr 面向 NovaSeq/NextSeq 平台，可评估移植 nextseq 式。
3. **poly-A 位置相关错误率（`poly_a_trim_index`）**：cutadapt 按"已扫到的尾巴长度"限制错误率（`errors*5 <= i+1` / `n-i`），比"整条按错误数"更精细；pgr `clean` 的 `--trim-poly-a` 走 BBDuk `trimpolyA`，语义不同。
4. **`--max-ee` 期望错误数过滤未落地为独立选项**：pgr 的 `expected_errors`（`trim_adapter.rs:907`、`clump.rs:607`）是 BBTools 式（`sum(10^(-Q/10))`，`prob[128]` 查表），用于 dedup/clump 排序，**不是** cutadapt 的 Edgar 2015 用户级过滤。长读场景可把 `--max-ee` 作为 `fq filter` 的过滤选项参考。
5. **`--quality-base` 可配置**：pgr 的 `--quality-base`（`33/64/auto`）比 cutadapt 更强——除显式覆盖外，还通过 BBDuk flip-flop 启发式（`detect_quality_base`, `trim.rs:88`）自动探测编码。

### 5.4 与 `fq clean` 的关系（重要澄清）

pgr 目前有**两条质量修剪路径，算法来源不同**：

- `pgr fq trim-qual`：滑窗（sickle）+ Mott（cutadapt 本文算法），纯质量修剪，单/双端，可配 `--outfile-2/--outfile-single`。
- `pgr fq clean`：整体是 **BBDuk（bbduk.sh）** 移植（`libs/fq/trim_adapter.rs`），其 `--qtrim` 取 `r/l/rl/w/f`——这是 **BBDuk 的质量修剪模式**（`r`=右侧逐碱基、`w`=滑窗等，对应 `trimq`/`qtrim-window`），**不是** cutadapt 的 Mott。即 `clean` 的质量修剪不来自 cutadapt。

因此 cutadapt 对 pgr 的算法价值**已集中在 `trim-qual --method mott`**；接头/适配器去除部分 pgr 走 BBDuk k-mer 路线（`clean --ref` + `--k`），**未采用** cutadapt 的 `-a/-g/-b` 适配器匹配（banded 半全局 DP，见 §3）。

### 5.5 边界与实现注意

- `quality_trim_index` 对质量越界**不校验**（Cython 直接 `qual[i] - base`，可能得负值），cutadapt 只在 `expected_errors` 里校验。pgr 的 `validate_quality` 已补齐（零 panic 硬约束）。
- 左闭右开切片 `read[start:stop]` 与 pgr `SeqRecord::sequence()[start..end]` 语义一致。
- cutadapt `-q 0` 禁用修剪（`cli.py:1059` `cutoff != "0"`）；pgr 以 `--no-fiveprime` / 阈值 0 近似，语义可对齐。
- pgr 用 `f64` 阈值而非整数，注意与 cutadapt 精确逐碱基复现时，`f64` 累积与整数累积在极端大 read 上可能有浮点尾差，但常规阈值下等价。

---

*参考来源: [cutadapt GitHub](https://github.com/marcelm/cutadapt) | 本项目源码 `cutadapt-main/`（v5.2） | [sickle.md](./sickle.md)*
