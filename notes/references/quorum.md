# quorum-1.1.2：基于 k-mer 计数的 read 纠错（源码分析）

> 2026-08 整理，纯源码分析（`quorum-1.1.2/`）。quorum 是 Guillaume Marçais
> 的 read 纠错工具，依赖 Jellyfish 2.0 的 k-mer 哈希计数；项目较老（源码
> Copyright 2012，NEWS 仅 "Initial release 0.1.0"，无后续版本记录），但
> **算法基础扎实**，是 pgr 未来做 read 纠错的主要参考。

## 1. 概况

- **定位**：对 Illumina 双端/单端 reads 做**基于 k-mer 计数的纠错**——把
  低计数的错误 k-mer 修正为可信碱基，或截断不可信区域。
- **构建**：autotools（configure/make），依赖 Jellyfish 2.0
  （`pkg-config --exists jellyfish-2.0`）；需 gcc 4.4-4.7 时代代码风格。
- **输入输出**：输入 FASTQ；默认输出 `quorum_corrected.fa`（FASTA），
  `--paired-files` 时输出 `_1.fa`/`_2.fa`。FASTA 头带纠错日志：
  `>1204 86:sub:T-C 91:3_trunc 62:5_trunc`（坐标 0-based）。

## 2. 架构：入口脚本 + 四个工具

`src/quorum.in`（Perl）串起四个二进制（CMake 前时代的命令名约定）：

| 工具 | 作用 |
|---|---|
| `quorum_create_database` | 读 FASTQ，建 k-mer 计数库（`hash_with_quality`），写二进制 `.jf` 文件 |
| `quorum_error_correct_reads` | 读计数库 + FASTQ，逐 read 纠错，输出 FASTA |
| `merge_mate_pairs` | 双端两个文件 → 交错（配对保序） |
| `split_mate_pairs` | 交错 → `_1.fa`/`_2.fa` 两个文件 |
| `query_mer_database` / `histo_mer_database` | 调试工具（`check_PROGRAMS`，make check 构建、不随安装）：查单个 mer 的 (count, quality)、输出计数直方图 |

入口脚本额外做：质量编码自动检测（读前 ~1000 条 read，取最小 quality
char；遇 35/66 特殊 −2，校验须为 33/59/64）、k-mer 长度（默认 24，README
注明上限 31）、哈希大小（`-s` 默认 200M，README 给估式 `(G + k·n) / 0.8`，
估小会报 "Failed: Increase the size parameter"）。

## 3. 核心数据结构（`mer_database.hpp`）

### 3.1 `hash_with_quality`：带质量偏置的 k-mer 计数

关键设计：**每个 k-mer 的计数带 1 位"高质量"标记**，值编码为
`(count << 1) | quality`，存于 `atomic_bits_array`（bits 由 `-b` 指定，
1-63；计数溢出即封顶）。`add(key, quality)` 的更新规则：

```
首次出现：quality=1 → nval=3（count=1, 高质量）
           quality=0 → nval=2（count=1, 低质量）
已有值：  新 quality=1 且旧为低质量 → 重置为 3（提升为高质量）
          新 quality=0 且旧为高质量 → 不计数（低质量证据不污染高质量计数）
          同质量 → nval += 2（即 value 增 2，count 实际增 1，质量位不变）
```

> 精确实现见 `mer_database.hpp` `hash_with_quality::add()`：`if ((nval&1) < quality) nval=3;`
> `else if ((nval>>1)==max_val_ || (nval&1) > quality) return; else nval += 2;`
> 注意 `nval += 2` 是对**编码值**加 2（低 1 位是质量位），等价于 count 增 1、质量位保留。

即**至少一次高质量出现就把该 k-mer 记为高质量，且权重从 3 起步**；低质量
reads 的错误 k-mer 不会膨胀计数。这是 quorum 与朴素计数最大的区别。

**quality 位如何判定**（`create_database.cc` 的 `quality_mer_counter`）：
逐碱基维护 `low_len`/`high_len`，凡碱基质量 ≥ `qual_thresh`（= `min-q-char +
min-quality`，`quorum.in` 里 `-q` 传入）则 `high_len++`、否则清零；当
`low_len >= k`（凑满一个 k-mer）时，`quality = (high_len >= k)`——即**只有
连续 k 个高质量碱基的 k-mer 才标记为高质量**。窗口内任一个低质量碱基都会
清零 `high_len`，从而取消该位置 k-mer 的高质量标记。

> `quorum.in` 建库时把 bits 硬编码为 `-b 7`（计数上限 2^7-1=127，足够覆盖
> 一般覆盖度），并默认 `-s 200M`（Jellyfish 哈希槽位，估小会
> "Failed: Increase the size parameter"）。

哈希满时自动翻倍扩容（`handle_full_ary`，多线程屏障同步，旧表迁移到
2× 新表）。

### 3.2 `database_query`：只读 mmap 查询

纠错阶段把 `.jf` 库 mmap（或整读），`operator[](mer)` 返回
`(count, quality)`；`get_best_alternatives(m, counts[4], ucode, level)`：
把 `m` 第 0 位替换为 A/C/G/T，逐个查询 canonical 计数，返回 4 个计数、
最高质量等级（level）与命中数——**纠错的核心查询原语**。

## 4. 纠错算法（`error_correct_reads.cc`）

逐 read 处理，**从可信 anchor 向两端扩展**（`extend` 对 forward/backward
对称，`forward_ptr`/`backward_ptr` 模板统一方向语义）：

### 4.1 `find_starting_mer`：找 anchor

- 从 5' 端滑窗（k-mer 窗口），跳过 N；
- 连续 `good` 个 k-mer 计数 ≥ `anchor`（anchor-count）即认为找到可信
  起点；找不到 → 整条丢弃（`--no-discard` 时输出单碱基 N）。

### 4.2 `extend`：逐碱基扩展与纠错

每个位置把新碱基移入 k-mer，然后 `get_best_alternatives` 检查 4 种可能：

| 情况 | 处理 |
|---|---|
| `count==0`（无任何延续） | `truncation`，截断该端 |
| `count==1`（唯一延续） | 若与当前碱基不同 → `substitution`（替换） |
| `count>1` 且 `counts[ori]` > min-count 且（≥ cutoff 或质量够） | 保留原碱基 |
| `count>1` 且 `counts[ori]` > min-count 但（< cutoff 且质量不足） | **Poisson 碰撞检验**：`p = Σcounts × (先验错误率/3)`，`poisson_term(p, counts[ori]) < 阈值` → 视为随机碰撞、保留；否则落入候选替换 |
| `count>1` 且 `counts[ori]` ≤ min-count，且 `level==0 && counts[ori]==0`（原碱基为错误、无高质量候选） | `truncation` |
| 其余情况（含 N 碱基） | 进入候选替换：对每个计数 > min-count 的候选检查**延续性**（替换后移一位，下一个 k-mer 的 level ≥ 当前）；选计数最接近 `prev_count` 的候选；平局时用 read 下一个碱基仲裁；仍多个候选 → 不纠 |
| 原碱基为 N 且无候选 | `truncation` |

### 4.3 `err_log`：窗口错误数限制（防过度纠错）

所有 sub/trunc 事件按位置记录；**滑动窗口（`-w`，quorum.in 不传时默认
10）内错误数 ≥ `-e`（默认 3）即触发回退截断**（`remove_last_window` 丢弃
窗口内事件并截断到窗口起点；`window()`/`error()` 仅当显式传 0 时才回退到
k 与 k/2）。输出日志：`pos:sub:from-to`、`pos:3_trunc`（3' 端）、
`pos:5_trunc`（5' 端）。

### 4.4 其他

- `compute_poisson_cutoff`：未显式给 `-p` 时，从计数分布自动估计 cutoff
  ——只统计**高质量** mer（编码值奇数且 ≥2），`coverage = total/distinct`
  （总 k-mer / 去重数，即平均计数），`lambda = coverage × 先验错误率/3`，
  取首个满足 `poisson_term < 阈值` 的 x 并**返回 x+1**。
- `homo_trim`（`--homo-trim`）：从 corrected read 末端向前扫描，逐位累计
  homopolymer 评分（`(same<<1)-1`：同碱基 +1、异 −1），记录累计最大值
  位置；**仅当最大评分 ≥ 阈值才在该位置截断**（否则不截）。
- `contaminant`：可选 Jellyfish 污染库（需与主库同 k），命中即丢弃
  （或 `--trim-contaminant` 截断）。

## 5. 参数语义（`quorum.in`）

| 参数 | 含义 |
|---|---|
| `-s` | Jellyfish 哈希槽位大小（必须容下全部 k-mer） |
| `-k` | k-mer 长度（默认 24） |
| `-q` / `-m` | 质量下限字符 / 高质量阈值偏移（高质量 = char ≥ q+m） |
| `-w` / `-e` | 错误窗口大小（默认 10）/ 窗口内最大错误数（默认 3）；显式传 0 时回退到 k / k/2 |
| `--min-count` / `--skip` | 好 k-mer 最小计数 / 找 anchor 前跳过的碱基数 |
| `--anchor` / `--anchor-count` | anchor 连续个数 / anchor k-mer 最小计数 |
| `-p` | cutoff（显式覆盖 Poisson 自动估计） |
| `--contaminant` / `--trim-contaminant` | 污染库 / 命中即截断 |
| `-d` / `-P` / `--homo-trim` | 不丢弃（输出 N）/ 双端分文件 / homopolymer 截断 |

## 6. 与 pgr 的关联

- **pgr 现状**：`libs/kmer`（`KmerTable`：canonical 2-bit u128 key、精确
  计数、radix sort、rayon 并行）是**精确计数**路线；quorum/Jellyfish 是
  **哈希近似**路线（内存可控、自动扩容、带质量偏置）。pgr 目前**没有
  read 纠错功能**，但已有 `pgr fq norm`（`src/libs/fq/norm.rs` +
  `src/cmd_pgr/fq/norm.rs`）：基于精确 KmerTable + minq 高质量过滤的
  **低深度 read 过滤**，正是 §6.1 所述"判定器"的现成雏形。
- **可借鉴的算法点**（若 pgr 做纠错）：
  1. **质量加权计数**（高质量 k-mer 权重 3 起步、低质量不污染）——直接
     可移植到 KmerTable 的计数语义；
  2. **anchor + 双向扩展**的纠错框架（`find_starting_mer` + `extend` +
     `get_best_alternatives`）；
  3. **Poisson 碰撞检验**与**窗口错误数限制**（`err_log`）——防止过度
     纠错的两个关键机制；
  4. 纠错日志输出格式（`pos:sub:X-Y`/`pos:N_trunc`）便于审计。
- **pgr 的差异优势**：精确计数（无哈希碰撞）、SIMD 能力（canonical_keys
  是滚动 2-bit 编码）、`.pkt` 缓存；quorum 的哈希自动扩容/近似计数在
  大内存场景值得对照。

## 6.1 应用解读：anchr 用法 = 错误 read 直接丢弃（2026-08 补充）

用户确认：**根本不想纠错，只要"检测有错误的 reads 并直接丢弃"**。anchr 的
`quorum.tera.sh` 之所以"完整跑纠错再丢弃"，只是因为脚本层无法改 quorum
行为——quorum 只输出修正后的序列，脚本只能靠 header 里的 `:sub:`/`trunc`
标注反推"这条 read 被判定有错"，再用 `hnsm some` 剔除。

因此 pgr 若实现该功能，**不需要修正器，只需要判定器**：

- **需要保留的信号**（quorum 触发 sub/trunc 的条件，即"有错"证据）：
  1. `find_starting_mer` 失败（找不到连续 `good` 个高计数 anchor k-mer）；
  2. `extend` 中 `count==0`（当前 k-mer 无任何碱基延续）；
  3. `count==1` 且唯一延续碱基 ≠ 当前碱基（sub 信号）；
  4. 多候选时 Poisson 碰撞检验判定为错误；
  5. 窗口内错误数超限（`err_log`）。
- **可以砍掉**：替换/截断的输出序列、`err_log` 的坐标日志格式、
  `homo_trim`、修正结果的 FASTA 输出。
- **输出的形态**：保留判定通过的 reads（原样 FASTQ），丢弃有错 reads
  （可输出 discard 名单或直接过滤）。
- 简化空间：只做"read 内是否存在低计数且无延续的 k-mer 区域"判定，可能
  比完整 anchor+extend 更简单；但 anchor+Poisson+窗口限制的判定质量是
  用户满意的部分，移植时应保留其判定语义。
- **pgr 现成落点**：`pgr fq norm` 已实现"按 k-mer 深度判定并过滤 read"
  同类能力（精确计数表 + minq 高质量过滤 + 按 min-depth 丢 read），与本
  场景高度契合；若需保留 quorum 的 anchor+extend+Poisson+窗口限制判定
  语义，可在此基础上扩展而非重写。

## 7. 局限

- 依赖 Jellyfish 2.0（外部 C++ 库），构建链 autotools + yaggo，较旧；
- 输出是 FASTA（丢质量），且纠错后 read 可能变短/截断；
- 单机内存模型（`-s` 手动指定哈希大小，估小会失败）；
- 无更新（2012-2014 年项目）。
- **source quirk**：`quorum.in` 在 `--paired-files` 模式会**强制置
  `--no-discard=1`**（`$opts{"no-discard"} = 1 if $paired_files`），即双端模式
  从不会丢 read、错误 read 一律输出单碱基 N——这是为了保持 `_1.fa`/`_2.fa`
  的 mate 配对结构。若 pgr 要做"判定并丢弃"，单端模式才是默认丢弃语义。

---

*参考来源: 本项目源码 `quorum-1.1.2/`（src/ 全部 + README + quorum.in）*
