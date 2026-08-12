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
| `merge_mate_pairs` | 读**偶数个**文件，把偶/奇索引文件按位交错（配对保序），写 FASTQ 到 stdout；文件个数为**奇数**报 "Must give a even number files"、两流配对数不一致报 "Input files are not paired reads."（`merge_mate_pairs.cc:65-84`） |
| `split_mate_pairs` | 从 stdin 读 FASTA，把相邻两行（`>header`+序列）交替写到 `<prefix>_1.fa`/`<prefix>_2.fa` |
| `query_mer_database` / `histo_mer_database` | 调试工具（`check_PROGRAMS`，make check 构建、不随安装）：`query` 查单个 mer 的 (count, quality)，`histo` 输出 (count, 高质量/低质量) 双通道计数直方图 |

`Makefile.am` 里 4 个 `bin_PROGRAMS`；`all_tests`（`unit_tests/test_mer_database.cc`，gtest
参数化测试 bits 与计数语义）、`query_mer_database`、`histo_mer_database` 归 `check_PROGRAMS`。
`data/adapter.jf` 由 Makefile 规则 `.fa.jf` 用 `jellyfish count -m 24 -s 5k -C` 从适配器
FASTA 生成（`dist_data_DATA`）。

**双端模式是三条管道的流水线**（`quorum.in:169-231`）：`--paired-files` 时不直接调 EC，
而是 `merge_mate_pairs @ARGV | quorum_error_correct_reads db /dev/fd/0 | split_mate_pairs prefix`，
用 Perl `pipe`+`fork`+`exec` 串起三个子进程（stderr 重定向到 `<prefix>.log`）。EC 从 `/dev/fd/0`
读 stdin、输出 FASTA 到 stdout，最后由 `split_mate_pairs` 按相邻两行切回 `_1.fa`/`_2.fa`。
注意单端模式才是 EC 的默认丢弃语义；双端模式强制 `--no-discard`（见 §7）。

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

> 实现细节：底层 `vals_` 用 `val_array(bits + 1, ...)`——多出的 1 位正是
> 质量位；`max_val_ = (2^bits) - 1` 即计数封顶（`create_database -b` 校验
> 1-63）。写入 `.jf` 头时 `header->bits(vals_->bits() - 1)` 存回用户 bits。
> `create_database` 的 `-q/-Q`（min-qual-value/char）二选一互斥、且必须给出
> 其一；`-Q` 是单 ASCII 字符、`-q` 是整数值；`-p/--reprobe` 是独立选项（默认
> 126，Jellyfish 哈希最大 reprobe 次数）。

> **双链 canonical 计数**（`create_database.cc:77-85` 的 `quality_mer_counter`）：
> 每个位置同时维护正向 `m` 与反向互补 `rm` 两个滚动 mer
> （`m.shift_left(code)` / `rm.shift_right(complement(code))`），入库时取
> `m < rm ? m : rm`（canonical），因此**正反链 k-mer 合并成一个计数**。这与
> pgr `libs/kmer` 的 canonical 2-bit 思路一致（一条 read 及其反向互补只贡献
> 一个计数）。非 ACGT 碱基直接清零两个长度游标、跳过该位置（`not_dna`）。

> **无锁并发计数**：`add` 走 `keys_->set()`（Jellyfish `large_hash::array`
> 的原子探针 + reprobe）与 `vals_`（`atomic_bits_array<uint64_t>`）的 CAS 循环
> （`mer_database.hpp:94-113`），多线程对同一共享哈希做原子自增、无需互斥锁。
> 哈希满时 `handle_full_ary` 经 pthread barrier 同步后翻倍扩容迁移到新表
> （`mer_database.hpp:137-187`）——这正是 pgr 精确计数（radix sort + rayon）
> 之外的另一条并发路线，可作为对照。

> `quorum.in` 建库时把 bits 硬编码为 `-b 7`（计数上限 2^7-1=127，足够覆盖
> 一般覆盖度），并默认 `-s 200M`（Jellyfish 哈希槽位，估小会
> "Failed: Increase the size parameter"）。

哈希满时自动翻倍扩容（`handle_full_ary`，多线程屏障同步，旧表迁移到
2× 新表）。

### 3.2 `database_query`：只读 mmap 查询

纠错阶段把 `.jf` 库 mmap（或整读），`operator[](mer)` 返回
`(count, quality)`；`get_best_alternatives(m, counts[4], ucode, level)`：把 `m` 第 0 位替换为 A/C/G/T，逐个查询 canonical 计数，返回 4 个计数、最高质量等级（level）与命中数——**纠错的核心查询原语**。

> 精确语义（`mer_database.hpp:303-329`）：`counts[]` **只含最高质量等级**的替代——遍历
> A/C/G/T 时若遇到 quality 更高的替代，会把之前记录的低质量位置的 `counts[j]`（j<i）清零、
> 命中数 `count` 归零（`if(v.second > level && count>0) { for(j<i) counts[j]=0; count=0; }`）。
> 因此 `counts[ori]==0` 有两种成因：原碱基是错误（低质量、无高质量替代），或**原碱基低质量而
> 存在高质量替代**（其低质量计数被清零）——这正是 §4.2 `extend` 中"原碱基为错误、无高质量
> 候选"截断分支与候选替换分支的判定依据。

> **k 从头部恢复**：`.jf` 以每碱基 2 bit 存 mer，头部的 `key_len` 实为 `2k`；纠错端读库后
> 取 `mer_dna::k(mer_database.header().key_len() / 2)` 恢复 k（`error_correct_reads.cc:688`），
> 故建库与纠错的 k 必须一致（不一致时的 k 校验仅对污染库显式检查）。

## 4. 纠错算法（`error_correct_reads.cc`）

逐 read 处理，**从可信 anchor 向两端扩展**（`extend` 对 forward/backward
对称，`forward_ptr`/`backward_ptr` 模板统一方向语义）：

### 4.1 `find_starting_mer`：找 anchor

- 从 **`skip` 偏移**（`--skip`，默认 1）起的 5' 端滑窗（k-mer 窗口），
  遇 N 则重置并跳到下一个 k-mer（`shift_left` 返回 false 即重新装配）；
- 每个位置先查污染（命中且非 trim → 直接判 "Contaminated read"）；
- 非污染时取 **`get_val`（仅高质量计数）**，连续 `good`（`--anchor`→`-g`，
  默认 2）个 k-mer 计数 ≥ `anchor-count`（默认 3）即认为找到可信起点；
  找不到 → 整条丢弃，`--no-discard` 时输出单碱基 N。
- 跳过/丢弃的 read 会在 `.log` 里记一行 `Skipped <header>: <error>`，
  错误信息三种：`Contaminated read` / `No high quality mer` /
  `Entire read is an homopolymer`。

> anchor 判定细节（`error_correct_reads.cc:607-641`）：`found` 在 `get_val >= anchor_count`
> 时 +1、否则清零，累计 `found >= good` 即命中（`found = (int)val >= _ec.anchor() ? found+1 : 0;
> if(found >= _ec.good())`）。注意这里取的是**仅高质量**计数（`get_val` 对低质量 mer 返回 0），
> 因此 anchor 必须是高质量 k-mer。

### 4.2 `extend`：逐碱基扩展与纠错

每个位置把新碱基移入 k-mer，然后 `get_best_alternatives` 检查 4 种可能：

| 情况 | 处理 |
|---|---|
| `count==0`（无任何延续） | `truncation`，截断该端 |
| `count==1`（唯一延续） | 若与当前碱基不同 → `substitution`（替换）；相同则原样保留（`log_substitution` 内 `from==to` 直接返回 OK，`error_correct_reads.cc:361`） |
| `count>1` 且 `counts[ori]` > min-count 且（≥ cutoff 或质量够） | 保留原碱基 |
| `count>1` 且 `counts[ori]` > min-count 但（< cutoff 且质量不足） | **Poisson 碰撞检验**：`p = Σcounts × (先验错误率/3)`，`poisson_term(p, counts[ori]) < 阈值` → 视为随机碰撞、保留；否则落入候选替换 |
| `count>1` 且 `counts[ori]` ≤ min-count，且 `level==0 && counts[ori]==0`（原碱基为错误、无高质量候选） | `truncation` |
| 其余情况（含 N 碱基） | 进入候选替换：对每个计数 > min-count 的候选检查**延续性**（替换后移一位，下一个 k-mer 的 level ≥ 当前）；选计数最接近 `prev_count` 的候选；平局时用 read 下一个碱基仲裁；仍多个候选 → 不纠 |
| 原碱基为 N 且 `level==0`（所有候选均低质量） | `truncation`（`error_correct_reads.cc:457-460`） |
| 原碱基为 N 且替换后仍无延续（候选替换失败、`check_code<0`） | `truncation`（`error_correct_reads.cc:554-557`） |

> 表格中"质量够"即 `*qual >= qual_cutoff`（EC 的 `-q/-Q`）；其**默认值是
> `char` 最大值 127**——即默认情况下该分支几乎不因质量直接保留原碱基，除非
> 显式给 `-q`/`-Q` 压低阈值。另外 `prev_count` 用 `get_val`（仅高质量计数）
> 初始化，每步随 `count==1` 分支更新；候选替换里"选最接近 prev_count"
> 在 `prev_count <= min_count` 时退化为选**计数最大**的候选
> （`_prev_count = prev_count<=min_count ? UINT32_MAX : prev_count`，
> `error_correct_reads.cc:514`）。

### 4.3 `err_log`：窗口错误数限制（防过度纠错）

所有 sub/trunc 事件按位置记录；**滑动窗口（`-w`，quorum.in 不传时默认
10）内错误数 ≥ `-e`（默认 3）即触发回退截断**（`remove_last_window` 丢弃
窗口内事件并截断到窗口起点；`window()`/`error()` 仅当显式传 0 时才回退到
k 与 k/2）。输出日志：`pos:sub:from-to`、`pos:3_trunc`（3' 端）、
`pos:5_trunc`（5' 端）。

> `err_log` 语义（`err_log.hpp`）：事件按位置保序，`_lwin` 维护滑动窗口左边界；
> `check_nb_error` 在 `pos > _lwin.pos + window` 时推进 `_lwin`，返回
> `_log.size() - _lwin - 1 >= error`（窗口内错误数是否超限）。`substitution`/`truncation`
> 都返回该布尔值；`remove_last_window` 返回 `last.pos - lwin.pos`（即回退的碱基数），
> `extend` 据此把输出指针回退并截断到窗口起点。backward 方向经 `backward_log::truncation`
> 做 `pos-1` 修正（`error_correct_reads.hpp:170-172`）。

### 4.4 其他

- `compute_poisson_cutoff`：未显式给 `-p` 时，从计数分布自动估计 cutoff
  ——只统计**高质量** mer（编码值奇数且 ≥2），`coverage = total/distinct`
  （总 k-mer / 去重数，即平均计数），`lambda = coverage × 先验错误率/3`，
  从 `x=2` 起取首个满足 `poisson_term(lambda,x) < 阈值` 的 x 并**返回 x+1**。
  注意这里传给它的**阈值是 `poisson-threshold / apriori-error-rate`**（默认
  `1e-6 / 0.01 = 1e-4`），与 `extend` 内碰撞检验直接用的 `poisson-threshold`
  （1e-6）**不是同一个值**；若自动估计失败（返回 0）且未给 `-p`，程序
  `err::die("Cutoff computation failed. Pass it explicitly with -p switch.")`。
  另外 `poisson_term(λ,i)` 对 `i<11` 查阶乘表、`i≥11` 用 Stirling 近似。
- `homo_trim`（`--homo-trim`）：从 corrected read 末端向前扫描，逐位累计
  homopolymer 评分（`(same<<1)-1`：同碱基 +1、异 −1），记录累计最大值
  位置；**仅当最大评分 ≥ 阈值才在该位置截断**（否则不截）。
- `contaminant`：可选 Jellyfish 污染库（需与主库同 k），命中即丢弃
  （或 `--trim-contaminant` 截断）。

### 4.5 k-mer 表示与方向抽象（`kmer.hpp` / `error_correct_reads.hpp`）

- **dual-mer**：`kmer_t` 同时维护正向 `_fmer` 与反向互补 `_rmer` 两个滚动
  mer；`shift_left(c)` 同时做 `_fmer.shift_left` + `_rmer.shift_right(complement)`，
  `replace(i,x)` 同步写正反两条（`kmer.hpp:24-50`）；canonical 取
  `_fmer < _rmer ? _fmer : _rmer`（`kmer.hpp:43`）——滚动过程中正反链同步维护，
  与 pgr `libs/kmer` 的滚动 canonical 2-bit 编码**同构**（pgr 无需同时维护两条
  序列，因为 2-bit 编码天然自带互补对称）。
- **方向统一**：`forward_mer`/`backward_mer` 把 "shift" 抽象成同一 `shift()`
  接口，方向差异由适配器承担（`forward_mer::shift→shift_left`、
  `backward_mer::shift→shift_right`），再叠加 `forward_ptr`/`backward_ptr`
  （指针方向反转）与 `forward_counter`/`backward_counter`（坐标方向反转），
  让 `extend` **只用一套模板代码表达两个方向的扩展**（`error_correct_reads.hpp:16-149`）。
  这是 "方向无关" 的工程范本——pgr 若实现双向扩展可借鉴，但 Rust 下更自然的
  做法是 `forward`/`backward` 两个闭包或显式参数化，避免 C++ 的模板指针体操。

## 5. 参数语义（`quorum.in`）

| 参数 | 含义 |
|---|---|
| `-s` | Jellyfish 哈希槽位大小（默认 200M，必须容下全部 k-mer；估小报 "Failed: Increase the size parameter"） |
| `-t` | 线程数（默认自动检测 CPU 数） |
| `-p` | **输出前缀**（默认 `quorum_corrected`） |
| `-k` | k-mer 长度（默认 24，README 注明上限 31） |
| `-q` / `-m` | 质量下限字符（默认自动检测）/ 高质量阈值偏移（默认 5；高质量 = char ≥ q+m） |
| `-w` / `-e` | 错误窗口大小（默认 10）/ 窗口内最大错误数（默认 3）；显式传 0 时回退到 k / k/2 |
| `--min-count` / `--skip` | 好 k-mer 最小计数（EC 默认 1）/ 找 anchor 前跳过的碱基数（默认 1） |
| `--anchor` / `--anchor-count` | anchor 连续个数 / anchor k-mer 最小计数（EC 默认 2 / 3） |
| `--contaminant` / `--trim-contaminant` | 污染库（Jellyfish `.jf`）/ 命中即截断（而非丢弃） |
| `-d` / `-P` / `--homo-trim` | 不丢弃（输出单碱基 N）/ 双端分文件 / homopolymer 截断（需传整数值阈值） |
| `--debug` / `--version` / `-h` | 调试（回显执行的命令行）/ 版本 / 帮助 |

> **quorum.in 参数名与 `error_correct_reads` 内部名并不一一对应**：
> `--anchor` → EC 的 `-g/--good`（连续 good 数，默认 2）；`--anchor-count` →
> EC 的 `-a/--anchor-count`（默认 3）。EC 二进制还有若干 **quorum.in 未暴露** 的
> 开关：`-p/--cutoff`（Poisson cutoff，未给则自动估计）、`--apriori-error-rate`
> （默认 0.01）、`--poisson-threshold`（默认 1e-6）、`-q/-Q`（质量 cutoff 值/字符）、
> `--gzip`、`-M/--no-mmap`、`-v/--verbose`。因此上一版把 `-p` 记为 cutoff 是**错的**——
> 那是 EC 二进制内部的 `-p`，而 quorum.in 的 `-p` 是输出前缀；cutoff 在 quorum.in 流程里
> 只能靠 Poisson 自动估计，无法从脚本层显式指定。

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
  5. **canonical 双链计数**（正反链合并为一个计数）与 pgr 的 canonical
     2-bit 思路天然契合，直接沿用即可。
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
