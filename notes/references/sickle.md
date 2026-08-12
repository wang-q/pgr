# sickle: 基于滑动窗口的 FASTQ 自适应质量修剪

> 整理于 2026-08，源自对 `sickle-master/`（v1.33）源码的分析。目的：sickle 的滑动窗口质量修剪算法是 pgr `fq trim-qual`（`src/libs/fq/trim.rs`）中默认 `--method sliding` 的直接算法来源。本文分析 sickle 的滑动窗口自适应修剪算法、CLI 与配对逻辑，并记录 pgr 移植时的取舍与两者差异。`kseq.h` 等 I/O 基础设施仅作背景（pgr 用 noodles 已覆盖）。

## 1. 简介

`sickle`（Nikhil Joshi, UC Davis Bioinformatics Core, 2011）是一个对 FASTQ 读段做**质量自适应修剪**的工具。其核心思想：高通量测序读段质量在 3' 端（部分也在 5' 端）逐渐退化，错误碱基会污染后续组装/比对。sickle 用滑动窗口配合质量阈值与长度阈值，确定该在何处修剪 3' 端，以及（可选）何处修剪 5' 端，并按长度阈值丢弃过短读段。

- **窗口大小自适应**：窗口长度 = `(int)(0.1 × 读段长度)`（int 截断）；若为 0（读长 <10bp）则取读段全长。
- **5' 修剪**：从 5' 端滑动，当窗口平均质量**首次超过阈值**时，找到窗口内首个达阈值的碱基位置作为 5' 切点。可禁用（`-x`）。
- **3' 修剪**：当窗口平均质量**低于阈值**（或到达读段末端）时，找到窗口内首个低于阈值的碱基位置作为 3' 切点。
- **丢弃**：修剪后剩余长度 < 长度阈值则整条丢弃（或在 `-M` 模式下输出单碱基 `N` 记录以保持配对）。
- **质量编码**：支持 Sanger / Illumina / Solexa（Solexa 为线性近似）。Illumina 1.3–1.7 用 `+64`，CASAVA ≥ 1.8 即 Sanger（`+33`）。
- **格式细节**：输出时可把 `+` 行后的重复 header 替换为单个 `+`（CASAVA ≥ 1.8 默认格式）；支持 gzip 输入/可选 gzip 输出；`-n` 可在首个 N 处截断。

> **范围说明**：pgr 已在 `fq trim-qual` 中移植了 sickle 的滑动窗口修剪算法（`--method sliding`，见 §4），无需复现其 C 实现或 `kseq.h`/`getopt` 层。本文同时对比 pgr 与 sickle 的实现差异，并记录 pgr 在滑动窗口之外新增的算法（Mott、poly-G）与 CLI 取舍。

## 2. 核心概念 (Key Concepts)

### 2.1 质量编码表（`sickle.h:82-88`）

| 类型 | offset | min | max | 说明 |
| :--- | :--- | :--- | :--- | :--- |
| PHRED | 0 | 4 | 60 | 未在 CLI 暴露，仅代码内定义 |
| SANGER | 33 | 33 | 126 | CASAVA ≥ 1.8 |
| SOLEXA | 64 | 58 | 112 | **近似**，真实转换为非线性 |
| ILLUMINA | 64 | 64 | 110 | CASAVA 1.3–1.7 |

质量值 = `ASCII 码 - offset`。`get_quality_num`（`sliding.c:10`）读取时校验质量字符是否落在 `[min, max]` 区间，越界即 `fprintf` 报错并 `exit(1)`（`sliding.c:21-28`）——**这是硬校验**。pgr 移植时保留为友好报错而非 panic（见 §4.4）。

> 注意：`sickle se` / `pe` 的 `-t` 只接受 `solexa|illumina|sanger` 三值，`PHRED` 类型未在命令行暴露。

### 2.2 滑动窗口修剪算法（`sliding.c:35` `sliding_window`）

单条读段的处理流程：

```
1. 若 seq.l < length_threshold → 返回 five=three=-1（丢弃整条）
2. window_size = (int)(0.1 * seq.l)；若为 0（读长 <10bp）则取 seq.l（全长）
3. 初始化窗口 [0, window_size) 的质量和 window_total
4. 对每个窗口起点 i（0..=qual.l - window_size）：
   a. window_avg = window_total / window_size
   b. 5' 修剪（未禁用且尚未找到）：若 window_avg >= qual_threshold，
      在窗口内找首个 >= threshold 的碱基位置作为 five_prime_cut
   c. 3' 修剪（已找到 5' 或禁用 5'）：若 window_avg < threshold
      或到达读段末端，在窗口内找首个 < threshold 的碱基位置作为
      three_prime_cut，然后 break
   d. 滑动：window_total -= 窗口首碱基；window_total += 新进入碱基
5. 若 -n 启用且序列含 N/n，three_prime_cut = 首个 N 的位置
6. 若 (始终未找到 5' 且未禁用) 或 (three - five < length_threshold)
   → 两者置 -1（丢弃）
```

**关键实现细节**（pgr 移植时需注意）：

- **数据结构 `cutsites`**（`sickle.h:90-93`）：仅含两个 `int` 字段 `five_prime_cut`、`three_prime_cut`，`sliding_window` 通过 `malloc` 返回指针（成功与丢弃两种路径都分配），调用方（`trim_single.c:217` / `trim_paired.c:471-472`）用完 `free`。
- **切点初始值**：`three_prime_cut = seq.l`（默认不切 3'）、`five_prime_cut = 0`、`found_five_prime = 0`（`sliding.c:41-43`）。若整个循环从未触发 3' 分支，`three_prime_cut` 保持 `seq.l`，即 3' 端不做任何修剪。
- **3' 分支的守卫条件** `(found_five_prime == 1 || no_fiveprime)`（`sliding.c:93`）：当 `-x` 未开且 5' 从未找到时，3' 分支永不进入，循环滑完全长后由第 6 步 `found_five_prime == 0 && !no_fiveprime` 兜底丢弃——这正是"5' 未找到即丢弃"的实现机制。
- **滑动用差分而非重算**：`window_total -= 首碱基; window_total += 新碱基`，O(1) 滑动整个窗口（`sliding.c:107-111`）。pgr 用相同的差分滑动实现。
- **差分滑动的越界守卫**：加新碱基前先判 `window_start+window_size < qual.l`（`sliding.c:108`），最后一个窗口迭代对 `window_total` 的修改因无下一轮而实际未使用。pgr 的 `sliding_cut` 沿用了同一守卫（`trim.rs:175`）。
- **`seq.l` 与 `qual.l` 混用**：窗口大小与长度阈值取 `fqrec->seq.l`，而滑动循环上界用 `fqrec->qual.l`（`sliding.c:37,49,64`）。对规范 FASTQ 二者相等；pgr 的 `sliding_cut` 统一以 `qual.len()` 为准，无此混用。
- **`window_start+window_size > qual.l` 作为"最后窗口"判定**：注意循环条件已是 `i <= qual.l - window_size`，故迭代内该条件实际恒为假（`qual.l - window_size + window_size > qual.l` 为假）。这是一段**冗余/死代码**（`sliding.c:92`），3' 修剪实际只由 `window_avg < qual_threshold` 触发。pgr 的 `sliding_cut` 移植时确实丢弃了该死条件（`trim.rs` 只保留 `avg < threshold`）。
- **5' 切点是窗口内首个达阈值碱基，3' 切点是窗口内首个低于阈值碱基**：两者都返回**绝对位置**（基于原始序列），不是相对窗口偏移。
- **边界行为**：5' 与 3' 切点可能重叠（如读段质量整体很差时 five > three），由第 6 步长度检查兜底丢弃。
- **`-n` 截断在滑动窗口之后**：截断到首个 N 处，再走第 6 步长度检查。

### 2.3 修剪记录的输出（`print_record.c`）

```
@name [comment]
seq[five..three]
+
qual[five..three]
```

- 输出长度 = `three_prime_cut - five_prime_cut`（已按 5' 切点偏移）。
- `+` 行固定输出单个 `+`（丢弃原 header 内容）。
- `print_record_N`（`-M` 模式）：整个记录替换为单碱基 `N`，质量取 `quality_constants[qualtype][Q_MIN]`（该编码最小值），header 保留原 name/comment。**用于保持配对**。

## 3. 双端（`pe`）与单端（`se`）模式

### 3.1 单端 `sickle se`

单一输入 FASTQ → 单一修剪输出。参数：`-f` 输入、`-t` 质量类型、`-o` 输出、`-q` 质量阈值（默认 20）、`-l` 长度阈值（默认 20）、`-x` 禁用 5' 修剪、`-n` 截断 N、`-g` gzip 输出、`-z` 静默，另有 `-d` 调试输出（逐窗口打印 avg 与切点，`trim_single.c:143-145`）。

- **`-f/-t/-o` 三者皆必需**，缺一即 usage 报错（`trim_single.c:161`）；`-f` 与 `-o` 同名会被拒绝（`trim_single.c:165`）。
- 输入经 `gzopen` 读取，透明支持 gzip 与普通文件；`-g` 仅控制输出是否 gzip。
- 默认在 stdout 打印统计（Total FastQ records / kept / discarded，`trim_single.c:220`），`-z` 静默关闭。`pe` 的统计更细：paired kept/discarded、singles kept/discarded（分离文件形态还区分来自 PE1/PE2 各多少，`trim_paired.c:482-494`）。

### 3.2 双端 `sickle pe`

三种输入/输出形态：

| 形态 | 输入 | 输出 | 说明 |
| :--- | :--- | :--- | :--- |
| 分离文件 | `-f` 正向 + `-r` 反向 | `-o` 正向 + `-p` 反向 + `-s` singles | 两文件记录数须相同 |
| 交错文件 | `-c` 交错输入 | `-m` 交错输出 + `-s` singles | 交错成对 |
| 交错全保留 | `-c` 交错输入 | `-M` 单交错输出 | 丢弃的读段输出单碱基 N，保持配对记录数 |

**配对逻辑**（`trim_paired.c:369` 主循环）：

```
对每对 (fqrec1, fqrec2)：
  若两端都通过 → 输出两记录（分离文件分别写，交错文件依次写）
  若仅一端通过 → 通过端写 singles；-M 模式则通过端原样 + 失败端 N
  若两端都失败 → 丢弃；-M 模式则输出两条 N
```

**"singles" 文件**：只在一个方向通过过滤的读段，配对关系被打破后单独存放。

**参数校验较严格**（pgr 移植时值得参考的错误处理）：
- 分离文件形态要求 `-f/-r/-o/-p/-s` 五个参数**同时给出**，缺一即报错（`trim_paired.c:287-289`）；`-s` 在本形态下为必填（singles 文件）。
- `-c` 与 `-f/-r/-o/-p` 互斥（`trim_paired.c:246`）。
- `-m` 与 `-M` 二选一（`trim_paired.c:250`），且 `-m`/`-M` 共用同一输出参数 `outfnc`，都写交错输出文件（`trim_paired.c:168-178`）。
- `-m` 必须配 `-s`，`-M` 不能配 `-s`（`trim_paired.c:254`）。
- 输入/输出文件名**全对去重**，禁止覆盖（`trim_paired.c:259,295`）——与 pgr 的 `ensure_outfile_distinct` 约束一致。
- 双端文件长度不匹配时给出警告并截断（`trim_paired.c:372,475`），且 `pe` 同样有 `-d` 调试选项。

## 4. 与 pgr 的关联性

### 4.1 现状对比

pgr 的 `fq` 命令目前已拥有 `to-fa`、`interleave`、`norm`、`range`、`sample`、`split`、`clean`、`filter`、`clump`、`trim-qual` 等子命令（`src/cmd_pgr/fq/mod.rs`）。其中 **`trim-qual`**（`src/cmd_pgr/fq/trim_qual.rs`，实现位于 `src/libs/fq/trim.rs`）正是本笔记主题：**按质量分数修剪**，且默认算法即 sickle 式滑动窗口。

即：笔记最初撰写时"pgr 没有质量修剪功能"的现状已改变，sickle 填补的空白已由 `fq trim-qual` 实现。sickle 分析的价值从"未来移植参考"转为"对照 pgr 实际实现、澄清差异与取舍"。

### 4.2 pgr `fq trim-qual` 对 sickle 的移植与差异

pgr `trim-qual` 核心 `sliding_cut`（`trim.rs:143`）忠实对齐了 sickle 的 `sliding_window`，但有以下差异：

| 维度 | sickle | pgr `trim-qual` |
| :--- | :--- | :--- |
| **窗口大小** | `(int)(0.1*l)`，为 0 则取全长 | `max(1, n/10)` |
| **5'/3' 切点** | 窗口内首个达/低于阈值的绝对位置 | 相同 |
| **差分滑动** | `window_total ±=` | 相同 |
| **死条件** `window_start+window_size>l` | 存在（死代码） | 已删除 |
| **5' 未找到 → 丢弃** | `found_five_prime==0 && !no_fiveprime` | `sliding_cut` 返回 `None` → `trim_interval` 丢弃 |
| **质量编码** | 4 类型表（SANGER/SOLEXA/ILLUMINA + 未暴露 PHRED），`-t` 必填 | `--quality-base` 33/64/`auto`（默认 auto，BBDuk flip-flop 自动检测），无 SOLEXA 近似 |
| **质量校验** | 越界 `exit(1)` | `validate_quality` 校验 Phred ∈ [0,93]，越界返回 `anyhow` 错误（不 panic） |
| **N 截断 `-n`** | 有 | **无**（pgr 未移植） |
| **额外修剪** | 无 | 新增 `--method mott`（Mott 累积质量）、`--polyg-right`（3' poly-G） |
| **输出 `+` 行** | 单个 `+`，保留 name/comment | 相同 |
| **配对输出** | `-M` 输出 N 保配对 / `-m -s` singles | 双端可 `--outfile-2` 分离或省略为交错；`--outfile-single` 收 singles；**无 `-M` 保配对模式**，失败端直接丢弃 |

**值得注意的行为差异**：

1. **短读窗口大小**：对读长 `<10bp`，sickle 把窗口取为**读段全长**，pgr 取 `max(1, n/10)=1`。对默认长度阈值 20 而言此类读段通常早已被丢弃，故实际影响有限，但语义并不完全一致。
2. **配对保配策略**：sickle 的 `-M` 模式会把失败读段输出为单碱基 `N` 以保持记录数/配对。pgr `trim-qual` **未实现**此模式——双端仅一端通过时通过端写入 `--outfile-single`，两端都失败则直接丢弃，记录数不保持。pgr 的 `interleave` 等命令不与 `trim-qual` 联动补齐该行为。
3. **质量校验时机（顺序）**：sickle 的"读长 < 长度阈值即丢弃"发生在任何质量读取之前（`sliding.c:49-54` 先于窗口初始化），短读根本不调用 `get_quality_num`，故无效质量的短读会被静默丢弃、不报错；pgr 的 `trim_interval` 却是先 `validate_quality`（`trim.rs:270`）再查长度阈值（`trim.rs:271`），即对即将被丢弃的短读也做质量校验并可能 `bail`。这是"零 panic + 友好报错"取向带来的严格化差异：pgr 对坏质量短读会报错，而 sickle 直接丢弃。
4. **5' 切点内层回退**：sickle 的 5' 内层循环必能找到首个达阈碱基（`sliding.c:77-82`）；pgr 用 `.find(...).unwrap_or(window_start)`（`trim.rs:163-165`）补了一个**实际不可达**的默认回退（`avg >= threshold` 已保证窗口内至少一个达阈碱基）。两者结果等价，pgr 多出的 `unwrap_or` 是防御性写法。

### 4.3 值得借鉴的健壮性设计（pgr 已落实）

1. **质量值越界即报错**（`sliding.c:21`）：不是静默裁剪。pgr 项目硬约束是"任何用户输入都不应 panic"，此处已落实为 `validate_quality` 返回 `anyhow` 错误（`trim.rs:245`），而非 sickle 的 `exit(1)`。
2. **文件名去重校验**（`trim_paired.c:259,295`）：输入输出不可相同，防止覆盖输入——pgr `trim_qual::execute` 用 `ensure_outfile_distinct`（`trim_qual.rs:125`）及输出文件两两去重落实。
3. **双端不匹配的容错**：文件长短不一时警告并只处理公共部分，不崩——pgr `run_paired` 用 `warn_pair_mismatch` 实现，只警告一次（定义于 `trim.rs:457`）。
4. **`+` 行 header 丢弃**：输出统一为单个 `+`，避免嵌套 header 不一致——pgr `write_record` 保留原 name/comment 但 `+` 行固定输出（`trim.rs:296-313`）。

### 4.4 与 pgr 现有约束的差异点（移植取舍）

- sickle 用 `exit(1)` 处理错误；pgr 一律用 `anyhow::Result` + 友好错误信息，且质量校验选择"报错终止整条读段所在流程"（`validate_quality` 对整个记录 `bail`），而非"记录错误并跳过/继续"——这是对零 panic 原则的落实，但与早期笔记"倾向于跳过"的预判不同。
- sickle 的 SOLEXA 是线性近似；pgr 未实现 SOLEXA，只区分 Phred33/64，避免引入近似误差。
- sickle 的质量类型由 `-t` 强制指定；pgr 默认 `auto` 自动检测（BBDuk flip-flop 启发式），更省心但也引入检测不确定性（`--quality-base` 可显式覆盖）。
- pgr 未移植 sickle 的 `-n` 首个 N 截断，而是以 `--polyg-right` 填补 3' 端杂质（poly-G）场景——两者不重叠，是功能取舍而非替代。

### 4.5 现代质量修剪算法对比（2026 视角）与 pgr 现状

> 调研补充（2026-08）：**按质量分数修剪**（非去接头）在现代工具中仍是这几类算法，滑窗并未过时。sickle 的"旧"不在算法，而在缺乏多线程/自动接头/报告等工程能力；但这些对**纯质量修剪**不重要。真正不同的算法方向是 cutadapt 的 Mott 累积质量法（更精细）。

| 算法 | 代表工具 | 原理 | 特点 |
| :--- | :--- | :--- | :--- |
| **滑动窗口** | sickle、Trimmomatic `SLIDINGWINDOW`、fastp | 窗口内平均质量低于阈值即切 3' 端 | 默认/主流；Trimmomatic 的同类滑窗是事实标准 |
| **Mott 算法** | cutadapt `-q` | 从 3' 端逐碱基累加质量，`累积质量 - 位置惩罚` 决定切点，可修中间低质量区 | 每碱基阈值，比滑窗精细；是唯一真正不同的算法方向 |
| **Leading/Trailing** | Trimmomatic、fastp | 只切两端低于阈值的连续碱基，不处理中间 | 最保守，常配合滑窗使用 |
| **逐读动态窗口** | fastp | 滑窗 + 基于读长的动态窗口，单遍 O(n) | 对短读优化，作者将复杂度从 O(n²) 降为 O(n) |

**pgr 落地方案**：
- 滑窗作为默认 `--method sliding`（忠实对齐 sickle / Trimmomatic `SLIDINGWINDOW` 语义，最通用）——已实现。
- **Mott 算法**作为 `--method mott` 与滑窗并存，提供更精细的修剪路径——已实现（`trim.rs` `mott_cut`，移植自 cutadapt `qualtrim.pyx` 的 `quality_trim_index`），正是当初建议的"两种方法并存"。
- 现代工具差异主要在**多线程**与**单遍多操作**（质控统计+过滤+修剪一次扫描）；pgr 用 noodles + rayon 天然契合。`fq clean`/`fq filter`（bbduk 风格）已实现"adapter 修剪 + 质量过滤"的组合，与 `trim-qual`（纯质量）分工互补（见 `docs/fq.md`）。

---

*参考来源: [sickle GitHub](https://github.com/najoshi/sickle) | 本项目源码 `sickle-master/`（v1.33） | pgr `src/libs/fq/trim.rs` + `src/cmd_pgr/fq/trim_qual.rs`*
