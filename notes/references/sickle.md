# sickle: 基于滑动窗口的 FASTQ 自适应质量修剪

> 整理于 2026-08，源自对 `sickle-master/`（v1.33）源码的分析。目的：pgr 的 `fq` 命令目前只有 `to_fa`（FASTQ→FASTA）与 `interleave`（双端交错）两个子命令，**没有质量修剪功能**。本文分析 sickle 的滑动窗口自适应修剪算法，为 pgr 未来可能的 `fq trim` 提供算法参考。本文聚焦算法与 CLI 设计，`kseq.h` 等 I/O 基础设施仅作背景（pgr 用 noodles 已覆盖）。

## 1. 简介

`sickle`（Nikhil Joshi, UC Davis Bioinformatics Core, 2011）是一个对 FASTQ 读段做**质量自适应修剪**的工具。其核心思想：高通量测序读段质量在 3' 端（部分也在 5' 端）逐渐退化，错误碱基会污染后续组装/比对。sickle 用滑动窗口配合质量阈值与长度阈值，确定该在何处修剪 3' 端，以及（可选）何处修剪 5' 端，并按长度阈值丢弃过短读段。

- **窗口大小自适应**：窗口长度 = `0.1 × 读段长度`；若小于 1 则取读段全长。
- **5' 修剪**：从 5' 端滑动，当窗口平均质量**首次超过阈值**时，找到窗口内首个达阈值的碱基位置作为 5' 切点。可禁用（`-x`）。
- **3' 修剪**：当窗口平均质量**低于阈值**（或到达读段末端）时，找到窗口内首个低于阈值的碱基位置作为 3' 切点。
- **丢弃**：修剪后剩余长度 < 长度阈值则整条丢弃（或在 `-M` 模式下输出单碱基 `N` 记录以保持配对）。
- **质量编码**：支持 Sanger / Illumina / Solexa（Solexa 为线性近似）。Illumina 1.3–1.7 用 `+64`，CASAVA ≥ 1.8 即 Sanger（`+33`）。
- **格式细节**：输出时可把 `+` 行后的重复 header 替换为单个 `+`（CASAVA ≥ 1.8 默认格式）；支持 gzip 输入/可选 gzip 输出；`-n` 可在首个 N 处截断。

> **范围说明**：pgr 不需要复现 sickle 的 C 实现。若实现 **`fq trim-q`**（按质量分数修剪，见 §4），应移植其**滑动窗口修剪算法**（`sliding.c`，约 100 行 C）。`kseq.h` 的流式读取由 noodles 替代，`getopt` CLI 由 clap 替代。

## 2. 核心概念 (Key Concepts)

### 2.1 质量编码表（`sickle.h:82-88`）

| 类型 | offset | min | max | 说明 |
| :--- | :--- | :--- | :--- | :--- |
| PHRED | 0 | 4 | 60 | 未在 CLI 暴露，仅代码内定义 |
| SANGER | 33 | 33 | 126 | CASAVA ≥ 1.8 |
| SOLEXA | 64 | 58 | 112 | **近似**，真实转换为非线性 |
| ILLUMINA | 64 | 64 | 110 | CASAVA 1.3–1.7 |

质量值 = `ASCII 码 - offset`。`get_quality_num`（`sliding.c:10`）读取时校验质量字符是否落在 `[min, max]` 区间，越界即报错退出（**这是硬校验，pgr 移植时应保留友好报错而非 panic**）。

> 注意：`sickle se` / `pe` 的 `-t` 只接受 `solexa|illumina|sanger` 三值，`PHRED` 类型未在命令行暴露。

### 2.2 滑动窗口修剪算法（`sliding.c:35` `sliding_window`）

单条读段的处理流程：

```
1. 若 seq.l < length_threshold → 返回 five=three=-1（丢弃整条）
2. window_size = max(1, 0.1 * seq.l)   // 窗口长度自适应
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

- **滑动用差分而非重算**：`window_total -= 首碱基; window_total += 新碱基`，O(1) 滑动整个窗口（`sliding.c:107-111`）。pgr 可用滑动窗口和值或前缀和实现。
- **`window_start+window_size > qual.l` 作为"最后窗口"判定**：注意循环条件已是 `i <= qual.l - window_size`，故迭代内该条件实际恒为假（`qual.l - window_size + window_size > qual.l` 为假）。这是一段**冗余/死代码**（`sliding.c:92`），3' 修剪实际只由 `window_avg < qual_threshold` 触发。移植时只需等价于"平均质量低于阈值"。
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

单一输入 FASTQ → 单一修剪输出。参数：`-f` 输入、`-t` 质量类型、`-o` 输出、`-q` 质量阈值（默认 20）、`-l` 长度阈值（默认 20）、`-x` 禁用 5' 修剪、`-n` 截断 N、`-g` gzip 输出、`-z` 静默。

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
- `-c` 与 `-f/-r/-o/-p` 互斥（`trim_paired.c:246`）。
- `-m` 与 `-M` 二选一（`trim_paired.c:250`）。
- `-m` 必须配 `-s`，`-M` 不能配 `-s`（`trim_paired.c:254`）。
- 输入/输出文件名**全对去重**，禁止覆盖（`trim_paired.c:259,295`）——与 pgr 的 `ensure_outfile_distinct` 约束一致。
- 双端文件长度不匹配时给出警告并截断（`trim_paired.c:372,475`）。

## 4. 与 pgr 的关联性

### 4.1 现状对比

pgr 的 `fq` 命令目前只有 `to_fa`（FASTQ→FASTA）与 `interleave`（双端交错），**没有质量修剪子命令**。sickle 填补的正是这个空白。

### 4.2 潜在移植点：`pgr fq trim-q`

> **命名说明**：质量修剪子命令命名为 `fq trim-q`（而非 `fq trim`），用 `-q`（quality）明确表示"按质量分数修剪"，避免与"去接头"(`trim`/trimming) 混淆。

若未来实现，可参考但不直接照搬：

| sickle 概念 | pgr 移植建议 |
| :--- | :--- |
| 窗口 = 0.1×读长 | 可保留自适应窗口，或提供固定 `--window` 选项 |
| 质量阈值 / 长度阈值 | 对应 `--qual-threshold` / `--length-threshold`（默认 20/20） |
| 5' 修剪（可禁用） | `--no-fiveprime` 开关 |
| 质量编码表 | 用 noodles 的 Phred 解码，或按 offset 表实现 |
| 3' 修剪 + 丢弃 | 核心算法，`three - five < len` 则丢弃 |
| `-n` 首个 N 截断 | 可选开关 |
| `-M` 单碱基 N 保配对 | 双端交错模式特有，`interleave` 命令可复用 |
| "singles" 文件 | 双端配对打破时单独输出 |

**算法体量**：核心 `sliding_window` 约 100 行 C（含质量校验）。Rust 移植后预计 60–100 行，无复杂数据结构，是标准"分隔 + 差分滑动窗口"。

### 4.3 值得借鉴的健壮性设计

1. **质量值越界即报错**（`sliding.c:21`）：不是静默裁剪。pgr 的项目硬约束是"任何用户输入都不应 panic"，此处应返回 `anyhow` 错误而非直接 `exit(1)`。
2. **文件名去重校验**（`trim_paired.c:259,295`）：输入输出不可相同，防止覆盖输入——与 pgr 的 `ensure_outfile_distinct` 硬约束精神一致。
3. **双端不匹配的容错**：文件长短不一时警告并只处理公共部分，不崩。
4. **`+` 行 header 丢弃**：输出统一为单个 `+`，避免嵌套 header 不一致。

### 4.4 与 pgr 现有约束的冲突点

- sickle 用 `exit(1)` 处理错误，pgr 必须改为 `anyhow::Result` + 友好错误信息。
- sickle 的质量校验是"越界即终止整条程序"，pgr 的零 panic 原则会更倾向于"记录错误并跳过/继续"，需在移植时决策。
- sickle 的 SOLEXA 是线性近似，若 pgr 需要严格正确应实现非线性转换或明确标注近似。

### 4.5 现代质量修剪算法对比（2026 视角）

> 调研补充（2026-08）：**按质量分数修剪**（非去接头）在现代工具中仍是这几类算法，滑窗并未过时。sickle 的"旧"不在算法，而在缺乏多线程/自动接头/报告等工程能力；但这些对**纯质量修剪**不重要。真正不同的算法方向是 cutadapt 的 Mott 累积质量法（更精细）。

| 算法 | 代表工具 | 原理 | 特点 |
| :--- | :--- | :--- | :--- |
| **滑动窗口** | sickle、Trimmomatic `SLIDINGWINDOW`、fastp | 窗口内平均质量低于阈值即切 3' 端 | 默认/主流；Trimmomatic 的同类滑窗是事实标准 |
| **Mott 算法** | cutadapt `-q` | 从 3' 端逐碱基累加质量，`累积质量 - 位置惩罚` 决定切点，可修中间低质量区 | 每碱基阈值，比滑窗精细；是唯一真正不同的算法方向 |
| **Leading/Trailing** | Trimmomatic、fastp | 只切两端低于阈值的连续碱基，不处理中间 | 最保守，常配合滑窗使用 |
| **逐读动态窗口** | fastp | 滑窗 + 基于读长的动态窗口，单遍 O(n) | 对短读优化，作者将复杂度从 O(n²) 降为 O(n) |

**对 `fq trim-q` 的建议**：
- 滑窗仍是默认选择（忠实对齐 Trimmomatic `SLIDINGWINDOW` 语义，最通用）。
- 可将 **Mott 算法**作为 `--method mott` 选项与滑窗（`--method sliding`）并存，提供更精细的修剪路径。
- 现代工具差异主要在**多线程**与**单遍多操作**（质控统计+过滤+修剪一次扫描），这些 pgr 用 noodles + rayon 天然契合，不作为算法选型依据。

---

*参考来源: [sickle GitHub](https://github.com/najoshi/sickle) | 本项目源码 `sickle-master/`（v1.33）*