# pgr fq trim-q：按质量分数修剪（设计稿）

> 定位：`fq` 子命令，按质量分数修剪读段两端（不做去接头）。与
> [sickle.md](../references/sickle.md)、[cutadapt.md](../references/cutadapt.md)
> 两篇参考笔记配套。设计输入来自业务侧真实用法——anchr 的
> `templates/trim.tera.sh#L135`。
>
> 状态：**设计稿（未实现）**。本文件先定设计，落地前需用户确认。

## 0. 需求来源：anchr 中 sickle 的真实用法

anchr（`wang-q/anchr`，commit `04f827a`）的 read 清洗流程里，sickle 被用作**第二道
质量修剪 + 参数扫描**。原始调用（`templates/trim.tera.sh#L135`）：

```bash
parallel --no-run-if-empty --linebuffer -k -j 2 "\
  mkdir -p Q{1}L{2}; cd Q{1}L{2}; \
  sickle pe -t sanger -q {1} -l {2} \
    -f ../opt.prefix1.fq.gz -r ../opt.prefix2.fq.gz \
    -o prefix1.fq -p prefix2.fq -s prefix.s.fq" \
  ::: opt.qual ::: opt.len
```

拆解出的真实需求：

| 需求 | 来源细节 |
|---|---|
| **双端修剪 + singles 分离** | `sickle pe` 产出 `1.fq`/`2.fq`/`s.fq`（singles） |
| **singles 二次修剪再合并** | 对 `s.fq.gz` 再跑 `sickle se`，append 回 `s.fq` |
| **质量/长度阈值可配** | `-q {1}`（质量阈值）、`-l {2}`（长度阈值） |
| **参数扫描** | `parallel` 遍历 qual×len 组合，各建 `Q{qual}L{len}/` 目录 |
| **sanger 编码** | `-t sanger`（质量偏移 33） |

> **要点**：anchr 的 hone 是"同批数据用不同质量/长度阈值各修剪一份，供下游组装
> 参数寻优"。故 `fq trim-q` 的核心价值不只是单次修剪，而是**低成本、可并列的
> 多阈值批量修剪**。

## 1. CLI 设计

### 1.1 命名

`pgr fq trim-q`（用户已确认）。`-q`（quality）显式表明"按质量分数修剪"，与
"去接头"（`trim`/trimming）区分。子命令用连字符风格，与 `pgr fq to-fa`、
`pgr fq interleave` 一致。

### 1.2 参数草案

```
pgr fq trim-q [options] <infiles...>

Input:
  <infiles...>  单端 1 个文件；双端 2 个文件（分别对应 R1/R2）

Options:
  -o, --outfile        输出文件（单端）
  --outfile-2           双端 R2 输出文件
  --outfile-single      双端 singles 输出文件
  -q, --qual-threshold  质量阈值（默认 20）
  -l, --length-threshold 长度阈值（默认 20）
  --method <sliding|mott>  修剪算法（默认 sliding）
  --no-fiveprime        禁用 5' 端修剪（仅 sliding 有效）
  --quality-base N      Phred 偏移（默认 33）
  --nextseq             以 NextSeq 变体修剪 polyG 尾巴（默认关）
```

### 1.3 与参数扫描的衔接

anchr 的场景是"一批阈值各跑一遍"。两个可选方案：

- **A（推荐，最小）**：`fq trim-q` 单次处理一组阈值。参数扫描由外层
  `parallel`/shell 循环完成（与 anchr 现状一致），pgr 不内置。
- **B**：内置多阈值扫描（类似 `--qual 15,20,25 --len 50,60,70` 笛卡尔积）。
  与 anchr 的 `Q{qual}L{len}` 目录结构对应，但属"便利功能"，需求不足不先做。

> 遵循 AGENTS.md「简洁优先」，先实现 A。B 留作未来方向（§5）。

## 2. 算法设计

### 2.1 双算法并存（`--method`）

| 算法 | 原理 | 来源 | 用途 |
|---|---|---|---|
| `sliding`（默认） | 窗口平均质量低于阈值即切 3' 端 | sickle / Trimmomatic `SLIDINGWINDOW` | 直观通用，默认 |
| `mott` | 累积质量取局部最大切点，可修中部低质量 | cutadapt `-q` / BWA `bwa_trim_read` | 更精细，备选 |

两算法都以相似方式结合 `--qual-threshold`/`--length-threshold` 与长度过滤。

### 2.2 滑窗核心（移植自 sickle `sliding_window`）

- 窗口大小 = `max(1, 0.1 × 读长)`（自适应）。
- 5' 端：窗口平均质量首次 ≥ 阈值时，取窗口内首个达阈值碱基为切点。
- 3' 端：窗口平均质量 < 阈值时，取窗口内首个低于阈值碱基为切点。
- 修剪后 `three - five < length_threshold` 则丢弃。
- 窗口滑动用差分更新（`window_total -= 首碱基; += 新碱基`），O(1)。

> 移植时注意 [sickle.md](../references/sickle.md) §2.2 记录的**冗余死代码**：
> `window_start+window_size > qual.l` 的"最后窗口"判定恒为假，忽略即可。

### 2.3 Mott 核心（移植自 cutadapt `quality_trim_index`）

- 5'/3' 端独立计算：单遍累积 `cutoff - q`，累积和首次转负停止，切点在累积和
  局部最大处。
- 返回 `(start, stop)` 左闭右开区间；`start >= stop` 时置 `(0,0)`。
- 复杂度 O(n)。参考 Rust 实现见 [cutadapt.md](../references/cutadapt.md) §4.1。

### 2.4 双端 + singles（对齐 sickle pe 语义）

```
对每对 (R1, R2)：
  两端都通过 → 写 R1 输出、R2 输出
  仅一端通过 → 通过端写 singles
  两端都失败 → 丢弃
```

- 双端文件记录数不匹配：警告并只处理公共部分，不 panic。
- singles 再修剪：这是 anchr 的流程（先用 pe 拿 singles，再 se 修剪 append），
  属于**外层编排**，pgr 不内置，用户用管道/脚本完成。

## 3. 与 pgr 现有约束对接

| 参考项 | 落地要求 |
|---|---|
| 零 panic | 质量字符越界（`qual - base` 得负或超范围）须返回 `anyhow` 错误，不静默处理 |
| 文件名去重 | 输入/输出不可相同，防覆盖（对齐 `ensure_outfile_distinct` 硬约束） |
| 质量编码 | 支持 Sanger(+33)/Illumina(+64)/Solexa(+64 近似)，`--quality-base` 可配 |
| 分层 | 算法（滑窗/Mott）放 `libs/fq/`，`cmd_pgr/fq/` 仅做 clap 编排 |
| 依赖 | 用 noodles 读 FASTQ，不引入新依赖 |

## 4. 测试计划

- **单元**：滑窗与 Mott 对已知质量串的切点断言（两端质量高/低、中部低质量区、
  空序列、长度阈值边界）。
- **集成**：`tests/cli_fq_trim_q.rs`，覆盖单端/双端/singles、`--method` 两种、
  `-q`/`-l` 阈值、质量越界报错、输入输出同名报错。
- **对照**：与 `sickle` 同参数跑同一对文件，比对修剪后序列（滑窗语义应一致）。

## 5. 未来方向（暂不做）

- 内置多阈值参数扫描（§1.3 方案 B）。
- `--nextseq` polyG 变体默认开启（当前默认关，待证据）。
- `expected_errors`（Edgar 2015）整条读错误数过滤，作为独立 `fq` 功能。

---

*参考来源: [sickle.md](../references/sickle.md) | [cutadapt.md](../references/cutadapt.md) | [anchr trim.tera.sh](https://github.com/wang-q/anchr/blob/04f827afe37d5f40f12cd0602d54086cf8b0078c/templates/trim.tera.sh)*