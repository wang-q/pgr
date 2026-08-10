# pgr fq trim-qual：按质量分数修剪（设计稿）

> 定位：`fq` 子命令，按质量分数修剪读段两端（不做去接头）。与
> [sickle.md](../references/sickle.md)、[cutadapt.md](../references/cutadapt.md)
> 两篇参考笔记配套。设计输入来自业务侧真实用法——anchr 的
> `templates/trim.tera.sh#L135`。
>
> 状态：**已实现（2026-08）**，设计细节见 [fq-trim-qual.md](fq-trim-qual.md) 全文与 §7。

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
> 参数寻优"。故 `fq trim-qual` 的核心价值不只是单次修剪，而是**低成本、可并列的
> 多阈值批量修剪**。

## 1. CLI 设计

### 1.1 命名

`pgr fq trim-qual`（用户已确认）。`-q`（quality）显式表明"按质量分数修剪"，与
"去接头"（`trim`/trimming）区分。子命令用连字符风格，与 `pgr fq to-fa`、
`pgr fq interleave` 一致。

### 1.2 参数（已定稿）

```
pgr fq trim-qual [options] <infiles...>

Input:
  <infiles...>  单端 1 个文件；双端 2 个文件（分别对应 R1/R2）

Options:
  -o, --outfile        输出文件（单端 / 双端 interleaved）
  --outfile-2           双端 R2 输出文件（不给出则 R2 与 R1 interleaved 写入 -o）
  --outfile-single      双端 singles 输出文件
  -q, --qual-threshold  质量阈值（默认 20，允许浮点）
  -l, --length-threshold 长度阈值（默认 20）
  --method <sliding|mott>  修剪算法（默认 sliding）
  --no-fiveprime        禁用 5' 端修剪
  --quality-base <33|64|auto>  输入质量编码（默认 auto，自动检测）
  --polyg-right N       修剪 ≥N 的 3' polyG 尾巴（默认 0=关）
```

**双端输出形态（BBTools 式，隐式）**：给了 `--outfile-2` 就分离输出 R1/R2；不给
则 R2 与 R1 interleaved 写入 `-o`。`--outfile-single` 可选，给则接收幸存单端。
双端多输出均要求文件路径，只有单端允许 `-o stdout`。

### 1.3 与参数扫描的衔接

anchr 的场景是"一批阈值各跑一遍"。两个可选方案：

- **A（推荐，最小）**：`fq trim-qual` 单次处理一组阈值。参数扫描由外层
  `parallel`/shell 循环完成（与 anchr 现状一致），pgr 不内置。
- **B**：内置多阈值扫描（类似 `--qual 15,20,25 --len 50,60,70` 笛卡尔积）。
  与 anchr 的 `Q{qual}L{len}` 目录结构对应，但属"便利功能"，需求不足不先做。

> 遵循 AGENTS.md「简洁优先」，先实现 A。B 留作未来方向（§5）。

### 1.4 质量编码与 auto 检测

- `--quality-base 33|64|auto`，默认 `auto`。`auto` 逐条移植 BBDuk 的检测算法
  （BBTools-40.01 `stream/FASTQ.java` `testQuality`）。它比"min/max 阈值"式规则
  更稳：翻转阈值 87 保证不把高质量 +33（Q42-54）误判为 +64；长读强制 +33 保护
  ONT 数据。
  - 采样：第一个输入文件的前 N 条记录（或前 ~1 MB 质量字符，先到者）。BBDuk
    只采前 2 条记录（8 行），我们放大采样量以更稳，判定逻辑与其一致。
  - 初始假设 +33，逐字符判定：
    - 字符 > 87（`33 + QUAL_THRESH`，`QUAL_THRESH=54`）→ 翻转假设为 +64；
    - 当前假设 +64 且字符 < 59（Q < -5，Solexa 下限）→ 翻回 +33；
    - 假设 +33 且碱基为 N、质量字符为 64（`@`）或 66（`B`）→ 翻转假设为 +64
      （+64 数据中 N 常带 Q0/Q2）；
    - 采样中任一条读长 ≥ 200 bp → 强制 +33（长读的 +33 质量可超 Q41，会误触发
      +64）。
  - 两次翻转 → 放弃检测，按 +33 处理（与 BBDuk 一致：无法判定时默认 +33，
    不报错，只警告）。
- 无质量字符（FASTA 输入）→ 报错"FASTQ quality required"。
- 输出与输入同编码：只修剪、不改写质量值，`--quality-base` 只影响输入解析。
- Solexa（sickle `-t solexa`）按 +64 近似，不单独设值。
- 双端输入只在第一个文件上检测，另一文件沿用（同批数据编码一致）。

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
| 质量编码 | `--quality-base 33|64|auto`，默认 auto（规则见 §1.4）；Solexa 按 +64 近似；输出与输入同编码 |
| 分层 | 算法（滑窗/Mott）放 `libs/fq/`，`cmd_pgr/fq/` 仅做 clap 编排 |
| 依赖 | 用自研 `SeqReader`（`fmt/seq.rs`）读 FAFQ，不引入新依赖 |

## 4. 测试计划

- **单元**：滑窗与 Mott 对已知质量串的切点断言（两端质量高/低、中部低质量区、
  空序列、长度阈值边界）。
- **集成**：`tests/cli_fq_trim_q.rs`，覆盖单端/双端/singles、`--method` 两种、
  `-q`/`-l` 阈值、质量越界报错、输入输出同名报错。
- **质量编码检测**：普通 +33（字符 ≤74 不触发翻转）、+64（字符 >87 触发）、
  N+质量 `@`/`B` 触发、读长 ≥200 bp 强制 +33、两次翻转回退 +33、
  显式 `--quality-base` 覆盖、FASTA 输入报错。
- **对照**：与 `sickle` 同参数跑同一对文件，比对修剪后序列（滑窗语义应一致）。

## 5. 未来方向（暂不做）

- 内置多阈值参数扫描（§1.3 方案 B）。
- `--nextseq` polyG 变体默认开启（当前默认关，待证据）。
- `expected_errors`（Edgar 2015）整条读错误数过滤，作为独立 `fq` 功能。

## 6. 参考对照：BBTools BBDuk（2026-08 补充）

用户常用 BBTools，其 `bbduk.sh`（BBTools-40.01）的质量修剪相关选项与本文设计的对照如下。
原则：**只采纳真正低成本且对需求（anchr 的 sickle 语义）无害的选项；其余记入未来方向**。

| 主题 | BBDuk 做法 | 本文设计 | 采纳建议 |
|---|---|---|---|
| 修剪方向 | `qtrim=f\|r\|l\|rl\|w`：方向是显式枚举，默认 `f`（关） | 默认两端修剪，`--no-fiveprime` 禁 5' | 保留本文模型（对应 sickle `-x`；BBDuk 把"方向"和"算法"耦合在一个枚举里，不采用） |
| 修剪算法 | 默认（`qtrim=r/l/rl`）是 BWA 式"末尾需 ≥2（`minGoodInterval`）个连续达标碱基"（`TrimRead.testLeft/testRight`）；`qtrim=w` 是滑窗（固定窗口 4，只修 3'）；`optitrim` 才是累积最大（`testOptimal`，与我们的 mott 同族） | `--method sliding\|mott`，sliding 默认 | 保持本文双算法；BBDuk 三种算法均不等于 sickle 自适应窗口，仅作对照 |
| 质量阈值 | `trimq=6`，**支持浮点**（源码确认 `Float.parseFloat`，转错误率参与判定） | `-q` 浮点，默认 20 | 已采纳：允许浮点；默认 20 保持（anchr 用整数阈值扫描） |
| 长度过滤 | `minlength`（默认 10）、`mlf`（原始长度分数）、`maxlength` | `-l` 绝对长度 | 只保留 `-l`；mlf/maxlength 无需求 |
| 额外过滤 | `maq`（修剪后平均质量）、`mbq`（单碱基最低质量）、`maxns`（N 上限）、`mcb` | 无 | 本轮不引入；记入未来方向 |
| 质量编码 | `qin=auto`（33/64 自动检测，检测算法见 §1.4）、`qout=auto`（输出保持输入编码） | `--quality-base 33|64|auto`，默认 auto | 已采纳：逐条移植 BBDuk 检测算法（翻转判定 + 长读强制 +33），显式值兜底；输出与输入同编码（只修剪不改写质量） |
| polyG / NextSeq | 无 `--nextseq` 开关；用显式 `trimpolygleft/right`、`trimpolyg=N`（修剪 ≥N 的 polyG 前缀/尾），另有 `filterpolyg` | `--polyg-right N`（默认 0=关） | 已采纳：对齐 BBDuk `trimpolygright=N`，N=0 关 |
| 双端过滤语义 | `outs` = 幸存单端；`removeifeitherbad=t`（默认，源码确认 `shouldRemove`）：任一端 bad → pair 移出主输出，幸存端写 `outs`；`rieb=f` 才要求两端都 bad | 仅一端过 → 幸存端写 singles | 与 sickle 完全一致：默认（rieb=t）+ 可选 `outs` 即是我们 §2.4 的语义；`rieb=f` 是额外选项，不引入 |
| 输出顺序 | `ordered=f` 默认乱序（多线程） | 单线程天然有序 | 无参数 |
| 输出编码/折行/压缩 | `qout=auto`、`fastawrap=70`、`ziplevel=2` | 输出保持输入编码；FASTQ 不折行；压缩级别用默认 | 不引入对应参数 |

> 源码依据（BBTools-40.01）：`stream/FASTQ.java` `testQuality`（质量编码检测）、
> `parse/Parser.java` `parseTrim`/`parseTrimq`/`parsePoly`（qtrim/trimq/polyG 解析）、
> `shared/TrimRead.java` `testLeft`/`testRight`/`testRightWindow`/`testOptimal`
> （修剪算法）、`jgi/BBDuk.java` `shouldRemove` 与 outs 路由（双端语义）。

**结论**：已采纳 `-q` 允许浮点、质量编码 33/64/auto 默认 auto（BBDuk 检测算法）、
`--polyg-right N` 替代布尔开关、BBTools 式隐式 interleaved 双端输出；保持现状的是
方向/算法模型、`-l` 绝对长度、singles 语义；明确拒绝 maq/mbq/maxns/mlf/maxlength
等额外过滤与 `rieb=f` 选项。

## 7. 实现记录（2026-08）

**落地**：`pgr fq trim-qual` 已实现并验证。

- 代码：`src/libs/fq/trim.rs`（滑窗/Mott 算法、质量编码检测、单/双端编排）、
  `src/cmd_pgr/fq/trim_q.rs`（clap 薄壳）、`src/libs/fmt/seq.rs`（`SeqRecord`
  加 `Clone`，供 auto 检测的采样缓冲）。
- 与设计稿的差异：
  - `--no-fiveprime` 对 mott 同样有效：5' cutoff 置 0（cutadapt `-q 0,N` 语义），
    不再是"仅 sliding 有效"。
  - `--polyg-right N` 只修剪**连续** G 尾（BBDuk `trimpolygright` 的
    `maxNonPoly=0` 简化），中间有非 G 即不剪。
  - 质量校验范围 [0, 93]（对齐 BBDuk `MAX_CALLED_QUALITY` 上限）；越界报错
    含记录名与位置。
  - auto 检测采样放大为前 1000 条记录或 1 MB 质量字符（BBDuk 只采 2 条），
    判定逻辑与常量（`QUAL_THRESH=54`、长读 200 bp 强制 +33、N+`@`/`B`、
    两次翻转回退 +33）完全一致。
- 双端：`--outfile-2` 给出分离输出，否则 interleaved 写入 `-o`；幸存单端写
  `--outfile-single`（可选，不写则丢弃）；记录数不匹配警告并只处理公共前缀。
- 验证：lib 单元 13 个（滑窗/Mott 切点、编码检测、polyG、越界报错）+ 集成 14 个
  （单端/双端/交错/singles、mott、auto 与显式编码、gzip、错误路径）；
  `cargo clippy --all-targets -- -D warnings` clean；全量 lib 685 + 集成套件全绿。
- 吞吐 sanity（release，103 MB / 50 万条 100 bp 读）：trim-qual 0.23 s；同文件
  `fq to-fa` 0.13 s。修剪相对纯解析+写出的开销约 +0.1 s，无异常掉速。
- 基准测试：按约定本命令不引入 criterion 基准（全新命令、O(n) 算法、瓶颈在既有
  FAFQ 解析）；未来若对修剪算法做 SIMD 化，再按"先写基准"原则补。

---

*参考来源: [sickle.md](../references/sickle.md) | [cutadapt.md](../references/cutadapt.md) | [anchr trim.tera.sh](https://github.com/wang-q/anchr/blob/04f827afe37d5f40f12cd0602d54086cf8b0078c/templates/trim.tera.sh)*
