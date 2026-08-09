# FAFQ（FASTA/FASTQ）读取与 noodles 切换决策（2026-08-09）

> 评估是否用自研 **FAFQ 统一顺序读取**替代 noodles 的 fasta/fastq 两个
> reader（兼顾 FASTA 与 FASTQ，借鉴 seqtk kseq.h），以及 noodles 依赖的
> 边界。基准与 BGZF 补充清单见下文。

## 1. noodles 在 pgr 的用途

`Cargo.toml`：noodles 0.104（core/fasta/fastq/bgzf/gff）。实际使用：

| 用途 | 位置 |
|---|---|
| FASTA 顺序读取/写入 | `libs/fmt/fa.rs`（`fmt::fa::reader` → `fasta::io::Reader`） |
| FASTQ 读取 | `libs/fmt/fq.rs` |
| BGZF 随机访问 | `libs/loc.rs`（FastaStore，`bgzf::io::IndexedReader` + `.loc`/`.gzi`） |
| GFF | gff 命令 |

## 2. 读取性能基准（50 MB 内存 FASTA，`benches/fa_read_benchmark.rs`）

| 实现 | 耗时 | 说明 |
|---|---:|---|
| noodles 完整记录（`records()`） | 34.2 ms | 序列行拼接 + 定义解析 + Record 分配 |
| noodles 底层读序列（无构造） | 12.4 ms | 接近内存带宽（I/O 主导） |
| 自研 naive 逐行 | 18.0 ms | 实现低效（反证：自研要 memchr/批量读才快） |

**记录构造占读取 CPU 开销的 ~64%**（34.2 − 12.4 = 22 ms/50 MB）。磁盘
场景（fa size 50 MB 纯文本 87 ms）里还有页清零（memset 54%）+ 文件 I/O；
gz 场景 inflate 主导，记录构造占比被稀释。

## 3. noodles_bgzf 功能不完整，pgr 后加的清单

noodles_bgzf 只提供 IndexedReader（读 .gzi + 虚拟偏移）与 multithreaded_writer
（压缩），**不生成 .gzi 索引**。pgr 补充：

| 功能 | 位置 | 内容 |
|---|---|---|
| `.gzi` 索引生成 | `libs/fmt/fa.rs` `build_gzi_index` | 手写 BGZF 块解析：12B 头、FLG/FEXTRA、BC 子块拿 bsize、块尾 ISIZE、构建偏移对 |
| BGZF 头部检测 | `libs/io.rs` `is_bgzf` | `1f 8b 08 04` + XLEN=6 + "BC" + SLEN=2 |
| `fa gz` 块读写 | `cmd_pgr/fa/gz.rs` | 64KB 块手动读写（绕开 std::io::copy 8KB，降低 multithreaded_writer channel 开销） |

## 4. 切换决策的边界

- **BGZF 随机访问是独立基础设施**，与 FAFQ 顺序读取两条路：loc 走
  noodles_bgzf + pgr gzi；`fmt/fa.rs`/`fmt/fq.rs` 顺序读走 noodles。
  切换 FAFQ 顺序读取**不碰 BGZF/gff**——部分切换可行。
- **"自补功能"是现成先例**：gzi 索引证明 noodles 覆盖不了的可以自己补；
  自研 FASTA reader 的边界（CRLF/空行/非法定义）同理，但需测试覆盖到
  noodles 的成熟度。
- **候选方案**：`fmt/fa.rs` 加一条**借用式零拷贝 reader**（行缓冲复用，
  `sequence()` 返回借用 slice），保留 noodles 兜底；按"先基准后修改"
  原则先做原型 benchmark。

## 5. seqtk kseq.h 参考（2026-08-09）

seqtk-1.5 的 `kseq.h` 是经典 FASTA/FASTQ 读取器（2008 年起，被大量工具
沿用）。

### 机制

- **kstream_t**：自实现 16 KB 缓冲流（read 回调），绕过 std::io；
- **ks_getuntil2**：`memchr` 批量找分隔符（KS_SEP_LINE 用 `memchr('\n')`），
  行读取是批量不是逐字符；
- **kseq_read**：`last_char` 状态机（记住已读到的下一个头字符，避免重复
  读）+ `is_fastq` 标记。

### FAFQ 兼容（kseq 与 noodles 的关键差异）

单个 `kseq_read` 同时处理两种格式：跳到第一个 `>` **或** `@` 开头，
统一读 `name`/`comment`/`seq`；读到 `+` 则继续读质量行（FASTQ），否则
直接返回（FASTA）。`qual.l > 0` 即 FASTQ。noodles 则 `fasta::Reader` 与
`fastq::Reader` 分离，格式必须预先知道，不能自动检测。

### 缓冲复用（解释 noodles 记录构造 64% 的一个来源）

kseq 的 `name/comment/seq/qual` 是 `kstring_t`（**保留容量 m**），每次
`kseq_read` 只 `l=0` 重置、不重新分配，仅在增长时 realloc。noodles
`records::next` 每记录 `Vec::new()` + `read_sequence` **从 0 增长
（多次 realloc + 复制）**。这是自研 reader 必须学的：**缓冲复用**。

### 对自研 reader 的借鉴清单

1. 缓冲复用（kstring_t 式保留容量，避免每记录 realloc）；
2. `memchr` 批量行读（noodles 已有，自研照做）；
3. FAFQ 统一 reader（可选：pgr 有 fa/fq 通用场景，如 fq to-fa）；
4. kseq 也非零拷贝（seq.s 持有缓冲）——借用式 slice 是更进一步，但
   消费方要按借用模式改造。

## 6. tva 的 SIMD 分隔符搜索参考（2026-08-09）

`~/Scripts/tva/src/libs/tsv/simd/`（另一项目）有现成的 **SSE2/NEON
128-bit 分隔符搜索器**：单趟找 `\t` + `\n`（`_mm_cmpeq_epi8` 判等 →
`_mm_or_si128` 合并 → `_mm_movemask_epi8` 掩码 → `trailing_zeros` 遍历位），
CR 单独查 newline 前一字节。基准 ~6.5 GiB/s（比两遍法快 ~6.7×）。

**与 kseq 的类似点**：都是在缓冲里快速找分隔符——kseq 用 libc `memchr`
（本身已 SIMD），tva 自写。对自研 FASTA reader 的意义：

- FASTA 解析需要找 `\n`（行尾）和 `>`（定义行），可用同样的 SIMD 搜索
  原语（多分隔符：`>` + `\n`，或 `@` + `+` 兼容 FASTQ）；
- tva 的实现模式与我们的 SIMD 原则一致（SSE2/NEON 128-bit、架构门控、
  trait 抽象），是现成样板；`simd-minimizers` 的 packed_seq 同理。
- 结合 kseq 的缓冲复用 + tva 的 SIMD 搜索 = 自研 reader 的两个核心
  优化原语；若 tva 代码可共享（同作者/同仓库族），直接移植模式。

## 7. FAFQ 读取基准与实现（2026-08-09）

### 实现：`libs/fmt/seq.rs`

`SeqReader`/`SeqRecord`：kseq 式缓冲复用 + `memchr`/`memchr2` 批量行读 +
FAFQ 自动检测（`>`/`@` 头、`+` 质量）+ 边界处理（CRLF、空行、质量长度
校验）；`SeqReader::new` 走 `io::reader`（支持 stdin/gz），
`from_reader` 支持借用/owned 缓冲。`memchr` 2.8.3 已入正式依赖。

### bstr 引入（2026-08-09）

`SeqRecord` 的 `name`/`comment` 用 `BString`（bstr 1.13，传递依赖已有、
现入直接依赖），`name()`/`comment()` 返回 `&BStr`——**名称是字节字符串，
读取层不强制 UTF-8**（新增 `non_utf8_name_is_byte_clean` 测试：非 UTF-8
名称可正常读取，消费方自行决定解码）。`description()` 保持 `Option<&[u8]>`
兼容写回/签名路径。`io.rs` 的 `read_lines`/`read_names` 评估后**不改**：
处理的是用户文本列表（UTF-8 语义合理），bstr 无实质收益。

### 基准（50 MB 合成数据）

50 MB 合成数据（FASTA 80 bp 多行；FASTQ 单行——noodles_fastq
`records()` 只读单行序列/质量，多行 FASTQ 不支持，是其不完整性的又一例）：

| 实现 | FASTA | FASTQ |
|---|---:|---:|
| noodles `records()` | 34.2 ms | 7.8 ms |
| kseq 式原型（缓冲复用 + memchr） | **4.5 ms** | **4.4 ms** |
| **产品级 `SeqReader`** | **5.6 ms** | **4.3 ms** |
| 产品级提升 vs noodles | **~6.1×** | ~1.75× |

### 已接入命令与端到端（2026-08-09）

- `fa size`：纯文本 50 MB 87 → **25 ms（3.5×）**；gz 207 → **148 ms（28%）**
- `fa count`：纯文本 25 ms；gz 114 ms（读取型，同 size 级提升）
- `fa masked`：输出型（1234 万行），输出 I/O 主导，读取优化被淹没
- `fq to-fa`：**多行 FASTQ 50 MB × 2（序列+质量）54 ms**——noodles_fastq
  `records()` 不支持多行（只读单行），这是 SeqReader 独有能力

### 全部顺序读取已切换（2026-08-09 晚）

`rg 'fmt::fa::reader|records()'` 清理完毕（除 `fmt/fa.rs` 定义、`loc.rs`
BGZF 随机访问、`plot/histogram.rs` CSV 外零残留）：

- **fa 命令**：size / masked / count / one / rc / order / six_frame /
  replace / some / n50 / dedup / mask / filter / to_2bit / split（15 个）
- **fq 命令**：`to-fa`、`interleave`（含 `interleave_read` 手动双循环）
- **库层**：`pgi::build::read_fasta`、`hv`×2、`hash`×3、`pbit/compressor`、
  `pl/repeat`×2、`alignment/msa`（外部 aligner 输出）
- 写回 noodles 记录的命令经 `fmt::fa::new_record_with_desc`（新增，
  `bstr` 1.13 入依赖）

### 行为差异（kseq 式宽容 vs noodles 严格）

- noodles_fastq：多行 FASTQ 报 "invalid description prefix"；缺序列行的
  畸形记录报错。SeqReader：两者都宽容（多行拼接、空记录），零 panic。
  `command_fq_to_fa_ucsc_malformed_no_panic` 已更新为预期 success。
- noodles_fasta：尾部 `>c`（定义后 EOF）返回空记录——SeqReader 一致。

### noodles 去留

**保留**：全部顺序读取已走 SeqReader；noodles 仅剩 BGZF 随机访问
（`loc.rs`）、GFF、写入路径（`fmt/fa.rs` writer）。顺序读取侧 noodles
路径已无生产消费方。

解读：

- **多行拼接是 noodles 慢的主因**：FASTA 80 bp 多行时每记录逐行 append
  （34 ms）；FASTQ 单行时 noodles 也快（7.8 ms）——kseq 式对两者都
  ~4.5 ms。
- 原型已含：缓冲复用（kstring_t 式）、memchr 批量行读、FAFQ 统一
  （`> ` 或 `@` 自动检测）、质量读取；未含完整边界（CRLF/空行/非法
  定义）与 SIMD 分隔符搜索（tva 式，可再降）。
- **决策依据就绪**：若消费方接受借用式/缓冲复用 API，自研 FAFQ reader
  的收益上限明确（FASTA 7.6×、FASTQ 1.8×，内存读取层面）。
