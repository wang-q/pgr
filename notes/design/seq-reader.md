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

## 4. BGZF 随机访问：inflate 瓶颈与方案搜索（2026-08-09）

### 现状基准（50 MB FASTA → 24 MB BGZF）

- 随机 1000 区间：**2.38 s**；perf 分布：inflate 57%、memset 30%、
  FASTA 解析 <2%；gzi 构建仅 8 ms。
- 同块密集访问（记录级缓存命中）：**6 ms** —— 证明 inflate 是纯重复
  劳动，命中缓存可数量级下降。
- 顺序 FASTA 的 6× 优化空间在 BGZF 上**不存在**：瓶颈是解压固有成本，
  自研解析无收益。

### noodles-bgzf 0.45.0 既有能力（源码核实）

- 默认后端 `flate2` + `zlib-rs`（`default-features = false`），已是 Rust
  侧较快的选择。
- `libdeflate` feature 存在（`["dep:libdeflater"]`）：`src/deflate.rs`
  `#[cfg]` 编译期切换 decode/encode，**零源码改动**。
- `MultithreadedReader`（io 线程 + inflate 线程池）实现了
  `crate::io::Seek`（`seek_to_virtual_position`/`seek_with_index`），但
  seek 内部 `get_mut()` 会 `pause()`（join 全部线程）+ `resume()`
  （重建线程池）——**随机 seek 每次重启管线，实际不可用**。
- `IndexedReader` = `Reader` + `.gzi`，seek = 底层 seek + 读新块 +
  inflate，**无块缓存**；每块新建 `DeflateDecoder`（窗口缓冲
  分配+清零，memset 30% 的主要来源）。

### 方案矩阵

| 方案 | 原理 | 单块收益 | 成本 | 适用 |
|---|---|---:|---|---|
| A. 块 LRU 缓存 | 同块命中不重复 inflate | 同块数量级；混合看命中率 | 中（pgr 层包一层，~200 行） | **loc 随机访问（热点）** |
| B. libdeflate feature | 换更快后端 | ~+18%（~650 vs ~550 MB/s） | 低（1 行 feature），但引入 C 依赖（cmake） | 全部 BGZF 读 |
| C. linflate | 纯 Rust 更快后端 | ~+27%（~700 MB/s） | 高（noodles 无支持，需 fork/自研块解析） | 全部 BGZF 读 |
| D. MultithreadedReader | 并行 inflate | 顺序读近核数倍 | 低（现成） | 顺序读；与"fa 单线程"约束冲突，且随机 seek 不可用 |
| E. Intel ISA-L | 最快 C 后端 | ~14×（7.9 GB/s） | 很高（C+cmake） | 大批量解压，单块 64 KB 启动开销吃收益 |

### 结论与建议

- **A 优先**：正好打在 loc 瓶颈（同块重复解压），零新依赖、保持单线程；
  但 0.45 `IndexedReader` 不可插拔，需 pgr 自研带缓存的 Bgzf 读取层
  （或 fork noodles-bgzf）。
- **B 作为低成本可选**：一个 feature 开关 +18%，代价是 C 依赖，需用户
  批准。
- C/D/E 记录在案，现阶段不推荐：C 工程量大、D 与单线程约束冲突且随机
  seek 不可用、E 依赖过重。

## 5. 方案 A 落地：CachedBgzfReader（2026-08-09）

### 实现

- `src/libs/bgzf.rs`：`CachedBgzfReader`，块级 LRU（key = 压缩偏移，
  value = 解压块 + 块大小），复用单个 `flate2::Decompress`（消掉每块
  新建解码器的分配/清零），`seek(uncompressed)` 走 `.gzi` 索引转虚拟
  位置，命中缓存不碰文件；空块（EOF 标记）跳过，文件尾当 EOF。
- 接入：`Input::Bgzf(CachedBgzfReader)`（默认缓存 16 块 = 1 MB），
  `paf/fasta.rs` 的 FastaStore 同步切换；`create_loc` 顺序建索引仍用
  noodles IndexedReader。
- 块解析复用 `fmt/fa.rs build_gzi_index` 的经验（BC 子字段 bsize）；
  CRC32 校验保留（flate2::Crc）。

### 基准（50 MB FASTA → BGZF，764 块，criterion）

| 实现 | cold 764 块读取 | warm 同块 1000 读 |
|---|---:|---:|
| noodles IndexedReader（现状） | 41.2 ms | 51.2 ms |
| CachedBgzfReader cap=1 | 37.1 ms | 16.7 ms |
| CachedBgzfReader cap=4 | 37.3 ms | 149 µs |
| CachedBgzfReader cap=16 | 37.0 ms | 149 µs |

- **warm 51.2 ms → 149 µs（343×）**：块缓存消除重复 inflate；
  cap=1 时 warm 只有 16.7 ms——warm 探针有 ~2% 跨块，单块缓存抖动
  （两块反复驱逐），cap≥4 稳定。
- **cold 41.2 → 37.0 ms（1.11×）**：缓存无命中收益，纯来自 inflater
  复用（noodles 每块新建 DeflateDecoder）。
- 集成测试：`cli_fa_index` 12 个 + `cli_paf_graph`/`cli_paf_stat`
  全过。

## 6. 候选 inflate 后端清单与系统测试计划（2026-08-09）

用户批准：全部候选装为 **dev-dependencies**（不进主线产物），系统对比
后再决定是否引入。

安装：`cargo add --dev linflate libdeflater isal-rs libz-ng-sys rust-htslib`

| crate | 版本 | 类型 | 宣称 | 备注 |
|---|---|---|---|---|
| linflate | 0.1.x | 纯 Rust | ~700 MB/s | 方案 C；full-buffer 匹配 BGZF 块；原项目已迁 nordisk/znippy |
| libdeflater | 1.25.x | C | ~650 MB/s | noodles `libdeflate` feature 底层 |
| isal-rs | 0.5.3 | C | ~7.9 GB/s | 注意 crates.io 另有只做 erasure code 的 `isa-l` |
| libz-ng-sys | 1.1.x | C | ~1 GB/s+ | zlib-ng SIMD |
| rust-htslib | 0.50 | C | — | 完整 BGZF 对照（多线程/索引） |
| flate2/zlib-rs | 1.1.9 | Rust | ~550 MB/s | 现状基线（pgr 在用） |
| zlib-rs 直接 | 0.6.x | Rust | ~550 MB/s | flate2 底层直调，少一层包装 |
| miniz_oxide | 0.8 | Rust | 慢 | 对照用（已在树中） |

注意事项：

1. C 后端需构建工具（cc / cmake / htslib vendored）。
2. **不要**通过 flate2 feature 开 zlib-ng/libz-ng：backend feature 全局
   合并，会切换整个依赖树（含 noodles 与 CachedBgzfReader）的后端；
   C 后端一律直接调 C API。
3. linflate 无 CRC 校验，需 crc32fast 自算（已在树中，正式使用需显式
   声明）。

测试矩阵（装好后执行）：

1. 微基准：单 64 KB deflate 块解压吞吐（GB/s）各后端对比。
2. 宏观基准：CachedBgzfReader 换后端 × 随机访问 cold/warm + 顺序读。
3. rust-htslib 整体方案对照（含并行解压选项）。

## 7. 后端系统测试结果（2026-08-09）

### 环境与实现

- 后端抽象：`src/libs/bgzf.rs` 增加 `BlockInflater` trait + 默认
  `Flate2Inflater`（flate2/zlib-rs，复用解码器）；`open_with_inflater`
  供基准注入。产品默认路径不变。
- dev-deps：linflate 0.1.11、libdeflater 1.25.2、isal-rs 0.5.3
  （+isal-sys）、libz-ng-sys 1.1.29、miniz_oxide 0.8。
- rust-htslib 1.0.1 放弃：hts-sys 需要 clang/bindgen 生成绑定，
  环境不支持；且其代表的是并行解压路线（已排除），价值有限。
- isal-sys 从源码构建 ISA-L 需要 autotools（autoreconf）与 nasm，
  本机已补装（apt install autoconf automake libtool nasm）。
- isal/libz-ng 的基准包装需要少量 unsafe（sys crate 本质），仅限
  bench；libz-ng 复用 strm（inflateReset）在本项目场景报
  Z_STREAM_ERROR，改为每次 init/end（结果仍无优势，未深挖）。

### 微基准：单 64 KB deflate 块解压（DNA 数据，criterion）

| 后端 | 时间/块 | 吞吐 |
|---|---:|---:|
| **libdeflater** | **32.9 µs** | **1.99 GB/s** |
| flate2/zlib-rs reuse（现状） | 43.7 µs | 1.50 GB/s |
| flate2/zlib-rs fresh（noodles） | 44.4 µs | 1.48 GB/s |
| libz-ng（zlib-ng C） | 43.1 µs | 1.52 GB/s |
| isal stateless（ISA-L） | 44.2 µs | 1.48 GB/s |
| linflate | 55.5 µs | 1.18 GB/s |
| miniz_oxide | 55.5 µs | 1.18 GB/s |

### 宏观基准：CachedBgzfReader 换后端（cap=16，50 MB FASTA → BGZF）

| 后端 | cold 764 块 | warm 同块 1000 读 |
|---|---:|---:|
| IndexedReader（noodles 现状） | 40.1 ms | 50.9 ms |
| **libdeflater** | **28.8 ms** | **125.2 µs** |
| flate2/zlib-rs（默认） | 36.6 ms | 147.4 µs |
| isal | 36.1 ms | 148.4 µs |
| libz-ng（fresh init/end） | 36.5 ms | 157.8 µs |
| linflate | 45.2 ms | 173.6 µs |

### 结论

- **libdeflater 全面领先**：微基准 +25%、宏观 cold +21%、warm +15%，
  且 C 后端复用轻量。**方案 B 验证通过，是唯一值得产品化的后端**。
- **linflate（方案 C）证伪**：本机/本数据比 zlib-rs 慢 27%（宣传
  ~700 MB/s 未兑现，可能与数据形态/无 AVX512 有关）；宏观 cold 甚至
  比 noodles IndexedReader 还慢。不再考虑。
- isal / libz-ng 与 zlib-rs 持平：宣传的 7.9 GB/s（ISA-L）在单块
  64 KB 场景被启动开销吃掉，无引入价值。

### 落地（2026-08-09，用户批准）

- libdeflater 1.25.2 提升为正式依赖（C 编译，cc 即可）。
- **不引入 feature 开关**（用户裁定：不要开 features，直接切换）：
  `CachedBgzfReader::open` 默认后端 = `LibdeflaterInflater`；
  `Flate2Inflater` 保留供基准对比（`open_with_inflater` 注入）。
- 切换后宏观：cold 33.6 → 29.0 ms、warm 142 → 136 µs
  （与 flate2 变体同轮对比；数据有系统波动，libdeflater 稳定领先）。
- 集成测试 `cli_fa_index`（12）/`cli_paf_graph`/`cli_paf_stat` 全过；
  fmt/clippy clean。
- 其余 dev-deps（linflate/isal/libz-ng/miniz_oxide）保留，仅用于
  基准复现；rust-htslib 已从 dev-deps 移除。

## 8. 去除 noodles-bgzf：自研 BGZF 全套（2026-08-09）

用户裁定：**去除 noodles-bgzf 依赖，从头实现**。`src/` 里
noodles_bgzf 引用清零（仅 benches/tests 的 dev 对照保留）；
`Cargo.toml` 移除 noodles-bgzf 直接依赖与 noodles 主 crate 的
`bgzf` feature。

### 自研组件（`src/libs/bgzf/`）

- `index.rs`：`GziIndex`（bgzip 兼容 .gzi 读写 + 二分 query）。
- `mod.rs`：`VirtualPos`（48+16 位打包）、通用 gzip 头解析
  （FEXTRA/FNAME/FCOMMENT/FHCRC + BC 子字段）、整块
  `gzip_compress`/`gzip_decompress`（libdeflater + crc32fast，
  ISIZE 预分配 + bomb 上限）、`CachedBgzfReader`（随机访问块缓存，
  现支持无索引的虚拟位置 seek + `BufRead`）。
- `writer.rs`：`BgzfWriter`（单线程）+ `ParallelBgzfWriter`
  （worker 池压缩、按序输出；块上限 0xff00 留 stored 兜底余量）。
- `reader.rs`：`GzReader`（流式多成员 gzip，zlib-ng raw inflate +
  手写 gzip 头/trailer/CRC 校验，替代 MultiGzDecoder）+
  `ParallelBgzfReader`（读线程解析块 + ≤4 worker libdeflater 解压 +
  按序输出，替代 bgzf::io::Reader 的并行版）。

### 关键 bug 与教训

1. **z_stream 地址移动**：zlib 的 inflate_state 保存 init 时的 strm
   地址，Rust struct 移动（入 Box/struct）后 `state->strm != strm` →
   Z_STREAM_ERROR。**修复：`Box<z_stream>` 固定地址**。这解释了之前
   基准中 libz-ng 复用失败的现象。
2. **gzip 模式（window_bits=31）在此 zlib-ng 版本不可用**（首次
   inflate 即 Z_STREAM_ERROR），改用手写 gzip 头 + raw inflate(-15) +
   手动 CRC/ISIZE 校验。
3. **pbit 段边界既有 bug**：collection 段长用 `footer_start`（文件尾）
   而非 `paf_data_offset`，把 PAF recovery 数据并入 collection；
   flate2 GzDecoder 只解第一个成员而侥幸通过，自研
   gzip_decompress（尾部 ISIZE 预分配）暴露。**修复：
   `collection_end = paf_data_offset.min(footer_start)`**。
4. CachedBgzfReader 初始 current=None 时顺序读直接 EOF（随机访问先
   seek 未暴露）——初始化为从块 0 开始。

### 依赖变更

- 新增正式：libdeflater（已有）、crc32fast、libz-ng-sys（流式 gzip）。
- 移除正式：noodles-bgzf、flate2（pgr 直接依赖；仅 dev 测试/基准用）。
- flate2 使用点迁移：pbit 4 处解压 + 3 处压缩 → gzip_compress/
  gzip_decompress（保留 bomb 上限）；sd/search_lastz、fq is_fq →
  自研；io.rs MultiGzDecoder → GzReader/ParallelBgzfReader。

### 性能（50 MB FASTA → BGZF，release）

| 顺序读实现 | 耗时 |
|---|---:|
| GzReader 单线程（zlib-ng 流式） | 42.3 ms |
| ParallelBgzfReader 1 worker（libdeflater） | 33.0 ms |
| ParallelBgzfReader 2 workers | 17.7 ms |
| **ParallelBgzfReader 4 workers** | **9.5 ms** |

- 4 线程 vs 单线程流式 **4.5×**；`pgr fa size` 端到端
  148 ms → **20 ms（7.4×，含解析与输出）**。
- 随机访问（CachedBgzfReader + libdeflater）与写入侧迁移后，
  `cli_fa_index`/`cli_paf_bgzf` 等全部集成测试通过。

### 剩余

- dev-deps 的 noodles-bgzf 仅作基准对照（IndexedReader）；正式代码
  已完全自研。

解读：

- **多行拼接是 noodles 慢的主因**：FASTA 80 bp 多行时每记录逐行 append
  （34 ms）；FASTQ 单行时 noodles 也快（7.8 ms）——kseq 式对两者都
  ~4.5 ms。
- 原型已含：缓冲复用（kstring_t 式）、memchr 批量行读、FAFQ 统一
  （`> ` 或 `@` 自动检测）、质量读取；未含完整边界（CRLF/空行/非法
  定义）与 SIMD 分隔符搜索（tva 式，可再降）。
- **决策依据就绪**：若消费方接受借用式/缓冲复用 API，自研 FAFQ reader
  的收益上限明确（FASTA 7.6×、FASTQ 1.8×，内存读取层面）。
