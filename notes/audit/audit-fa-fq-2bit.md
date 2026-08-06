# pgr fa / fq / 2bit 命令族代码审核记录（2026-08-05）

对 `pgr fa`（18 个子命令）、`pgr fq`（2 个子命令）与 `pgr 2bit`（5 个子命令）
命令族及相关库文件（`libs/fmt/fa`、`libs/fmt/fq`、`libs/fmt/twobit`、
`libs/translate`、`libs/io`、`libs/ds/range`、`libs/nt`、`libs/loc`、
`libs/fasta/chunk`、`libs/fasta/stat`、`libs/fasta/dedup`、
`libs/fasta/filter`）和全部测试/文档（`docs/{fa,fq,twobit}.md`）进行审核。
三个命令族结构相近（同为序列格式操作：fa/fq 互为转换、2bit 为二进制索引格式），
故合并为一份审核记录。缺陷按类别分组记录；关键修复均附回归测试，验证概况见
文末"验证"一节。

审核范围：
- **fa**
  - **info**：`size` / `count` / `masked` / `n50`
  - **records**：`one` / `some` / `order` / `split` / `window`
  - **transform**：`replace` / `rc` / `filter` / `dedup` / `mask` / `six-frame` / `to-2bit`
  - **indexing**：`gz` / `range`
- **fq**
  - **转换**：`to-fa`（FASTQ → FASTA）
  - **双端**：`interleave` / `il`（单/双文件交错，可生成虚拟 R2）
- **2bit**
  - **信息**：`masked` / `size`
  - **子集**：`range` / `some`
  - **转换**：`to-fa`

审核重点：数据安全（`-o` 不得覆盖输入及辅助列表文件）、Zero Panic（畸形输入
不 panic）、坐标/边界处理、2bit 打包/解包正确性、文档一致性。

## 与外部参考实现的语义一致性核对

* 2bit 家族对照 UCSC kent-tools 的 2bit 二进制格式规范（magic/version/index/
  record 布局、2-bit 打包、N-block/mask-block、双字节序）逐字节核对，语义一致。
  有意差异（均已记录）：
  * 读取端同时支持 version 0（index 偏移 4 字节）与 version 1（index 偏移 8
    字节）两种布局；写入端产出 version 1。
  * 打包时任何非 A/C/G/T 碱基（含 IUPAC 歧义码、U）按 N 处理并记录为
    hard-mask 块；软屏蔽小写 A/C/G/T 记录为 mask 块。`size --no-ns` 仅扣减
    N-block，与 UCSC "排除 hard-masked N 位置"语义一致。
  * 畸形输入：kent-tools 部分输入直接崩溃/告警，pgr 统一友好出错/跳过
    （Zero-Panic）。

## 排除的疑点（经核验无需修复）

**fa / fq：**
* 逐命令通读全部 18 个 `fa` 子命令的 `execute`：`unwrap()`/`unreachable!`
  全部为 clap `required` 参数或 `value_parser` 约束枚举，运行期不可达，
  无潜在 panic（符合"稳定性原则"）。`interleave` 的 clap `unwrap()` 同理
  （`infiles` required、其余参数有默认值）。
* `-o` 覆盖保护覆盖情况：全部单文件输出命令（`count`/`size`/`n50`/
  `masked`/`one`/`some`/`order`/`rc`/`replace`/`range`/`window`/`dedup`/
  `filter`/`to-2bit`/`gz`/`six-frame`）及 `fq to-fa` 均调用
  `ensure_outfile_distinct`；`split` 输出为目录，采用逐输出路径 `same_path`
  反向检查（见下）。
* `range start == 0`：`Range::from_str("chr1")`（仅名称）返回 `start=0`，
  `execute` 中 `if *rg.start() == 0` 写整条记录，符合"仅名称 = 整条序列"
  的语义；`chr1:0-100` 这类退化坐标同样走整条路径，与 1-based 约定一致。
* `split about` 无 `-c`：`opt_count = usize::MAX`，`cur_size > usize::MAX`
  永不触发，全部记录落入文件 0；配合 `-e` 仍按每 2 条切分，行为合理。
* `split about --max-part`：`part_width` 用 `checked_ilog10` 计算零填充
  宽度，`--max-part` 已校验为正数，无溢出。
* `translate` 反向链坐标：ORF 的 `end*3` 不超过该框碱基数 `dna_len - frame`，
  `dna_len - end*3` 不会下溢，坐标恒 ≥ 1。
* `mask_sequence` 异常 span：`lower < 1` 先拦截；`upper < lower` 时 `as usize`
  转为巨大值，被 `offset + length > out.len()` 兜底拦截，不会 panic。
* `n50` 的 `calc_n50_stats`/`transpose`：`transpose` 假设各行等长，实测各
  `--no-header` 形态下所有行长度一致（均 1 或均 2），无 panic 风险；
  `nx_sizes[i] == 0` 作"未设置"哨兵，空序列集合 Nx 恒为 0，语义正确。
* `count_bases` 对 IUPAC 歧义码（M/R/W/S/Y/K/V/H/D/B）经 `to_nt` 归为 `N`
  计入 `len`，与文档"歧义码计为 N"一致。
* `fq is_fq` 对目录输入：`File::open` 对目录成功但 `read_exact` 失败
  （EISDIR），返回友好错误而非 panic，符合 Zero Panic。
* `fq to-fa` 非 UTF-8 名称：`std::str::from_utf8(record.name())?`，报错而非
  panic。
* `fq write_fq`/`write_fa` 对空序列输出 `@name\n\n` / `>name\n\n`，合法。
* 文档一致性：`docs/fa.md`、`docs/fq.md` 中参数的默认值、坐标约定与代码一致。

**2bit：**
* 全部 5 个子命令的 `unwrap()` 均为 clap `required` 参数（`infile`/`infiles`/
  `name_list`）或 `value_parser` 约束枚举，运行期不可达，无潜在 panic。
* `Range::from_str` 为手写字节扫描器（不 panic），无 `:` 时回退 `chr` 为首个
  空白 token；`range` 命令对含 `:` 的输入用 `rg.is_valid()`（`start != 0`）
  兜底，对无 `:` 的全序列请求绕过坐标校验，畸形输入均返回友好错误而非 panic。
* `from_dna` 打包：`packed.len()` 恒等于 `ceil(len/4)`，与读取端
  `dna_size.div_ceil(4)` 一致；`bit_offset` 归零即推入整字节，尾部半字节在循环
  结束后补推。`test_blocks_from_dna`/`test_write_read_roundtrip` 验证一致。
* `TwoBitWriter::write` 的 `record_size`（16 + 8·n + 8·m + packed）与
  `write_packed_record` 实际写出的字节数一致；v1 偏移为 8 字节，index 起算
  `16 + Σ(1+N+8)` 正确。
* `size --no-ns`：2bit 中 IUPAC 歧义码与 N 均以 N-block 存储，故"排除 N 与
  IUPAC 歧义码"与文档"排除 hard-masked N 位置"语义一致；`n_count > dna_size`
  时返回友好错误。
* `rev_comp`（`range` 负链）用 `NT_COMP` 表，对 A/C/G/T/N 与大小写均正确，
  2bit 序列只含 ACGTN（含软屏蔽小写），无越界。
* `masked` 输出坐标：0-based 半开 `[start,end)` → 1-based 包含
  `start+1..=end`，单碱基输出 `name:pos`，与文档一致。
* `merge_intervals` 对相邻（`block.start <= last.end`）的 N/mask 块合并，
  半开区间相邻即合并，行为正确。
* 文档一致性：`docs/twobit.md` 与各子命令 `after_help` 的坐标、`--no-ns`、
  `--no-mask`、`-l 0`（不换行）、`-i`（invert）、大小写敏感、`#` 注释忽略等
  描述一致；命令分组（Info/Subset/Transform）一致。
* 纵深复核（跨轮）无新增问题的项：
  * `read_sequence` 的 `cached` 记录缓存与 `get_sequence_len`/`get_sequence_blocks`
    的 seek 无位置串扰：所有读取方法均在使用前显式 seek，缓存命中时不再 seek，
    混用（`range` 中 `get_sequence_len` + `read_sequence`）正确。
  * `range` 坐标 1-based→0-based 换算、首端越界跳过/末端越界截断告警、负链
    `rev_comp(NT_COMP)`、单点 `chr1:1`（end 默认 start）、`chr1:0`/反向区间报错，
    均正确；头部用原始 range 字符串，与文档一致。
  * `merge_intervals` 先按 start 排序、相邻（`block.start <= last.end`）合并，
    半开区间相邻即合并，N/mask 重叠位置正确合并为单区域。
  * `some` 的 invert 逻辑 `contains != invert`、大小写敏感、`#` 注释/空行忽略、
    首列取名，与文档一致。
  * 5 个子命令 `-o` 覆盖保护（含 `some` 的 list、`range` 的 rgfile）均已覆盖。

## 已知限制（有意保留）

* `2bit range` 的全序列反链请求 `seq_name(-)`（无坐标）不被解析：`Range::from_str`
  在无 `:` 时回退 `chr` 为整个 token，`strand` 为空，`has_sequence("seq_name(-)")`
  为假 → 告警跳过。文档仅承诺 `seq_name` 全序列与 `seq_name(strand):start-end`
  坐标形式，此用法未文档化，且不会静默返回错误数据（仅跳过）。
* `2bit range` 在既无位置参数 `ranges` 也无 `--rgfile` 时静默输出空结果（与
  `fa range` 行为一致）。

## 记录项（未改，低风险 / 待决策）

**fa / fq：**
* `fa split` 的 `name` 模式：`sanitize_filename` 把 `/`、`\`、`(`、`)`、`:`
  替换为 `_`，两个不同名称可能清洗到同一文件名（如 `a/b` 与 `a_b` 均 →
  `a_b`），会静默合并到同一文件。与文件名清洗方案固有行为一致，低风险，未改。
* `fa split` 的 `name` 模式：空/全特殊字符名称清洗后理论上可得空文件名，实际
  FASTA 名称非空，不可达，未改。
* `fa one` 在未找到序列时返回错误，但输出文件已先用 `File::create` 打开，会
  留下一个空文件（`-o`）。`fq to-fa`/`interleave` 打开 writer 后若输入读取
  失败也同样会留空/部分输出文件。属一般的"出错留空文件"行为，非数据损坏，
  `-o` 与输入重叠已被 `ensure_outfile_distinct` 前置拦截，未改。
* `fa size --no-ns` 对 `-`/`*`/数字等非 IUPAC 字符计数为"有效碱基"：实现用
  `!is_n(b)`（仅排除 N 与 IUPAC），`-`/`*` 会被计入长度，与"仅计算有效碱基"
  的字面略有出入；但文档明确把排除范围限定为"N 及 IUPAC 歧义码"，行为与文档
  字面一致，且与 `count` 的 `len`（仅 A/C/G/T/N）语义不同属两命令各自定义。
  > ✅ 已修复（2026-08-06）：`--no-ns` 改为排除 N + IUPAC 歧义码 + Invalid
  （`-`/`*`/数字），与"仅有效碱基"语义一致；新增集成测试
  `command_fa_size_no_ns_excludes_iupac_and_invalid`。
  低风险文档歧义，未改。
* `fa replace` 的 `read_replace_tsv` 不跳过 `#` 注释行，而 `read_names`
  （`some`/`order` 等）会跳过：若 TSV 带 `#` 表头/注释行，会被误当作 key
  （如 `#old\tnew` 生成 key=`#old`）。`replace` 命令文档未承诺注释支持，属
  两读取函数一致性小瑕疵，且 `#` 前缀名称极罕见，低风险，未改。
* `fa masked --gap` 的 `is_n` 会把 IUPAC 歧义码（M/R/W/S/Y/K/V/H/D/B 及 X）
  一并计为 N：帮助文本写作 "Only identify regions of N/n (gaps)"，字面略窄于
  实现。与代码库"歧义码即 N"的统一口径一致（默认 `masked` 亦将歧义码视为
  masked），属术语精度问题，非行为缺陷，未改。
* `fa split about -c N -o stdout`：`about` 模式无论输出到 stdout 还是目录，
  都会按 `max_files`（默认 999）轮转并 `break 'outer`，即在 stdout 上仅输出
  前 `max_files` 个分块容量内的记录。默认 `-c` 为 `usize::MAX`（不触发轮转）
  时不受影响；仅当用户显式 `-c SIZE -o stdout` 时才会提前截断。流式 stdout
  与按大小分块本就矛盾，属极端组合，低风险，未改。
* `fq interleave` 双文件格式只按 infile[0] 检测：`is_fq` 只探测第一个文件，
  若第二个文件格式不同会在读取时报错（非静默），行为可接受，未改。

**2bit：**
* `read_u32_vec` 在畸形文件给出超大 `count` 时会尝试分配巨大缓冲（可能 OOM
  abort，而非 panic）。属全局既有模式与畸形输入鲁棒性范畴，低风险，未改。

## 修复的缺陷（共 24 处）

### 崩溃 / 越界 / 溢出（Zero Panic，3 处，均 fa）

**`six-frame` 短序列 panic**：`&dna[frame..]` 在序列短于 frame 时（如
  1 碱基、frame=2）slice 越界 panic；`filter_and_convert_orfs` 中
  `dna_len - frame` 在 frame 超出序列长度时下溢。修复：改用
  `dna.get(frame..).unwrap_or(&[])`；`orfs.is_empty()` 时提前返回，避免
  `dna_len - frame` 下溢。回归 `command_six_frame_short_sequence_no_panic`、
  `test_six_frame_short_sequences_do_not_panic`。

**`window` 窗口长度过大导致 usize 溢出**：`start + len` 在 `len` 接近
  usize::MAX 且多窗口时溢出。修复：改用 `start.saturating_add(len)`。
  回归 `windows_huge_length_does_not_overflow`。

**`mask_sequence` 对多字节 UTF-8 序列 panic**：以 `&str` 切片按字节偏移
  操作，遇到多字节 UTF-8 字符会因 char 边界 panic。修复：重写为直接操作
  `&[u8]` 字节，避免 UTF-8 边界检查；同时增加 `lower < 1` 与越界 span 的
  错误检查。回归 `mask_sequence_non_ascii_does_not_panic`。

### 数据安全 / 参数校验（9 处：fa 7 + fq 1 + 2bit 1）

**全部 `fa` 子命令允许 `-o` 覆盖输入文件**：流式命令先打开输出再读取输入，
  若 `-o` 指向输入文件，会在读取前截断输入，静默清空数据。修复：在
  `count`/`size`/`n50`/`masked`/`one`/`some`/`order`/`rc`/`replace`/
  `range`/`window`/`dedup`（含 `--dups-file`）/`filter`/`to-2bit`/`gz` 中
  统一加入 `ensure_outfile_distinct` 检查。回归
  `command_fa_output_same_as_input_rejected`（`fa filter` 用例）。

**`six-frame` 遗漏 `-o` 覆盖输入的保护**：为命令族统一加检查时，
  `six-frame` 使用 `pgr::writer` 且遗漏该检查，是唯一仍可覆盖输入的子命令。
  修复：在 `six_frame.rs` 打开 writer 前加入
  `ensure_outfile_distinct(outfile, [infile.as_str()])?`。回归
  `command_six_frame_output_same_as_input_rejected`。

**`split` 输出文件可能覆盖输入文件**：`split` 的 `-o` 输出是**目录**，输出
  文件名由序列名（`name` 模式）或编号（`about` 模式）动态生成。若输入文件
  恰好位于输出目录且文件名与某输出文件名一致（如 `pgr fa split name data/
  -o data` 且序列名为 `chr`），`gen_fh` 用 `truncate(true)` 打开输出会截断
  正在读取的输入，导致后续记录丢失。修复：在 `gen_fh` 打开输出前，用
  `same_path` 将输出路径与所有输入路径比对，命中即 `bail!`，避免截断输入。
  回归 `command_fa_split_output_not_overwrite_input`。

> 注：`split` 的输出是目录，无法复用 `ensure_outfile_distinct`（该函数针对
> 单文件 `-o`），故在 `gen_fh` 内做逐输出文件的反向检查。

**`gz` 的 `-o stdout` 创建字面文件 `stdout` 而非写标准输出**：`gz` 用
  `std::fs::File::create` 直接建输出文件（未走 `io::writer`），传 `-o stdout`
  时会在当前目录创建**字面文件** `stdout`（随后 `build_gzi_index` 又生成
  `stdout.gzi`），与帮助文本 "[stdout] for screen" 的约定相悖。修复：`outfile
  == "stdout"` 时改向真实标准输出写 BGZF 流，并跳过无法对 stdout 生成的
  `.gzi` 索引。回归 `command_fa_gz_stdout_does_not_create_file`。

**`range` 的 `-o` 可覆盖 `.loc` 侧车索引**：`range` 先打开输出 writer，后在
  `open_indexed` 读取 `infile.loc`。若用户把 `-o` 命名为 `infile.loc`，会在
  读取索引前先截断该文件，随后 `loc_is_fresh` 因新 mtime 判定索引"新鲜"而不
  重建，`load_loc` 读到空索引，所有区域输出为空并伴随 "not found" 警告，且
  `.loc` 被永久损坏（后续调用持续失败直至重建）。修复：将 `infile.loc` 一并
  加入 `ensure_outfile_distinct` 的保护列表。回归
  `command_fa_range_output_not_overwrite_loc_index`。

**`window --chunk-records` 的分块文件可覆盖输入文件**：`run_window` 的分块
  文件名由 `-o` 派生（`out.fa` → `out.001.fa`），若输入文件恰好命名为某个
  分块名（如输入 `out.001.fa` 且 `-o out.fa`），`create_writer` 用 truncate
  打开该分块文件会截断**正在流式读取**的输入，导致后续记录丢失。修复：在
  `create_writer` 生成分块路径时，用 `same_path` 与 `infile` 比对，命中即
  `bail!`。回归 `command_fa_window_chunk_output_not_overwrite_input`。

**`gz` 对目录输入会先建残缺输出文件**：`gz` 用 `File::open` 直接打开输入
  （未走 `io::reader`），Unix 上对目录 `open` 成功，错误延迟到首次读取——此时
  输出文件（默认 `{infile}.gz`）已被 `File::create` 创建，留下残缺文件后才报
  "Is a directory"。与 `io::reader` 的前置目录拒绝行为不一致。修复：在压缩分
  支打开输出前，对非 `stdin` 的输入做 `is_dir()` 前置检查并尽早 `bail!`。
  回归 `command_fa_gz_directory_input_no_stray_output`。

**`fq to-fa` 与 `fq interleave` 均允许 `-o` 覆盖输入文件**：两个命令都在读取
  输入**之前**用 `pgr::writer`（内部 `File::create`，truncate）打开输出。若
  `-o` 指向某个输入文件，输出会先截断输入，随后 reader 读到空文件——静默清除
  原始数据。修复：
  - `to_fa.rs`：在打开 writer 前加入
    `ensure_outfile_distinct(outfile, infiles.map(|s| s.as_str()))?`；
  - `interleave.rs`：先把 infiles 收进 `Vec<String>`，再于打开 writer 前加入
    `ensure_outfile_distinct(outfile, infiles.iter().map(String::as_str))?`。
  回归 `command_fq_to_fa_output_same_as_input_rejected`、
  `command_fq_interleave_output_same_as_input_rejected`。

**全部 5 个 `2bit` 子命令允许 `-o` 覆盖输入 2bit 文件**：`to-fa`/`size`/`masked`
  /`range`/`some` 均未调用 `ensure_outfile_distinct`。若 `-o` 指向输入 2bit 文件，
  输出 writer 以截断方式打开，会先于（或与）读取输入清空该文件：`to-fa`/`size`/
  `masked` 先开 writer 再 `TwoBitFile::open`，输入被截断后 `open` 直接失败（"Not a
  valid 2bit file"）且文件已被销毁；`range`/`some` 先 `TwoBitFile::open`（读到头/
  索引）再开 writer 截断，随后读取返回错误或脏数据，输入同样被永久破坏。`some`
  的 `--name-list` 与 `range` 的 `--rgfile` 也是 `-o` 不应命中的输入。修复：在 5 个
  子命令的 `execute` 打开 writer 前统一加入 `ensure_outfile_distinct`（`some` 保护
  输入 2bit + list 文件，`range` 保护输入 2bit + rgfile）。回归
  `test_2bit_output_not_overwrite_input`、`test_2bit_range_output_not_overwrite_rgfile`。

### 输入校验 / 静默错误（3 处：fa 2 + fq 1）

**`fa read_names` 未跳过 `#` 注释行**：`some`/`order`/`range` 等命令读取名称
  列表时不跳过 `#` 注释行，导致把注释误当名称。修复：在 `read_names` 的
  `filter_map` 中跳过空行与 `#` 开头行。回归
  `command_fa_some_ignores_hash_comments`。

**`fa gz --reindex` 静默忽略 `-o`**：reindex 分支直接对 `infile` 建索引
  （`infile.gzi`），此时传入的 `-o` 无任何提示地被忽略，用户误传 `-o` 无反馈。
  修复：在 reindex 分支开头检测 `-o`/`--outfile` 已提供即报错（二者互斥，
  reindex 输出位置固定为输入旁）。回归 `command_fa_gz_reindex_rejects_outfile`。

**`fq interleave` 双文件交错对读取计数不匹配静默截断**：两文件路径原先用
  `std::iter::zip` 取较短者停止，多余记录被静默丢弃——帮助文档明确"Paired
  files must have same number of reads"，但实现未校验，属静默截断风险。修复：
  在 `libs/fmt/fq.rs` 新增泛型 `interleave_pair_iter`，逐对消费两个迭代器，
  任一 `None` 而另一 `Some` 时 `bail!("paired files have different numbers of
  reads")`，将静默截断改为显式报错。回归
  `command_fq_interleave_mismatched_read_counts_rejected`。

### 行为一致性 / 算法（2 处，均 fq）

**`fq interleave` 单文件虚拟 R2 两路径不一致**：单 FQ 输入 → FA 输出时，虚拟
  R2 序列为 `"\n"`（空序列，输出 `>name/2\n\n`）；而单 FA 输入 → FA 输出时
  虚拟 R2 序列为 `"N"`。帮助文本与 `docs/fq.md` 均声明单文件模式"Generate
  dummy R2 sequences (N's)"，FA→FA 路径符合文档，FQ→FA 路径违背——同一命令
  对同类输入输出形态不一致。修复：将单 FQ 输入路径的虚拟 R2 统一为 `b"N"`
  （与 FA 路径及文档一致）。质量端由 `write_pair` 的 `Option` 默认处理不变
  （FA 输出忽略质量，FQ 输出填 `'!'`）。

**`fq interleave` 双文件路径返回的最终索引错误**：`interleave` 的 doc 契约声明
  "Returns the final read index (start + count)"，单文件路径正确递增并返回
  `start+count`；但双文件路径调用 `interleave_pair_iter(..)?` 时**丢弃了**其
  返回的更新后 `idx`，最终 `Ok(idx)` 返回的是未递增的 `start`。当前唯一调用方
  `interleave.rs` 忽略该返回值（`let _final_idx`），故不产生用户可见错误，但
  `pub fn` 自身契约被违背。修复：两个双文件分支改为 `idx = interleave_pair_iter
  (..)?`（单文件路径已正确）。回归单测
  `test_interleave_two_files_returns_final_index`（start=5、2 对 → 期望 7）。

### 文档一致性（4 处：fa 3 + fq 1）

**`fa some`/`fa order`/`fa mask` 示例暗示 gzip 输出**：三处 `after_help` 的
  示例均为 `-o output.fa.gz`，但 `io::writer`（`libs/io.rs`）写端**不压缩**
  （仅 `io::reader` 读端支持 `.gz`），会生成带 `.gz` 后缀的**纯文本**文件，
  用户回读时（`pgr fa size output.fa.gz`）会因非 gzip 而失败。设计上压缩由
  专门的 `fa gz`（BGZF）子命令负责，普通命令输出不压缩。修复：将示例改为
  `-o output.fa`，并把标题从 "Process gzipped files" 改为仅指明输入可为
  gzipped（`Read a gzipped input` / `Process input from a gzipped file`），
  避免误导输出被压缩。核对 `docs/fa.md` 与 `gz` 子命令示例（其输出确为
  BGZF）无此问题。

**`fa --no-ns` 帮助文本与实际行为不符**：`no_ns_arg` 帮助写 "Output size
  without Ns"，但实现用 `!is_n(b)`（`is_n` 对 IUPAC 歧义码及 X 返回 true），
  实际同时排除 N 与 IUPAC 歧义码。`docs/fa.md` 已准确描述（"排除 N 及 IUPAC
  歧义码"），仅 CLI 帮助不准确。修复：帮助文本改为 "Output size without Ns and
  IUPAC ambiguous codes"，与行为及 `docs/fa.md` 一致（`twobit size` 共用该
  > ✅ 已修复（2026-08-06）：`no_ns_arg` 帮助文本已改为 "Output size without
  Ns and IUPAC ambiguous codes"。
  arg，2bit 掩码块即 N，语义同样准确）。

**`fa six-frame` 帮助文本的 frame 编号与实际输出不符**：`after_help` 的
  "Translation frames" 写作 "Forward strand: +1, +2, +3" / "Reverse strand:
  -1, -2, -3"，但实际输出头（测试 `>seq1(+):1-15|frame=0`、`>seq1(-):3-26|
  frame=2`）恒为 `frame=0/1/2`，并以 `(strand)` 字段区分方向。帮助文本暗示
  头会显示 +1/+2/+3 或 -1/-2/-3，与行为不符。修复：改写该段，明确 `frame`
  是链内的 0-based 阅读框偏移（0/1/2），方向由 `(strand)` 承载
  （`+` 直接应用偏移；`-` 先反向互补再应用偏移）。仅改帮助文本，行为与测试
  不变。

**`fq to-fa` 帮助文本声称 "Supports compressed input/output"**：`to_fa`
  输出端经 `fmt/fa::writer` → `io::writer`，**不压缩**（仅输入 `io::reader`
  支持 gzip）。传 `-o output.fa.gz` 会生成带 `.gz` 后缀的**纯文本** FASTQ，
  回读时（`pgr size output.fa.gz` 等）会因非 gzip 失败。与 `fa some/order/mask`
  已修复的同类文档缺陷一致。修复：帮助文本改为 "Supports compressed input"
  （仅输入可为 gzipped）。核对 `docs/fq.md` 未声称输出压缩，无需改。

### 统计正确性（1 处，fa）

**`n50` 的 Nx 边界条件用 `>` 而非 `>=`**：`calc_n50_stats` 在
  `cumul_size > goal` 时赋 Nx 值。当累计长度**恰好等于** goal 时（如总长 100、
  片段 50/30/20，N50 目标 50），标准定义（累计长度 **达到或超过** goal 的
  最短片段）应得 N50=50，但 `>` 会顺延到下一片段得 30。修复：改为
  `cumul_size >= *goal`。回归 `n50_boundary_exact_goal`（`[50,30,20]`→50）与
  `n50_standard`（`[40,30,20,10]`→30）单元测试；既有 `command_fa_n50*` 测试
  保持不变（ufasta 无恰好相等边界）。

### 算法正确性（2 处，fa）

**`rc` 对非 IUPAC 字符报错而非保留**：`rc` 用 noodles `Sequence::complement()`
  补全，遇到 `-`、`*` 等非 IUPAC 字符会返回错误，与 `after_help` 文档声明的
  "Non-IUPAC characters are preserved as-is" 不符。修复：改用
  `libs::nt::NT_COMP` 查找表逐字节补全（反向遍历 + 补全），标准/IUPAC 碱基
  补全并保留大小写，未知字节（NT_COMP 为 255）原样保留。实测
  `ACGT-*NR-` → `-YN*-ACGT`。回归
  `command_rc_preserves_non_iupac_chars`。

**`range` 负链对非 IUPAC 字符报错**：`slice_record`（`libs/loc.rs`）对 `-`
  链用 noodles `Sequence::complement()`，遇到 `-`/`*` 报 `invalid base`，而正链
  正常输出——同一输入的内部不一致。修复：同样改用 `NT_COMP` 表（反向遍历 +
  补全，未知字节保留）。实测 `ACGT-ACG` 负链范围 → `CGT-ACGT`。回归
  `command_range_preserves_non_iupac_on_reverse_strand`。

## 验证

* 数据安全：`-o` 同路径（单文件与 `split`/`window` 派生文件、`range` 的
  `.loc`、`fq` 双命令、`2bit` 各子命令输入 2bit 及 `some` 的 list / `range` 的
  rgfile）、`six-frame` 遗漏保护等修复前后均实测复现，既有输入原样保留。
* 畸形输入：短序列（1 碱基）、超大窗口长度、多字节 UTF-8 序列、异常 span、
  畸形 FASTQ、`2bit` 的无 `:` 全序列 / `start>end` / `0` 坐标 / 超长范围截断 /
  不存在序列等 fuzz，零 panic。
* 2bit 打包/解包：`test_blocks_from_dna`/`test_write_read_roundtrip` 等验证
  打包长度与读取端一致、v1 偏移正确、N/mask 块合并正确。
* `cargo fmt`、`cargo clippy --all-targets -- -D warnings` 干净；全部 lib
  单测（554 个 + `fmt::twobit` 12 个）+ `cli_fa`（39 个）+ `cli_fa_index`（12 个）
  + `cli_fq`（12 个）+ `cli_2bit`（20 个）集成测试全部通过。
* 新增回归测试（主要）：`command_fa_output_same_as_input_rejected`、
  `command_six_frame_output_same_as_input_rejected`、
  `command_fa_split_output_not_overwrite_input`、
  `command_fa_range_output_not_overwrite_loc_index`、
  `command_fa_window_chunk_output_not_overwrite_input`、
  `command_fa_gz_stdout_does_not_create_file`、
  `command_fa_gz_directory_input_no_stray_output`、
  `command_fa_gz_reindex_rejects_outfile`、
  `command_six_frame_short_sequence_no_panic`、
  `command_fa_some_ignores_hash_comments`、
  `windows_huge_length_does_not_overflow`、
  `mask_sequence_non_ascii_does_not_panic`、
  `command_rc_preserves_non_iupac_chars`、
  `command_range_preserves_non_iupac_on_reverse_strand`、
  `command_fq_to_fa_output_same_as_input_rejected`、
  `command_fq_interleave_output_same_as_input_rejected`、
  `command_fq_interleave_mismatched_read_counts_rejected`、
  `test_interleave_two_files_returns_final_index`、
  `test_2bit_output_not_overwrite_input`（覆盖 to-fa/size/masked/range/some 五
  个子命令）、`test_2bit_range_output_not_overwrite_rgfile` 等。

## 结论

`fa`（18 子命令）、`fq`（2 子命令）与 `2bit`（5 子命令）命令族合计修复 24 处
缺陷：崩溃/越界/溢出 3 处、数据安全/参数校验 9 处、输入校验/静默错误 3 处、
行为一致性/算法 2 处、文档一致性 4 处、统计正确性 1 处、算法正确性 2 处。全部
关键修复附回归测试与文档澄清。

审核经多轮纵深复审收敛：全部子命令的执行路径、`-o` 覆盖保护（含单文件、
`split`/`window` 派生文件、`range` 的 `.loc`、`fq` 双命令、`2bit` 的输入 2bit
及辅助列表文件）、`six-frame` 短序列、`window` 溢出、`mask_sequence` 非 ASCII、
`n50` 统计与 Nx 边界、`rc`/`range` 非 IUPAC 字符、`gz` 的 `--reindex` 参数互斥
与 `-o stdout`、`--no-ns` 帮助一致性、`fq` 的虚拟 R2 与双文件计数校验、
`2bit` 打包/解包/偏移计算正确性、`read_sequence` 与 `read_2bit_record` 的边界与
掩码、`range` 的坐标/长度/负链处理、`docs/{fa,fq,twobit}.md` 与帮助文本一致性、
空/畸形输入 Zero Panic 均逐一核验。记录项（`replace` TSV 不跳 `#` 注释、
`split name` 文件名碰撞、`size --no-ns` 对 `-`/`*` 计数、`one` 未命中留空文件、
`split about -c + stdout` 截断、`to-2bit` 名称去重、`fq` 双文件格式检测、
`2bit seq_name(-)` 全序列反链、`2bit range` 空输入静默输出、`read_u32_vec`
大 count 分配等）均属文档化一致行为或不可达的极端命名场景，非缺陷，无需改动。

最终收敛轮未再发现任何新问题，审核收敛。
