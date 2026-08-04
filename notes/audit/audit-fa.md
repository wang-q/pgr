# pgr fa 命令族代码审核记录（2026-08-05）

对 `pgr fa` 命令族（全部 18 个子命令）以及相关库文件（`libs/fmt/fa`、
`libs/translate`、`libs/io`、`libs/fasta/chunk`、`libs/fasta/stat`、
`libs/fasta/dedup`、`libs/fasta/filter`）和全部测试/文档进行审核。缺陷按
类别分组记录；关键修复均附回归测试（见文末），验证概况见文末"验证"一节。

审核范围：
- **info**：`size` / `count` / `masked` / `n50`
- **records**：`one` / `some` / `order` / `split` / `window`
- **transform**：`replace` / `rc` / `filter` / `dedup` / `mask` / `six-frame` / `to-2bit`
- **indexing**：`gz` / `range`

审核重点：数据安全（`-o` 不得覆盖输入）、Zero Panic（畸形输入不 panic）、
坐标/边界处理、文档一致性。

## 排除的疑点（经核验无需修复）

* 逐命令通读全部 18 个 `fa` 子命令的 `execute`：`unwrap()`/`unreachable!`
  全部为 clap `required` 参数或 `value_parser` 约束枚举，运行期不可达，
  无潜在 panic（符合"稳定性原则"）。
* `-o` 覆盖保护覆盖情况：全部单文件输出命令（`count`/`size`/`n50`/
  `masked`/`one`/`some`/`order`/`rc`/`replace`/`range`/`window`/`dedup`/
  `filter`/`to-2bit`/`gz`/`six-frame`）均调用 `ensure_outfile_distinct`；
  `split` 输出为目录，采用逐输出路径 `same_path` 反向检查（见下）。
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
* 文档一致性：`docs/fa.md` 中各子命令的参数、默认值、坐标约定与代码一致。

## 记录项（未改，低风险 / 待决策）

* `split` 的 `name` 模式：`sanitize_filename` 把 `/`、`\`、`(`、`)`、`:` 替换
  为 `_`，两个不同名称可能清洗到同一文件名（如 `a/b` 与 `a_b` 均 → `a_b`），
  会静默合并到同一文件。与文件名清洗方案固有行为一致，低风险，未改。
* `split` 的 `name` 模式：空/全特殊字符名称清洗后理论上可得空文件名，实际
  FASTA 名称非空，不可达，未改。
* `one` 在未找到序列时返回错误，但输出文件已先用 `File::create` 打开，会留下
  一个空文件（`-o`）。属一般的"出错留空文件"行为，非数据损坏，未改。
* `size --no-ns` 对 `-`/`*`/数字等非 IUPAC 字符计数为"有效碱基"：实现用
  `!is_n(b)`（仅排除 N 与 IUPAC），`-`/`*` 会被计入长度，与"仅计算有效碱基"
  的字面略有出入；但文档明确把排除范围限定为"N 及 IUPAC 歧义码"，行为与文档
  字面一致，且与 `count` 的 `len`（仅 A/C/G/T/N）语义不同属两命令各自定义。
  低风险文档歧义，未改。
* `read_replace_tsv` 不跳过 `#` 注释行，而 `read_names`（`some`/`order` 等）
  会跳过：若 TSV 带 `#` 表头/注释行，会被误当作 key（如 `#old\tnew` 生成
  key=`#old`）。`replace` 命令文档未承诺注释支持，属两读取函数一致性小瑕疵，
  且 `#` 前缀名称极罕见，低风险，未改。
* `masked --gap` 的 `is_n` 会把 IUPAC 歧义码（M/R/W/S/Y/K/V/H/D/B 及 X）一并
  计为 N：帮助文本写作 "Only identify regions of N/n (gaps)"，字面略窄于实现。
  与代码库"歧义码即 N"的统一口径一致（默认 `masked` 亦将歧义码视为 masked），
  属术语精度问题，非行为缺陷，未改。
* `split about -c N -o stdout`：`about` 模式无论输出到 stdout 还是目录，都会按
  `max_files`（默认 999）轮转并 `break 'outer`，即在 stdout 上仅输出前
  `max_files` 个分块容量内的记录。默认 `-c` 为 `usize::MAX`（不触发轮转）时
  不受影响；仅当用户显式 `-c SIZE -o stdout` 时才会提前截断。流式 stdout 与
  按大小分块本就矛盾，属极端组合，低风险，未改。

## 修复的缺陷（共 16 处）

### 崩溃 / 越界 / 溢出（Zero Panic，3 处）

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

### 数据安全 / 参数校验（`-o` 同输入保护，6 处）

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

**`range` 的 `-o` 可覆盖 `.loc` 侧车索引**：`range` 先在行 79 打开输出
  writer，后在 `open_indexed`（行 89）读取 `infile.loc`。若用户把 `-o` 命名
  为 `infile.loc`，会在读取索引前先截断该文件，随后 `loc_is_fresh` 因新 mtime
  判定索引"新鲜"而不重建，`load_loc` 读到空索引，所有区域输出为空并伴随
  "not found" 警告，且 `.loc` 被永久损坏（后续调用持续失败直至重建）。修复：
  将 `infile.loc` 一并加入 `ensure_outfile_distinct` 的保护列表。回归
  `command_fa_range_output_not_overwrite_loc_index`。

**`window --chunk-records` 的分块文件可覆盖输入文件**：`run_window` 的分块
  文件名由 `-o` 派生（`out.fa` → `out.001.fa`），若输入文件恰好命名为某个
  分块名（如输入 `out.001.fa` 且 `-o out.fa`），`create_writer` 用 truncate
  打开该分块文件会截断**正在流式读取**的输入，导致后续记录丢失。修复：在
  `create_writer` 生成分块路径时，用 `same_path` 与 `infile` 比对，命中即
  `bail!`。回归 `command_fa_window_chunk_output_not_overwrite_input`。

### 输入校验 / 静默错误（2 处）

**`read_names` 未跳过 `#` 注释行**：`some`/`order`/`range` 等命令读取名称
  列表时不跳过 `#` 注释行，导致把注释误当名称。修复：在 `read_names` 的
  `filter_map` 中跳过空行与 `#` 开头行。回归
  `command_fa_some_ignores_hash_comments`。

**`gz --reindex` 静默忽略 `-o`**：reindex 分支直接对 `infile` 建索引
  （`infile.gzi`），此时传入的 `-o` 无任何提示地被忽略，用户误传 `-o` 无反馈。
  修复：在 reindex 分支开头检测 `-o`/`--outfile` 已提供即报错（二者互斥，
  reindex 输出位置固定为输入旁）。回归 `command_fa_gz_reindex_rejects_outfile`。

### 文档一致性（2 处）

**`fa some`/`fa order`/`fa mask` 示例暗示 gzip 输出**：三处 `after_help` 的
  示例均为 `-o output.fa.gz`，但 `io::writer`（`libs/io.rs`）写端**不压缩**
  （仅 `io::reader` 读端支持 `.gz`），会生成带 `.gz` 后缀的**纯文本**文件，
  用户回读时（`pgr fa size output.fa.gz`）会因非 gzip 而失败。设计上压缩由
  专门的 `fa gz`（BGZF）子命令负责，普通命令输出不压缩。修复：将示例改为
  `-o output.fa`，并把标题从 "Process gzipped files" 改为仅指明输入可为
  gzipped（`Read a gzipped input` / `Process input from a gzipped file`），
  避免误导输出被压缩。核对 `docs/fa.md` 与 `gz` 子命令示例（其输出确为
  BGZF）无此问题。

**`--no-ns` 帮助文本与实际行为不符**：`no_ns_arg` 帮助写 "Output size without
  Ns"，但实现用 `!is_n(b)`（`is_n` 对 IUPAC 歧义码及 X 返回 true），实际同时
  排除 N 与 IUPAC 歧义码。`docs/fa.md` 已准确描述（"排除 N 及 IUPAC 歧义码"），
  仅 CLI 帮助不准确。修复：帮助文本改为 "Output size without Ns and IUPAC
  ambiguous codes"，与行为及 `docs/fa.md` 一致（`twobit size` 共用该 arg，
  2bit 掩码块即 N，语义同样准确）。

### 统计正确性（1 处）

**`n50` 的 Nx 边界条件用 `>` 而非 `>=`**：`calc_n50_stats` 在
  `cumul_size > goal` 时赋 Nx 值。当累计长度**恰好等于** goal 时（如总长 100、
  片段 50/30/20，N50 目标 50），标准定义（累计长度 **达到或超过** goal 的
  最短片段）应得 N50=50，但 `>` 会顺延到下一片段得 30。修复：改为
  `cumul_size >= *goal`。回归 `n50_boundary_exact_goal`（`[50,30,20]`→50）与
  `n50_standard`（`[40,30,20,10]`→30）单元测试；既有 `command_fa_n50*` 测试
  保持不变（ufasta 无恰好相等边界）。

### 算法正确性（2 处）

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

* 数据安全：`-o` 同路径（单文件与 `split` 目录）、`six-frame` 遗漏保护等
  修复前后均实测复现，既有输入原样保留。
* 畸形输入：短序列（1 碱基）、超大窗口长度、多字节 UTF-8 序列、异常 span
  等 fuzz，零 panic。
* `cargo fmt`、`cargo clippy --all-targets -- -D warnings` 干净；全部测试
  + doctest 通过（含 38 个 `cli_fa` 集成测试）。
* 新增回归测试（主要）：`command_fa_output_same_as_input_rejected`、
  `command_six_frame_output_same_as_input_rejected`、
  `command_fa_split_output_not_overwrite_input`、
  `command_fa_range_output_not_overwrite_loc_index`、
  `command_fa_window_chunk_output_not_overwrite_input`、
  `command_six_frame_short_sequence_no_panic`、
  `command_fa_some_ignores_hash_comments`、
  `windows_huge_length_does_not_overflow`、
  `mask_sequence_non_ascii_does_not_panic` 等。

## 结论

`fa` 命令族审核完成（累计修复 16 处缺陷、补回归测试与文档澄清），并经多轮
纵深复审（全部 18 个子命令的执行路径、`-o` 覆盖保护含 `range` 的 `.loc` 与
`window` 分块文件、`six-frame` 短序列、`window` 溢出、`mask_sequence` 非
ASCII、`split about` 边界、`n50` 统计与 Nx 边界、`rc`/`range` 非 IUPAC 字符、
`gz --reindex` 参数互斥、`gz -o stdout` 写标准输出、`--no-ns` 帮助一致性、
坐标/文档一致性、空/畸形输入 Zero Panic）复核。

**最终收敛轮**：对全部记录项（`replace` TSV 不跳 `#` 注释、`split name`
文件名碰撞合并、`size --no-ns` 对 `-`/`*` 计数、`one` 未命中留空文件、
`split about -c + stdout` 截断、`n50` `-N` 默认与 Append、`count` 跨文件
total 聚合、`to-2bit` 名称去重顺序等）逐一重新核验，均属文档化一致行为或
不可达的极端命名场景，非缺陷，无需改动。此轮未再发现任何新问题，审核收敛。

**追加复审轮（2026-08-05）**：重新通读全部 18 个子命令的 `execute` 与核心库
（`fmt/fa`、`nt`、`io`、`loc`、`translate`、`fasta/{stat,filter,chunk,dedup}`），
并复核 `docs/fa.md` 一致性。重点核验：
- `mask_sequence` 的 `upper < lower` 兜底：确认 `IntSpan::try_from` /
  `runlist_to_ranges` 对反序范围（如 `5-3`）直接判非法并报错，`mask` 经
  `read_runlist` 只会拿到 `lower <= upper` 的合法 span，故 `offset+length`
  溢出路径在 `fa` 命令上下文不可达（此前记录的安全兜底逻辑描述有误，但
  结论——实际不可达——成立）。
- `window --chunk-records` 分块路径：`create_writer` 的 `same_path` 检查在
  `truncate` 打开前执行，首个分块与流式中途分块对输入文件的碰撞均在截断前
  拦截。
- `rc`/`range` 负链的 `NT_COMP` 非 IUPAC 保留、`n50` 的 `>=` 边界、`count`
  的 total 聚合、`to-2bit` 的 `u32` 长度上限等均正确。
- 全部 `fa` 文件 `cargo fmt --check` 干净；`cargo clippy --all-targets -- -D
  warnings` 干净；551 个 lib 单测 + 38 个 `cli_fa` + 12 个 `cli_fa_index`
  集成测试全部通过。
- 注：`cargo fmt --check` 在 `src/cmd_pgr/fas/to_xlsx.rs` 有一处既有格式差异，
  属 `fas` 命令族（与 `fa` 无关，处于他人进行中的改动），不在本次 `fa` 审核
  范围内。

此追加轮未发现任何新的 `fa` 缺陷，审核保持收敛。