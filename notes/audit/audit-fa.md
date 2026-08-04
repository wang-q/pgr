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

## 修复的缺陷（共 6 类）

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

### 数据安全 / 参数校验（`-o` 同输入保护，3 处）

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

### 输入校验 / 静默错误（1 处）

**`read_names` 未跳过 `#` 注释行**：`some`/`order`/`range` 等命令读取名称
  列表时不跳过 `#` 注释行，导致把注释误当名称。修复：在 `read_names` 的
  `filter_map` 中跳过空行与 `#` 开头行。回归
  `command_fa_some_ignores_hash_comments`。

## 验证

* 数据安全：`-o` 同路径（单文件与 `split` 目录）、`six-frame` 遗漏保护等
  修复前后均实测复现，既有输入原样保留。
* 畸形输入：短序列（1 碱基）、超大窗口长度、多字节 UTF-8 序列、异常 span
  等 fuzz，零 panic。
* `cargo fmt`、`cargo clippy --all-targets -- -D warnings` 干净；全部测试
  + doctest 通过（含 34 个 `cli_fa` 集成测试）。
* 新增回归测试（主要）：`command_fa_output_same_as_input_rejected`、
  `command_six_frame_output_same_as_input_rejected`、
  `command_fa_split_output_not_overwrite_input`、
  `command_six_frame_short_sequence_no_panic`、
  `command_fa_some_ignores_hash_comments`、
  `windows_huge_length_does_not_overflow`、
  `mask_sequence_non_ascii_does_not_panic` 等。

## 结论

`fa` 命令族审核完成（累计修复 6 类缺陷、补回归测试与文档澄清），并经多轮
纵深复审（全部 18 个子命令的执行路径、`-o` 覆盖保护、`six-frame` 短序列、
`window` 溢出、`mask_sequence` 非 ASCII、`split about` 边界、`n50` 统计、
坐标/文档一致性）复核，未再发现新问题，审核收敛。