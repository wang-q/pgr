# pgr 2bit 命令族代码审核记录（2026-08-05）

对 `pgr 2bit` 命令族（全部 5 个子命令：`masked`/`size`/`range`/`some`/`to-fa`）
以及相关库文件（`libs/fmt/twobit`、`libs/ds/range`、`libs/nt`、`libs/loc`、
`libs/io`）、全部测试与文档（`docs/twobit.md`）进行审核。缺陷按类别分组记录；
关键修复均附回归测试（见文末），验证概况见文末"验证"一节。

审核范围：
- **信息**：`masked` / `size`
- **子集**：`range` / `some`
- **转换**：`to-fa`

审核重点：数据安全（`-o` 不得覆盖输入 2bit / 辅助列表文件）、Zero Panic（畸形
输入不 panic）、坐标/长度边界处理、2bit 打包/解包正确性、`docs/twobit.md` 与
帮助文本一致性。

## 与外部参考实现的语义一致性核对

2bit 家族对照 UCSC kent-tools 的 2bit 二进制格式规范（magic/version/index/record
布局、2-bit 打包、N-block/mask-block、双字节序）逐字节核对，语义一致。有意差异
（均已记录）：

* 读取端同时支持 version 0（index 偏移 4 字节）与 version 1（index 偏移 8 字节）
  两种布局；写入端产出 version 1。
* 打包时任何非 A/C/G/T 碱基（含 IUPAC 歧义码、U）按 N 处理并记录为 hard-mask
  块；软屏蔽小写 A/C/G/T 记录为 mask 块。`size --no-ns` 仅扣减 N-block，与 UCSC
  "排除 hard-masked N 位置"语义一致。
* 畸形输入：kent-tools 部分输入直接崩溃/告警，pgr 统一友好出错/跳过
  （Zero-Panic）。

## 排除的疑点（经核验无需修复）

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

* `read_u32_vec` 在畸形文件给出超大 `count` 时会尝试分配巨大缓冲（可能 OOM
  abort，而非 panic）。属全局既有模式与畸形输入鲁棒性范畴，低风险，未改。

## 修复的缺陷（共 1 处）

### 数据安全（`-o` 覆盖保护，1 处）

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

## 验证

* 数据安全：`-o` 同路径（各子命令输入 2bit、`some` 的 list、`range` 的 rgfile）
  修复前后均实测复现，输入字节原样保留，命令返回友好错误。
* 畸形输入：无 `:` 全序列、`start>end`、`0` 坐标、超长范围截断、不存在序列等
  fuzz，零 panic。
* `cargo build`、`cargo clippy --all-targets -- -D warnings`、`cargo fmt --check`
  全部干净；`cli_2bit` 20 个集成测试 + `libs::fmt::twobit` 12 个 lib 单测全通过。
* 新增回归测试：`test_2bit_output_not_overwrite_input`（覆盖 to-fa/size/masked/
  range/some 五个子命令）、`test_2bit_range_output_not_overwrite_rgfile`。
* 构建解阻后的完整验证（pbit 未提交改动曾阻塞整个 crate 编译，E0502，来自
  align/pgi 审核；等待其修复后重验）：
  * `cargo build`：通过。
  * `cargo test --test cli_2bit`：**20 个集成测试全部通过**，含新增回归
    `test_2bit_output_not_overwrite_input`、
    `test_2bit_range_output_not_overwrite_rgfile`，以及坐标/掩码/负链/正则兼容/
    真实 UCSC 文件等既有用例。
  * `cargo test --lib fmt::twobit`：**12 个库单测全部通过**（打包/解包/切片/掩码/
    版本/溢出等）。
  * `cargo clippy --all-targets -- -D warnings`：干净无警告。
  * `cargo fmt --check`：`twobit/*.rs`、`tests/cli_2bit.rs` 全部 `OK`（`cargo fmt`
    报告的 diff 全部来自他人 pbit 审核的未格式化文件，与本命令无关）。

## 结论

`2bit` 命令族审核完成（累计修复 1 处数据安全缺陷并补回归测试），并经多轮纵深
复审（全部 5 个子命令的执行路径、`-o` 覆盖保护（含 `some` 的 list 与 `range` 的
rgfile）、`libs/fmt/twobit` 的打包/解包/偏移计算正确性、`read_sequence` 与
`read_2bit_record` 的边界与掩码、`range` 的坐标/长度/负链处理、
`docs/twobit.md` 与帮助文本一致性）复核，未再发现新问题，审核收敛。

对全部记录项与已知限制（`seq_name(-)` 全序列反链、`range` 空输入静默输出、
`read_u32_vec` 大 count 分配）逐一核验，均属文档化一致行为或非崩溃的极端场景、
或超出本范畴的全局既有模式，非缺陷或另行拆分范畴，无需改动。`2bit` 命令族
累计修复 1 处数据安全缺陷，回归测试覆盖，构建/clippy/测试/fmt 全绿。