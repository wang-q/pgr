# pbit 命令族代码审核记录（2026-08-05）

对 `pgr pbit` 命令族（`create` / `append` / `append-ref` / `stat` / `range` /
`some` / `to-fa`）及核心库（`libs/pbit`：`compressor` / `decompressor` /
`format` / `collection` / `segment` / `lz_diff` / `cigar_delta` / `paf_index`）
与全部测试、文档进行多轮审核。缺陷按类别分组记录；关键修复均附回归测试，
验证概况见文末"验证"一节。

## 与格式规范的一致性核对

pbit 格式（`notes/design/pbit.md` §文件格式规范 v1004）：

- Header 36 字节、Footer 24 字节、固定大小小端整数、字符串 u32 长度前缀，
  `format.rs` 逐字段读写与规范一致（含 36/24/10 字节长度断言测试）。
- `Decompressor::new` 就地校验 `ref_group_count`（header/index/delta 三处）、
  `sample_count`（header/collection）、`sample_index_offset <= footer_start`，
  损坏文件不会越界。
- 参考层复用标准 2bit 记录（`read_2bit_record`），保留 N/mask blocks，
  与 twobit 共享代码。

## 排除的疑点（经核验无需修复）

* `format.rs` 读取：magic/version 校验、`read_string` 的 `MAX_STRING_LEN`
  （16 MB）防护、`PbitFooter::read_at_end` 的 `file_size < 24` 检查，均防
  损坏文件。
* `decompressor.rs` `decode_delta`：`ref_group_id`/`delta_id` 越界、CIGAR
  `ref_start >= ref_end`、`ref_end > ref_dna.len()` 均显式报错，无 panic。
* `lz_diff.rs` `decode`：N-run 长度超 ref_len、`ref_pos` 计算溢出、match
  区间越界、畸形分隔符（缺 `N_CODE`/`.`/`,`）均返回 `Err`，抗崩溃完善。
* `cigar_delta.rs` `apply_cigar`：`=`/`X`/`D` 超参考长、`X`/`I` 超碱基流、
  `M` 意外出现、`xi`/`rt` 消费不一致均显式报错，无 panic。
* `get_contig`：`s >= e` 跳过（反向坐标在 CLI 层已被拒绝，见记录项 9）；
  `start`/`end` clamp 到 `[0, total_len]`，整数安全。
* 多参考路由：CIGAR 路径（`try_encode_segment_cigar`）与 LZ-diff 回退
  （`append_sample_with_paf`）均按当前参考的 `group_start/group_count`
  过滤 `ref_group_ids`，共享 contig 名的样本不会路由到错误参考。
* `append`/`append-ref` 原地更新：临时文件 + 原子重命名，失败不损坏原归档；
  `TempFileGuard` 未 `disarm` 时随 drop 删除临时文件，无残留。
* `args.rs` pbit 参数构建器（`--samples`/`--refs`/`--contigs`/`-s`/
  `-l` 等）帮助文本与命令一致。
* `to-fa`/`some` 样本名路径穿越防护：拒绝 `/`、`\`、`.`、`..`。

## 记录项（未改，低风险 / 待决策）

* `create` 记录的 `cmd_line` 不含样本 `-i`/`--name` 信息（仅含 `-r`/`-o`/
  `-s`/`-k`/`-l`），属于归档内的溯源元数据不完整，不影响正确性与数据。
* `to-fa` 对空样本名会生成 `{outdir}/.fa`；实际样本名来自 basename 或
  TSV，非空，运行期不可达。

## 已修复缺陷（按类别分组）

### 数据安全（`-o` 覆盖输入）

1. `create`：`-o` 可能截断参考 FASTA、样本 FASTA、PAF 或 `--name` TSV
   输入（`Compressor::create_multi` 用 `File::create` 截断输出先于读输入）。
   修复：`ensure_outfile_distinct(outfile, inputs)` 覆盖全部输入。
2. `append`：指定 `-o` 时 `stage_work_path` 先把归档复制到 `-o`，覆盖样本/
   PAF/TSV 输入。修复：`-o` 存在时对样本输入做 `ensure_outfile_distinct`。
3. `append-ref`：`-o` 可能覆盖新参考 FASTA。修复：`-o` 存在时对 `-r` 列表
   做 `ensure_outfile_distinct`。
4. `range`：`-o` 可能覆盖输入归档或 `-r/--rgfile` 列文件（writer 截断先于
   `get_contig` 惰性读归档）。修复：`ensure_outfile_distinct` 覆盖归档与
   rgfile。
5. `some`：`-o` 可能覆盖输入归档或名单文件。修复：`ensure_outfile_distinct`
   覆盖归档与 name_list。
6. `stat`：`-o` 可能覆盖输入归档。修复：`ensure_outfile_distinct`。
7. `to-fa`：生成的 `{outdir}/{sample_name}.fa` 可能等于输入归档路径（writer
   截断先于 `get_sample` 惰性读归档）。修复：写前对每个输出路径
   `same_path` 检查，冲突即报错。

### 功能正确性

8. 多参考 PAF 回退路由：`append_sample_with_paf` 的 LZ-diff 回退未按当前
   参考过滤 `ref_group_ids`，共享 contig 名在多参考归档中会落到错误参考的
   段。修复：按当前参考 `group_start/group_count` 过滤后路由。
9. `range` 反向坐标（如 `chr1:100-50`）：此前 `is_valid` 仅查 `start != 0`，
   反向坐标被 `get_contig` 的 `s >= e` 静默跳过、无输出无警告。修复：在
   CLI 层将 `start <= end` 纳入校验，显式报错。

### 文档一致性

10. `docs/pbit.md` 「输入输出格式 · `--name` TSV」段只写了 3 列，而实现
    （及 Options 段）支持可选第 4 列 `ref_name`（多参考路由）。修复：补记
    第 4 列及缺省路由到参考 0 的说明。

## 验证

- `create`/`append`/`append-ref`/`range`/`some`/`stat` 的 `-o` 覆盖保护
  均先于 writer/`File::create` 打开，顺序正确。
- 新增回归测试：
  - `test_pbit_to_fa_output_not_overwrite_input`（记录项 7）
  - `test_pbit_range_reversed_coordinates_rejected`（记录项 9）
  - 既有多参考路由测试 `test_pbit_multi_reference_routing` 覆盖记录项 8。
- `cargo build`、`cargo fmt`、`cargo clippy --all-targets -- -D warnings`
  均 clean；全量 `cargo test` 通过（含 `cli_pbit` 39 项、pbit 库 102 项）。

## 结论

共修复 10 处缺陷（数据安全 7 + 功能 2 + 文档 1）。多轮复核未再发现新问题，
pbit 命令族审核收敛。