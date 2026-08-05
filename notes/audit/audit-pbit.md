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
* `get_contig`：`s >= e` 跳过（反向坐标在 CLI 层已被拒绝，见修复项）；
  `start`/`end` clamp 到 `[0, total_len]`，整数安全。
* 多参考路由：CIGAR 路径（`try_encode_segment_cigar`）与 LZ-diff 回退
  （`append_sample_with_paf`）均按当前参考的 `group_start/group_count`
  过滤 `ref_group_ids`，共享 contig 名的样本不会路由到错误参考。
* `append`/`append-ref` 原地更新：临时文件 + 原子重命名，失败不损坏原归档；
  `stage_work_path` 拒绝 `-o` 与输入归档相同（符文档"省略 -o 原地更新"）；
  `TempFileGuard` 未 `disarm` 时随 drop 删除临时文件，无残留。
* `args.rs` pbit 参数构建器（`--samples`/`--refs`/`--contigs`/`-s`/
  `-l` 等）帮助文本与命令一致；`-s` 短选项在 `create`（segment-size）与
  `stat`/`to-fa`（sample）分属不同子命令，无冲突。
* `to-fa`/`some` 样本名路径穿越防护：拒绝 `/`、`\`、`.`、`..`（另见修复项
  空名拒绝）。
* `paf_index.rs`：`coord_to_i32` 对 >i32::MAX 的 PAF 坐标拒绝并 skip；
  `query_start > query_end` / `target_start > target_end` 拒绝；`query_id`
  分配用 `or_insert(next_id)` 处理重复名；全部行解析失败时 bail（非 PAF）。
  内存分配量由 PAF 行数决定（用户输入），无攻击面。
* CIGAR 切片（`slice_cigar_by_query`）边界 D 处理经完整推演正确：边界 D 不入
  `sliced_ops` 但仍推进 `cur_t`，使 `target_start`/`target_end` 投影包含其
  目标跨度，ref_slice 与 ops 消费一致；`split_m_to_eqx` 的 `=`/`X`/`I`/`D`
  越界与 `rt/si` 完整消费校验齐全；反向链 `forward_to_rc_coords` 投影正确。
* `decompressor.rs` `delta_cache` 键含 `ref_start`/`ref_end`，正确处理 CIGAR
  delta 去重后不同 ref 切片；三层解码（LZ-diff/CIGAR）长度均对
  `delta_meta.raw_length` 校验；`get_contig` 智能切片的 `offset` 累计与
  `saturating_sub` 防下溢。
* 命令层 `create`/`append`/`append-ref`：`collect_samples_from_args` 拒绝
  `--name` 与 `-i`/`-p` 混用、`--paf` 与 `-i` 数量不匹配、重复样本名；
  参考索引/名称解析越界均报错。
* `range`/`some`/`stat`/`to-fa`：`-o` 覆盖输入防护（`ensure_outfile_distinct`
  / `same_path`）、路径穿越防护、`-s` 过滤样本不存在时报错、坐标 1-based
  含端点 → 0-based 半开转换正确、反向坐标在 CLI 层拒绝。

## 已知限制（有意保留）

* `Decompressor` 参考层 `SequenceReader::read_sequence` 在多参考归档中按
  `contig_groups[name]`（跨所有参考的段）拼接同名 contig，语义上是"跨参考
  拼接"，与 `stat --refs` 的展示歧义同源。当前无任何 `pbit` 命令走该路径
  （`range`/`some`/`to-fa` 均用 `get_contig`/`get_sample`），仅内部测试调用，
  故暂记为 API 限制而非缺陷。

## 修复的缺陷（共 31 处）

### 数据安全（`-o` 覆盖输入，7 处）

1. **`create` `-o` 截断输入**：`Compressor::create_multi` 用 `File::create`
   截断输出先于读参考/样本 FASTA、PAF、`--name` TSV。修复：`ensure_outfile_
   distinct(outfile, inputs)` 覆盖全部输入。
2. **`append` `-o` 截断样本输入**：`stage_work_path` 先把归档复制到 `-o`，
   覆盖样本/PAF/TSV 输入。修复：`-o` 存在时对样本输入做
   `ensure_outfile_distinct`。
3. **`append-ref` `-o` 截断参考**：`-o` 可能覆盖新参考 FASTA。修复：`-o`
   存在时对 `-r` 列表做 `ensure_outfile_distinct`。
4. **`range` `-o` 截断归档/rgfile**：writer 截断先于 `get_contig` 惰性读归档。
   修复：`ensure_outfile_distinct` 覆盖归档与 rgfile。
5. **`some` `-o` 截断归档/名单**：修复：`ensure_outfile_distinct` 覆盖归档与
   name_list。
6. **`stat` `-o` 截断归档**：修复：`ensure_outfile_distinct`。
7. **`to-fa` 输出等于输入归档**：`{outdir}/{sample_name}.fa` 可能等于输入归档
   路径（writer 截断先于 `get_sample` 惰性读归档）。修复：写前对每个输出路径
   `same_path` 检查，冲突即报错。

### 数据损坏 / 样本名冲突（1 处）

8. **样本名冲突静默数据损坏**：样本名来自 `-i` basename 或 TSV，名冲突（如
   `dup.1.fa`/`dup.2.fa` 都成 `dup`；或 `append` 追加已存在样本名）时
   `append_sample` 把两段样本 segment 合并进同一名字，提取返回错误序列且无
   报错。修复三处：`collect_samples_from_args` 单次命令内重复名 → 报错（提示
   用 `--name`）；`compressor` 新增 `has_sample` 查询归档内已存在名；`append`
   追加前逐样本 `has_sample` 检查，已存在 → 报错。修复顺序正确：`create` 在
   打开输出前拦截（不残留输出文件）；`append` 在 `stage_work_path` 临时副本
   上检查，bail 后 guard 随 drop 删除临时文件，原归档不受影响。

### 功能正确性（2 处）

9. **多参考 PAF 回退路由**：`append_sample_with_paf` 的 LZ-diff 回退未按当前
   参考过滤 `ref_group_ids`，共享 contig 名在多参考归档中落到错误参考的段。
   修复：按当前参考 `group_start/group_count` 过滤后路由。
10. **`range` 反向坐标静默空输出**：`chr1:100-50` 此前仅 `is_valid` 查
    `start != 0`，反向坐标被 `get_contig` 的 `s >= e` 静默跳过、无输出无警告。
    修复：CLI 层将 `start <= end` 纳入校验，显式报错。

### 展示歧义 / 超界 UX（2 处）

11. **`stat --refs` 多参考展示歧义**：多参考归档中同名 contig 原按
    `contig_name` 跨参考聚合，无法区分所属参考。修复：参考数 > 1 时每行前缀
    参考名（`ref_name<TAB>contig<TAB>count`），单参考输出不变。
12. **`range` 超界坐标静默空输出**：如 500 bp 的 contig 查询 `chr1:1000-2000`，
    `get_contig` clamp 后 `s >= e` 静默跳过、无输出无警告。修复：`get_contig`
    返回实际写入的 FASTA 条目数（`Result<usize>`），`range.rs` 对切片请求写
    0 条时 `log::warn!`。返回类型变更对所有调用方兼容。

### 溯源元数据（2 处）

13. **`create` 溯源元数据不完整**：`cmd_line` 不含样本输入信息。修复：在
    `cmd_line` 中追加每个样本的 `-i name:path`、可选 `-p paf`、可选 `@ref
    ref_spec`，完整记录输入来源。
14. **`append`/`append-ref` 溯源元数据不完整**：`append` 未记录追加样本的
    `-i name:path`/`-p paf`/`@ref ref_spec`；`append-ref` 未记录追加的 `-r`。
    修复：两命令在 `set_cmd_line` 前追加各自输入来源。

### 死代码 / 警告措辞（3 处）

15. **`append_sample` 重复空检查死代码**：`ref_group_ids.is_empty()` 检查重复，
    且注释 "empty reference (0 bp)" 误导（过滤后为空实为"该 contig 在当前
    参考无段"）。修复：删除冗余块。
16. **`append` 不可达 `None` 分支**：`match ref_spec` 含重复 `None => 0` 分支
    不可达（reachable 的 `None` 已在前处理 `num_refs > 1` 警告并返回 0），trigger
    `unreachable_patterns` 警告。修复：删除重复分支，行为不变。
17. **contig 无段警告措辞误导**：原 "empty reference (0 bp)" 令人误以为参考
    为空。修复：改 "has no segments in the current reference"，并删除与之重复
    的第二处 `ref_group_ids.is_empty()` 死代码块。

### 文档一致性（3 处）

18. **`docs/pbit.md` `--name` TSV 只写 3 列**：实现（及 Options 段）支持可选
    第 4 列 `ref_name`（多参考路由）。修复：补记第 4 列及缺省路由到参考 0 的
    说明。
19. **CLI 帮助文本 `--name` TSV 仍写 3 列**：`pbit_name_arg`、`create.rs`
    after_help 未提可选第 4 列。修复：两处补记 `[<TAB>ref_name]`（与修复项 18
    一致）。
20. **`docs/pbit.md` `to-fa` Notes 未提空名拒绝**：只写 `/`、`\`、`.`/`..`。
    修复：补记"样本名不能为空"。

### 内存 DoS（5 处）

21. **`packed_size` 无界分配**：`decode_delta`/`DeltaEntry::read_from` 以
    `vec![0u8; meta.packed_size]` 按归档声明分配，`packed_size`（u32，最大
    ≈4GB）无上限；恶意归档可在小文件中写膨胀 `packed_size`，`new` 的 delta
    扫描 seek 越过 EOF 仍成功，随后先分配 4GB 再读失败。修复：`format.rs` 新增
    `MAX_PACKED_SIZE`（256 MB），在 `Decompressor::new` 扫描与
    `DeltaEntry::read_from` 两处校验报错。
22. **delta 解压 gzip bomb**：`decode_delta` 的 LZ-diff 路径与 `unpack_cigar`
    以 `read_to_end` 分配解压后缓冲；几 KB 压缩流可解压出数 GB 原始 delta。
    修复：新增 `MAX_DELTA_UNCOMPRESSED`（256 MB），两处解压路径用
    `decoder.take(limit + 1)` 限长读取，超出即报错（`limit + 1` 保证恰好等于
    上限的合法 delta 不被误拒）。
23. **sample index gzip bomb**：`Collection::deserialize` 以 `read_to_end`
    一次性解压 sample index，在解析计数前分配整个原始缓冲；各计数虽有上限但
    解压发生在计数校验前。修复：`MAX_COLLECTION_UNCOMPRESSED`（256 MB）+
    `decoder.take(limit + 1)` 限长解压。`open_for_append` 经 `open` 复用。
24. **`min_match_len` 无界 → 段 padding DoS**：`Decompressor::new` 对
    `min_match_len` 仅校验 `<= segment_size`，而 `segment_size` 无上限；
    `decode_delta` 每解一段 `Segment::new(prepare)` 触发
    `reference.resize(ref_dna.len() + key_len)`（`key_len ≈ min_match_len - 3`）。
    恶意归档把两者都设 ~2GB 使校验相等通过，解码单段即触发 ~2GB 分配。修复：
    `Decompressor::new` 增加绝对上限（复用 `MAX_PACKED_SIZE` 256 MB）。
25. **CIGAR 段计数无界分配**：`unpack_cigar` 的 `op_count`/`xi_count` 是攻击者
    可控的 u32（可达 ~40 亿），直接用于 `Vec::with_capacity(op_count)`（~16GB）
    与 `vec![0u8; xi_count]`（~4GB），分配先于读取循环命中 EOF。对照：
    `collection.rs` 用 `with_capacity(x.min(1024))` 封顶、`read_string` 有
    `MAX_STRING_LEN=16MB` 上限、`lz_diff.rs` decode 均有界，唯独 `unpack_cigar`
    缺失。修复：分配前按实际 payload 大小校验 —— `op_count` 不超过
    `(raw_len - 8) / 4`，`xi_count` 不超过剩余 `raw_len - 8 - op_count*4`。

### 健壮性（1 处）

26. **`to-fa` 空样本名生成 `.fa` 点文件**：样本名来自归档 collection（不可信
    `.pbit` 输入），构造的恶意归档可嵌入空样本名，经 `to-fa`（无 `-s` 时遍历
    全部样本）产生散落的隐藏点文件。修复：在路径穿越守卫中补
    `sample.is_empty()` 检查，一并拒绝。

### 一致性 / 报告去重（2 处）

27. **`append` 多参考缺警告**：未指定 `ref_spec`（TSV 无第 4 列）时静默路由到
    参考 0，与 `create`/`resolve_ref_id` 的多参考警告不一致。修复：参考数 > 1
    且未指定参考时对每个样本 `log::warn!`。
28. **报告重复记录项**：第 8 轮 `packed_size` 条目在报告中重复出现（两处同标题
    与重复条目）。修复：删除重复块。

### 反向链压缩率（1 处）

29. **反向链 LZ-diff 段路由**：`encode_segment_lzdiff` 对反向互补 contig 按
    `seg_idx → ref_group_ids[seg_idx]` 正向一一对应，未按反向倒序匹配；解码用
    同一映射故正确，但样本段 i 实为参考段 N-1-i 的反向互补，正向路由让每段都
    错配 → 折叠后 delta 偏大。修复：`contig_is_rev_comp` 为真时改用
    `ref_group_ids[N-1-seg_idx]`（clamp 到首段）。因 `ref_group_id` 逐段写入并
    随归档存储、解码按同一映射还原，改动仅影响压缩率不影响正确性。

### 零 panic 校验一致性（2 处）

30. **`append`/`append-ref` 对畸形归档 panic**：`Decompressor::new` 对
    `segment_size`/`kmer_len` 仅校验正数，而 `open_for_append` 复用头部的值对
    样本/参考重新分段，`segment_sequence` 调 `chunks(0)`、`detect_rev_comp` 调
    `windows(0)`，参数为 0 时 panic。构造的归档（头部 `segment_size=0` 或
    `kmer_len=0`）可让两命令崩溃。修复：`Decompressor::new` 显式拒绝
    `segment_size == 0` 与 `kmer_len == 0`。
31. **`create` 的 `min_match_len` 校验与 decompressor 不一致**：`create` 只校验
    `<= segment_size`（上限 ≈2GB），可创建 `min_match_len` 达 2GB 的归档；
    decompressor 对 > `MAX_PACKED_SIZE`（256MB）绝对上限拒绝，导致 `create`
    生成 `stat`/`range`/`to-fa` 都拒绝的非法归档。修复：`create` 增加一致的上限
    校验（`min_match_len <= MAX_PACKED_SIZE`），并将 `MAX_PACKED_SIZE` 由
    `pub(crate)` 提升为 `pub`（二进制 crate 的 `cmd_pgr` 无法访问库 crate 的
    `pub(crate)` 项）。

## 验证

- 数据安全：`create`/`append`/`append-ref`/`range`/`some`/`stat`/`to-fa` 的
  `-o` 覆盖保护均先于 writer/`File::create` 打开，顺序正确。
- 新增回归测试（主要）：
  - 数据安全：`test_pbit_to_fa_output_not_overwrite_input`（修复项 7）。
  - 功能正确性：`test_pbit_range_reversed_coordinates_rejected`（10）；
    既有多参考路由 `test_pbit_multi_reference_routing`（9）。
  - 数据损坏：`test_pbit_create_duplicate_sample_name_rejected` 与
    `test_pbit_append_existing_sample_name_rejected`（8）。
  - 超界 UX：`test_pbit_range_out_of_bounds_warns`（12）。
  - 溯源元数据：`test_pbit_create_cmd_line_includes_samples`（13）、
    `test_pbit_append_cmd_line_includes_samples` 与
    `test_pbit_append_ref_cmd_line_includes_refs`（14）。
  - 内存 DoS：`test_decompressor_rejects_huge_packed_size`（21）、
    `test_unpack_cigar_rejects_gzip_bomb`（22）、
    `test_deserialize_rejects_gzip_bomb`（23）、
    `test_decompressor_rejects_huge_min_match_len`（24）、
    `test_unpack_cigar_rejects_huge_op_count` 与
    `test_unpack_cigar_rejects_huge_xi_count`（25）。
  - 零 panic 校验一致性：`test_decompressor_rejects_zero_segment_size` 与
    `test_decompressor_rejects_zero_kmer_len`（30）、
    `test_pbit_create_invalid_params_rejected` 追加绝对上限用例（31）。
  - 反向链压缩率：`test_append_rev_comp_sample_multi_segment`（29）。
- 死代码/警告措辞/文档/一致性类修复行为不变或仅文本，无需新增测试；既有测试
  验证无回归。
- `cargo test --lib pbit` 112 全绿；`cargo test --test cli_pbit pbit_` 46 全绿。
- `cargo build`、`cargo fmt`、`cargo clippy --all-targets -- -D warnings` 均
  clean。

## 结论

`pbit` 命令族审核完成（累计修复 31 处缺陷：数据安全 7 + 数据损坏 1 + 功能
正确性 2 + 展示歧义/UX 2 + 溯源元数据 2 + 死代码/警告措辞 3 + 文档一致性 3 +
内存 DoS 5 + 健壮性 1 + 一致性/报告去重 2 + 反向链压缩率 1 + 零 panic 校验
一致性 2），补回归测试与文档澄清，并经多轮纵深复审（首轮对 `compressor`/
`decompressor`/`lz_diff`/`format`/`collection`/`segment`/`paf_index` 与全部命令
层逐行深审；第 16 轮复核命令层与 PAF 索引；第 17 轮对核心库全文重读并核对文档
与 CLI 实现一致性）均未再发现新缺陷，审核收敛。

剩余记录项仅参考层 `read_sequence` 跨参考拼接（非命令可达，见"已知限制"）。