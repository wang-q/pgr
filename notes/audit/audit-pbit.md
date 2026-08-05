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

原暂记的两处低风险待决策项均已在此前轮次或第 13 轮修复，见记录项 22（`to-fa`
空样本名）与 28（反向链 LZ-diff 段路由）。

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

### 第 2 轮追加

11. `compressor.rs` `append_sample` 中存在重复的 `ref_group_ids.is_empty()`
    检查（第 552 行已处理，其后的第二处为死代码，且注释"empty reference (0 bp)"
    易误导——过滤后为空实为"该 contig 在当前参考无段"）。修复：删除冗余块。
12. CLI 帮助文本（`args.rs` `pbit_name_arg`、`create.rs` after_help）中
    `--name` TSV 格式仍只写 3 列，未提可选第 4 列 `ref_name`（与已修复的
    文档项 10 不一致）。修复：两处帮助文本补记 `[<TAB>ref_name]`。

## 验证

- `create`/`append`/`append-ref`/`range`/`some`/`stat` 的 `-o` 覆盖保护
  均先于 writer/`File::create` 打开，顺序正确。
- 新增回归测试：
  - `test_pbit_to_fa_output_not_overwrite_input`（记录项 7）
  - `test_pbit_range_reversed_coordinates_rejected`（记录项 9）
  - 既有多参考路由测试 `test_pbit_multi_reference_routing` 覆盖记录项 8。
- 第 2 轮改动（项 11、12）为删死代码与帮助文本，无需新增测试；既有
  `cargo test pbit`（库 102 + `cli_pbit` 39）全绿验证无回归。
- `cargo build`、`cargo fmt`、`cargo clippy --all-targets -- -D warnings`
  均 clean。

## 第 3 轮深审（无新缺陷）

逐行复核 `compressor`（`append_sample`/`append_reference`/`open_for_append`/
`finish`、`slice_cigar_by_query`/`split_m_to_eqx`/`try_encode_segment_cigar`
的 ± 链投影与 ref 投影）、`decompressor`（`decode_delta` 的 LZ-diff/CIGAR
双路径、`get_contig`/`get_sample` 的智能切片与长度校验、`SequenceReader`）、
`lz_diff`（V2 编解码、N-run、match-to-end、`!` 回引、溢出防护）、`format`、
`collection`、`segment`、`paf_index`，以及全部命令层与 `args.rs` 参数，未发现
新的正确性缺陷。

- 确认 `append_reference` 在截断点追加 2bit 记录、`finish` 重算
  `ref_index_offset`、`deltas.resize` 保留旧段，原地追加参考正确。
- 确认 ± 链 CIGAR 编码：`try_encode_segment_cigar` 对 `-` 链用
  `forward_to_rc_coords` 切片、`sample_slice = rev_comp(seg)`、
  `is_rev_comp = strand == '-'`，解码端 `apply_cigar` 后按 `meta.is_rev_comp`
  还原，正/反向链样本 roundtrip 均正确。
- 确认 `lz_diff` N-run 不推进 `pred_pos` 与编码端一致，roundtrip 含 N-run 通过。
- 确认 `decompressor` 各索引越界、长度不匹配、损坏 `raw_length` 均显式报错，
  无 panic（`Zero Panic` 满足）。

## 第 4 轮（样本名冲突 → 数据损坏）

13. 样本名冲突导致静默数据损坏：样本名来自 `-i` 文件 basename 或 TSV，若
    名冲突（如 `dup.1.fa`/`dup.2.fa` 都成 `dup`；或 `append` 追加已存在
    的样本名），`append_sample` 会把两段样本的 segment 合并进同一名字，
    提取时返回错误序列且无任何报错。修复三处：
    - `mod.rs` `collect_samples_from_args`：单次命令内重复样本名 → 报错
      （提示用 `--name` 显式命名）。
    - `compressor.rs` 新增 `has_sample`：查询归档内是否已存在指定样本名。
    - `append.rs`：追加前对每个样本名 `has_sample` 检查，已存在 → 报错。
  修复顺序正确：`create` 在 `Compressor::create_multi` 打开输出前拦截重复名
  （不残留输出文件）；`append` 在 `stage_work_path` 复制出的临时副本上检查，
  bail 后 guard 随 drop 删除临时文件，原归档不受影响。

新增回归测试：
- `test_pbit_create_duplicate_sample_name_rejected`（记录项 13 第一点）
- `test_pbit_append_existing_sample_name_rejected`（记录项 13 第二、三点）

## 第 5 轮（展示歧义与超界 UX）

14. `stat --refs` 多参考展示歧义（此前暂记为待决策）：多参考归档中同名
    contig（如两个参考都有 `chr1`）原按 `contig_name` 跨参考聚合，无法区分
    所属参考。修复：`stat.rs` 检测到参考数 > 1 时，每行前缀参考名
    （`ref_name<TAB>contig<TAB>count`），单参考输出 `contig<TAB>count` 不变。
    既有 `test_pbit_stat_refs`（单参考 `chr1\t2`）验证无回归。
15. `range` 超界坐标静默空输出（此前暂记为待决策）：如 500 bp 的 contig 查询
    `chr1:1000-2000`，`get_contig` clamp 后 `s >= e` 静默跳过、无输出无警告。
    修复：`get_contig` 返回实际写入的 FASTA 条目数（`Result<usize>`），
    `range.rs` 对切片请求写入 0 条时报 `log::warn!`（"nothing extracted"）。
    返回类型变更对所有调用方兼容（`some.rs`/基准/测试均用 `?` 或丢弃）。
    新增回归测试 `test_pbit_range_out_of_bounds_warns`。

## 第 6 轮（溯源元数据与警告措辞）

16. `create` 记录的 `cmd_line` 不含样本输入信息（此前暂记为待决策）：归档
    内溯源元数据不完整，无法追溯每个样本的 FASTA/PAF/参考来源。修复：
    `create.rs` 在 `cmd_line` 中追加每个样本的
    `-i name:path`、可选 `-p paf`、可选 `@ref ref_spec`，完整记录输入来源。
    新增回归测试 `test_pbit_create_cmd_line_includes_samples`（验证 `-i s1:`、
    `-p` PAF 路径、`@ref ref_2000`、`-i s2:` 均写入）。
17. `compressor.rs` contig 无段警告措辞误导：原 "empty reference (0 bp)"
    令人误以为参考为空，实为"该 contig 在当前参考无段（可能存在于其他参考）"。
    修复：改warning文案为 "has no segments in the current reference"，并删除
    与之重复的第二处 `ref_group_ids.is_empty()` 死代码块。

## 第 7 轮（append / append-ref 溯源元数据）

18. `append` / `append-ref` 的 `cmd_line` 溯源不完整（与第 6 轮 `create` 修复
    不一致）：`append` 只记录 `pgr pbit append infile [-o out]`，未记录追加的
    样本 `-i name:path`/`-p paf`/`@ref ref_spec`；`append-ref` 未记录追加的
    `-r ref`。`append-ref` 还可能在 `cmd_line` 中丢失 `-r` 信息（它本就只记
    归档名与 `-o`）。修复：两命令在 `set_cmd_line` 前追加各自的输入来源。
    新增回归测试 `test_pbit_append_cmd_line_includes_samples`（验证 `-i s2:<path>`）
    与 `test_pbit_append_ref_cmd_line_includes_refs`（验证 `-r <ref>`）。

## 第 8 轮（delta packed_size 无界 → 内存 DoS）

19. `Decompressor::decode_delta` 与 `DeltaEntry::read_from` 均以
    `vec![0u8; meta.packed_size]` 按归档声明的 `packed_size` 分配内存，而
    `packed_size`（u32，最大 ≈4GB）无任何上限。恶意归档可在小文件中写入
    膨胀的 `packed_size`，使 `new` 的 delta 扫描 seek 越过 EOF 仍成功，随后
    `decode_delta` 先分配 4GB 缓冲区再读失败 → 内存耗尽 DoS。修复：
    `format.rs` 新增 `MAX_PACKED_SIZE`（256 MB，远超任何真实单段 delta），
    在 `Decompressor::new` 的 delta 扫描与 `DeltaEntry::read_from` 两处校验并
    报错。新增回归测试 `test_decompressor_rejects_huge_packed_size`。

## 第 9 轮（delta 解压无界 → gzip bomb 内存 DoS）

20. `decode_delta` 的 LZ-diff 路径以 `decoder.read_to_end(&mut delta)` 按
    解压后大小分配内存，`unpack_cigar` 同样 `read_to_end`。`packed_size`
    虽已设上限（记录项 19），但压缩包可解压出远超 `packed_size` 的膨胀
    数据（gzip bomb）：恶意归档用几 KB 的压缩流解压出数 GB 的原始 delta，
    `decode_delta` 先分配该缓冲再被 `Segment::get` / `apply_cigar` 消费 →
    内存耗尽 DoS。修复：新增 `MAX_DELTA_UNCOMPRESSED`（256 MB），在
    LZ-diff 与 CIGAR 两处解压路径用 `decoder.take(limit + 1)` 限长读取，
    超出即报错（`limit + 1` 保证恰好等于上限的合法 delta 不被误拒）。
    新增回归测试 `test_unpack_cigar_rejects_gzip_bomb`（压缩 ~256MB 零，
    断言小粒度压缩流被拒绝）。

## 第 10 轮（sample index 解压无界 → gzip bomb 内存 DoS）

21. `Collection::deserialize` 以 `decoder.read_to_end(&mut raw)` 一次性解压
    sample index（flate2 压缩），在解析任何计数之前就分配整个原始缓冲。
    各计数虽有上限（`MAX_SAMPLE_COUNT` 等），但解压发生在计数校验之前，
    恶意归档仍可用几 KB 的压缩流解压出数 GB → `Decompressor::new` 内存耗尽
    DoS。修复：新增 `MAX_COLLECTION_UNCOMPRESSED`（256 MB，覆盖 ~20 个
    全基因组样本；一致于 delta 256 MB 上限），用 `decoder.take(limit + 1)`
    限长解压，超出即拒。`open_for_append` 经 `Decompressor::open` 复用该
    防护。新增回归测试 `test_deserialize_rejects_gzip_bomb`。

## 第 11 轮（to-fa 空样本名）

22. `to-fa` 对空样本名会生成 `{outdir}/.fa`（此前暂记为低风险不可达）。样本名
    来自归档 collection（不可信的 `.pbit` 输入），构造的恶意归档可嵌入空样本名，
    经 `to-fa`（无 `-s` 时遍历全部样本）产生散落的隐藏点文件。修复：在既有
    路径穿越守卫中补 `sample.is_empty()` 检查，一并拒绝。属健壮性改进，无新增
    测试（既有 CLI 测试验证无回归）。

## 第 12 轮（min_match_len 无界 → 段 padding 内存 DoS）

23. `Decompressor::new` 对 `min_match_len` 仅校验 `<= segment_size`，而
    `segment_size` 本身无上限。`decode_delta` 每解一段都新建
    `Segment::new(header.min_match_len)` 并 `prepare`，后者
    `reference.resize(ref_dna.len() + key_len)`（`key_len ≈ min_match_len - 3`）。
    恶意归档把 `segment_size` 与 `min_match_len` 都设为 ~2GB（使第一处校验
    相等通过），解码单段即触发 ~2GB 分配 → 内存 DoS。修复：在
    `Decompressor::new` 增加绝对上限（复用 `MAX_PACKED_SIZE` 256 MB），
    覆盖 `min_match_len` 大于该值的归档。新增回归测试
    `test_decompressor_rejects_huge_min_match_len`（`segment_size` 与
    `min_match_len` 均 = `MAX_PACKED_SIZE + 1`，断言被拒）。

## 第 12 轮补充（一致性 / 文档）

24. `append` 在未指定 `ref_spec`（TSV 无第 4 列）时静默路由到参考 0，与
    `create`/`resolve_ref_id` 的多参考警告不一致。修复：`append.rs` 在
    参考数 > 1 且未指定参考时对每个样本 `log::warn!`。
25. `docs/pbit.md` `to-fa` 的 Notes 只写了样本名不能含 `/`、`\`、`.`/`..`，
    未提第 11 轮新增的空名拒绝。修复：补记"样本名不能为空"。
26. 第 8 轮记录项 19 在报告中重复出现（两处 `## 第 8 轮` 标题与重复的
    `packed_size` 条目）。修复：删除重复块。

## 第 13 轮（死代码 / 反向链压缩率）

27. `append.rs` 的 `match ref_spec` 含重复的 `None => 0` 分支：第二个 `None`
    为不可达模式（reachable 的 `None` 分支已在前处理 `num_refs > 1` 警告并返回
    0），是死代码，且会触发 `unreachable_patterns` 警告（`clippy -D warnings`
    下为错误）。修复：删除重复分支。行为不变（`None` 仍路由到参考 0 并告警）。
28. 反向链 LZ-diff 段路由（此前暂记为低风险待决策）：`encode_segment_lzdiff`
    对反向互补 contig 按 `seg_idx → ref_group_ids[seg_idx]`（正向顺序）一一对应，
    未按反转向量倒序匹配。解码用同一映射故正确，但样本段 i 实为参考段 N-1-i 的
    反向互补，正向路由会让每个段都匹配到错误的参考段 → 折叠后 delta 偏大，压缩率
    低于最优。修复：`contig_is_rev_comp` 为真时改为 `ref_group_ids[N-1-seg_idx]`
    （clamp 到首段）。因 `ref_group_id` 逐段写入并随归档存储、解码按同一映射还原，
    改动仅影响压缩率不影响正确性。新增回归测试
    `test_append_rev_comp_sample_multi_segment`（5000 bp 参考 2 段 + 反向互补样本，
    断言 roundtrip 精确还原）。

### 验证（第 13 轮）

- `cargo test --lib pbit` 108 全绿（含新增 `test_append_rev_comp_sample_multi_segment`）。
- `cargo test --test cli_pbit pbit_` 46 全绿。
- `cargo clippy --all-targets -- -D warnings` clean（记录项 27 消除后无警告）。

## 第 14 轮（append/append-ref 零 panic 与 create 校验一致性）

29. `append` / `append-ref` 处理畸形归档可能 panic：`Decompressor::new` 对
    `segment_size` / `kmer_len` 仅校验正数，而 `open_for_append` 复用头部的
    `segment_size`/`kmer_len` 对样本/参考 FASTA 重新分段，
    `segment_sequence` 会调用 `chunks(0)`、`detect_rev_comp` 调用
    `windows(0)`，二者在参数为 0 时 panic。构造的归档（头部 `segment_size=0`
    或 `kmer_len=0`）可让 `append` / `append-ref` 崩溃（违反 Zero Panic）。
    修复：`Decompressor::new` 显式拒绝 `segment_size == 0` 与 `kmer_len == 0`。
30. `create` 的 `min_match_len` 校验与 decompressor 不一致：`create` 只校验
    `min_match_len <= segment_size`，而 `segment_size` 上限为 `i32::MAX`
    （≈2GB），可创建 `min_match_len` 达 2GB 的归档；decompressor（记录项 23）
    对 `min_match_len > MAX_PACKED_SIZE`（256MB）绝对上限拒绝，导致 `create`
    生成 `stat`/`range`/`to-fa` 都会拒绝的非法归档。修复：`create` 增加与
    decompressor 一致的绝对上限校验（`min_match_len <= MAX_PACKED_SIZE`）。
    为此将 `format.rs` 的 `MAX_PACKED_SIZE` 由 `pub(crate)` 提升为 `pub`
    （二进制 crate 的 `cmd_pgr` 无法访问库 crate 的 `pub(crate)` 项）。

新增回归测试：
- `test_decompressor_rejects_zero_segment_size` / 
  `test_decompressor_rejects_zero_kmer_len`（记录项 29，库内测试）
- `test_pbit_create_invalid_params_rejected` 追加绝对上限用例
  （`-s/-l 300000000`，断言 "must not exceed the per-segment bound"）（记录项 30）

### 验证（第 14 轮）

- `cargo test --lib pbit` 108 全绿。
- `cargo test --test cli_pbit pbit_` 46 全绿。
- `cargo build`、`cargo fmt`、`cargo clippy --all-targets -- -D warnings` clean。

## 第 15 轮（CIGAR 段计数无界 → 内存 DoS）

31. `cigar_delta.rs` `unpack_cigar` 的 `op_count` / `xi_count` 无界分配：
    解压后的 payload 虽受 `MAX_DELTA_UNCOMPRESSED`（256MB）约束，但这两个
    计数字段是攻击者可控的 u32（可达 ~40 亿）。`op_count` 直接用于
    `Vec::with_capacity(op_count)`（每个 op 4 字节，~40 亿 → ~16GB 分配），
    `xi_count` 直接用于 `vec![0u8; xi_count]`（~4GB 分配），且分配先于读取
    循环命中 EOF，恶意归档可凭极小的 CIGAR delta 触发 OOM/abort（gzip-bomb
    变体）。对照：`collection.rs` 用 `with_capacity(x.min(1024))` 封顶、
    `format.rs` 的 `read_string` 有 `MAX_STRING_LEN=16MB` 上限、`lz_diff.rs`
    `decode` 对 N-run / `ref_pos+len` / overflow 均有界，唯独 `unpack_cigar`
    缺失。修复：在分配前按实际 payload 大小校验 —— `op_count` 不超过
    `(raw_len - 8) / 4`，`xi_count` 不超过剩余 `raw_len - 8 - op_count*4`
    （用 `saturating_sub`/`saturating_mul` 防下溢/溢出）。

新增回归测试：
- `test_unpack_cigar_rejects_huge_op_count`（`op_count=u32::MAX`，无 op 数据）
- `test_unpack_cigar_rejects_huge_xi_count`（`xi_count=u32::MAX`，无 xi 数据）

### 验证（第 15 轮）

- `cargo test --lib pbit` 112 全绿（含两个新增测试）。
- `cargo fmt`、`cargo clippy --lib -- -D warnings` clean。

## 第 16 轮（命令层与 PAF 索引复核，无新缺陷）

对前几轮未逐段复核的剩余模块与命令层做了一轮深审，均未发现需修复的缺陷：
* `paf_index.rs`：`coord_to_i32` 对 >i32::MAX 的 PAF 坐标拒绝并 skip；
  `query_start > query_end` / `target_start > target_end` 拒绝；`query_id`
  分配用 `or_insert(next_id)` 处理重复名；全部行解析失败时 bail（非 PAF）。
  内存分配量由 PAF 行数决定（用户输入），无攻击面。
* 命令层 `create`/`append`/`append-ref`：`collect_samples_from_args` 拒绝
  `--name` 与 `-i`/`-p` 混用、`--paf` 与 `-i` 数量不匹配、重复样本名；
  `append` 用 `comp.has_sample` 拒绝归档内已存在的样本名；参考索引/名称
  解析越界均报错。`stage_work_path` 原地更新走临时文件 + 原子重命名。
* `range`/`some`/`stat`/`to-fa`：`-o` 覆盖输入防护（`ensure_outfile_distinct`
  / `same_path`）、路径穿越防护、`-s` 过滤样本不存在时报错、坐标 1-based
  含端点 → 0-based 半开转换正确、反向坐标在 CLI 层拒绝。
* `args.rs` pbit 参数：段大小/`kmer`/`min_match_len` 默认值、`-s` 短选项在
  `create`（segment-size）与 `stat`/`to-fa`（sample）分属不同子命令，无冲突。
* CIGAR 切片（`slice_cigar_by_query`）边界 D 处理、反向链 `forward_to_rc_coords`
  投影、`split_m_to_eqx` 的 M→=/X 拆分与 `=` 越界保护，均经既有测试覆盖。

### 验证（第 16 轮）

- `cargo test --test cli_pbit pbit_` 46 全绿。
- `cargo clippy --all-targets -- -D warnings` clean。

## 第 17 轮（核心库全文重读 + 文档一致性核对，无新缺陷）

对核心库做了第 3 轮之外的又一次全文重读，并交叉核对文档与 CLI 实现：
* `compressor.rs`：`slice_cigar_by_query` 边界 D 处理经完整推演正确（边界 D
  不入 sliced_ops 但仍推进 `cur_t`，使 `target_start`/`target_end` 投影包含
  其目标跨度，ref_slice 与 ops 消费一致）；`split_m_to_eqx` 的 `=`/`X`/`I`/`D`
  越界与 `rt/si` 完整消费校验齐全；`encode_segment_lzdiff` 反向段路由
  `last-saturating_sub` 防下溢；`append_reference` 的 `deltas.resize`、
  `ref_group_count` 更新、`finish` 回写 header 均正确。
* `decompressor.rs`：`delta_cache` 键含 `ref_start`/`ref_end` 正确处理 CIGAR
  delta 去重后不同 ref 切片；三层解码（LZ-diff/CIGAR）长度均对
  `delta_meta.raw_length` 校验；`decode_delta` 对 `gid`/`did` 越界、CIGAR
  `ref_start>=ref_end`、`ref_end>ref_dna.len()` 均 bail；`get_contig` 智能切片
  的 `offset` 累计与 `saturating_sub` 防下溢。
* `format.rs`：`read_ref_index`/`read_ref_table`/`read_string`/`DeltaEntry::read_from`
  对 count/len/packed_size 均有界；`PbitHeader::read_from` 校验 magic/version。
* `cigar_delta.rs`：`apply_cigar` 对 `=`/`X`/`I`/`D` 越界、`xi` 消费、`rt`
  消费均校验；第 15 轮的 op_count/xi_count 修复完好。
* 文档一致性：`docs/pbit.md` 各子命令的选项名与 `make_subcommand` 完全一致；
  `append` 的 `-o` 与归档相同时被 `stage_work_path` 拒绝（符文档"省略 -o 原地
  更新"）；`to-fa` 样本名净化（非空/无路径分隔符/非 `.`/`..`）与文档一致；
  `create` 的 `--name` 与 `-i`/`--paf` 互斥由 `collect_samples_from_args` 强制。

### 验证（第 17 轮）

- `cargo test --lib pbit` 112 全绿。
- `cargo test --test cli_pbit pbit_` 46 全绿。
- `cargo fmt`、`cargo clippy --all-targets -- -D warnings` clean。

## 已知限制（暂不改，非命令可达）

* `Decompressor` 参考层 `SequenceReader::read_sequence` 在多参考归档中按
  `contig_groups[name]`（跨所有参考的段）拼接同名 contig，语义上是"跨参考
  拼接"，与 `stat --refs` 的展示歧义同源。当前无任何 `pbit` 命令走该路径
  （`range`/`some`/`to-fa` 均用 `get_contig`/`get_sample`），仅内部测试调用，
  故暂记为 API 限制而非缺陷。

## 结论

共修复 31 处缺陷（数据安全 7 + 功能 2 + 文档 2 + 死代码 2 + 样本名冲突 1 +
展示歧义 1 + 超界 UX 1 + 溯源元数据 2 + 警告措辞 1 + 内存 DoS 5 + 健壮性 1 +
一致性/文档 1 + 报告去重 1 + 反向链压缩率 1 + 零 panic 校验一致性 2）。
第 3 轮对核心库与命令层逐行深审未再发现正确性缺陷；第 4 轮修复样本名冲突数据损坏；
第 5 轮核对了此前暂记的 `stat --refs` 展示歧义与 `range` 超界静默空输出两处记录项
并修复；第 6 轮修复 `create` 溯源元数据缺失与 contig 无段警告措辞误导；
第 7 轮补齐 `append` / `append-ref` 的溯源元数据；第 8 轮封堵 `packed_size`
无界分配的内存 DoS；第 9 轮封堵 delta 解压 gzip bomb 的内存 DoS；第 10 轮封堵
sample index 解压 gzip bomb 的内存 DoS；第 11 轮补 `to-fa` 空样本名守卫；
第 12 轮封堵 `min_match_len` 无界的段 padding 内存 DoS，并补齐 `append` 多参考
警告、`to-fa` 文档空名说明与报告去重；第 13 轮删除 `append` 中重复的不可达
`None` 分支（死代码），并修复反向链 LZ-diff 段路由以提升压缩率（解码一致故原
实现正确，改动仅为压缩率优化）；第 14 轮封堵 `append`/`append-ref` 对零
`segment_size`/`kmer_len` 畸形归档的 panic，并使 `create` 的 `min_match_len`
校验与 decompressor 绝对上限一致（`MAX_PACKED_SIZE` 提升为 `pub`）；
第 15 轮封堵 `unpack_cigar` 的 `op_count`/`xi_count` 无界分配内存 DoS；
第 16 轮对命令层与 PAF 索引复核未再发现新缺陷；第 17 轮对核心库（compressor/
decompressor/format/cigar_delta）全文重读并核对文档与 CLI 实现一致性，亦未
发现新缺陷（审核收敛）。
剩余记录项仅参考层 `read_sequence` 跨参考拼接（非命令可达）。pbit 命令族审核收敛。

验证：`cargo test --lib pbit`（112）+ `cargo test --test cli_pbit pbit_`（46）
全绿；`cargo fmt`、`cargo clippy --all-targets -- -D warnings` clean。