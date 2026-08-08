# pbit 命令族代码审核记录（2026-08-05）

对 `pgr pbit` 命令族（`create` / `append` / `append-ref` / `stat` / `range` /
`some` / `to-fa`）及核心库（`libs/pbit`：`compressor` / `decompressor` /
`format` / `collection` / `segment` / `lz_diff` / `cigar_delta` / `paf_index`）
与全部测试、文档进行多轮审核。以下仅保留有借鉴意义的结论；验证过程已精简。

## 与格式规范的一致性核对

pbit 格式（`notes/design/pbit.md` §文件格式规范 v1004）：
- Header 36 字节、Footer 24 字节、固定大小小端整数、字符串 u32 长度前缀，
  `format.rs` 逐字段读写与规范一致（含 36/24/10 字节长度断言测试）。
- `Decompressor::new` 就地校验 `ref_group_count`（header/index/delta 三处）、
  `sample_count`（header/collection）、`sample_index_offset <= footer_start`，
  损坏文件不会越界。
- 参考层复用标准 2bit 记录（`read_2bit_record`），保留 N/mask blocks，
  与 twobit 共享代码。

## 排除的疑点（安全不变量，经核验无需修复）

- `decompressor.rs` `decode_delta`：`ref_group_id`/`delta_id` 越界、CIGAR
  `ref_start >= ref_end`、`ref_end > ref_dna.len()` 均显式报错。
- `lz_diff.rs` `decode`：N-run 长度超 ref_len、`ref_pos` 计算溢出、match 区间越界、
  畸形分隔符均返回 `Err`。
- `cigar_delta.rs` `apply_cigar`：`=`/`X`/`D` 超参考长、`X`/`I` 超碱基流、
  `M` 意外出现、`xi`/`rt` 消费不一致均显式报错。
- 多参考路由：CIGAR 路径与 LZ-diff 回退均按当前参考的 `group_start/group_count`
  过滤 `ref_group_ids`，共享 contig 名的样本不会路由到错误参考。
- `append`/`append-ref` 原地更新：临时文件 + 原子重命名，失败不损坏原归档；
  `stage_work_path` 拒绝 `-o` 与输入归档相同；`TempFileGuard` 未 `disarm` 时随
  drop 删除临时文件。
- `paf_index.rs`：`coord_to_i32` 对 >i32::MAX 拒绝并 skip；`query_start > query_end`
  / `target_start > target_end` 拒绝；内存分配量由 PAF 行数决定（用户输入）。
- CIGAR 切片边界 D 处理经完整推演正确：边界 D 不入 `sliced_ops` 但仍推进 `cur_t`，
  使 `target_start`/`target_end` 投影包含其目标跨度。
- `decompressor.rs` `delta_cache` 键含 `ref_start`/`ref_end`，正确处理 CIGAR delta
  去重后不同 ref 切片。

## 已知限制（有意保留）

- `Decompressor` 参考层 `SequenceReader::read_sequence` 在多参考归档中按
  `contig_groups[name]`（跨所有参考的段）拼接同名 contig，语义上是"跨参考拼接"，
  与 `stat --refs` 的展示歧义同源。当前无任何 `pbit` 命令走该路径（`range`/
  `some`/`to-fa` 均用 `get_contig`/`get_sample`），仅内部测试调用，暂记为 API 限制。

## 修复的缺陷（根因模式）

### 数据安全（`-o` 覆盖输入）

- **各子命令 `-o` 截断输入**：`create` 截断先于读参考/样本 FASTA、PAF、`--name`
  TSV；`append`/`append-ref` 经 `stage_work_path` 复制归档时覆盖样本/参考输入；
  `range`/`some`/`stat`/`to-fa` 截断先于惰性读归档。修复：`create` 覆盖全部输入、
  `append` 对样本输入、`append-ref` 对 `-r` 列表、`range` 覆盖归档与 rgfile、
  `some` 覆盖归档与 name_list、`stat` 覆盖归档、`to-fa` 写前对每个输出路径
  `same_path` 检查。

### 数据损坏 / 样本名冲突

- **样本名冲突静默数据损坏**：样本名来自 `-i` basename 或 TSV，名冲突时
  `append_sample` 把两段样本 segment 合并进同一名字，提取返回错误序列且无报错。
  修复：`collect_samples_from_args` 单次命令内重复名报错；`compressor` 新增
  `has_sample`；`append` 追加前逐样本检查已存在名。`create` 在打开输出前拦截（不
  残留输出），`append` 在临时副本上检查、bail 后 guard 删临时文件（原归档不受影响）。

### 功能正确性（多参考 / 反向坐标 / 遮蔽还原）

- **多参考 PAF 回退路由**：`append_sample_with_paf` 的 LZ-diff 回退未按当前参考
  过滤 `ref_group_ids`，共享 contig 名在多参考归档中落到错误参考的段。修复：按
  当前参考 `group_start/group_count` 过滤后路由。
- **`range` 反向坐标 / 超界坐标静默空输出**：`chr1:100-50`、`chr1:1000-2000`
  （500bp contig）被 `get_contig` 的 `s >= e` 静默跳过，无输出无警告。修复：CLI
  层将 `start <= end` 纳入校验显式报错；`get_contig` 返回实际写入条目数
  （`Result<usize>`），`range.rs` 对写 0 条时 `log::warn!`。
- **`some`/`range` 丢失 soft-mask 小写还原**：`get_sample`（`to-fa`）已还原小写遮蔽，
  但 `get_contig`（`some`/`range` 共用）未还原，`to-fa` 与 `some`/`range` 输出不
  一致（`docs/pbit.md` 声明三者均"原样还原小写"）。修复：`get_contig` 收集样本段时
  克隆 `mask_blocks`，对正向切片按 `slice_start=s` 偏移应用 `apply_mask_blocks_at`，
  再按需反向互补（`rev_comp` 保留大小写）。

### 内存 DoS（恶意归档）

- **`packed_size` 无界分配**（`vec![0u8; meta.packed_size]`，u32 最大 ~4GB）。
  修复：`MAX_PACKED_SIZE`（256 MB），`Decompressor::new` 扫描与 `DeltaEntry::
  read_from` 两处校验。
- **delta 解压 gzip bomb**：LZ-diff 路径与 `unpack_cigar` 以 `read_to_end` 分配，
  几 KB 压缩流可解压出数 GB。修复：`MAX_DELTA_UNCOMPRESSED`（256 MB）+ `decoder.
  take(limit + 1)` 限长读取。
- **sample index gzip bomb**：`Collection::deserialize` 一次性解压 sample index。
  修复：`MAX_COLLECTION_UNCOMPRESSED`（256 MB）+ `take(limit+1)`。
- **`min_match_len` 无界 → 段 padding DoS**：`decode_delta` 每解一段触发
  `reference.resize(... + key_len)`（`key_len ≈ min_match_len - 3`）。修复：
  `Decompressor::new` 增加绝对上限（复用 `MAX_PACKED_SIZE` 256 MB）。
- **CIGAR 段计数无界分配**：`unpack_cigar` 的 `op_count`/`xi_count` 是攻击者可控
  u32，直接 `Vec::with_capacity`。修复：分配前按实际 payload 校验 —— `op_count`
  不超过 `(raw_len - 8) / 4`，`xi_count` 不超过剩余。

### 健壮性 / 零 panic 校验一致性

- **`to-fa` 空样本名生成 `.fa` 点文件**（样本名来自不可信 `.pbit` collection）。
  修复：路径穿越守卫补 `sample.is_empty()` 检查。
- **`append`/`append-ref` 对畸形归档 panic**：`Decompressor::new` 仅校验
  `segment_size`/`kmer_len` 为正数，而 `open_for_append` 复用头部值重新分段，
  `segment_sequence` 调 `chunks(0)`、`detect_rev_comp` 调 `windows(0)`，参数为 0
  时 panic。修复：显式拒绝 `segment_size == 0` 与 `kmer_len == 0`。
- **`create` 的 `min_match_len` 校验与 decompressor 不一致**：`create` 只校验
  `<= segment_size`（上限 ~2GB），可创建 decompressor 拒绝的非法归档。修复：`create`
  增加一致上限（`min_match_len <= MAX_PACKED_SIZE`），并将 `MAX_PACKED_SIZE` 提升
  为 `pub`。

### 反向链压缩率

- **反向链 LZ-diff 段路由**：`encode_segment_lzdiff` 对反向互补 contig 按正向
  `seg_idx → ref_group_ids[seg_idx]` 一一对应，未按反向倒序匹配；解码用同一映射
  故正确，但每段都错配 → delta 偏大。修复：`contig_is_rev_comp` 为真时改用
  `ref_group_ids[N-1-seg_idx]`。因 `ref_group_id` 逐段写入并随归档存储、解码按
  同一映射还原，改动仅影响压缩率不影响正确性。

### 展示歧义 / 一致性 / 溯源元数据

- **`stat --refs` 多参考展示歧义**：同名 contig 跨参考聚合无法区分所属参考。修复：
  参考数 > 1 时每行前缀参考名（`ref_name<TAB>contig<TAB>count`），单参考不变。
- **`create`/`append`/`append-ref` 溯源元数据不完整**：`cmd_line` 不含样本/参考
  输入来源。修复：追加每个样本的 `-i name:path`、可选 `-p paf`、可选 `@ref
  ref_spec`、`append-ref` 的 `-r`。
- **`append` 多参考缺警告**：未指定 `ref_spec` 时静默路由到参考 0，与 `create`/
  `resolve_ref_id` 不一致。修复：参考数 > 1 且未指定时 `log::warn!`。

### 死代码 / 警告措辞 / 文档一致性（一次性小修，已精简）

`append_sample` 重复空检查死代码、`append` 不可达 `None` 分支、contig 无段警告
措辞误导均修正；`docs/pbit.md` 与 CLI 帮助的 `--name` TSV 第 4 列（可选
`ref_name`）、`to-fa` 空名拒绝说明补齐。

## 结论

`pbit` 命令族审核完成（累计修复 32 处缺陷：数据安全 7、数据损坏 1、功能正确性 2、
遮蔽还原一致性 1、展示歧义/UX 2、溯源元数据 2、死代码/警告措辞 3、文档一致性 3、
内存 DoS 5、健壮性 1、一致性/报告去重 2、反向链压缩率 1、零 panic 校验一致性 2），
补回归测试与文档澄清，并经多轮纵深复审（`compressor`/`decompressor`/`lz_diff`/
`format`/`collection`/`segment`/`paf_index`/`cigar_delta` 全部库文件与命令层）
均未再发现新缺陷，审核收敛。剩余记录项仅参考层 `read_sequence` 跨参考拼接（非命令
可达，见"已知限制"）。
