# sd / rept / align 命令族代码审核记录（2026-08-04）

对新增命令族 sd（8 命令 + libs/sd）、rept（6 命令 + libs/pl）、
align（pgi/lastz + libs/pgi、libs/lastz、libs/fmt/lav、alignment DP）
约 8000 行代码及全部文档（docs/{sd,rept,align-pgi,align-lastz}.md）进行
审核。缺陷按发现顺序全局编号 #1–#37，按类别分组记录；#26a/26b/31a 为
外部源码复核，见「与外部参考实现的语义一致性核对」。验证概况见文末。

## 与外部参考实现的语义一致性核对

两处关键修复均对照官方源码复核，方向一致：

* #26 `psl lift` 负链坐标提升：与 UCSC kent `pslLiftSubrangeBlat.c` 的
  `liftSide` 行为一致（子范围命名约定 pgr 1-based vs kent 0-based 为
  记录在案的有意差异）。
* #31 greedy 链合并：与 FastGA `align_contigs` / `ALNchain.c` 的链化
  语义一致——同对角线纯间隔是两条独立链，仅对角线平移才缝合（pgr 的
  自身扩展）。

## 排除的疑点（经核验无需修复）

* `sd run` 的 cluster set_id 重编号值域各簇两两不相交，不可能碰撞。
* 60,423 → 75,413 数据差异来自 tncentral 库更新与编译时序，非代码 bug
  （repeat-masking.md §2.3.5 已勘误）。
* sd cluster minus 链序列提取：按 pgr PAF 正向坐标约定提取，逐碱基一致。
* wave 初始 trim 越界经几何推演与约 20 万次 fuzz 均不可达，不加防御。
* `spanr fill -n 0` 为 no-op，与设计一致，仅多一次冗余进程。
* LAV `s`/`h` stanza 含空格文件名解析错位，与 UCSC lavToPsl 一致，记录不修。

## 记录项（未改，低风险 / 待决策）

* tube 工作流对"库 vs 基因组"的结构性失效：根因是跨对角桶链被切断，结论
  基于修复前代码，syncmer/排序键修复后待真实数据重测。
* `decompose.rs` 负链投影依赖 header 与序列长度一致（cluster 内部保证）。
* cluster/cover 的 u32→i32 坐标转换（仅 >2.1 Gb 染色体溢出）。
* `run_lastz` self 模式仍构建 n×n job 列表，大目录可提前过滤。
* `syncmer.rs` 参考实现与 `collect_one_contig` 重复发射同一位置，消费方
  已 HashSet 去重，可后续合并。
* wave.rs 的 `unreachable!`/`panic!` 均为算法不变量，有测试兜底。
* s_align / sd search --engine pgi 传不支持类型时报错可读性差，不 panic。
* `fa split name` 名称碰撞（`chr(1)` 与 `chr_1`）概率极低，记录不修。

## 已知限制（有意保留）

* 子范围命名 pgr 1-based vs kent 0-based（#26b 决策）：pgr 生成端/消费端
  自洽，直接消费 UCSC/blat 生态子范围名时需先确认语义。
* s-kmer 对染色体尾部重复保守丢弃：Profex `-z` 不输出末 run 深度，有阈值
  时无法区分唯一尾与重复尾（与 anchr 参考管线一致）。
* 单 contig > 4.3 Gb 的 pgi 索引：pos 为 u32，超长单 contig 不被支持。

## 修复的缺陷（共 37 处）

### 崩溃 / 越界 / 溢出（Zero Panic，10 处）

4. **sd/run.rs 解析 elem.bed 短行越界**：直接取 `f[4]`。修复：加
   `f.len() < 8` 检查（与 cover.rs 一致）。
5. **sd decompose 负链投影 usize 下溢**（畸形 header）。修复：拒绝
   end < start，投影 saturating。回归 `malformed_header_does_not_panic`。
6. **lav d stanza 边界差一越界**。修复：守卫改 `+ 6`。回归
   `truncated_d_stanza_errors_not_panics`。
7. **构造 .pgi/.hv 头容量溢出 panic/OOM**（未校验 n_records/n_contigs）。
   修复：头解析校验 + `try_reserve_exact`。回归 3 个 crafted 测试。
8. **e-align span 过滤 `(t_end - t_start) as usize` 回绕**。修复：i64
   运算 `.max(0)` 再转 usize。
9. **lav `l` 行负跨度回绕成超大 block**。修复：t_end < t_start 等报
   InvalidData。回归 `negative_span_l_line_rejected`。
10. **pgi build `positions.len() as u32` 静默截断**（>42 亿记录）。修复：
    `payloads.len() <= u32::MAX` 防御检查。
29. **非 UTF-8 临时目录路径 `to_str().unwrap()` panic**。修复：
    `io::path_to_str` 友好报错。
35. **`align_banded_local` 序列长度悬殊时 DP 数组越界**。修复：j_lo/j_hi
    与对角带求交、空交集跳行。回归 `unbalanced_lengths_do_not_panic`。
36. **lav `l` 行极值坐标 `-1` 下溢/跨度比较溢出**。修复：`checked_sub`。
    回归 `extreme_l_line_values_do_not_panic`。

### 功能正确性 / 算法（13 处；#1-3 为重大索引/链算法缺陷）

1. **（重大）pgi 索引 k-mer key 与位置错配**（2 Mb 随机基因组 39% 错配、
   self 比对 101 条伪块）。修复：pending 去重、flush 按位置重算 key、RC
   用 `rc_key`。回归 `index_records_match_sequence_positions`。
2. **（重大）tube 排序键 anti/bucket 溢出**（>8 Mb 基因组失效）。修复：
   anti/bucket 扩到 32 位。回归
   `tube_sort_key_supports_large_anti_coordinates`。
3. **（重大）tube 排序键负对角线回绕**（>64 Mb 间距失效）。修复：
   `BUCK_OFF = 1 << 26`。回归深负对角线两个测试。
12. **cluster 重叠 union 漏连嵌套区间**。修复：扫描时跟踪最大右端。回归
    `nested_overlapping_intervals_form_one_cluster`。
13. **sd cluster 去重键忽略链向/物种**（回文倒位拷贝被折叠）。修复：键加
    strand。回归 `same_coordinates_on_opposite_strands_are_distinct_copies`。
14. **s-align 漏做带点 contig 名映射**（spanr 截断，`fa mask` 失配）。
    修复：复用 chrom.sizes 映射。回归 `command_rept_s_align_dotted_name`。
15. **Profex `-z` 坐标右端多 +1 + e-kmer 染色体尾部丢失**。修复：end 不再
    +1；无阈值时用染色体长度闭合尾 run。回归
    `command_rept_e_kmer_tandem_coordinates`。
16. **sd cluster/run 不支持普通 gzip**（生成垃圾 `.loc`）。修复：非 BGZF
    先解压到临时文件。回归 `command_sd_run_gzipped_genome`。
17. **align lastz 省略 query 未启用 self 模式**。修复：传 `self_mode`。
    回归 `command_align_lastz_omitted_query_is_self`。
26. **`psl lift` 负链外层坐标提升错误（违反 UCSC 约定）**。修复：
    `qStart/qEnd += start_0`、`qStarts += (size - end_0)`，夹具修正。
    回归 `test_lift_minus_strand_forward_coordinates`。
27. **s-align/e-align soft-mask 警告误报 N gap**。修复：`has_soft_mask`
    只扫 lowercase。回归 `soft_mask_detection_ignores_n_gaps`。
31. **greedy 链合并导致倒位 SD 漏检**。修复：合并条件加
    `|diagA − diagB| > 0`。回归 `command_sd_search_pgi_inverted_repeat`。
37. **pgi merge 频率过滤两侧边界不一致**（`== freq` 处理与 FastGA 不符）。
    修复：A/B 侧统一 `>= freq` 跳过、`< freq` 计入。回归
    `freq_boundary_drops_exact_freq_on_reference_side`、
    `exact_freq_query_entries_are_absent_not_range_killers`。

### 输入校验 / 静默错误（3 处）

11. **repeat.rs 两处 `map_while(Result::ok)` 吞 IO 错误**。修复：
    `let line = line?;` 传播错误。
32. **e-align PSL 过滤静默跳过畸形行**。修复：补 `log::warn!`。
34. **decompose 对解析失败的 FASTA 头静默丢弃**。修复：补 `log::warn!`。

### 外部工具与参数 / CLI / 文档（11 处）

18. **lastz 静默失败**（只打日志返回 Ok）。修复：统计失败数并 bail。
19. **lastz 失败原因被吞**（status 丢 stderr）。修复：`cmd.output()` 记录
    首个失败的 stderr。
20. **参数校验缺失/不一致**（`--min-identity` 范围、kmer/window/parallel/
    minscore 正值有限性）。修复：统一校验，帮助同步 "(0, 1]"。
21. **trf 特殊字符文件名找不到**。修复：`sanitize_filename(chr)`。回归
    `command_rept_trf_special_chars`。
22. **sd search/cross `--preset` 默认值未注册**。修复：
    `.default_value("set01")`。回归 `command_sd_search_lastz_default_preset`。
23. **sd run --engine lastz --preset 拼装错误**。修复：`Vec<String>` +
    `$[preset_args]` 展开。回归 `command_sd_run_lastz_preset_parses`。
24. **噪音与帮助文本多处小修**：lav mask stanza 静默、`#` 元数据行跳过、
    lastz `[multiple]`/`-s` 修正、align.md 示例输出修正、pgi 帮助默认
    syncmer 修正。
25. **文档一致性**：rept.md 补 e-align；soft-mask 说明；`.pgi` 命名；
    sd.md lastz 单序列/纯文本约束；TnCentral 路径。
28. **lastz 单序列约束帮助/文档未同步**。修复：四处补齐。
30. **e-align identity 定义未说明**（gap-compressed）。修复：补文档。
33. **主帮助 rept 子命令列表漏 e-align/s-align**。修复：补齐。

## 验证

* 引擎交叉验证：pgi 与 lastz 检出同一对 1200 bp 拷贝（边界修剪 2 bp）；
  倒位重复经 chainnet 后两条 `-` 链保留，`sd run` 输出坐标正确的
  elementary SD。
* 端到端坐标：`rept trf` 输出 "101-1100"、`rept s-align` 输出 "501-2900"
  （1-based 全覆盖）；多 contig s-kmer 编号与 chr.sizes 一致；e-kmer
  `--keep-index` 缓存复用/失效正确；普通 gz 与 BGZF 全流程通过。
* 鲁棒性：截断/负跨度 lav、越界 PAF、垃圾 BED/.pgi、构造头、空输入、
  全 N、极值参数、短行等畸形输入全部友好报错或空输出，零 panic。
* `cargo test` 全量 995 通过；`cargo fmt --check` 与 `cargo clippy
  --all-targets -- -D warnings` 干净；sd/rept/align 四组 CLI 端到端
  （38 测试）全绿。
