# sd / rept / align 命令族代码审核记录（2026-08-04）

对新增命令族 sd（8 命令 + libs/sd）、rept（6 命令 + libs/pl）、
align（pgi/lastz + libs/pgi、libs/lastz、libs/fmt/lav、alignment DP）
约 8000 行代码及全部文档（docs/{sd,rept,align-pgi,align-lastz}.md）进行
审核。缺陷按类别分组记录；`psl lift` 与 greedy 链合并两处关键修复对照
官方源码复核（见「与外部参考实现的语义一致性核对」）。验证概况见文末。

## 与外部参考实现的语义一致性核对

两处关键修复均对照官方源码复核，方向一致：

* `psl lift` 负链坐标提升：与 UCSC kent `pslLiftSubrangeBlat.c` 的
  `liftSide` 行为一致（子范围命名约定 pgr 1-based vs kent 0-based 为
  记录在案的有意差异）。
* greedy 链合并：与 FastGA `align_contigs` / `ALNchain.c` 的链化
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

* 子范围命名 pgr 1-based vs kent 0-based（记录在案）：pgr 生成端/消费端
  自洽，直接消费 UCSC/blat 生态子范围名时需先确认语义。
* s-kmer 对染色体尾部重复保守丢弃：Profex `-z` 不输出末 run 深度，有阈值
  时无法区分唯一尾与重复尾（与 anchr 参考管线一致）。
* 单 contig > 4.3 Gb 的 pgi 索引：pos 为 u32，超长单 contig 不被支持。

## 修复的缺陷（共 37 处）

### 崩溃 / 越界 / 溢出（Zero Panic，10 处）

**sd/run.rs 解析 elem.bed 短行越界**：直接取 `f[4]`。修复：加
   `f.len() < 8` 检查（与 cover.rs 一致）。
**sd decompose 负链投影 usize 下溢**（畸形 header）。修复：拒绝
   end < start，投影 saturating。回归 `malformed_header_does_not_panic`。
**lav d stanza 边界差一越界**。修复：守卫改 `+ 6`。回归
   `truncated_d_stanza_errors_not_panics`。
**构造 .pgi/.hv 头容量溢出 panic/OOM**（未校验 n_records/n_contigs）。
   修复：头解析校验 + `try_reserve_exact`。回归 3 个 crafted 测试。
**e-align span 过滤 `(t_end - t_start) as usize` 回绕**。修复：i64
   运算 `.max(0)` 再转 usize。
**lav `l` 行负跨度回绕成超大 block**。修复：t_end < t_start 等报
   InvalidData。回归 `negative_span_l_line_rejected`。
**pgi build `positions.len() as u32` 静默截断**（>42 亿记录）。修复：
    `payloads.len() <= u32::MAX` 防御检查。
**非 UTF-8 临时目录路径 `to_str().unwrap()` panic**。修复：
    `io::path_to_str` 友好报错。
**`align_banded_local` 序列长度悬殊时 DP 数组越界**。修复：j_lo/j_hi
    与对角带求交、空交集跳行。回归 `unbalanced_lengths_do_not_panic`。
**lav `l` 行极值坐标 `-1` 下溢/跨度比较溢出**。修复：`checked_sub`。
    回归 `extreme_l_line_values_do_not_panic`。

### 功能正确性 / 算法（13 处，含 3 处重大索引/链算法缺陷）

**（重大）pgi 索引 k-mer key 与位置错配**（2 Mb 随机基因组 39% 错配、
   self 比对 101 条伪块）。修复：pending 去重、flush 按位置重算 key、RC
   用 `rc_key`。回归 `index_records_match_sequence_positions`。
**（重大）tube 排序键 anti/bucket 溢出**（>8 Mb 基因组失效）。修复：
   anti/bucket 扩到 32 位。回归
   `tube_sort_key_supports_large_anti_coordinates`。
**（重大）tube 排序键负对角线回绕**（>64 Mb 间距失效）。修复：
   `BUCK_OFF = 1 << 26`。回归深负对角线两个测试。
**cluster 重叠 union 漏连嵌套区间**。修复：扫描时跟踪最大右端。回归
    `nested_overlapping_intervals_form_one_cluster`。
**sd cluster 去重键忽略链向/物种**（回文倒位拷贝被折叠）。修复：键加
    strand。回归 `same_coordinates_on_opposite_strands_are_distinct_copies`。
**s-align 漏做带点 contig 名映射**（spanr 截断，`fa mask` 失配）。
    修复：复用 chrom.sizes 映射。回归 `command_rept_s_align_dotted_name`。
**Profex `-z` 坐标右端多 +1 + e-kmer 染色体尾部丢失**。修复：end 不再
    +1；无阈值时用染色体长度闭合尾 run。回归
    `command_rept_e_kmer_tandem_coordinates`。
**sd cluster/run 不支持普通 gzip**（生成垃圾 `.loc`）。修复：非 BGZF
    先解压到临时文件。回归 `command_sd_run_gzipped_genome`。
**align lastz 省略 query 未启用 self 模式**。修复：传 `self_mode`。
    回归 `command_align_lastz_omitted_query_is_self`。
**`psl lift` 负链外层坐标提升错误（违反 UCSC 约定）**。修复：
    `qStart/qEnd += start_0`、`qStarts += (size - end_0)`，夹具修正。
    回归 `test_lift_minus_strand_forward_coordinates`。
**s-align/e-align soft-mask 警告误报 N gap**。修复：`has_soft_mask`
    只扫 lowercase。回归 `soft_mask_detection_ignores_n_gaps`。
**greedy 链合并导致倒位 SD 漏检**。修复：合并条件加
    `|diagA − diagB| > 0`。回归 `command_sd_search_pgi_inverted_repeat`。
**pgi merge 频率过滤两侧边界不一致**（`== freq` 处理与 FastGA 不符）。
    修复：A/B 侧统一 `>= freq` 跳过、`< freq` 计入。回归
    `freq_boundary_drops_exact_freq_on_reference_side`、
    `exact_freq_query_entries_are_absent_not_range_killers`。

### 输入校验 / 静默错误（3 处）

**repeat.rs 两处 `map_while(Result::ok)` 吞 IO 错误**。修复：
    `let line = line?;` 传播错误。
**e-align PSL 过滤静默跳过畸形行**。修复：补 `log::warn!`。
**decompose 对解析失败的 FASTA 头静默丢弃**。修复：补 `log::warn!`。

### 外部工具与参数 / CLI / 文档（11 处）

**lastz 静默失败**（只打日志返回 Ok）。修复：统计失败数并 bail。
**lastz 失败原因被吞**（status 丢 stderr）。修复：`cmd.output()` 记录
    首个失败的 stderr。
**参数校验缺失/不一致**（`--min-identity` 范围、kmer/window/parallel/
    minscore 正值有限性）。修复：统一校验，帮助同步 "(0, 1]"。
**trf 特殊字符文件名找不到**。修复：`sanitize_filename(chr)`。回归
    `command_rept_trf_special_chars`。
**sd search/cross `--preset` 默认值未注册**。修复：
    `.default_value("set01")`。回归 `command_sd_search_lastz_default_preset`。
**sd run --engine lastz --preset 拼装错误**。修复：`Vec<String>` +
    `$[preset_args]` 展开。回归 `command_sd_run_lastz_preset_parses`。
**噪音与帮助文本多处小修**：lav mask stanza 静默、`#` 元数据行跳过、
    lastz `[multiple]`/`-s` 修正、align.md 示例输出修正、pgi 帮助默认
    syncmer 修正。
**文档一致性**：rept.md 补 e-align；soft-mask 说明；`.pgi` 命名；
    sd.md lastz 单序列/纯文本约束；TnCentral 路径。
**lastz 单序列约束帮助/文档未同步**。修复：四处补齐。
**e-align identity 定义未说明**（gap-compressed）。修复：补文档。
**主帮助 rept 子命令列表漏 e-align/s-align**。修复：补齐。

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

## 复核 2（2026-08-04 后续轮次）

### 修复的缺陷（2 处）

**（崩溃）crafted .pgi 记录 contig id 越界 panic**：构造索引的
occurrence 记录携带超出 contig 表的 cid 时，`emit_entry_hits` 的
`b.contigs()[bc]` 直接越界 panic（Zero Panic 违反）。三个读取路径此前
只校验头部不校验记录体。修复：`PgiIndex::read` / `PgiStream::next_batch`
逐记录校验 `cid < n_contigs` 且 `pos + k <= contig len`；
`PgiMmap`（惰性解码）在 `emit_entry_hits` 解码命中记录时同步校验，
报友好错误。回归 3 个测试：`crafted_record_contig_rejected_not_panic`、
`mmap_merge_rejects_out_of_range_contig`、
`command_align_pgi_crafted_index_errors_not_panics`。

**（功能）相邻链合并把两条独立同源对缝成嵌合链，SD 命中丢失**：
多拷贝家族中两条拷贝对的对角线差可在 band 内（如 56 bp）且两轴间隔
均在 merge_gap 内（如 3.6 kb），纯几何 merge 会将其缝成一条跨两段真实
匹配 + 随机中段的嵌合链；扩展出的嵌合块身份 ~72% 被 SD 过滤丢弃，
两条真实命中（及下游一个拷贝的 CORE 标记）随之丢失。几何条件无法区分
"同源块种子缺口"与"两条独立块"（两者形状完全一致），必须用序列判定。
修复：`merge_adjacent_chains` 增加可选序列参数；两侧间隔均非空时要求
中段 banded 对齐身份 ≥ 0.9 且 query 覆盖 ≥ 0.9 才合并（单轴插入的
IS 缝合场景一侧间隔为空，跳过检查）。回归：单元测试
`merge_requires_homologous_middle_with_sequences`（随机中段不合并、
同源中段合并）、CLI 测试
`command_sd_search_pgi_multi_copy_close_diagonals`（4 拷贝 × 2 方向
共 12 条命中全保留）。随机化 6 组端到端 trial 全部拷贝 CORE 覆盖。

### 文档/帮助一致性（4 处）

* rept.md 仍写 "`align` variants are planned"（e-align/s-align 早已
  存在）→ 改为现况描述。
* rept.md "All four emit runlist JSON" → "All five"（trf 在内共 5 命令）。
* rept.md e-align 空 "### Dependencies" 章节删除。
* align-pgi.md `--freq` 语义写 "more than this many times"，代码与帮助
  均为 "at least this many times"（`>= freq`）→ 文档修正。
* `sd run` 帮助/文档补齐：`--preset` 默认 set01（与 search/cross 一致，
  值透传给 search 的默认）、`--min-identity (0, 1]`、lastz 引擎需
  单序列 FASTA 的说明。

## 复核 2 验证

* 修复后 trial 场景 `sd search` 输出 12/12 条命中（此前 8/12），
  `sd run` elementary BED 四拷贝均有 CORE 行。
* MG1655 真实基因组 `sd search --engine pgi`：229 条 ≥1 kb/≥90% 命中、
  136 个互惠对（43 条单侧为边界修剪所致），全部长度/身份自洽；无新增
  嵌合块。修复仅收紧合并，不可能减少合法 ≥1 kb 命中。
* 极端参数（k/smer/window/band/merge-gap/freq/query-depth 0、
  k=65、max-period 0、pm 0、min-depth 0）与全 N/小序列/带 N 基因组
  全部友好报错或空输出，零 panic。
* `cargo test` 全量 1201 通过；`cargo fmt --check` 与 `cargo clippy
  --all-targets -- -D warnings` 干净。

## 复核 3（2026-08-04 后续轮次）

### 修复的缺陷（4 处）

* `sd run` 合并 elementary BED 时按 `read_dir` 顺序枚举 cluster 文件，
  set_id 全局重编号与输出行序依赖文件系统枚举顺序，跨运行不确定。修复：
  按 cluster 文件名的**数值**编号排序（词法排序会把 cluster_10 排在
  cluster_2 前）。回归验证：10 家族基因组输出 set_id 严格 1..10 升序。
* `sd align`（`chainnet_to_paf`）按 `read_dir` 顺序迭代 MAF 文件，
  PAF 输出行序不确定。修复：排序后再合并。
* `sd search --engine lastz` 按 `read_dir` 顺序迭代 LAV 文件，输出 PSL
  行序不确定。修复：排序后再转换。
* `search_lastz::decompress_if_gz` 把解压输出统一命名为
  `{base}.plain.fa`（base 取首个 `.` 段）：嵌套目录中同名 `.fa.gz`
  会解压到同一路径，后写的静默覆盖先写的，两个 job 比对同一份错误序列。
  修复：同一次调用内用 HashSet 去重，重复 basename 追加输入序号后缀；
  序号确定性保证 self 模式下 target/query 两次调用生成相同路径。
  回归 `decompress_colliding_basenames_stay_distinct`。

### 记录项（未改，低风险）

* `rept e-kmer`/`s-kmer` 的 `--fill-kmer` 以 `usize as i32` 传入
  `IntSpan::fill`；超 i32 值静默截断为负 → fill 变 no-op（无 panic，
  且 `excise` 同理安全）。极端参数属用户误用，行为安全，记录不修。
* `s-align` 的 `--min-depth` 以 `usize as u32` 传入深度阈值；超 u32 值
  截断。同上，记录不修。
* `sd cluster` 的同染色体重叠合并只按 chrom 名（物种前缀已剥离）分组：
  跨基因组 PAF 在两端基因组 contig 名与文件 stem 均相同时会把不同基因
  组的同名区间并簇。`sd cluster` 文档仅面向自比对 PAF，记录不修。

## 复核 3 验证

* `sd run` 两次运行输出逐字节一致；`sd search --engine lastz` 与
  `sd align` 输出逐字节一致（多次运行）。
* `sd run` 10 个重复家族端到端：20 条 elementary 行、set_id 1..10
  升序、全部 CORE。
* 全量 1202 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 4（2026-08-04 后续轮次）

### 修复的缺陷（3 处）

**（功能）`psl lift` 的 `parse_subrange` 误切含 `.`/`:` 的 contig 名**：
窗口名 `{contig}:{start}-{end}` 经共享 `Range` 解析器时，`NC_000913.1:1-200`
被读成 name="NC_000913" + chr="1"、`chr1:alt:1-200` 被读成 chr="alt"，
`lift_query` 在 sizes 表里查错键、静默跳过提升（仅 warn）。s-align 的
coverage 于是来自未提升的窗口坐标（带点名用例恰好在 runlist 还原后"看起来
对"，坐标实际是窗口空间）。修复：`parse_subrange` 改为取最后一个
`:`+数字后缀切分，前缀整体作为 contig 名；回归
`parse_subrange_keeps_dotted_and_colon_contigs`。

**（功能）`s-align` 安全名改写仍用 `split_once(':')`**：coverage.rg 行
`chr1:alt:1-200` 被切成 name="chr1"，带 `:` 的 contig 输出键被截断成
"alt"。修复：改用 `parse_subrange` 解析再写占位名。回归
`command_rept_s_align_colon_name`；`command_rept_s_align_dotted_name`
断言从"含 '-' 即可"加强为精确坐标（301-800,1001-1500，即真实重复区）。

**（帮助文本）`pgr align` 的 about 写 "into PSL blocks"，但 lastz 子命令
输出 LAV**。改为 "Aligns genomes or .pgi indexes"。

## 复核 4 验证

* `s-align` 对 `>chr1:alt` 与 `>NC_000913.1` 的输出键完整、区间与真实
  重复拷贝坐标逐碱基一致。
* `psl lift` 既有 4 个 CLI 测试与全部单元测试通过。
* 全量 1205 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 6（2026-08-04 后续轮次）

### 修复的缺陷（1 处）

**（静默错误）`sd search`/`sd cross`/`sd run` 传入 `.pgi` 索引**：pgi
引擎对 `.pgi` 输入不做扩展（无序列），输出块全部 0 分，SD 过滤后静默返回
空结果。修复：`pgi_to_hits`/`lastz_to_hits` 前置拒绝 `.pgi` 输入（magic
或扩展名），报友好错误。回归 `command_sd_search_rejects_pgi_input`。

### 记录项（未改）

* 顶层路径为 `.pgi` 扩展名的目录会被 `is_pgi_input` 误判拒绝（目录名恰好
  以 .pgi 结尾）。概率极低，记录不修。

## 复核 6 验证

* `sd search`/`sd run` 传 `.pgi` 均报 "sd search/cross needs genome
  FASTA" 而非空输出。
* 全量 1206 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 7（2026-08-04 后续轮次）

### 测试补充（未发现新缺陷）

* `merge_checks_minus_strand_middle_in_rc_space`：覆盖链合并中段同源检查的
  负链分支（b 坐标为 RC 空间，需先对正向中段做 rev_comp），含同源中段
  合并与随机中段分离两个断言。此前该分支无直接测试。

### 交叉验证

* MG1655 `sd search --engine pgi`：229 条命中（min block 1007 bp、min
  identity 0.959），136 个拷贝对（93 双方向、43 单侧为边界修剪），与
  历轮一致。
* `sd cross` 与 `sd run` 多次运行输出逐字节一致（确定性）。
* 全量 1207 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 5（2026-08-04 后续轮次）

### 文档（1 处）

* `docs/align-pgi.md` 的 `--merge-gap` 说明补齐序列校验语义：两侧间隔均
  非空时合并还要求中段同源（banded 对齐验证），近对角线的独立拷贝对保持
  分离（对应复核 2 的嵌合链修复）。

### 专项复核（未发现新问题）

* tube/greedy 双工作流在 5 拷贝密集重复、>20 kb 长程重复上均输出一致的
  全部互惠块（20/20、2/2），无 panic。
* 随机化 5 组 4 拷贝端到端 `sd run`：全部拷贝 CORE 覆盖。
* `e-align` 对库的负链（RC）拷贝输出正确靶区间；含 `#` 的 cluster 头、
  含 `:` 的 PAF 名、带未知 tag 的 PAF 均友好处理。
* `sd run --engine lastz --preset set03` 端到端通过（preset 转发与默认
  一致）。
* 全量 1205 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 8（2026-08-04 后续轮次）

### 测试补充（未发现新缺陷）

* `randomized_single_pass_matches_reference`：30 组随机序列（随机位置/长度
  N 段，默认 k=40/syncmer 8/5）上 `collect_one_contig` 与参考
  （syncmer_dna + rolling keys）逐记录一致。

### 交叉验证

* `e-kmer --keep-index` 库文件被替换（mtime 更新）后缓存正确失效并重建。
* `docs/sd.md` 全文复查与实现一致。
* 全量 1208 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 9（2026-08-04 后续轮次）

### 专项复核（未发现新问题）

* `wave_extend`/`forward_wave` 为无生产调用方的导出 API（仅自身测试），
  扩展路径实际使用 `local_alignment`；`extend_tube` 的 `alow` 推进（None
  时 +BUCK_ANTI、Some 时到 eant）无死循环；`forward_wave_mid` 对越界
  mid-line 的 x/y 有 dead-cell 守卫。
* MG1655 `e-kmer`（TnCentral 库）输出 48 个区间，与 docs 记录一致；
  `s-align` 输出 1457 区间 / 244,460 bp（广谱自比对语义，无参考数值冲突）。
* 全量 1208 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 10（2026-08-04 后续轮次）

### 修复的缺陷（1 处，性能语义）

**`align pgi --parallel` 未约束自动索引构建的 rayon 并行度**：`resolve_side`
（内部 `build_from_seqs` → `radix_sort_u128_par`）在自定义线程池创建前
执行，索引构建走全局 rayon 池，`--parallel N` 只约束 merge/扩展阶段。
文档承诺 "--parallel: rayon thread count"，行为不一致。修复：把从
`resolve_side`（索引构建）到 merge/扩展的整个流程移入 `pool.install`，
`--parallel` 现约束整个命令的 rayon 用量（`sd search --engine pgi` 与
`rept e-align` 经由 `align pgi` 同步受益）。

## 复核 10 验证

* `align pgi` 15 个 CLI 测试全通过；`-p 2` 端到端输出与修复前一致。
* 全量 1208 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 11（2026-08-04 后续轮次）

### 专项复核（未发现新问题）

* `--query-depth 0`：lastz 引擎正常执行（视为无深度限制），`sd search` 与
  `align lastz` 均找到全部命中。
* gz 基因组走 `s-align`/`trf` 全流程正常；`trf` 对 4×300 bp 串联重复输出
  `1-1200`（全覆盖）。
* 空输出、`--min-len` 过滤、`fa size` 对 2bit 报错的链路均为友好错误。

### 记录项（未改）

* `sd search --engine pgi` 接受 `.2bit` 输入（`align pgi` 原生支持），但
  下游 `sd align`/`sd run` 的 chainnet 需要 FASTA，2bit 会在
  `fa size` 步骤报错（外层 run_cmd 只显示失败命令、不含根因）。文档仅
  承诺 FASTA；2bit 部分支持是既有行为，记录不修。

## 复核 11 验证

* 全量 1208 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 12（2026-08-04 后续轮次）

### 专项复核（未发现新问题）

* 真实 Profex `-z` 输出（"Read 1:" 头、`start - end (depth)` 行、末行裸
  start）与 `run_profex_per_chr` 的正则/尾部处理逐行匹配；s-kmer 在
  min_depth=2 下正确保留深度 2 区间、丢弃尾 run。
* `sd search --engine pgi` 与 `align pgi` 多次运行输出逐字节一致。
* `--min-shared 1`（极值）正常输出大量块、`--min-shared 40`（精确）正常，
  无 panic。

### 记录项（未改）

* `align lastz --lastz-args` 的值以 `-` 开头时需用 `--lastz-args=<val>`
  形式（clap 对空格形式的值为标准行为）；帮助文本未提示该写法。

## 复核 12 验证

* 全量 1208 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 13（2026-08-04 后续轮次）

### 专项复核（未发现新问题）

* `extend_chain` 为仅测试使用的公开 API（生产路径走 `psls_from_hits` 的
  窗口扩展），无死代码风险。
* 多 contig `e-kmer`/`s-kmer`：每个 contig 正确出现在 runlist 中，重复
  拷贝坐标精确（如 c1 1201-1500/1701-2000、c2 801-1100）。
* 10 拷贝 cluster 的 `sd decompose`：10 行输出、正/负链家族正确分组为
  2 个 set，11 ms 完成。
* 多 contig `trf`：c1/c2 各自报告串联重复区间（1-600、1-1200）。

## 复核 13 验证

* 全量 1208 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 14（2026-08-04 后续轮次）

### 专项复核（未发现新问题）

* 三个家族的 CLI 参数审计：`get_one`/`get_flag` 全部对应已定义的 Arg
  （类型、默认值、possible values 一致），无缺失参数路径。
* 1 Mb 人工基因组（8 个重复家族 × 2 拷贝）`sd run` 端到端 2.2 s
  （debug 构建），输出 16 行 elementary、全部 CORE，无 panic。

## 复核 14 验证

* 全量 1208 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 15（2026-08-04 后续轮次）

### 专项复核（未发现新问题）

* pgi 与 lastz 引擎在同一合成基因组上检出一致：两引擎各 4 条命中，覆盖
  相同两个重复家族，坐标差异仅边界修剪 4–8 bp。
* `e-align --keep-index`：首次运行生成 `genome.pgi`/`lib.pgi`，后续运行
  （按 mtime 验证）复用缓存索引且输出逐字节一致。

## 复核 15 验证

* 全量 1208 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 16（2026-08-04 后续轮次）

### 专项复核（未发现本族问题）

* 三个家族在 **release 模式**下全部通过：pgi 52、sd 13、pl 1、
  alignment 23 个 lib 测试 + 42 个 CLI 集成测试全绿。

### 记录项（范围外观察，未改）

* `cargo test --release` 全量有 1 个既有失败：
  `libs::paf::cigar::tests::test_invalid_op_panics`。`CigarOp::new` 对非法
  op 用 `debug_assert!` + release 回退为 'M'（doc comment 已声明该设计），
  `#[should_panic]` 测试只在 debug 下成立。属 paf 家族、非本审计范围，
  且非本次改动引入（git diff 为空）；建议 paf 审计轮次将测试加
  `#[cfg(debug_assertions)]` 或改为断言 release 回退行为。

## 复核 16 验证

* 本族 debug 全量 1208 测试通过；`cargo fmt --check` 与 `cargo clippy
  --all-targets -- -D warnings` 干净。

## 复核 17（2026-08-04 后续轮次）

### 专项复核（未发现新问题）

* `--keep-index` 遇同名目录：报 "could not open genome.pgi: is a
  directory"（友好错误）。
* `sd cover` 贪心平局（两个相同覆盖集合）：确定性选最后一个，输出正确
  CORE 标记。
* `sd run -o` 指向已存在文件：报 "File exists"；非法 `--preset`：clap
  拒绝。
* 低复杂度基因组（长同聚物 + 真实重复）：freq 过滤正确丢弃低复杂度
  k-mer，`sd run` 仅检出真实重复拷贝（15000-16189 / 17200-18389）。

## 复核 17 验证

* 全量 1208 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 18（2026-08-04 后续轮次）

### 综合回归验证（未发现新问题）

* 全部关键回归测试批量通过：11 个 lib 回归（crafted pgi、mmap merge、
  链合并正/负链、decompress 碰撞、parse_subrange、randomized syncmer）
  + cli_sd 15 + cli_rept 15 + cli_align_pgi 15 + cli_align_lastz 1。
* 复核 10 的 `--parallel` 重构后日志行为正常（built reference/query
  index、wrote N blocks 均按预期输出）。

## 复核 18 验证

* 全量 1208 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 19（2026-08-04 后续轮次）

### 回归核对（未发现新问题）

* 首轮审核的 37 处修复回归测试抽查全通过：pgi 索引 key/位置、
  tube 排序键（大 anti/深负对角线/跨 contig）、lav 负跨度、
  `sd search` 倒位重复/gz 基因组/lastz preset、`rept` 特殊字符/带点名等。
* MG1655 `sd run` 端到端：41 s（debug 构建）完成，117 行 elementary、
  111 行 CORE，管线正常。

## 复核 19 验证

* 全量 1208 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 20（2026-08-05 后续轮次）

### 收尾（未发现新问题）

* 全量改动 diff 复查：36 个文件（本族 25 个 + rg/runlist 会话 11 个），
  本族改动均为历轮记录项的对应实现，无遗留调试代码或未完成片段。
* `sd search` 默认 `--parallel 4`（文档一致）与 `align pgi` 默认 8 为
  各自文档化的既有差异。

## 复核 20 验证

* 全量 1208 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 21（2026-08-05 后续轮次）

### 修复的缺陷（1 处，功能）

**（功能）倒位拷贝间隔 < max_gap 时 greedy 链循环把两条互惠链并成嵌合
链，SD 对完全漏检**：两条互惠链（同一倒位对的两个方向）落在同一对角线上，
拷贝间隔小于 max_gap 时其种子在 greedy 循环内直接连成一条链（绕过
复核 2 修复的 merge_adjacent_chains 的 `|diag|>0` 守卫），扩展出一个横跨
两个拷贝 + 间隙的嵌合块，身份被稀释到 SD 阈值以下被过滤，整对丢失。
最小复现：1200 bp 倒位对间隔 800 bp → 修复前 `sd search` 输出 0 条命中、
一条 3190 bp / 身份 0.879 的嵌合块；修复后输出 2 条干净命中（1183 bp、
100% 身份）。修复：greedy 循环在"双侧种子间隙 ≥ 200 bp"时用中段同源
检查门控（不通过则闭合当前链、以该种子起新链）；间隙 < 200 bp 不检查
（对 ≥1000 bp 块，200 bp 随机间隙的嵌合身份 ≥ 0.909，高于 SD 阈值，
不会静默漏检，且保持稠密种子流廉价）。回归
`command_sd_search_pgi_close_inverted_repeat`。

### 验证

* 随机化 6 组 4 拷贝（间隙 300-900 bp、混合链向）`sd run`：全部拷贝
  CORE 覆盖（此前间隙 < max_gap 的倒位对会漏检）。
* tube workflow 同类基因组天然分开两条互惠链（anti 间隙 > BREAK），
  未受影响。
* 全量 1209 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 22（2026-08-05 后续轮次）

### 复核 21 修复的真实数据验证（无回退）

* MG1655 `sd search --engine pgi`：229 → 232 条命中（+8 新增、-5 移除）。
  8 条新增为 3233 bp 嵌合块的正确拆分（各拆成两条 ~1700 bp 干净块）；
  5 条移除的 q/t 区域在新输出中仍各有 ≥1000 bp 命中覆盖（区域仍被标注
  为 SD），无覆盖丢失。耗时 12.0 s vs 修复前 11.8 s（开销可忽略）。
* tube workflow 在间隙 300-900 bp 的多组倒位基因组上均检出真实配对。

### 记录项（未改，极端参数）

* `--max-gap` 调大（如 10000）时，greedy 循环的 off-band 忽略规则会把
  后续不同对角线的种子全部忽略（属于"另一条管"），远距重复家族可能整体
  丢失。该行为是 off-band-ignore 的既有设计（默认 1000 下正确），非本
  轮修复引入，记录不修。

## 复核 22 验证

* 全量 1209 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 23（2026-08-05 后续轮次）

### 复核 21 修复对非自比对路径的回归确认（未发现新问题）

* `rept e-align` 输出与修复前逐字节一致（2007-2298, 2607-2898）；
  `rept s-align`（走 lastz，不涉及 pgi）不变；`sd cross` 正常检出 2 条
  跨基因组命中。
* 含 300 bp 分歧中段的库重复（真实块，非种子缺口）：e-align 仍检出两个
  完整拷贝（2006-3098 / 3606-4698），greedy 门控未误伤。

## 复核 23 验证

* 全量 1209 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 24（2026-08-05 后续轮次）

### 复核 21 修复的参数空间模糊验证（未发现新问题）

* 12 组随机基因组（3-6 拷贝、拷贝间隔 150-1800 bp 跨越 max_gap 边界、
  混合链向）`sd run`：全部拷贝均出现在至少一条 elementary 行中。
* 观察到的"某拷贝仅 non-core"是贪心 cover 的合法冗余消除：当某集合
  （含该拷贝同源代表）已覆盖全部命中时，冗余集合标 non-core，与 BISER
  文档语义一致。
* 复核脚本曾误报漏检，根因是分析脚本种子不一致（rng 与 rep 种子混淆），
  非代码问题。

## 复核 24 验证

* 全量 1209 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 25（2026-08-05 后续轮次）

### 复核 21 门控的边界行为确认（未发现新问题）

* 重读修改后的 greedy 循环控制流：门控失败 → 闭合当前链 + 以该种子起新
  链；off-band 忽略；跨组闭合——三种分支均正确，`saturating_add(k)` 对
  近 u32::MAX 的 last_a 安全（范围反转为空时视为可延伸）。
* 倒位对尺寸/间隙组合（L=1100-2500、gap=120-800）行为符合阈值设计：
  gap ≥ 200 → 两条干净命中；gap < 200 → 一条高身份合并块（如 L=1300、
  gap=120 → 2708 bp / 身份 0.98 ≥ 0.9），区域始终被标注，无检测丢失。
  L < 1000（如 800）被 SD 阈值正确过滤。

## 复核 25 验证

* 全量 1209 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 26（2026-08-05 后续轮次）

### 收尾（未发现新问题）

* `docs/align-pgi.md` 的 `--max-gap` 说明补齐复核 21 门控语义（双侧间隙
  ≥ 200 bp 时仅同源中段才桥接，近距离倒位对保持分离）。
* 近距离倒位基因组 `sd run` 端到端：两条拷贝均输出 elementary 行（各
  含 CORE 方向），配对不再丢失。
* 全量 1209 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 27（2026-08-05 后续轮次）

### 文档与深度语义复查（未发现新问题）

* `docs/rept.md` 全文重读：e-kmer 缓存命名/失效、e-align 身份定义、
  s-align 深度语义（50% 窗口基线 2、`--min-depth 4` = 至少 2 拷贝）、
  各参数表与实现一致。
* 4 拷贝串联重复的 `s-align` 深度：默认 min-depth 4 → 全区间 1-1200；
  min-depth 8 → 101-1200（起点边界深度不足被切除），语义正确。

## 复核 27 验证

* 全量 1209 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。
