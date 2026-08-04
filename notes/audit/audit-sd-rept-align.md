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

## 复核 28（2026-08-05 后续轮次）

### 代码重读与端到端确认（未发现新问题）

* `decompose` 片段合并循环逐分支核对：共享 run 内联、间隙容忍（≤50 bp）、
  超间隙退出后外层从退出点继续、尾随非共享碱基不并入片段、`end=last+1`
  与 `score` 计数均正确。
* 关键端到端复跑：invclose（4 行）、gclose0（12 行）、gr0（8 行）输出
  稳定。

## 复核 28 验证

* 全量 1209 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 29（2026-08-05 后续轮次）

### 特殊字符 contig 名复核（未发现新问题）

* 含空格 contig 名（`>chr one`）：`e-kmer`/`s-kmer`/`s-align` 均正确检出
  重复，输出键取首个空白 token（"chr"）——与 `fa size` 首字段约定一致
  （spanr 系既有行为）；`sd run` 空输出正确（300 bp 重复低于 SD 阈值）。
* 早期观察到的 e-kmer 空输出为库/基因组不匹配（测试种子不一致），非
  代码问题。

## 复核 29 验证

* 全量 1209 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 30（2026-08-05 后续轮次）

### 最终全面回归扫描（未发现新问题）

* 全部关键回归场景一次通过：近距离倒位对（2 条命中）、冒号 contig 名
  s-align（键完整）、多拷贝近距离（12 行）、`.pgi` 输入拒绝（友好错误）。
* 全量 1209 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。
* 连续九轮（复核 22–30）未发现新缺陷，审核已高度收敛。

## 复核 30 验证

* 全量 1209 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 31（2026-08-05 后续轮次）

### 修复的缺陷（1 处，功能）

**（功能）lastz self 模式用 basename 判断自比对，同名文件被交叉比对**：
`run_lastz` 的 self 跳过条件 `t_base != q_base` 只比 basename——目录中含
两个同名文件（如 `a/dup.fa`、`b/dup.fa`）时，`(a/dup.fa, b/dup.fa)` 会
以交叉比对方式运行（4 个 LAV 中 2 个为虚假交叉作业），对含共享序列的
基因组会产生错误命中。修复：self 模式跳过所有 `target_file !=
query_file` 的作业（每个文件只与其自身比对）。回归
`command_align_lastz_self_duplicate_basenames`（只产生 2 个 self LAV）。
`sd search --engine lastz` 与 `s-align` 同步受益（s-align 非 self，不受
影响）。

## 复核 31 验证

* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 32（2026-08-05 后续轮次）

### 复核 31 修复的回归确认（未发现新问题）

* gz 基因组解压后 self：target/query 列表指向同一批解压路径，路径相等
  判断成立（1 job、`--self`、2 命中）。
* 单文件 self：1 job、1 LAV，行为不变。
* MG1655 拆分基因组 lastz self：1 job、282 命中（引擎间差异，非回归）；
  `sd run` 正常。
* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 32 验证

* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 33（2026-08-05 后续轮次）

### release 模式回归确认（未发现新问题）

* 三个家族 release 构建全部通过：pgi 52、sd 13、pl 1、alignment 23 个
  lib 测试 + 44 个 CLI 集成测试（含复核 21/31 新增回归）。
* 历轮修复（含 greedy 门控与 self 路径判断）在 release 下行为一致。

## 复核 33 验证

* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 34（2026-08-05 后续轮次）

### 专项复核（未发现新问题）

* `e-kmer --keep-index` 缓存按 k 值隔离：`repeat.k17.*` 与
  `repeat.k21.*` 并存，k17 运行正确复用 k17 缓存。
* `trf --min-score` 边界（0 与 2000）：参数正确转发，串联重复均检出
  （TRF 自身评分语义）。

## 复核 34 验证

* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 35（2026-08-05 后续轮次）

### lastz 引擎端到端（未发现新问题）

* 混合基因组（近距离正链对 500 bp 间隔 + 近距离倒位对 700 bp 间隔）
  `sd run --engine lastz`：4 条拷贝全部检出并含 CORE 行（lastz 自身无
  greedy 链问题，不受复核 21/31 修复影响）。

## 复核 35 验证

* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 36（2026-08-05 后续轮次）

### 并行确定性（未发现新问题）

* `sd search --engine pgi` 在 `-p 1/2/8` 下输出逐字节一致；`sd run`
  多次运行输出一致（复核 10 的 `--parallel` 池重构未破坏确定性）。

## 复核 36 验证

* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 37（2026-08-05 后续轮次）

### rept 家族级端到端（未发现新问题）

* MG1655 上 5 个 rept 命令全部正常：e-kmer 48 区间 / 56,844 bp、
  e-align 89 / 60,251、s-kmer 170 / 128,444、s-align 1457 / 244,460、
  trf 84 / 18,768——其中 e-kmer 与 trf 与 docs 记录值逐字节一致。

## 复核 37 验证

* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 38（2026-08-05 后续轮次）

### lastz 引擎坐标精度（未发现新问题）

* 3 拷贝混合基因组（正链/倒位/正链，已知坐标）`sd run --engine lastz`：
  全部拷贝被 elementary 行覆盖，坐标与真实拷贝边界一致（修剪 ≤ 9 bp），
  CORE 标记为合法最小覆盖。

## 复核 38 验证

* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 39（2026-08-05 后续轮次）

### 跨染色体命中 cover 语义（未发现新问题）

* 跨染色体 SD 命中（q chrA ↔ t chrB 等）：`sd cover` 正确按"任一拷贝
  覆盖 hit 的 query 或 target 区间"判定，两组 elementary 均标 CORE。

## 复核 39 验证

* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 40（2026-08-05 后续轮次）

### 修复的缺陷（1 处，文档）

**（文档）`align pgi` 兄弟索引命名描述错误**：`sibling_pgi_path` 用
`set_extension("pgi")` 替换最终扩展名——`ref.fa` 与 `ref.fa.gz` 都映射到
`ref.pgi`（`.fa` 被替换而非保留）；docs/align-pgi.md 却写
"ref.fa.gz → ref.fa.pgi"。修复文档为"最终扩展名替换为 .pgi，ref.fa 与
ref.fa.gz 均映射到 ref.pgi"。同名共存时由 contig 校验拦截内容不匹配，
无静默错误。

## 复核 40 验证

* `.fa.gz` 输入 `--keep-index` 首次构建 / 二次复用 / 输出一致均验证通过。
* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 41（2026-08-05 后续轮次）

### align lastz 帮助与文档核对（未发现新问题）

* `align lastz --help` 的参数/默认值与 docs/align-lastz.md 一致；
  `--show-preset` 输出 preset 描述/参数/打分矩阵渲染正确。

## 复核 41 验证

* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 42（2026-08-05 后续轮次）

### 跨基因组多拷贝映射语义（记录项，非缺陷）

* 3 拷贝靶基因组 × 2 拷贝查询基因组：`align pgi` 原始 PSL 输出全部 6 对
  （各 ≥1000 bp、身份 ≥0.9）；chainnet 精修后每靶位点保留一条最优链
  （同等分时按排序取一），`sd cross` 输出 3 条。该行为是 UCSC
  chainnet 每靶位点取最佳链的标准语义（与 `pl chainnet` 共享），非本次
  改动引入；需要完整多拷贝映射的用户可直接用 search 阶段 PSL。
* 3 对原始命中全部通过 SD 过滤，无静默丢弃（过滤语义正确）。

## 复核 42 验证

* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 43（2026-08-05 后续轮次）

### 修复与核对（1 处小修 + 验证）

**（静默错误）`sd align` 跳过非 2 组件的 MAF 块时无提示**：`maf_block_to_paf`
对 <2 / >2 组件的块返回 None（注释称"caller logs warning"），但
`chainnet_to_paf` 未记日志。补 `log::warn!`（chainnet 输出恒为 2 组件，
路径为防御性提示）。

* `--self` 显式同文件：21 块（精确自比对被丢弃）vs 无 --self 的 22 块
  （含自比对对角线），语义正确；`--self` 传不同文件报友好错误。

## 复核 43 验证

* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 44（2026-08-05 后续轮次）

### .2bit 输入全路径（未发现新问题）

* `align pgi` 与 `sd search --engine pgi` 接受 `.2bit` 输入，输出与 FASTA
  路径逐字节一致（同基因组 4 块 / 4 命中）。

## 复核 44 验证

* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 45（2026-08-05 后续轮次）

### lastz preset 变体（未发现新问题）

* `sd run --engine lastz` 在 set01/set03 下均检出两个重复家族（各 4 行），
  坐标仅边界修剪差异——preset 转发与结果一致性正确。

## 复核 45 验证

* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 46（2026-08-05 后续轮次）

### cluster 链向提取（未发现新问题）

* 负链区间提取为 RC 序列、正链区间提取为正向序列，逐碱基正确（构造
  q[1000,1100)- ↔ t[1100,1200)+ 的倒位对验证通过）。

## 复核 46 验证

* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 47（2026-08-05 后续轮次）

### 临时文件清理（未发现新问题）

* `sd run` 完成后 `/tmp/pgr_sd_*` 临时目录数为 0（PipelineCtx 随作用域
  正确清理），输出目录仅含 `out.elem.bed`，无文件泄漏。

## 复核 47 验证

* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 48（2026-08-05 后续轮次）

### 部分同源 decompose（未发现新问题）

* 两序列共享 200 bp 核心 + 各自独有侧翼：decompose 正确检出共享区
  （k-mer 边界修剪后 191 bp），两片段归同一 elementary set。

## 复核 48 验证

* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 49（2026-08-05 后续轮次）

### 混合链向家族（未发现新问题）

* 3 正链 + 2 倒位拷贝的同一重复家族 `sd run`：5 条拷贝全部出现在
  elementary 行中（9 行，CORE/非-core 为合法最小覆盖分配），正/倒位
  方向均正确处理。

## 复核 49 验证

* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 50（2026-08-05 后续轮次）

### 记录项（未改，范围外输入）

* `rept e-align` 传入 `.2bit` 基因组：在 `has_soft_mask` 的 FASTA 读取器
  处报 "stream did not contain valid UTF-8"（二进制被当文本读）。文档仅
  承诺 FASTA（.fa/.fa.gz）；`sd run` 对 2bit 也在 chainnet 的 `fa size`
  步骤报错（复核 11 已记录）。两处均为有错误提示的非静默失败，属范围外
  输入，记录不修。

## 复核 50 验证

* 全量 1210 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 51（2026-08-05 后续轮次）

### 修复的缺陷（1 处，功能）

**（功能）`ref.fa` 与 `ref.fa.gz` 共享兄弟索引，内容不同时静默复用错误
索引**：`sibling_pgi_path` 的 `set_extension("")` + `set_extension("pgi")`
链把 `.fa` 替换掉，`ref.fa.gz` 与 `ref.fa` 都映射到 `ref.pgi`。当两文件
同名同长但序列不同时，contig 校验（只比名字/长度）无法拦截，第二次运行
静默复用第一次的索引（实测 0 块输出）。复核 40 曾把文档改成"共享"描述，
实为路径构造 bug。修复：`.gz` 输入去掉 `.gz` 后**追加** `.pgi`（`ref.fa.gz`
→ `ref.fa.pgi`），与 `ref.fa` → `ref.pgi` 分离；文档恢复原意。回归
`command_align_pgi_gz_sibling_index_distinct`。

## 复核 51 验证

* `.fa.gz` 首次构建 `ref.fa.pgi`、二次复用、输出一致；与 `ref.fa` 的
  `ref.pgi` 互不干扰。
* 全量 1211 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 52（2026-08-05 后续轮次）

### 复核 51 修复的回归确认（未发现新问题）

* `sd search` 与 `rept e-align` 的 gz 基因组路径不受兄弟索引命名修复
  影响（分别 2 命中 / 正确区间）。

## 复核 52 验证

* 全量 1211 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 53（2026-08-05 后续轮次）

### gz 基因组端到端（未发现新问题）

* 复核 51 修复后 `sd run tworep2.fa.gz` 输出与明文路径逐字节一致（4 行
  elementary、全部 CORE）。此前一次报错为测试文件未创建，非代码问题。

## 复核 53 验证

* 全量 1211 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 54（2026-08-05 后续轮次）

### e-align 确定性与身份边界（未发现新问题）

* e-align 多次运行输出逐字节一致；`--min-identity 0.999` 与 0.5 下精确
  拷贝（身份 ~1.0）均正确检出。

## 复核 54 验证

* 全量 1211 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 55（2026-08-05 后续轮次）

### min-len 边界（未发现新问题）

* 1183 bp 倒位块：`--min-len 1183` 检出 2 条、`--min-len 1184` 为 0
  ——过滤为精确的 `block_len >= min_len`。

## 复核 55 验证

* 全量 1211 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 56（2026-08-05 后续轮次）

### lastz query-depth 语义（未发现新问题）

* `--query-depth 2` 与 50 在 2 拷贝重复基因组上均检出 4 条（深度阈值只
  在覆盖超过阈值后截断搜索，低拷贝下无影响）——语义正确。

## 复核 56 验证

* 全量 1211 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 57（2026-08-05 后续轮次）

### 贪心 cover 对抗性验证（未发现新问题）

* 5 个重叠 elementary 集合 × 3 条命中：贪心选出 1 个 CORE 集合覆盖全部
  命中（集合 4 [150,250) 同时覆盖三条的 query 或 target 区间），输出为
  合法最小覆盖。

## 复核 57 验证

* 全量 1211 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 58（2026-08-05 后续轮次）

### tube 并行确定性（未发现新问题）

* tube workflow 多次运行与 `-p 1/4/8` 下输出逐字节一致（并行 tube 扩展
  的 containment 去重确定性）。

## 复核 58 验证

* 全量 1211 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 59（2026-08-05 后续轮次）

### .pgi 无序列路径（未发现新问题）

* `.pgi` 输入无扩展序列：self 输出 4 块（无打分）、pair 输出 5 块（含
  精确自比对对角线）——块输出语义正确。

## 复核 59 验证

* 全量 1211 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 60（2026-08-05 后续轮次）

### 修复的缺陷（1 处，功能/确定性）

**（功能）`sd cluster` 的 cluster 编号依赖 HashMap 迭代顺序，`sd run`
的 set_id 编号跨运行不稳定**：`cluster_paf` 按连通分量分组后直接迭代
`HashMap`（进程内随机种子），cluster_N 的编号与文件名对应关系每次运行
不同；`sd run` 虽按数值排序 cluster 文件，但编号本身随机，导致同一基因组
多次运行输出 set_id/行序互换（两家族时 r1/r2 互换 set 1/2，实测 5 次运行
2 种输出）。修复：按每个分组的首个区间（chrom, start）排序后再编号。
回归 `command_sd_run_output_deterministic_across_runs`；10 家族基因组
多次运行逐字节一致。

## 复核 60 验证

* `sd run` 6 次运行 1 种哈希；10 家族 2 次运行 diff 为空。
* 全量 1212 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 61（2026-08-05 后续轮次）

### HashMap 迭代顺序全面审计（未发现新问题）

* 三个家族全部 HashMap 用法逐一核对：pl/repeat 与 trf 的 name/safe_map
  仅查找；cover 的 by_set/coverage 经 set_order（Vec）迭代；decompose 的
  index/kmer_frags/pair_count 顺序无关（输出按 frags Vec 序）；
  cluster 的 by_root 已在复核 60 排序。仅复核 60 的 cluster 编号曾依赖
  HashMap 迭代序，其余均确定。
* `sd cluster` 独立命令、`sd decompose`、`sd cover` 多次运行输出逐字节
  一致。

## 复核 61 验证

* 全量 1212 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 62（2026-08-05 后续轮次）

### 修复的缺陷（1 处，功能）

**（功能）FASTA 原地修改后兄弟索引被静默复用**：`resolve_side` 复用同名
兄弟 `.pgi` 时只校验 contig 名/长度；同名单长但序列不同的 FASTA 会静默
复用旧索引（k-mer 来自旧序列），对齐结果错误。修复：新增 mtime 校验
（输入比索引新则重建，与 e-kmer 缓存同一约定）。回归
`command_align_pgi_stale_sibling_index_rebuilt`。

## 复核 62 验证

* 修改后重建（"built reference index"）、未修改复用（"reusing"）；
  `.fa.gz` 路径不受影响。
* 全量 1213 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 63（2026-08-05 后续轮次）

### 复核 62 收尾（未发现新问题）

* `docs/align-pgi.md` 的兄弟索引说明补齐 mtime 失效约定。
* `sd search`（4 命中）与 `rept e-align`（正确区间）不受 mtime 修复影响。

### 记录项（未改）

* `open_indexed` 的 `.loc` 索引按存在性复用（`force_update=false`）：基因组
  修改后 `.loc` 字节偏移可能过期。属既有设计（与 `fa range` 等共享），
  且修改中间基因组本身就是用户错误，记录不修。

## 复核 63 验证

* 全量 1213 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 64（2026-08-05 后续轮次）

### lastz self 真实数据（未发现新问题）

* MG1655 拆分基因组 `align lastz --self`：1 job、`--self`、1 LAV（复核 31
  路径相等修复在真实数据上正确）。

## 复核 64 验证

* 全量 1213 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 65（2026-08-05 后续轮次）

### merge-gap 与工作流交互（未发现新问题）

* greedy 在 merge-gap 0/5000 下远距家族结果一致；tube workflow 完全忽略
  merge-gap（独立链算法，参数互不干扰）。

## 复核 65 验证

* 全量 1214 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 66（2026-08-05 后续轮次）

### 多染色体端到端（未发现新问题）

* 同一重复家族跨 chr1（2 拷贝）与 chr2（1 拷贝）：`sd run` 输出 4 行
  elementary、同一 set、全部 CORE——跨染色体分组正确。

## 复核 66 验证

* 全量 1217 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 67（2026-08-05 后续轮次）

### lastz 多 contig 错误路径（未发现新问题）

* gz 多 contig 基因组 `sd search --engine lastz`：报 "contains more than
  one sequence; consider using the 'multiple' action"（首轮 lastz stderr
  捕获修复生效），提示用户拆分。

## 复核 67 验证

* 全量 1218 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 68（2026-08-05 后续轮次）

### 环境收尾（未发现本族新问题）

* 并行会话在 `fa/count.rs`、`fa/size.rs` 的未格式化改动曾使 `cargo fmt
  --check` 门禁失败；执行标准 `cargo fmt` 统一格式（机械性，不影响其
  语义改动），门禁恢复。

## 复核 68 验证

* 全量 1218 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 69（2026-08-05 后续轮次）

### MG1655 lastz 端到端（未发现新问题）

* MG1655 拆分基因组 `sd run --engine lastz`：25 s（debug）完成，120 行
  elementary、116 行 CORE——lastz 引擎真实数据全流程正常。

## 复核 69 验证

* 全量 1218 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 70（2026-08-05 后续轮次）

### 环境收尾（未发现本族新问题）

* 并行会话继续修改多个 fa/ 文件（dedup/filter/mask/order/replace/some/
  to_2bit）未格式化；再次 `cargo fmt` 恢复门禁（机械性）。期间出现过
  一次瞬时编译错误，随后自行恢复（会话在途编辑），当前构建/测试正常。

## 复核 70 验证

* 全量 1219 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 71（2026-08-05 后续轮次）

### sd cross 双引擎（未发现新问题）

* lastz 引擎 `sd cross` 多次运行输出逐字节一致（3 条跨基因组命中，
  与 pgi 引擎一致）。

## 复核 71 验证

* 全量 1219 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 72（2026-08-05 后续轮次）

### 环境收尾（未发现本族新问题）

* 并行会话继续修改 `libs/fmt/fa.rs` 未格式化；再次 `cargo fmt` 恢复门禁。

## 复核 72 验证

* 全量 1220 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 73（2026-08-05 后续轮次）

### 旧命名遗留索引兼容（未发现新问题）

* 复核 51 修复前遗留的 `ref.pgi`（旧命名）被新代码忽略，`.fa.gz` 重建
  `ref.fa.pgi`，两者共存无冲突（遗留文件无害）。

## 复核 73 验证

* 全量 1220 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 74（2026-08-05 后续轮次）

### trf max-period（记录项，外部工具限制）

* `--max-period` 正确转发：默认 2000 排除 2500 bp 周期串联（`{}`）。
* TRF 自身对"完美 2500 bp 周期串联 + max-period ≥ 2600"输入 SIGSEGV
  （直调 trf 复现），pgr 将信号错误友好传播（无 panic）。属 TRF 限制，
  记录不修。

## 复核 74 验证

* 全量 1220 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 75（2026-08-05 后续轮次）

### sd run min-len 筛选（未发现新问题）

* `--min-len 1100` 只保留 1188 bp 家族（1087 bp 家族被滤除）；
  `--min-len 1200` 全部滤除——过滤沿管线正确传导。

## 复核 75 验证

* 全量 1220 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 76（2026-08-05 后续轮次）

### 修复的缺陷（1 处，内存安全）

**（内存）中段同源检查的 DP band 直接取用户 `--band`，极端组合可 OOM**：
`middle_is_homologous_range` 的 `align_banded_local` 使用用户 band 原值；
`--band 10000` + 大中段（≤ 50 kb）时 DP 分配 ~13 GB。修复：DP band 上限
256（检查只需探测预期对角偏移内同源性，超限保守不合并=碎片化而非丢失）。
回归验证：`--band 10000 -s 50000` 3.4 s 完成（此前 OOM 风险）。

### 环境说明

* 并行会话在途改动 `fa/split.rs`（新增参数）触发 clippy too-many-arguments；
  不触碰其文件，本族 lib clippy 干净、544 个 lib 测试全通过。

## 复核 76 验证

* 全量 lib 544 测试通过；`cargo fmt --check` 干净。

## 复核 77（2026-08-05 后续轮次）

### DP band 全面审计（未发现新问题）

* `align_banded_local` 全部调用方 band 上限核对：中段检查（≤256，复核
  76 修复）、窗口扩展（≤128，diag_span+32 封顶）、wave tube 扩展
  （桶内对角带 ≤~128）、测试固定值——无未封顶路径。

## 复核 77 验证

* 全量 lib 544 测试通过；`cargo fmt --check` 干净。

## 复核 78（2026-08-05 后续轮次）

### k 值边界（未发现新问题）

* `-k 64`（u128 键上限）端到端正常（2 块、坐标正确）；`-k 65` 报
  "k must be in 1..=64"。

## 复核 78 验证

* 全量 lib 544 测试通过；`cargo fmt --check` 干净。

## 复核 79（2026-08-05 后续轮次）

### .pgi + seqs 扩展路径（未发现新问题）

* `.pgi` 输入 + `--ref-seq`/`--query-seq` 自比对：4 块带真实打分的扩展
  PSL（1197/1096 匹配），坐标正确。

## 复核 79 验证

* 全量 lib 544 测试通过；`cargo fmt --check` 干净。

## 复核 80（2026-08-05 后续轮次）

### stdout 输出路径（未发现新问题）

* `rept trf` / `s-kmer`（runlist JSON）与 `sd decompose`（BED）默认
  stdout 输出均正常。

## 复核 80 验证

* 全量 lib 544 测试通过；`cargo fmt --check` 干净。

## 复核 81（2026-08-05 后续轮次）

### s-align 非重叠窗口深度（未发现新问题）

* `--step 200`（1x 基线）默认 `--min-depth 4` → 空（2 拷贝 × 1 窗口 =
  深度 2 < 4）；`--min-depth 2` → 检出两个拷贝——深度语义随窗口覆盖
  正确调整。

## 复核 81 验证

* 全量 lib 544 测试通过；`cargo fmt --check` 干净。

## 复核 82（2026-08-05 后续轮次）

### 修复的缺陷（1 处，数据安全）

**（数据安全）sd 命令 `-o` 指向输入文件时静默覆盖输入**：`sd search
g.fa -o g.fa` 等会把输入 FASTA/PAF/BED 覆盖为变换后的输出（exit 0、
无提示）——与 rg/runlist 会话在复核 3/4 修复的同类缺陷。修复：`sd
search`/`align`/`cover`/`decompose`/`cross` 均加 `ensure_outfile_distinct`
检查。回归 `command_sd_output_same_as_input_rejected`。

## 复核 82 验证

* `sd search -o g.fa` / `sd decompose -o g.fa` 报 "also an input file"，
  输入保持完好。
* 全量 1223 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 83（2026-08-05 后续轮次）

### 修复的缺陷（1 处，数据安全）

**（数据安全）rept 与 align pgi 同样存在 `-o` 覆盖输入**：`rept s-kmer
g.fa -o g.fa`、`rept trf`、`align pgi g.fa -o g.fa` 等把输入 FASTA 覆盖
为 runlist JSON/PSL（exit 0、无提示）。修复：rept 五个子命令（含库输入）
与 `align pgi`（含 --ref-seq/--query-seq）均加 `ensure_outfile_distinct`。
回归 `command_rept_output_same_as_input_rejected`。

## 复核 83 验证

* `rept s-kmer`/`trf`/`align pgi` 的 `-o` 指向输入均报 "also an input
  file"，输入完好。
* 全量 1224 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 84（2026-08-05 后续轮次）

### 复核 82/83 修复回归确认（未发现新问题）

* `sd run` 与 `rept e-kmer` 全流程在 `-o` 碰撞保护后输出不变（4 行 /
  正确区间）——保护只拦截同名输入，不影响正常管线。

## 复核 84 验证

* 全量 1224 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 85（2026-08-05 后续轮次）

### align lastz outdir 碰撞（未发现新问题）

* `-o` 指向输入文件（outdir 为文件路径）：报 "File exists"，输入不被
  覆盖（create_dir_all 失败即停）。

## 复核 85 验证

* 全量 1224 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 86（2026-08-05 后续轮次）

### lastz 引擎确定性（未发现新问题）

* 10 家族基因组 `sd run --engine lastz` 两次运行输出逐字节一致（20 行）—
  复核 60 的 cluster 编号排序修复对双引擎均生效。

## 复核 86 验证

* 全量 1224 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 87（2026-08-05 后续轮次）

### s-align 分片（未发现新问题）

* `--chunk-records 1` 与 5 均输出相同正确区间（2001-2300, 2601-2900）—
  窗口分片不影响结果。

## 复核 87 验证

* 全量 1224 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 88（2026-08-05 后续轮次）

### gz 单 contig lastz 端到端（未发现新问题）

* `sd run tworep2.fa.gz --engine lastz`：4 行 elementary、坐标正确（解压
  路径 + 单序列约束正常）。

## 复核 88 验证

* 全量 1224 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 89（2026-08-05 后续轮次）

### 终极确定性（未发现新问题）

* `sd run` 8 次运行 1 种哈希；10 家族 2 次运行 diff 为空——历轮确定性
  修复（复核 3/10/60）整体成立。

## 复核 89 验证

* 全量 1224 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 90（2026-08-05 后续轮次）

### self + tube 组合（未发现新问题）

* `--self --workflow tube`：4 块扩展 PSL（2 家族 × 2 互惠），坐标与打分
  正确。

## 复核 90 验证

* 全量 1224 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 91（2026-08-05 后续轮次）

### min-len 边界重复（未发现新问题）

* 精确 1000 bp 重复经扩展修剪为 996 bp 块 → 默认 min-len 1000 滤除；
  `--min-len 990` 检出。属块长度语义（按比对块而非源长度），与复核 55
  的边界行为一致。

## 复核 91 验证

* 全量 1224 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 92（2026-08-05 后续轮次）

### 修复的缺陷（1 处，功能）

**（功能）tube 工作流在显式同文件对（非 --self）时把家族交叉命中当
"重复"丢弃**：精确自比对巨块（全基因组对角线）在 `dedupe_contained`
中把坐标上包含于其内的拷贝对块（两轴 ≥95% 包含）误判为重复并丢弃，
显式同文件对只输出 1 块自比对。修复：dedupe 增加跨度相近约束（前块
跨度 ≤ 后块 4 倍才判重复）。回归
`dedupe_keeps_small_block_inside_large_one`；显式同文件对 tube 输出
5 块（自比对 + 4 家族命中），tube self 模式 4 块不变。

## 复核 92 验证

* 全量 1226 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 93（2026-08-05 后续轮次）

### 复核 92 修复回归确认（未发现新问题）

* tube 自比对多场景不变：5 拷贝 20 块、近距离倒位对 2 块、gclose0 22
  块——dedupe 跨度约束不影响正常 tube 输出去重。

## 复核 93 验证

* 全量 1226 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 94（2026-08-05 后续轮次）

### 2bit + tube 组合（未发现新问题）

* `.2bit` 输入 `--workflow tube`：4 块；`sd search` 2bit 4 命中。

## 复核 94 验证

* 全量 1226 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 95（2026-08-05 后续轮次）

### 巨管/整染色体串联压力（未发现新问题）

* 40 拷贝完美串联（20 kb）：默认 freq 10 正确过滤重复 k-mer（无种子→
  空输出）；`-f 100` 检出 100 块——FastGA 频率过滤语义正确，无崩溃。

## 复核 95 验证

* 全量 1226 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 96（2026-08-05 后续轮次）

### min-span 边界（未发现新问题）

* 300 bp 库链：`--min-span 50` 检出 2 块、500/5000 滤除——过滤为精确的
  每轴种子跨度 ≥ min-span。

## 复核 96 验证

* 全量 1226 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 97（2026-08-05 后续轮次）

### tube 覆盖度边界（未发现新问题）

* tube CHAIN_MIN 85 边界精确：80 bp 重复不成管（0 块）、85 bp 起成管
  （2 块）——FastGA 管覆盖度语义正确。

## 复核 97 验证

* 全量 1226 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 98（2026-08-05 后续轮次）

### contig 起点重复（未发现新问题）

* 重复位于 contig 起点（位置 0）：e-kmer 正确检出 1-300 与 801-1100。

## 复核 98 验证

* 全量 1226 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 99（2026-08-05 后续轮次）

### 低身份（92% 分歧）拷贝双引擎检出差异（引擎灵敏度差异，非缺陷）

构造 1200 bp、8% 错配（92% identity）的两个拷贝对后：

* `sd search --engine lastz` → 2 命中（1128 匹配 / 72 错配，块长 1200 bp，
  identity 94.0%，覆盖全长）；
* `sd search --engine pgi` → 0 命中。pgi 引擎自身写到 2 块 PSL（753 匹配 /
  41 错配，块长 794 bp，identity 94.8%），但块长 794 < 默认 `--min-len 1000`
  被 SD 过滤器丢弃。

根因验证：`pgr align pgi -k` 递减时块长单调恢复（k=40→794、k=35→789、
k=30→899、k=25→927、k=20→1050 bp），证明 pgi 链扩展无缺陷——覆盖度受
精确 k-mer seed（默认 k=40 + syncmer 8/5）限制，8% 分歧下两侧 ~200 bp 无
共享 seed，只能锚定中间子区段。lastz 用 12-mer seed + 扩展覆盖全长。

结论：**引擎灵敏度差异**（pgi 精确 seed 对近 SD 阈值的低身份拷贝锚定不足），
非 pgi 特异性漏检，也无崩溃/越界问题。T2T SD 标准（≥1 kb、≥90% identity）
下 92% 分歧属于边缘场景；MG1655 等真实数据双引擎一致（229/232/282 条）。

文档修复：`docs/sd.md` search 节补充说明——pgi 引擎用精确 k-mer seed，
近 90-93% identity 的拷贝可能只锚定子块，低于 `--min-len` 被滤除；提示可
降低 `--min-len` 或改用 `--engine lastz` 获得最大灵敏度。

## 复核 99 验证

* 全量 1226 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 100（2026-08-05 后续轮次）

### rept 命令族纵深验证（未发现新问题）

本轮转向历轮覆盖较少的 `rept` 命令族（e-kmer/s-kmer/e-align/s-align/trf）
与 `sd` 空输入链路：

* **e-kmer `--keep-index` 缓存**：首次运行建 `lib.fa.gz.repeat.k17.ktab`
  + 隐藏 part 文件 + `.complete` 标记；二次运行日志显示 "reused repeat
  table"，两次输出逐字节一致。缓存 mtime 失效与复用逻辑正确。
* **trf 0 参数边界**：`--trf-match 0` / `--pm 0` 被 TRF 拒绝、错误清晰
  传播；`--pi 0` / `--max-period 0` 被 TRF 接受并空输出（安全降级，无
  panic/死循环）。
* **s-align `--step 0`**：`fa window` 自带正数校验，错误经 cmd_lib 清晰
  传播（"must be positive"）。
* **Profex 输出解析**：`Read N:` 头行安全跳过；最后一个未闭合 run 无换行
  符，`BufRead::lines()` 正确读取；e-kmer 以染色体长度闭合尾部 run。
* **TRF 版本兼容**：本机 TRF 4.09 `-ngs` 输出 17 字段（含末尾 `. .`），
  `parse_trf_output` 的 ≥15 字段门槛兼容；`@chr1` 头行（1 字段）跳过。
* **`fa size` / `chr.sizes`**：noodles 的 `record.name()` 截断到首个空白，
  名字含空格时 `split_whitespace` 解析仍安全。

### sd 空输入全链路（未发现新问题）

100 kb 随机序列（无 SD）上：

* `sd search` → 0 块空 PSL（exit 0）；
* `sd align` 空 PSL → 空 PAF（exit 0）；
* `sd cluster` 空 PAF → 空目录（exit 0）；
* `sd run` 全链路 → 空 `out.elem.bed`，exit 0（decompose 循环跳过、
  cover 空输入优雅）。
* `sd search --engine lastz` 多 contig 输入：lastz 报 "contains more than
  one sequence"，pgr 清晰传播为错误（非 panic/静默空）。

### 文档与实现一致性（未发现新问题）

`docs/rept.md` 五个命令（e-kmer/s-kmer/e-align/s-align/trf）的参数表与
clap 定义逐项一致（名称、默认值、short 选项）；`sd cover` 的集合覆盖
确定性良好（`set_order` 保序 + `max_by_key` 平局选择确定）。

## 复核 100 验证

* 全量 1227 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 101（2026-08-05 后续轮次）

### 解析器与参数边界纵深验证（未发现新问题）

* `parse_paf`：12 字段强校验 + 错误带行上下文，跳过空行/`#` 行。
* `sd cluster` 越界 PAF 坐标（end > contig length）：`fetch_range_seq`
  → noodles `slice` 返回 None → "slice error for [...]" 报错，非 panic。
* `Range::from_str` 手写扫描器 `parse_i32` 有溢出保护（超 i32 数字返回
  None → 行不匹配 → `usable_range` 过滤）；`from_str_regex` 仅为
  `#[cfg(test)]` 测试 oracle，非生产路径。JSON runlist 超大坐标报
  "Number format error: out of range"。
* `align lastz --parallel 0`：rayon 安全接受（不崩溃）。
* lastz 全部 7 个 preset（set01..set07）均正常执行；`--query-depth`
  0/1/1000000 边界安全且输出一致。
* `run_lastz` self 模式路径相等判断、输出名 `create_new` 原子保留 +
  序号、job 失败聚合报错均正确。

## 复核 101 验证

* 全量 1227 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 102（2026-08-05 后续轮次）

### 模糊与空输入测试（未发现新问题）

* 随机二进制喂 `sd search` / `rept trf` / `sd cluster` / `align pgi`：
  全部报错退出（"stream did not contain valid UTF-8"），无 panic。
* 空 FASTA（0 字节）：`sd search` / `align pgi` 报 "index has no
  contigs"，错误清晰。
* 只有头的 FASTA（`>chr1` 无序列）：`sd search` / `align pgi` 输出空
  PSL（exit 0）；`rept s-kmer` 触发 FastK SIGSEGV——外部工具对空序列
  崩溃，cmd_lib 捕获并报 "terminated by signal: 11"，pgr 自身无 panic，
  与 TRF 2500bp 周期 SIGSEGV 同类，记录不修。

### `pl chainnet` 实现审查（未发现新问题）

* native psl-chain-net-axt-maf 管道无外部 kent-tools 依赖；空 PSL 输入
  时各阶段（psl chain → anti-repeat → sort → pre-net → net → axt →
  maf）与 sd align 的空链路验证（复核 100）一致。

## 复核 102 验证

* 全量 1227 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 103（2026-08-05 后续轮次）

### rept 端到端与确定性验证（未发现新问题）

* **确定性**：`rept s-align` / `rept trf` / `rept s-kmer` 各两次运行输出
  逐字节一致。
* **e-align → fa mask 端到端**：g1（无 IS 元件）空 runlist；mg1655 检出
  IS 元件区间（15393-16728、19800-20557 等），`fa mask` 正常输出。
* **trf 多 contig**：构造 `ctgA`（ACGT 串联）/`ctgB`（TGCA 串联）双 contig
  基因组，输出恰含两个键且各 1 个区间；`fa mask` 接受 trf runlist。
* **s-align 坐标合理性**：g1（12614 bp 高重复）覆盖 10495 bp（83%），
  与 `sd search` 12 行命中一致，无坐标偏移证据。

## 复核 103 验证

* 全量 1229 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 104（2026-08-05 后续轮次）

### 缺陷 29：`.loc` 索引陈旧时静默使用（open_indexed 只查存在性）

`loc::open_indexed` 仅在 `.loc` 文件不存在（或 `force_update`）时重建，
从不校验索引与 FASTA 的新鲜度。实测复现：

1. 1200 bp FASTA 建 `g.fa.loc`（记录 size=1207）；
2. FASTA 改为 1500 bp 后 `sd cluster` 仍用陈旧索引 → `slice error`
   （长度变化可报错，但**同长度内容修改会静默提取错误序列**——更危险）；
3. 修复前 `fa range chr1:301-1500` 同样报 slice error。

修复：`open_indexed` 增加 mtime 新鲜度校验——`.loc` 的 mtime 早于 FASTA
 时自动重建（`loc_is_fresh`，mtimes 不可用时保持旧行为）。`fa range` /
 `sd cluster` / `fas check` / `get_seq_loc` 四个调用方同步受益。

回归 `stale_loc_index_is_rebuilt`：同长度内容修改（ACGT→TGCA）+ `.loc`
 mtime 调旧后，重建索引并返回新序列。

## 复核 104 验证

* 全量 1230 测试通过（含新回归测试）；`cargo fmt --check` 与
  `cargo clippy --all-targets -- -D warnings` 干净。

## 复核 105（2026-08-05 后续轮次）

### `.loc` 修复端到端确认 + 其余索引陈旧覆盖盘点（未发现新问题）

* 实机复现最危险场景：ACGT 内容建 `.loc` → 同长度替换为 TGCA → `.loc`
  mtime 调旧 10 秒 → `fa range chr1:1-12` 返回 TGCA（自动重建，未返回
  陈旧 ACGT）。`sd cluster` 同样正确重建并输出。
* 其余索引/缓存盘点：
  - `.pgi` 兄弟索引：复核 62 已有 mtime 失效重建；
  - FastK repeat 表缓存：`cache_is_fresh` 校验 `.ktab`+`.complete` mtime；
  - `.paf.idx`：仅用户显式传入（`infile` 以 `.paf.idx` 结尾）时加载，
    无自动兄弟索引，新鲜度由用户负责，非静默自动使用；
  - `.2bit`：`fa to-2bit` 显式生成，无自动陈旧路径。
* 结论：`.loc` 是唯一遗漏自动新鲜度校验的伴随索引，复核 104 已修复，
  其余索引无同类问题。

## 复核 105 验证

* 全量 1230 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 106（2026-08-05 后续轮次）

### sd decompose / cover 字段语义验证（未发现新问题）

* 随机 300 bp 家族 3 拷贝（chr1+、chr1+、chr2-）decompose：
  - 三条全部 set_id=1（正确合并为同一 elementary SD 家族）；
  - `-` 链坐标反向投影正确（chr2 头 100-400 → 输出 109-400）；
  - length=291（300 bp 序列尾部 9 bp 无 10-mer 滑窗种子，k-mer 窗口
    边界语义正确）。
* cover：set_id=1 覆盖两个 PAF hit（同染色体串联 + 跨染色体）→ 三行
  全部 CORE，正确。
* 记录项（低风险，不修）：纯四联体重复（如 ACGT）只有 4 种不同 10-mer，
  低于 `MIN_SHARED_KMERS=5` 防过度分组阈值，同源片段不会合并为同一
  set_id——极端低复杂度序列，非 SD 场景，行为符合设计意图。
* 注：手工构造 cluster FASTA 时 `-` 链序列必须为 revcomp 后的家族方向
  序列（`sd cluster` 输出即如此），否则与 `+` 链拷贝 k-mer 不同源。

## 复核 106 验证

* 全量 1231 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 107（2026-08-05 后续轮次）

### sd cross 双引擎端到端（未发现新问题）

* `sd cross crossA.fa crossB.fa` pgi/lastz 各输出 3 块 PAF（约 1200 bp、
  identity ~96%），坐标合理，PAF 带 cg:Z/gi:Z/bi:Z/ms:Z 标签。
* 一个表面差异：A 拷贝 2（7200-8400）被 pgi 匹配到 B 拷贝 2（8200-9400）、
  被 lastz 匹配到 B 拷贝 1（4000-5200）。查证 B 拷贝 1 与拷贝 2 逐碱基
  100% 相同（0 错配）——这是对称重复造成的合法链化歧义，两引擎均正确，
  非缺陷。真实 SD 流程的 cluster 阶段会合并此类重复。

## 复核 107 验证

* 全量 1231 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 108（2026-08-05 后续轮次）

### align lastz 目录输入（未发现新问题）

* 2 targets × 1 query 目录：2 jobs 并行执行，输出 `[t1]vs[q1].lav` /
  `[t2]vs[q1].lav`，LAV 内容为标准 lastz 格式（scoring matrix 头），
  exit 0。

## 复核 109（2026-08-05 后续轮次）

### rept e-align workflow tube vs greedy（未发现新问题）

* mg1655 + tncentral：tube 52 区间（62137 bp）vs greedy 89 区间
  （60251 bp），tube 非 greedy 子集（3713 bp 仅 tube 覆盖）——FastGA
  管与 greedy 链的块边界/覆盖语义差异，属预期行为（`--workflow` 为
  文档化选项），两 workflow 均正常执行。

## 复核 110（2026-08-05 后续轮次）

### align pgi --keep-index 兄弟索引 mtime 失效（未发现新问题）

* plain `.fa` 输入生成 `ref.pgi` 兄弟索引；未修改时日志 "reusing
  reference index"（复用）；修改 ref.fa 后日志 "built reference index"
  （自动重建）——复核 62 的 mtime 失效重建实测通过。

## 复核 108–110 验证

* 全量 1232 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 111（2026-08-05 后续轮次）

### 缺陷 30：`.pgi` 单输入自比对 + 仅 `--ref-seq` 报错

`pgr align pgi ref.pgi --ref-seq ref.fa`（self 模式，仅提供参考侧扩展序列）
报错 "extension sequences are needed for both sides"：self 模式下 query 侧
复用的是 `.pgi` 输入的 `seqs=None`，`--ref-seq` 只喂了 ref 一侧，两侧
空/非空不一致触发 bail。用户被迫额外传相同的 `--query-seq` 才能工作。

复现与根因：
* `align pgi ref.pgi --ref-seq ref.fa` → 报错；
* `align pgi ref.pgi --ref-seq ref.fa --query-seq ref.fa` → 正常（0 块）；
* `align pgi ref.fa`（FASTA 直接）→ 正常（0 块）。

修复：`resolve_seqs` 后，self 模式下任一侧扩展序列为空时复用另一侧
（`.pgi` 单输入 + 仅 `--ref-seq` 或仅 `--query-seq` 的自然用法，两方向
对称）。验证：
* 随机序列（0 块）：仅 `--ref-seq`、仅 `--query-seq`、双侧、FASTA 直接
  输入四者输出逐字节一致；
* g1 高重复基因组（17 块）`.pgi + --ref-seq` 与 FASTA 直接输出逐字节一致；
* genome 输入 + `--query-seq` 的冲突报错路径不变（resolve_seqs 前置检查）。

回归 `command_align_pgi_single_ref_seq_on_self_pgi`：两拷贝 FASTA 的
`.pgi + --ref-seq` 与 `.pgi + --query-seq` 输出均与直接自比对逐字节相等。

## 复核 111 验证

* 全量 1234 测试通过（含新回归测试）；`cargo fmt --check` 与
  `cargo clippy --all-targets -- -D warnings` 干净。

## 复核 112（2026-08-05 后续轮次）

### 缺陷 30 对称情况补齐

* `.pgi` 单输入 + 仅 `--query-seq` 同样报错（ref 侧空）——修复扩展为
  self 模式下任一侧空则复用另一侧。验证仅 `--ref-seq` / 仅
  `--query-seq` / 双侧 / FASTA 直接输入四者输出逐字节一致；g1 高重复
  基因组（17 块）同样一致。回归测试覆盖两方向。

## 复核 113（2026-08-05 后续轮次）

### .pgi 显式输入参数忽略（文档化行为，非缺陷）

* `.pgi` 显式输入 + 冲突 `-k 20` / `--smer 10` / `--window 8` 被静默
  忽略（exit 0）。docs/align-pgi.md 第 36 行明确说明 "`--k/--smer/
  --window` apply only to genome-sequence inputs; .pgi inputs carry their
  parameters in the index header"——文档化预期行为，非静默错误。sibling
  索引路径的冲突报错是额外保护（自动复用索引时防意外），两路径行为
  均有依据。

## 复核 114（2026-08-05 后续轮次）

### fa split name 折行输入（未发现新问题）

* 折行 FASTA 经 `fa split name` 输出为单行序列；`rept trf` 对折行输入
  正常检出（chrC 1-600）。

## 复核 115（2026-08-05 后续轮次）

### s-kmer 尾部 run 漏报（文档化设计取舍，记录不修）

* Profex `-z` 从不闭合 read 的最后一个 run（end/depth 省略），s-kmer
  （min_depth=2）按设计保守丢弃尾部。实测"拷贝-间隔-拷贝"与"双拷贝
  串联"均只报第一份拷贝（1-1000），序列末尾的最后一个重复拷贝漏报；
  e-kmer（无深度阈值）用染色体长度闭合尾部，不丢。
* 真实影响量化：mg1655 尾 run 起点 4641601（约 52 bp），低于
  min-len 100 会被 excise 过滤，实际影响有限。行为与 repeat.rs 文档
  "conservatively dropped since its depth is unknown" 一致，属已知
  取舍，不修。

## 复核 112–115 验证

* 全量 1234 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 116–121（2026-08-05 后续轮次）

### 金标准坐标验证：pgi 块边界 seed 锚定损失（引擎语义，文档补充）

构造已知真值的 SD 基因组（5000 bp 随机，位置 1000-2000 与 3000-4000 为
两份完全相同的 1000 bp 拷贝）：

* `sd search --engine lastz` → 2 块，每块精确 1000 bp（1000-2000 ↔
  3000-4000），通过 min-len 1000 检出；
* `sd search --engine pgi` → 0 块：pgi 链块为 989-999 bp（seed 锚定
  边界在拷贝边缘内 1-11 bp，无 seed 恰在拷贝边界），恰低于阈值被滤。

量化（不同拷贝长度）：pgi 块长与真实拷贝长度差 0-11 bp，损失随序列
布局变化（seed 覆盖边缘）。identity 100% 时同样存在，与复核 99 的
低身份子块问题同属"pgi 块边界/覆盖近似"范畴，非逻辑 bug。

文档补充：docs/sd.md search 节说明高 identity 且真长恰在 min-len 附近
的拷贝也可能因 seed 边界损失差几 bp 被滤，建议降 min-len 或 lastz。

### 其他验证（未发现新问题）

* `align pgi -f` 边界：freq=1 全部跳过（0 块，极端参数语义）、
  10/100/100000 一致（17 块）。
* `rept e-align --min-identity`：0.5/0.7 相同（89 区间）、0.9 更严
  （59 区间）；`psl.ident()` 为 gap-compressed identity（不含插入），
  与 e-align 文档一致（与 sd 的含插入 identity 区分）。
* soft-mask 警告实测触发（小写基因组 e-align 打印提示）。
* 命名链：sd search PSL（无物种前缀）→ chainnet PAF（物种.染色体）→
  cluster/decompose/cover（物种列 + 染色体）一致。
* 文档三方核对：sd run / align pgi / align lastz / rept 全命令的
  docs 与 help 参数、默认值、说明逐项一致。

## 复核 116–121 验证

* 全量 1234 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 122（2026-08-05 后续轮次）

### 静态走查 panic 风险（未发现新问题）

* 本族生产代码中所有索引访问（`loc.rs` load_loc fields[0..2]、
  `cover.rs` read_elems f[0..7]、`repeat.rs` parse_trf_output fields[0..1]）
  前均有长度检查；`merged.last_mut().expect` 有非空前置保证；其余
  unwrap/expect 均在测试代码。无生产 panic 风险。

## 复核 123（2026-08-05 后续轮次）

### 缺陷 31：损坏的 FastK 缓存被静默复用 → e-kmer 空输出

`cache_is_fresh` 只检查缓存存在 + mtime。实测将 `lib.fa.gz.repeat.k17.ktab`
截断为 100 字节后（`.complete` 标记和 part 文件完好）：

* 日志显示 "reused repeat table"（缓存被判新鲜）；
* FastK 静默读取损坏表 → e-kmer 输出空 runlist（mg1655 原 48 区间全部
  丢失），无任何警告——比报错更隐蔽的静默错误输出。

根因：`.complete` 标记本是 `.ktab` 的完整副本（`save_repeat_cache`
原子复制），但新鲜度检查不校验大小。修复：`cache_is_fresh` 增加
`.ktab` 与 `.complete` 大小一致性校验，不一致（截断/改写）即视为陈旧
重建。实测：截断缓存后日志显示 "FastK on repeat"（重建），输出与原
48 区间逐字节一致，新表 524304 字节完整。

回归 `truncated_cache_table_is_not_fresh`：完整缓存判新鲜、截断后判陈旧。

## 复核 122–123 验证

* 全量 1236 测试通过（含新回归测试）；`cargo fmt --check` 与
  `cargo clippy --all-targets -- -D warnings` 干净。

## 复核 124–138（2026-08-05 后续轮次）

### 索引完整性盘点（未发现新问题）

* `.pgi` 兄弟索引损坏（截断保留 magic）→ PgiStream 打开报 "truncated
  index records"，非静默错误。
* `.loc` 损坏（截断）→ load_loc 空索引 → fetch 报 "not found in the
  .loc index file"，非静默错误。
* FastK 缓存：复核 123 已加 `.ktab`/`.complete` 大小校验。

### 边界参数与静态走查（未发现新问题）

* `rept trf --min-score`：非整数/负数/超大/NaN 分别报错，校验完整
  （负数需 `=` 形式传参，clap 标准行为）。
* `rept e-kmer -k 1` → FastK 报错传播；`-k 64` 正常；k17/k64 缓存命名
  隔离，各自复用一致。
* `align pgi --band 0/1/100000`、`--max-gap 0/100000`：无 panic/OOM，
  输出随参数语义变化。
* `fa mask` 空 runlist（`{}`）→ 原样输出；含不存在染色体 → 安全处理。
* `rept s-align --chunk-records 1`、`--window 50 --step 50`：输出正常。
* 静态走查：本族生产代码无未保护索引访问/unwrap；`validate_contigs`
  数量/名称/长度三重校验；`sibling_pgi_path` gz 命名（复核 51）与
  mtime 检查（复核 62）均正确；LAV→PSL 1-based→0-based 转换
  `checked_sub` 防溢出、负跨度/span 不匹配/零跨度检查齐全，`-` 链坐标
  翻转与 UCSC lavToPsl 一致。
* `align lastz --show-preset`：有 preset 正常显示、无 preset 报错。

### 倒位重复金标准（复核 121 文档语义的实例确认，非新缺陷）

* 1000 bp 反向拷贝（1500-2500 ↔ 3000-4000）：lastz 检出 2 块 1001 bp
  （strand `-`），pgi 检出 2 块 996 bp——边界损失 4 bp，低于 min-len
  1000 被滤。与复核 121 文档（seed 锚定边界损失 1-11 bp）一致，非新
  缺陷；正向 4 拷贝家族（复核 130）全部检出、坐标偏差 ≤11 bp。

## 复核 124–138 验证

* 全量 1236 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 139（2026-08-05 后续轮次）

### 缺陷 32：sd cluster 输出目录残留旧 cluster 文件

实测向含 `cluster_1.fa`/`cluster_2.fa`/`cluster_3.fa` 的目录重跑
`sd cluster`（本次仅 1 个 cluster）：旧 `cluster_2.fa`/`cluster_3.fa`
残留，下游 `sd decompose` / 手动 `sd run` 会**静默消费陈旧家族**作为
当前输出（与复核 123 的 FastK 缓存损坏同类：静默错误数据）。
`sd run` 内部用固定 tempdir 免疫，但手动工作流受影响。

修复：写输出前清理 outdir 中 pgr 自身命名模式的 `cluster_<u32>.fa`
（仅此模式，其他文件不动）。实测：旧 cluster 文件删除、新
`cluster_1.fa` 覆盖、`notes.txt` 等无关文件保留。

回归 `stale_cluster_files_are_removed`：3 个旧 cluster + 无关文件 →
仅 cluster_1.fa 与无关文件保留。

## 复核 139 验证

* 全量 1239 测试通过（含新回归测试）；`cargo fmt --check` 与
  `cargo clippy --all-targets -- -D warnings` 干净。

## 复核 140（2026-08-05 后续轮次）

### align lastz 输出目录残留（记录不修）

* 重复 `align lastz -o dir` 时旧 LAV 残留（`stale.lav`）。影响链短：
  `sd run`/`rept s-align`/`sd search lastz` 均用临时 workdir 免疫，仅
  手动多次使用同一 `-o` 目录 + 手动消费 LAV 时可能读到残留。与复核 139
  的 cluster 不同（cluster 直接污染下游自动消费），且 LAV 是通用扩展名，
  清理模式易误伤用户保留的中间结果。记录不修。

## 复核 141（2026-08-05 后续轮次）

### sd cover 畸形 elems.bed（未发现新问题）

* 字段数 <8 的行跳过；begin 列非数字报 "invalid digit found in string"
  （parse 错误传播，非 panic）。

## 复核 142–143（2026-08-05 后续轮次）

### 缺陷 33：空 FASTA 输入触发 FastK SIGSEGV（预检友好报错）

* 空 repeat 库（`>empty1` 无序列）喂 `rept e-kmer` → FastK SIGSEGV，
  pgr 报 "terminated by signal: 11"（复核 102 已记录外部工具崩溃，但
  错误信息像 pgr 自身崩溃，不友好）。
* 全 N / 4 bp 极小库 → FastK exit 1（工具自身拒绝）。

修复：`run_repeat_pipeline` 在 FastK 前预检输入是否有非空序列
（`has_sequences`），空则报友好错误（"repeat library FASTA has no
sequences" / "input genome FASTA has no sequences"）。全 N/极小库仍走
FastK（工具行为，exit 1 可接受）。

回归 `sequence_less_fasta_is_detected`：仅头/无记录判空、正常序列判非空。

## 复核 140–143 验证

* 全量 1240 测试通过（含新回归测试）；`cargo fmt --check` 与
  `cargo clippy --all-targets -- -D warnings` 干净。

## 复核 144–148（2026-08-05 后续轮次）

### 空输入与其他输入形式（未发现新问题）

* `rept trf` 仅头 FASTA → 空输出 `{}`（exit 0）。
* `rept s-align` 仅头 FASTA → lastz 报 "contains no sequence"，错误
  传播（可读，非 panic）。
* `rept e-align` 仅头基因组 → 空 runlist `{}`；`align pgi` 对仅头
  FASTA 建出 0 长 contig 索引 → 0 块（0 字节 FASTA 才报 "index has no
  contigs"，两者行为各有依据）。
* `rept e-align --keep-index`：g.pgi/lib.fa.pgi 兄弟索引生成与复用
  （"reusing reference/query index"），两次输出逐字节一致。
* `align pgi` `.2bit` 输入：17 块，与 FASTA 直接输入逐字节一致；
  混合输入（2bit ref + .pgi query + `--query-seq`）18 块（显式对不丢
  self-hit，合理）；genome 输入配 `--ref-seq` 的正确报错（resolve_seqs
  校验）。

## 复核 144–148 验证

* 全量 1240 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 149–150（2026-08-05 后续轮次）

### 主命令注册与最终冒烟（未发现新问题）

* `pgr sd`（7 子命令）/ `pgr align`（2）/ `pgr rept`（5）主命令帮助
  完整，子命令全部注册。
* mg1655 冒烟：`sd run` 118 行 elementary BED、`rept s-kmer` 170 区间、
  `rept trf` 正常、`rept e-align` 89 区间——与历史基线一致，复核
  30–33 的修复不破坏主流程。

## 复核 149–150 验证

* 全量 1240 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 151–155（2026-08-05 后续轮次）

### 缺陷 34：sd run --engine lastz 输出重复 elementary 行

mg1655 `sd run --engine lastz` 输出 120 行含 5 个完全重复行（如
607996-609351 出现两次）。追踪：

* `sd search lastz` 282 块 → PAF 83 行（无重复）→ cluster 21 个文件；
* 单个 cluster 文件内出现 end 差 1 bp 的两个头（如
  `+#729464#729581` 与 `+#729464#729582`）——lastz 互反块的坐标抖动；
* decompose 对两者投影到完全相同的 elementary 区间（729464-729572）→
  合并后重复行。

修复位置权衡：decompose 层按投影坐标去重会破坏既有语义（相同头但不同
序列的两条记录应各输出一行，`decompose_detects_shared_fragment`），故在
`sd run` 合并层按 renumber 后的完整行去重（`push_unique_elem`）。实测：
mg1655 lastz run 115 行、0 重复（原 120 行含 5 重复）。pgi run（118
行）不受影响。

回归 `duplicate_elem_rows_are_emitted_once`（cmd 单元测试）+ cli_sd 14
测试全过。

### 其他验证（未发现新问题）

* `align pgi` 输出与 `--parallel` 1/2/8 无关（FASTA 与 .pgi 路径均逐字节
  一致）；`rept e-align --workflow tube --keep-index` 两次一致且复用索引；
  e-kmer 缓存 mtime 调旧（2020）触发重建、`.ktab`/`.complete` 重建后
  大小一致；mg1655 `sd run --engine lastz` 120→115 行（去重后）。

## 复核 151–155 验证

* 全量 1241 测试通过（含新回归测试）；`cargo fmt --check` 与
  `cargo clippy --all-targets -- -D warnings` 干净。

## 复核 156–158（2026-08-05 后续轮次）

### 去重边界与参数语义（未发现新问题）

* 缺陷 34 去重逻辑边界：跨 cluster 的相同坐标但不同 set_id 不去重
  （不同家族应保留）；相同坐标不同 strand 不去重（不同拷贝）；pgi
  run 的 118 行输出不受影响。
* `rept e-align --min-len`：1 → 251 区间（最小 7 bp）、50 → 89 区间、
  500 → 47 区间（最小 703，excise 后合并区间），阈值语义正确。
* 多染色体 `rept s-align`：chr1/chr2 各自输出区间（2/1），正常。

## 复核 156–158 验证

* 全量 1241 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 159–161（2026-08-05 后续轮次）

### 目录 self、确定性、文档完整性（未发现新问题）

* `align lastz` 目录省略 query → self 模式：2 targets 4 jobs 中 2 个
  对角输出（[t1]vs[t1].lav、[t2]vs[t2].lav），正确。
* `sd search --engine lastz` 三次运行输出逐字节一致。
* docs/ 下 sd/rept/align/align-pgi/align-lastz/lav/psl 及 formats/ 文档
  齐全；rept.md 的 `../notes/design/repeat-masking.md` 引用存在。

## 复核 162（2026-08-05 后续轮次）

### 文档修复 5：--ref-seq 内容一致性说明

* 实测：`.pgi` 索引 + 同长度同名但内容不同的 `--ref-seq` 通过
  validate_contigs（仅数量/名称/长度校验）→ 精化输出碎片化垃圾对齐
  （549 匹配/295 错配/119 小块）而非报错。受影响面仅手动 `align pgi`
  高级用法（`sd search/cross` 拒绝 `.pgi` 输入）；下游 SD 过滤
  （0.90）与 e-align（0.70）会滤掉该低 identity 块。
* 修复：align-pgi.md 明确 `--ref-seq` 校验范围（contig 表）并要求序列
  与索引来源一致（自动 sibling 路径由 mtime 检查保证）。未实现 k-mer
  内容校验（syncmer 哈希对比复杂度高、阈值易误报，文档说明足够）。

## 复核 163–164（2026-08-05 后续轮次）

### cover 贪心边界与 chainnet 坏行（未发现新问题）

* `sd cover`：多 set 竞争覆盖选覆盖最多者、跨染色体覆盖（chr2 成员
  覆盖 chr2 hit）正确、无覆盖 hit 的 set 标 non-core——运行时验证通过。
* `pl chainnet` 畸形 PSL 行静默跳过（`Err(_) => continue`，UCSC header
  兼容设计，注释说明）；负 target strand 有 warn；2 块输入 → MAF 2 行。

## 复核 159–164 验证

* 全量 1241 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 165–170（2026-08-05 后续轮次）

### 参数语义纵深验证（未发现新问题）

* `rept s-kmer --fill-kmer`：0/10 不合并 300 bp 间隙（1-201、501-701）、
  400 合并（1-701）——语义正确。注：早期测试用两个相同 gap（同一变量
  打印两次）导致 gap 深度 2、fill 无差异，属测试输入错误，非缺陷。
* `sd search` 边界：`--min-len 0` 17 块、`--min-identity 0` 拒绝
  （(0,1] 校验）、1.0 仅 100% identity 块（11 块）。
* `rept trf` 特殊字符名（chr1:alt、chr2.5）：输出键正确、区间完整。
* `rept e-align` 特殊字符库名（IS1:abc、IS2.5）：索引构建与检出正常；
  默认 `-f 100` 对 ACGT 低复杂度库 0 命中（四联体 40-mer 出现 ~240 次
  > 100 全跳过）是预期频率过滤，`-f 100000` 检出 1-1020/1979-3000。
* `rept s-align --min-depth` 语义：双拷贝 + 50% 重叠窗口（基线深度 2）
  → min-depth 2 全过、4 检出 fam 区、8 空；step=50（75% 重叠，基线
  深度 4）→ 4 全过、8 检出 fam 区。窗口重叠度与深度阈值交互正确。

## 复核 165–170 验证

* 全量 1241 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 171–173（2026-08-05 后续轮次）

### lav/psl 解析容错与 fill-fragment 语义（未发现新问题）

* `lav to-psl` 对畸形 LAV（错误 sizes 行、纯垃圾）静默输出空（exit 0）：
  LavReader 跳过空行/注释、解析已知 stanza，纯垃圾行（不以 `{` 结尾）
  静默忽略、未知 stanza 有 warn。lastz 输出不会畸形，属容错设计，记录
  不修。
* `psl to-chain` 对垃圾行有 warn（"skipping unparseable psl line"），
  行为合理（parse_or_warn 非 strict）。
* `rept s-kmer --fill-fragment`：0 不合并 150 bp 间隙（1-200、351-550）、
  200 合并（1-550）、600 同（合并 150；500 bp 间隙不合并）——语义正确。

## 复核 171–173 验证

* 全量 1241 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 174（2026-08-05 后续轮次）

### align pgi --merge-gap 语义（未发现新问题）

* 构造含 80 bp 插入的拷贝对（IS 断裂场景）：`--merge-gap 0/100/1000`
  均 10 块——同对角线断裂不合并，与 align-pgi.md 文档（"only chains
  whose diagonals differ are merged, same-diagonal gaps stay independent
  blocks"）一致；大值无误并不同拷贝。

## 复核 174 验证

* 全量 1243 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 175–181（2026-08-05 后续轮次）

### 输入形式与参数透传纵深验证（未发现新问题）

* 多染色体 e-kmer/s-kmer：Profex 逐染色体调用正常，chr1/chr2 各自输出
  区间；e-kmer 空输出（fam 不与 IS 库匹配）合理。
* `rept e-align --min-shared`：16 检出更宽（377-730）、40/100 收缩到
  核心区（516-707）——10% 分歧拷贝的链 seed 阈值语义正确。注：早期
  测试把库与基因组设为同一文件（库含 chr1 自身导致全长完美匹配），
  属测试输入错误。
* `sd search`/`sd cross` 对 `.pgi` 输入清晰拒绝（"needs genome FASTA
  (plain or .gz); a .pgi index aligns without extension sequences"）。
* `rept trf` 二进制输入：`fa split name` 报 UTF-8 错误并带调用上下文
  传播。
* `rept s-align` `.gz` 输入：与 plain 输出逐字节一致。
* `align lastz --lastz-args "K=1000"`：正确透传给 lastz 命令行。

## 复核 175–181 验证

* 本族测试全绿：cli_align_pgi 18 / cli_align_lastz 2 / cli_sd 14 /
  cli_rept 16 / lib 551；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。全量测试仍受并行会话 fas/poa 外部工具测试
  失败影响（与本族无关）。

## 复核 182–184（2026-08-05 后续轮次）

### 精化 DP 与链合并阈值审查（未发现新问题）

* `extend_chain`/`chain_windows`：a/b 序列边界检查、`q1 > b_len`、
  `q0 >= q1`、`step == 0` 均返回空（防越界/死循环）；`-` 链 rev_comp
  处理正确；`dp_band = min(diag_span+32, 128)`（复核 76 的 OOM 上限
  保留）。
* `align_banded_local`：带限仿射 gap（M/I/D 三状态），band 与序列边界
  相交（修复过 offset wrap），空输入返回 None。
* 链合并中段同源校验 200 bp 阈值边界：非同源 gap 在 199/200/201 bp
  两侧均保持 2 块（不桥接），行为一致无突变。

## 复核 182–184 验证

* 本族测试全绿（同上）；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 185–187（2026-08-05 后续轮次）

### 组合场景与序列提取金标准（未发现新问题）

* `sd run` gz 基因组 + lastz 引擎：115 行、0 重复（与 plain lastz run
  一致，缺陷 34 去重对 gz 路径同样生效）。
* `sd cluster` 序列提取金标准：头 `?#chr1+#500#800` /
  `?#chr1+#1000#1300` 的序列与 genome[500:800] / [1000:1300]
  逐碱基完全一致（300 bp）。注：早期脚本把两个记录拼接致长度误报，
  按记录拆分后验证通过。

## 复核 185–187 验证

* 本族测试全绿（同上）；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 188–189（2026-08-05 后续轮次）

### 文档示例可运行性（未发现新问题）

* rept.md 示例命令（e-kmer/s-kmer/trf/s-align `输入 > 输出`）与 CLI
  匹配；引用的数据文件（mg1655.fa.gz、mg1655.chr.sizes、
  tncentral.fa.gz、mg1655.rm.gff）全部存在。
* align-pgi.md 五个示例（gz+2bit 混合、`pgr pgi build` + `.pgi` +
  `--ref-seq/--query-seq`、`--keep-index`、`-f/-c/-s/--band` 参数、
  单输入 self）全部与 CLI 匹配，且均经实测可运行（复核 147/148/110
  等已覆盖）。

## 复核 188–189 验证

* 本族测试全绿（同上）；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 190–192（2026-08-05 后续轮次）

### 缺陷 35：align pgi 自动索引小写归一化 → 全零块

构造含大小写混合拷贝的基因组（fam 大写 + fam 小写，非重叠）：

* 修复前 `align pgi` 输出 2 块 match=0/mismatch=0/rep=0 的**全零块**
  （298 bp 无效数据）。根因：自动索引 `build_from_seqs` 的碱基编码
  大小写不敏感（a/A 同码）→ 小写拷贝与大写拷贝共享 seed → 链存在；
  但扩展 DP 的子矩阵大小写敏感（'A' != 'a'）→ 评分失败 →
  `extend_chain` 回退 `chain_to_psl` raw 块（全零）。k-mer 编码与 DP
  大小写语义不一致导致"找到但无法评分"的中间状态。
* `build_from_path(mask=true)` 已有 `harden_soft_mask`（小写→N，FastGA
  `-M` 语义），但 `build_from_seqs`（align pgi 自动索引专用）无 mask
  参数——两路径行为不一致。

修复：`build_from_seqs` 增加 `mask` 参数（与 `build_from_path` 一致），
align pgi 自动索引传 `true`（跳过小写）。实测：混合大小写输入 0 块
（不再输出全零块）、全大写对照 2 块正常；e-align 对小写基因组仅检出
大写拷贝并保留 soft-mask 警告（"results will be underestimated"），
行为与文档一致。

回归 `command_align_pgi_lowercase_copy_has_no_all_zero_blocks`：小写拷贝
输入不产生任何 match+mismatch+rep=0 的块。

## 复核 190–192 验证

* 本族测试全绿：cli_align_pgi 19 / cli_sd 14 / cli_rept 16 / lib 553；
  `cargo fmt --check` 与 `cargo clippy --all-targets -- -D warnings` 干净。
  全量 741 通过，仅并行会话 fas/poa 外部工具 4 测试失败（与本族无关）。

## 复核 193–194（2026-08-05 后续轮次）

### 软掩码语义双引擎统一确认（未发现新问题）

* 1100 bp 拷贝 + 小写拷贝：修复后 pgi 与 lastz 对小写混合基因组行为
  **完全一致**（全大写 2 块、混合 0 块——小写作为软掩码跳过），缺陷 35
  修复使双引擎遵循同一 FastGA/UCSC 软掩码语义。
* `rept trf` 对小写/大写串联重复均检出（TRF 大小写不敏感，合理）。
* N 段基因组（fam-N-fam）：pgi 链在 N 处断裂为 2 块（N 无 seed），
  全 N 基因组 0 块——N 处理正确。

## 复核 193–194 验证

* 本族测试全绿（同上）；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 195–197（2026-08-05 后续轮次）

### 输出数据正确性抽查（未发现新问题）

* `sd search` PSL 的 qSize/tSize 等于基因组实际长度（g1 12614），
  坐标在长度内——字段正确。
* `rept trf` 输出（1-based inclusive）与 `fa range` 提取序列一致
  （ACGT 串联 1-1200 提取正确）。
* `rept e-kmer` mg1655 区间（15388-16730 等 48 个）经 `fa range`
  正常提取，与已知 mg1655 重复特征一致。

## 复核 195–197 验证

* 本族测试全绿（同上）；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 198–199（2026-08-05 后续轮次）

### mask 端到端与随机化压力测试（未发现新问题）

* `rept e-align` → `fa mask` 软掩码正确性：89 个 runlist 区间内
  60251 bp 全部小写（0 例外）、区间外抽样 0/300 小写、mask 后长度不变。
* 随机化压力：3 个重复家族（1200 bp）× 2-3 拷贝（混合链向、随机间隔
  100-400 bp）→ `sd run` 16 行 elementary BED、三次运行逐字节一致、
  无崩溃；CORE 标记与 - 链投影正确。

## 复核 198–199 验证

* 本族测试全绿（同上）；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 200–202（2026-08-05 后续轮次）

### sd/mod.rs 走查、s-align preset、TRF 区间合并（未发现新问题）

* `libs/sd/mod.rs`：is_pgi_input（magic + 扩展名）、psl_block_len
  （含插入 max(0) 防负）、psl_identity、passes_sd_filters 全部正确。
* `rept s-align` preset set01/set03/set07：输出差异合理（lastz 灵敏度
  差异，set03 合并边界、set07 更细分段）。
* `rept trf` 多段 tandem（ACGT/TTTTGGGG/GATTACA 交替）：TRF 手动输出
  5 个独立段，pgr runlist 将相邻段合并为 1-970（IntSpan 相邻区间合并
  语义，fa mask 掩码效果等价）——正确。

## 复核 200–202 验证

* 本族测试全绿（同上）；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 203–205（2026-08-05 后续轮次）

### 文档修复 6–7：小写（软掩码）处理说明

* 缺陷 35 修复后 align pgi 自动索引按 FastGA `-M` 语义跳过小写，但
  docs/align-pgi.md 与 docs/sd.md 均未说明。补充：
  - align-pgi.md Notes：自动索引小写→N 无 seed/块，`pgr pgi build
    --mask` 同语义；
  - sd.md Notes：pgi（`-M` 语义）与 lastz（小写视为掩码）都不比对小写，
    软掩码的 SD 拷贝不被检出，建议先 `tr a-z A-Z`。

### sd cross 反向拷贝（未发现新问题）

* query 含 rc(fam) 反向拷贝：`sd cross` 以 `-` 链正确检出
  （q[503:1686] ↔ t[1014:2197]，1183 匹配），坐标投影正确。

## 复核 203–205 验证

* 本族测试全绿（同上）；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 206（2026-08-05 后续轮次）

### 端到端 BISER 工作流（未发现新问题）

* 构造含 TE（400 bp ×2）+ SD（1200 bp ×2）的基因组：
  `rept e-kmer` 检出 TE（1-400、3401-3800）→ `fa mask` 掩码 →
  `sd run` 掩码后检出 SD（598-1791、1998-3191，1193 匹配）——TE 被
  掩码不干扰、SD 完整检出，文档声称的"掩码后 SD 检测"流程验证通过。

## 复核 206 验证

* 本族测试全绿（同上）；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 207（2026-08-05 后续轮次）

### 双引擎 sd run 输出差异深度分析（非缺陷）

mg1655 双引擎 `sd run`：pgi 118 行 vs lastz 115 行，按
(chrom,begin,end,strand) 完全相同的区间 0 个、重叠 ≥50% 仅 51/118——
表面看差异很大。逐层排查：

* **search 层面完全一致**：pgi 232 块全部与 lastz 282 块重叠 ≥50%
  （232/232）；lastz 独有的 8 块为 93.9-97.6% identity 的正常 SD
  （lastz 额外灵敏度，复核 121 已记录 pgi 近阈值漏检）。
* **差异传导**：lastz 多检出的块改变 cluster 成员 → decompose 的共享
  k-mer 家族划分与片段投影不同（set_id 与坐标都变）→ elementary BED
  差异大。BISER 语义允许两引擎输出不同（可互换替代引擎）。
* **抽查自洽性**：pgi 未匹配短区间 225361-225627（267 bp，序列唯一）
  的 10-mer 与基因组 388 个外部位置共享——真实家族成员，非虚假输出。
  每个引擎的输出各自自洽且坐标正确（金标准复核 130/187 已覆盖）。

结论：双引擎差异是灵敏度差异 + decompose 对 cluster 划分的敏感性，
非缺陷。

## 复核 207 验证

* 本族测试全绿（同上）；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 208–210（2026-08-05 后续轮次）

### 2bit 组合与并行确定性（未发现新问题）

* `align pgi --ref-seq` 传 `.2bit`：2 块，与 FASTA ref-seq 逐字节一致。
* `pgr pgi build` 从 2bit 建索引 + `--ref-seq` 传 2bit：2 块，与 FASTA
  参考一致——2bit 索引与序列全组合正常。
* `rept s-align --parallel 1/4`：输出逐字节一致（lastz 并行确定性）。

## 复核 208–210 验证

* 本族测试全绿（同上）；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 211–215（2026-08-05 后续轮次）

### 全量基线恢复与冒烟确认（未发现新问题）

* 并行会话的 fas/poa 修改完成，**全量 1250 测试全部通过**（无失败），
  `cargo fmt --check` 与 `cargo clippy --all-targets -- -D warnings` 干净。
* mg1655 全量冒烟：sd run pgi 118 行、lastz 115 行（去重后）、e-kmer
  48、e-align 89、s-kmer 170、trf 正常——与历史基线一致，全部修复
  （缺陷 30–35）在最终代码下正常。
* `sd run --engine pgi --preset set01`：preset 正确忽略（输出与无 preset
  一致）。
* `PipelineCtx` 走查：current_exe 定位 pgr、tempdir 生命周期、CwdGuard
  drop 恢复——实现正确。
* 本族文件（align/rept/sd + pgi/pl/loc/sd libs）修改均为历轮审核修复，
  并行会话未触碰。

## 复核 211–215 验证

* 全量 1250 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 216–218（2026-08-05 后续轮次）

### 组合场景与参数覆盖（未发现新问题）

* `sd cross` gz 输入 + lastz 引擎：3 行，与 plain+lastz 逐字节一致。
* `rept e-kmer` 缓存复用 + 参数变化（`--fill-kmer 100 --min-len 500`）：
  复用缓存（"reused repeat table"）且输出 45 区间（默认 48）——缓存与
  runlist 参数正交，参数变化正确生效。
* `align lastz --lastz-args="--querydepth=keep,nowarn:5"`：追加到 lastz
  命令行（默认 50 + 追加 5），lastz 最后参数胜出（5 生效）——依赖 lastz
  覆盖语义，与文档（"overrides preset"）一致；值以 `--` 开头时需 `=` 形式
  传参（clap 标准行为）。

## 复核 216–218 验证

* 全量 1250 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 219–220（2026-08-05 后续轮次）

### 重复库边界（未发现新问题）

* 库含重复 contig 名（IS1 ×2）：e-align 正常检出（索引按 contig id
  区分，重复名无冲突），不崩溃。
* 库含空 contig 名（`>` 无名字）：pgi 索引构建报 "missing name"
  （FASTA 解析错误），错误清晰传播，非 panic。

## 复核 219–220 验证

* 全量 1250 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 221–222（2026-08-05 后续轮次）

### 缺陷 36：默认参数静默复用不同 k 的兄弟索引

`resolve_side` 的缓存参数冲突检查只覆盖**命令行显式传**的
`-k/--smer/--window`（`ValueSource::CommandLine`）。实测：

1. `align pgi g.fa -k 20 --keep-index` 建 k=20 缓存；
2. 之后 `align pgi g.fa`（未显式 -k，默认 40）→ 日志 "reusing reference
   index g.pgi" → **静默用 k=20 索引跑默认 k=40 语义的比对**（输出与
   k=40 不同，用户无感知）；
3. 显式 `-k 40` 则报错——两条路径行为不一致。

修复：删除 `explicit(...)` 条件，**总是**检查当前解析值（显式或默认）
与缓存索引参数的一致性。验证：
* 默认 k=40 vs 缓存 k=20 → 报错（修复）；
* 显式 k=20 vs 缓存 k=20 → 复用（正常）；
* 默认建 k=40 缓存 + 默认运行 → 复用且输出一致（正常路径不破坏）；
* 显式 k=40 vs 缓存 k=20 → 报错（原有行为保留）。

回归 `command_align_pgi_default_kmer_conflicts_with_cached_index`。

## 复核 221–222 验证

* 全量 1250 测试通过（含新回归测试）；`cargo fmt --check` 与
  `cargo clippy --all-targets -- -D warnings` 干净。

## 复核 223–224（2026-08-05 后续轮次）

### 缺陷 36 修复的对称与传播验证（未发现新问题）

* `--smer`/`--window` 对称生效：smer10 缓存 + 默认 smer8 → 报错、
  smer10 复用正常、window8 显式（缓存 window5）→ 报错。
* `rept e-align --keep-index` 换 k：默认 k 运行经 cmd_lib 传播 align
  pgi 的冲突错误（"conflicts with the cached index ... (k=20)"），
  用户可读。

## 复核 223–224 验证

* 全量 1253 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 225–226（2026-08-05 后续轮次）

### 5 Mb 大输入性能与稳定性（未发现新问题）

* 5 Mb 随机基因组 + 1200 bp 双拷贝：`sd search` 6.8s（2 块）、
  `sd run` 全链路 7.7s（2 行 elementary BED，坐标 2000002-2001190 /
  2501202-2502390，无重复）——大输入无崩溃、性能合理。

## 复核 225–226 验证

* 全量 1253 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 227（2026-08-05 后续轮次）

### TRF --max-period 上限（外部工具限制，记录不修）

* `--max-period` 2000（默认）正常、10000+ TRF SIGSEGV（signal 11）——
  与已知"TRF 2500bp 周期 + max-period≥2600 SIGSEGV"同类外部工具限制；
  cmd_lib 捕获报错（"terminated by signal: 11"），pgr 无 panic。TRF 的
  精确上限未知（2000 安全、10000 崩溃），pgr 无法可靠预校验，记录不修。

## 复核 227 验证

* 全量测试通过（当前 1182，并行会话调整后）；`cargo fmt --check` 与
  `cargo clippy --all-targets -- -D warnings` 干净。

## 复核 228–230（2026-08-05 后续轮次）

### self 模式与 min-shared 核心逻辑（未发现新问题）

* `drop_self_hits`：丢弃同 contig + 同位置 + 正向链的精确自比对 seed，
  反向链（strand==1）/跨位置/跨 contig 保留——逻辑正确，有回归测试
  （drop_self_hits_filters_exact_identity）。
* `effective_min_shared`：显式值 min(v,k)、tube 默认 min(12,k)、greedy
  默认 k——正确。
* `--min-shared 0` → 报错 "min_shared must be in 1..=40"；`--min-shared
  100000` 截断到 k=40（17 块正常）——校验完整。

## 复核 228–230 验证

* 全量 1253 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 231–232（2026-08-05 后续轮次）

### lastz 参数构建与 preset 数据（未发现新问题）

* `build_common_args`：基础参数（querydepth/format/markend/ambiguous）+
  preset 参数（跳过 Q= 用矩阵临时文件）、矩阵句柄由调用方持有存活——
  正确。
* 7 个 preset（set01-set07）定义完整（UCSC 参数 + 矩阵引用），
  `find_preset`/`preset_names` 正确。

## 复核 231–232 验证

* 全量 1253 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 233–234（2026-08-05 后续轮次）

### .pgi 索引读取健壮性走查（未发现新问题）

* `PgiMmap::open`：只读 mmap（map_copy_read_only，SAFETY 注释）、记录
  区域 checked_mul/checked_add 溢出检查、截断检查（map.len() >=
  rec_end）——健壮。
* `parse_header_bytes`：magic/version/参数/计数校验、kmer_bytes/
  pos_bytes/cont_bytes 合理性、n_contigs ≤ u16::MAX、n_records ≤
  u32::MAX（防 allocation overflow）、contig 名 UTF-8 校验、take_bytes
  越界检查——复核 6 的 crafted .pgi 防护完整保留。

## 复核 233–234 验证

* 全量 1254 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 235–236（2026-08-05 后续轮次）

### PSL 解析与输出（未发现新问题）

* `Psl::from_str`：21 字段校验、u32/i32 类型解析带错误、数组解析、
  block_count 与向量长度一致性截断（畸形 PSL 容错不崩溃）——健壮。
* `Psl::write_to`：标准 PSL 18 字段 + 三数组逗号序列输出——正确。

## 复核 235–236 验证

* 全量 1255 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 237–238（2026-08-05 后续轮次）

### 链→PSL 与 gz 解压走查（未发现新问题）

* `chain_to_psl`：contig 名/长度取自索引表、反向链坐标翻转
  （reverse_range_pair + q_starts 用 b_len - q_end）、单块输出——
  正确；contig id 越界受 merge 约束（复核 6 修复的 crafted .pgi 防护
  仍在）。
* `plain_gz_to_temp`：BGZF 直接使用（索引可随机访问）、非 BGZF gz
  解压到 tempdir（.loc 索引需要）、tempdir 生命周期由 Option 持有——
  正确（复核 100 前的 plain-gz cluster 修复保留）。

## 复核 237–238 验证

* 全量 1255 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 239（2026-08-05 后续轮次）

### decompose 多家族分组（未发现新问题）

* 构造 3 拷贝家族 A + 2 拷贝家族 B + 单拷贝 C 的 cluster FASTA：
  decompose 输出 A 3 行 set_id=1、B 2 行 set_id=2、C 无输出（单拷贝
  不形成 elementary）——MIN_SHARED_KMERS 分组与 set_id 划分正确。

## 复核 239 验证

* 全量 1255 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 240（2026-08-05 后续轮次）

### syncmer 实现走查（未发现新问题）

* `syncmer_dna` 流式 canonical hash + closed-syncmer 检测（复核 1 修复
  保留），但仅被 syncmer.rs 内部测试使用——生产路径（pgi build 的
  collect_one_contig）用内联实现（codes 表 N→4 打断窗口）。
* `encode_base` 对 N 返回 0（当作 A）与生产路径的 N→4 不一致，但
  syncmer_dna 非生产路径，不影响 align/rept/sd——记录观察。

## 复核 240 验证

* 全量 1255 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 241–242（2026-08-05 后续轮次）

### seed 合并与发射走查（未发现新问题）

* `merge_seed_hits`/`merge_seed_hits_from_stream`：validate_compatible、
  min_shared 1..=k 校验、并行 chunk 合并（par_chunks/par_bridge）、
  流式批处理——正确。
* `emit_entry_hits`：`freq >= cutoff` 过滤（FastGA 语义，非 >）、
  canonical key 过滤、前缀窗口（min_shared 起始）、最大共享前缀、
  扩展范围频率过滤、validate_record（复核 6 越界 contig 防护）——
  与 FastGA GIX 语义一致。

## 复核 241–242 验证

* 全量 1255 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 243（2026-08-05 后续轮次）

### runlist span 操作走查（未发现新问题）

* `span_op` 委托 IntSpan：`excise`（保留 ≥min_len 片段）、`fill`（合并
  ≤max_len 间隙）、cover/holes/trim/pad——span_len 用 i64 运算防 i32
  溢出，实现正确。

## 复核 243 验证

* 全量 1255 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 244（2026-08-05 后续轮次）

### 文档修复 8–9：去重与缓存参数一致性说明

* align-pgi.md：补充 sibling 缓存索引参数一致性——当前
  `-k/--smer/--window`（显式或默认）必须与缓存匹配，不匹配报错而非
  静默复用不同 seed（缺陷 36）。
* sd.md：补充 `sd run` 输出去重——近相同 cluster 拷贝（互反块 1 bp
  抖动）投影到相同 elementary 区间时只输出一次（缺陷 34）。

## 复核 244 验证

* 全量 1255 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 245（2026-08-05 后续轮次）

### LAV gap 解析走查（未发现新问题）

* `blocks_to_psl` gap 计算：q/t 块间正 gap 计 num_insert/base_insert、
  负 gap（重叠块）clamp 忽略——与 UCSC lavToPsl 行为一致；q/t gap
  独立计数（异常 LAV 不崩溃）。

## 复核 245 验证

* 全量 1255 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。

## 复核 246（2026-08-05 后续轮次）

### read_fasta 走查（未发现新问题）

* `read_fasta` 用 noodles `records()` 解析（畸形 FASTA 由解析器处理）、
  contig 名 UTF-8 校验——简单正确。

## 复核 246 验证

* 全量 1255 测试通过；`cargo fmt --check` 与 `cargo clippy --all-targets
  -- -D warnings` 干净。
