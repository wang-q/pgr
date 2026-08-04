# sd / rept / align 命令族代码审核记录（2026-08-04 / 2026-08-05）

对 `pgr sd`（8 命令 + `libs/sd`）、`pgr rept`（6 命令 + `libs/pl`）、
`pgr align`（pgi/lastz + `libs/pgi`、`libs/lastz`、`libs/fmt/lav`、
`libs/fmt/psl`、`alignment` DP）三个命令族约 8000 行代码及全部文档
（`docs/{sd,rept,align-pgi,align-lastz,lav,psl}.md`）进行审核。缺陷按类别
分组记录；关键修复均附回归测试，验证概况见文末"验证"一节。

审核范围：
- **sd**：`search` / `cross` / `align` / `cluster` / `decompose` / `cover` /
  `run`（+ `libs/sd`、`libs/pgi` 索引、tube/greedy 双工作流）
- **rept**：`e-kmer` / `s-kmer` / `e-align` / `s-align` / `trf`
  （+ `libs/pl/repeat`、FastK/Profex/TRF 外部工具封装）
- **align**：`align pgi` / `align lastz`（+ `libs/lastz`、`libs/fmt/lav`、
  `align_banded_local` 等 alignment DP）

审核重点：Zero Panic（畸形输入不 panic）、数据安全（`-o` 不得覆盖输入、
陈旧/损坏索引不得静默复用）、确定性（跨运行输出逐字节一致）、与外部参考
实现（FastGA / UCSC kent / lastz / TRF）的语义一致。

## 与外部参考实现的语义一致性核对

关键修复均对照官方源码复核，方向一致：

* `psl lift` 负链坐标提升：与 UCSC kent `pslLiftSubrangeBlat.c` 的
  `liftSide` 行为一致（子范围命名约定 pgr 1-based vs kent 0-based 为
  记录在案的有意差异）。
* greedy 链合并：与 FastGA `align_contigs` / `ALNchain.c` 的链化
  语义一致——同对角线纯间隔是两条独立链，仅对角线平移才缝合（pgr 的
  自身扩展）。
* pgi 索引 k-mer 频率过滤：`emit_entry_hits` 的 canonical key 过滤、
  `freq >= cutoff`（FastGA 语义，非 `>`）、前缀窗口 / 最大共享前缀 /
  扩展范围过滤，与 FastGA GIX 语义一致。
* 软掩码语义：`build_from_seqs` 新增 `mask`（小写→N，FastGA `-M` 语义），
  与 `build_from_path` 一致；pgi 与 lastz 对小写拷贝行为完全一致。
* LAV→PSL 转换：`blocks_to_psl` 的 q/t 正 gap 计 insert、负 gap clamp 忽略、
  1-based→0-based `checked_sub` 防溢出、`-` 链坐标翻转，与 UCSC
  `lavToPsl` 一致。
* `sd run` 的 chainnet：每靶位点保留一条最优链（同等分按序取一）是 UCSC
  chainnet 每靶位点取最佳链的标准语义（与 `pl chainnet` 共享）。
* 有意差异（已记录）：子范围命名 pgr 1-based vs kent 0-based；`stat`/
  `statop` 类输出格式差异不涉及本族。

## 排除的疑点（经核验无需修复）

* `sd run` 的 cluster set_id 重编号值域各簇两两不相交，不可能碰撞。
* 60,423 → 75,413 数据差异来自 tncentral 库更新与编译时序，非代码 bug
  （repeat-masking.md §2.3.5 已勘误）。
* sd cluster minus 链序列提取：按 pgr PAF 正向坐标约定提取，逐碱基一致。
* wave 初始 trim 越界经几何推演与约 20 万次 fuzz 均不可达，不加防御。
* `spanr fill -n 0` 为 no-op，与设计一致，仅多一次冗余进程。
* LAV `s`/`h` stanza 含空格文件名解析错位，与 UCSC lavToPsl 一致，记录不修。
* wave.rs 的 `unreachable!`/`panic!` 均为算法不变量，有测试兜底。
* tube/greedy 双工作流块边界为 FastGA 管与 greedy 链的语义差异（`--workflow`
  为文档化选项），非缺陷。
* 双引擎差异逐层归因：search 层面 pgi 232 块全部与 lastz 282 块重叠 ≥50%；
  差异传导至 decompose 的 cluster 划分（set_id 与坐标都变）。BISER 语义允许
  两引擎输出不同（可互换替代引擎），各自自洽且坐标正确，非缺陷。
* 对称重复链化歧义（A→B1 与 A→B2 等价）：两引擎均正确，真实 SD 流程的
  cluster 阶段会合并，非缺陷。
* 全量扫描家族生产代码 `unwrap()`/`unreachable!`：`loc.rs` fields[0..2]、
  `cover.rs` f[0..7]、`repeat.rs` fields[0..1] 等索引访问前均有长度检查；
  `merged.last_mut().expect` 有非空前置保证；其余 unwrap/expect 均在测试
  代码。无生产 panic 风险。
* 全部 HashMap 用法逐一核对：pl/repeat 与 trf 的 name/safe_map 仅查找；
  cover 的 by_set 经 set_order（Vec）迭代；decompose 的 index/kmer_frags
  顺序无关；仅 cluster 的 by_root 曾依赖 HashMap 迭代序（复核 60 已修复、
  按首个区间排序）。
* 静态走查：`wave_extend`/`forward_wave`/`extend_chain` 为仅测试使用的导出
  API（生产路径走 `local_alignment`/`psls_from_hits`），无死代码风险。
* `syncmer_dna` 的 `encode_base` 对 N 返回 0（当作 A）与生产路径 N→4 不一致，
  但 syncmer_dna 非生产路径（仅内部测试），不影响 align/rept/sd——记录观察。
* 含空格 contig 名（`>chr one`）：输出键取首个空白 token（"chr"），与
  `fa size` 首字段约定一致（spanr 系既有行为），非缺陷。
* `lav to-psl` 对畸形 LAV 静默输出空：LavReader 跳过空行/注释、未知 stanza 有
  warn，lastz 输出不会畸形，属容错设计，记录不修。

## 记录项（未改，低风险 / 待决策）

* tube 工作流对"库 vs 基因组"的结构性失效：根因是跨对角桶链被切断，结论
  基于修复前代码，syncmer/排序键修复后待真实数据重测。
* `decompose.rs` 负链投影依赖 header 与序列长度一致（cluster 内部保证）。
* cluster/cover 的 u32→i32 坐标转换（仅 >2.1 Gb 染色体溢出）。
* `run_lastz` self 模式仍构建 n×n job 列表，大目录可提前过滤。
* `syncmer.rs` 参考实现与 `collect_one_contig` 重复发射同一位置，消费方
  已 HashSet 去重，可后续合并。
* s_align / sd search --engine pgi 传不支持类型时报错可读性差，不 panic。
* `fa split name` 名称碰撞（`chr(1)` 与 `chr_1`）概率极低，记录不修。
* `rept e-kmer`/`s-kmer` 的 `--fill-kmer` 以 `usize as i32` 传入 `IntSpan::fill`；
  超 i32 值静默截断为负 → fill 变 no-op（无 panic，`excise` 同理安全）。
  极端参数属用户误用，记录不修。
* `s-align` 的 `--min-depth` 以 `usize as u32` 传入深度阈值；超 u32 值截断。
  同上，记录不修。
* `sd cluster` 的同染色体重叠合并只按 chrom 名（物种前缀已剥离）分组：跨
  基因组 PAF 在两端基因组 contig 名与文件 stem 均相同时会把不同基因组的
  同名区间并簇。`sd cluster` 文档仅面向自比对 PAF，记录不修。
* 顶层路径为 `.pgi` 扩展名的目录会被 `is_pgi_input` 误判拒绝（目录名恰好
  以 .pgi 结尾）。概率极低，记录不修。
* `sd search --engine pgi` 接受 `.2bit` 输入（`align pgi` 原生支持），但
  下游 `sd align`/`sd run` 的 chainnet 需要 FASTA，2bit 在 `fa size` 步骤
  报错（外层 run_cmd 只显示失败命令、不含根因）。文档仅承诺 FASTA；2bit
  部分支持是既有行为，记录不修。
* `align lastz --lastz-args` 的值以 `-` 开头时需用 `--lastz-args=<val>` 形式
  （clap 对空格形式的值为标准行为）；帮助文本未提示该写法。
* `--max-gap` 调大（如 10000）时，greedy 循环的 off-band 忽略规则会把后续
  不同对角线的种子全部忽略（"另一条管"），远距重复家族可能整体丢失。该
  行为是 off-band-ignore 的既有设计（默认 1000 下正确），记录不修。
* `rept e-align` 传入 `.2bit` 基因组：在 `has_soft_mask` 的 FASTA 读取器处报
  "stream did not contain valid UTF-8"（二进制被当文本读）。文档仅承诺 FASTA；
  有错误提示的非静默失败，记录不修。
* TRF 外部工具限制：完美 2500 bp 周期 + max-period ≥ 2600、`--max-period`
  10000+ 均触发 TRF SIGSEGV（signal 11），pgr 将信号错误友好传播（无 panic）。
  精确上限未知，pgr 无法可靠预校验，记录不修。
* TRF 版本兼容：本机 TRF 4.09 `-ngs` 输出 17 字段（含末尾 `. .`），
  `parse_trf_output` 的 ≥15 字段门槛兼容；`@chr1` 头行（1 字段）跳过。
* 只有头的 FASTA（`>chr1` 无序列）：`rept s-kmer` 触发 FastK SIGSEGV——
  外部工具对空序列崩溃，cmd_lib 捕获报 "terminated by signal: 11"，pgr 自身
  无 panic，记录不修。
* 纯四联体重复（如 ACGT）只有 4 种不同 10-mer，低于 `MIN_SHARED_KMERS=5`
  防过度分组阈值，同源片段不会合并为同一 set_id——极端低复杂度序列，非 SD
  场景，行为符合设计意图。
* `open_indexed` 的 `.loc` 索引按存在性复用（`force_update=false`）：基因组
  修改后 `.loc` 字节偏移可能过期。复核 104 已加 mtime 新鲜度重建；但修改
  中间基因组本身是用户错误，其余索引（`.paf.idx` 仅显式传入、`.2bit` 显式
  生成）无自动陈旧路径。
* `align lastz -o dir` 重复使用旧 LAV 残留：影响链短（`sd run`/`s-align`/`sd
  search lastz` 均用临时 workdir 免疫），且 LAV 是通用扩展名清理易误伤，
  记录不修。
* s-kmer 尾 run 保守丢弃：Profex `-z` 从不闭合 read 的最后一个 run，s-kmer
  （min_depth=2）按设计保守丢弃尾部（mg1655 尾 run 起点 4641601、约 52 bp，
  低于 min-len 100 会被 excise 过滤，实际影响有限）。行为与 repeat.rs 文档
  "conservatively dropped since its depth is unknown" 一致，记录不修。
* `.pgi` 显式输入 + 冲突 `-k`/`--smer`/`--window` 被静默忽略（exit 0）：
  docs/align-pgi.md 明确说明 "apply only to genome-sequence inputs; .pgi
  inputs carry their parameters in the index header"——文档化预期行为，
  记录不修（sibling 索引路径的冲突报错是额外保护）。

## 已知限制（有意保留）

* 子范围命名 pgr 1-based vs kent 0-based（记录在案）：pgr 生成端/消费端
  自洽，直接消费 UCSC/blat 生态子范围名时需先确认语义。
* s-kmer 对染色体尾部重复保守丢弃：Profex `-z` 不输出末 run 深度，有阈值
  时无法区分唯一尾与重复尾（与 anchr 参考管线一致）。
* 单 contig > 4.3 Gb 的 pgi 索引：pos 为 u32，超长单 contig 不被支持。
* pgi 引擎灵敏度限制：精确 k-mer seed（默认 k=40 + syncmer 8/5）对近
  90–93% identity 或真长恰在 `--min-len` 附近的拷贝可能只锚定子块、边界
  损失 1–11 bp，低于阈值被滤（lastz 用 12-mer seed + 扩展覆盖全长）。属
  引擎语义差异而非逻辑 bug，文档已提示降 `--min-len` 或用 lastz（见
  docs/sd.md）。

## 修复的缺陷（共 81 处：60 处代码/行为 + 21 处 CLI/帮助/文档）

### 崩溃 / 越界 / 溢出（Zero Panic，11 处）

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
**crafted .pgi 记录 contig id 越界 panic**：构造索引的 occurrence 记录携带
    超出 contig 表的 cid 时，`emit_entry_hits` 的 `b.contigs()[bc]` 直接越界
    panic。三个读取路径此前只校验头部不校验记录体。修复：
    `PgiIndex::read` / `PgiStream::next_batch` 逐记录校验 `cid < n_contigs`
    且 `pos + k <= contig len`；`PgiMmap`（惰性解码）在 `emit_entry_hits`
    解码命中记录时同步校验，报友好错误。回归 `crafted_record_contig_rejected_
    not_panic`、`mmap_merge_rejects_out_of_range_contig`、
    `command_align_pgi_crafted_index_errors_not_panics`。

### 内存安全（1 处）

**中段同源检查的 DP band 直接取用户 `--band`，极端组合可 OOM**：
    `middle_is_homologous_range` 的 `align_banded_local` 使用用户 band 原值；
    `--band 10000` + 大中段（≤ 50 kb）时 DP 分配 ~13 GB。修复：DP band 上限
    256（检查只需探测预期对角偏移内同源性，超限保守不合并=碎片化而非丢失）。
    回归验证：`--band 10000 -s 50000` 3.4 s 完成（此前 OOM 风险）。后续
    `chain_windows`/`wave_extend` 的 band 上限统一核对：中段（≤256）、窗口
    扩展（`min(diag_span+32, 128)`）、wave tube（桶内对角带 ≤~128）——无
    未封顶路径。

### 功能正确性 / 算法（31 处，含 3 处重大索引/链算法缺陷）

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
**相邻链合并把两条独立同源对缝成嵌合链，SD 命中丢失**：多拷贝家族中两条
    拷贝对的对角线差可在 band 内（如 56 bp）且两轴间隔均在 merge_gap 内
    （如 3.6 kb），纯几何 merge 会将其缝成一条跨两段真实匹配 + 随机中段的
    嵌合链；扩展出的嵌合块身份 ~72% 被 SD 过滤丢弃。几何条件无法区分
    "同源块种子缺口"与"两条独立块"（形状完全一致），必须用序列判定。
    修复：`merge_adjacent_chains` 增加可选序列参数；两侧间隔均非空时要求
    中段 banded 对齐身份 ≥ 0.9 且 query 覆盖 ≥ 0.9 才合并（单轴插入的 IS
    缝合场景一侧间隔为空，跳过检查）。回归 `merge_requires_homologous_middle_
    with_sequences`、`command_sd_search_pgi_multi_copy_close_diagonals`。
    随机化 6 组端到端 trial 全部拷贝 CORE 覆盖。
**`sd run` 合并 elementary BED 按 `read_dir` 顺序枚举 cluster 文件**：
    set_id 全局重编号与输出行序依赖文件系统枚举顺序，跨运行不确定。修复：
    按 cluster 文件名的**数值**编号排序（词法排序会把 cluster_10 排在
    cluster_2 前）。回归验证：10 家族基因组输出 set_id 严格 1..10 升序。
**`sd align`（`chainnet_to_paf`）按 `read_dir` 顺序迭代 MAF 文件**：
    PAF 输出行序不确定。修复：排序后再合并。
**`sd search --engine lastz` 按 `read_dir` 顺序迭代 LAV 文件**：输出 PSL
    行序不确定。修复：排序后再转换。
**`search_lastz::decompress_if_gz` 解压命名碰撞**：统一命名为 `{base}.plain.
    fa`（base 取首个 `.` 段），嵌套目录中同名 `.fa.gz` 会解压到同一路径，
    后写的静默覆盖先写的，两个 job 比对同一份错误序列。修复：同一次调用内
    用 HashSet 去重，重复 basename 追加输入序号后缀；序号确定性保证 self
    模式下 target/query 两次调用生成相同路径。回归
    `decompress_colliding_basenames_stay_distinct`。
**`psl lift` 的 `parse_subrange` 误切含 `.`/`:` 的 contig 名**：窗口名
    `{contig}:{start}-{end}` 经共享 `Range` 解析器时，`NC_000913.1:1-200`
    被读成 name="NC_000913" + chr="1"、`chr1:alt:1-200` 被读成 chr="alt"，
    `lift_query` 在 sizes 表里查错键、静默跳过提升（仅 warn）。修复：
    `parse_subrange` 改为取最后一个 `:`+数字后缀切分，前缀整体作为 contig
    名；回归 `parse_subrange_keeps_dotted_and_colon_contigs`。
**`s-align` 安全名改写仍用 `split_once(':')`**：coverage.rg 行
    `chr1:alt:1-200` 被切成 name="chr1"，带 `:` 的 contig 输出键被截断成
    "alt"。修复：改用 `parse_subrange` 解析再写占位名。回归
    `command_rept_s_align_colon_name`；`command_rept_s_align_dotted_name`
    断言从"含 '-' 即可"加强为精确坐标（301-800,1001-1500）。
**倒位拷贝间隔 < max_gap 时 greedy 链循环把两条互惠链并成嵌合链，SD 对
    完全漏检**：两条互惠链（同一倒位对的两个方向）落在同一对角线上，拷贝
    间隔小于 max_gap 时其种子在 greedy 循环内直接连成一条链（绕过
    `|diag|>0` 守卫），扩展出横跨两个拷贝 + 间隙的嵌合块，身份被稀释到 SD
    阈值以下被过滤。最小复现：1200 bp 倒位对间隔 800 bp → 修复前 0 条命中、
    一条 3190 bp / 身份 0.879 的嵌合块；修复后 2 条干净命中（1183 bp、100%
    身份）。修复：greedy 循环在"双侧种子间隙 ≥ 200 bp"时用中段同源检查
    门控（不通过则闭合当前链、以该种子起新链）；间隙 < 200 bp 不检查（对
    ≥1000 bp 块，200 bp 随机间隙的嵌合身份 ≥ 0.909，高于 SD 阈值，不会静默
    漏检，且保持稠密种子流廉价）。回归
    `command_sd_search_pgi_close_inverted_repeat`。
**lastz self 模式用 basename 判断自比对，同名文件被交叉比对**：`run_lastz`
    的 self 跳过条件 `t_base != q_base` 只比 basename——目录中含两个同名文件
    （如 `a/dup.fa`、`b/dup.fa`）时，`(a/dup.fa, b/dup.fa)` 会以交叉比对
    方式运行（4 个 LAV 中 2 个虚假交叉），对含共享序列的基因组产生错误命中。
    修复：self 模式跳过所有 `target_file != query_file` 的作业（每个文件只
    与其自身比对）。回归 `command_align_lastz_self_duplicate_basenames`。
**`ref.fa` 与 `ref.fa.gz` 共享兄弟索引，内容不同时静默复用错误索引**：
    `sibling_pgi_path` 的 `set_extension("")` + `set_extension("pgi")` 链把
    `.fa` 替换掉，两文件都映射到 `ref.pgi`；同名同长但序列不同时 contig
    校验（只比名字/长度）无法拦截，第二次运行静默复用第一次的索引（实测 0
    块输出）。修复：`.gz` 输入去掉 `.gz` 后**追加** `.pgi`（`ref.fa.gz` →
    `ref.fa.pgi`），与 `ref.fa` → `ref.pgi` 分离。回归
    `command_align_pgi_gz_sibling_index_distinct`。
**`sd cluster` 的 cluster 编号依赖 HashMap 迭代顺序，`sd run` 的 set_id
    编号跨运行不稳定**：`cluster_paf` 按连通分量分组后直接迭代 `HashMap`
    （进程内随机种子），cluster_N 的编号与文件名对应关系每次运行不同；
    `sd run` 虽按数值排序 cluster 文件，但编号本身随机，导致同一基因组多次
    运行输出 set_id/行序互换（两家族时 r1/r2 互换 set 1/2，实测 5 次运行
    2 种输出）。修复：按每个分组的首个区间（chrom, start）排序后再编号。
    回归 `command_sd_run_output_deterministic_across_runs`。
**FASTA 原地修改后兄弟索引被静默复用**：`resolve_side` 复用同名兄弟 `.pgi`
    时只校验 contig 名/长度；同名单长但序列不同的 FASTA 会静默复用旧索引
    （k-mer 来自旧序列），对齐结果错误。修复：新增 mtime 校验（输入比索引
    新则重建，与 e-kmer 缓存同一约定）。回归
    `command_align_pgi_stale_sibling_index_rebuilt`。
**tube 工作流在显式同文件对（非 --self）时把家族交叉命中当"重复"丢弃**：
    精确自比对巨块（全基因组对角线）在 `dedupe_contained` 中把坐标上包含
    于其内的拷贝对块（两轴 ≥95% 包含）误判为重复并丢弃，显式同文件对只
    输出 1 块自比对。修复：dedupe 增加跨度相近约束（前块跨度 ≤ 后块 4 倍
    才判重复）。回归 `dedupe_keeps_small_block_inside_large_one`；显式同
    文件对 tube 输出 5 块（自比对 + 4 家族命中），tube self 模式 4 块不变。
**`.loc` 索引陈旧时静默使用（open_indexed 只查存在性）**：`loc::open_indexed`
    仅在 `.loc` 不存在（或 force_update）时重建，从不校验新鲜度。实测：1200
    bp FASTA 建 `g.fa.loc` → 改为 1500 bp 后 `sd cluster` 用陈旧索引 →
    `slice error`（长度变化可报错，但**同长度内容修改会静默提取错误序列**）。
    修复：`open_indexed` 增加 mtime 新鲜度校验（`.loc` 的 mtime 早于 FASTA
    时自动重建，`loc_is_fresh`，mtimes 不可用时保持旧行为）。`fa range` /
    `sd cluster` / `fas check` / `get_seq_loc` 四个调用方同步受益。回归
    `stale_loc_index_is_rebuilt`（同长度内容 ACGT→TGCA + mtime 调旧）。
**`.pgi` 单输入自比对 + 仅 `--ref-seq` 报错**：`align pgi ref.pgi --ref-seq
    ref.fa` 报 "extension sequences are needed for both sides"——self 模式下
    query 侧复用 `.pgi` 输入的 `seqs=None`，两侧空/非空不一致触发 bail。
    修复：`resolve_seqs` 后 self 模式下任一侧扩展序列为空时复用另一侧（两
    方向对称）。验证：仅 `--ref-seq` / 仅 `--query-seq` / 双侧 / FASTA 直接
    输入四者输出逐字节一致。回归
    `command_align_pgi_single_ref_seq_on_self_pgi`。
**`sd run --engine lastz` 输出重复 elementary 行**：mg1655 输出 120 行含 5
    个完全重复行（如 607996-609351 出现两次）。追踪：lastz 互反块坐标抖动
    使单个 cluster 文件出现 end 差 1 bp 的两个头，decompose 对两者投影到完全
    相同的 elementary 区间 → 合并后重复。修复：在 `sd run` 合并层按 renumber
    后的完整行去重（`push_unique_elem`），decompose 层按投影坐标去重会破坏
    "相同头不同序列各输出一行"的既有语义。实测：mg1655 lastz run 115 行、
    0 重复（原 120 行含 5 重复），pgi run 不受影响。回归
    `duplicate_elem_rows_are_emitted_once`。
**align pgi 自动索引小写归一化 → 全零块**：构造含大小写混合拷贝的基因组
    （fam 大写 + fam 小写）：修复前输出 match=0/mismatch=0/rep=0 的全零块。
    根因：`build_from_seqs` 的碱基编码大小写不敏感 → 小写与大写拷贝共享 seed
    → 链存在；但扩展 DP 大小写敏感 → 评分失败 → `extend_chain` 回退 raw 块
    （全零）。修复：`build_from_seqs` 增加 `mask` 参数（与 `build_from_path`
    一致），align pgi 自动索引传 `true`（跳过小写）。实测：混合大小写 0 块
    （不再输出全零块）、全大写对照 2 块正常；小写作为软掩码跳过，pgi 与
    lastz 双引擎语义统一。回归
    `command_align_pgi_lowercase_copy_has_no_all_zero_blocks`。
**默认参数静默复用不同 k 的兄弟索引**：`resolve_side` 的缓存参数冲突检查
    只覆盖命令行显式传的 `-k/--smer/--window`（`ValueSource::CommandLine`）。
    实测：`-k 20 --keep-index` 建 k=20 缓存后，`align pgi g.fa`（默认 40）
    静默用 k=20 索引跑 k=40 语义的比对（输出不同，用户无感知）；显式 `-k
    40` 则报错——两条路径行为不一致。修复：删除 `explicit(...)` 条件，**总是**
    检查当前解析值（显式或默认）与缓存索引参数的一致性（smer/window 对称
    生效）。回归 `command_align_pgi_default_kmer_conflicts_with_cached_index`。

### 输入校验 / 静默错误（6 处）

**repeat.rs 两处 `map_while(Result::ok)` 吞 IO 错误**。修复：
    `let line = line?;` 传播错误。
**e-align PSL 过滤静默跳过畸形行**。修复：补 `log::warn!`。
**decompose 对解析失败的 FASTA 头静默丢弃**。修复：补 `log::warn!`。
**`sd search`/`sd cross`/`sd run` 传入 `.pgi` 索引**：pgi 引擎对 `.pgi`
    输入不做扩展（无序列），输出块全部 0 分，SD 过滤后静默返回空结果。修复：
    `pgi_to_hits`/`lastz_to_hits` 前置拒绝 `.pgi` 输入（magic 或扩展名），报
    友好错误。回归 `command_sd_search_rejects_pgi_input`。
**`sd align` 跳过非 2 组件的 MAF 块时无提示**：`maf_block_to_paf` 对 <2 /
    >2 组件的块返回 None（注释称 "caller logs warning"），但
    `chainnet_to_paf` 未记日志。补 `log::warn!`（chainnet 输出恒为 2 组件，
    路径为防御性提示）。
**空 FASTA 输入触发 FastK SIGSEGV（预检友好报错）**：空 repeat 库（`>empty1`
    无序列）喂 `rept e-kmer` → FastK SIGSEGV，pgr 报 "terminated by signal:
    11"（错误信息像 pgr 自身崩溃，不友好）。修复：`run_repeat_pipeline` 在
    FastK 前预检输入是否有非空序列（`has_sequences`），空则报友好错误；全
    N/4 bp 极小库仍走 FastK（工具 exit 1 可接受）。回归
    `sequence_less_fasta_is_detected`。

### 数据安全（`-o` 同输入保护 / 陈旧索引 / 静默数据丢失，4 处）

**sd 命令 `-o` 指向输入文件时静默覆盖输入**：`sd search g.fa -o g.fa` 等把
    输入 FASTA/PAF/BED 覆盖为变换后的输出（exit 0、无提示）。修复：`sd
    search`/`align`/`cover`/`decompose`/`cross` 均加 `ensure_outfile_distinct`
    检查。回归 `command_sd_output_same_as_input_rejected`。
**rept 与 align pgi 同样存在 `-o` 覆盖输入**：`rept s-kmer g.fa -o g.fa`、
    `rept trf`、`align pgi g.fa -o g.fa` 等把输入覆盖为 runlist JSON/PSL。
    修复：rept 五个子命令（含库输入）与 `align pgi`（含 `--ref-seq`/
    `--query-seq`）均加 `ensure_outfile_distinct`。回归
    `command_rept_output_same_as_input_rejected`。
**损坏的 FastK 缓存被静默复用 → e-kmer 空输出**：`cache_is_fresh` 只检查缓存
    存在 + mtime。实测将 `lib.fa.gz.repeat.k17.ktab` 截断为 100 字节后（
    `.complete` 标记和 part 文件完好）：日志显示 "reused repeat table"，
    FastK 静默读取损坏表 → e-kmer 输出空 runlist（mg1655 原 48 区间全部丢失），
    比报错更隐蔽。修复：`cache_is_fresh` 增加 `.ktab` 与 `.complete` 大小一致
    性校验，不一致即视为陈旧重建。实测重建后输出与原 48 区间逐字节一致。
    回归 `truncated_cache_table_is_not_fresh`。
**`sd cluster` 输出目录残留旧 cluster 文件**：向含 `cluster_1.fa`/
    `cluster_2.fa`/`cluster_3.fa` 的目录重跑 `sd cluster`（本次仅 1 个
    cluster）：旧 `cluster_2.fa`/`cluster_3.fa` 残留，下游 `sd decompose` /
    手动 `sd run` 会**静默消费陈旧家族**（`sd run` 内部用固定 tempdir 免疫，
    但手动工作流受影响）。修复：写输出前清理 outdir 中 pgr 自身命名模式的
    `cluster_<u32>.fa`（仅此模式）。回归 `stale_cluster_files_are_removed`。

### 性能（1 处）

**`align pgi --parallel` 未约束自动索引构建的 rayon 并行度**：`resolve_side`
    （内部 `build_from_seqs` → `radix_sort_u128_par`）在自定义线程池创建前
    执行，索引构建走全局 rayon 池，`--parallel N` 只约束 merge/扩展阶段。
    文档承诺 "--parallel: rayon thread count"，行为不一致。修复：把从
    `resolve_side`（索引构建）到 merge/扩展的整个流程移入 `pool.install`，
    `--parallel` 现约束整个命令的 rayon 用量（`sd search --engine pgi` 与
    `rept e-align` 经由 `align pgi` 同步受益）。`-p 1/2/8` 输出逐字节一致
    （确定性未破坏）。

### 外部工具与参数 / CLI（6 处）

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

### CLI / 文档（21 处）

**噪音与帮助文本多处小修**：lav mask stanza 静默、`#` 元数据行跳过、lastz
    `[multiple]`/`-s` 修正、align.md 示例输出修正、pgi 帮助默认 syncmer 修正。
**文档一致性**：rept.md 补 e-align；soft-mask 说明；`.pgi` 命名；sd.md lastz
    单序列/纯文本约束；TnCentral 路径。
**lastz 单序列约束帮助/文档未同步**。修复：四处补齐。
**e-align identity 定义未说明**（gap-compressed）。修复：补文档。
**主帮助 rept 子命令列表漏 e-align/s-align**。修复：补齐。
**rept.md 仍写 "`align` variants are planned"**（e-align/s-align 早已存在）
    → 改为现况描述。
**rept.md "All four emit runlist JSON" → "All five"**（trf 在内共 5 命令）。
**rept.md e-align 空 "### Dependencies" 章节删除**。
**align-pgi.md `--freq` 语义错误**：写 "more than this many times"，代码与
    帮助均为 "at least this many times"（`>= freq`）→ 文档修正。
**`sd run` 帮助/文档补齐**：`--preset` 默认 set01、`--min-identity (0, 1]`、
    lastz 引擎需单序列 FASTA。
**`pgr align` 的 about 写 "into PSL blocks"，但 lastz 子命令输出 LAV**。
    改为 "Aligns genomes or .pgi indexes"。
**`align pgi` 兄弟索引命名描述错误**：`sibling_pgi_path` 用 `set_extension`
    替换最终扩展名，docs/align-pgi.md 却写 "ref.fa.gz → ref.fa.pgi"。修复
    文档为"最终扩展名替换为 .pgi，ref.fa 与 ref.fa.gz 均映射到 ref.pgi"
    （复核 51 后改为"追加 .pgi"的分离命名，文档恢复原意）。
**align-pgi.md `--merge-gap` 说明补齐序列校验语义**：两侧间隔均非空时合并
    还要求中段同源（banded 对齐验证），近对角线的独立拷贝对保持分离。
**align-pgi.md `--max-gap` 说明补齐 greedy 门控语义**：双侧间隙 ≥ 200 bp 时
    仅同源中段才桥接，近距离倒位对保持分离。
**align-pgi.md 兄弟索引说明补齐 mtime 失效约定**。
**align-pgi.md 补充 sibling 缓存索引参数一致性**：当前 `-k/--smer/--window`
    （显式或默认）必须与缓存匹配，不匹配报错而非静默复用不同 seed（缺陷 36）。
**sd.md 补充 `sd run` 输出去重**：近相同 cluster 拷贝（互反块 1 bp 抖动）投影
    到相同 elementary 区间时只输出一次（缺陷 34）。
**sd.md search 节补充 pgi 引擎灵敏度限制**：精确 k-mer seed 对近 90–93%
    identity 拷贝可能只锚定子块，低于 `--min-len` 被滤；提示降 `--min-len` 或
    用 `--engine lastz`（复核 99）。高 identity 且真长恰在 min-len 附近的拷贝
    也可能因 seed 边界损失差几 bp 被滤（复核 121）。
**align-pgi.md 明确 `--ref-seq` 校验范围**（contig 表）并要求序列与索引来源
    一致（自动 sibling 路径由 mtime 检查保证）；未实现 k-mer 内容校验
    （syncmer 哈希对比复杂度高、阈值易误报，文档说明足够）。
**align-pgi.md Notes 补充小写（软掩码）处理**：自动索引小写→N 无 seed/块，
    `pgr pgi build --mask` 同语义。
**sd.md Notes 补充软掩码语义**：pgi（`-M` 语义）与 lastz（小写视为掩码）
    都不比对小写，软掩码的 SD 拷贝不被检出，建议先 `tr a-z A-Z`。

## 验证

* 引擎交叉验证：pgi 与 lastz 检出同一对 1200 bp 拷贝（边界修剪 2 bp）；
  倒位重复经 chainnet 后两条 `-` 链保留，`sd run` 输出坐标正确的 elementary
  SD；合成基因组上两引擎各 4 条命中覆盖相同两个重复家族（坐标差异仅边界
  修剪 4–8 bp）。
* 端到端坐标：`rept trf` 输出 "101-1100"、`rept s-align` 输出 "501-2900"
  （1-based 全覆盖）；多 contig s-kmer 编号与 chr.sizes 一致；e-kmer
  `--keep-index` 缓存复用/失效正确；普通 gz 与 BGZF 全流程通过；`.2bit`
  输入与 FASTA 路径逐字节一致。
* 鲁棒性：截断/负跨度 lav、越界 PAF、垃圾 BED/.pgi、构造头、空输入、全 N、
  极值参数、短行、随机二进制喂各命令等畸形输入全部友好报错或空输出，零
  panic；`sd`/`rept`/`align` 空输入全链路（search 0 块 → align 空 PAF →
  cluster 空目录 → run 空 BED，exit 0）通过。
* 确定性：`sd run`（8 次 1 种哈希、10 家族 diff 为空）、`sd search --engine
  lastz`/`sd align`/`sd cross`/`rept s-align`/`trf`/`s-kmer` 多次运行输出
  逐字节一致；`--parallel 1/2/8` 下输出逐字节一致；tube workflow 多次运行
  逐字节一致。
* 数据安全：`-o` 同输入（sd/rept/align pgi）均报 "also an input file" 且输入
  完好；`.loc` 陈旧重建、`.pgi` sibling mtime 重建、FastK 截断缓存重建、
  `sd cluster` 旧文件清理均实测复现；`sd run` 完成后 `/tmp/pgr_sd_*` 临时
  文件数为 0，无文件泄漏。
* 真实数据：MG1655 `sd search --engine pgi` 229→232 条命中（复核 21 后 8 条
  新增为嵌合块拆分、无覆盖丢失）、136 个拷贝对；`sd run` pgi 118 行 / lastz
  115 行（去重后）；`rept` 五命令 e-kmer 48 区间、e-align 89、s-kmer 170、
  s-align 1457、trf 84——其中 e-kmer 与 trf 与 docs 记录值逐字节一致。
* 性能：5 Mb 随机基因组 + 双拷贝 `sd search` 6.8 s、`sd run` 全链路 7.7 s；
  `align pgi --band 10000 -s 50000` 3.4 s（修复前 OOM 风险）；MG1655 `sd run`
  debug 41 s。
* 新增回归测试（主要）：`command_sd_search_pgi_inverted_repeat`、
  `command_sd_search_pgi_close_inverted_repeat`、
  `command_sd_search_pgi_multi_copy_close_diagonals`、
  `command_sd_run_output_deterministic_across_runs`、
  `command_sd_output_same_as_input_rejected`、
  `command_rept_output_same_as_input_rejected`、
  `command_align_pgi_crafted_index_errors_not_panics`、
  `command_align_pgi_gz_sibling_index_distinct`、
  `command_align_pgi_stale_sibling_index_rebuilt`、
  `command_align_pgi_default_kmer_conflicts_with_cached_index`、
  `command_align_pgi_lowercase_copy_has_no_all_zero_blocks`、
  `command_align_lastz_self_duplicate_basenames`、
  `stale_loc_index_is_rebuilt`、`truncated_cache_table_is_not_fresh`、
  `stale_cluster_files_are_removed`、`sequence_less_fasta_is_detected`、
  `duplicate_elem_rows_are_emitted_once`、
  `merge_checks_minus_strand_middle_in_rc_space`、
  `dedupe_keeps_small_block_inside_large_one`、
  `randomized_single_pass_matches_reference` 等。
* `cargo test` 全量 1255 通过（历轮 995→1201→1209→1210→1211→1212→1213→
  1218→1219→1220→1223→1224→1226→1227→1229→1230→1231→1232→1234→1236→
  1239→1240→1241→1243→1249→1250→1253→1254→1255 递增）；本族 release 模式
  全绿（pgi 52、sd 13、pl 1、alignment 23 个 lib + 44 个 CLI）；`cargo fmt
  --check` 与 `cargo clippy --all-targets -- -D warnings` 干净。

## 结论

`sd`/`rept`/`align` 三个命令族审核完成（累计修复 81 处缺陷：60 处代码/行为 +
21 处 CLI/帮助/文档），并经多轮纵深复核（`libs/sd`、`libs/pgi`、`libs/lastz`、
`libs/fmt/lav`、`libs/fmt/psl`、`libs/pl/repeat`、`libs/alignment` DP 与全部
命令执行路径、tube/greedy 双工作流、索引/缓存新鲜度、确定性、`-o` 覆盖保护、
HashMap 迭代序、外部工具封装）复核，未再发现新问题，审核收敛。