# sd / rept 命令族代码审核记录（2026-08-05）

对 `pgr sd`（8 命令 + `libs/sd`）与 `pgr rept`（6 命令 + `libs/pl`）两个命令族
约 5000 行代码及全部文档（`docs/{sd,rept}.md`）进行审核。两个命令族结构相近
（同为重复序列分析工作流，共享 `libs/pgi` 索引消费端与 `libs/loc` 索引），故
合并为一份审核记录。缺陷按类别分组记录；关键修复均附回归测试，验证概况见文末
"验证"一节。

> 注：`pgr align` 命令族的审核记录见 `notes/audit/audit-pgi-align.md`。`libs/pgi`
> 索引的**构建**缺陷（k-mer key、构造头、记录越界、sibling 索引、`--parallel`
> 等）记录在 audit-pgi-align.md；本文件记录 sd 对 pgi 的**消费**缺陷（`sd search`
> 传 `.pgi` 拒绝、pgi merge 频率过滤、greedy/tube 链逻辑等）。

审核范围：
- **sd**：`search` / `cross` / `align` / `cluster` / `decompose` / `cover` /
  `run`（+ `libs/sd`、tube/greedy 双工作流、pgi 消费端）
- **rept**：`e-kmer` / `s-kmer` / `e-align` / `s-align` / `trf`
  （+ `libs/pl/repeat`、FastK/Profex/TRF 外部工具封装）

审核重点：Zero Panic（畸形输入不 panic）、数据安全（`-o` 不得覆盖输入、陈旧/
损坏索引不得静默复用）、确定性（跨运行输出逐字节一致）、与外部参考实现
（FastGA / UCSC kent / lastz / TRF）的语义一致。

## 与外部参考实现的语义一致性核对

关键修复均对照官方源码复核，方向一致：

* greedy 链合并：与 FastGA `align_contigs` / `ALNchain.c` 的链化
  语义一致——同对角线纯间隔是两条独立链，仅对角线平移才缝合（pgr 的
  自身扩展）。
* pgi 索引 k-mer 频率过滤：`emit_entry_hits` 的 canonical key 过滤、
  `freq >= cutoff`（FastGA 语义，非 `>`）、前缀窗口 / 最大共享前缀 /
  扩展范围过滤，与 FastGA GIX 语义一致。
* 软掩码语义：`build_from_seqs` 新增 `mask`（小写→N，FastGA `-M` 语义），
  与 `build_from_path` 一致；pgi 与 lastz 对小写拷贝行为完全一致。
* `sd run` 的 chainnet：每靶位点保留一条最优链（同等分按序取一）是 UCSC
  chainnet 每靶位点取最佳链的标准语义（与 `pl chainnet` 共享）。

## 排除的疑点（经核验无需修复）

* `sd run` 的 cluster set_id 重编号值域各簇两两不相交，不可能碰撞。
* 60,423 → 75,413 数据差异来自 tncentral 库更新与编译时序，非代码 bug
  （repeat-masking.md §2.3.5 已勘误）。
* sd cluster minus 链序列提取：按 pgr PAF 正向坐标约定提取，逐碱基一致。
* wave 初始 trim 越界经几何推演与约 20 万次 fuzz 均不可达，不加防御。
* `spanr fill -n 0` 为 no-op，与设计一致，仅多一次冗余进程。
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
* 含空格 contig 名（`>chr one`）：输出键取首个空白 token（"chr"），与
  `fa size` 首字段约定一致（spanr 系既有行为），非缺陷。
* `syncmer_dna` 的 `encode_base` 对 N 返回 0（当作 A）与生产路径 N→4 不一致，
  但 syncmer_dna 非生产路径（仅内部测试），不影响 align/rept/sd——记录观察。
* tube 工作流"库 vs 基因组"结构性失效（重测确认**无需修复**）：
  原疑"跨对角桶链被切断"；greedy 已移除、tube 为唯一流程后，MG1655 vs
  TnCentral 库 `rept e-align` 正常检出（71.6 kb，79% 与 e-kmer 重叠），
  失效随 syncmer/排序键修复消失。

## 记录项（未改，低风险 / 待决策）

* `decompose.rs` 负链投影依赖 header 与序列长度一致（cluster 内部保证）。
* cluster/cover 的 u32→i32 坐标转换（仅 >2.1 Gb 染色体溢出）。
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
* `save_repeat_cache` 中途失败（`.ktab` 已写、部分 part 文件已写、`.complete`
  未写）时，残留的 `.ktab`/part 文件无 `.complete` 标记，`cache_is_fresh`
  判陈旧会自动重建，不会复用；但若后续一次**成功**重建后 part 文件数变少，
  旧的高序号 part 文件（`.<base>.ktab.N`）可能残留被 FastK `-p:` 读到。
  需"保存中途失败 + 后续重建 part 数减少"同时成立，概率极低，记录不修。
* `sd search --engine pgi` 接受 `.2bit` 输入（`align pgi` 原生支持），但
  下游 `sd align`/`sd run` 的 chainnet 需要 FASTA，2bit 在 `fa size` 步骤
  报错（外层 run_cmd 只显示失败命令、不含根因）。文档仅承诺 FASTA；2bit
  部分支持是既有行为，记录不修。
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
* s-kmer 尾 run 保守丢弃：Profex `-z` 从不闭合 read 的最后一个 run，s-kmer
  （min_depth=2）按设计保守丢弃尾部（mg1655 尾 run 起点 4641601、约 52 bp，
  低于 min-len 100 会被 excise 过滤，实际影响有限）。行为与 repeat.rs 文档
  "conservatively dropped since its depth is unknown" 一致，记录不修。

## 已知限制（有意保留）

* s-kmer 对染色体尾部重复保守丢弃：Profex `-z` 不输出末 run 深度，有阈值
  时无法区分唯一尾与重复尾（与 anchr 参考管线一致）。

## 修复的缺陷（共 51 处：41 处代码/行为 + 10 处 CLI/帮助/文档）

### 崩溃 / 越界 / 溢出（Zero Panic，4 处）

**sd/run.rs 解析 elem.bed 短行越界**：直接取 `f[4]`。修复：加
   `f.len() < 8` 检查（与 cover.rs 一致）。
**sd decompose 负链投影 usize 下溢**（畸形 header）。修复：拒绝
   end < start，投影 saturating。回归 `malformed_header_does_not_panic`。
**非 UTF-8 临时目录路径 `to_str().unwrap()` panic**（sd run 临时目录）。
   修复：`io::path_to_str` 友好报错。
**e-align span 过滤 `(t_end - t_start) as usize` 回绕**。修复：i64
   运算 `.max(0)` 再转 usize。

### 功能正确性 / 算法（22 处，含 2 处重大链算法缺陷）

**（重大）tube 排序键 anti/bucket 溢出**（>8 Mb 基因组失效）。修复：
   anti/bucket 扩到 32 位。回归
   `tube_sort_key_supports_large_anti_coordinates`。
**（重大）tube 排序键负对角线回绕**（>64 Mb 间距失效）。修复：
   `BUCK_OFF = 1 << 26`。回归深负对角线两个测试。

**pgi 引擎灵敏度限制**（记录项升级）：精确 k-mer seed 对近
  90–93% identity 或真长恰在 `--min-len` 附近的拷贝可能只锚定子块被滤。
  已解决：`sd search` 默认 `freq=50/k=31` 后，10 个 E. coli 漏检率
  13.1%→0.26%，遮蔽流程后 pgi/lastz 互相漏检 3.2%/6.0%（详见
  `design/sd.md` §4.9/§4.10）。
**cluster 重叠 union 漏连嵌套区间**。修复：扫描时跟踪最大右端。回归
    `nested_overlapping_intervals_form_one_cluster`。
**sd cluster 去重键忽略链向/物种**（回文倒位拷贝被折叠）。修复：键加
    strand。回归 `same_coordinates_on_opposite_strands_are_distinct_copies`。
**sd cluster/run 不支持普通 gzip**（生成垃圾 `.loc`）。修复：非 BGZF
    先解压到临时文件。回归 `command_sd_run_gzipped_genome`。
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
**`sd cluster` 的 cluster 编号依赖 HashMap 迭代顺序，`sd run` 的 set_id
    编号跨运行不稳定**：`cluster_paf` 按连通分量分组后直接迭代 `HashMap`
    （进程内随机种子），cluster_N 的编号与文件名对应关系每次运行不同；
    `sd run` 虽按数值排序 cluster 文件，但编号本身随机，导致同一基因组多次
    运行输出 set_id/行序互换（两家族时 r1/r2 互换 set 1/2，实测 5 次运行
    2 种输出）。修复：按每个分组的首个区间（chrom, start）排序后再编号。
    回归 `command_sd_run_output_deterministic_across_runs`。
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
**`sd run --engine lastz` 输出重复 elementary 行**：mg1655 输出 120 行含 5
    个完全重复行（如 607996-609351 出现两次）。追踪：lastz 互反块坐标抖动
    使单个 cluster 文件出现 end 差 1 bp 的两个头，decompose 对两者投影到完全
    相同的 elementary 区间 → 合并后重复。修复：在 `sd run` 合并层按 renumber
    后的完整行去重（`push_unique_elem`），decompose 层按投影坐标去重会破坏
    "相同头不同序列各输出一行"的既有语义。实测：mg1655 lastz run 115 行、
    0 重复（原 120 行含 5 重复），pgi run 不受影响。回归
    `duplicate_elem_rows_are_emitted_once`。
**s-align 漏做带点 contig 名映射**（spanr 截断，`fa mask` 失配）。
    修复：复用 chrom.sizes 映射。回归 `command_rept_s_align_dotted_name`。
**Profex `-z` 坐标右端多 +1 + e-kmer 染色体尾部丢失**。修复：end 不再
    +1；无阈值时用染色体长度闭合尾 run。回归
    `command_rept_e_kmer_tandem_coordinates`。
**s-align/e-align soft-mask 警告误报 N gap**。修复：`has_soft_mask`
    只扫 lowercase。回归 `soft_mask_detection_ignores_n_gaps`。
**`s-align` 安全名改写仍用 `split_once(':')`**：coverage.rg 行
    `chr1:alt:1-200` 被切成 name="chr1"，带 `:` 的 contig 输出键被截断成
    "alt"。修复：改用 `parse_subrange` 解析再写占位名。回归
    `command_rept_s_align_colon_name`；`command_rept_s_align_dotted_name`
    断言从"含 '-' 即可"加强为精确坐标（301-800,1001-1500）。

### 输入校验 / 静默错误（6 处）

**decompose 对解析失败的 FASTA 头静默丢弃**。修复：补 `log::warn!`。
**`sd search`/`sd cross`/`sd run` 传入 `.pgi` 索引**：pgi 引擎对 `.pgi`
    输入不做扩展（无序列），输出块全部 0 分，SD 过滤后静默返回空结果。修复：
    `pgi_to_hits`/`lastz_to_hits` 前置拒绝 `.pgi` 输入（magic 或扩展名），报
    友好错误。回归 `command_sd_search_rejects_pgi_input`。
**`sd align` 跳过非 2 组件的 MAF 块时无提示**：`maf_block_to_paf` 对 <2 /
    >2 组件的块返回 None（注释称 "caller logs warning"），但
    `chainnet_to_paf` 未记日志。补 `log::warn!`（chainnet 输出恒为 2 组件，
    路径为防御性提示）。
**repeat.rs 两处 `map_while(Result::ok)` 吞 IO 错误**。修复：
    `let line = line?;` 传播错误。
**e-align PSL 过滤静默跳过畸形行**。修复：补 `log::warn!`。
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
**`sd cluster` 输出目录残留旧 cluster 文件**：向含 `cluster_1.fa`/
    `cluster_2.fa`/`cluster_3.fa` 的目录重跑 `sd cluster`（本次仅 1 个
    cluster）：旧 `cluster_2.fa`/`cluster_3.fa` 残留，下游 `sd decompose` /
    手动 `sd run` 会**静默消费陈旧家族**（`sd run` 内部用固定 tempdir 免疫，
    但手动工作流受影响）。修复：写输出前清理 outdir 中 pgr 自身命名模式的
    `cluster_<u32>.fa`（仅此模式）。回归 `stale_cluster_files_are_removed`。
**rept 命令 `-o` 指向输入文件时静默覆盖输入**：`rept s-kmer g.fa -o g.fa`、
    `rept trf` 等把输入覆盖为 runlist JSON。修复：rept 五个子命令（含库输入）
    均加 `ensure_outfile_distinct`。回归
    `command_rept_output_same_as_input_rejected`。
**损坏的 FastK 缓存被静默复用 → e-kmer 空输出**：`cache_is_fresh` 只检查缓存
    存在 + mtime。实测将 `lib.fa.gz.repeat.k17.ktab` 截断为 100 字节后（
    `.complete` 标记和 part 文件完好）：日志显示 "reused repeat table"，
    FastK 静默读取损坏表 → e-kmer 输出空 runlist（mg1655 原 48 区间全部丢失），
    比报错更隐蔽。修复：`cache_is_fresh` 增加 `.ktab` 与 `.complete` 大小一致
    性校验，不一致即视为陈旧重建。实测重建后输出与原 48 区间逐字节一致。
    回归 `truncated_cache_table_is_not_fresh`。

### 性能（1 处）

**`run_lastz` self 模式 n×n job 列表**（记录项移入）：self 模式
  只构建对角 n 个 job（不再生成 n²），执行期防御保留。

### 外部工具与参数 / CLI（4 处）

**参数校验缺失/不一致（sd 侧）**：`--min-identity` 范围、minscore 正值
    有限性。修复：统一校验，帮助同步 "(0, 1]"。
**sd search/cross `--preset` 默认值未注册**。修复：
    `.default_value("set01")`。回归 `command_sd_search_lastz_default_preset`。
**sd run --engine lastz --preset 拼装错误**。修复：`Vec<String>` +
    `$[preset_args]` 展开。回归 `command_sd_run_lastz_preset_parses`。
**trf 特殊字符文件名找不到**。修复：`sanitize_filename(chr)`。回归
    `command_rept_trf_special_chars`。

### CLI / 文档（10 处）

**文档一致性（rept 侧）**：rept.md 补 e-align；sd.md lastz 单序列/纯文本
    约束；TnCentral 路径。
**e-align identity 定义未说明**（gap-compressed）。修复：补文档。
**主帮助 rept 子命令列表漏 e-align/s-align**。修复：补齐。
**rept.md 仍写 "`align` variants are planned"**（e-align/s-align 早已存在）
    → 改为现况描述。
**rept.md "All four emit runlist JSON" → "All five"**（trf 在内共 5 命令）。
**rept.md e-align 空 "### Dependencies" 章节删除**。
**`sd run` 帮助/文档补齐**：`--preset` 默认 set01、`--min-identity (0, 1]`、
    lastz 引擎需单序列 FASTA。
**sd.md 补充 `sd run` 输出去重**：近相同 cluster 拷贝（互反块 1 bp 抖动）投影
    到相同 elementary 区间时只输出一次（缺陷 34）。
**sd.md search 节补充 pgi 引擎灵敏度限制**：精确 k-mer seed 对近 90–93%
    identity 拷贝可能只锚定子块，低于 `--min-len` 被滤；提示降 `--min-len` 或
    用 `--engine lastz`（复核 99）。高 identity 且真长恰在 min-len 附近的拷贝
    也可能因 seed 边界损失差几 bp 被滤（复核 121）。
**sd.md Notes 补充软掩码语义**：pgi（`-M` 语义）与 lastz（小写视为掩码）
    都不比对小写，软掩码的 SD 拷贝不被检出，建议先 `tr a-z A-Z`。
    > 修正：实测 lastz 大小写不敏感（小写仍参与匹配），仅 pgi
    感知小写；且 mask 仅影响 `sd search` 发现阶段，后续阶段读原始基因组。

## 验证

* 引擎交叉验证：pgi 与 lastz 检出同一对 1200 bp 拷贝（边界修剪 2 bp）；
  倒位重复经 chainnet 后两条 `-` 链保留，`sd run` 输出坐标正确的 elementary
  SD；合成基因组上两引擎各 4 条命中覆盖相同两个重复家族（坐标差异仅边界
  修剪 4–8 bp）。
* 端到端坐标：`rept trf` 输出 "101-1100"、`rept s-align` 输出 "501-2900"
  （1-based 全覆盖）；多 contig s-kmer 编号与 chr.sizes 一致；e-kmer
  `--keep-index` 缓存复用/失效正确；普通 gz 与 BGZF 全流程通过。
* 鲁棒性：截断/负跨度 lav、越界 PAF、垃圾 BED/.pgi、构造头、空输入、全 N、
  极值参数、短行、随机二进制喂各命令等畸形输入全部友好报错或空输出，零
  panic；`sd`/`rept` 空输入全链路（search 0 块 → align 空 PAF → cluster 空
  目录 → run 空 BED，exit 0）通过。
* 确定性：`sd run`（8 次 1 种哈希、10 家族 diff 为空）、`sd search --engine
  lastz`/`sd align`/`sd cross`/`rept s-align`/`trf`/`s-kmer` 多次运行输出
  逐字节一致；tube workflow 多次运行逐字节一致。
* 数据安全：`-o` 同输入（sd/rept）均报 "also an input file" 且输入完好；
  `.loc` 陈旧重建、FastK 截断缓存重建、`sd cluster` 旧文件清理均实测复现；
  `sd run` 完成后 `/tmp/pgr_sd_*` 临时文件数为 0，无文件泄漏。
* 真实数据：MG1655 `sd search --engine pgi` 229→232 条命中（复核 21 后 8 条
  新增为嵌合块拆分、无覆盖丢失）、136 个拷贝对；`sd run` pgi 118 行 / lastz
  115 行（去重后）；`rept` 五命令 e-kmer 48 区间、e-align 89、s-kmer 170、
  s-align 1457、trf 84——其中 e-kmer 与 trf 与 docs 记录值逐字节一致。
* 性能：5 Mb 随机基因组 + 双拷贝 `sd search` 6.8 s、`sd run` 全链路 7.7 s；
  MG1655 `sd run` debug 41 s。
* 新增回归测试（sd）：`command_sd_search_pgi_inverted_repeat`、
  `command_sd_search_pgi_close_inverted_repeat`、
  `command_sd_search_pgi_multi_copy_close_diagonals`、
  `command_sd_run_output_deterministic_across_runs`、
  `command_sd_output_same_as_input_rejected`、`stale_cluster_files_are_removed`、
  `duplicate_elem_rows_are_emitted_once`、`stale_loc_index_is_rebuilt`、
  `dedupe_keeps_small_block_inside_large_one`、
  `merge_requires_homologous_middle_with_sequences`、
  `command_sd_run_gzipped_genome`、`command_sd_search_rejects_pgi_input` 等。
* 新增回归测试（rept）：`command_rept_output_same_as_input_rejected`、
  `command_rept_s_align_colon_name`、`command_rept_s_align_dotted_name`、
  `command_rept_e_kmer_tandem_coordinates`、`command_rept_trf_special_chars`、
  `truncated_cache_table_is_not_fresh`、`sequence_less_fasta_is_detected`、
  `soft_mask_detection_ignores_n_gaps` 等。
* 共享库回归：`randomized_single_pass_matches_reference`、
  `merge_checks_minus_strand_middle_in_rc_space` 等。
* `cargo test` 全量 1255 通过（逐轮递增至收敛）；本族 sd 13、pl 1
  个 lib + 相关 CLI 测试 release 模式全绿；`cargo fmt --check` 与 `cargo
  clippy --all-targets -- -D warnings` 干净。
* 复查（cross 解压 + e-align 文档修复后）：`cargo test --lib` 568 通过、
  `cli_sd` 15 / `cli_rept` 17 通过、clippy 干净。新增的回归测试
  `decompress_colliding_basenames_stay_distinct` 触发 clippy
  `cloned_ref_to_slice_refs`（`&[a.clone()]`），已改为 `std::slice::from_ref(&a)`
  消除（测试代码，非生产逻辑）。复查 sd/rept 六命令 lib 与 repeat 管线
  （runlist JSON 恢复、`-` 空标记丢弃、safe 名双射、DSU 容量、`set_id` 全局
  重编号）未再发现新问题。

## 结论

`sd`/`rept` 两个命令族审核完成（累计修复 51 处缺陷：41 处代码/行为 + 10 处
CLI/帮助/文档），并经多轮纵深复核（`libs/sd`、`libs/pl/repeat`、tube/greedy
双工作流、索引/缓存新鲜度、确定性、`-o` 覆盖保护、HashMap 迭代序、外部工具
封装）复核，未再发现新问题，审核收敛。pgi 索引构建侧缺陷见 audit-pgi-align.md。
