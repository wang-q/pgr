# sd / rept / align 命令族代码审核记录（2026-08-04）

对新增命令族 sd / rept / align 的代码与文档进行多轮深入审核。范围：sd（8
命令 + libs/sd）、rept（6 命令 + libs/pl）、align（pgi/lastz + libs/pgi、
libs/lastz、libs/fmt/lav、alignment DP），约 8000 行代码逐文件审完；文档
docs/{sd,rept,align-pgi,align-lastz}.md 全部核对。

轮次：第一轮 #1-25；连续两轮无新问题后，应要求追加第二轮 #26-34（UCSC
PSL 负链约定、greedy 链合并等），对照 UCSC kent 与 FastGA 官方源码复核
（#26a/26b/31a）确认修复方向；第三轮 #35-37；第四轮 #37 精修；第五/六轮
确认。每轮发现问题后修复并进入下一轮复核，连续两轮无新问题后收束。最终
995 测试全绿、`cargo fmt` / `cargo clippy --all-targets -- -D warnings`
干净。

## 修复的缺陷（共 37 处）

修复按发现顺序编号 #1-37，按类别分组；#26a/26b/31a 为外部源码复核，
见「与外部参考实现的语义一致性核对」。

### 崩溃 / 越界 / 溢出（Zero Panic，10 处）

4. **sd/run.rs 解析 elem.bed 越界 panic**：短行直接取 `f[4]`。加
   `f.len() < 8` 检查（与 cover.rs 一致）。
5. **sd decompose 负链投影 usize 下溢 panic**：畸形 header（start > end 或
   跨度小于序列长度）使 `gend - end` 下溢，debug 构建直接崩溃。修复：
   `parse_header` 拒绝 end < start；投影改用 saturating 运算。回归测试
   `malformed_header_does_not_panic`。
6. **lav d stanza 边界检查差一 → 越界 panic**：守卫
   `header_idx + 5 > lines.len()` 对"matrix 头 + 4 行矩阵而无参数行"的场景
   放行，随后 `lines[header_idx + 5]` 越界。修复为 `+ 6`。回归测试
   `truncated_d_stanza_errors_not_panics`。
7. **构造 .pgi/.hv 头容量溢出 panic / OOM**：`PgiIndex::read` 用未校验的
   `n_records` 直接 `Vec::with_capacity`（u64::MAX 时 capacity overflow
   panic）；`read_header` 的 `buf.reserve(n_contigs * 16)` 在 u32::MAX 时
   预留 ~64 GiB（OOM abort）。修复：头解析增加与构建器一致的上下限校验
   （n_contigs ≤ u16::MAX、n_records ≤ u32::MAX），预分配改
   `try_reserve_exact`；顺带修复 `to_hv::read_hv` 的 name/dim 未校验分配。
   回归测试 `crafted_record_count_rejected_not_panic`、
   `crafted_contig_count_rejected_not_panic`、
   `crafted_hv_header_rejected_not_panic`。
8. **e-align span 过滤两处溢出**：`(t_end - t_start) as usize` 对畸形记录
   （t_end < t_start）回绕成超大 span 通过长度过滤；改为 i64 运算
   `.max(0)` 后再转 usize（极端 i32 坐标下减法本身也会溢出，i64 修复）。
9. **lav `l` 行负跨度静默产生回绕垃圾**：`l 5 5 1 1 95` 这类
   t_end < t_start 的行通过 `parse_a` 后，`blocks_to_psl` 的
   `(t_end - t_start) as u32` 回绕成超大 block。修复：
   `t_end < t_start || q_end < q_start` 直接报 InvalidData。回归测试
   `negative_span_l_line_rejected`。
10. **pgi build u32 溢出防护**：`pos_start: positions.len() as u32` 在
    >42 亿 k-mer 记录时静默截断索引。加 `payloads.len() <= u32::MAX` 防御
    检查（无法用真实规模测试，代码审查覆盖）。
29. **非 UTF-8 临时目录路径 panic（记录项落实）**：`PipelineCtx::new`/
    `enter` 与 `sd search`/`sd cross` 对 `tempfile::TempDir` 路径做
    `to_str().unwrap()`，`$TMPDIR` 含非 UTF-8 字节时 panic（Zero Panic
    违反）。统一改用 `io::path_to_str` 返回友好错误。
35. **`align_banded_local` 序列长度悬殊时 DP 数组越界 panic**：
    `banded.rs` 的行循环把 j_lo 钳到 `[1, m]`（`min(m)`）而不是与对角带
    求交；当一条序列远长于另一条时（如 q=2000、t=100、band=4），行尾的
    `off = j - i - diag0 + band` 为负，`as usize` 回绕成超大下标，
    `mscore[c]` 直接越界 panic（debug 实测 "attempt to add with overflow"）。
    修复：j_lo/j_hi 改为 `max(1, i+diag0-band)` / `min(m, i+diag0+band)`，
    空交集直接跳过该行；off 恒落在 `[0, width)`。当前 align 管线的调用
    （等长窗口）不触发，但作为公开 lib API 必须零 panic。回归测试
    `unbalanced_lengths_do_not_panic`。
36. **lav `l` 行极值坐标减法溢出 panic**：`parse_a` 用 `parse::<i64>()? - 1`
    解析 1-based 起点（i64::MIN 时下溢），跨度一致性检查
    `(q_end - q_start) != (t_end - t_start)` 在 `q_end=i64::MAX`、
    `q_start=-1` 时也溢出；构造输入在 debug 构建直接 panic。修复：起点
    `checked_sub(1)`、跨度用 `checked_sub` 比较，越界值报 InvalidData。
    回归测试 `extreme_l_line_values_do_not_panic`。

### 功能正确性 / 算法（13 处；#1-3 为重大索引/链算法缺陷）

1. **（重大）pgi 索引 k-mer key 与位置错配（所有 pgi 工作流受影响）**：
   `collect_one_contig` 两处叠加——closed-syncmer 会把同一位置选中两次
   （窗口 last-min 与 first-min 重合），且 flush 循环用当前滚动 key 给所有
   待处理位置配 key，多位置同迭代弹出时早期位置的 key 错配成后者的 k-mer。
   实测 2 Mb 随机基因组索引多出 39% 错配/重复记录，self 对齐产出 101 条
   伪自比对块。修复：pending 入队去重（HashSet）；flush 时仅当
   `pos + k - 1 == i` 才用滚动 key，否则按位置重算（含 N 防护）；RC 记录
   改用 `rc_key(key, k)`。回归测试
   `index_records_match_sequence_positions`。
2. **（重大）tube 排序键 anti/bucket 溢出（>8 Mb 基因组 tube 失效）**：
   `chain_tubes` 的 u128 排序键中 anti 只分配 24 位，a_pos + b_pos ≥ 2^24
   时高位污染 bucket 字段，tube 链被碎片化/丢失（8.9 Mb 串联重复修复前
   tube 输出 0 条正确块）。修复：anti 扩到 32 位、bucket 扩到 32 位，
   radix 键位不变。回归测试
   `tube_sort_key_supports_large_anti_coordinates`。
3. **（重大）tube 排序键负对角线回绕（>64 Mb 间距的重复失效）**：
   `BUCK_OFF = 1,000,000` 只覆盖到 -64 Mb；更深负对角线 `as u64` 回绕后
   高 32 位 0xFFFF 污染 strand/contig 字段，跨 contig 交错排序、链被切断
   （多 contig 深负对角线场景实测 0 条 tube）。修复：`BUCK_OFF` 改为
   `1 << 26`。回归测试 `tube_sort_key_supports_deeply_negative_diagonals`
   与 `tube_sort_key_does_not_mix_contigs_at_negative_diagonals`。
12. **cluster.rs 重叠 union 漏连嵌套区间**：按 start 排序后仅用相邻对做
    union-find，A 包含 B、C 且 B/C 不相交时 C 被拆成独立簇。修复：扫描时
    跟踪已见区间的最大右端，与先前区间重叠者必与最大右端区间重叠。回归
    测试 `nested_overlapping_intervals_form_one_cluster`。
13. **sd cluster 区间去重键忽略链向与物种**：同坐标异链向的两个拷贝
    （回文型倒位重复）被折叠成单个区间、丢失一个拷贝。修复：键改为
    `(species, chrom, start, end, strand)`。回归测试
    `same_coordinates_on_opposite_strands_are_distinct_copies`。
14. **rept s-align 带点 contig 名被 spanr 截断**：除 s_align 外所有 spanr
    管线都做了 `c1..cN` 名字映射，s_align 漏做，输出 key 被截断
    （`NC_000913.1` → `1`），下游 `fa mask` 无法命中。修复：复用
    chrom.sizes 建映射，改名后再跑 spanr，最后恢复真实名字。回归测试
    `command_rept_s_align_dotted_name`。
15. **Profex `-z` 坐标转换 + e-kmer 染色体尾部丢失**：Profex 输出
    `start` 是 0-based k-mer 起点（转 1-based 需 +1），`end` 已是
    1-based inclusive 末位坐标（不应再 +1）——原代码两端都 +1，每个区间
    右端偏一位；且 Profex 从不打印每个 read 最后一个 run 的 end/depth，
    原正则整段丢弃，e-kmer 丢失染色体尾部区间。修复：end 不再 +1；
    `run_profex_per_chr` 增加染色体长度参数，无深度阈值（e-kmer）时用
    染色体长度闭合尾 run 并输出（`-p:repeat` 只在命中重复库时非零，无假
    阳性风险），有阈值（s-kmer）时保守丢弃（深度未知，与 anchr 参考管线
    一致）。回归测试 `command_rept_e_kmer_tandem_coordinates`（2×2000 bp
    串联输出 "1-4000"，修复前 "1-2001"）与 s-kmer 坐标增强
    （"1-2000"）。
16. **sd cluster / sd run 不支持普通 gzip**：`sd search`/`sd align` 接受
    普通 `.gz`，但 cluster 的 `.loc` 随机访问索引只支持纯文本与 BGZF——
    普通 gzip 被当纯文本读，生成垃圾 `.loc` 并在 cluster 步骤中途失败。
    修复：`cluster_paf` 对非 BGZF 的 `.gz` 先解压到临时文件再建索引
    （TempDir 随用随清，不再污染用户目录）。回归测试
    `command_sd_run_gzipped_genome`。
17. **align lastz 省略 query 未真正启用 self 模式**：省略 query 或 `--self`
    时 `is_self` 仍为 false，目录输入跑 n×n 全交叉、单文件不传 `--self`
    （产生序列 vs 自身的平凡比对），与文档"省略 query 即 --self"相悖。
    修复：传 `self_mode`。回归测试 `command_align_lastz_omitted_query_is_self`。
26. **`psl lift` 负链记录的外层坐标提升错误（UCSC 约定违反）**：
    `Psl::lift_query`/`lift_target` 对 `-` 链记录把 `qStart/qEnd`
    （`tStart/tEnd`）按 `real_size - end_0` 平移。但按 UCSC PSL 规范（以及
    pgr 自己的两个产出方 `pgi align` 与 `lav_to_psl`），`qStart/qEnd` 恒为
    正向坐标，只有块级 `qStarts` 才在反向互补坐标系中。正确做法：
    `qStart/qEnd += start_0`，`qStarts += (real_size - end_0)`。
    实测影响：s-align 中间产物 `lifted.psl` 中 `-` 链记录
    `chr:1201-1400` 窗口的 qStart/qEnd 输出 `[4200,4301)`，正确值应为
    `[1200,1301)`。s-align 最终输出不受影响（`psl to-range` 只消费
    `qStarts`，双重反转恰好抵消），但任何直接消费 lifted PSL 的下游会拿到
    错误坐标。原 `test_lift_basic`/`test_lift_target` 夹具本身不符合 UCSC
    约定（`-` 记录的 qStarts 与 qStart/qEnd 同框），把错误行为固化进了
    测试。修复 `lift_query`/`lift_target`、夹具改为符合约定
    （`qStarts = qSize - qEnd`），并新增回归测试
    `test_lift_minus_strand_forward_coordinates`（含 to-range 往返校验）。
27. **s-align / e-align 的 soft-mask 警告误报 N gap**：检测用
    `pgr fa masked`（同时报告 lowercase 与 N 区域），纯 N gap 基因组也会
    触发"lowercase soft-mask"警告并误导用户 `tr a-z A-Z`。改为
    `has_soft_mask` 直接扫描 FASTA 的 lowercase 碱基（内存按记录有界），
    新增单元回归测试 `soft_mask_detection_ignores_n_gaps`。
31. **greedy 链合并导致倒位 SD 漏检（上一轮"待决策"项，已定方案 1 落实）**：
    `merge_adjacent_chains` 把同一对角线、间隔在 `--merge-gap` 内的两条链
    合并。倒位重复的两个互反链共享同一对角线（实测 diag 均为 500），间隔
    1.8 kb 无种子，被合并成一条嵌合链；单窗口扩展（16 kb）把两端 100%
    拷贝与中间非同源区连成一个 78% identity 大块，`sd search` identity
    过滤整块丢弃。同一输入 lastz 与 tube 均能检出。修复：合并条件增加
    `|diagA − diagB| > 0`——真实 IS 插入会平移对角线（与 merge 的设计初衷
    "对角线平移断链"一致），同对角线纯间隔是两个独立同源块。同步更新
    `merge_adjacent_chains_stitches_syntenic_blocks`（原"1.2 kb 插入"场景
    对角线差为 0，与真实插入不符，改为 4 bp 平移）并新增
    `command_sd_search_pgi_inverted_repeat` 回归测试（构造
    rand+dup+rand+rc(dup)+rand，检出两条 100% 的 1200 bp 倒位块；
    修复前 0 条）。`sd run` 端到端现在输出 2 set × 2 CORE 行。
37. **pgi merge 频率过滤两侧边界不一致**：A 侧 `ea_freq > freq` 保留
    `== freq` 的 k-mer，B 侧扩展范围 `occ >= freq` 与 FastGA GIXmake
    （"only k-mers whose count is less than the adaptamer frequency cutoff"）
    都丢弃 `>= freq`。交叉比对时（A 侧恰好 freq 次、B 侧罕见）pgr 会多出
    FastGA 不存在的种子。修复：A 侧改为 `ea_freq >= freq`；B 侧把
    `== freq` 的条目整体视为"不在索引中"（FastGA GIX 在构建期就排除
    `>= FREQ` 的 k-mer）——最大共享前缀扫描、occ 累计、发射循环三处统一
    用 `>= freq` 跳过 / `< freq` 计入。原实现里 `== freq` 的 B 侧条目会
    抬高 m 并触发整段丢弃，与 FastGA 不符（FastGA 中该条目不存在，同范围
    的稀有条目仍可命中）。回归测试
    `freq_boundary_drops_exact_freq_on_reference_side` 与
    `exact_freq_query_entries_are_absent_not_range_killers`。
    （第四轮精修）复查时发现还有两处残留不一致：B 侧"最大共享前缀 m
    扫描"与"扩展范围 occ 累计"仍把 `== freq` 的条目视为存在（`> freq` /
    `<= freq`），而 FastGA GIXmake 在构建期就排除 `>= FREQ` 的 k-mer——
    这些条目应整体视为"不在索引中"。统一改为：m 扫描 `>= freq` 跳过、
    occ 累计 `< freq` 计入、发射循环 `< freq`。原
    `merge_filters_extended_range_by_occurrences` 测试断言的是旧行为
    （`== freq` 条目丢弃整段），按 FastGA 语义改写为"多个低于阈值的条目
    总和 >= freq 才丢弃整段"，并新增
    `exact_freq_query_entries_are_absent_not_range_killers`
    （`== freq` 条目旁稀有条目仍以 25 bp 前缀命中）。

### 输入校验 / 静默错误（3 处）

11. **repeat.rs 两处吞 IO 错误**：`map_while(Result::ok)` 在读取中途出错时
    静默截断（e_align 的 PSL 过滤、run_profex_per_chr 的 Profex 输出）。
    改为 `let line = line?;` 传播错误。
32. **e-align PSL 过滤静默跳过畸形行（记录项落实）**：`run_align_repeat_
    pipeline` 对解析失败的 PSL 行 `let Ok(psl) = ... else { continue }`
    静默丢弃，与其余命令的 `parse_or_warn`（warn + 跳过）不一致，畸形输入
    难以诊断。补 `log::warn!`。
34. **decompose 对无法解析的 FASTA 头静默丢弃（记录项落实）**：
    `decompose::parse_fasta` 对解析失败的记录头静默跳过，畸形 cluster
    FASTA 的输入被悄悄忽略。补 `log::warn!`（与 `parse_or_warn` 行为一致）。

### 外部工具与参数 / CLI / 文档（11 处）

18. **lastz 静默失败**：`run_lastz` 对失败只打日志、返回 Ok，所有 job 失败
    时调用方拿到空结果无提示。改为统计失败数并 bail（实测损坏输入报
    `lastz failed for 1 of 1 jobs`）。
19. **lastz 失败原因被吞**：`cmd.status()` 丢弃 stderr，无法诊断。改为
    `cmd.output()` 并记录首个失败的 stderr（多 contig 报 "contains more
    than one sequence"，gz 报 "bad fasta character"）。
20. **参数校验缺失/不一致**：e_align 与 sd run/search/cross 的
    `--min-identity` 无范围校验（>1 全拒、<0 全过）或放行 0.0 与文案
    `(0,1]` 不符；e_align 的 kmer/smer/window/parallel、trf 的 minscore
    无正值/有限性校验。统一为 `> 0.0 && <= 1.0`、`> 0`、非负有限整数，
    帮助文本同步为 "(0, 1]"。
21. **trf 特殊字符文件名**：`fa split name` 生成 `sanitize(name).fa`，trf
    却用原始名拼 `${chr}.fa`——染色体名含 `/\():` 或双下划线时找不到文件。
    改用 `sanitize_filename(chr)`。回归测试 `command_rept_trf_special_chars`。
22. **sd search/cross 的 `--preset` 默认值未注册**：help 声称默认 set01 但
    clap 无 default_value，省略时完全不应用预设参数。补
    `.default_value("set01")`。回归测试 `command_sd_search_lastz_default_preset`。
23. **sd run --engine lastz --preset 参数拼装错误**：`--preset set01` 被拼
    成单个字符串经 `${var}` 传入，内层 clap 收到一个 argv 报
    "unexpected argument"。改用 `Vec<String>` + `$[preset_args]` 列表展开。
    回归测试 `command_sd_run_lastz_preset_parses`。
24. **噪音与帮助文本**（多处小修）：lav `m { ... }` mask stanza 每次转换
    warn 一行（新增 `LavStanza::Mask` 静默跳过）；sd search lastz 过滤不
    跳过 `#` 元数据行导致逐行 warn（改为跳过）；align lastz 帮助虚假声明
    `[multiple]` 修饰符（删除并注明单序列输入约束）；align-lastz 文档的
    `-s` 短选项实际不存在（改 `--preset`）；align.md 的 lastz 示例
    `-o out.psl` 错误（输出为 LAV 目录）；pgi build 帮助残留 "(12,8)
    canonical" 与默认 syncmer 8/5 不符（改正）。
25. **文档一致性**：docs/rept.md 补 e-align 章节；s_align soft-mask 警告
    行为补进帮助与用户文档；align-pgi 注明 `.gz` 输入的索引命名为
    `<name-without-.gz>.pgi`；docs/sd.md 同步 lastz 单序列/纯文本约束与
    SD 搜索前勿用自比对的提醒；TnCentral 示例路径修正。
28. **帮助文本与文档不一致（lastz 单序列约束）**：docs/sd.md 已注明 lastz
    引擎要求单序列 FASTA，但 `sd search`/`sd cross` 的命令内帮助漏掉该
    说明；`sd cross` 的文档段落也未注明。已同步四处。
30. **e-align identity 定义未说明（记录项落实）**：e-align 的
    `--min-identity` 用 gap-compressed identity（不含 insert 碱基），与 sd
    的 `(matches+repeats)/block_len` 不同，rept.md 与命令帮助均未说明。
    已补文档说明。
33. **主帮助文本 rept 子命令列表遗漏**：`src/pgr.rs` 的 about 帮助里
    `rept - Repeat detection: e-kmer, s-kmer, trf` 少了 `e-align` 与
    `s-align`，与实际注册的 5 个子命令不符。已补齐。

## 与外部参考实现的语义一致性核对

### #26 的 UCSC kent 源码复核

用 UCSC 官方 kent 仓库源码
`src/hg/utils/pslLiftSubrangeBlat/pslLiftSubrangeBlat.c` 的 `liftSide`
逐字段复核 #26 的修复，结论一致：

* kent 行为：`qStart/qEnd += subrange_start`；`-` 链时仅块级
  `qStarts += (seqSize - subrange_end)`（`reverseIntRange` 的结果）。
* kent 官方测试 `tests/input/qSubrange.psl` → `expected/qSubRangeTest.psl`
  的 `-` 链记录（chr21:33043300-33104451，qSize 61151，qStart 12，
  qEnd 61131，qStarts 首块 20）lift 后 qStart 33043312、qEnd 33104431、
  qStarts 首块 15025464，满足 `qStarts[0] = seqSize - qEnd =
  48129895 - 33104431 = 15025464`，且与同文件未子范围化的 chr21
  记录逐字段一致。原夹具/实现的"负链 qStarts 与 qStart/qEnd 同帧"
  在 kent 生态中不存在。

**子范围命名约定：pgr 1-based vs kent 0-based（记录在案，不改）**：kent
的 `liftSide` 把 `chr21:33043300-33104451` 的数字直接当 0-based 偏移用
（`+regStart` 无 `-1`，测试输出证实；也因此才能与未子范围的记录完全吻合）。
pgr 的 `parse_subrange`（intspan）按 1-based inclusive 解释
（`start - 1` 为偏移）。两者正链 qStart/qEnd 差 1 bp；负链块级 `qStarts`
偏移同为 `size - end`，不受影响。**决策：保持 pgr 1-based 约定**——生成端
`fa window`（`name:start+1-end`）与消费端 `psl lift` 自洽，`rept s-align`
端到端无 off-by-one，`psl lift` 帮助与 `docs/psl.md` 均已写明 1-based。
若未来直接消费 UCSC/blat 生态产出的子范围名，需先确认其语义（blat 文档
为 1-based，但 kent 工具按 0-based 处理）。

### #31 的 FastGA 源码复核

对照 `FASTGA-main/FastGA.c`（Gene Myers，V1.5）与 `ALNchain.c`，修复方向
与 FastGA 语义一致：

* `align_contigs`（FastGA.c:2973）链化是**单次 anti 扫描**：种子按
  anti 归并，`anti < ahgh + CHAIN_BREAK`（2000 anti，=2×-s 1000）
  并入当前链并累积 cov，否则**立即断链**；`cov >= CHAIN_MIN` 的链
  各自触发 tube wave 对齐。源码中**没有任何事后合并/缝合相邻链的
  步骤**——同对角线、间隔超阈值的链就是两条独立链。
* 链只跨两个相邻对角线桶（cdiag/cdiag+1，桶宽 `BUCK_WIDTH=64` →
  128 宽带）；对角线平移超带宽的种子落不同桶，不可能同链。FastGA
  同样不缝合（IS 断开的两个块由下游 chainnet 连接）。
* `ALNchain.c` 的 `localChain`/`KDRangeChain` 是比对块层面的 KD 树
  范围 DP 链化，连接条件是 X/Y 二维 gap ≤ maxGap 且得分增益为正、
  `backtrackLocal` 用 maxDrop 断链——也不缝合间隔超阈值的独立块。
* 数值换算：同对角线时 anti 间距 = 2×pos 间距，FastGA
  `CHAIN_BREAK=2000` anti 等价于 pgr greedy 的 `max_gap=1000`
  （pos）；`merge_gap=5000` 是 pgr 独有的一层缝合阈值，FastGA 无
  对应物。修复前 1.8 kb 间隔（>1000 bp，FastGA 会断链）被 pgr
  merge 合并；修复后同对角线不合并，与 FastGA 行为一致。
* 结论：#31 的"同对角线纯间隔 = 两个独立同源块、仅对角线平移才
  缝合"与 FastGA 链化语义一致；pgr 的对角线平移缝合是自身扩展
  （为 IS 场景），与 FastGA 无冲突（平移链在 FastGA 里也是独立
  tube，由下游 chainnet 连接）。

## 排除的疑点（经核验无需修复）

* **`sd run` 的 cluster set_id 合并**：`sid + set_offset`、`set_offset +=
  cluster_max` 的重编号值域 `[offset+1, offset+max]` 各 cluster 两两不相交，
  即使 set_id 不连续（实际 `sd decompose` 从 1 连续编号）也不可能碰撞；
  空 cluster / 短行跳过与 `sd cover` 的 `read_elems` 一致。无需修复。
* **60,423 → 75,413 数据勘误**：e_align 对 MG1655+tncentral 的历史记录
  差异来自 tncentral 库更新（6073 → 6093 条）与编译产物时序，非代码 bug，
  笔记 repeat-masking.md §2.3.5 已加勘误。
* **sd cluster minus 链序列提取误报**：一度怀疑提取序列错误，实为构造的
  测试 PAF 坐标语义问题；按 pgr 内部 PAF 约定（qstart/qend 恒为正向坐标）
  提取的序列与区间反向互补逐碱基一致。
* **wave `forward_wave_mid` 初始 trim 越界**：无 PATH_AVE 质量匹配时返回
  初始 trim 点，理论上一端可越界 1 bp 或超出序列；经几何推演（越界点
  必然镜像到空盒 → 反向 wave 返回 None，或 `at <= ab || bt <= bb` 守卫
  拦截）与 8×8 长度 × 全 band/amid 参数空间的 fuzz（约 20 万次调用）均
  未触发 panic，判定不可达，不加防御代码。
* **`spanr fill -n 0`（e-align 的 fk=0）**：实测 spanr 0.8.6 对 `-n 0`
  fill 为 no-op，行为与设计文档（e-align 不做 k-mer 填补）一致，仅多一次
  冗余进程，不改。
* **LAV `s`/`h` stanza 含空格文件名**（`[nameparse=darkspace]` 场景）：
  `split_whitespace` 解析会把后续字段错位；与 UCSC lavToPsl 的简单解析
  一致，pgr 自身产线文件名不含空格，记录不修。

## 记录项（未改，低风险 / 待决策）

* **tube 工作流对"库 vs 基因组"的结构性失效**：对照实验（酵母 + repbase）
  greedy 出 2220 个 PSL 块、tube 只有 4 个，根因是 tube 按相邻对角桶对
  独立 merge、库比对种子稀疏导致跨桶链被切断。该结论基于修复前代码，
  syncmer/排序键修复后 tube 行为可能显著改善，待真实数据重测。
* `decompose.rs` 负链投影依赖 header 与序列长度一致（cluster 内部保证）；
* cluster/cover 的 u32→i32 坐标转换（仅 >2.1 Gb 染色体才溢出）；
* `run_lastz` self 模式仍构建 n×n job 列表，运行时才按 basename 跳过
  （大目录 self 对齐可提前过滤）；
* `syncmer.rs` 参考实现与 `collect_one_contig` 是两套独立实现，前者也会
  重复发射同一位置，消费方用 HashSet 去重无实际影响，可后续合并；
* wave.rs 的 `unreachable!`/`panic!` 均为 Myers 算法不变量，有穷举与随机
  测试兜底，未发现用户输入可触发路径；
* s_align / sd search --engine pgi 传不支持的类型（2bit / 目录）时，内层
  命令错误被 cmd_lib 包装成不透明的 "Running [...] exited with error"——
  不 panic、有报错，仅可读性差；
* `fa split name` 名称碰撞（`chr(1)` 与 `chr_1` 同 sanitize 结果）：trf
  依赖 split 产物，但碰撞根因在 `fa split name`（共享命令），概率极低，
  记录不修。

## 已知限制（有意保留）

* **子范围命名约定：pgr 1-based vs kent 0-based**（#26b 决策）：pgr 的
  `psl lift` 按 1-based inclusive 解释子范围名，与生成端 `fa window` 和
  消费端 `psl to-range` 自洽；直接消费 UCSC/blat 生态的子范围名时需先
  确认语义。
* **s-kmer 对染色体尾部的重复保守丢弃**：Profex `-z` 不输出末 run 深度，
  无法区分唯一尾与重复尾，有阈值时保守丢弃；与 anchr 参考管线一致。
* **单 contig > 4.3 Gb 的 pgi 索引**：`pack_position` 的 pos 为 u32，
  `write` 会为超长 contig 计算 pos_bytes=5 而 reader 拒绝 `> 4`，属自相
  一致的上限（超长单 contig 不被支持）。

## 验证

### 引擎交叉验证

同一 tandem 基因组上 pgi 与 lastz 引擎检出同一对 1200 bp 拷贝（t
0-1200 ↔ q 1200-2400 及反向），pgi 边界修剪 2 bp（shared-k-mer 种子边界），
坐标一致；倒位重复（fwd + RC 拷贝）经 chainnet（无 --syn）后两条 `-` 链
PAF 保留，`sd run` 端到端输出两个坐标正确的 elementary SD（同 set_id、
CORE）。

### 端到端坐标正确性

* `rept trf` 对 100N + 1000 bp 串联输出精确 "101-1100"；
* `rept s-align` 对 500 + 2×1200 + 400 构造输出精确 "501-2900"
  （1-based inclusive 全覆盖）；
* 倒位重复 X + rc(X) 整段 [500,2900) 为回文，[513,2887) 100% 负链自比对
  成立，双链 CORE 行坐标正确；
* 多 contig s-kmer 的 read 编号与 chr.sizes 顺序一致；
* e-kmer `--keep-index` 二次运行复用缓存、输出逐字节一致，touch 库文件后
  失效重建；
* `sd run` 对普通 gz 与 BGZF fixture（mg1655.fa.gz）全流程成功；
* `sd cross` 双引擎对普通 gz 输入正常。

### 鲁棒性（Zero Panic）

畸形输入扫测——截断 lav / 负跨度 lav / 越界 PAF / 垃圾 BED / 垃圾 .pgi /
构造头 / 空 FASTA / 空 PAF / 空库 / 全 N / 极小基因组 / N 密集 / poly-A /
GC-rich / `/dev/urandom` / `-k 0` / `--step 0` / `--window 0` /
`--chunk-records 0` / `--min-shared 0` / 短 PAF 行——全部友好报错或空输出，
无 panic。索引/序列校验（contig 名与长度不匹配）友好报错；s-align /
s-kmer 无重复基因组输出 `"-"` 占位，下游 `fa mask` 正确消费。

### 最终状态

* 测试数演进：956 → 958 → 960 → 962 → 968 → 970 → 972 → 973 → 974 →
  975 → 976 → 977 → 995（每处修复均带回归测试或端到端验证；第四轮 #37
  精修后 995）。
* `cargo fmt --check` 与 `cargo clippy --all-targets -- -D warnings` 干净
  （plot/dot.rs 3 个与 pgi/align.rs 1 个测试目标既有告警，收尾轮已一并
  修掉）。
* 第一轮连续两轮无新问题后追加第二轮；两轮修复均经外部源码复核（UCSC
  kent `pslLiftSubrangeBlat.c`、FastGA `align_contigs`），无新增问题。
* 第五轮：重读本轮全部改动（banded 边界、lav checked 算术、pgi 频率过滤
  三处调用点）与剩余高风险路径（mmap 越界守卫、tube 合并循环、wave 端点
  几何推演），无新问题。
* 第六轮：`cargo test` 全量 995 通过、`cargo fmt --check` 与
  `cargo clippy --all-targets -- -D warnings` 干净；sd/rept/align 四组 CLI
  端到端（38 测试）全绿。
* 连续两轮无新问题，符合"修复后继续复核直到稳定"的审计约定。

## 提交状态

* 本文件所述修复已随 sd/rept/align 系列提交逐步落地（`53a8a78`、
  `b911869`、`bcc1a9e`、`801e21d`、`f752a3a` 等）。
* 当前工作区仅含本文档结构调整（与 `audit-runlist-rg.md` 统一骨架），
  无代码改动。
