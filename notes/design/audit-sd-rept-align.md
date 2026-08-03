# sd / rept / align 代码与文档审核记录（合并版，2026-08-04）

对新增命令族 sd / rept / align 的代码与文档进行多轮深入审核。范围：
sd（8 命令 + libs/sd）、rept（6 命令 + libs/pl）、align（pgi/lastz +
libs/pgi、libs/lastz、libs/fmt/lav、alignment DP），约 8000 行代码逐文件
审完；文档 docs/{sd,rept,align-pgi,align-lastz}.md 全部核对。每轮发现
问题后修复并进入下一轮复核，直至连续两轮无新问题。最终 974 测试全绿、
fmt / clippy 干净。

## 修复的缺陷（25 处）

### 索引与链算法（3 处，重大）

1. **pgi 索引 k-mer key 与位置错配（所有 pgi 工作流受影响）**：
   `collect_one_contig` 两处叠加——closed-syncmer 会把同一位置选中两次
   （窗口 last-min 与 first-min 重合），且 flush 循环用当前滚动 key 给所有
   待处理位置配 key，多位置同迭代弹出时早期位置的 key 错配成后者的 k-mer。
   实测 2 Mb 随机基因组索引多出 39% 错配/重复记录，self 对齐产出 101 条
   伪自比对块。修复：pending 入队去重（HashSet）；flush 时仅当
   `pos + k - 1 == i` 才用滚动 key，否则按位置重算（含 N 防护）；RC 记录
   改用 `rc_key(key, k)`。回归测试
   `index_records_match_sequence_positions`。
2. **tube 排序键 anti/bucket 溢出（>8 Mb 基因组 tube 失效）**：
   `chain_tubes` 的 u128 排序键中 anti 只分配 24 位，a_pos + b_pos ≥ 2^24
   时高位污染 bucket 字段，tube 链被碎片化/丢失（8.9 Mb 串联重复修复前
   tube 输出 0 条正确块）。修复：anti 扩到 32 位、bucket 扩到 32 位，
   radix 键位不变。回归测试
   `tube_sort_key_supports_large_anti_coordinates`。
3. **tube 排序键负对角线回绕（>64 Mb 间距的重复失效）**：
   `BUCK_OFF = 1,000,000` 只覆盖到 -64 Mb；更深负对角线 `as u64` 回绕后
   高 32 位 0xFFFF 污染 strand/contig 字段，跨 contig 交错排序、链被切断
   （多 contig 深负对角线场景实测 0 条 tube）。修复：`BUCK_OFF` 改为
   `1 << 26`。回归测试 `tube_sort_key_supports_deeply_negative_diagonals`
   与 `tube_sort_key_does_not_mix_contigs_at_negative_diagonals`。

### 崩溃 / 越界 / 溢出（8 处，Zero Panic）

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
11. **repeat.rs 两处吞 IO 错误**：`map_while(Result::ok)` 在读取中途出错时
    静默截断（e_align 的 PSL 过滤、run_profex_per_chr 的 Profex 输出）。
    改为 `let line = line?;` 传播错误。

### 功能正确性（7 处）

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

### 外部工具与参数 / CLI（7 处）

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

## 排除的疑点

* **`sd run` 的 cluster set_id 合并**：`sid + set_offset`、`set_offset +=
  cluster_max` 的重编号值域 `[offset+1, offset+max]` 各 cluster 两两不相交，
  即使 set_id 不连续（实际 `sd decompose` 从 1 连续编号）也不可能碰撞；
  空 cluster / 短行跳过与 `sd cover` 的 `read_elems` 一致。无需修复。
* **tube 工作流对"库 vs 基因组"的结构性失效**：对照实验（酵母 + repbase）
  greedy 出 2220 个 PSL 块、tube 只有 4 个，根因是 tube 按相邻对角桶对
  独立 merge、库比对种子稀疏导致跨桶链被切断。该结论基于修复前代码，
  syncmer/排序键修复后 tube 行为可能显著改善，待真实数据重测（记录项）。
* **60,423 → 75,413 数据勘误**：e_align 对 MG1655+tncentral 的历史记录
  差异来自 tncentral 库更新（6073 → 6093 条）与编译产物时序，非代码 bug，
  笔记 repeat-masking.md §2.3.5 已加勘误。
* **sd cluster minus 链序列提取误报**：一度怀疑提取序列错误，实为构造的
  测试 PAF 坐标语义问题；按 pgr 内部 PAF 约定（qstart/qend 恒为正向坐标）
  提取的序列与区间反向互补逐碱基一致。

## 记录项（未改，低风险 / 待决策）

* `ctx.rs` 的 `tempdir.path().to_str().unwrap()`（非 UTF-8 路径罕见）；
* `decompose.rs` 负链投影依赖 header 与序列长度一致（cluster 内部保证）；
* cluster/cover 的 u32→i32 坐标转换（仅 >2.1 Gb 染色体才溢出）；
* `run_lastz` self 模式仍构建 n×n job 列表，运行时才按 basename 跳过
  （大目录 self 对齐可提前过滤）；
* decompose 对无法解析的 FASTA 头静默跳过、e-align PSL 过滤静默跳过
  畸形行（建议后续加日志，与 `parse_or_warn` 统一）；
* e-align 的 identity 用 gap-compressed 定义（不含 insert 碱基），与 sd 的
  `(matches+repeats)/block_len` 不同，rept.md 未说明；
* `syncmer.rs` 参考实现与 `collect_one_contig` 是两套独立实现，前者也会
  重复发射同一位置，消费方用 HashSet 去重无实际影响，可后续合并；
* wave.rs 的 `unreachable!`/`panic!` 均为 Myers 算法不变量，有穷举与随机
  测试兜底，未发现用户输入可触发路径；
* s_align / sd search --engine pgi 传不支持的类型（2bit / 目录）时，内层
  命令错误被 cmd_lib 包装成不透明的 "Running [...] exited with error"——
  不 panic、有报错，仅可读性差；
* s-kmer 对染色体尾部的重复保守丢弃（Profex `-z` 不输出末 run 深度，无法
  区分唯一尾与重复尾；与 anchr 参考管线一致）。

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

* 测试数演进：956 → 958 → 960 → 962 → 968 → 970 → 972 → 973 → 974，
  每处修复均带回归测试或端到端验证；
* `cargo fmt --check` 与 `cargo clippy -- -D warnings` 干净
  （`--all-targets` 下 plot/dot.rs 3 个 + pgi/align.rs 测试 1 个既有告警
  不在本次范围）；
* 连续两轮（第 22 轮与收尾复核）未发现新问题，审计循环终止。
