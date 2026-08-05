# align 命令族代码审核记录（2026-08-04 / 2026-08-05）

对 `pgr align` 命令族（pgi/lastz + `libs/pgi`、`libs/lastz`、`libs/fmt/lav`、
`libs/fmt/psl`、`alignment` DP）约 3000 行代码及全部文档
（`docs/{align-pgi,align-lastz,lav,psl}.md`）进行审核。缺陷按类别分组记录；
关键修复均附回归测试，验证概况见文末"验证"一节。

> 注：`pgr sd` / `pgr rept` 命令族的审核记录见 `notes/audit/audit-sd-rept.md`。
> `libs/pgi` 索引的**构建**与读取缺陷（k-mer key、构造头、记录越界、sibling
> 索引、`--parallel` 等）记录在本文件；sd 对 pgi 的**消费**缺陷（`sd search`
> 传 `.pgi` 拒绝、greedy/tube 链逻辑等）记录在 audit-sd-rept.md。

审核范围：
- **align pgi**：`align pgi`（+ `libs/pgi` 索引构建/读取、`psl`/`sort` 输出）
- **align lastz**：`align lastz`（+ `libs/lastz` 封装、`libs/fmt/lav` 解析、
  `align_banded_local` 等 alignment DP）

审核重点：Zero Panic（畸形输入/构造索引不 panic）、数据安全（`-o` 不得覆盖
输入、陈旧/损坏索引不得静默复用）、确定性（`--parallel` 下输出可复现）、与
外部参考实现（FastGA GIX / UCSC kent pslLiftSubrangeBlat / lavToPsl / lastz）
的语义一致。

## 与外部参考实现的语义一致性核对

关键修复均对照官方源码复核，方向一致：

* `psl lift` 负链坐标提升：与 UCSC kent `pslLiftSubrangeBlat.c` 的
  `liftSide` 行为一致（子范围命名约定 pgr 1-based vs kent 0-based 为
  记录在案的有意差异）。
* pgi 索引 k-mer 频率过滤：`emit_entry_hits` 的 canonical key 过滤、
  `freq >= cutoff`（FastGA 语义，非 `>`）、前缀窗口 / 最大共享前缀 /
  扩展范围过滤，与 FastGA GIX 语义一致。
* 软掩码语义：`build_from_seqs` 新增 `mask`（小写→N，FastGA `-M` 语义），
  与 `build_from_path` 一致；pgi 与 lastz 对小写拷贝行为完全一致。
* LAV→PSL 转换：`blocks_to_psl` 的 q/t 正 gap 计 insert、负 gap clamp 忽略、
  1-based→0-based `checked_sub` 防溢出、`-` 链坐标翻转，与 UCSC
  `lavToPsl` 一致。
* 有意差异（已记录）：子范围命名 pgr 1-based vs kent 0-based；`stat`/
  `statop` 类输出格式差异不涉及本族。

## 排除的疑点（经核验无需修复）

* LAV `s`/`h` stanza 含空格文件名解析错位，与 UCSC lavToPsl 一致，记录不修。
* `lav to-psl` 对畸形 LAV 静默输出空：LavReader 跳过空行/注释、未知 stanza 有
  warn，lastz 输出不会畸形，属容错设计，记录不修。
* `pgi build` 的 `unwrap()`/`expect` 逐一核对：均在测试代码或 clap required/
  value_parser 约束下，运行期不可达，无潜在 panic。
* 全量扫描家族生产代码 `unwrap()`：`blocks_to_psl` 的 `last_mut().expect` 有
  非空前置保证；PgiReader 各字段解析均有长度/类型检查。无生产 panic 风险。
* `libs/pgi` 非测试代码 `unwrap`/`expect`/`unreachable!` 全量核对：均在定长
  切片取值、`Some` 前置守卫或 while 循环条件下不可达（如 `peak_rss_mb` 的
  `unwrap_or(0)`、`tubes_for_group` 的 `(None, None)` 由 `bi < len || mi < len`
  保证不可达），无生产 panic 风险。
* `chain_tubes` 的 u128 排序键字段布局：`a_contig`(16b)@89、`b_contig`(16b)@73、
  `strand`@72、`bucket`@40、`anti`@0。`anti`(≤2^33) 与 `bucket`(diag/64+2^26，
  ≤27b) 各占不足 40 位，字段间无重叠（bits 40–66 vs 72+），`a_contig` 上限 16
  位占 bits 89–104 < 128，无溢出/碰撞；`a_contig`/`b_contig` 的 `u16` 截断因
  `n_contigs <= u16::MAX` 约束（`build_from_seqs` 保证 + 头解析校验）而安全。
* reference 侧记录校验：`chain_to_psl` 直接以 `chain.a_contig` 索引
  `a.contigs[...]`；`validate_record` 只覆盖 query 侧，但 reference 侧逐路径
  安全——`align pgi` 命令只走 streaming 路径（`align_to_psl_streaming` /
  `align_to_psl_ext_streaming`），reference 经 `PgiStream::next_batch` 逐记录
  `validate_record(cid, pos, ...)`，`ac < contigs.len()`；非 streaming 的
  `merge_seed_hits` 其 `a` 来自 `PgiIndex::read`（逐记录校验）或
  `build_from_seqs`（构造即合法）。无 crafted reference 越界 panic 风险。
* `emit_entry_hits` 的 `a` 侧记录不重复校验但逐路径安全：resident 路径来自
  `build_from_seqs`（构造即合法）或 `PgiIndex::read`（逐记录 `validate_record`）；
  streaming 路径经 `PgiStream::next_batch` 逐记录校验；`b` 侧经 `emit_entry_hits`
  内 `validate_record` 校验。无越界 panic 风险。
* 负链坐标：`chain_to_psl` 对 minus 链用 `reverse_range_pair(b_start, b_end,
  b_len) = (b_len - b_end, b_len - b_start)` 与 `extend_window` 的 `reverse_range`
  一致；`b_end <= b_len` 由 `validate_record`（`pos + k <= len`）保证，无下溢；
  `q_starts` 内部块按 UCSC minus 约定置于 RC frame（`b_len - q_end`），回归
  `psl_block_coordinates` / `extend_chain_rc_query` 覆盖。
* `emit_entry_hits` 频率过滤两侧对称：`a` 侧 `ea_freq >= freq` 丢弃、`b` 侧
  最大前缀/扩展范围均按 `>= freq` 处理（FastGA GIX 语义），无误用 `>` 的残余。
* build `kmer_key_at` 切片越界排除：pending 位置均经 `start + k <= n` /
  `j + k <= n` 守卫后才入队（`build.rs`），`seq[pos..pos+k]` 恒在界内，无 panic。

## 记录项（未改，低风险 / 待决策）

* `align lastz --lastz-args` 的值以 `-` 开头时需用 `--lastz-args=<val>` 形式
  （clap 对空格形式的值为标准行为）；帮助文本未提示该写法。
* `align lastz -o dir` 重复使用旧 LAV 残留：影响链短（`sd run`/`s-align`/`sd
  search lastz` 均用临时 workdir 免疫），且 LAV 是通用扩展名清理易误伤，
  记录不修。
* `.pgi` 显式输入 + 冲突 `-k`/`--smer`/`--window` 被静默忽略（exit 0）：
  docs/align-pgi.md 明确说明 "apply only to genome-sequence inputs; .pgi
  inputs carry their parameters in the index header"——文档化预期行为，
  记录不修（sibling 索引路径的冲突报错是额外保护）。
* `--self` 只校验索引输入 `ref == query`，不校验 `--ref-seq` 与 `--query-seq`
  是否一致。`align pgi ref.pgi --self --ref-seq a.fa --query-seq b.fa`（两文件
  均匹配索引 contig 表但内容不同）会在 self 模式比对 a.fa vs b.fa，仅丢弃"精确
  自同一性"命中，对结果影响极小；属用户自相矛盾的请求，且 contig 校验在文件名/
  长度不符时已能拦截。按"简洁优先"记录不修。
* crafted 索引可携带超大 contig `len`（u64 无上限），`oriented = b_len - k -
  bpos` 超出 `u32` 时 `as u32` 截断、`chain_to_psl` 的 `*len as u32` 同样截
  断——仅产生错误坐标，不 panic。真实索引受"单 contig ≤ 4.3 Gb"的已知限制
  约束，故仅对畸形输入成立，按"简洁优先"记录不修。

## 已知限制（有意保留）

* 子范围命名 pgr 1-based vs kent 0-based（记录在案）：pgr 生成端/消费端
  自洽，直接消费 UCSC/blat 生态子范围名时需先确认语义。
* 单 contig > 4.3 Gb 的 pgi 索引：pos 为 u32，超长单 contig 不被支持。
* `ref.2bit` 与 `ref.fa` 同 stem 时共享 `ref.pgi` 兄弟索引（有意保留）：
  `.2bit` 是 `.fa` 的压缩变换，作为 drop-in 替换共享索引符合设计（docs
  明确 "2bit inputs are preferred"）。与 `.fa`/`.fa.gz` 的分离不同（两者
  可能是内容无关的独立文件），同目录下 `ref.fa` 与 `ref.2bit` 内容分歧属
  罕见用户错误，且 mtime 新鲜度检查在 `ref.2bit` 更新时自动重建，已部分
  缓解。若改为 `ref.2bit.pgi` 分离会破坏 drop-in 语义，故不修。

## 修复的缺陷（共 41 处：26 处代码/行为 + 15 处 CLI/帮助/文档）

### 崩溃 / 越界 / 溢出（Zero Panic，10 处）

**lav d stanza 边界差一越界**。修复：守卫改 `+ 6`。回归
   `truncated_d_stanza_errors_not_panics`。
**构造 .pgi/.hv 头容量溢出 panic/OOM**（未校验 n_records/n_contigs）。
   修复：头解析校验 + `try_reserve_exact`。回归 3 个 crafted 测试。
**lav `l` 行负跨度回绕成超大 block**。修复：t_end < t_start 等报
   InvalidData。回归 `negative_span_l_line_rejected`。
**pgi build `positions.len() as u32` 静默截断**（>42 亿记录）。修复：
    `payloads.len() <= u32::MAX` 防御检查。
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
**临时索引目录创建失败 panic**：`resolve_side` 的
    `TempDir::new().expect("creating temporary index directory")` 在系统临时
    目录不可写/磁盘满时 panic。修复：改为 `?` 传播友好错误（内部
    `tmp.as_ref().expect` 有刚赋值的前置保证，不可达）。
**`read_header` 按不受信任的 contig 名长度扩容**：`PgiIndex::read` /
    `PgiStream::open`（`align pgi` 的 reference 流式读取与 `read_index_params`
    都走它）逐 contig 读取 `nb` 后直接 `buf.resize(start + nb, 0)`；`nb` 来自
    头部、无上限（可到 `u32::MAX`），构造索引可迫使数 GB 分配（OOM abort）；
    而 mmap 路径经 `parse_header_bytes` 用实际缓冲区边界 `take_bytes` 校验
    `nb`，两路径不一致。此前只校验 `n_contigs`/`n_records`，未校验逐名 `nb`。
    修复：`read_header` 在扩容前校验 `nb <= MAX_CONTIG_NAME`（1 MiB），超限报
    友好错误。回归 `crafted_huge_contig_name_rejected_not_panic`。
**`emit_entry_hits` 前缀窗口在 k=64 时 `hi` 溢出**：`window(len)` 计算
    `hi = lo + r`；k=64 时 `k_bits=128`、`mask=u128::MAX`，若某 k-mer 前 `len`
    个碱基全为 T（位 104..128 置位），`lo = 2^128 - r`、`hi = lo + r = 2^128`
    溢出 `u128`（debug panic / release 静默回绕，前缀范围错乱）。k=64 经
    `-k 64` 可达（CLI 无上限、`build_from_seqs` 允许 k<=64）。mmap 路径对
    `2^(2k)` 哨兵有 `k<64` 的特判，但 `2^128` 无法表示。修复：`hi` 改
    `lo.saturating_add(r)`——饱和后仅排除全 T 的 `u128::MAX` 单键，非真实 seed，
    可接受。回归 `merge_k64_high_key_no_prefix_overflow`（构造 index 直接触发
    `window(12)` 的 `lo + r = 2^128` 路径）。

### 功能正确性 / 算法（10 处，含 1 处重大索引缺陷）

**（重大）pgi 索引 k-mer key 与位置错配**（2 Mb 随机基因组 39% 错配、
   self 比对 101 条伪块）。修复：pending 去重、flush 按位置重算 key、RC
   用 `rc_key`。回归 `index_records_match_sequence_positions`。
**align lastz 省略 query 未启用 self 模式**。修复：传 `self_mode`。
    回归 `command_align_lastz_omitted_query_is_self`。
**`psl lift` 负链外层坐标提升错误（违反 UCSC 约定）**。修复：
    `qStart/qEnd += start_0`、`qStarts += (size - end_0)`，夹具修正。
    回归 `test_lift_minus_strand_forward_coordinates`。
**`psl lift` 的 `parse_subrange` 误切含 `.`/`:` 的 contig 名**：窗口名
    `{contig}:{start}-{end}` 经共享 `Range` 解析器时，`NC_000913.1:1-200`
    被读成 name="NC_000913" + chr="1"、`chr1:alt:1-200` 被读成 chr="alt"，
    `lift_query` 在 sizes 表里查错键、静默跳过提升（仅 warn）。修复：
    `parse_subrange` 改为取最后一个 `:`+数字后缀切分，前缀整体作为 contig
    名；回归 `parse_subrange_keeps_dotted_and_colon_contigs`。
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
**FASTA 原地修改后兄弟索引被静默复用**：`resolve_side` 复用同名兄弟 `.pgi`
    时只校验 contig 名/长度；同名单长但序列不同的 FASTA 会静默复用旧索引
    （k-mer 来自旧序列），对齐结果错误。修复：新增 mtime 校验（输入比索引
    新则重建，与 e-kmer 缓存同一约定）。回归
    `command_align_pgi_stale_sibling_index_rebuilt`。
**`.pgi` 单输入自比对 + 仅 `--ref-seq` 报错**：`align pgi ref.pgi --ref-seq
    ref.fa` 报 "extension sequences are needed for both sides"——self 模式下
    query 侧复用 `.pgi` 输入的 `seqs=None`，两侧空/非空不一致触发 bail。
    修复：`resolve_seqs` 后 self 模式下任一侧扩展序列为空时复用另一侧（两
    方向对称）。验证：仅 `--ref-seq` / 仅 `--query-seq` / 双侧 / FASTA 直接
    输入四者输出逐字节一致。回归
    `command_align_pgi_single_ref_seq_on_self_pgi`。
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

### 数据安全（`-o` 同输入保护 / 陈旧索引 / 静默数据丢失，2 处）

**`align pgi` 的 `-o` 指向输入时静默覆盖输入**：`align pgi g.fa -o g.fa` 等把
    输入 FASTA/.pgi 覆盖为 PSL（exit 0、无提示）；`--ref-seq`/`--query-seq`
    同样可能被覆盖。修复：`align pgi` 对 `-o` 及 `--ref-seq`/`--query-seq`
    均加 `ensure_outfile_distinct`。实测覆盖输入报 "also an input file" 且
    输入完好。
**`align pgi` 的 `-o` 可静默覆盖兄弟索引**：`ensure_outfile_distinct` 只保护
    `[ref, query, ref_seq, query_seq]`，未包含基因组输入映射的兄弟索引路径
    （`ref.fa` → `ref.pgi`、`ref.fa.gz` → `ref.fa.pgi`）。实测 `align pgi
    ref.fa -o ref.pgi` 把 PSL 输出写到 `ref.pgi`，覆盖/破坏该兄弟索引，下一次
    运行时 `resolve_side` 把 PSL 当 pgi 读，报 "reading header / failed to fill
    whole buffer"。修复：在 `execute` 中把每个基因组输入的 `sibling_pgi_path`
    一并加入 `ensure_outfile_distinct`（跳过 `stdin`）。实测 `-o ref.pgi` /
    `-o ref.fa.pgi` / `--keep-index -o ref.pgi` 均报 "output file ... is also an
    input file" 且索引完好；正常 `-o out.psl` 不受影响。回归
    `command_align_pgi_output_not_overwrite_sibling_index`。

### 性能（1 处）

**`align pgi --parallel` 未约束自动索引构建的 rayon 并行度**：`resolve_side`
    （内部 `build_from_seqs` → `radix_sort_u128_par`）在自定义线程池创建前
    执行，索引构建走全局 rayon 池，`--parallel N` 只约束 merge/扩展阶段。
    文档承诺 "--parallel: rayon thread count"，行为不一致。修复：把从
    `resolve_side`（索引构建）到 merge/扩展的整个流程移入 `pool.install`，
    `--parallel` 现约束整个命令的 rayon 用量（`sd search --engine pgi` 与
    `rept e-align` 经由 `align pgi` 同步受益）。`-p 1/2/8` 输出逐字节一致
    （确定性未破坏）。

### 外部工具与参数 / CLI（3 处）

**lastz 静默失败**（只打日志返回 Ok）。修复：统计失败数并 bail。
**lastz 失败原因被吞**（status 丢 stderr）。修复：`cmd.output()` 记录
    首个失败的 stderr。
**参数校验缺失/不一致（align 侧）**：kmer/window/parallel 正值有限性。
    修复：统一校验，帮助同步。

### CLI / 文档（15 处）

**噪音与帮助文本多处小修**：lav mask stanza 静默、`#` 元数据行跳过、lastz
    `[multiple]`/`-s` 修正、align.md 示例输出修正、pgi 帮助默认 syncmer 修正。
**文档一致性（align 侧）**：`.pgi` 命名说明。
**lastz 单序列约束帮助/文档未同步**。修复：完整补齐。
**align-pgi.md `--freq` 语义错误**：写 "more than this many times"，代码与
    帮助均为 "at least this many times"（`>= freq`）→ 文档修正。
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
**align-pgi.md 明确 `--ref-seq` 校验范围**（contig 表）并要求序列与索引来源
    一致（自动 sibling 路径由 mtime 检查保证）；未实现 k-mer 内容校验
    （syncmer 哈希对比复杂度高、阈值易误报，文档说明足够）。
**align-pgi.md Notes 补充小写（软掩码）处理**：自动索引小写→N 无 seed/块，
    `pgr pgi build --mask` 同语义。
**`align pgi` after_help 的 `--freq` 语义写 "more than"**：帮助文本写
    "K-mers occurring more than --freq times on either side are skipped"，但
    代码（`emit_entry_hits`：`if ea_freq >= freq`）、`-f` 参数帮助与
    docs/align-pgi.md 均为 "at least"（`>=`，FastGA 语义）。"more than" 与代码/
    其余文档矛盾，会把用户引向 `> freq` 的错误预期。修复：改为 "at least
    --freq times"。
**`align pgi` after_help 使用不存在的 `--k` 长选项**：该参数实际注册为
    `.short('k').long("kmer")`（无 `--k` 别名），after_help 写 `--k/--smer/
    --window`，`--k` 会直接报 "unexpected argument"。修复：改用 `--kmer`（与
    `--smer`/`--window` 命名风格一致）。
**`sibling_pgi_path` 陈旧 doc 注释**：原写 "ref.fa / ref.fa.gz / ref.2bit all
    map to ref.pgi"，与 `.gz` 分离为 `ref.fa.pgi` 的实现不符。修复：更正注释
    （仅注释，无行为变化）。

## 验证

* 引擎交叉验证：pgi 与 lastz 检出同一对 1200 bp 拷贝（边界修剪 2 bp）；合成
  基因组上两引擎各 4 条命中覆盖相同两个重复家族（坐标差异仅边界修剪 4–8 bp）。
* 端到端坐标：`.2bit` 输入与 FASTA 路径逐字节一致；`--ref-seq` 单侧 / 双侧 /
  直接 FASTA 输入输出逐字节一致；`psl lift` 对 `>chr1:alt` 与 `>NC_000913.1`
  输出键完整、区间与真实坐标一致。
* 鲁棒性：截断/负跨度 lav、越界 PAF、垃圾 BED/.pgi、构造头、空输入、全 N、
  极值参数、短行、随机二进制喂各命令等畸形输入全部友好报错或空输出，零 panic。
* 确定性：`--parallel 1/2/8` 下输出逐字节一致；`align pgi` 反复运行逐字节
  一致。
* 数据安全：`align pgi -o` 同输入报 "also an input file" 且输入完好；`-o` 指向
  兄弟索引（`ref.pgi`/`ref.fa.pgi`）报 "output file ... is also an input file"
  且索引完好；`.pgi` sibling mtime 重建、`.gz` 与 `.fa` 兄弟索引分离、默认 k
  冲突报错均实测复现。
* 性能：`align pgi --parallel N` 现约束整个命令 rayon 用量；`-p 2` 端到端
  输出与修复前一致。
* 新增回归测试（主要）：`command_align_pgi_crafted_index_errors_not_panics`、
  `command_align_pgi_gz_sibling_index_distinct`、
  `command_align_pgi_stale_sibling_index_rebuilt`、
  `command_align_pgi_default_kmer_conflicts_with_cached_index`、
  `command_align_pgi_lowercase_copy_has_no_all_zero_blocks`、
  `command_align_pgi_single_ref_seq_on_self_pgi`、
  `command_align_pgi_output_not_overwrite_sibling_index`、
  `command_align_lastz_self_duplicate_basenames`、
  `command_align_lastz_omitted_query_is_self`、
  `index_records_match_sequence_positions`、
  `crafted_record_contig_rejected_not_panic`、
  `mmap_merge_rejects_out_of_range_contig`、
  `crafted_huge_contig_name_rejected_not_panic`、
  `merge_k64_high_key_no_prefix_overflow`、
  `parse_subrange_keeps_dotted_and_colon_contigs`、
  `test_lift_minus_strand_forward_coordinates`、
  `truncated_d_stanza_errors_not_panics`、`negative_span_l_line_rejected`、
  `extreme_l_line_values_do_not_panic`、`unbalanced_lengths_do_not_panic` 等。
* `align pgi` 22 个 CLI 测试与 `libs::pgi` 55 个单测全通过；`cargo test` 全量
  通过；本族 release 模式全绿（pgi、alignment 各 lib + 相关 CLI）；`cargo fmt
  --check` 与 `cargo clippy --all-targets -- -D warnings` 干净。

## 收尾复核（2026-08-05）

在前述收敛结论之上，对当前工作区代码与文档做最终回归核对，确认报告所述修复
均已落地且状态一致：

* 代码逐项核对：`cmd_pgr/align/pgi.rs`（`--self` 校验、`ensure_outfile_distinct`
  含 sibling 索引、`pool.install` 约束 `--parallel`、mtime / 参数一致性检查、
  `.gz` 分离命名）、`libs/pgi/build.rs`（`mask` 软掩码、pending 去重、按位置
  重算 key、`rc_key`）、`libs/pgi/align.rs`（`saturating_add` 防 k=64 前缀溢出、
  `freq >=` 双侧过滤、streaming 逐记录校验、负链 RC frame）、`libs/pgi/mmap.rs`
  （惰性解码 + 越界校验、截断拒绝）均与报告一致。
* 文档核对：`docs/align-pgi.md` 的 `--freq`("at least")、`--kmer`（非 `--k`）、
  `.gz` 分离命名、mtime 失效、缓存参数一致性、软掩码说明均与代码一致。
* 测试回归：`cargo test --lib libs::pgi` 64 通过；`cargo test --test
  cli_align_pgi` 23 通过（含 `crafted_index_errors_not_panics`、
  `output_not_overwrite_sibling_index`、`stale_sibling_index_rebuilt` 等）。

本轮未再发现新的代码/行为/CLI/文档问题，审核收敛状态得到确认。

## 结论

`align` 命令族审核完成（累计修复 41 处缺陷：26 处代码/行为 + 15 处 CLI/帮助/
文档），并经多轮纵深复核（`libs/pgi` 索引构建/读取、`libs/lastz`、
`libs/fmt/lav`、`libs/fmt/psl`、`alignment` DP、sibling/缓存索引新鲜度与
`-o` 覆盖保护、`--parallel` 确定性、`emit_entry_hits` 频率过滤与 k=64 前缀域、
`chain_tubes` 排序键布局、负链坐标、reference 侧记录校验）复核，未再发现新的
问题，审核收敛。