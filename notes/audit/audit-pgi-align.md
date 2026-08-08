# align 命令族代码审核记录（2026-08-07）

对 `pgr align` 命令族（pgi/lastz + `libs/pgi`、`libs/lastz`、`libs/fmt/lav`、
`libs/fmt/psl`、`alignment` DP）约 3000 行代码及全部文档进行审核。以下仅保留
有借鉴意义的结论；验证过程已精简。

> 注：`pgr sd` / `pgr rept` 命令族的审核记录见 `notes/audit/audit-sd-rept.md`。
> `libs/pgi` 索引的**构建**与读取缺陷记录在本文件；sd 对 pgi 的**消费**缺陷
> 记录在 audit-sd-rept.md。

## 与外部参考实现的语义一致性核对

- `psl lift` 负链坐标提升与 UCSC kent `pslLiftSubrangeBlat.c` 的 `liftSide` 行为
  一致（子范围命名 pgr 1-based vs kent 0-based 为记录在案的有意差异）。
- pgi 索引 k-mer 频率过滤：`emit_entry_hits` 的 canonical key 过滤、
  `freq >= cutoff`（FastGA 语义，非 `>`）、前缀窗口 / 最大共享前缀 / 扩展范围过滤，
  与 FastGA GIX 语义一致。
- 软掩码语义：`build_from_seqs` 的 `mask`（小写→N，FastGA `-M` 语义）与
  `build_from_path` 一致；pgi 与 lastz 对小写拷贝行为完全一致。
- LAV→PSL 转换：`blocks_to_psl` 的 q/t 正 gap 计 insert、负 gap clamp 忽略、
  1-based→0-based `checked_sub` 防溢出、`-` 链坐标翻转，与 UCSC `lavToPsl` 一致。

## 排除的疑点（安全不变量，经核验无需修复）

- `chain_tubes` 的 u128 排序键字段布局：`anti`(≤2^33) 与 `bucket`(≤27b) 各占不足
  40 位、字段间无重叠，`a_contig`/`b_contig` 的 `u16` 截断因 `n_contigs <=
  u16::MAX` 约束而安全。
- reference 侧记录校验：`align pgi` 只走 streaming 路径，reference 经
  `PgiStream::next_batch` 逐记录 `validate_record`，`ac < contigs.len()`；非
  streaming 的 `merge_seed_hits` 其 `a` 来自 `PgiIndex::read`（逐记录校验）或
  `build_from_seqs`（构造即合法）。无 crafted reference 越界 panic。
- `emit_entry_hits` 频率过滤两侧对称：`a` 侧 `ea_freq >= freq` 丢弃、`b` 侧最大
  前缀/扩展范围均按 `>= freq` 处理（FastGA GIX 语义）。
- build `kmer_key_at` 切片越界排除：pending 位置均经 `start + k <= n` 守卫后才入队。
- Myers 波前反向波镜像数学与 `rt/rq` 的镜像坐标逐一推导一致；self 模式对角线 0
  边界在正/反向波两侧均正确钳制。自比对 banded DP 对角带取并集，终点对角线恒在
  `[k_lo, k_hi]`。
- TrimSpec 数值域：`bias` 落入 `[3,9]`、`mscore∈[207,300]`、`dscore∈[700,793]`、
  `score∈[-10500,4500]` 均落于 i16；全 N 参考经 `total==0→0.5` 兜底无除零。
- `tubes_for_group`/`extend_tube` 循环与坐标域：`anti`/`diag`/`bucket` 均在 i64 内
  无溢出；`cov` 覆盖计账仅在 `cps > ahgh` 时转 `u64`；`alow` 恒严格递增无死循环；
  `(None, None) => unreachable!()` 由 while 条件保证不可达。
- `extend_tube` 的 `eant`：`LocalAlign.t_end`/`q_end` 是 `usize`，`(t_end+q_end)
  as i64` 在 64 位上无现实溢出（此前"i32 溢出"担忧不成立，字段并非 i32）。
- `dedupe_contained`：`overlap_frac` 的 `own <= 0` 返回 0、`(ov.max(0))` 防负。

## 记录项（未改，低风险 / 待决策）

- `align lastz -o dir` 重复使用旧 LAV 残留：影响链短（`sd run`/`s-align`/`sd
  search lastz` 均用临时 workdir 免疫），且 LAV 是通用扩展名清理易误伤，记录不修。
- `.pgi` 显式输入 + 冲突 `-k`/`--smer`/`--window` 被静默忽略（exit 0）：
  docs/align-pgi.md 明确说明 `.pgi` 输入在索引头携带参数——文档化预期行为。
- crafted 索引可携带超大 contig `len`（u64 无上限），`as u32` 截断仅产生错误坐标
  不 panic；真实索引受"单 contig ≤ 4.3 Gb"已知限制约束。按"简洁优先"记录不修。
- `-f/--freq 0` 会因 `ea_freq >= 0` 恒成立而丢弃全部 reference 条目，静默输出空
  PSL。FastGA 默认 10，`0` 无意义取值，非崩溃非错误，记录不修。

## 已知限制（有意保留）

- 子范围命名 pgr 1-based vs kent 0-based：pgr 生成端/消费端自洽，直接消费
  UCSC/blat 生态子范围名时需先确认语义。
- 单 contig 过大：pgi 索引 pos 为 u32（>4.3 Gb 不支持），且 PSL 坐标字段为 32 位
  有符号（>~2.1 Gb 回绕）。均为格式固有上限（UCSC 亦同），真实基因组最大 contig
  ~250 Mb 远达不到。
- `ref.2bit` 与 `ref.fa` 同 stem 时共享 `ref.pgi` 兄弟索引（有意保留）：`.2bit`
  是 `.fa` 的压缩变换，作为 drop-in 替换共享索引符合设计。与 `.fa`/`.fa.gz` 的
  分离不同（两者可能是内容无关的独立文件）。

## 修复的缺陷（根因模式）

### 崩溃 / 越界 / 溢出（Zero Panic）

- **LAV 解析多处越界/下溢/溢出**：`d` stanza 边界差一越界（守卫改 `+ 6`）；`l`
  行负跨度回绕成超大 block（`t_end < t_start` 报 InvalidData）；`l` 行极值坐标
  `-1` 下溢/跨度比较溢出（`checked_sub`）。
- **构造 .pgi/.hv 头容量溢出 panic/OOM**（未校验 n_records/n_contigs）。修复：头
  解析校验 + `try_reserve_exact`。
- **pgi build `positions.len() as u32` 静默截断**（>42 亿记录）。修复：
  `payloads.len() <= u32::MAX` 防御检查。
- **`align_banded_local` 序列长度悬殊时 DP 数组越界**。修复：j_lo/j_hi 与对角带
  求交、空交集跳行。
- **crafted .pgi 记录 contig id 越界 panic**：`emit_entry_hits` 的 `b.contigs()[bc]`
  直接越界（读取路径只校验头部不校验记录体）。修复：`PgiIndex::read`/
  `PgiStream::next_batch` 逐记录校验 `cid < n_contigs` 且 `pos + k <= contig len`；
  `PgiMmap`（惰性解码）在解码命中记录时同步校验。
- **临时索引目录创建失败 panic**：`TempDir::new().expect` 在系统临时目录不可写/盘
  满时 panic。修复：改 `?` 传播友好错误。
- **`read_header` 按不受信任的 contig 名长度扩容**：`nb` 来自头部、无上限，可迫使
  数 GB 分配。修复：扩容前校验 `nb <= MAX_CONTIG_NAME`（1 MiB）。
- **`emit_entry_hits` 前缀窗口在 k=64 时 `hi` 溢出**：k=64 时 `mask=u128::MAX`，
  若某 k-mer 前 `len` 个碱基全为 T，`lo + r = 2^128` 溢出（`2^128` 无法表示）。
  修复：`hi` 改 `lo.saturating_add(r)`（饱和后仅排除全 T 的 `u128::MAX` 单键，非
  真实 seed）。
- **`build_from_seqs` 的 `k * 2` 参数校验溢出 panic**（极端 `k`）。修复：直接校验
  `k <= 64`。
- **`collect_one_contig` 的 `window + 2` 溢出 panic/OOM**。修复：`SyncmerParams::
  validate` 增加 `window <= 1_000_000` 上限。
- **`--parallel` 无范围校验导致 rayon 线程风暴**：`--parallel 18446744073709551615`
  让 rayon 创建近似无限线程（load 1000+），`--parallel 0` 被静默当默认值。修复：
  `RangedU64ValueParser::<usize>::new().range(1..=1024)`，clap 在构造任何线程池前
  拒绝（共享 helper，同时覆盖 lastz/sd/prefilter 等所有 `-p` 消费方）。
- **构造 .pgi 头大 `k` 使 `pack_kmer`/`rc_key` 移位溢出 panic**：`parse_header_bytes`
  未校验 `k <= 64`。修复：增加 `k in 1..=64` 校验（与 build 侧一致）。

### 功能正确性 / 算法

- **（重大）pgi 索引 k-mer key 与位置错配**（2 Mb 随机基因组 39% 错配、self 比对
  101 条伪块）。修复：pending 去重、flush 按位置重算 key、RC 用 `rc_key`。
- **align lastz 省略 query 未启用 self 模式**。修复：传 `self_mode`。
- **`psl lift` 负链外层坐标提升错误（违反 UCSC 约定）**。修复：
  `qStart/qEnd += start_0`、`qStarts += (size - end_0)`。
- **`psl lift` 的 `parse_subrange` 误切含 `.`/`:` 的 contig 名**（`NC_000913.1:1-200`
  被读成 name="NC_000913"、`chr1:alt:1-200` 被读成 chr="alt"），`lift_query` 查错键
  静默跳过。修复：取最后一个 `:`+数字后缀切分，前缀整体作 contig 名。
- **lastz self 模式用 basename 判断自比对，同名文件被交叉比对**（`a/dup.fa`、
  `b/dup.fa`）。修复：self 模式跳过所有 `target_file != query_file` 的作业。
- **`ref.fa` 与 `ref.fa.gz` 共享兄弟索引，内容不同时静默复用错误索引**。修复：
  `.gz` 输入去掉 `.gz` 后**追加** `.pgi`（`ref.fa.gz` → `ref.fa.pgi`），与
  `ref.fa` → `ref.pgi` 分离。
- **FASTA 原地修改后兄弟索引被静默复用**（同名单长但序列不同）。修复：`resolve_side`
  增加 mtime 校验（输入比索引新则重建）。
- **`.pgi` 单输入自比对 + 仅 `--ref-seq` 报错**。修复：`resolve_seqs` 后 self 模式
  下任一侧扩展序列为空时复用另一侧。
- **align pgi 自动索引小写归一化 → 全零块**：`build_from_seqs` 碱基编码大小写不敏感
  → 小写与大写拷贝共享 seed，但扩展 DP 大小写敏感 → 评分失败 → 回退全零块。修复：
  `build_from_seqs` 增加 `mask` 参数，align pgi 自动索引传 `true`（跳过小写）。
- **默认参数静默复用不同 k 的兄弟索引**：`resolve_side` 缓存参数冲突检查只覆盖
  命令行显式传的 `-k/--smer/--window`。修复：删除 `explicit(...)` 条件，**总是**
  检查当前解析值（显式或默认）与缓存索引的一致性。
- **`build_from_seqs` 在小 k 下产生重复 (kmer, pos, strand) 记录**（`k <= smer +
  window - 1` 时，位置在第二次选中前已被 flush 出队随后重新入队发射），虚增频率
  计数使 `--freq` 过滤误丢真实 seed。修复：分组环节按 payload 精确去重。
- **`emit_entry_hits` 的 lcp 收窄窗口仅含高频条目时漏种子**：FastGA 的 GIX 构建期
  就剔除 `>= freq` 高频 k-mer，其"收窄窗口只含高频条目"等价于"为空"必然回退；pgr
  保留全部 k-mer、merge 期过滤，收窄窗口可非空却只含高频条目，此时直接返回空漏掉
  低频种子。修复：`m < min_shared` 且 `start > min_shared` 时回退重扫 floor 窗口。

### 数据安全（`-o` 同输入 / 陈旧 / 静默数据丢失）

- **`align pgi` 的 `-o` 指向输入时静默覆盖输入**。修复：对 `-o` 及 `--ref-seq`/
  `--query-seq` 均加 `ensure_outfile_distinct`。
- **`align pgi` 的 `-o` 可静默覆盖兄弟索引**（`ref.fa` → `ref.pgi`）。修复：把每个
  基因组输入的 `sibling_pgi_path` 一并加入 `ensure_outfile_distinct`（跳过
  `stdin`）。

### 性能

- **`align pgi --parallel` 未约束自动索引构建的 rayon 并行度**：`resolve_side` 在
  自定义线程池创建前执行，索引构建走全局 rayon 池，`--parallel` 只约束 merge/扩展
  阶段（文档承诺约束整个命令）。修复：把从 `resolve_side` 到 merge/扩展的整个流程
  移入 `pool.install`（`sd search --engine pgi` 与 `rept e-align` 同步受益）。

### 外部工具与参数 / CLI

- **`--self` / 省略 query 的 self 模式未校验 `--ref-seq`/`--query-seq` 一致**（可在
  self 模式把两个不同文件交叉比对）。修复：扩展序列一致性校验移到 `self_mode` 分支
  （单输入与 `--self` 统一生效）：两侧都给出必须相同；仅给一侧时省略 query 的 self
  模式复用该侧，显式 `--self` 仍要求两侧都给。
- **lastz 静默失败**（只打日志返回 Ok）。修复：统计失败数并 bail；记录首个失败
  stderr。
- **参数校验缺失/不一致（align 侧）**：kmer/window/parallel 正值有限性统一校验。

### CLI / 文档（一次性小修，已精简）

`--freq` 语义（"more than"→"at least"，`>=` FastGA 语义）、`--k`→`--kmer` 长选项、
`pgr psl to_chain`→`to-chain` 命令名、`--min-shared` 示例措辞、`--lastz-args` `=`
写法提示等均修正；`docs/align-pgi.md` 的 `.pgi` 命名、mtime 失效、缓存参数一致性、
软掩码说明、`--merge-gap`/`--max-gap` 语义补齐。

## 结论

`align` 命令族审核完成，累计修复 **52 处缺陷**（34 处代码/行为 + 18 处 CLI/帮助/
文档），经多轮纵深复核未再发现新问题，审核收敛。
