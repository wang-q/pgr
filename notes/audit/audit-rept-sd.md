# sd / rept 命令族代码审核记录（2026-08-05）

对 `pgr sd`（8 命令 + `libs/sd`）与 `pgr rept`（6 命令 + `libs/pl`）两个命令族
约 5000 行代码及全部文档进行审核。两个命令族结构相近（同为重复序列分析工作流，
共享 `libs/pgi` 索引消费端与 `libs/loc` 索引），合并为一份记录。以下仅保留有
借鉴意义的结论；验证过程已精简。

> 注：`pgr align` 命令族的审核记录见 `notes/audit/audit-pgi-align.md`。`libs/pgi`
> 索引的**构建**缺陷记录在 audit-pgi-align.md；本文件记录 sd 对 pgi 的**消费**
> 缺陷。

## 与外部参考实现的语义一致性核对

- greedy 链合并：与 FastGA `align_contigs` / `ALNchain.c` 的链化语义一致——同
  对角线纯间隔是两条独立链，仅对角线平移才缝合（pgr 的自身扩展）。
- `sd run` 的 chainnet：每靶位点保留一条最优链（同等分按序取一）是 UCSC chainnet
  的标准语义（与 `pl chainnet` 共享）。

## 排除的疑点（安全不变量，经核验无需修复）

- `sd run` 的 cluster set_id 重编号值域各簇两两不相交，不可能碰撞。
- sd cluster minus 链序列提取：按 pgr PAF 正向坐标约定提取，逐碱基一致。
- wave 初始 trim 越界经几何推演与约 20 万次 fuzz 均不可达，不加防御。
- 全量扫描家族生产代码 `unwrap()`/`unreachable!`：`loc.rs` fields[0..2]、
  `cover.rs` f[0..7]、`repeat.rs` fields[0..1] 等索引访问前均有长度检查；
  `merged.last_mut().expect` 有非空前置保证。
- 全部 HashMap 用法逐一核对：pl/repeat 与 trf 的 name/safe_map 仅查找；cover 的
  by_set 经 set_order（Vec）迭代；decompose 的 index/kmer_frags 顺序无关；仅
  cluster 的 by_root 曾依赖 HashMap 迭代序（已修复、按首个区间排序）。
- 双引擎差异逐层归因：search 层面 pgi 232 块全部与 lastz 282 块重叠 ≥50%；
  差异传导至 decompose 的 cluster 划分。BISER 语义允许两引擎输出不同（可互换替代
  引擎），各自自洽且坐标正确，非缺陷。
- tube 工作流"库 vs 基因组"结构性失效（重测确认**无需修复**）：greedy 已移除、
  tube 为唯一流程后，MG1655 vs TnCentral 库 `rept e-align` 正常检出，失效随
  syncmer/排序键修复消失。

## 记录项（未改，低风险 / 待决策）

- `syncmer.rs` 参考实现与 `collect_one_contig` 重复发射同一位置，消费方已
  HashSet 去重，可后续合并。
- `fa split name` 名称碰撞（`chr(1)` 与 `chr_1`）概率极低。
- `rept e-kmer`/`s-kmer` 的 `--fill-kmer` 以 `usize as i32` 传入 `IntSpan::fill`；
  超 i32 值静默截断为负 → fill 变 no-op（`excise` 同理安全）。`s-align` 的
  `--min-depth` 以 `usize as u32` 截断。极端参数属用户误用，记录不修。
- `sd cluster` 同染色体重叠合并只按 chrom 名分组：跨基因组 PAF 在两端 contig 名
  与文件 stem 均相同时会把不同基因组的同名区间并簇。`sd cluster` 文档仅面向自
  比对 PAF，记录不修。
- 顶层路径为 `.pgi` 扩展名的目录会被 `is_pgi_input` 误判拒绝（概率极低）。
- `save_repeat_cache` 中途失败（`.ktab` 已写、`.complete` 未写）时残留文件会被
  `cache_is_fresh` 判陈旧自动重建；但若"保存中途失败 + 后续重建 part 数减少"
  同时成立，旧高序号 part 文件可能残留被 FastK `-p:` 读到。概率极低，记录不修。
- `sd search --engine pgi` 接受 `.2bit` 输入但下游 `sd align`/`sd run` 需 FASTA，
  2bit 在 `fa size` 步骤报错。文档仅承诺 FASTA；2bit 部分支持是既有行为。
- `--max-gap` 调大（如 10000）时，greedy 循环的 off-band 忽略规则会把后续不同
  对角线种子全部忽略，远距重复家族可能整体丢失。off-band-ignore 是既有设计
  （默认 1000 下正确），记录不修。
- `rept e-align` 传入 `.2bit`：在 `has_soft_mask` 的 FASTA 读取器处报 "stream did
  not contain valid UTF-8"。文档仅承诺 FASTA；有错误提示的非静默失败。
- TRF 外部工具限制：完美 2500 bp 周期 + max-period ≥ 2600、`--max-period` 10000+
  均触发 TRF SIGSEGV，pgr 将信号错误友好传播。精确上限未知，无法可靠预校验。
- TRF 版本兼容：本机 TRF 4.09 `-ngs` 输出 17 字段，`parse_trf_output` 的 ≥15
  字段门槛兼容；`@chr1` 头行（1 字段）跳过。
- 只有头的 FASTA（`>chr1` 无序列）：`rept s-kmer` 触发 FastK SIGSEGV——外部工具对
  空序列崩溃，cmd_lib 捕获报 "terminated by signal: 11"，pgr 自身无 panic。
- 纯四联体重复（如 ACGT）只有 4 种不同 10-mer，低于 `MIN_SHARED_KMERS=5` 阈值，
  同源片段不会合并为同一 set_id——极端低复杂度序列，非 SD 场景，行为符合设计。
- s-kmer 尾 run 保守丢弃：Profex `-z` 从不闭合 read 的最后一个 run，s-kmer 按
  设计保守丢弃尾部（与 repeat.rs 文档一致）。

## 已知限制（有意保留）

- s-kmer 对染色体尾部重复保守丢弃：Profex `-z` 不输出末 run 深度，有阈值时无法
  区分唯一尾与重复尾（与 anchr 参考管线一致）。

## 修复的缺陷（根因模式）

### 崩溃 / 越界 / 溢出（Zero Panic）

- **sd/run.rs 解析 elem.bed 短行越界**（直接取 `f[4]`）。修复：加 `f.len() < 8`
  检查。
- **sd decompose 负链投影 usize 下溢**（畸形 header）。修复：拒绝 end < start，
  投影 saturating。
- **非 UTF-8 临时目录路径 `to_str().unwrap()` panic**（sd run 临时目录）。修复：
  `io::path_to_str` 友好报错。
- **e-align span 过滤 `(t_end - t_start) as usize` 回绕**。修复：i64 运算
  `.max(0)` 再转 usize。

### 功能正确性 / 算法（含 2 处重大链算法缺陷）

- **（重大）tube 排序键 anti/bucket 溢出**（>8 Mb 基因组失效）。修复：anti/bucket
  扩到 32 位。
- **（重大）tube 排序键负对角线回绕**（>64 Mb 间距失效）。修复：`BUCK_OFF = 1 << 26`。
- **pgi 引擎灵敏度限制**（记录项升级）：精确 k-mer seed 对近 90–93% identity 或
  真长恰在 `--min-len` 附近的拷贝可能只锚定子块被滤。已解决：`sd search` 默认
  `freq=50/k=31` 后，E. coli 漏检率 13.1%→0.26%。
- **cluster 重叠 union 漏连嵌套区间**。修复：扫描时跟踪最大右端。
- **sd cluster 去重键忽略链向/物种**（回文倒位拷贝被折叠）。修复：键加 strand。
- **sd cluster/run 不支持普通 gzip**（生成垃圾 `.loc`）。修复：非 BGZF 先解压到
  临时文件。
- **greedy 链合并导致倒位 SD 漏检**。修复：合并条件加 `|diagA − diagB| > 0`。
- **pgi merge 频率过滤两侧边界不一致**（`== freq` 处理与 FastGA 不符）。修复：
  A/B 侧统一 `>= freq` 跳过、`< freq` 计入。
- **相邻链合并把两条独立同源对缝成嵌合链，SD 命中丢失**：多拷贝家族两条拷贝对的
  对角线差与轴间隔都在 band/merge_gap 内时，纯几何 merge 无法区分"同源块种子缺口"
  与"两条独立块"（形状一致）。修复：`merge_adjacent_chains` 增加可选序列参数，
  两侧间隔均非空时要求中段 banded 对齐身份 ≥ 0.9 且 query 覆盖 ≥ 0.9 才合并。
- **`sd run` 合并 elementary BED 按 `read_dir` 顺序枚举 cluster 文件**：set_id 全局
  重编号与输出行序依赖文件系统枚举顺序，跨运行不确定。修复：按 cluster 文件名的
  **数值**编号排序（词法排序会把 cluster_10 排在 cluster_2 前）。
- **`sd align`（`chainnet_to_paf`）/ `sd search --engine lastz` 按 `read_dir`
  顺序迭代 MAF/LAV 文件**，输出行序不确定。修复：排序后再合并/转换。
- **`search_lastz::decompress_if_gz` 解压命名碰撞**：统一命名 `{base}.plain.fa`
  （base 取首个 `.` 段），嵌套目录中同名 `.fa.gz` 会解压到同一路径静默覆盖。
  修复：同一次调用内 HashSet 去重，重复 basename 追加输入序号后缀。
- **倒位拷贝间隔 < max_gap 时 greedy 链循环把两条互惠链并成嵌合链，SD 完全漏检**：
  两条互惠链（同一倒位对的两个方向）落在同一对角线上，种子在 greedy 循环内直接
  连成一条链（绕过 `|diag|>0` 守卫），扩展出嵌合块身份被稀释过滤。修复：greedy
  循环在"双侧种子间隙 ≥ 200 bp"时用中段同源检查门控（不通过则闭合当前链、以该
  种子起新链）；间隙 < 200 bp 不检查（对 ≥1000 bp 块，200 bp 随机间隙的嵌合身份
  ≥ 0.909 高于 SD 阈值，不会静默漏检）。
- **`sd cluster` 的 cluster 编号依赖 HashMap 迭代顺序，`sd run` 的 set_id 编号
  跨运行不稳定**（同一基因组多次运行输出 set_id/行序互换）。修复：按每个分组的
  首个区间（chrom, start）排序后再编号。
- **tube 工作流在显式同文件对（非 --self）时把家族交叉命中当"重复"丢弃**：精确自
  比对巨块在 `dedupe_contained` 中把坐标上包含于其内的拷贝对块误判为重复丢弃。
  修复：dedupe 增加跨度相近约束（前块跨度 ≤ 后块 4 倍才判重复）。
- **`.loc` 索引陈旧时静默使用（open_indexed 只查存在性）**：同长度内容修改会静默
  提取错误序列。修复：`open_indexed` 增加 mtime 新鲜度校验（`.loc` 的 mtime 早于
  FASTA 时自动重建）。`fa range`/`sd cluster`/`fas check`/`get_seq_loc` 四个调用方
  同步受益。
- **`sd run --engine lastz` 输出重复 elementary 行**：lastz 互反块坐标抖动使单个
  cluster 文件出现 end 差 1 bp 的两个头，decompose 投影到相同 elementary 区间 →
  合并后重复。修复：`sd run` 合并层按 renumber 后的完整行去重（`push_unique_elem`）。
- **s-align 漏做带点 contig 名映射**（spanr 截断，`fa mask` 失配）。修复：复用
  chrom.sizes 映射。
- **Profex `-z` 坐标右端多 +1 + e-kmer 染色体尾部丢失**。修复：end 不再 +1；无
  阈值时用染色体长度闭合尾 run。
- **s-align/e-align soft-mask 警告误报 N gap**。修复：`has_soft_mask` 只扫
  lowercase。
- **`s-align` 安全名改写仍用 `split_once(':')`**：`chr1:alt:1-200` 被切成
  name="chr1"，带 `:` 的 contig 输出键被截断成 "alt"。修复：改用 `parse_subrange`
  解析再写占位名。

### 输入校验 / 静默错误

- **decompose 对解析失败的 FASTA 头静默丢弃**。修复：补 `log::warn!`。
- **`sd search`/`sd cross`/`sd run` 传入 `.pgi` 索引**：pgi 引擎对 `.pgi` 不做
  扩展（无序列），输出块全部 0 分，SD 过滤后静默返回空结果。修复：`pgi_to_hits`/
  `lastz_to_hits` 前置拒绝 `.pgi` 输入。
- **`sd align` 跳过非 2 组件的 MAF 块时无提示**。修复：补 `log::warn!`。
- **repeat.rs 两处 `map_while(Result::ok)` 吞 IO 错误**。修复：`let line = line?;`
  传播错误。
- **e-align PSL 过滤静默跳过畸形行**。修复：补 `log::warn!`。
- **空 FASTA 输入触发 FastK SIGSEGV（预检友好报错）**：空 repeat 库喂 `rept
  e-kmer` → FastK SIGSEGV，报错信息像 pgr 自身崩溃。修复：`run_repeat_pipeline`
  在 FastK 前预检输入是否有非空序列，空则报友好错误。

### 数据安全（`-o` 同输入保护 / 陈旧索引 / 静默数据丢失）

- **sd 命令 `-o` 指向输入文件时静默覆盖输入**。修复：`sd search`/`align`/`cover`/
  `decompose`/`cross` 均加 `ensure_outfile_distinct`。
- **`sd cluster` 输出目录残留旧 cluster 文件**：向含 `cluster_1.fa`/`cluster_2.fa`
  的目录重跑 `sd cluster`（本次仅 1 个 cluster）时旧文件残留，下游会静默消费陈旧
  家族。修复：写输出前清理 outdir 中 pgr 自身命名模式的 `cluster_<u32>.fa`。
- **rept 命令 `-o` 指向输入文件时静默覆盖输入**。修复：rept 五个子命令（含库输入）
  均加 `ensure_outfile_distinct`。
- **损坏的 FastK 缓存被静默复用 → e-kmer 空输出**：`cache_is_fresh` 只检查存在 +
  mtime，截断的 `.ktab`（`.complete` 和 part 完好）被静默读取 → 空 runlist。
  修复：增加 `.ktab` 与 `.complete` 大小一致性校验，不一致即视为陈旧重建。

### 性能

- **`run_lastz` self 模式 n×n job 列表**：self 模式只构建对角 n 个 job（不再生成
  n²），执行期防御保留。

### 外部工具与参数 / CLI

- **参数校验缺失/不一致（sd 侧）**：`--min-identity` 范围、minscore 正值有限性
  统一校验，帮助同步 "(0, 1]"。
- **sd search/cross `--preset` 默认值未注册**。修复：`.default_value("set01")`。
- **sd run --engine lastz --preset 拼装错误**。修复：`Vec<String>` + `$[preset_args]`
  展开。
- **trf 特殊字符文件名找不到**。修复：`sanitize_filename(chr)`。

### CLI / 文档（一次性小修，已精简）

主帮助与文档补齐 e-align/s-align、`sd run` 的 `--preset` 默认值/`--min-identity`
范围、e-align identity 定义（gap-compressed）、sd.md 输出去重与 pgi 引擎灵敏度
限制、软掩码语义（修正：实测 lastz 大小写不敏感，仅 pgi 感知小写；mask 仅影响
`sd search` 发现阶段）等说明。

## 结论

`sd`/`rept` 两个命令族审核完成（累计修复 51 处缺陷：41 处代码/行为 + 10 处
CLI/帮助/文档），经多轮纵深复核收敛，未再发现新问题。pgi 索引构建侧缺陷见
audit-pgi-align.md。
