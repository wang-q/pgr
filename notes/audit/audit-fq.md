# pgr fq 命令族代码审核记录（2026-08-11）

对 `pgr fq`（10 个子命令）命令族及相关库文件（`libs/fq/{clump,norm,pairs,
sample,split,trim,trim_adapter}`、`libs/fmt/fq`、`libs/fmt/seq`、`libs/loc`、
`libs/par`）和全部测试/文档进行审核。本报告由 `audit-fa-fq-2bit.md` 中的 fq
部分拆分并扩写，覆盖 fq 全命令族。以下仅保留有借鉴意义的结论；验证过程已精简。

## 命令族

`fq` 现有 10 个子命令：`interleave`(`il`)、`to-fa`、`clump`、`split`、
`sample`、`clean`、`filter`、`norm`、`range`、`trim-qual`。

## 与外部参考实现的语义一致性核对

fq 家族以 BBTools 39.38 为主参考，逐字节核对（golden 数据见
`tests/bbtools/Lambda/README.md`），另有 `fa range` 复用的 `.loc` 索引与
`sickle`/`cutadapt` 风格的质量修剪。均已核验一致：

- `clump` → `clumpify.sh`（k=31, seed=1 默认输出顺序，逐字节）。
- `split` → `repair.sh` `rp` 模式（R1/R2/singles 逐字节）。
- `sample` → `reformat.sh samplebasestarget sampleseed`（精确模式，去重下无
  上行采样）。
- `clean`/`filter` → `bbduk.sh`（anchr 管线两趟调用，`ordered=t` 确定性）。
- `norm` → `bbnorm.sh passes=1 bits=16 min=<n> target=9999999` 的读决策逻辑
  （计数用精确 canonical 表，非近似 `bits=16` 哈希表）。
- `trim-qual` 的 sliding/mott → `sickle`/`cutadapt` 的质量修剪算法。
- `range` → 复用 `fa range` 的 `.loc` 索引语义（name、plain offset、record
  size），BGZF/plain 支持一致。

有意差异（均已记录）：
- `norm` 计数用精确表而非 `bits=16` 近似计数，深度边界附近的读可能因计数精度
  不同而被取舍（文档已注明）。
- `clump` 外部 bucket 路径输出为"按桶拼接"序，与 BBTools 大数据行为一致但与
  全局内存序不同（文档已注明）。

## 排除的疑点（安全不变量，经核验无需修复）

- `fq is_fq` 对目录输入：`File::open` 成功但 `read_exact` 失败（EISDIR），返回
  友好错误而非 panic（Zero-Panic）。
- `sample` 的 `score_cov`/`read_stats` 索引：`cov` 非空后，`above_limit` 回退到
  -1 时分支短路、`depth_al` 保持 -1，不会越界索引。
- `norm` 外部路径 `load_table` 置 `k: 0`：`read_stats`/`canonical_keys` 均用
  `opts.k`，不依赖表的 `k` 字段。
- `trim.rs` `sliding_cut`/`mott_cut`：窗口尺寸至少 1、空序列提前返回，无除零/
  越界。
- `trim_adapter.rs` 的 `make_codes`/`antialias`/`JavaRandom`：与 JVM 序列化
  逐项有单元测试锁定。
- `fq range` 的 LRU 缓存：`LruCache<String, Vec<u8>>` 用 `&str` 借钥命中，
  `Borrow` 关系成立。

## 记录项（未改，低风险 / 待决策）

- `fq range` 的 `name:0-N`（start=0）被当作"整条记录"而非子序列（与 `fa range`
  的 `start==0` 约定一致）；未文档化该边界，但不返回错误数据（仅取整条）。
- `norm` 输出会对保留的读施加 `changequality` 归一化（N 质 0、ACGT 质最低 2），
  与 bbnorm 默认一致；文档已说明 changequality 被应用。

## 已知限制（有意保留）

- `interleave` 双文件格式只按 `infiles[0]` 探测（`is_fq` 只看第一个文件），若
  第二个文件格式不同会在实际读取时报错（非静默）。
- `trim-qual` 配对质量编码只从 R1 自动检测（文档已注明）。

## 修复的缺陷（根因模式）

按根因模式分组；标注"承袭"者自 `audit-fa-fq-2bit.md` 沿袭，其余为本次审核新发现。

### Zero-Panic / clap 参数缺失

- **`split` 缺少 `--outfile-2` 时 panic**：`args.get_one::<String>("outfile_2")
  .unwrap()` 对 `None` 解包崩溃。修复：`outfile_2` 参数加 `.required(true)`，
  clap 先于执行校验，缺失时输出友好用法错误。新增回归测试
  `command_fq_split_missing_outfile2_is_clap_error`。
- **`sample` 缺少 `--bases` 时 panic**：`args.get_one::<i64>("bases").unwrap()`
  对 `None` 解包崩溃。修复：`bases` 参数加 `.required(true)`。新增回归测试
  `command_fq_sample_missing_bases_is_clap_error`。
- **`clean`/`filter`/`norm`/`clump` 的 `--parallel` 未做 1..=1024 范围校验**：
  四个命令的 `--parallel` 均接受 `auto` 或整数，但只 `parse::<usize>()` 未限上
  界，数值直接进入 `rayon::ThreadPoolBuilder::num_threads` 创建线程池
  （clean/filter 经 `par::ordered_map`、norm/clump 直接建池），越界值（如
  `--parallel 1000000`）会尝试创建海量线程，导致系统资源耗尽。违反全局硬约束
  "`--parallel` 必须经 clap 校验 1..=1024"。修复：新增共享助手
  `cmd_pgr::args::parse_parallel_auto`（`auto` 取逻辑 CPU 数，整数须在
  1..=1024，越界友好报错），四处统一改用。新增回归测试
  `command_fq_{clean,filter,norm,clump}_parallel_out_of_range_*`。
- **`clump` 的 `--buckets` 未校验范围，`--buckets 0` 触发除零 panic**：
  `clump_buckets` 中 `key.kmer as u64 % buckets as u64`（`libs/fq/clump.rs`）
  对 `--buckets 0` 除零崩溃。修复：`execute` 中校验 `--buckets` 须在
  1..=4096，越界友好报错（与 clump 内部 `MAX_BUCKETS=4096` 一致）。新增回归
  测试 `command_fq_clump_buckets_out_of_range_is_friendly_error`。
- **`sample` 输入以空（0 碱基）记录结尾时除零 panic**：`sample` 循环里
  `remaining` 每轮按记录碱基数递减，当输入末尾存在空记录时 `remaining` 先减到
  0，下一轮 `target / remaining`（`libs/fq/sample.rs`）除零崩溃。修复：在除法
  前若 `remaining == 0`（说明剩余全为空记录）直接 `break`。新增回归测试
  `command_fq_sample_trailing_empty_records_do_not_panic`。

### 数据安全（`-o` 同输入保护）

- **流式命令允许 `-o` 覆盖输入文件**：`fq to-fa`/`fq interleave` 已统一加入
  `ensure_outfile_distinct`。

### 输入校验 / 静默错误

- **`interleave` 双文件交错对读取计数不匹配静默截断**（`zip` 取较短者）。修复：
  `interleave_read` 中任一文件先读完而另一未读完即 `bail!`。

### 行为一致性 / 算法

- **`interleave` 单文件虚拟 R2 两路径不一致**：单 FQ→FA 为 `"\n"`（空序列）、单
  FA→FA 为 `"N"`；帮助与 `docs/fq.md` 均声明 "N's"。修复：统一为 `b"N"`。
- **`interleave` 双文件路径返回的最终索引错误**：更新后的 `idx` 被丢弃，最终返回
  未递增的 `start`，违背 pub fn 契约。修复：两文件分支改为
  `idx = interleave_read(..)?`。

### 文档一致性

- **`trim-qual` 命令名错写为 `trim-q`**：`trim_qual.rs` 帮助文本/示例、
  `libs/fq/trim.rs` 的 `TrimOptions` doc、`docs/fq.md` 均改为 `trim-qual`。
- **`docs/fq.md` 子命令清单不完整**：补全 10 个子命令列表。
- **`docs/fq.md` 缺失 `trim-qual` 小节**：新增完整 Options/Examples。
- **`clean` 文档 gzip 输出示例误导**（`io::writer` 写端不压缩）：
  `-o unmerged.trim.fq.gz` 改为 `-o unmerged.trim.fq` 并仅指明输入可为 gzipped。
- **`to-fa` 文档误置于 `norm` 小节**：移回其所属小节。

### 死代码 / 功能不可达

- **`clean` 的 `--mask-kmers` 静默失效（死代码）**：`trim_adapter.rs` 中按
  `ktrim_right` 分派"ktrim / kmask / filter"三个分支，但 `clean.rs` 把
  `ktrim_right` 硬编码为 `true`，使 kmask 分支永远不可达。于是文档化的
  `--mask-kmers`、`--mask-fully-covered`、`--trim-pad` 三个选项被静默忽略，
  掩码功能完全失效（`filter` 用 `ktrim_right: false` 走 filter 分支，不受影响）。
  修复：
  - `ktrim_right` 改为 `kmask.is_none()`——默认保持 ktrim=right（与 bbduk
    逐字节 golden 行为不变），指定 `--mask-kmers` 时切到 kmask 掩码分支。
  - 新增守卫：`--mask-fully-covered` / `--trim-pad`（仅掩码语义）在未给
    `--mask-kmers` 时报友好错误，避免静默无效；`--mask-kmers` 在未给 `--ref`
    时报友好错误（与文档"requires --ref"一致），避免静默无操作。
  - 新增回归测试 `command_fq_clean_kmask_masks_instead_of_trims`（掩码为 N、
    全长保留 vs 默认 ktrim 截短）、
    `command_fq_clean_kmask_mask_only_options_require_mask_kmers`（友好错误）、
    `command_fq_clean_kmask_requires_ref`（缺 --ref 友好错误）。
  `docs/fq.md` 中 `--mask-kmers`/`--mask-fully-covered`/`--trim-pad` 的描述
  现与行为一致。

## 第二轮：merge 家族（merge / ec-kmer / ec-overlap / extend / s-filter）审核

首轮报告完成后，代码库新增了 anchr merge 管线的 5 个子命令及其库（`libs/fq/{merge,
overlap,bbnet}`、`libs/asm/{tadpole,assemble}`、`libs/kmer/{quality,qcheck}`）。
本报告拆分后这部分未被首轮覆盖，故新增第二轮审核。逐命令/库核对 Zero-Panic、clap
参数、数据安全（`-o` 同输入保护）、算法边界与文档一致性。

### 排除的疑点（安全不变量，经核验无需修复）

- `merge`/`ec-overlap` 的 `process_pair` 对 `seq<2` 短读提前返回 `Ambiguous`，
  `join_reads` 的 `overlap==0` 分支切片起点恒 ≤ 终点（`insert` 已由
  `min_insert>0` 保证），`corrected_pair` 的 `rc2.len()-c2len` 由
  `c2len ≤ b.len()` 保证不欠位；`expected_tip_errors` 中 `bases`/`quals` 同长。
- `bbnet::parse` 对 23 维输入层、非空层、逐层权重维度做了校验，
  `feed_forward` 的最后 `last()/first()` 均有保障，畸形 `.bbnet` 返回友好错误
  而非越界。
- `tadpole::run` 对 `k∈[1,128]` 校验；`count_read_kmers` 对 phred 索引
  `min(127)` 限位；`extend_read` 的 rollback 用 `saturating_add(1)` 防 `%0` 与
  `usize::MAX+1` 溢出。
- `libs/kmer/quality.rs` 对 `k > MAX_K` 提前返回空表；`qcheck::check_read` 对
  短读/无锚点返回 `NoAnchor` 而非 panic。

### 修复的缺陷

- **`s-filter` `-k ≥ 65` 触发 u128 移位越界 panic**：`qcheck` 的锚点/延伸扫描用
  u128 滚动 k-mer（每碱基 2 bit，最多 64 碱基），`masks` 返回 `rc_top=2k-2`，对
  `k≥65` 时 `1u128 << rc_top` 等移位在 debug 下 panic
  （实测 `pgr fq s-filter -k 65` 崩）。修复：`execute` 校验 `k∈1..=64` 并报友好
  错误。新增回归测试 `command_fq_s_filter_kmer_out_of_range_is_friendly_error`。
- **`s-filter` 的 `--discard-file` 可覆盖输入文件**（数据安全）：仅 `-o` 受
  `ensure_outfile_distinct` 保护，`--discard-file` 指向输入会破坏原文件。修复：
  对 `discard_file` 也加 `ensure_outfile_distinct`。新增回归测试
  `command_fq_s_filter_discard_file_same_as_input_rejected`。
- **`ec-overlap` 的 `--efilter 0` 未真正禁用**（行为不一致）：`merge` 用
  `(*x>0.0).then_some(*x)`（0 置 `None` 禁用），`ec-overlap` 却 `Some(*x)`，
  使 `--efilter 0` 在 `process_pair` 中 `(expected+offset)*0 < bad` 恒真，意外
  抑制 pfilter。修复：`ec-overlap` 与 `merge` 对齐。
- **`bbnet` 单元测试与新增 23 维校验冲突**：`parse_and_forward_minimal` 用旧的
  2 维输入网，被 `dims[0]==23` 校验拒绝后 `unwrap()` panic。修复：测试改用 23 维
  输入（隐藏层补足 23 权重），并新增 `parse_rejects_wrong_input_dims` 用例。

### 文档一致性

- **`docs/fq.md` 的 `ec-kmer`/`extend` `-k` 标注 "no upper bound" 不实**：实际
  `tadpole::run` 限制 `k≤128`。改为 "up to 128, the k-mer key table limit"。
- **`docs/fq.md` 的 `s-filter` `-k` 未注明上限**：标注 `1..=64, the u128
  rolling-key limit`，与新校验一致。
- **`docs/fq.md` `ec-overlap` 与 `extend` 之间缺 `---` 分隔**：补全，保持小节间
  分隔一致。

## 第三轮：迁移重构后复核 + 未覆盖路径补查

期间外部并行完成了 anchr 迁移阶段 1（`notes/design/fq-asm-migrate.md`）：
`detect_quality_base`/`PHRED33`/`PHRED64` 从 `libs/fq/trim.rs` 抽到基础层
`libs/fq/qual.rs`，`kmer qhist`/`fq s-filter`/`trim` 改引用新位置；kmer 的
`base_codes`、`count::count_keys` 由 `pub(crate)` 改 `pub`；新增
`tests/migrate_api.rs`（外部 crate 视角编译证明 anchr 依赖的基础符号可达）。

复核结论（均绿）：
- 重构未引入新问题：`cargo build`/`clippy --all-targets -D warnings`/`fmt --check`
  clean，全量 `cargo test` 通过（含 `migrate_api` 与既有 fq/asm 测试）。
- 无残留旧路径引用：`grep trim::detect_quality_base|trim::PHRED` 无命中；
  `detect_quality_base` 的 `seq[i]` 因 `SeqReader` 强制 seq/qual 等长而安全。
- 补查此前未逐行覆盖的路径，确认校验齐全：
  - `clump`：`-k` 2..=31（`fill_max` 的 `x2 << (2k-2)` i64 移位上限）、
    `--buckets` 1..=4096、`--parallel` 走 `parse_parallel_auto`；
  - `norm`：`-k` 2..=31、`--parallel` 走 `parse_parallel_auto`；
  - `sample`：`remaining==0` 提前 break 防除零，`--bases` 必填；
  - `s-filter`：`-k` 1..=64（第二轮修复）、`--discard-file` 防覆盖输入。
  未发现新缺陷。

## 结论

首轮：`fq`（10 子命令）命令族合计修复本次缺陷（Zero-Panic / clap 参数缺失 4：
`split --outfile-2`、`sample --bases`、`clump --buckets 0` 除零、
`sample` 空记录结尾除零；`--parallel` 范围校验 1 [clean/filter/norm/clump
四处]；死代码/功能不可达 1 [clean `--mask-kmers`]），均含回归测试；另有承袭自
`audit-fa-fq-2bit.md` 的数据安全 1、输入校验 1、行为一致性/算法 2、文档一致性 5。
全部 fq CLI 集成测试与 fq 库单元测试通过，`cargo clippy -- -D warnings` 与
`cargo fmt --check` clean。

第二轮：merge 家族（`merge`/`ec-kmer`/`ec-overlap`/`extend`/`s-filter` 5 个子命令
及 `libs/fq/{merge,overlap,bbnet}`、`libs/asm/tadpole`、`libs/kmer/{quality,qcheck}`
）修复缺陷 4（Zero-Panic 1 [s-filter `-k≥65` 移位越界]、数据安全 1 [s-filter
`--discard-file` 覆盖输入]、行为一致性 1 [ec-overlap `--efilter 0`]、测试修复 1
[bbnet 23 维校验与单元测试对齐]），文档一致性 3，均含回归测试或已修正文档。

第三轮：anchr 迁移阶段 1 重构（`detect_quality_base` 抽到 `qual.rs`、kmer pub 化）
复核全绿，无残留旧引用；补查 `clump`/`norm`/`sample`/`s-filter` 校验齐全，未发现
新缺陷。

经多轮纵深复审收敛，未再发现新问题。
