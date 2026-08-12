# asm 命令族代码审核记录（2026-08-12）

对 `pgr asm` 全部 7 个子命令（contig/unitig/ovlp/layout/cns/olc/map）及相关库
文件（`libs/asm/assemble.rs`、`libs/asm/tadpole.rs`、`libs/olc/overlap.rs`、
`libs/olc/layout.rs`、`libs/olc/consensus.rs`、`libs/map.rs`，以及
`cmd_pgr/asm/common.rs`）和全部测试/文档进行审核。以下仅保留有借鉴意义的结论；
逐轮验证过程已精简（第一轮发现并修复 5 类问题，第二、三轮复核未再发现新问题）。

| 子命令 | 命令文件 | 主要库 |
| ------ | -------- | ------ |
| contig | `src/cmd_pgr/asm/contig.rs` | `libs/asm/assemble.rs`, `libs/asm/tadpole.rs` |
| unitig | `src/cmd_pgr/asm/unitig.rs` | `libs/asm/assemble.rs`, `libs/asm/tadpole.rs` |
| ovlp   | `src/cmd_pgr/asm/ovlp.rs`   | `libs/olc/overlap.rs` |
| layout | `src/cmd_pgr/asm/layout.rs` | `libs/olc/layout.rs` |
| cns    | `src/cmd_pgr/asm/cns.rs`    | `libs/olc/consensus.rs` |
| olc    | `src/cmd_pgr/asm/olc.rs`    | 上述全部 |
| map    | `src/cmd_pgr/asm/map.rs`    | `libs/map.rs` |

共享：`src/cmd_pgr/asm/common.rs`。

## 与外部参考实现的语义一致性核对

asm 家族对照 BBTools tadpole（contig）、BCALM2/GATB `ograph.cpp graph3`
（unitig）、BBTools bbmap `perfectmode`（map）语义，OLC 三段各自对照参考流程，
并逐字节/逐行对拍验证。已确认的一致/有意行为：
- `contig`：k-mer 计数质量门控、种子深度阈值、双向贪心扩展、泡泡消除，复现
  tadpole 默认 contig 模式（`popbubbles=t`）。
- `unitig`：仅当唯一后继且其前驱唯一时扩展，停于分支/接点/覆盖缺口/环，与 k
  扫描顺序无关（BCALM2 压缩语义）。
- `map`：种子-验证、完美匹配（无错配无缺口）、`ambiguous=all`（单端）、
  `--paired` 携带配对标志/伴侣坐标/TLEN，作为 bbmap `perfectmode maxindel=0`
  的替代。
- OLC：unitig → ovlp（精确重叠）→ layout（贪心互惠接点）→ cns（精确重叠拼接）
  三段闭环。
- 确定性与复现：种子扫描按规范化 k-mer 排序快照线性迭代（不依赖哈希表内存序），
  多线程输出保序。

## 排除的疑点（安全不变量，经核验无需修复）

- `find_overlaps`：seed_k 先 `min(最短 unitig)` 再 `ensure!(seed_k >= 1)`，
  空/超短 unitig 先归 0 被友好错误拦截，不会在 `query_unitig` 的
  `seq[0..seed_k]` 越界 panic；seed_k 已钳到最短长度，所有 unitig 窗口切片有效。
- `verify_seed`/`extend`：四个方向扩展均有边界守卫，`q[p..p+k]`/`t[tpos..tpos+k]`
  不越界；`query_unitig` 的 tid/tpos 来自排序后 payload，索引安全。
- `build_layouts`/`extend`：`placed` 保证不死循环；`2*u+end` 索引在
  `qid/tid < n` 保证下安全；`is_repeat` 的 `length*10` 在 64 位 usize 下不溢出。
- `consensus`/`dedup_contained`：`overlap <= seq.len() && overlap <= piece.len()`
  防切片越界；`contains` 对空 needle 有守卫。
- `make_contig`：`bb.len() == k` 守卫在种子处即拒绝 contig；扩展循环
  （`extend_to_right`/`right_counts_of`/`left_counts_of`/`calc_coverage`）与
  tadpole 移植一致；`resolved_min_contig_len` 的 `2*k` 因 k ≤ 128 不溢出。
- `ovlp` 的 `--overlap-k 0`/`--min-overlap 0`：seed_k 经 `.max(1)` 归 1；
  `min_overlap=0` 仅放行所有 found overlap，均无 panic。
- `map`：单端/配对统计、SAM FLAG/TLEN/RNEXT 计算正确；`--max-reads` 单端跨文件
  与配对成对语义均正确。
- `common.rs`：`to_paf`/`read_unitigs`/`format_cov` 边界安全。

## 记录项（未改，低风险 / 待决策）

- `contig`/`unitig` 对**单个文件**按"交错配对"解释（1 交错或 2 配对），单端 reads
  需用 2 个文件或显式交错文件；奇数条 reads 时报 `unpaired trailing read`。与文档
  一致，属设计语义，未改。
- `map` 配对模式仅报告两端**首个**命中（`ambiguous=all` 语义仅单端成立），属配对
  模式固有约束，未改。

## 已知限制（有意保留）

- `map --max-reads` 奇数时按"成对计两条"取整到偶数（配对原子性）。
- `contig`/`unitig` 的 `-p/--parallel` 仅为兼容而接受，实际单遍确定性处理
  （文档已注明）。

## 修复的缺陷（根因模式）

### Zero-Panic / 溢出

- **`-k` 超上限 panic**（contig/unitig/olc）：只校验 `k >= 1` 未校验上界，
  `-k 129`（> `key::Kmer::MAX_K = 128`）在 `count_read_kmers` → `Kmer::new(k)
  .expect()` 处 panic。修复：`assemble()` 与 `assemble_unitigs_core()` 加
  `ensure!(k <= MAX_K)` 友好报错。回归 `command_asm_contig_rejects_kmer_above_
  limit` / `command_asm_unitig_rejects_kmer_above_limit`。
- **`cns` 布局 `contig_0` 下溢 panic**：`parse_index(...)? - 1` 对 1-based 布局
  id 直接减一，`contig_0` 使 usize 下溢（debug panic / release 环绕成
  `usize::MAX` 触发 OOM）。修复：先校验 `>= 1` 再减一。回归
  `command_asm_cns_rejects_contig_zero`。

### 数据安全（`-o` 同输入保护 / 输出截断）

- **`contig`/`unitig` 的 `-o` 可覆盖输入**：writer 在读取 reads 前打开，`-o`
  指向输入会先截断 reads。修复：打开 writer 前 `ensure_outfile_distinct`。回归
  `command_asm_contig_outfile_not_input` / `command_asm_unitig_outfile_not_input`。
- **`map` 的 `--outm`/`--outu` 可覆盖输入**：SAM writer 在流式读取 reads
  **之前**打开并截断，`--outm`/`--outu` 指向 reads 文件会丢数据。修复：对
  `outm`/`outu` 与 ref+reads 全部 `ensure_outfile_distinct`。回归
  `command_asm_map_outm_not_input`。
- **一致性补充**：`ovlp`/`layout`/`cns`/`olc` 四个写 `-o` 的子命令虽在打开
  writer 前已把输入读入内存（无中途丢失），为对齐项目统一约定也补上
  `ensure_outfile_distinct`（分别对 unitig FASTA / PAF+FASTA / TSV+FASTA /
  reads）。

### 文档一致性 / 帮助文本

- **`-k` 上界描述错误**：`docs/asm.md` contig/unitig 段写 "no upper bound"，
  map 段与 `map` help 写 "1..=64"，与代码实际上限 128 不符。修复：contig/unitig
  段改 "up to 128, the k-mer key table limit — k > 64 uses multi-word k-mers"；
  map 段与 help 改 "1..=128"。

## 实机边界验证（均无 panic，返回友好错误或正常退出）

- `contig` reads 长度 < k → 正常退出，0 contigs
- `unitig` 空输入文件 → 正常退出，0 unitigs
- `ovlp` 空 unitigs → 友好错误 `cannot overlap: unitigs are empty...`
- `map -k 129`（超上限） → 友好错误 `k-mer length must be in 1..=128`
- `map --paired` 两端 reads 数目不等 → 友好错误 `paired files have different read counts`
- `map` 空 reads → 正常退出，0 mapped

## 结论

`asm` 命令族审核完成（累计修复 5 类问题 + 4 处 `-o` 防护统一），经纵深复核收敛；
与 BBTools/BCALM 语义对拍、边界输入验证零 panic，`cargo fmt`/`clippy` 干净，
相关集成测试与 `cargo test --lib` 全部通过。
