# fairy（fairy-prime）：FracMinHash 稀疏采样 + 宏基因组 coverage（源码分析）

> 2026-08 整理，基于本地 `fairy-prime/`（v0.5.8，2024 Microbiome，
> bluenote-1577/fairy，从 sylph fork）。功能：多样本宏基因组 MAG binning
> 的 contig coverage 计算，替代 all-to-all read alignment（BWA/minimap2），
> 声称快 100-1000×。对应 pgr 语境：`pgr fq norm` 大数据量方案调研中
> "reads 采样 vs k-mer 采样"的参考——fairy 是典型的 **k-mer 采样**路线。

## 1. 概况

- **定位**：`fairy sketch`（reads → 每样本 `.bcsp` 索引）+
  `fairy coverage`（contigs × 样本 → coverage 矩阵，MetaBAT2/MaxBin2/
  SemiBin2 兼容）。作者明确 caveat：不适于单样本 binning、PacBio HiFi、
  菌株分辨组装。
- **依赖**：needletail（FASTA/Q 解析）、rayon（并行）、fxhash、
  serde+bincode（`.bcsp` 落盘）、scalable_cuckoo_filter（近似去重）、
  statrs（泊松/伽马）、fastrand（bootstrap）、memory-stats（ram-barrier）；
  musl 静态编译时切 jemalloc。
- **入口**：`main.rs` 把 `coverage` 恒以 pseudotax=true 调用
  （`contain(args, true)`）→ 实际走 pseudotax 分支；sketch 的 dedup
  默认关闭（见 §8 quirk）。

## 2. 架构

| 模块 | 行数 | 内容 |
|---|---|---|
| `sketch.rs` | 853 | sketch 主流程、read dedup、序列化 |
| `contain.rs` | 1202 | coverage 主流程、ANI/λ 估计、矩阵输出 |
| `seeding.rs` | 209 | FracMinHash 标量实现（滚动 2-bit） |
| `avx2_seeding.rs` | 266 | 同上 AVX2 4 通道版（仅 x86_64） |
| `types.rs` | 189 | SequencesSketch / GenomeSketch / AniResult |
| `inference.rs` | 121 | λ 的矩估计/二分（nb 路径） |
| `cmdline.rs` | 126 | clap 参数 |

## 3. FracMinHash 采样（seeding.rs）

- **2-bit 编码**：A=0、C=1、G=2、T=3；反向互补 `nuc_r = 3 - nuc_f`
  （A↔T、C↔G）；canonical = `min(f, r)`。注意编码与 khmer
  （A=0,T=1,C=2,G=3）不同。**非 ACGT 碱基查 `BYTE_TO_SEQ` 表一律得 0(A)**
  ——含 N 的 read 会生成人为的 A-kmer。
- **哈希**：`mm_hash64` = murmur64 最终化变体（源码带
  `//TODO this is bugged. Fix after release` 注释，scalar 与 AVX2 版实际
  等价；bug 应指"对 canonical min(f,r) 直接哈希"这一层）。
- **采样**：`threshold = u64::MAX / c`；`hash < threshold` 才保留
  → 采样率 ≈ 1/c。默认 `c=50`（约 1/50，sylph 为 1/200）。
- **滚动**：f 左移 2 位累积、r 右移 + 顶部补补链，与 pgr 现有滚动同构；
  AVX2 版把 read 切成 4 段重叠窗口并行滚动（`extract_markers_avx2`），
  只支持 k=21/31（`2(k-1)=40/60` 硬编码）。

## 4. sketch：read → .bcsp

- **存储**：`SequencesSketch.kmer_counts: FxHashMap<u64, u32>` 全内存；
  落盘前转 `Vec<(u64, u32)>`（`SequencesSketchEncode`，序列化快一个量级）
  + 元数据（c、k、file_name、sample_name、paired、mean_read_length），
  bincode 写 `.bcsp`。每样本一张表，`threads`（默认 3）个样本并行。
- **去重（关键）**：
  - **pair marker**：固定 k=16（`Marker=u32`），把一对 read（或单端 read
    前半/后半）的奇偶位交错拼成两个 16-mer（`pair_kmer`/`pair_kmer_single`）。
    长度不足（双端 <33bp、单端 <66bp）或单端 >400bp → 无 marker、不去重。
  - **规则**：对每个采样的 kmer，若 `(km, marker)` 已见过 → 该 kmer
    不计数（`num_dup_removed++`）；否则插入并 `c += 1`。
  - **效果**：`c` ≈ 该 kmer 出现过的**不同 read-pair 数**——完全相同的
    read pair 只贡献一次；部分重叠的 pair 只对未见过 marker 的 kmer 计数。
  - 精确模式用 `FxHashSet`；`--fpr>0`（默认 0.0001，隐藏参数）切
    `ScalableCuckooFilter`。
- **单端上限**：`MAX_DEDUP_COUNT=4`，`c < 4` 才查去重（高拷贝序列的
  计数不再被门控）。双端无上限。
- 其余：`mean_read_length` 逐条 moving average；`MAX_DEDUP_LEN` 常量
  未使用。

## 5. coverage：contig × 样本

- **contig sketch**（`sketch_genome_individual`，每 contig 独立）：
  FracMinHash 采样 → 按 kmer 去重（重复出现的 kmer 直接丢弃，不保留）
  → `min_spacing=30`（间距 <30nt 的相邻 kmer 丢弃）→ `genome_kmers:
  Vec<u64>`；`gn_size` = contig 长度。可预存 `.bcdb`。
- **查询**：对每个 contig 遍历 `genome_kmers` 在样本表中查 multiplicity
  → `covs` 向量；`contain_count` = 命中数。过滤：`gn_kmers.len() ≥
  min_number_kmers`（默认 8）且 `covs` 非空。
- **ANI**：`naive_ani = (contain/total)^(1/k)`；随后用 λ 校正
  （`ani_from_lambda`：`contain/(1-e^-λ)/total` 再开 k 次方）。
  **输出阈值 0.95**（pseudotax 分支，main 恒走；普通分支为 0.9）
  ——与 wiki 的说法一致，但 0.9 分支实际不可达。
- **覆盖度估计**（`get_stats`）：
  1. `median_cov` = covs 中位数；median<30 时按 `Poisson(median)` CDF
     < 0.9999999999 剪掉高倍噪声（`max_cov`）；
  2. `full_covs` = 未命中补 0 + 命中的 covs（≤max_cov）；
  3. λ：默认 `ratio_lambda` = `count(mode+1)/count(mode) × (mode+1)`，
     要求 ≥25 个命中、mode+1 存在、两侧计数 ≥ `min_count_correct`（默认 3）；
     备选 `mme` / `nb`（矩估计二分）/ `mle`（零位 + Newton-Raphson）；
  4. `final_est_cov` = λ（可估）| median<15 时 `geq1_mean_cov` | median。
- **方差**：对 `full_covs` 前 95% 窗口算（`VAR_CUTOFF=10` 以下不剪）。
- **CI**：100 次 bootstrap（`fastrand::seed(7)` 固定），5-95 分位，
  <50 次成功则输出 NA。
- **pseudotax**（恒生效）：`winner_table` 把共享 kmer 按 ANI 最高的 genome
  重新分配，二次 `get_stats`；`estimate_true_cov` 再乘
  `read_length/(read_length-k+1)` 与 `1/(seq_id^k)` 校正。
- **输出**：默认 MetaBAT2 格式（contigName/contigLen/totalAvgDepth +
  每样本 cov/var），`--concoct-format`/`--aemb-format` 去方差列。

## 6. 内存与并行

- sketch：每样本一张 FxHashMap（论文：土壤样本约 4GB/样本；1/50 采样下
  内存 ≈ 采样 kmer 数 × (8+4) 字节 + 哈希表开销）；`--ram-barrier`
  （隐藏）是**软限制**：虚拟内存超限时 `sleep` 阻塞等回收，不保证上限。
- coverage：默认 `step=1` → 同一时刻只有 1 个样本 sketch 驻留
  （genome sketches 全量常驻）；pseudotax 时 `step = threads/2 + 1` 个样本
  并行。论文数字：10 样本索引 9min + 查询 7min vs BWA 40h。

## 7. 对 pgr 的启示

1. **fairy = "k-mer 采样"，不是 "reads 采样"**：它保留约 1/c 的 kmer 并
   只统计这些 kmer 的 multiplicity。pgr norm 的 bbnorm 语义（highpass
   filter，按 read 内 kmer 深度分位判定）**不需要采样**——采样会直接改变
   深度分位的语义；此前讨论的"对 reads 采样"也不是 fairy 路线。
2. **dedup 思路**（pair marker 指纹 + 按 kmer 门控计数）与 pgr 现有
   `fq clump` 精确整对去重是不同抽象层级；pgr 不需要引入。
3. 若 pgr 未来做 coverage/丰度类工具，FracMinHash + FxHashMap +
   `Vec<(u64,u32)>` bincode 是现成的最小实现模板；`c` 与内存线性反比。
4. 与 pgr 现有设施对照：

| 项 | pgr 现有 | khmer | fairy |
|---|---|---|---|
| 精确计数 | `KmerTable`（u128+u32）、`count.rs` .pkt sort-merge | 无 | 无 |
| 近似计数 | 无（bbnorm 精确表语义） | CMS（u8/u16 饱和） | FracMinHash 稀疏采样 + u32 multiplicity |
| 判定 | truedepth/depthAL 分位 + toss | median ≥ cutoff（在线） | 中位数 + 泊松剪枝 + λ 校正 |
| 内存 | 外部桶（mem_cap 约束） | 固定但随装载率失真 | 与采样率反比（1/c） |

## 8. 源码 quirks

- `--no-dedup` 默认 `true` 且无反开关（clap SetTrue 语义）→ **去重实际
  默认关闭**，与 README/论文描述的 illumina 去重不符，疑似 0.5.8 有意/无意
  改动。
- `mm_hash64` 带 "TODO this is bugged" 注释；AVX2 与 scalar 等价，
  问题在 canonical 编码层。
- 非 ACGT 一律编码为 A，含 N 的序列会产生人为 kmer。
- `pair_kmer` 的 k=16 与主 k=21/31 无关，仅作 read 指纹。
- coverage 经 main.rs 恒走 pseudotax 分支 → 默认 ANI 阈值实际是 0.95；
  `contain()` 里 0.9 的分支在 CLI 下不可达。
