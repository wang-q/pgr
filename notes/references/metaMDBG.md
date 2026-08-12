# metaMDBG（1.4）：minimizer-space de Bruijn 图宏基因组组装器（源码分析）

> 2026-08 整理，纯源码分析（`metaMDBG-metaMDBG-1.4/`，版本 `1.4`）。metaMDBG 是
> 面向**长而准的宏基因组 reads**（PacBio HiFi、Nanopore R10.4+）的组装器，论文
> [High-quality metagenome assembly from long accurate reads with metaMDBG]
> (Nature Biotechnology 2023)，作者 Gaëtan Benoit、Rayan Chikhi、Christopher
> Quince 等。**与 pgr 的关系**：它是 rust-mdbg（minimizer-space DBG）的宏基因组
> 工程化版本，与 pgr 的 `asm unitig`（bcalm 移植）同属 k-mer/unitig 路线，但其
> 核心创新——**local progressive abundance filter（用丰度替代气泡解析处理菌株
> 多样性）**——正好回应 pgr 讨论中"气泡不如不处理"的直觉。
> **与 OLC 的连接（2026-08-12）**：pgr `asm olc` 已落地，metaMDBG 的
> 渐进丰度过滤与 RepeatRemover 直接映射其 v1 的"覆盖度证据 repeat
> breaking"（见 §9）。

## 1. 概况

- **定位**：从长 reads（HiFi/ONT R10）构建宏基因组 contigs，输出带
  `length= coverage= circular=` 头信息的 FASTA；可选输出 GFA 组装图。
- **两种数据模式**：`getParamsHifi()`（ONT 模式默认开 read correction、HiFi
  默认关）；nanoMDBG（ONT R10 simplex）方法已集成进 1.4。
- **与 BCALM 的关系**：论文作者含 Rayan Chikhi（bcalm 作者），但实现**不基于
  BCALM/GATB**——BCALM 的 k-mer 是碱基 k-mer，metaMDBG 的节点是 **k′-min-mer
  （连续 minimizer 序列）**，走 rust-mdbg 的 minimizer-space DBG 路线。
- **语言/构建**：C++20 + OpenMP，CMake；内嵌 `ext/`（minimap2、spoa、htslib、
  TurboPFor 整数压缩），**全部 vendored，无外部运行时依赖**。
- **命令形态**：单一二进制 + 子命令（`asm` 为主，`graph`/`contig`/`toMinspace`/
  `readSelection`/`readCorrection`/`toBasespace`/`gfa` 等为底层步，见 §3）。
- **checkpoint 断点续跑**：每步写 `tmp/checkpoints/<step>.checkpoint` 文件，
  重跑同一命令自动跳过已完成步骤（README 明示，`AssemblyPipeline.hpp` 里
  `createCheckpoint`/`isCheckpoint` 实现）。

## 2. 仓库结构

```
metaMDBG-metaMDBG-1.4/
├── src/
│   ├── MdbgAssembler.cpp        # 主入口：子命令分发
│   ├── Commons.hpp              # 全局类型/常量/工具（8378 行）
│   ├── pipeline/AssemblyPipeline.hpp   # ★ asm 主流水线（multi-k′ 调度）
│   ├── graph/
│   │   ├── CreateMdbg.{cpp,hpp}        # ★ k′-min-mer 计数 + MDBG 建图
│   │   ├── ProgressiveAbundanceFilter.hpp  # ★ 图简化 + 渐进丰度过滤
│   │   ├── Graph.hpp              # UnitigGraph2 数据结构
│   │   ├── GenerateContigGraph.hpp # final 轮的 repeat solver（嵌合清理）
│   │   ├── GraphPOA.hpp / GraphSimplify.hpp / GenerateGfa.hpp / GfaParser.hpp
│   ├── assembly/GenerateContigs.hpp   # ★ unitig 路径 → contigs.nodepath
│   ├── toBasespace/
│   │   ├── ToMinspace.hpp        # contig nodepath → minimizer 序列（反馈下一轮）
│   │   ├── ToBasespace2.hpp      # ★ minimizer contig → 碱基序列（minimap2+POA 抛光）
│   │   ├── RepeatRemover.hpp / OverlapRemover2.hpp / DerepSmallContigs.hpp
│   │   ├── ContigPolisher.hpp / ContigTrimmer.hpp / ReadVsContigMapper.hpp
│   ├── readSelection/
│   │   ├── ReadSelection.hpp     # reads → minimizer 表示（含 MinimizerParser 入口）
│   │   ├── ReadCorrection.hpp    # ONT 的 reads 纠错（minimap2 内嵌）
│   ├── contigFeatures/           # KminmerCounter/KmerCounter（当前未参与主流程）
│   └── utils/                    # args.hxx / BooPHF / BloomFilter / edlib / MurmurHash3
└── ext/                          # vendored: minimap2, spoa, htslib, TurboPFor
```

> 源码总量约 8.3 万行（`wc -l src/*`），大头是 `Commons.hpp`、`CreateMdbg.cpp/hpp`、
> `ProgressiveAbundanceFilter.hpp`、`ToBasespace2.hpp`、`ReadCorrection.hpp`。
> 文件多、注释少（大量被注释的旧代码），核心算法集中在 ProgressiveAbundanceFilter
> 与 ToBasespace2。

## 3. 命令入口（`MdbgAssembler.cpp`）

| 命令 | 实现 | 作用 |
|---|---|---|
| `asm` | `AssemblyPipeline` | 完整组装（唯一面向用户的命令） |
| `graph` | `CreateMdbg` | 计数 k′-min-mer，建 unitig 图 |
| `contig` | `GenerateContigs` | 图简化（superbubble/tip/丰度过滤），导出 nodepath |
| `toMinspace` | `ToMinspace` | nodepath → minimizer 序列 |
| `readSelection` | `ReadSelection` | reads → minimizer 表示 |
| `readCorrection` | `ReadCorrection` | reads 纠错（ONT） |
| `toBasespace` | `ToBasespace2` | minimizer contig → 碱基序列 + 抛光 |
| `toBasespaceGfa` / `gfa` | `ToBasespaceGfa` / `GenerateGfa` | GFA 输出 |
| `derepSmall` / `removeOverlaps` / `removeRepeats` | 对应类 | 最终后处理 |
| `map` | `MappingContigToGraph` | contig → 图映射（调试用） |

## 4. asm 主流水线（`AssemblyPipeline::execute_pipeline`）

`execute()` 先跑 `convertReadsToMinimizerSpace()`，再跑多轮
`executePass(k, prevK, pass)`；**最终后处理（`derepSmallContigs → removeOverlaps →
removeRepeats → toBasespace`）不是循环结束后的独立步骤，而是在最后一轮
`executePass(_lastK, ...)` 的 `isFinalPass` 分支内联执行**（`AssemblyPipeline.hpp:1111`
起），输出 `contigs.fasta.gz`。首轮/中间轮只生成 unitig 序列反馈下一轮，不做碱基重建。

### 4.1 multi-k′ 迭代

```cpp
_firstK = 4;
_lastK = Commons::computeLastK(_minimizerDensityAssembly, readStats._n50ReadLength, _firstK, _maxK);
// lastK = N50ReadLength * density * 2.0（10 kb 读长、density 0.005 → lastK ≈ 100）
// 每轮 k += Commons::getMultikStep(k)  // 1.4 里恒为 1
// 最后一轮再显式补一次 executePass(_lastK, ...)
```

每轮 `executePass(k)`：

1. `createGraph(k, pass)`：`graph` 子命令。**首轮**（pass 0）从
   `read_data_corrected.txt` 计数全部 k′-min-mer（`--min-abundance`，默认 0 =
   rescue 模式）；**后续轮**加载上一轮 refined abundance
   （`loadRefinedAbundances`），并从 `read_data_corrected.txt` **和**
   `unitig_data.txt`（上一轮的 unitig 序列）一起计数——这是
   "unitig 反馈进下一轮"的实现。
2. `generateContigs(k, pass)`：`contig` 子命令。加载 unitig 图 →
   `ProgressiveAbundanceFilter::execute`（图简化，见 §6）→ `generateContigs3`
   从各 cutoff 快照生成 `contigs.nodepath`。
3. `toMinspaceContigs(...)`：把 nodepath 转回 minimizer 序列，写入
   `unitig_data.txt`（非 final）或 `contig_data_init.txt`（final）。
4. 每轮结束时 `savePassData(k)`；非首轮 `dumpUnitigAbundances` 备份
   `unitigGraph_prev.*` 并写 refined abundance。

**k′ 的语义与长度换算**（`AssemblyPipeline::writeParameters`,
`AssemblyPipeline.hpp:1479`）：multi-k 里的 `k` 是**一个 k′-min-mer 包含的
minimizer 个数**，不是碱基长度。换算关系写在 `parameters.gz` 里，供各子命令
`Parameters::load` 复读：

```cpp
minimizerSpacingMean = 1 / assemblyDensity;   // 相邻 minimizer 的平均间距(碱基)
kminmerLengthMean   = minimizerSpacingMean * (k-1);   // 一个 k-min-mer 的期望碱基跨度
kminmerOverlapMean  = kminmerLengthMean - minimizerSpacingMean; // 相邻 k-min-mer 重叠
```

故 `k` 每 +1，k′-min-mer 期望碱基长度约 +`1/density`（assembly density 0.005 →
每轮约 +200 bp）。这就是"unitig 反馈 + k 递增"实现**渐进长单元化**
（longer k-min-mer → 更多直链、更少分支）的机制。

**assembly graph 导出节奏**：`--gen-graph` 默认在第 11 轮（`_nextGenGraphIteration=11`，
之后每 +10）导出一次 GFA（`doesGenerateAssemblyGraph`，
`AssemblyPipeline.hpp:831`；首轮过大不导出）——用"隔轮导出"控制磁盘/内存，
pgr 的 `pl` 管道若多轮组装可参考。

### 4.2 最终后处理（isFinalPass）

- `derepSmallContigs`：去小 contig 重复（`DerepSmallContigs`）。
- `removeOverlaps`：去除 contig 间 overlap（`OverlapRemover2`）。
- `removeRepeats`：`ReadVsContigMapper` 把 reads 映射回 contig，找未桥接的
  重复位点并断开（`RepeatRemover`）。
- `toBasespace`：`ToBasespace2`，见 §7。

## 5. minimizer 提取（`ReadSelection` + `Kmer.hpp::MinimizerParser`）

### 5.1 编码与同聚体压缩

- 2-bit 编码（`DnaBitset`），`KmerModel` 滚动 k-mer；`EncoderRLE` 对 HiFi 做
  homopolymer 压缩（ONT 关）。
- `MinimizerParser(_minimizerSize, _minimizerDensity, ...)`：默认
  `--kmer-size 15`（cap ≤ 16），assembly density 0.005、correction 0.025。

### 5.2 采样规则（FracMinHash 式"通用 minimizer"）

```cpp
_minimizerBound = minimizerDensity * maxHashValue;   // u_int64
// 对每个 k-mer：
u_int64_t kmerHash = MurmurHash3_x64_128(&kmerValue, sizeof(kmerValue), 42);
if(kmerHash < _minimizerBound){ minimizers.push_back(kmerValue); ... }
```

即**不是窗口内取最小 k-mer**，而是对每个 k-mer 哈希、保留低于阈值者——等价于
rust-mdbg 的 universal minimizer / FracMinHash 采样。这样 minimizer 是无序空间
均匀采样，密度 ≈ 0.5%（assembly）/2.5%（correction）。

> 与 pgr 的对照：pgr `asm map` 用全量 k-mer 索引（讨论过 minimizer/syncmer 但
> 结论是不优化）；metaMDBG 的 minimizer 密度极低（0.5%），是**组装图节点**，
> 不是比对种子，两者目的不同。

## 6. 建图 + 图简化（`CreateMdbg` + `ProgressiveAbundanceFilter`）

### 6.1 k′-min-mer 计数（`CreateMdbg::createMDBG`）

- 分区数 `_nbPartitions = nbBases / 20Gb`，clamp 到 `[nbCores, 5000]`。
- `KminmerCounter`：把每条 read 的连续 minimizer 切成 k′-min-mer（`KmerVec`，
  canonical normalize，`hash128()` = MurmurHash3 128-bit），按 `hash % nbPartitions`
  写分区文件 → 分区内去重计数 → 合并写 `kminmerData_min.txt` +
  `kminmerData_abundance.txt`。
- **rescue 机制**（`--min-abundance 0` 默认）：首轮计数后把 `abundance==1` 的
  singleton k′-min-mer 单独"营救"一遍——对每条 read，若其 k-min-mer 中**多数是
  solid（丰度 >1）、少数是 singleton**，则这些 singleton 很可能是低覆盖基因组
  的真实 k-mer 而非测序错误，予以保留（`RescueKminmerFunctor`：read 内丰度中位
  数的 10% 为阈值，`median*0.1 <= 1` 才营救）。
- 节点表 `MdbgNodeMapLight`：`phmap::parallel_flat_hash_map<u_int128_t,
  DbgNodeLight>`（10 分区 + mutex）。
- 之后 `indexEdges → computeUnitigNodes → computeDeterministicUnitigs →
  indexUnitigEdges`，输出 `unitigGraph.nodes.bin` / `edges.successors.bin` /
  `nodes.abundances.bin` / `stats.bin`。

### 6.2 UnitigGraph2（`Graph.hpp`）

- `UnitigNode`：`_unitigName`、`_successors_forward/reverse`、`_nbMinimizers`、
  `_abundance`（float）、`_abundances`（每个组成 k-min-mer 的丰度向量）、
  `_unitigMerge`。
- **丰度语义**：unitig 丰度 = 组成 k-min-mer 丰度向量的**中位数**
  （`computeMedianAbundance`）；`recompact` 合并两个 unitig 时把两个丰度向量
  merge 后再取中位。
- 长度估计 `getLength = (nbMinimizers-1) * _minimizerSpacingMean`（minimizer
  空间里没有碱基坐标，长度是期望值）。

### 6.3 渐进丰度过滤（`ProgressiveAbundanceFilter`）★

`execute` → `simplifyProgressive(functor)` 主循环：

```cpp
maxAbundance = min(图内最大丰度, 10000);
currentCutoff = 0;
while(true){
    isModification = simplify();          // superbubble + tip + repeat solver
    checkSaveState(currentCutoff);        // 每个新 cutoff 存一张图快照
    nbErrorRemoved = removeAbundanceNoQueue(maxAbundance, currentCutoff);
    if(!isModification && !nbErrorRemoved) break;
}
```

**`simplify()`**（图结构简化）：

- `SuperbubbleRemoverOld`：找 `nbSuccessors>1` 的节点做 BFS 找出口，`isSuperbubble`
  判定，`collapseSuperbubble2` 收集并删除低丰度分支（丰度高于
  `currentCutoff/0.25` 的受 repeat solver 保护不移除），邻接 unitig `recompact`。
- `TipRemover`：按 `_nbMinimizers` 升序队列删 tip。
- final 轮额外挂 `_repeatSolver`（`GenerateContigGraph`）做嵌合 unitig 清理。

**`removeAbundanceNoQueue`**（丰度渐进过滤，核心）：

```cpp
float t = 1.1;
while(t < abundanceCutoff_min){
    currentCutoff = t;
    移除所有 abundance < t 的 unitig（记录前驱/后继）；
    对受影响的邻接 unitig 排序后 recompact（合并成更长 unitig）；
    t = t * (1 + 0.1);            // 每次约 +10%
    increaseStep = min(t_new - t, 10);
}
```

即**从丰度 1.1x 起步、按 ~10% 的步长逐步抬升阈值，每次移除低于阈值的 unitig
并重新压实图**。这是论文里 "local progressive abundance filter" 的实现：它不
试图区分菌株气泡里的正确路径，而是用丰度把低覆盖分支逐级删掉，让高丰度物种
的主路径自然收敛——**与 pgr 讨论中"气泡经常引入不确定性、不如不处理"的直觉
一致**。

**cutoff 快照（`dumpUnitigs`）**：每个新 cutoff 把当前 unitig 图导出到
`filter/unitigs_<idx>.bin`（node path + circular/repeat 标记 + 丰度），
`_cutoffIndexes` 记录 `{idx, cutoff}`。供 `generateContigs3` 从高 cutoff 到低
cutoff 倒序消费（见下）。

## 7. contig 生成 + 碱基空间重建

### 7.1 `GenerateContigs::generateContigs3`

- 从 `_cutoffIndexes` **最后一个（最高 cutoff）倒序**读快照：
  `_minUnitigAbundance = cutoff / 0.5`，跳过 `contigAbundance < _minUnitigAbundance`
  的路径、已组装过的 unitig（`_processedNodeNames`）、final 轮的重复 unitig。
- 圆形 contig 特殊处理（`isCircular` 标记，`nbMinimizers += 1` 以闭合）。
- 输出 `contigs.nodepath`（unitig index 路径）+ `_nodeNameAbundances` →
  `unitigGraph.nodes.refined_abundances.bin`（供下一轮复用）。

### 7.2 `ToBasespace2`（minimizer → 碱基）

- 内嵌 minimap2（`mm_dbg_flag |= MM_DBG_NO_KALLOC`），HiFi `map-hifi`/`ava-pb`、
  ONT `map-ont`/`ava-ont`。
- 流程：reads 映射到 minimizer contigs → `partitionReads`（按内存分片，
  `minimapBatchSize = peakMemory/8`）→ `createBaseContigs` 每片读回 reads
  用 POA（内嵌 spoa）抛光 → 输出 `contigs.fasta.gz`（header
  `ctg<id> length= coverage= circular=yes|no`）。
- `--skip-correction` / `--min-contig-length 50` / `--min-contig-coverage 1`
  可过滤输出。

## 8. 与 pgr 的对应/借鉴点

1. **丰度过滤替代气泡解析**（最值得借鉴）：`ProgressiveAbundanceFilter` 的
   "1.1x 起步、~10% 步长、边删边压实"策略可以直接映射到 pgr `asm contig`/
   `unitig` 的 `--min-coverage` 语义：不是单阈值一刀切，而是**多轮渐进 + 每轮
   重压实**，低覆盖菌株分支逐步被吞并。pgr 目前只有全局 `--min-coverage`，
   可考虑加"渐进模式"。
2. **cutoff 快照倒序输出**：metaMDBG 存多个 cutoff 的图快照、生成 contig 时
   从高丰度往低丰度补——天然适合"先出高置信 contig、再补低丰度"的宏基因组
   输出策略，pgr 的 `asm contig --min-coverage` 是单值，可借鉴快照思路。
3. **unitig 丰度 = 中位数向量**：与 pgr `asm unitig`（bcalm 移植）的
   `km:f:` 平均丰度不同，metaMDBG 保留每个 unitig 的丰度向量并取中位数，
   merge 时合并向量。若 pgr unitig 要输出稳健丰度（宏基因组场景），可参考
   中位数语义。
4. **k′-min-mer = minimizer 序列**：pgr 的 kmer 表是碱基 k-mer（u128 ≤ 64），
   metaMDBG 的节点是"minimizer 序列"（k′ 个 minimizer 的向量，hash128）。
   两者维度不同：minimizer 空间天然支持长读（HiFi/ONT），pgr 目前是短读工具，
   这一块暂不对齐，但知道差距在哪。
5. **内嵌 minimap2+spoa 抛光**：pgr 的 `asm map` 是完美匹配、无 gap；metaMDBG
   的 toBasespace 用 minimap2 容忍错误 + POA 抛光。若 pgr 未来要支持长读纠错
   或容错比对，metaMDBG 是"内嵌依赖"的参考，但 pgr 目前不引新依赖（用户约束）。
6. **断点续跑**：checkpoint 文件机制简单实用，pgr 的 `pl` 管道若做多步任务
   可参考（不过 pgr 目前坚持原语路线，优先级低）。
7. **multi-k 反馈 = 迭代式参数精化**（新增，2026-08-12）：把"组装"重构成
   **"参数化子命令 + 磁盘中间文件 + 循环调度"**——同一 `graph`/`contig`/
   `toMinspace` 子命令被 `AssemblyPipeline` 以不同 `k` 反复调用，跨进程只通过
   `parameters.gz`（gzip 二进制参数 blob，`Parameters::load/save`）和
   `unitig_data.txt`/`refined_abundances` 传递状态。pgr 的单进程 `libs/` 路线
   不必照搬子进程，但**"迭代长度参数 + 反馈 unitig"** 的骨架可直接映射到
   `asm` 的多趟 OLC/unitig 循环；`parameters.gz` 可类比 pgr 用 struct 传参。
8. **外部分区计数（scale-out）**：k′-min-mer 计数不把全量 k-mer 塞内存，而是
   `hash128 % nbPartitions`（`nbBases/20Gb`，clamp `[nbCores, 5000]`）写分区文件
   → 分区内去重计数 → 合并（`KminmerCounter::partitionKminmers`，
   `CreateMdbg.hpp:3652`）。这是典型的"外排序式"大数据手法，pgr 若做超大
   数据集（如 `kmer count` 溢出内存）可参考分区+归并，而非一味加大内存。
9. **内存驱动的批量分片**：toBasespace 用 `--max-memory`（`peakMemory/8`，
   clamp `[1,100]`）决定 minimap2 一次读入多少 reads（`ToBasespace2.hpp:337`）——
   峰值内存预算显式控制批大小。pgr 若加长读抛光，可把内存预算作为一等参数。

> 结论：metaMDBG 对 pgr 的最大价值是**§6.3 的渐进丰度过滤**——它给出了
> "不解析菌株气泡、用丰度逐级简化"的成熟实现，正好验证用户对气泡处理的直觉；
> 其次是 unitig 丰度中位数语义与多 cutoff 快照输出。其余（minimizer-space、
> minimap2 抛光）与 pgr 短读+完美匹配的路线距离较远，暂不借鉴。

## 9. OLC v1 借鉴映射（2026-08-12）

承接 `design/olc.md` 的 v1 待决项：

1. **渐进丰度过滤 → unitig 覆盖度驱动的布局前过滤**：
   `ProgressiveAbundanceFilter::removeAbundanceNoQueue`
   （`ProgressiveAbundanceFilter.hpp:2181`）：`t=1.1` 起步、`~10%` 步长、
   每轮删 `abundance < t` 的 unitig 并 recompact 邻接——不是单阈值一刀切。
   pgr `asm unitig` 头部已带 `cov=`，v1 可在 `asm olc` 布局前按 unitig
   丰度多轮剔除（或给 `asm unitig` 加渐进 `--min-coverage` 模式）。
2. **RepeatRemover 的桥接 reads 证据 → OLC repeat breaking 的实现路径**：
   `RepeatRemover.hpp:283` 把 reads 映射回 contig（`ReadVsContigMapper`，
   minimizer 索引）→ 按比对边界分片算覆盖度与 `_nbBridgingReads`
   （`:1195`）→ 无桥接的片段边界断开（判定在 `RepeatRemover.hpp:1257`，
   `nbContigsFinal>1` 即 split）。这正是 `canu.md` §8.3 预言的"pgr
   `asm map` + `sam to-rg` + `rg coverage` 回放"的成熟实现——pgr 设施
   齐全，v1 可直接照搬语义
   （桥接 reads = 覆盖度证据，对应 Celera 的 6/15 阈值）。
3. **cutoff 快照倒序输出 → 宏基因组 contig 分级输出**：多个 cutoff 的图
   快照、从高丰度往低丰度补——适合"先出高置信 contig、再补低丰度"策略。
4. **unitig 丰度中位数语义**：pgr `asm unitig` 的 `cov=` 是平均丰度；
   宏基因组场景若要稳健丰度，可参考"组成 k-min-mer 丰度向量取中位数、
   merge 时合并向量再取中位"（`Graph.hpp` `computeMedianAbundance`）。
