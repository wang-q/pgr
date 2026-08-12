# SKESA / skesa-rs：de Bruijn 图短读组装器（源码分析）

> 2026-08 整理，纯源码分析。SKESA（`SKESA-master/`，C++）是 NCBI 的微生物
> 基因组 de-novo 短读组装器，论文 [SKESA: strategic k-mer extension for
> scrupulous assemblies](https://doi.org/10.1186/s13059-018-1540-z)（Genome
> Biology 2018）；`skesa-rs/` 是基于 SKESA v2.4.0 快照（commit
> `27caba2`，2024-10-11）的**逐位忠实 Rust 移植**（henriksson-lab/rustification
> 项目），并追求与 C++ 输出字节级一致。两者共享同一算法族，是 pgr 做
> k-mer 计数 / de Bruijn 图遍历的**首选参考**。
> **与 OLC 的连接（2026-08-12）**：pgr 的 `asm olc`（多 k unitig 层 OLC，
> `design/olc.md`）已落地，SKESA 的 fork 过滤 / 可逆性 / 迭代多 k 语义
> 直接映射其 v1 待决项（见 §7.1）。

## 1. 概况

- **定位**：SKESA 面向 Illumina 短读（单/双端），用**保守启发式**在重复区
  断裂，换取序列质量；k 从 mate 长度一直增到 insert size，兼顾 N50。
- **确定性**：同输入（含 read 顺序）下输出 contig 的顺序/方向**确定**——
  依赖排序 + 稳定启发式，不依赖多线程调度（这点对 pgr 做字节级一致性很关键）。
- **语言/构建**：
  - C++：Boost + gcc，`make`（NGS 版）/ `make -f Makefile.nongs`（文件版）；
    variant（`boost::variant<LargeInt<1>..LargeInt<16>>`）做运行时多精度 k-mer。
  - Rust：Cargo，依赖仅 `clap`（optional）、`noodles`、`flate2`、`rayon`；
    默认不编译 CLI（库形式），`--features cli` 才引 clap。
- **两者关系**：skesa-rs 是**逐函数移植**（含复刻 bug 以求可复现），C++ 与
  Rust 文件基本一一对应；Rust 在 `--cores 4` 基准上与 C++ 时间/RSS **基本持平**
  （README 快照 wall 0.972x、RSS 1.000x，contigs 输出 SHA-256 一致）。
- **未移植**：SAUTE / saute-prot / gfa-connector 完整对等、SRA 输入（Rust 显式拒绝）。

## 2. C++ 仓库结构（SKESA-master/）

| 文件 | 作用 |
|---|---|
| `skesa.cpp` | 主入口，参数解析（boost::program_options）、流程编排 |
| `Integer.hpp` / `LargeInt.hpp` / `LargeInt1/2.hpp` | 大整数 / 变长 k-mer 编码 |
| `KmerInit.hpp` / `Model.hpp` | k-mer 初始化、2-bit 编码表 |
| `concurrenthash.hpp` | **并发 k-mer 计数**：blocked bloom filter + 分块哈希表 |
| `counter.hpp` | **排序计数** `CKmerCount`（variant of 有序 vector） |
| `kmercounter.cpp` | 独立 k-mer 计数工具（`kmercounter` 子命令） |
| `DBGraph.hpp` | de Bruijn 图：`CDBGraph`（排序计数版）、`CDBHashGraph` |
| `graphdigger.hpp` | **图遍历 + contig 组装**（保守扩展、fork 解析） |
| `assembler.hpp` | **迭代组装编排**（多轮 k、paired-end 连接、insert 估计） |
| `guidedassembler.hpp` / `guidedgraph.hpp` / `guidedpath_naa.hpp` | SAUTE 目标富集组装 |
| `saute.cpp` / `saute_prot.cpp` | SAUTE 入口（核酸/蛋白） |
| `gfa.hpp` / `gfa_connector.cpp` | GFA 图输出 / contig 连接成图 |
| `glb_align.cpp/hpp` / `nuc_prot_align.hpp` | 全局/核酸蛋白比对 |
| `readsgetter.hpp` | FASTA/FASTQ/gzip/SRA 读入 + 适配器裁剪 |

Rust `src/` 同名对应：`concurrent_hash.rs`、`sorted_counter.rs`/
`counter.rs`/`flat_counter.rs`、`db_graph.rs`、`graph_digger.rs`、
`assembler.rs`、`guided_*.rs`、`snp_discovery.rs`、`linked_contig.rs`、
`paired_reads.rs`、`clean_reads.rs`、`glb_align.rs`、`gfa.rs`。

## 3. 核心数据结构

### 3.1 大整数 k-mer（`LargeInt.hpp` / Rust `large_int.rs` + `kmer.rs`）

- 2-bit 编码，k-mer 长度 → 精度 `precision = (kmer_len+31)/32`（用几个 u64
  存储），最大 16×64=1024 bit = 512 nt。
- C++ 用 `boost::variant<LargeInt<1>..LargeInt<16>>` 做运行时多精度；
  Rust 用 `enum Kmer { K1(LargeInt<1>), ..., K16(LargeInt<16>) }` +
  `macro_rules! define_kmer_enum` 展开 16 个变体的全操作分派（`kmer.rs`）。
- 关键操作：`revcomp`、`shl/shr`、`oahash`（SKESA 自定义哈希，`KmerOaHasher`
  复刻为 Rust `Hasher`）、`resize`（换精度，左截断/补零 + 顶字掩码）。
- 变长 k-mer 用 `u64` 字数组；`<=32nt` 走 `Flat((u64,u64))`，多字走内联数组
  保持缓存局部性。

### 3.2 并发 k-mer 计数（`concurrenthash.hpp` / Rust `concurrent_hash.rs`）

两条计数路线（`--hash_count` 切换）：

**A. blocked counting bloom filter（`CConcurrentBlockedBloomFilter<128>`）**
- 每 `SBloomBlock` 128 字节、`alignas(64)` 缓存行对齐（`concurrenthash.hpp:145`）；
  计数元素按块内位偏移打包（每计数器 2/4/8 bit）。
- 由**两个**哈希值 `(hashp, hashm)` 生成 `hash_num` 个哈希位：`hashp += hashm`
  迭代（`concurrenthash.hpp:77-85`），块内取 `hashp & (elements_in_block-1)`。
- 每块一个 `SAtomic<uint8_t>` 自旋锁；`Insert` 返回
  `eNewKmer / eAboveThresholdKmer / eExistingKmer`，用于**只把达到 min_count
  的 k-mer 灌入真实哈希表**——bloom 过滤是内存控制的关键。
- 计数封顶 `m_max_element = (1<<counter_bit_size)-1`（饱和计数）。

**B. 分块并发哈希表（`SHashBlock<Key,V,BucketBlock=32>`）**
- 每桶一个**小数组（≤32 槽）+ 溢出 forward list**（`concurrenthash.hpp:395-445`）：
  先试哈希指定位置，再线性扫小数组，最后溢出链表——`BucketBlock=32` 折中缓存
  命中与溢出量。
- 每槽状态原子 `eAssigned / eKeyExists`，`Lock/Wait` 实现无锁读 + 自旋写；
  `CDeque` 分块并行初始化大表。
- k-mer 的 (key, count) 打包存储；`count` 在计数期低 32 位 total、高 32 位
  plus-strand（见 §3.3），`(plusf<<48)+(branches<<32)+total`（`concurrenthash.hpp:1380`）。

### 3.3 排序计数（`counter.hpp` / Rust `counter.rs`/`sorted_counter.rs`/`flat_counter.rs`）

- `CKmerCount` 是 `vector<pair<LargeInt<N>, size_t>>` 的 variant，**只存 canonical**
  （kmer 与其 revcomp 中较小的）；排序后二分查找（`lower_bound`）。
- **count 打包**（`counter.hpp:42-44` + `DBGraph.hpp:180`）：
  ```
  低 32 bit: total count（self+revcomp）
  高 32 bit: 计数期=plus-strand count；进 CDBGraph 后重排为:
             [0:31]=total | [32:39]=8bit 分支信息 | [40:47]=未用 | [48:63]=16bit plus-fraction
  ```
- `--memory`（GB）决定**多轮外部归并**：内存预算 → 每轮可装多少元素 →
  分块排序落盘再归并（`counter.hpp` 的 `SortAndExtractUniq`/`MergeTwoSorted`）。
- Rust `KmerCount` 用 `enum Storage` 区分 `Flat/Words2..8/General`，1..8 字内联、
  其余 boxed 兜底（`counter.rs`）；排序用 `rayon::par_sort_unstable_by_key`
  （>10000 并行、否则串行，保证小规模确定性）。`flat_counter.rs`/`sorted_counter.rs`
  是不同精度/计数策略的变体。

### 3.4 de Bruijn 图（`DBGraph.hpp` / Rust `db_graph.rs`）

- **节点编码**：`Node` 包一个 `size_t m_node`；偶数=正链、奇数=负链、0=无效；
  `Index() = m_node/2 - 1` 映射回数组下标（`DBGraph.hpp:102-124`）。
  图里只存 canonical k-mer，`GetNode` 对 `kmer<revcomp` 返回正链节点、否则
  revcomp 的负链节点（`DBGraph.hpp:155-164`）。
- **分支信息**：`GetNodeSuccessors` 读 count 高 32 位的 8bit 分支掩码，负链取
  高 4 位、正链取低 4 位（`DBGraph.hpp:249-267`）；`shifted=(kmer<<2)&max_kmer`
  + nt 查后继。**用打包位快速跳过无后继的碱基**，避免逐个探测。
- **strand 信息**：`PlusFraction = (count>>48)/65535`（`DBGraph.hpp:185-190`），
  供图遍历区分正负链计数。
- visited 用原子 uint8（1=永久占用、2=临时、3=多 contig），多线程安全标色
  （`DBGraph.hpp:204-222`）。
- **哈希版 visited 复用了 count 的高位**（`CDBHashGraph`，`DBGraph.hpp:432,454`）：
  排序版 `CDBGraph` 用独立 `m_visited` 数组（`DBGraph.hpp:301`），而
  `CDBHashGraph` 把 count 的 `[40:47]`（对应 `DBGraph.hpp:180` 的"未用" 8 bit）
  改作 visited 控制位（`eVisited/eTemp/eMulti = 1<<40 / 1<<41 / 1<<42`，
  `SetColor` 用 `mask<<40` 打色、`GetColor` 取 `(count>>40)&0xFF`）——省下
  一个 per-node 数组，是"状态打包进计数"的另一处实例。
- Rust `SortedDbGraph` 完全对应（`db_graph.rs`），另有 `HashNode` 对应哈希计数版 `CDBHashGraph`；两个具体图共用 `DBGraph` trait。

## 4. 图遍历与 contig 组装（`graphdigger.hpp` / Rust `graph_digger.rs`）

核心是**保守扩展 + 在重复区断裂**，用"只沿唯一、可信路径延伸"换质量：

- **fork 类型**（`graphdigger.hpp:93`）：`eNoFork/eLeftFork/eRightFork/
  eLeftBranch/eRightBranch/eSecondaryKmer`——记录左右分支与次生 k-mer。
- **后继过滤**（`FilterNeighbors` / `FilterLowAbundanceNeighbors`，
  `graphdigger.hpp:1769-1887`），按序：
  1. **低丰度 fork 剔除**：`abundance(后继) <= fraction × Σabundance` 的删除
     （`fraction` 即 `--fraction` 默认 0.1，噪音/信号比上限）；`LowCount()==1`
     且首后继丰度>5 时，把丰度==1 的尾巴删掉。
  2. **strand 特异的 Illumina 噪音**（`GGT→GG[ACG]` 现象）：对以 `GGT` 结尾的
     后继，用 `abundance×(1-PlusFraction)` 与 `fraction×am` 比较剔噪——正负链
     两处分别处理（`graphdigger.hpp:1793-1815, 1837-1859`）。
  3. **不可扩展 fork**：首后继丰度>5 时，剔除 `ExtendableSuccessor` 为假的。
  4. **strand 平衡问题**：存在 `min(plusf,minusf)>0.25` 的双链好节点时，剔掉
     `min(plusf,minusf) < 0.1×fraction×max(...)` 的偏链后继
     （`graphdigger.hpp:1861-1884`）。
- **可逆性检查**（`GetReversibleNodeSuccessors`）：扩展前验证每个后继再回退
  （对后继的 revcomp 求后继）能回到原节点，否则该 fork 不可逆、断裂
  （`graphdigger.hpp:1739-1762`）。
- **conservative 扩展**：沿唯一可逆路径延伸；遇到 fork/低置信即停，宁可断。
- **`jump`/`max_snp_len`**（`--max_snp_len` 默认 150）：扩展时允许跨过一个
  SNP 的"跳"以桥接多态区；`--allow_snps` 开 `check_repeats` 额外做 SNP 发现
  （Rust `snp_discovery.rs`）。

## 5. 迭代组装编排（`assembler.hpp` / Rust `assembler.rs`）

`CDBGAssembler` 多轮迭代，k 从小到大逐步解重复：

1. **建首图 @ min_kmer**（默认 21）；算 read 平均长、genome size 估计
   （从 k-mer histogram 的 `CalculateGenomeSize`）。
2. **自动抬阈值**（`assembler.hpp:963-981`，Rust `assembler.rs:178-199`）：
   若 coverage 过高，`new_min_count = coverage/50`、`new_max_kmer_count =
   coverage/10`，并 `remove_low_count` 剪枝。
3. **GenerateNewSeeds → ImproveContigs**：`graph_digger` 保守组装出 seed contig
   （jump=0 的保守版）；有 `--seeds` 则从种子扩展；`mark_previous_contigs`
   标已用 k-mer，`assemble_contigs_with_visited` 找新种子（避开已组装区）。
4. **max_kmer 估计**：`max_kmer = read_len+1 - (max_kmer_count/avg_count)×(read_len-min_kmer+1)`，
   clamp 到奇数（`assembler.rs:315-328`）。
5. **paired-end 连接**：
   - `estimate_insert_size`（抽样 10000 对，用首轮图估 insert N50，clamp 到
     `MAX_KMER`）；`paired_insert_limit = 3×N50`。
   - 若 `N50 > 1.5×max_kmer` 才启用**长 insert 双端迭代**（`use_long_paired_iterations`），
     否则直接 `connect_pairs` 在首图连 mate。
   - Rust `paired_reads.rs`：`connect_pairs` / `estimate_insert_size_full`。
6. **clean reads**：`clean_reads` 把完全落在已组装 contig 内的 read 剔除
   （`cleanup_min_contig_len = max(max_kmer, paired_insert_n50)`，Rust
   `assembler.rs:407`），防止陈旧 k-mer 污染下一轮 histogram。
7. **后续轮**：`max_kmer` 往上的每轮，用上一轮 contig 做"引导"，把未解重复区
   用更长 k 重连；`linked_contig.rs` 的 `ConnectFragments` 走连接链。
8. 输出 contigs（`--min-contig` 默认 200 过滤），可选 GFA。

> 关键工程点：**clean 的阈值语义**——min_contig_len 用的是 `max_kmer` 与
> `paired_insert_n50`（连接 mate 的 N50），**不是** insert 的 3×N50 上限；
> 用错会把 (N50, 3×N50) 区间已组装的 contig 排除在 kmer→contig 映射外，导致
> 陈旧 k-mer 残留（Rust `assembler.rs:400-407` 注释专门说明）。

## 6. SAUTE / 引导组装（`guidedassembler.hpp`、`saute.cpp`）

- SAUTE 用**目标序列（参考）引导**：`guidedgraph`/`guidedpath_naa` 对目标区域
  做目标富集 de Bruijn 组装，输出 GFA + 两条 FASTA。
- Rust `guided_assembly.rs`/`guided_graph.rs`/`guided_path.rs` 仅为**简化版辅助**
  （README 明确 "full SAUTE parity is not yet implemented"）；`spider_graph.rs`
  对应 gfa-connector 的连接路径辅助，同样未完全对等。

## 7. 与 pgr 的关联 / 借鉴点

pgr 已有 `libs/kmer`（KmerTable：canonical 2-bit u128、精确计数、radix sort、
rayon 并行），是**精确计数路线**；SKESA/skesa-rs 提供两条互补路线 + 全套
de Bruijn 图遍历启发式：

1. **count 打包布局**（`DBGraph.hpp:180`）：total(32bit)+branch(8bit)+plus-fraction(16bit)
   一个 u64 同时承载计数/分支/链向——pgr 若扩展 k-mer 表，可参考这种打包省内存。
2. **canonical 存储 + Node 偶数/奇数编码**（`DBGraph.hpp:102-164`）：图只存
   canonical，用奇偶位表达链向、`Index()=m_node/2-1` 映射数组——pgr `KmerTable`
   已是 canonical key，若做 de Bruijn 图可直接套用该节点编码。
3. **分支位快速找后继**：8bit 分支掩码 + `(kmer<<2)&max` + 打包位跳过，避免
   逐碱基探测——对 pgr 未来 `asm` 类功能的图遍历是高性能范式。
4. **保守扩展 + fork 过滤启发式**（`graphdigger.hpp:1769-1887`）：低丰度 fork
   剔除、strand 特异 GGT 噪音剔除、strand 平衡检查、可逆性检查——四层过滤是
   "在重复区断裂"的实现核心，pgr 若实现 de Bruijn 组装应移植这层语义。
5. **迭代多 k + paired-end 连接 + read 清理**：从 min k 到 max k 渐进解重复，
   clean 阈值用 `max(max_kmer, paired_insert_n50)` 的细节值得照搬。
6. **确定性**：排序 + 稳定启发式保证输出确定（多线程不改变结果）——与 pgr
   "字节级一致"的硬约束一致，是排序计数优于哈希计数的理由之一。
7. **Rust 移植经验**（skesa-rs 对 pgr 的直接价值）：
   - `enum Kmer` + 宏展开 16 精度变体，替代 C++ boost::variant；
   - `Storage` enum 按精度选内联数组 vs boxed，兼顾缓存与正确性；
   - 排序阈值（>10000 并行 else 串行）保小规模确定性；
   - 用 `noodles`+`flate2` 替代 C++ NGS 库；把 CLI 做成 optional feature（库优先）。

### 7.1 OLC v1 借鉴映射（2026-08-12）

承接 `design/olc.md` 的 v1 待决项，SKESA 提供三块直接素材：

1. **覆盖度/丰度驱动的 fork 过滤 → OLC repeat breaking 参数**：
   `FilterLowAbundanceNeighbors`（`graphdigger.hpp:1770`）的
   `abundance <= fraction × Σabundance`（`--fraction` 默认 0.1）与
   "`LowCount()==1` 且首后继丰度>5 时删丰度==1 尾巴"；`FilterNeighbors`
   的不可扩展 fork 剔除（`:1827`）与 strand 平衡检查（`:1863`，
   `min(plusf,minusf) < 0.1×fraction×max` 剔偏链）。这是"在重复区断裂"的
   成熟多层阈值语义——pgr layout 的 v0 repeat 检测只有 top2 近等边近似
   （`canu.md` §8.5 记录了 6× 低覆盖漏检案例），v1 应移植这层丰度阈值
   （pgr `asm unitig` 头部已带 `cov=`，可直接取用）。
2. **可逆性检查（`GetReversibleNodeSuccessors`，`graphdigger.hpp:1740`）→
   layout 延伸的"回得来"保证**：SKESA 扩展前验证后继的 revcomp 能回到原
   节点、否则断裂；与 pgr layout 的互惠 best edge（`canu.md` §8.5 连接端
   语义）是同一思想的两个实现——SKESA 在 k-mer 图、pgr 在 unitig 重叠图。
3. **迭代多 k + read 清理 → `asm olc` 的多 k 反馈**：pgr `asm olc` 目前
   各 k 独立出 unitig 再合并（无反馈）；SKESA 的"每轮用上一轮 contig 引导
   + `clean_reads`（阈值 `max(max_kmer, paired_insert_n50)`，
   `assembler.rs:407`）"与 metaMDBG 的 unitig 反馈（`metaMDBG.md` §4.1）
   同族，是 v2 候选。
4. **count 打包 + Node 奇偶编码（`DBGraph.hpp:180/102`）**：若 pgr 从
   bcalm 式哈希表切换到排序 DBG（kmer 表已走 radix 排序路线），
   total(32)+branch(8)+plus-fraction(16) 打包与 `Index()=m_node/2-1`
   是现成模式。

## 8. 局限

- C++ 依赖 Boost（variant）+ NGS 库（SRA），构建链较重；Rust 版已剥离。
- SAUTE / gfa-connector / SRA 在 Rust 版未对等移植。
- 哈希计数（bloom+分块表）有哈希碰撞/近似性，确定性弱于排序计数；pgr 的精确
  路线与哈希路线各有所长（内存 vs 精确）。
- SKESA 定位微生物短读组装；长读（HiFi/ONT）场景不在其路线（参照 metaMDBG）。
- skesa-rs 是 LLM 中介的"忠实翻译"，README 明示**不可完全信任**、复刻 bug、
  需自行验证——参考其工程手法而非当作权威实现。

---

*参考来源: `SKESA-master/`（skesa.cpp、concurrenthash.hpp、counter.hpp、
DBGraph.hpp、graphdigger.hpp、assembler.hpp、guidedassembler.hpp、kmercounter.cpp、
gfa.hpp）+ `skesa-rs-main/`（src/kmer.rs、counter.rs、db_graph.rs、assembler.rs、
graph_digger.rs、cli.rs、Cargo.toml、README.md）*
