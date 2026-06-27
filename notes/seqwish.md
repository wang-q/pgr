# seqwish 分析笔记

> 整理于 2026-06-28，源自对 `seqwish-master/src/` 全部 22 个源文件的通读。 目的：理解
> seqwish 如何从 PAF 诱导出 GFA 变异图，并与 pgr 的隐式图路线对照，提取可借鉴的工程细节。

## 0. 项目定位

`seqwish`（"sequence wish"）是 PGGB 流水线中的**图物化器**（variation graph inducer）：
输入一组序列 + 它们的 all-vs-all PAF 比对，输出 GFA v1.0 变异图。它在 PGGB 中的位置是
`wfmash`（比对）→ **`seqwish`**（诱导图）→ `smoothxg`（归一化）。

一句话概括其本质：**把 pairwise 比对蕴含的"同源等价类"通过传递闭包物化成图节点，再沿输入
序列的邻接关系派生出图边。** 与 pgr/impg 的"隐式图"路线（不物化、按需 BFS）相对，seqwish
是"显式物化"路线的代表。本文档既是对其算法的拆解，也是 pgr V4（GFA 物化阶段）的直接参考。

## 1. 整体流程（6 阶段）

[main.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/main.rs) 把整个流程串成 6 个阶段，
每阶段对应一个模块：

1. **序列索引** — [seqindex.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/seqindex.rs)
   `SeqIndex::build_index`：读 FASTA/FASTQ，拼接成单一字节流，建 FM-index 索引序列名，
   用 `SparseBitVec` 记录序列边界，mmap 拼接文件做 O(1) 随机访问。
2. **比对索引** — [alignments.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/alignments.rs)
   `unpack_paf_alignments`：多线程读 PAF，按 CIGAR 解析成逐碱基的 match 段，双向写入
   `aln_iitree`（query→target 与 target→query 各一份）。
3. **传递闭包** — [transclosure.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/transclosure.rs)
   `compute_transitive_closures`：核心算法，把对齐位置划分成等价类，输出图序列 `seq_v` +
   `node_iitree`（图位置→输入位置）+ `path_iitree`（输入位置→图位置）。
4. **节点压缩** — [compact.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/compact.rs)
   `compact_nodes`：在图序列上标记节点边界（分叉/汇合点），把无分叉的线性段压成单节点。
5. **边派生** — [links.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/links.rs)
   `derive_links`：沿每条输入序列走一遍，记录相邻节点对，去重后得到边集。
6. **GFA 输出** — [gfa.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/gfa.rs)
   `emit_gfa`：写 S（节点序列）、L（边）、P（路径）三段，并校验路径与图序列碱基一致。

[main.rs#L200-L420](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/main.rs#L200-L420) 给出了
完整的编排，包含进度日志、tempfile 管理、Rayon 线程池配置等工程细节。

## 2. 关键数据结构

### 2.1 PosT：offset + 方向的单 u64 编码

[pos.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/pos.rs) 定义了贯穿全项目的位置类型：

```rust
pub type PosT = u64;
// 低位（bit 0）：方向（0=正链，1=反链）
// 高 63 位：在拼接序列中的 offset
pub fn make_pos_t(offset: u64, is_rev: bool) -> PosT { (offset << 1) | (is_rev as u64) }
pub fn offset(pos: PosT) -> u64 { pos >> 1 }
pub fn is_rev(pos: PosT) -> bool { (pos & 1) != 0 }
```

**亮点**：方向和 offset 打包进一个 u64，`incr_pos` / `decr_pos` 用 `±2` 步进，反链时反向步进。
所有区间树都以 `PosT` 为 value，单棵树同时表达"位置 + 链方向"。pgr 的 `pgr paf` 模块若要支持
反链投影，可直接借鉴此编码。

### 2.2 SeqIndex：FM-index + SparseBitVec

[seqindex.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/seqindex.rs) 把所有输入序列拼接
成单一字节流（写临时文件后 mmap），用两个结构索引：

- `name_index: FMIndexWithLocate<u8>` — 对 `">name1 >name2 ..."` 文本建 FM-index，
  支持 O(m) 时间按名字查找序列 rank。空间 O(n log σ) bit。
- `seq_boundaries: SparseBitVec` — 只存 1-bit 的位置数组，`select1` 是 O(1) 数组访问，
  `rank1` 是 O(log m) 二分。注释明确说这比 RsVec 在稀疏数据上更快。

```rust
struct SparseBitVec {
    positions: Vec<usize>, // 排序的 1-bit 位置
    size: usize,
}
fn select1(&self, i: usize) -> usize { self.positions[i] } // O(1)
fn rank1(&self, i: usize) -> usize { self.positions.binary_search(&i).unwrap_or_else(|x| x) }
```

**对 pgr 的启示**：pgr 现有 `src/libs/io.rs` 用 HashMap 存序列名→offset，对 4 万大肠杆菌尚可，
但若扩到 HPRC 规模（数百单倍型、Gb 级），FM-index + SparseBitVec 是更省内存的方案。
"稀疏 bitvec 只存 1-bit 位置"这个朴素思路在小 m / 大 N 场景下完胜压缩位向量。

### 2.3 AdaptiveTree：磁盘/内存双后端区间树

[intervaltree.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/intervaltree.rs) 定义了
`IntervalTree` trait 与两个后端：

- `disk::DiskBackedTree` — 封装 `iitree-rs`，数据落盘 mmap，适合超 RAM 的数据集。
- `memory::InMemoryTree` — 自实现 IAITree（隐式堆式二叉树 + max 增广），`overlap_iaitree`
  用 64 槽显式栈做迭代遍历，避免递归。小数据集快得多。

`AdaptiveTree` 是二者的枚举，由 `-M` / `--in-memory` 开关选择。主流程对三棵树
（`aln_iitree`、`node_iitree`、`path_iitree`）都用同一接口，写阶段 `open_writer` /
`add`，查阶段 `index` / `overlap`。

**对 pgr 的启示**：pgr 现用 `coitrees`（Crimson 对 iitree 的 Rust 实现），与 seqwish 的
`InMemoryTree` 同源同思想。seqwish 的"磁盘后端兜底大数据集"是 pgr 处理 4 万大肠杆菌时
可考虑的兜底方案——内存吃紧时切到 disk-backed，牺牲速度换可跑性。

### 2.4 DisjointSets：无锁并查集（3 个实现）

seqwish 为传递闭包的 union-find 提供了三个等价接口的实现：

- [dset64.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/dset64.rs) —
  `DisjointSets`，用 `Vec<AtomicU128>`，parent 和 rank 打包进 128 位原子。便携、跨平台。
- [dset64_asm.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/dset64_asm.rs) —
  `DisjointSetsAsm`，手写 16 字节对齐分配 + `CMPXCHG16B`，注释明确说"匹配 C++ 行为"。
- [dset64_unsafe.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/dset64_unsafe.rs) —
  `DisjointSetsUnsafe`，去掉边界检查的极致版本。

三者都实现了 Anderson & Woll 1991 的 wait-free union-find：`find` 做路径压缩，
`unite` 做 union-by-rank，CAS 失败时重试。128 位原子把 (parent, rank) 当一个字段，
避免双 CAS 的 ABA 问题。

```rust
// dset64.rs 核心循环
pub fn unite(&self, mut id1: usize, mut id2: usize) -> usize {
    loop {
        id1 = self.find(id1); id2 = self.find(id2);
        if id1 == id2 { return id1; }
        // union-by-rank: 把小 rank 的挂到大 rank 的下面
        let old = ((r1 as u128) << 64) | (id1 as u128);
        let new = ((r1 as u128) << 64) | (id2 as u128);
        if self.data[id1].compare_exchange(old, new, SeqCst, SeqCst).is_err() { continue; }
        // rank 相等时给新根加 1
        ...
    }
}
```

**对 pgr 的启示**：pgr 现在的 BFS 传递闭包是**查询时**按需做，规模小，普通 `HashSet` 就够。
但若 pgr V4a 要做"粗全局 GFA"（4 万大肠杆菌全图物化），等价类规模会到 Gbp 级，
此时无锁并查集是必备。`portable_atomic::AtomicU128` 在 x86_64 上自动用 `CMPXCHG16B`，
Apple Silicon 上用 LDXR/STXR，无需手写汇编——优先选 `DisjointSets`（dset64.rs）这个版本。

## 3. 传递闭包：seqwish 的算法核心

[transclosure.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/transclosure.rs) 是全项目
最复杂、也最有价值的模块。它解决的问题是：**给定 `aln_iitree`（pairwise 对齐区间），
把所有"对齐过"的输入位置划进同一个等价类，每个等价类成为图序列中的一个碱基。**

朴素做法是 N² 级的"对每个位置查全部对齐"。seqwish 的工程优化分四步：

### 3.1 第 0 步：最大权生成树剪枝

`compute_spanning_tree` 先扫一遍 `aln_iitree`，统计每对 (seq_i, seq_j) 的对齐碱基数作为权重，
跑 Kruskal 算法得到一棵最大权生成树。后续 BFS 只沿生成树边走，把 N(N-1)/2 对序列对齐
压缩到 N-1 对。

```rust
// 关键日志：显示压缩比
eprintln!("[transclosure] Spanning tree: {} edges from {} total pairs ({}x reduction)",
    tree_edges, edges.len(), edges.len() / tree_edges);
```

**这是 seqwish 最聪明的优化。** 生成树覆盖所有序列且连通，BFS 沿树边能发现所有连通分量，
不必扫全部对齐。后续 phase 2 再补全非树边的等价类合并。

### 3.2 Phase 1：BFS 发现（仅标记，不合并不收集）

主循环按 `transclose_batch`（默认 1Mb）切 chunk，每个 chunk 内：

1. 用 `for_each_fresh_range` 把未访问的种子位置标进 `q_curr_bv`（AtomicBitVec）。
2. 启动 `num_threads * 2` 个 worker + 1 个 manager，从 `todo_out` 队列取任务，
   调 `explore_overlaps_discovery` 查 `aln_iitree`。
3. **关键过滤**：`explore_overlaps_discovery` 内只追 `spanning_adj.contains(source, target)`
   的对齐，非树边直接跳过。
4. 发现新位置时 `curr_bv.set` 用 `fetch_or`（编译成 `LOCK OR`），返回旧值判断是否新，
   新的才 push 回 `todo_in`。

Phase 1 **只做位置发现**，不做 union-find，不收集 overlap 列表。这是相对 C++ 原版的关键改进——
C++ 版 phase 1 同时收集 ovlp_q 用于 phase 2，内存压力大；Rust 版用 spanning tree 把 phase 1
瘦身为纯 BFS 标记，phase 2 用 per-sequence 查询独立完成。

### 3.3 Phase 1b：孤儿恢复

生成树 BFS 可能漏掉"只能通过非树边到达"的位置。`orphan_recovery` 循环：

1. `find_component_sequences` 找出当前 chunk 已标记位置涉及的序列。
2. 对每条序列的 [offset, offset+len) 区间查 `aln_iitree`，把对端的未标记位置标进 `q_curr_bv`。
3. 若本轮无新位置则收敛，否则继续。

日志显示典型情况 1-3 轮就收敛。这是个"补漏"机制，保证等价类完整性。

### 3.4 Phase 2：并查集合并

phase 1b 收敛后，`q_curr_bv` 标记了当前 chunk 的全部相关位置。phase 2 对这些位置做 union-find：

1. 并行收集 `q_curr_positions: Vec<u64>`（所有 1-bit 的位置）。
2. 建 `rank_table: Vec<u32>`（位置→rank 的直接查表，O(1) 查找，比 sdsl rank 快）。
3. 对每条 component sequence 查 `aln_iitree`，对每个对齐区间内同时被 `q_curr_bv` 标记的
   (j, t) 位置对，调 `dsets.unite(rank_table[j], rank_table[t])`。
4. 输出 `(dset_id, position)` 列表，按 dset_id 排序、压缩、按最小 position 重命名，
   再排序一次。最终每个 dset 成为图序列中的一个碱基。

phase 2 用 `DisjointSetsAsm`（无锁 CAS），`component_seqs.par_iter()` 并行，是性能热点。

### 3.5 图序列写入：write_graph_chunk

等价类排好序后，`write_graph_chunk` 顺序遍历 `(dset_id, position)`：

- 每遇新 dset_id，从 `seqidx.at(offset)` 取该位置碱基，push 进 `seq_v_out`。
- 对该 dset 的每个输入位置 `curr_q_pos`，调 `extend_range` 把 (图位置, 输入位置) 对
  写进 `range_buffer`，满一段后 flush 到 `node_iitree` 和 `path_iitree`。
- `repeat_max` / `min_repeat_dist` 参数控制重复区过滤：同一序列在图同一位置的拷贝数
  超阈值时，不再写进 range_buffer，避免高拷贝重复把图吹爆。

`write_graph_chunk` 在独立线程跑，主线程同时处理下一个 chunk，实现流水线。

## 4. 节点压缩与边派生

### 4.1 compact_nodes：标记节点边界

[compact.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/compact.rs) 的逻辑很直接：
**对每个输入碱基，查 `path_iitree` 看它映射到图序列的哪一段，该段的起止就是节点边界。**

```rust
// 对每条序列并行
(1..=num_seqs).into_par_iter().for_each(|i| {
    let mut j = j_start;
    while j < k {
        path_iitree.overlap(j, j+1, |_, start, end, pos| {
            // 每个输入碱基应映射到唯一图位置
            // 在该图位置段的起止标记 1-bit
            seq_id_abv.set(offset(pos_start_in_s));
            seq_id_abv.set(offset(pos_end_in_s));
        });
        j = ovlp_end_in_q;
    }
});
```

输出 `seq_id_bv: BitVec`，1-bit 表示节点边界。再用 `RankSelectBitVector::from_bitvec`
转成只存 1-bit 位置数组，供后续 select/rank。

注意代码里有个 `panic!("Overlap count mismatch")` —— 这是 [compact.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/compact.rs)
里少见的硬 panic，因为这种情况说明前序算法出错（图断裂），不该被用户输入触发。
pgr 的"零 panic"原则在此处需要换成 `bail!` 并附诊断信息。

### 4.2 derive_links：沿输入序列走边

[links.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/links.rs) 的思路：
**节点 N1 到 N2 有边 ⟺ 存在某条输入序列，它在 N1 的末端紧接 N2 的起端。**

```rust
// 对每个节点 id（并行）
(1..=n_nodes).into_par_iter().map(|id| {
    let node_start = seq_id_cbv.select(id);
    let node_end = seq_id_cbv.select(id + 1);
    // 查 node_iitree：图 [start,end) 映射到哪些输入位置段
    node_iitree.overlap(node_start, node_end, |_, _, _, pos_in_q| {
        let end_in_q = ...; // 该输入段在节点末端的位置
        // 查 path_iitree：end_in_q 的下一个碱基映射到哪个图位置
        path_iitree.overlap(end_in_q, end_in_q+1, |_, _, _, pos_in_s| {
            let next_id = seq_id_cbv.rank(offset(pos_in_s) + 1);
            local_links.push((curr_node, next_node));
        });
    });
}).collect();
```

最后 `sort + dedup` 得到去重边集。这里用 `RwLock::read` 而非 `Mutex`，允许并发读两棵区间树。

### 4.3 emit_gfa：路径校验

[gfa.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/gfa.rs) 写 GFA 时有一段重要校验：
对每条输入序列的每个 overlap 段，逐碱基比对 `seqidx.at(q)` 与 `seq_v_slice[p]`（反链时取
complement），不一致则 `bail!` 报"GRAPH BROKEN"。这保证 P 路径与 S 节点序列严格一致，
是 seqwish 正确性的最后兜底。

## 5. 与 pgr 隐式图路线的对照

| 维度 | seqwish（显式物化） | pgr（隐式图） |
|------|---------------------|---------------|
| 输入 | PAF + 序列 | MAF→PAF + 序列 |
| 传递闭包 | **一次性全图** DSU | **查询时** BFS，按需局部 |
| 等价类表达 | 图序列的一个碱基 | 不物化，对齐区间即隐式等价类 |
| 数据结构 | iitree + DSU + 位向量 | coitrees + HashSet/Vec |
| 输出 | GFA（S+L+P） | BED/PAF（V1），GFA（V4） |
| 适用场景 | 全图分析、归一化、可视化 | 单 locus 查询、区域 MSA |
| 规模上限 | 受图序列长度限制（Gbp 级） | 受对齐索引大小限制（可分片） |
| 重复处理 | `--repeat-max` / `--min-repeat-dist` | `--min-len` / `--merge-distance` |

**核心差异**：seqwish 的传递闭包是**全局、一次性**的——算出全部等价类再写图。
pgr 的传递闭包是**局部、按需**的——每次查询从一个区间出发 BFS，只算相关等价类。
这两种粒度对应不同的应用场景：全图统计 vs 单点查询。

**seqwish 的 spanning tree 优化对 pgr 有直接借鉴价值。** pgr 现在的隐式图查询对每个
起点都做 BFS，如果对齐网络稠密（如 4 万大肠杆菌 K=50 的稀疏对齐仍有 135k 条边），
单次 BFS 可能横跳很多序列。若 pgr 在加载 PAF 阶段预计算一棵最大权生成树，
查询时优先沿树边走，可显著减少 BFS 的边遍历数。这是 [[paf-route.md]] 可考虑的优化项。

## 6. 对 pgr 各版本的启示

### 6.1 V1（坐标输出 bed/paf）

- **PosT 编码**：pgr 的 `pgr paf query` 若要支持反链投影，可借鉴 `make_pos_t` 把方向
  打包进 u64，单棵区间树同时存正反链对齐。
- **SparseBitVec**：pgr 处理 4 万大肠杆菌时，序列边界用 `SparseBitVec`（只存 1-bit 位置）
  比位向量省内存且 select O(1)。

### 6.2 V4a（粗全局 GFA）

- **直接复用 seqwish 算法骨架**：spanning tree → BFS discovery → DSU union-find →
  compact → links → GFA，这套流程对 pgr V4a 完全适用。
- **`--min-var-len 100` 过滤**：对应 seqwish 的 `--repeat-max` 思路，但维度不同——
  seqwish 过滤的是重复拷贝数，pgr 要过滤的是 SV 长度。可在 `write_graph_chunk` 阶段
  加一层"长度 < 100bp 的变异不写进图"的过滤，类似 minigraph 的粗框架哲学（见 [[minigraph.md]]）。
- **磁盘后端兜底**：4 万大肠杆菌全图可能超 RAM，`AdaptiveTree` 的 disk-backed 模式
  是现成的兜底方案。

### 6.3 V4b（局部精细 GFA）

- **phase 1b orphan recovery 的思路**：pgr 的局部 GFA 从一个 region 出发，BFS 发现的
  等价类可能不完整（只覆盖部分序列）。seqwish 的 orphan recovery 循环（按序列查 iitree
  补漏）可直接用于 pgr 的局部 GFA 完整性保证。
- **流水线写入**：`write_graph_chunk` 在独立线程跑、主线程算下一个 chunk 的模式，
  对 pgr 处理大区域 MSA 时同样适用。

### 6.4 工程细节

- **零 panic 原则**：seqwish 的 [compact.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/compact.rs)
  有 `panic!("Overlap count mismatch")`，pgr 应改为 `bail!` + 诊断信息，符合 CLAUDE.md 的稳定性要求。
- **进度日志**：seqwish 每阶段都打 `%` 进度，pgr 处理大规模数据时应沿用此风格。
- **tempfile 管理**：[tempfile.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/tempfile.rs)
  提供 `set_dir` / `set_keep_temp`，pgr 的 `pgr paf` 若产生中间文件可借鉴。

## 7. 不打算借鉴的部分

- **SXS 格式**（[sxs.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/sxs.rs)）：
  seqwish 私有的对齐格式，pgr 用 PAF/MAF，不需要。
- **dset64_asm.rs / dset64_unsafe.rs**：手写汇编和 unsafe 优化，`portable_atomic` 的
  `AtomicU128` 已足够快，pgr 用便携版 `DisjointSets` 即可。
- **`--sparse-factor` 哈希稀疏化**（[alignments.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/alignments.rs)
  的 `keep_sparse`）：seqwish 用哈希函数随机丢弃对齐，pgr 走 Mash KNN 稀疏化（见
  [[ecoli-cohort.md]]），质量更高，不用哈希稀疏化。

## 8. 参考链接

- 源码：[seqwish-master/src/](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/)
- 关联文档：[[pangenome-tools.md]] §3.2（PGGB 流水线中 seqwish 的位置）、
  [[impg.md]] §1.1.2（隐式图 vs 物化图适用边界）、[[minigraph.md]]（粗框架过滤哲学）、
  [[graph-design.md]]（pgr V4a/V4b 路线）、[[paf-route.md]]（pgr 隐式图核心原则）。
