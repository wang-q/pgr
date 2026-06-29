# seqwish 分析笔记

> 整理于 2026-06-28，源自对 `seqwish-master/src/` 全部 22 个源文件的通读。 目的：理解 seqwish
> 如何从 PAF 诱导出 GFA 变异图，并与 pgr 的隐式图路线对照，提取可借鉴的工程细节。

## 0. 项目定位

`seqwish`（"sequence wish"）是 PGGB 流水线中的**图物化器**（variation graph inducer）：
输入一组序列 + 它们的 all-vs-all PAF 比对，输出 GFA v1.0 变异图。它在 PGGB 中的位置是
`wfmash`（比对）→ **`seqwish`**（诱导图）→ `smoothxg`（归一化，详见 [[smoothxg.md]]）。

一句话概括其本质：
**把 pairwise 比对蕴含的"同源等价类"通过传递闭包物化成图节点，再沿输入 序列的邻接关系派生出图边。**
与 pgr/impg 的"隐式图"路线（不物化、按需 BFS）相对，seqwish 是"显式物化"路线的代表。
本文档既是对其算法的拆解，也是 pgr graph / to-gfa（GFA 物化阶段）的直接参考。

## 1. 整体流程（6 阶段）

[main.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/main.rs) 把整个流程串成 6 个阶段，
每阶段对应一个模块：

1. **序列索引** — [seqindex.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/seqindex.rs)
   `SeqIndex::build_index`：读 FASTA/FASTQ，拼接成单一字节流，建 FM-index 索引序列名，用
   `SparseBitVec` 记录序列边界，mmap 拼接文件做 O(1) 随机访问。
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

[main.rs#L200-L420](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/main.rs#L200-L420)
给出了 完整的编排，包含进度日志、tempfile 管理、Rayon 线程池配置等工程细节。

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

- `name_index: FMIndexWithLocate<u8>` — 对 `">name1 >name2 ..."` 文本建 FM-index， 支持 O(m)
  时间按名字查找序列 rank。空间 O(n log σ) bit。
- `seq_boundaries: SparseBitVec` — 只存 1-bit 的位置数组，`select1` 是 O(1) 数组访问， `rank1` 是
  O(log m) 二分。注释明确说这比 RsVec 在稀疏数据上更快。

```rust
struct SparseBitVec {
    positions: Vec<usize>, // 排序的 1-bit 位置
    size: usize,
}
fn select1(&self, i: usize) -> usize { self.positions[i] } // O(1)
fn rank1(&self, i: usize) -> usize { self.positions.binary_search(&i).unwrap_or_else(|x| x) }
```

**对 pgr 的启示**：pgr 现有 `src/libs/io.rs` 用 HashMap 存序列名→offset，对 4 万大肠杆菌尚可，
但若扩到 HPRC 规模（数百单倍型、Gb 级），FM-index + SparseBitVec 是更省内存的方案。"稀疏 bitvec
只存 1-bit 位置"这个朴素思路在小 m / 大 N 场景下完胜压缩位向量。

### 2.3 AdaptiveTree：磁盘/内存双后端区间树

[intervaltree.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/intervaltree.rs) 定义了
`IntervalTree` trait 与两个后端：

- `disk::DiskBackedTree` — 封装 `iitree-rs`，数据落盘 mmap，适合超 RAM 的数据集。
- `memory::InMemoryTree` — 自实现 IAITree（隐式堆式二叉树 + max 增广），`overlap_iaitree` 用 64
  槽显式栈做迭代遍历，避免递归。小数据集快得多。

`AdaptiveTree` 是二者的枚举，由 `-M` / `--in-memory` 开关选择。主流程对三棵树 （`aln_iitree`、
`node_iitree`、`path_iitree`）都用同一接口，写阶段 `open_writer` /`add`，查阶段 `index` /
`overlap`。

**InMemoryTree 的 IAITree 算法**：`finalize()` 先按 `(start, end)` 排序，再自底向上算 每个节点的
`max`（子树最大 end），是 iitree-rs `index_core()` 的直接移植。查询时 `overlap_iaitree` 用
64 槽显式栈做迭代（不递归），核心剪枝：节点的 `max <= query_start`时整棵子树跳过。小 subtree
（`k <= 3`，即 ≤8 个节点）退化为线性扫描，避免栈操作开销。

```rust
// overlap_iaitree 的核心分支
if z.k <= 3 { /* 小 subtree: 线性扫描 */ }
else if z.w == 0 { /* 左子未处理: push 自身+左子 */ }
else if z.x < n && start < query_end { /* 处理当前+push 右子 */ }
```

**对 pgr 的启示**：pgr 现用 `coitrees`（Crimson 对 iitree 的 Rust 实现），与 seqwish 的
`InMemoryTree` 同源同思想。seqwish 的"磁盘后端兜底大数据集"是 pgr 处理 4 万大肠杆菌时
可考虑的兜底方案——内存吃紧时切到 disk-backed，牺牲速度换可跑性。

### 2.4 DisjointSets：无锁并查集（3 个实现）

seqwish 为传递闭包的 union-find 提供了三个**接口等价、存储表示不同**的实现。三者都基于 Anderson &
Woll 1991 的 wait-free 算法：`find` 做路径压缩（CAS 失败容忍，不重试），`unite` 做 union-by-rank
（CAS 失败重试）。核心技巧是把 (parent, rank) 打包进单个 128 位 原子，避免双 CAS 的 ABA 问题——parent
在低 64 位，rank 在高 64 位。

**存储表示差异**（决定性能与可移植性）：

- [dset64.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/dset64.rs) — `DisjointSets`，
  `data: Vec<AtomicU128>`。便携版，依赖 `portable_atomic::AtomicU128`的 `is_lock_free()` 断言。
  x86_64 上 `AtomicU128` 由 `portable_atomic` 在内部用 `CMPXCHG16B` 实现（需要 `target-cpu=native`
  或 `+cmpxchg16b` feature），ARM64 用 LDXR/STXR。索引访问有 Rust 标准库的边界检查。
- [dset64_asm.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/dset64_asm.rs) —
  `DisjointSetsAsm`，`data: *mut AlignedU128`，其中 `AlignedU128` 是 `#[repr(C, align(16))]`
  的包装。手写 16 字节对齐分配（`Layout::array::<AlignedU128>`），注释明确说"匹配 C++ 行为"——C++
  原版用 `alignas(16)` + `__sync_bool_compare_and_swap`，对齐的 128-bit load 编译成原子 SSE 指令。
  比 `Vec` 版快，但需要 `unsafe impl Send/Sync`。
- [dset64_unsafe.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/dset64_unsafe.rs) —
  `DisjointSetsUnsafe`，`data: *mut AtomicU128`。和 `DisjointSets` 用同一个 `AtomicU128` 类型，
  但用裸指针 + `unsafe { self.data.add(id) }` 跳过边界检查，所有访问走 `*_unchecked` 方法。`Drop`
  里手动 `drop_in_place` 每个元素再 `dealloc`。

**实际选用**：
[transclosure.rs#L12](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/transclosure.rs#L12)
`use crate::dset64_asm::DisjointSetsAsm;` —— phase 2 性能热点用的是**asm 版**，便携版和 unsafe
版主要作为对照/备份保留。三者 `unite` 的核心循环结构一致：

```rust
// unite 的关键路径（三者仅在 CAS 调用方式上不同）
loop {
    id1 = self.find(id1); id2 = self.find(id2);
    if id1 == id2 { return id1; }
    // union-by-rank: r1 < r2 或 (r1 == r2 且 id1 < id2) 时交换，保证小 rank 挂大 rank
    let old = ((r1 as u128) << 64) | (id1 as u128);  // parent=id1, rank=r1
    let new = ((r1 as u128) << 64) | (id2 as u128);  // parent=id2, rank=r1
    // CAS 失败说明并发改写，重试整个 find+unite
    if self.data[id1].compare_exchange(old, new, SeqCst, SeqCst).is_err() { continue; }
    if r1 == r2 {
        // rank 相等时给新根 id2 的 rank +1（CAS 失败容忍，因为 rank 只增不减）
        let _ = self.data[id2].compare_exchange_weak(
            ((r2 as u128) << 64) | (id2 as u128),
            (((r2 + 1) as u128) << 64) | (id2 as u128),
            SeqCst, SeqCst);
    }
    return id2;
}
```

**内存序**：全部 `SeqCst`，注释说"to match C++ `__sync_bool_compare_and_swap`"。 这是 C++
原版遗留的保守选择——`__sync_*` 内建函数全是 full barrier，Rust 移植时 直接照搬语义。理论上 `unite`
的 CAS 用 `AcqRel`/`Acquire`+`Release` 就够，`find` 的路径压缩 CAS 用 `Relaxed` 即可，但 seqwish
没做这个优化。

**对 pgr 的启示**：pgr 现在的 BFS 传递闭包是**查询时**按需做，规模小，普通 `HashSet` 就够。 但若
pgr graph 要做"粗全局 GFA"（4 万大肠杆菌全图物化），等价类规模会到 Gbp 级，此时无锁并查集是必备。
选型建议：

- **首选 `DisjointSets`（dset64.rs）**：`portable_atomic::AtomicU128` 在 x86_64 上自动用
  `CMPXCHG16B`，Apple Silicon 上用 LDXR/STXR，无需手写汇编也无需 `unsafe impl Send/Sync`。只要
  `Cargo.toml` 加 `portable_atomic` 依赖 + build.rs 设 `target-cpu=native` 即可。
- **不选 asm/unsafe 版**：asm 版的 16 字节对齐分配在 `portable_atomic` 内部已做（它会对齐到 16
  字节以满足 `CMPXCHG16B` 的对齐要求），手写一遍没有收益；unsafe 版的"跳边界检查"在现代 Rust
  编译器下收益有限（LLVM 能根据上下文消除已知安全的检查），却牺牲了 UB 安全性。

## 3. 传递闭包：seqwish 的算法核心

[transclosure.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/transclosure.rs)
是全项目 最复杂、也最有价值的模块。它解决的问题是：
**给定 `aln_iitree`（pairwise 对齐区间）， 把所有"对齐过"的输入位置划进同一个等价类，每个等价类成为图序列中的一个碱基。**

朴素做法是 N² 级的"对每个位置查全部对齐"。seqwish 的工程优化分四步：

### 3.1 第 0 步：最大权生成树剪枝

`compute_spanning_tree` 先扫一遍 `aln_iitree`，统计每对 (seq_i, seq_j) 的对齐碱基数作为权重， 跑
Kruskal 算法得到一棵最大权生成树。后续 BFS 只沿生成树边走，把 N(N-1)/2 对序列对齐 压缩到 N-1 对。

```rust
// 关键日志：显示压缩比
eprintln!("[transclosure] Spanning tree: {} edges from {} total pairs ({}x reduction)",
    tree_edges, edges.len(), edges.len() / tree_edges);
```

**这是 seqwish 最聪明的优化。**生成树覆盖所有序列且连通，BFS 沿树边能发现所有连通分量，
不必扫全部对齐。后续 phase 2 再补全非树边的等价类合并。

### 3.2 Phase 1：BFS 发现（仅标记，不合并不收集）

主循环按 `transclose_batch`（默认 1Mb）切 chunk，每个 chunk 内：

1. 用 `for_each_fresh_range` 把未访问的种子位置标进 `q_curr_bv`（AtomicBitVec）。
2. 启动 `num_threads * 2` 个 worker + 1 个 manager，从 `todo_out` 队列取任务， 调
   `explore_overlaps_discovery` 查 `aln_iitree`。
3. **关键过滤**：`explore_overlaps_discovery` 内只追 `spanning_adj.contains(source, target)`
   的对齐，非树边直接跳过。
4. 发现新位置时 `curr_bv.set` 用 `fetch_or`（编译成 `LOCK OR`），返回旧值判断是否新， 新的才 push
   回 `todo_in`。

Phase 1 **只做位置发现**，不做 union-find，不收集 overlap 列表。这是相对 C++ 原版的关键改进—— C++
版 phase 1 同时收集 ovlp_q 用于 phase 2，内存压力大；Rust 版用 spanning tree 把 phase 1 瘦身为纯
BFS 标记，phase 2 用 per-sequence 查询独立完成。

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
3. 对每条 component sequence 查 `aln_iitree`，对每个对齐区间内同时被 `q_curr_bv` 标记的 (j, t)
   位置对，调 `dsets.unite(rank_table[j], rank_table[t])`。
4. 输出 `(dset_id, position)` 列表，按 dset_id 排序、压缩、按最小 position 重命名， 再排序一次。
   最终每个 dset 成为图序列中的一个碱基。

phase 2 用 `DisjointSetsAsm`（无锁 CAS），`component_seqs.par_iter()` 并行，是性能热点。

### 3.5 图序列写入：write_graph_chunk

等价类排好序后，`write_graph_chunk` 顺序遍历 `(dset_id, position)`：

- 每遇新 dset_id，从 `seqidx.at(offset)` 取该位置碱基，push 进 `seq_v_out`。
- 对该 dset 的每个输入位置 `curr_q_pos`，调 `extend_range` 把 (图位置, 输入位置) 对 写进
  `range_buffer`，满一段后 flush 到 `node_iitree` 和 `path_iitree`。
- `repeat_max` / `min_repeat_dist` 参数控制重复区过滤：同一序列在图同一位置的拷贝数 超阈值时，
  不再写进 range_buffer，避免高拷贝重复把图吹爆。

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

输出 `seq_id_bv: BitVec`，1-bit 表示节点边界。再用 `RankSelectBitVector::from_bitvec` 转成只存
1-bit 位置数组，供后续 select/rank。

#### 4.1.1 两套 AtomicBitVec 实现

seqwish 在两处独立实现了 `AtomicBitVec`，**接口不同、语义不同**，没有复用：

- [transclosure.rs#L72-L110](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/transclosure.rs#L72-L110)
  的版本：`fn set(&self, index, value, ordering) -> bool`，用 `fetch_or` 返回**旧值**，
  调用方据此判断"这个位置是不是我刚标的"（新位置才 push 进 BFS 队列）。这是 phase 1 BFS 的去重关键——
  `fetch_or` 在 x86 上编译成 `LOCK OR`，原子地"标记并返回是否新"。
- [compact.rs#L19-L77](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/compact.rs#L19-L77)
  的版本：`fn set(&self, index) -> ()`（fire-and-forget，不返回旧值），额外提供 `iter_ones()`
  迭代器。compact 阶段不需要判断新旧（每个边界只需标一次，重复标记无害），省掉返回值开销；
  `iter_ones` 在最后把原子字 load 出来扫一遍，用于输出 `seq_id_bv`。

**对 pgr 的启示**：这两套实现本质是同一思路（`Vec<AtomicU64>` + `fetch_or`），分裂成两份 是历史包袱。
pgr 若引入原子位向量，应统一成一个带 `set -> bool`（返回旧值）的 API，`iter_ones` 作为 trait
method 提供，避免重复。`fetch_or` 的"返回旧值判新"是并行 BFS 去重的标准技巧，比 `Mutex<HashSet>`
快一个数量级。

#### 4.1.2 硬 panic 的诊断价值

[compact.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/compact.rs) 里有一处少见的 硬
panic：

```rust
if overlap_count != 1 {
    eprintln!("[compact] error: found {overlap_count} overlaps for seq {seq_name} idx {i} at j={j} of {k}");
    // ... 打印所有 overlap 的详细信息 ...
    panic!("Overlap count mismatch");
}
```

每个输入碱基应映射到**恰好一个**图位置——0 个说明图断裂，> 1 个说明 `path_iitree` 写重了。
两种都是前序算法（transclosure/write_graph_chunk）出 bug，不该被用户输入触发，所以用 panic 而非
`bail!`。panic 前先 `eprintln!` 打印所有 overlap 的 `(start, end, pos)` 供调试。

pgr 的"零 panic"原则（CLAUDE.md）要求此处改成 `bail!` + 诊断信息，但 seqwish 的"panic 前
打印完整上下文"做法值得借鉴——pgr 的 `bail!` 也应包含 `(seq_name, j, k, overlap_count)` +所有
overlap 详情，否则用户报错时无法定位是哪个 chunk 的哪个序列出了问题。

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

### 4.3 emit_gfa：多级路径校验

[gfa.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/gfa.rs) 在写 P 行之前对每条
输入序列做**4 级校验**，任一级失败都 `return Err` 带详细诊断。这是 seqwish 正确性的最后兜底，保证
P 路径与 S 节点序列严格一致。

**Level 1 — overlap 唯一性**（每碱基）： 对输入序列的每个碱基 `j` 查 `path_iitree.overlap(j, j+1)`，
`overlap_count` 必须恰好为 1。0 个说明图断裂（transclosure 漏标），> 1 说明 `path_iitree` 写重了。
失败信息：`"found {overlap_count} overlaps for seq {seq_name} idx {i} at j={j} of {k}"`。与
[compact.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/compact.rs) §4.1.2 的检查 同源，
但 compact 用 panic，gfa 用 Err——因为 gfa 是最后输出阶段，错误可恢复（跳过该序列）。

**Level 2 — 逐碱基一致性**（每 overlap 段）： 对每个
overlap 段的 `(q, p)` 对，取 `input_char = seqidx.at(q)` 与
`graph_char = seq_v_slice[offset(p)]`（反链时取 `complement`），逐碱基比对。不一致即报：
`"[gfa] GRAPH BROKEN @ {seq_name} pos {q} -> graph pos {p_offset}: expected {input_char} got {graph_char}"`。
这直接定位是哪个输入碱基与图序列对不上，是调试 transclosure/compact bug 的关键线索。

**Level 3 — 总覆盖长度**（每序列）： `seen_bp`（所有 overlap 段长度之和）必须等于 `seq_len`。
失败信息：`"length mismatch for {seq_name}, expected {seq_len} but got {seen_bp}"`。
说明输入序列有未被图覆盖的区段（可能有 novel 段未补全）。

**Level 4 — 路径步长之和**（每序列）： 把 `path_v` 里每个 node_id 的长度
（`select(id+1) - select(id)`）求和，必须等于 `seq_len`。失败信息最详细：
`"path step length mismatch for {seq_name}: expected {seq_len} bp but path steps sum to {path_step_len} bp ({n_nodes} nodes, seen_bp={seen_bp}, nodes: {first3}...{last3})"`。
`nodes` 字段在节点数 ≤10 时全列，否则只列前 3 + 后 3，避免错误信息过长。

**校验顺序的设计**：Level 1→2 在遍历 overlap 段时同步做（一次循环），Level 3→4 在遍历结束后 做。
Level 2 最贵（逐碱基比对），但因为 Level 1 已保证 overlap 唯一，Level 2 的内层循环不会 重复执行。
这 4 级校验共同保证：**输出的 GFA P 路径 walk 过的节点序列，逐碱基 reconstruct 回输入序列**，是
seqwish 作为"图物化器"的正确性契约。

**对 pgr 的启示**：pgr 的 `pgr paf graph` 目前只在 `emit_gfa` 里做 Level 2 的等价检查
（节点序列与输入序列比对），Level 1/3/4 缺失。建议补齐：

- Level 1：在 `compact_nodes` 阶段加 overlap 唯一性断言（用 `bail!` 不用 panic），早于 Level 2
  暴露问题。
- Level 3：在路径写入前检查 `seen_bp == seq_len`，捕获 novel 段补全逻辑的遗漏。
- Level 4：路径步长求和校验，捕获节点边界标记错误（如 `seq_id_bv` 漏标导致节点合并）。
- 错误信息格式：照搬 seqwish 的 `{seq_name} pos {q} -> graph pos {p}` 风格，用户能直接 grep
  到出错位点。

## 5. 与 pgr 隐式图路线的对照

| 维度       | seqwish（显式物化）                  | pgr（隐式图）                    |
|------------|--------------------------------------|----------------------------------|
| 输入       | PAF + 序列                           | MAF→PAF + 序列                   |
| 传递闭包   | **一次性全图** DSU                   | **查询时** BFS，按需局部         |
| 等价类表达 | 图序列的一个碱基                     | 不物化，对齐区间即隐式等价类     |
| 数据结构   | iitree + DSU + 位向量                | coitrees + HashSet/Vec           |
| 输出       | GFA（S+L+P）                         | BED/PAF（query），GFA（graph / to-gfa） |
| 适用场景   | 全图分析、归一化、可视化             | 单 locus 查询、区域 MSA          |
| 规模上限   | 受图序列长度限制（Gbp 级）           | 受对齐索引大小限制（可分片）     |
| 重复处理   | `--repeat-max` / `--min-repeat-dist` | `--min-len` / `--merge-distance` |

**核心差异**：seqwish 的传递闭包是**全局、一次性**的——算出全部等价类再写图。 pgr 的传递闭包是
**局部、按需**的——每次查询从一个区间出发 BFS，只算相关等价类。这两种粒度对应不同的应用场景：
全图统计 vs 单点查询。

**seqwish 的 spanning tree 优化对 pgr 有直接借鉴价值。** pgr 现在的隐式图查询对每个 起点都做 BFS，
如果对齐网络稠密（如 4 万大肠杆菌 K=50 的稀疏对齐仍有 135k 条边），单次 BFS 可能横跳很多序列。
若 pgr 在加载 PAF 阶段预计算一棵最大权生成树，查询时优先沿树边走，可显著减少 BFS 的边遍历数。
这是 [[paf-pangenome.md]] 可考虑的优化项。

## 6. 对 pgr 各版本的启示

### 6.1 query / to-bed（坐标输出）

- **PosT 编码**：pgr 的 `pgr paf query` 若要支持反链投影，可借鉴 `make_pos_t` 把方向 打包进 u64，
  单棵区间树同时存正反链对齐。
- **SparseBitVec**：pgr 处理 4 万大肠杆菌时，序列边界用 `SparseBitVec`（只存 1-bit 位置）
  比位向量省内存且 select O(1)。

### 6.2 graph（粗全局 GFA，✅ 已实现）

- **已实现**：`pgr paf graph [-f refs.fa] --min-var-len 100`，输出 GFA v1.0（S/L/P）；`-f` 可选，拓扑模式零序列依赖。
  `src/libs/paf/graph.rs` 470 行引擎 + `src/cmd_pgr/paf/graph.rs` CLI 包装，5 单元 + 7 集成测试。
- **算法骨架**：seqwish 风格段级 DSU（CIGAR 切分 → 段对 → DSU 传递闭包 → 节点序列 → 路径
    - novel 段补全 → 边派生 → GFA 输出），简化版（无 spanning tree 优化，等价类规模小）。
- **`--min-var-len 100` 过滤**：在 CIGAR 切分阶段即过滤（indel < 阈值不切分）， 比 seqwish 在
  `write_graph_chunk` 后过滤更早，避免无效段产生。对应 minigraph 的粗框架哲学。
- **简化项**（相对 seqwish）：无 disk-backed interval tree / SparseBitVec / lock-free DSU，
  路径方向恒 `+`（反向已翻转坐标到正链），rGFA SN/SO/SR tag 已补全（见 [[paf-pangenome.md]] §3.3）。
- **与 seqwish 的关键差异——零序列依赖**：seqwish 的 GFA 输出中每个节点序列（S 行）对应传递闭包
  的一个碱基等价类，必须从 `seqidx.at(offset)` 取原始碱基，因此**必须**有序列索引。pgr graph 的节点
  是段级（segment-level，CIGAR `=`/`X`/`M` 的累计长度），拓扑（边界、边、路径）完全从 PAF 坐标推断；
  S 行序列可填 `*` 并标注 `LN:i:` 长度，实现**拓扑模式零序列依赖**。这是 pgr 粗图能快速构建（无需
  加载 GB 级 FASTA）而 seqwish 必须先建序列索引的根本原因。
- **磁盘后端兜底（未启用）**：4 万大肠杆菌全图可能超 RAM，`AdaptiveTree` 的 disk-backed 模式
  是现成的兜底方案，待规模验证后再引入。

### 6.3 to-gfa（局部精细 GFA）

- **phase 1b orphan recovery 的思路**：pgr 的局部 GFA 从一个 region 出发，BFS 发现的 等价类可能不完整
  （只覆盖部分序列）。seqwish 的 orphan recovery 循环（按序列查 iitree 补漏）可直接用于 pgr 的局部
  GFA 完整性保证。
- **流水线写入**：`write_graph_chunk` 在独立线程跑、主线程算下一个 chunk 的模式， 对 pgr 处理大区域
  MSA 时同样适用。

### 6.4 工程细节

- **零 panic 原则**：seqwish 的
  [compact.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/compact.rs)有
  `panic!("Overlap count mismatch")`，pgr 应改为 `bail!` + 诊断信息，符合 CLAUDE.md 的稳定性要求。
- **进度日志**：seqwish 每阶段都打 `%` 进度，pgr 处理大规模数据时应沿用此风格。
- **tempfile 管理**：[tempfile.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/tempfile.rs)
  提供 `set_dir` / `set_keep_temp`，pgr 的 `pgr paf` 若产生中间文件可借鉴。

### 6.5 lib.rs：FFI 边界与 C/C++ 互操作

[lib.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/lib.rs) 共 1020 行，其中约 70% 是
`#[no_mangle] pub extern "C"` 的 FFI 包装函数——这是 seqwish 从 C++ 原版重写到 Rust 的 产物：Rust
实现作为库被 C++ 主程序调用，FFI 层保持与原 C++ 头文件兼容。

#### 6.5.1 lib.rs 的两层结构

- **顶层模块声明**（21 个 `pub mod`）— 对应 §1 的 6 阶段实现： `seqindex`（阶段 1）、`alignments`/
  `paf`/`cigar`（阶段 2）、`transclosure`/`dset64`/`dset64_asm`/`dset64_unsafe`（阶段 3）、
  `compact`（阶段 4）、`links`（阶段 5）、`gfa`（阶段 6），加上 `pos`/`dna`/`intervaltree`/
  `mmap`/`tempfile`/`time`/`utils`/`sxs`/`version` 等基础设施模块。模块划分与算法阶段一一对应，
  是阅读源码的天然地图。
- **FFI 包装层**（约 700 行）— 把 Rust API 翻译成 C ABI，按模块分组： `tempfile_*`（5 个）、
  `pos_*`（8 个）、`dna_*`（3 个）、`cigar_*`（5 个）、`mmap_*`（2 个）、`paf_row_*`（13 个）、
  `sxs_*`（9 个）、`alignments::match_hash`/`keep_sparse`、`utils::file_exists`/`handy_parameter`、
  `time::time_since_epoch_ms`。

#### 6.5.2 opaque handle 模式

FFI 层用 6 个 opaque struct 把 Rust 资源暴露给 C++，确保生命周期由 C++ 侧控制：

- `CigarHandle { cigar: Vec<CigarOp> }` — CIGAR 向量句柄
- `SeqIndexHandle { seqidx: Arc<SeqIndex> }` — 序列索引句柄（`Arc` 共享所有权）
- `IITreeHandle { iitree: Arc<RwLock<AdaptiveTree<u64, PosT>>> }` — 节点/路径区间树句柄（读写锁）
- `AlnIITreeHandle { iitree: Arc<Mutex<AdaptiveTree<u64, PosT>>> }` — 比对区间树句柄（互斥锁）
- `PafRowHandle { row: PafRow }` — PAF 行句柄
- `SxsHandle { aln: SxsAlignment }` — SXS 比对句柄

约定：`*_parse`/`*_new` 返回 `*mut Handle`（堆分配，`Box::into_raw`），`*_free` 用 `Box::from_raw`
回收。字段访问器（如 `paf_row_query_start`）通过 `unsafe { &(*handle).row }`借用，返回原始类型或
`*mut c_char`（C 字符串，调用方负责 `free`）。

**关键认识**：`IITreeHandle` 用 `RwLock` 而 `AlnIITreeHandle` 用 `Mutex`，反映了 §2.3
的读写模式差异——节点/路径树构建后只读（多线程并发查询），比对树写入阶段需独占。`SeqIndexHandle` 用
`Arc` 是因为序列索引在多线程查询时共享。

#### 6.5.3 FFI 的安全边界

FFI 层是 `unsafe` 的集中地，但遵循三条纪律：

1. **入口校验空指针** — 所有函数开头 `if ptr.is_null() { return ...; }`，避免 C++ 传 NULL 崩溃。
2. **`CStr::from_ptr` 包裹** — 字符串转换用 `to_str()` + `match` 处理非 UTF-8，不直接 `unwrap`。
3. **所有权显式转移** — `mmap_open_rust` 用 `mem::forget(handle)` 把 `MmapHandle` 所有权转给 C++，
   避免 Rust 析构关闭 fd；`mmap_close_rust` 重建 `MmapHandle` 后调用 `mmap_close`。

但仍存在风险：`dna_reverse_complement` 用 `std::ptr::copy_nonoverlapping` 直接写 C++ 缓冲区， 若 C++
传入的 `out` 缓冲区不足 `len` 字节会 UB——这是 FFI 的固有代价，文档化由调用方保证。

#### 6.5.4 对 pgr 的启示

1. **pgr 不需要 FFI 层** — pgr 是纯 Rust 项目，无 C++ 主程序，不应模仿 seqwish 的 FFI 包装。 pgr 的
   `lib.rs` 应保持简洁的 `pub mod` 声明 + 公共类型，不引入 `extern "C"`。
2. **模块划分与算法阶段对应** — seqwish 的 21 个 `pub mod` 与 §1 的 6 阶段一一对应，
   是阅读源码的天然地图。pgr 的 `libs/paf/` 已有类似实践（`index.rs`/`query.rs`/`graph.rs`/
   `to_gfa.rs` 按处理阶段划分），可继续保持。
3. **opaque handle 模式的 Rust 纯净版** — 若 pgr 未来需要把图构建引擎抽象成可替换后端 （如 seqwish
   风格 vs impg 风格），可借鉴 handle 模式的"所有权显式转移"思想，但用 `Arc<dyn GraphEngine>` +
   trait object 而非裸指针。
4. **`RwLock` vs `Mutex` 的读写模式区分** — seqwish 对只读树用 `RwLock`、写树用 `Mutex`
   是并行区间树的标准实践。pgr 的 `libs/paf/index.rs` 若未来引入并行查询，可借鉴此区分。

### 6.6 基础设施模块速览

除前述核心算法模块外，seqwish 还有 9 个基础设施模块。它们功能单一但各有工程细节值得记录：

#### 6.6.1 tempfile.rs：全局临时文件管理

[tempfile.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/tempfile.rs) 用
`once_cell::sync::Lazy<Mutex<TempFileState>>` 维护全局状态，提供 `create`/`remove`/`set_dir`/
`set_keep_temp`/`cleanup` 接口。关键设计：

- **目录选择优先级**：`set_dir` 显式设置 > `TMPDIR` 环境变量 > `/dev/shm`（Linux tmpfs，
  优先选以利用内存文件系统）> 当前工作目录。`/dev/shm` 的偏好是为了让大中间文件（如拼接 序列、
  `seq_v`）走 RAM-disk，避免磁盘 I/O 瓶颈。
- **mkdtemp + mkstemps**：先 `libc::mkdtemp` 建唯一父目录（避免多进程冲突），再 `libc::mkstemps`
  在其中建带后缀的文件。fd 建后立即 `File::from_raw_fd` + `drop` 关闭，只保留路径——后续按需
  reopen。
- **Drop 清理 + 显式 cleanup**：`TempFileState::Drop` 在程序退出时扫 `filenames` 集合 +
  父目录残余文件，`remove_file` + `remove_dir`。`cleanup()` 供多次构建之间显式调用 （main.rs 在每次
  graph build 后调）。
- **`keep_temp` 调试开关**：`--keep-temp` 设 true 后 Drop 不删文件，便于调试中间产物。

**对 pgr 的启示**：pgr 的 `pgr paf` 若产生大中间文件（如 graph 的 `seq_v`、`node_iitree`），
可借鉴这个模式——全局 `Lazy<Mutex<State>>` + `mkdtemp` 隔离 + `/dev/shm` 优先。但 pgr 现用
`tempfile` crate（标准库生态）已够用，seqwish 自己造轮子是 C++ 遗产，不必照搬。

#### 6.6.2 mmap.rs：手写内存映射

[mmap.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/mmap.rs) 直接封装 `libc::mmap`，
提供 `MmapHandle { ptr, fd, size }` + `mmap_open`/`mmap_close`。关键参数：

- **PROT_READ | PROT_WRITE + MAP_SHARED**：读写共享映射，对映射区的写会写回文件 （`msync` 后落盘）。
  用于 `seq_v` 这种需要边写边读的大文件。
- **madvise(MADV_WILLNEED | MADV_SEQUENTIAL)**：提示内核"我会顺序访问且即将需要"， 内核可预读 +
  释放已读页。对 `seq_v` 的顺序扫描场景很重要。
- **Drop 自动 munmap + close**：`MmapHandle` 的 `Drop` 调 `mmap_close`，幂等设计 （`ptr.is_null()`
  判空），多次 close 安全。

**对 pgr 的启示**：pgr 现用 `memmap2` crate（跨平台、安全封装），不需要手写 `libc::mmap`。 但
`MADV_SEQUENTIAL` 的 hint 思路可借鉴——`memmap2` 的 `MmapOptions::advise` 支持设置。

#### 6.6.3 dna.rs：256 字节查表的碱基操作

[dna.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/dna.rs) 用一个 256 字节的
`DNA_COMPLEMENT` 表覆盖所有单字节值的互补：

- **IUPAC 码全覆盖**：A/T/C/G/N + R/Y/S/W/K/M/B/D/H/V 双向，大小写各一份。
- **GCSA 特殊字符**：`$` ↔ `#`（GCSA 索引的 stop/start sentinel），`-` ↔ `-`（gap 自反）。 这是
  seqwish 与 GCSA/xg 互操作的需要，pgr 不需要这两个。
- **`reverse_complement_in_place`**：单遍扫描，`swap + complement` 同步做，奇数长度中间 元素单独
  complement。比"先 reverse 再 complement"少一次遍历。

**对 pgr 的启示**：pgr 的 `src/libs/dna.rs` 已有类似 256 表实现，思路一致。seqwish 的 in-place
版本对大序列反补（如 `seq_v` 处理反链 query 时）省内存，pgr 可借鉴。

#### 6.6.4 cigar.rs：极简 CIGAR 解析

[cigar.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/cigar.rs) 定义
`CigarOp { len: u64, op: u8 }`（`#[repr(C)]` 为 FFI 兼容），提供 `cigar_from_string`/
`cigar_to_string` 两个函数。解析是手写状态机（digit 累积 → 遇 alpha 切换 → flush），不依赖 regex。
支持所有 CIGAR op（M/I/D/N/S/H/P/X/=），空串返回空 vec。

**对 pgr 的启示**：pgr 已有 `src/libs/alignment.rs` 处理 CIGAR，功能更全（含逐碱基 walk）。 seqwish
的极简版只切 op 不解释语义，因为对齐解析在 `alignments.rs` 里做。pgr 不需要借鉴。

#### 6.6.5 paf.rs：PAF 行 + 多文件 spec 解析

[paf.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/paf.rs) 定义 `PafRow`（12
必填 字段 + CIGAR 可选字段），`from_line`/`to_string` 往返。额外提供 `parse_paf_spec`：解析
`"file1:weight1,file2:weight2,..."` 格式，支持单文件无权重（默认 0）、混合权重、空字段跳过、
无效权重跳过。这个 spec 用于 `-p` 参数指定多个 PAF 文件及其权重。

**对 pgr 的启示**：pgr 的 `src/libs/paf/` 已有更完整的 PAF 解析（含 MAF→PAF 转换）。 seqwish 的
`parse_paf_spec` 多文件权重思路可参考——pgr 若支持"多 PAF 加权合并"可借鉴此格式。

#### 6.6.6 sxs.rs：seqwish 私有对齐格式（不借鉴）

[sxs.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/sxs.rs) 定义 `SxsAlignment`， 解析
6 行格式的对齐记录（A/I/M/C/Q 各一行 + 空行分隔）。这是 seqwish C++ 原版的私有输入 格式，Rust
移植保留以兼容旧数据。pgr 用 PAF/MAF，**不需要这个模块**（见 §7）。

#### 6.6.7 utils.rs：file_exists + handy_parameter

[utils.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/utils.rs) 两个函数：

- `file_exists(filename) -> bool`：用 `libc::stat` 检查文件存在（比 `std::fs::metadata` 少一层抽象，
  C++ 遗产）。pgr 应用 `std::fs::metadata` 或 `Path::exists`。
- `handy_parameter(value, default) -> f64`：解析 `"10k"/"2.5m"/"1g"` 这类带 k/m/G 后缀的
  数字，失败返回 default。这是 seqwish CLI 参数解析的辅助（如 `-t 10k`）。pgr 用 `clap` 的
  `value_parser` + 自定义 parser 更规范，但 k/m/G 后缀的需求可借鉴此实现。

#### 6.6.8 time.rs：单函数时间戳

[time.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/time.rs) 只有
`time_since_epoch_ms() -> u64`，封装 `SystemTime::now().duration_since(UNIX_EPOCH).as_millis()`。
用于进度日志的耗时统计。`expect("System time before Unix epoch")` 在时钟异常时 panic，这是合理的
（时钟倒流说明系统严重故障），pgr 可照此处理。

#### 6.6.9 version.rs：构建时注入版本号

[version.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/version.rs) 用
`option_env!("SEQWISH_GIT_VERSION")` 在编译时注入 git 版本（build.rs 或 `cargo:rustc-env`），
`CODENAMES: HashMap` 维护版本号→代号映射（如 `v0.7.6` → `Temporaneo`）。提供 `get_version`/
`get_release`/`get_codename`/`get_short` 四个查询函数 + 对应 FFI 导出。

**对 pgr 的启示**：pgr 现用 `clap` 的 `crate_version!` + Cargo.toml 版本号，比 seqwish 的
`option_env!` + build.rs 更简洁。seqwish 的代号映射是 PGGB 项目的版本文化（意大利语形容词），pgr
不需要。

### 6.7 工程模式总结

通读 22 个源文件后，提炼出 seqwish 在工程层面的 6 个通用模式，对 pgr 大规模数据处理有直接参考价值：

#### 6.7.1 原子操作三件套

seqwish 把 `std::sync::atomic` 用到了极致，三个层次的原子操作覆盖不同场景：

- **`AtomicU64` + `fetch_or`** — `AtomicBitVec`（§4.1.1）的并行位标记。x86 编译成 `LOCK OR`，
  单指令原子地"置位并返回旧值"，用于 BFS 去重（transclosure）和并行边界标记 （compact）。比
  `Mutex<HashSet>` 快 10-100x。
- **`AtomicU128` + `compare_exchange`** — `DisjointSets`（§2.4）的无锁 union-find。 x86_64 编译成
  `CMPXCHG16B`（需 16 字节对齐），ARM64 用 LDXR/STXR pair。
- **`AtomicBool` + `AtomicU64` 计数器** — `alignments.rs` 的 `more: AtomicBool` 控制
  worker 循环（EOF 时 `store(false)`），`transclosure.rs` 的 `active_workers: AtomicU64`
  统计活跃线程数做收尾判断。比 `Mutex<bool>` + `Condvar` 轻量。

**对 pgr 的启示**：pgr 的并行 BFS（`pgr paf query`）目前用 `Mutex<HashSet>` 去重， 规模大时是瓶颈。
改用 `AtomicBitVec` + `fetch_or` 是低成本高收益的优化——只需把 `HashSet::insert` 换成
`bv.set(pos, true, Ordering::Relaxed) && !was_set`。

#### 6.7.2 流水线 + 生产者-消费者

seqwish 的 transclosure 主循环是经典的 manager-worker 流水线：

- **`crossbeam_queue::ArrayQueue`**作为无锁有界队列，`todo_in`/`todo_out` 两个队列
  分别承载"待处理"和"已发现"任务。
- **1 manager + N workers**：manager 从 `todo_in` 取任务分发给 workers，workers 调
  `explore_overlaps_discovery` 处理后把新发现的位置 push 回 `todo_in`。manager 统计
  `active_workers`，全部 worker idle 且 `todo_in` 空时收敛。
- **`write_graph_chunk` 独立线程**：phase 2 算完等价类后，`write_graph_chunk` 在独立 线程写图序列 +
  区间树，主线程同时处理下一个 chunk，实现 I/O 与计算重叠。

**对 pgr 的启示**：pgr 的 to-gfa 局部 GFA 若处理大区域 MSA，可借鉴此模式——BFS 发现和
等价类写出用双线程流水线，避免单线程的 I/O 等待。`crossbeam_queue::ArrayQueue` 是 Rust
生态成熟的无锁队列，比手写 `Mutex<VecDeque>` + `Condvar` 好。

#### 6.7.3 进度日志风格

seqwish 每个阶段都打 `%` 进度，格式统一：

```
[seqwish] seqindex: 100% (40/40 sequences) ...
[transclosure] phase 1: 45% (4500000/10000000 positions) ...
[compact] 78% (31200/40000 sequences) ...
```

关键要素：(1) 模块名前缀 `[xxx]` 便于 grep；(2) 百分比 + 绝对值 `(done/total)` 双指标； (3) 单位明示
（sequences/positions/bytes）。用 `eprintln!` 走 stderr，不污染 stdout 的 GFA 输出。

**对 pgr 的启示**：pgr 处理 4 万大肠杆菌时应有同风格日志。建议在 `src/libs/io.rs` 加
`progress_bar(total, module_name)` 辅助函数，每 N 条记录打一次（避免日志刷屏）。

#### 6.7.4 Rayon par_iter + RwLock 并行读

seqwish 的 compact/links/gfa 三个阶段都用 `(1..=n).into_par_iter().for_each(...)` 并行
处理每条序列/每个节点，共享数据用 `Arc<RwLock<AdaptiveTree>>` 包裹：

- **写阶段独占**（`Mutex`）：`alignments.rs` 写 `aln_iitree` 用 `Arc<Mutex<...>>`， 多 worker
  串行写入。
- **读阶段共享**（`RwLock`）：compact/links/gfa 读 `path_iitree`/`node_iitree` 用
  `Arc<RwLock<...>>`，多线程 `read()` 并发查。
- **不可变共享**（`Arc` 无锁）：`SeqIndex` 构建后只读，直接 `Arc<SeqIndex>` clone 到各线程，
  无锁开销。

这个"写用 Mutex、读用 RwLock、不可变用 Arc"的三层选择是 Rust 并发数据结构的最佳实践。

**对 pgr 的启示**：pgr 的 `libs/paf/index.rs` 若引入并行查询（多 region 同时 BFS），
应照此分层：索引构建后用 `Arc<Index>` 共享，BFS 队列用 `crossbeam_queue`，结果收集用
`par_iter().map().collect()`。

#### 6.7.5 全局状态 + Lazy 初始化

seqwish 用 `once_cell::sync::Lazy` 管理两类全局状态：

- **tempfile.rs 的 `TEMP_STATE`**：`Lazy<Mutex<TempFileState>>`，首次 `create()` 时初始化 临时目录，
  程序退出时 Drop 清理。
- **version.rs 的 `CODENAMES`**：`Lazy<HashMap<&str, &str>>`，首次 `get_codename()` 时建表。

`once_cell::Lazy` 比 `std::sync::OnceLock` 多了"返回值的解引用"语法糖，且稳定可用。 Rust 1.70+ 的
`OnceLock<T>` 是标准库替代，pgr 可任选。

**对 pgr 的启示**：pgr 若需要全局配置（如日志级别、临时目录），可用 `OnceLock<T>` 避免
`lazy_static!` 宏依赖。

#### 6.7.6 错误处理：io::Result + eprintln 诊断 + panic 兜底

seqwish 的错误处理分三层，对应不同严重性：

- **`io::Result<()>` + `return Err`**：可恢复错误（文件不存在、GFA 校验失败），用
  `io::Error::new(ErrorKind::Other, format!(...))` 构造带详细信息的错误。
- **`eprintln!` + 继续**：非致命异常（如某条 PAF 行解析失败），打日志后跳过该行继续处理。
- **`panic!` + 诊断信息**：不可恢复的算法 bug（如 compact 的 overlap 不唯一、time 的时钟 倒流），
  panic 前先 `eprintln!` 打印完整上下文。

**对 pgr 的启示**：pgr 的 CLAUDE.md 要求"零 panic"，前两层照搬，第三层改为 `anyhow::bail!` +
诊断信息。seqwish 的"panic 前打印完整上下文"做法（如 compact 的 overlap 详情、gfa 的 first3...last3
节点列表）是错误信息设计的范本，pgr 的 `bail!` 应 照此提供足够的定位信息。

## 7. 不打算借鉴的部分

- **SXS 格式**（[sxs.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/sxs.rs)）： seqwish
  私有的对齐格式，pgr 用 PAF/MAF，不需要。
- **dset64_asm.rs / dset64_unsafe.rs**：手写汇编和 unsafe 优化，`portable_atomic` 的 `AtomicU128`
  已足够快，pgr 用便携版 `DisjointSets` 即可。
- **`--sparse-factor` 哈希稀疏化**
  （[alignments.rs](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/alignments.rs)的
  `keep_sparse`）：seqwish 用哈希函数随机丢弃对齐，pgr 走 Mash KNN 稀疏化（见 [[ecoli-cohort.md]]），
  质量更高，不用哈希稀疏化。

## 8. 参考链接

- 源码：[seqwish-master/src/](file:///Volumes/ExtHome/Scripts/pgr/seqwish-master/src/)
- 关联文档：[[pangenome-tools.md]] §3.2（PGGB 流水线中 seqwish 的位置）、 [[impg.md]] §1.1.2
  （隐式图 vs 物化图适用边界）、[[minigraph.md]]（粗框架过滤哲学）、[[paf-pangenome.md]]（pgr
  graph / to-gfa 路线）、[[paf-pangenome.md]]（pgr 隐式图核心原则）。

