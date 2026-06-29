# smoothxg 分析笔记

> 整理于 2026-06-29，源自对 `smoothxg-master/src/` 全部 31 个源文件的通读。目的：理解 smoothxg
> 如何把 GFA 图切分成 POA 块、用 SPOA/abPOA 平滑后拼接回完整图，并与 pgr 的 `paf to-gfa`
> 局部 POA 路线对照，提取可借鉴的工程细节。

## 0. 项目定位

`smoothxg` 是 PGGB 流水线中的**图归一化器**（graph normalizer）：输入一张 GFA 变异图，
输出一张"平滑后"的 GFA 图，其中 collinear 的路径段被 POA 多序列比对重新对齐，消除
seqwish 诱导阶段留下的冗余 bubble 与碎片节点。它在 PGGB 中的位置是
`wfmash`（比对）→ `seqwish`（诱导图）→ **`smoothxg`**（归一化）→ `odgi`（统计/可视化）。

一句话概括其本质：
**把整张图按路径共线性切成可 POA 平滑的块，每块独立跑 SPOA/abPOA 生成共识子图，
再用 pathfragment 映射把所有子图拼接回一张保路径的完整图。**
与 pgr/impg 的"隐式图 + 局部 POA"路线（查询时按需 POA）相对，smoothxg 是
"全图物化 + 分块 POA"路线的代表。本文档既是对其算法的拆解，也是 pgr 后续若做全图
归一化时的直接参考。

**关键认识**：smoothxg 不改变路径的碱基序列（main.cpp 末尾有严格的 path 序列校验，
不一致直接 `exit(1)`），它改变的是**图的拓扑表示**——把 seqwish 留下的"同一 locus 上的
多个等价节点"重新用 POA 对齐成更紧凑的偏序结构。换句话说，平滑前后路径序列恒等，
但图的节点数/边数/bubble 结构被重新归一化。

## 1. 整体流程（6 阶段）

[main.cpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/main.cpp) 把整个流程串成 6 个阶段，
每阶段对应一个模块。注意 smoothxg 支持**多次迭代**（`-l` 可逗号分隔多个 target POA length），
每轮迭代都完整跑一遍 prep → blocks → breaks → smooth_and_lace，前一轮的输出 GFA 作为下一轮的输入：

1. **预处理** — [prep.cpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/prep.cpp)
   `prep`：GFA 加载 → path-guided SGD 排序 → groom（消除 path 翻转）→ toposort → chop（节点切到
   `max_node_length`，默认 100 bp）。等价于 `odgi chop` + `odgi sort -p Ygs`。
2. **XG 索引** — [xg.cpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/xg.cpp)
   `XG::from_gfa`：把预处理后的 GFA 加载进 XG 索引（succinct de Bruijn graph 表示），后续
   `smoothable_blocks` 与 `smooth_and_lace` 都基于 XG 做 O(1) step 查询。
3. **块分解** — [blocks.cpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/blocks.cpp)
   `smoothable_blocks`：沿每条 path 走 step，按 `max_block_weight`（块内总碱基数）、
   `max_path_jump`（path 上相邻 step 在图排序中的跳距）、`max_edge_jump`（边的跳距）切分出
   可 POA 平滑的 block。`toposplit_block` 用并查集按弱连通分量切分跨连通区的块。
4. **块再切** — [breaks.cpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/breaks.cpp)
   `break_blocks`：在 VNTR 边界进一步切分 block，避免把高拷贝重复塞进同一个 POA。用
   mash 距离聚类 + 自相关（sautocorr）检测重复周期，参数包括 `min_copy_length`/`max_copy_length`/
   `min_autocorr_z`/`autocorr_stride`。
5. **平滑与拼接** — [smooth.cpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/smooth.cpp)
   `smooth_and_lace`：对每个 block 跑 `smooth_abpoa` 或 `smooth_spoa` 生成 POA 子图（含 consensus
   path），把子图 zstd 压缩存进 `block_graphs` 向量；同时记录 `path_mapping`（path fragment →
   block_id）供拼接阶段使用。拼接阶段遍历 `path_mapping`，把每个 block 的节点/边按 `id_mapping`
   平移后插入输出图，再按 path 顺序 append step，最后 `odgi::algorithms::unchop` 合并线性段。
6. **共识图（可选）** — [consensus_graph.cpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/consensus_graph.cpp)
   `create_consensus_graph`：若指定 `-C/--consensus-spec`，在平滑后的图上构建共识图——保留
   consensus paths 与连接它们的 link paths，丢弃中间的变异 bubble。用于生成"参考路径风格"的
   简化图。

[main.cpp#L374-L1045](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/main.cpp#L374-L1045)
给出了完整的迭代编排，包含进度日志、tempfile 管理、path 序列校验（`exit(1)` 兜底）、
consensus path 嵌入与合并（`merged_block_id_intervals_tree_vector` 处理跨 block 的 consensus
合并）等工程细节。

**未在主流程调用的模块**：[chain.hpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/chain.hpp)
定义了 `collinear_blocks`/`chains`/`superchains`（锚点链式化算法，类似 minimap2 的链式化），
但 main.cpp 没有调用它们。这应该是 smoothxg 早期版本用锚点链找共线块的遗留代码，现已被
`smoothable_blocks` 的"path step 遍历 + 跳距切分"方案取代。文档中该模块一律从略。

## 2. 关键数据结构

### 2.1 blockset_t：外存块集（mmmulti::map）

[blocks.hpp#L70-L120](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/blocks.hpp#L70-L120)
定义了贯穿全项目的块集类型：

```cpp
class blockset_t {
private:
    uint64_t _num_blocks = 0;
    mmmulti::map<uint64_t, ranked_path_range_t>* _blocks;  // key=block_id+1, value=ranked_path_range_t
    std::string _path_tmp_blocks;  // tempfile 路径
public:
    void add_block(uint64_t block_id, block_t& block);  // append 每个 path_range
    void index(uint64_t num_threads);                    // mmmulti::map::index
    block_t get_block(uint64_t block_id) const;          // values(block_id+1)
};
```

**亮点**：块集不全部驻留内存，而是用 `mmmulti::map`（memory-mapped multimap）落盘，`index`
后通过 mmap 按块 id 随机读取。这对超大图（HPRC 规模）是必备——block 数可达百万级，
全驻留会 OOM。pgr 当前 `paf to-gfa` 的 POA 块是查询时临时构造、用完即弃，无需外存；
但若未来做全图归一化，`mmmulti::map` 是值得借鉴的外存方案。

### 2.2 block_t / path_range_t：块的物理表示

```cpp
struct path_range_t {
    step_handle_t begin = {0, 0};  // [begin, end) 左闭右开
    step_handle_t end = {0, 0};
    uint64_t length = 0;
};
struct block_t {
    std::vector<path_range_t> path_ranges;  // 块内的所有 path 片段
};
```

`step_handle_t` 是 odgi/handlegraph 的标准类型，是 `(path_rank, step_rank)` 的 128 位打包。
`path_range_t` 描述"某条 path 上从 step A 到 step B 的一段"，多个 path_range 组成一个 block。
**关键认识**：block 不存节点序列，只存 path step 区间——序列通过 XG 的
`get_sequence(get_handle_of_step(step))` 按需查询。这种"引用而非复制"的设计让块分解
几乎零内存开销。

### 2.3 path_position_range_t：拼接阶段的映射记录

[smooth.hpp#L41-L61](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/smooth.hpp#L41-L61)
定义了拼接阶段的核心记录类型：

```cpp
using path_position_range_t = std::tuple<path_handle_t, uint64_t, uint64_t, path_handle_t, uint64_t>;
// get<0>: base_path（原 path）
// get<1>: start_pos（原 path 上的起始碱基坐标）
// get<2>: end_pos（原 path 上的结束碱基坐标）
// get<3>: target_path（block 内的 path，可能是 consensus path）
// get<4>: block_id
```

`smooth_and_lace` 把每个 path fragment 的映射写进 `mmmulti::set<path_position_range_t>`，
拼接阶段按 `(base_path, start_pos)` 排序后顺序遍历，把 block 内的 step 翻译回输出图。
**这是 smoothxg 的"账本"**——它记录了"原 path 的哪一段被映射到哪个 block 的哪条 path 上"，
保证平滑前后路径序列可追溯。

### 2.4 pos_t：offset + 方向的单 u64 编码（低位方向位）

[pos.cpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/pos.cpp) 定义了图位置类型：

```cpp
typedef uint64_t pos_t;
pos_t make_pos_t(uint64_t offset, bool is_rev) {
    pos_t pos = offset << 1;
    pos = (pos & ~1) | (-is_rev & 1);  // 低位存 is_rev
    return pos;
}
uint64_t offset(const pos_t& pos) { return pos >> 1; }
bool is_rev(const pos_t& pos) { return pos & 1; }
void incr_pos(pos_t& pos) { is_rev(pos) ? pos -= 2 : pos += 2; }  // ±2 步进
```

**与 seqwish 的对比**：seqwish 的 `PosT` 把方向位放低位（bit 0），offset 放高位；
smoothxg 的 `pos_t` 也是方向位放低位，但用 `<<1` 而非 `<<1 | is_rev`，实现略有差异但语义
相同。两者都支持 `±2` 步进（反链时反向步进）。pgr 若要统一反链位置编码，建议直接采用
seqwish 的 `make_pos_t` 版本（更简洁）。

注意 chain.hpp 里有另一个 `seq_pos_t`，方向位在**最高位**（MSB），与 `pos_t` 不兼容——
这是 chain 模块遗留的独立编码，主流程不用。

### 2.5 seqindex_t：SDSL 序列索引（拼接校验用）

[seqindex.hpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/seqindex.hpp) 把所有 path
序列拼接成单一字节流，用两个 SDSL 结构索引：

- `seq_begin_cbv: sdsl::sd_vector<>` — 标记每条序列的起始 offset，`rank`/`select` 支持
  O(1) 序列边界查询。
- `seq_name_csa: sdsl::csa_wt<>` — 对 `">name1 >name2 ..."` 文本建压缩后缀数组，支持
  O(m) 按名字查找序列 rank。

`seqindex_t::build_index` 从 XG graph 提取所有 path 序列拼接，`seq(name)` / `subseq(...)`
通过 mmap 拼接文件做 O(1) 随机访问。**用途单一**：仅在 main.cpp 末尾的 path 校验阶段
调用 `seqidx.seq(path_name)` 获取原序列，与平滑后图的 path 序列比对。不参与块分解或 POA。

**对 pgr 的启示**：pgr 的 `paf to-gfa` 已有 BGZF FASTA + TSV 作为序列仓库，校验时直接
从 BGZF 随机读即可，无需 SDSL 索引。但若 pgr 未来要做全图归一化且图内 path 序列需高频
随机访问，`sd_vector` + `csa_wt` 的组合比 `HashMap<String, Vec<u8>>` 更省内存。

### 2.6 consensus_spec_t：共识图规格 DSL

[consensus_graph.hpp#L48-L56](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/consensus_graph.hpp#L48-L56)
定义了 `-C` 参数的解析结果：

```cpp
struct consensus_spec_t {
    std::string basename;          // 输出文件基名
    int min_allele_len = 0;        // 保留的偏离 consensus 的最小长度
    int max_allele_len = 1e6;      // 最大 allele 长度
    std::string ref_file;          // 参考路径列表文件
    bool keep_consensus_paths;     // 是否保留 POA consensus paths
    double min_consensus_path_cov = 0;  // consensus path 最低覆盖度
};
```

`parse_consensus_spec` 解析形如 `cons,100,1000:refs1.txt:n,1000:refs2.txt:y:2.3:1000000,10000`
的 DSL（逗号分隔 spec，冒号分隔字段）。**对 pgr 的启示**：这种"单字符串编码多参数"的
DSL 在 impg 的 stage 化管道里也出现过（`gfa:cut-n=100:pggb:crush:sort`）。pgr 若要支持
复杂的图构建流水线，可借鉴此模式，但要注意 impg 的反思——DSL 解析逻辑会膨胀 main.rs，
应下沉到独立模块。

### 2.7 link_path_t / path_part_t：共识图的边

```cpp
enum path_part_t : char { begin = 'b', middle = 'm', end = 'e' };
struct link_path_t {
    std::string* from_cons_name;  // 起点 consensus path 名
    std::string* to_cons_name;    // 终点 consensus path 名
    path_handle_t from_cons_path, to_cons_path;
    path_part_t from_cons_part, to_cons_part;  // 起止在 consensus path 的哪一段
    uint64_t length;      // 核苷酸数
    uint64_t jump_length; // 在偏序中的跳距
    step_handle_t begin, end;  // off-consensus 的 step 区间
    path_handle_t path;
    uint64_t rank;
};
```

`create_consensus_graph` 用 `link_path_t` 描述"两条 consensus path 之间通过某条普通 path
的某段连接"。`path_part_t` 把每条 path 分成 begin/middle/end 三段，consensus 图保留
consensus paths 的 middle 段 + 连接它们的 link paths。这是 smoothxg 共识图的核心抽象。

## 3. 各模块详解

### 3.1 prep：图预处理流水线

[prep.cpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/prep.cpp) 的 `prep` 函数是
一条紧凑的 4 步流水线：

1. **GFA 加载** — `odgi::gfa_to_handle(gfa_in, &graph, true, num_threads, true)` 直接加载到
   odgi 内存图。
2. **path-guided SGD 排序** — `odgi::algorithms::path_linear_sgd_order` 计算"让 path 上相邻
   step 在图排序中也相邻"的节点顺序。参数：`path_sgd_iter_max=100`、`path_sgd_zipf_theta=0.99`、
   `path_sgd_eps=0.01`、`path_sgd_cooling=0.5`。`p_sgd_min_term_updates` 控制每轮迭代的最小
   更新数（默认 `1 * sum_path_step_count`）。
3. **groom + toposort** — `odgi::algorithms::groom` 消除 path 翻转（让 path 尽量走正链），
   再 `topological_order` 做拓扑排序。
4. **chop** — `odgi::algorithms::chop(graph, max_node_length, ...)` 把长节点切到 ≤100 bp，
   保证后续 POA 的块内节点粒度统一。

**关键认识**：prep 是 smoothxg 正确性的前提——path-guided SGD 让 collinear 的节点在图
排序中相邻，`smoothable_blocks` 才能沿 path step 遍历切出共线块。没有 prep，块分解会
把不相关的节点塞进同一个 POA，产出垃圾比对。pgr 的 `paf to-gfa` 走隐式图路线不需要
全局排序，但若未来做全图归一化，path-guided SGD 是必经一步。

### 3.2 blocks：块分解算法

[blocks.cpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/blocks.cpp) 的
`smoothable_blocks` 是 smoothxg 的核心算法之一。逻辑：

1. **遍历 path step** — 对每条 path，按 step 顺序遍历，累积当前 block 的 `block_handles`。
2. **切分条件** — 当 `max_block_weight`（块内总碱基数）或 `max_block_path_length`（块内
   最长 path 长度）超阈值时，finalize 当前 block 并开新 block。
3. **toposplit_block** — 用并查集（`odgi::DisjointSets`）按弱连通分量切分跨连通区的块。
   对 block 内每对相邻 step，`unite` 它们的节点 id，最后按 `find` 结果把 block 拆成子块。
4. **finalize_block** — 收集 block 内所有 handle 的 traversals（step），按
   `(path_rank, step_rank)` 排序，按 `max_path_jump`（path 上相邻 step 在图排序中的跳距）
   切分 path_range，最终输出 `block.path_ranges`。

```cpp
// finalize_block 的切分逻辑
if (path_rank(last) != path_rank(step)
    || (get_position_of_step(step) - (get_position_of_step(last) + get_length(...))
        > max_path_jump)) {
    path_ranges.push_back({step, step, 0});  // 新 range
} else {
    last = step;  // 扩展 range
}
```

**关键认识**：`toposplit_block` 是 smoothxg 处理"图排序不完美"的兜底——即使 path-guided
SGD 把不相关的节点排到一起，toposplit 也能按连通性把它们拆开。pgr 的 `paf to-gfa` 当前
用 BFS 传递闭包找同源片段，本质等价于"隐式的 toposplit"——BFS 自然只连通同源节点。

### 3.3 breaks：VNTR 边界切分

[breaks.cpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/breaks.cpp) 的 `break_blocks`
在 `smoothable_blocks` 的基础上进一步切分 block，主要解决两个问题：

1. **VNTR（串联重复）** — 高拷贝重复会把 POA 吹爆（同一序列在 block 内出现几十次）。
   `break_blocks` 用自相关（`sautocorr` 库）检测重复周期，`min_autocorr_z=5`、
   `autocorr_stride=50` 控制检测灵敏度。检测到 VNTR 就在重复边界切分 block。
2. **mash 聚类** — 对长序列（`> min_length_mash_based_clustering`，默认 200 bp）用 mash
   距离聚类，把相似度低于 `block_group_est_identity` 的序列拆到不同 block。避免把
   异源序列塞进同一个 POA。

参数 `min_copy_length`（默认 1000）/`max_copy_length`（默认 20000）限定检测的重复单元
长度范围。**对 pgr 的启示**：pgr 的 `paf to-gfa` 当前不做 VNTR 检测，对高拷贝重复区
可能产出过大的 POA 块。若实测有问题，可借鉴 smoothxg 的自相关 + mash 聚类方案。

### 3.4 smooth：POA 平滑与拼接

[smooth.cpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/smooth.cpp) 是 smoothxg
最大的源文件，包含三个核心函数：

- **`smooth_spoa`** — 用 SPOA 库对 block 跑 POA，输出 odgi 子图。支持 local/global 对齐模式
  （`-Z` 切换）、adaptive POA 参数（`-a` 按 block 内序列相似度调整打分）、padding
  （`-O` 在每端加 flanking 序列，默认 0.001 * 平均序列长度）。
- **`smooth_abpoa`** — 用 abPOA 库（ SIMD 加速）替代 SPOA，接口与 `smooth_spoa` 对称。
  `-A` 开启。
- **`smooth_and_lace`** — 编排函数：并行对每个 block 调 `smooth_spoa`/`smooth_abpoa`，
  把结果子图 zstd 压缩存进 `block_graphs[block_id]`，同时写 `path_mapping` 记录 path
  fragment → block_id 映射。若 `merge_blocks` 开启，还处理跨 block 的 consensus 合并
  （`merged_block_id_intervals_tree_vector` + `inverted_merged_block_id_intervals_ranks`）。

`build_odgi_SPOA` / `build_odgi_abPOA` 把 POA 图（spoa::Graph 或 abpoa_t）转换成 odgi
graph_t，保留 consensus path（`consensus_name`）。**关键认识**：POA 子图不是孤立的——
每条原 path 在子图里都有一条对应的 path，consensus 是额外加的一条 path。拼接阶段通过
`path_mapping` 知道"原 path 的哪段在哪个子图的哪条 path 上"，从而把子图的 step 翻译回
输出图。

### 3.5 consensus_graph：共识图构建

[consensus_graph.cpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/consensus_graph.cpp)
的 `create_consensus_graph` 在平滑后的图上构建共识图：

1. **收集 consensus paths** — 从 `consensus_path_names`（smooth_and_lace 阶段生成）+ ref_file
   指定的参考路径合并。
2. **找 link paths** — 遍历每条普通 path，找它在两条 consensus path 之间的中间段
   （`path_part_t::middle`），记录为 `link_path_t`。
3. **构建输出图** — 保留 consensus paths 的 middle 段 + 连接它们的 link paths，丢弃
   consensus paths 的 begin/end 段和变异 bubble。

`min_allele_len` 控制保留的偏离 consensus 的最小长度（小于此长度的变异被丢弃），
`max_allele_len` 控制最大 allele 长度。**对 pgr 的启示**：共识图是 smoothxg 的"参考路径
投影"——把泛基因组图简化成参考路径 + 大变异。pgr 的 `paf to-vcf` 走类似思路（POA
consensus → VCF），但输出格式不同。

### 3.6 工具模块

- **[utils.cpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/utils.cpp)** —
  `graph_deep_copy`（节点/边/path 全复制）、`handy_parameter`（解析 `1k/1m/1g` 后缀）、
  `get_block_graph`/`save_block_graph`（zstd 压缩块图的序列化/反序列化）、`modulo`（位运算
  取模，要求 d 是 2 的幂）。
- **[tempfile.cpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/tempfile.cpp)** —
  线程安全临时文件管理，`mkdtemp` + `mkstemp`，程序退出时 `atexit` 自动清理。
- **[progress.hpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/progress.hpp)** —
  `ProgressMeter` 异步进度条，原子计数 + 后台线程每 500ms 打印。
- **[zstdutil.cpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/zstdutil.cpp)** —
  ZSTD 压缩/解压字符串，含流式版本。用于 `block_graphs` 的压缩存储。
- **[dna.cpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/dna.cpp)** — 256 字节
  查表反向互补，含 IUPAC + 大小写保留。
- **[maf.hpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/maf.hpp)** —
  `maf_partial_row_t`/`maf_t`/`write_maf_rows`，MAF 输出含列宽对齐。
- **[cleanup.cpp](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/cleanup.cpp)** —
  类似 prep 但只做排序，不加载 GFA。主流程未调用。

## 4. 关键算法深入

### 4.1 块分解：path step 遍历 + 跳距切分

smoothxg 的块分解策略与 seqwish 的"传递闭包 + 节点压缩"截然不同：

| 维度     | seqwish                          | smoothxg                              |
|----------|----------------------------------|---------------------------------------|
| 输入     | PAF pairwise 比对                | GFA 图（已有节点/边/path）            |
| 单位     | 单碱基等价类                     | path step 区间（path_range_t）        |
| 切分依据 | 比对坐标的传递闭包               | path step 在图排序中的跳距            |
| 输出     | 图序列 + 节点边界                | block 集合（每块含若干 path_range）   |
| 连通性   | 并查集（DSU）合并等价类          | toposplit_block 按 DSU 拆分跨连通块   |

**关键认识**：smoothxg 不重新发现同源关系（那是 seqwish 的职责），它在已有的图拓扑上
找"可 POA 平滑的共线段"。`max_path_jump`（默认 100）是核心参数——它允许 path 上相邻
step 在图排序中跳 100 个位置仍被视为同一块，超过则切分。这等价于"允许 block 跨过小的
图排序波动，但不能跨过大段不相关节点"。

### 4.2 POA 平滑：双引擎 + padding

smoothxg 同时支持 SPOA 和 abPOA 两个 POA 引擎，通过 `-A` 切换：

- **SPOA** — 单线程，精度高，适合小 block。
- **abPOA** — SIMD 加速，多线程，适合大 block。`-T/--poa-threads` 可独立设置 POA 线程数
  （与主线程数解耦，便于控制大 block 的内存）。

`poa_padding_fraction`（默认 0.001）是个精巧的设计：在每条序列两端加 `avg_seq_len * 0.001`
的 flanking 序列，让 POA 在 block 边界有重叠，事后修剪。这避免了"block 边界处的 POA
对齐不完整"问题。`max_block_depth_for_padding_more`（默认 1000）控制深 block 不再加
padding（节省算力）。

**对 pgr 的启示**：pgr 的 `libs/poa/` 已有 SPOA 移植（参见 `docs/spoa_port.md`）。smoothxg
的 padding 机制是 pgr 当前 POA 没有的——若 pgr 的 `paf to-gfa` 在 block 边界产出不完整
对齐，可借鉴此方案。

### 4.3 path 拼接：id_mapping + path_mapping

拼接阶段（main.cpp#L614-L754）的核心是两个映射：

- **`id_mapping: Vec<uint64_t>`** — `id_mapping[block_id]` 记录该 block 的节点 id 在输出图
  中的起始偏移。拼接时 `smoothed.get_handle(block.get_id(h) + id_trans, ...)` 把 block 内
  handle 翻译到输出图。
- **`path_mapping: mmmulti::set<path_position_range_t>`** — 记录原 path 的每段映射到哪个
  block 的哪条 path。拼接时按 `(base_path, start_pos)` 排序，顺序遍历，对每个 fragment
  调 `block.for_each_step_in_path(target_path, ...)` 把 step 翻译到输出图。

**关键认识**：拼接阶段有个 `assert(false)` 兜底——如果原 path 的某段没被任何 block 覆盖
（`start_pos - last_end_pos > 0`），直接 assert 失败。这说明 smoothxg 要求**全图覆盖**：
每个 path step 必须属于某个 block。pgr 的 `paf to-gfa` 走隐式图路线，没有这个约束——
未被比对覆盖的区域不出现在图里。

### 4.4 path 序列校验：硬 exit 兜底

[main.cpp#L762-L800](file:///Volumes/ExtHome/Scripts/pgr/smoothxg-master/src/main.cpp#L762-L800)
有严格的校验：

```cpp
// 对每条 path，比对原序列与平滑后序列
std::string orig_seq = seqidx.seq(smoothed->get_path_name(path));
std::string smoothed_seq;
smoothed->for_each_step_in_path(path, [&](const step_handle_t &step) {
    smoothed_seq.append(smoothed->get_sequence(smoothed->get_handle_of_step(step)));
});
if (orig_seq != smoothed_seq) {
    std::cerr << "] error! path " << smoothed->get_path_name(path)
              << " was corrupted in the smoothed graph" << std::endl
              << "original\t" << orig_seq << std::endl
              << "smoothed\t" << smoothed_seq << std::endl;
    exit(1);
}
```

**与 seqwish 的对比**：seqwish 的 `emit_gfa` 也有 4 级校验（overlap 唯一性、path 连续性、
边完备性、序列一致性），失败 `return Err`。smoothxg 用 `exit(1)` 硬退出，且打印完整序列
对比。pgr 的"零 panic"原则（CLAUDE.md）要求此处改成 `bail!` + 诊断信息，但 smoothxg 的
"打印完整序列对比"做法值得借鉴——pgr 的 `bail!` 也应包含 `(path_name, orig_seq, smoothed_seq)`
便于定位问题。

## 5. 对 pgr 的启示

### 5.1 路线对照：隐式图 vs 全图归一化

| 维度       | pgr `paf to-gfa`（隐式图 + 局部 POA）   | smoothxg（全图归一化）              |
|------------|-----------------------------------------|-------------------------------------|
| 输入       | PAF + BGZF FASTA                        | GFA 图                              |
| 触发       | 查询时按需                              | 离线全图处理                        |
| POA 范围   | 查询区间内的同源片段                    | 整图的 collinear block              |
| 输出       | 局部 GFA                                | 全图 GFA                            |
| 排序       | 不需要（BFS 自然有序）                  | path-guided SGD（必备）             |
| 内存       | 查询区间大小                            | 全图 + block_graphs（zstd 压缩）    |
| 适用场景   | 单 locus 查询、稀疏查询                 | 全图统计、genotyping、可视化        |

**关键认识**：pgr 当前路线（to-gfa 局部 GFA）与 smoothxg 不是替代关系，而是**互补**——
单 locus 查询用 pgr 的隐式图，全图分析用 smoothxg 风格的归一化。pgr 的 `paf graph`/
`paf stat` 已经走"粗粒度全图物化"路线（仅 ≥100 bp SV 切分节点），若未来需要"精细全图
归一化"，smoothxg 的块分解 + POA 平滑是直接参考。

### 5.2 可直接借鉴的工程细节

1. **zstd 压缩块图存储** — smoothxg 用 `zstdutil::CompressString` 把每个 block 子图序列化
   后压缩存进 `Vec<string*>`，拼接时按需解压。pgr 的 `paf graph` 当前全图驻留内存，
   若图变大可借鉴此方案分块压缩。
2. **`handy_parameter` 的 k/m/g 后缀解析** — smoothxg 所有大小参数都支持 `1k/1m/1g` 后缀，
   `handy_parameter` 统一解析。pgr 的 CLI 参数若涉及大文件大小，可借鉴此工具函数。
3. **`ProgressMeter` 异步进度条** — 原子计数 + 后台线程每 500ms 打印，不阻塞主线程。
   pgr 的长任务（如 `paf index`）可借鉴此模式。
4. **`tempfile` 的 atexit 自动清理** — `mkdtemp` + `atexit` 注册清理函数，程序退出自动
   删临时目录。pgr 当前用 `tempfile::TempDir`（Rust 生态）已有类似能力，但 smoothxg 的
   C++ 实现思路值得了解。
5. **padding 机制** — POA 块边界加 flanking 序列避免对齐不完整。pgr 的 `paf to-gfa`
   若遇边界问题可借鉴。

### 5.3 不建议借鉴的设计

1. **`exit(1)` 硬退出** — smoothxg 多处用 `exit(1)` 处理校验失败，违反 pgr 的"零 panic"
   原则。pgr 应改用 `bail!` + 诊断信息。
2. **单文件超大 main.cpp** — main.cpp 1134 行，参数解析 + 流程编排 + 校验全塞一起。
   pgr 的 `paf mod.rs` 已有更清晰的子命令分发结构。
3. **chain 模块遗留代码** — `collinear_blocks`/`chains`/`superchains` 定义了但主流程不调用，
   是历史包袱。pgr 应避免保留未使用的"备选实现"。
4. **`assert(false)` 兜底** — 拼接阶段用 `assert(false)` 处理"path 未被 block 覆盖"，
   在 release build 会被编译掉。pgr 应改用 `bail!`。

### 5.4 pgr 后续路线建议

基于 smoothxg 的分析，pgr 若要做全图归一化（V7+），建议路线：

1. **复用现有 `paf graph` 的物化图** — pgr 已有 seqwish 风格的 segment-level DSU 物化
   （`libs/paf/graph.rs`），无需重新诱导。
2. **引入 path-guided SGD 排序** — 这是 smoothxg prep 的核心，pgr 可移植 odgi 的
   `path_linear_sgd_order` 算法（或用 Rust 重写）。
3. **块分解复用 BFS 传递闭包** — pgr 的 `paf query -t` 已有 BFS 同源片段发现，可直接
   作为块分解的输入，无需 smoothable_blocks 的"path step 遍历"。
4. **POA 复用 `libs/poa/`** — pgr 已有 SPOA 移植，无需引入 abPOA。
5. **拼接复用 `path_mapping` 思路** — 记录 path fragment → block 映射，拼接时按 path
   顺序翻译 step。

**关键认识**：pgr 不需要完整复制 smoothxg 的流水线——pgr 的隐式图基础设施（BFS 传递
闭包、BGZF FASTA、区间树索引）已经覆盖了 smoothxg 的部分能力。全图归一化的增量工作
主要是 path-guided SGD 排序 + 块分解策略 + 拼接映射，这三块是 smoothxg 的核心价值。

## 6. 与其他参考项目的关系

- **seqwish** — smoothxg 的上游，诱导出 smoothxg 的输入 GFA。smoothxg 不重新发现同源关系，
  只在 seqwish 产出的图上做归一化。详见 [[seqwish.md]]。
- **impg** — impg 的 `smooth` 模块（`src/smooth.rs`）直接调用了 smoothxg 风格的块分解 +
  SPOA 平滑，是 smoothxg 算法的 Rust 移植版。impg 的 stage 化管道（`gfa:pggb:crush`）里
  `pggb` 阶段就是 smoothxg 的等价物。详见 [[impg.md]] §6.2。
- **odgi** — smoothxg 的图操作全部基于 odgi 库（`odgi::graph_t`、`odgi::algorithms::*`）。
  pgr 若引入 odgi 等价物，需评估是否用 handlegraph 抽象（impg 用了 `handlegraph` crate）。
- **PGGB 流水线** — `wfmash → seqwish → smoothxg → odgi`，smoothxg 是第三步。pgr 的
  `paf` 模块覆盖了 wfmash 输出处理 + seqwish 风格诱导 + impg 风格查询，但缺 smoothxg 的
  全图归一化能力。

## 7. 关键认识汇总

1. **smoothxg 不改变路径序列，只改变图拓扑** — 平滑前后 path 序列恒等，校验失败直接
   `exit(1)`。这是 smoothxg 正确性的硬约束。
2. **块分解依赖图排序** — path-guided SGD 是 smoothxg 的前提，没有它块分解会产出垃圾。
   pgr 的隐式图路线不需要全局排序，但全图归一化必须引入。
3. **block 不存序列，只存 path step 区间** — 序列通过 XG 按需查询，块分解零内存开销。
4. **POA 子图 + path_mapping 拼接** — 每块独立 POA，拼接阶段用映射表把子图 step 翻译回
   输出图。这是 smoothxg 处理大图的核心架构。
5. **chain 模块是遗留代码** — 主流程不调用，是早期版本的锚点链式化方案，现已被
   smoothable_blocks 取代。
6. **zstd 压缩块图** — block_graphs 用 zstd 压缩存内存，拼接时按需解压。这是 smoothxg
   处理百万级 block 的内存优化关键。
7. **padding 机制** — POA 块边界加 flanking 序列避免对齐不完整，`poa_padding_fraction`
   默认 0.001。
8. **path 序列校验是硬兜底** — `exit(1)` + 打印完整序列对比，pgr 应改用 `bail!` 但保留
   诊断信息。
