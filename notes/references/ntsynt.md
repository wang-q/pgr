# ntSynt：多基因组宏观共线性检测（minimizer 图）

> 整理于 2026-08-09，源自对 `ntSynt-1.0.8/` 目录源码的分析（ntSynt v1.0.8，Birol Lab / BC Cancer）。
> 目的：理解其"动态 minimizer 图 + 路径 → 多基因组共线性块"算法，评估用 pgr 的 `.pgi`
> （syncmer 有序排列）实现同类的多基因组 macro-synteny 检测。本文聚焦算法与数据流；
> ntJoin 依赖（`ntjoin_utils.py` / `ntjoin.py` 在 `subprojects/ntJoin/bin/`）**本 checkout 已检出**，
> 图构建/路径/去重等核心逻辑实际都在 ntJoin 侧实现，本文据此实际源码分析。
>
> **与 pgr 的关联**：`pgi build` 已产出**按 k-mer 排序的 syncmer 位置表**（`PgiQuery`），
> 而 ntSynt 的第一步（indexlr）正是"按基因组产出一张有序 k-mer 位置表"。因此 pgr 具备了
> 实现 ntSynt 式共享 minimizer 图的"原料"，缺的是跨基因组图构建与路径/共线性块逻辑。
> 详见 `design/pgi-ntsynt.md`。

## 1. 简介

`ntSynt`（**n**? **T**eam **Synt**eny）是多基因组 macro-synteny（宏观共线性）检测工具：
输入任意数量（≥2）的基因组 FASTA，输出这些基因组之间**共同存在**的共线性块（synteny block）。
它建立在 [ntJoin](https://github.com/BirolLab/ntJoin)（scaffolding 工具）的代码库之上，
核心是把"跨基因组共享的 minimizer 邻接关系"建成一张**无向图**，然后从图中找路径，
每条路径就是一个候选共线性块。

- **为什么用 minimizer**：minimizer 采样把每个基因组的 k-mer 稀疏化到 ~1/w 密度，
  使"跨基因组共享 k-mer"的比较可在内存/时间内完成；对近缘（低分歧）基因组尤其高效。
- **为什么是图而非 pairwise 比对**：天然支持 ≥2 基因组，把"多对一"的复杂关系折叠成
  一张共享邻接图，路径上的每个 minimizer 同时在所有基因组里出现。
- **性能参考**（README）：人类基因组（0.1% 模拟分歧，2 基因组，3 Gbp）26 min / 34 GB；
  大猿基因组（4 基因组，3 Gbp）48 min / 32 GB；蜜蜂（11 基因组，0.44 Gbp）15 min / 4 GB。
- **发表**：Coombe et al. 2025, BMC Biology（DOI 10.1186/s12915-025-02455-w）。
- **依赖**：btllib（Bloom filter / NtHash）、python intervaltree/pybedtools/ncls/igraph、
  snakemake、seqtk、samtools。

## 2. 核心概念 (Key Concepts)

### 2.1 参数与默认值

| 参数 | 默认 | 含义 |
| :--- | :--- | :--- |
| `-k` | 24 | minimizer k-mer 长度 |
| `-w` | 1000 | minimizer 窗口大小 |
| `-t` | 12 | 线程 |
| `--fpr` | 0.025 | common Bloom filter 假阳性率 |
| `--hashes` | 3 | Bloom filter 哈希函数数 |
| `-d/--divergence` | 必填 | 基因组间近似最大分歧（%），用于推默认块参数 |
| `-b/--block_size` | 依分歧 | 最小共线性块长度 (bp) |
| `--indel` | 依分歧 | indel 检测阈值 (bp) |
| `--merge` | 依分歧 | 共线块合并最大间距 (bp，或 `Nw` 表示 N 倍 w) |
| `--w_rounds` | 依分歧 | 递减的细化窗口序列 |

**分歧 → 默认参数**（`ntSynt:101-111`）：

| 分歧范围 | block_size | indel | merge | w_rounds |
|---|---|---|---|---|
| < 1% | 500 | 10000 | 10000 | 100 10 |
| 1%–10% | 1000 | 50000 | 100000 | 250 100 |
| > 10% | 10000 | 100000 | 1000000 | 500 250 |

> **顶层 CLI 约束**（`ntSynt` 顶部，L101-116）：`--w_rounds` 的**每个值都必须 < `-w`**
> （否则 `parser.error`，L114-116），因为细化窗口只能比初始窗口小。另：手动显式传
> `--indel/--merge/--w_rounds/--block_size` 会**覆盖**按分歧推导的默认值
> （源码用 `args.indel or 默认` 的 or 短路实现）。
>
> **`-d/--divergence` 边界怪癖**：只有 **`> 100`** 才触发 `parser.error`（L107-111）；
> 而**负值**会静默落入 `< 1%` 分支（L101 `if args.divergence < 1` 同时捕获负值与 0~1），
> 不报错。pgr 移植时应显式校验 `0 <= d <= 100`，堵住这个不报错的缺口。

> **双层 CLI 架构（易混淆点）**：ntSynt 实际是**两层命令**，上面这张表属于**顶层驱动
> `bin/ntSynt`**（snakemake 驱动器：吃 FASTA + `-d/--divergence`，据此推导 indel/merge/block_size/
> w_rounds，再拼 snakemake 命令），而真正跑图算法的是**底层 `bin/ntsynt_run.py`**（图阶段，
> `-k`/`-w` 反而**必填**、`-z`=block_size 默认 500、`--bp`=indel 阈值默认 500、
> `--collinear-merge` 默认 `1w`、`--w-rounds` 默认 `[100,10]`、`-m` 方向阈值默认 90、
> `-n`=最小边权重默认 0→按基因组数）。顶层把参数**改名后**传给底层（snakemake config→
> `ntsynt_run.py` 的 argparse）：`--indel→--bp`、`--merge→--collinear-merge`、`--block_size→-z`、
> `--w_rounds→--w-rounds`、`--no-common→common=False`。pgr 移植时读 `ntsynt_run.py` 的参数即可，
> 无需复刻顶层 snakemake 编排。
>
> **参数名冲突怪癖**：顶层 `-n` 是 `--dry-run`，底层 `-n` 是最小边权重——同名不同义，跨层引用时
> 极易踩坑。

### 2.2 数据流（snakemake 管线，`ntsynt_run_pipeline.smk`）

```
每个基因组.fa
  ├─ rule faidx            → .fai（samtools faidx，供后续坐标/掩膜）
  ├─ rule make_common_bf   → <prefix>.common.bf（C++ 级联 Bloom filter）
  ├─ rule make_repeat_bf   → <prefix>.repeat.bf（实验性，默认关，repeat=False）
  ├─ rule indexlr          → <genome>.k<k>.w<w>.tsv（每基因组有序 minimizer 位置表）
  └─ rule ntsynt_synteny   → <prefix>.synteny_blocks.tsv（最终共线性块）
```

只有最后一步（`ntsynt_run.py`）是"算法"，其余是索引/过滤预处理：
1. **faidx**：建 .fai，供坐标与 pybedtools 掩膜用。gzipped 输入先用 `gunzip -c | samtools faidx -o -`。
2. **make_common_bf**（C++，见 §4.1）：找**所有基因组共有**的 k-mer，建成一个 Bloom filter，
   供 indexlr 只保留共有 minimizer。
3. **make_repeat_bf**（`ntsynt_make_repeat_bfs.py`，实验性，见 §4.2）：默认 `repeat=False` 关。
4. **indexlr**（btllib）：对每个基因组做 minimizer 采样，命令行 `--long --seq --pos`；common 开时
   加 `-s <common.bf>`、repeat 开时加 `-r <repeat.bf>`。输出 TSV 每行 `contig\t<mx:pos:seq> ...`，
   即一行一个 contig、各 minimizer 以空格分隔、每个以 `哈希:位置:序列` 三段式编码（按位置有序）。
5. **ntsynt_run.py**：读所有基因组的 minimizer TSV → 建共享 minimizer 图 → 图简化/过滤 →
   找路径 → 找共线性块 → 细化（`w_rounds` 递减）→ 合并共线块 → 输出 TSV。

> **`--no-common`**：可跳过 common BF（`common=False`），此时图包含所有 minimizer，
> 靠 `-n`（最小边权重）过滤。common BF 默认开，是省内存/加速的手段。

### 2.3 共享 minimizer 图模型（核心思想）

把"多基因组共线性"编码成一张**无向图**，来自 ntJoin：

- **顶点 (vertex)** = 一个 minimizer **k-mer 哈希**（`ntjoin_utils.vertex_index(graph, mx)`）。
- **边 (edge)** = 两个 minimizer 哈希在**同一基因组的有序 minimizer 列表里相邻**
  （consecutive），`ntjoin_utils.edge_index(graph, mx_i, mx_i_next)`。
- **边权重 (weight)** = 该"相邻关系"在多少个基因组里同时出现
  （`weights = [1]*len(FILES)`，每基因组贡献 1）。
- **过滤**：删掉 `weight < n` 的边（`filter_graph_global`，`-n` 默认 = 基因组数）。
  剩下的边只在**所有基因组**里都相邻出现——构成"共线主干"。

直觉：如果两个共享 minimizer 在所有基因组的同一条 contig 上都是紧邻的，那么这段
"双 minimizer 邻接"跨基因组保守，可视为共线锚点；把这些邻接首尾相连成路径，
路径覆盖的区间就是跨基因组的共线性块。

### 2.4 输出格式（`<prefix>.synteny_blocks.tsv`）

8 列 TSV，`block_id` 相同的行属于同一共线性块（8 列 TSV，最终输出 `get_block_string(verbose=True)`）：

```
block_id  genome  contig  start  end  strand  num_minimizers  broken_reason
```

- **列数注意**：第 8 列 `broken_reason` 只在最终 `<prefix>.synteny_blocks.tsv`（
  `refine_block_coordinates` 末轮 `verbose=True`）输出；中间产物
  `<prefix>.pre-collinear-merge.synteny_blocks.tsv` 与初轮块是 7 列（无 broken_reason）。

- `start/end`：该基因组上的块坐标（start 为 0-based 第一个 minimizer 位置，
  end = 最后一个 minimizer 位置 + k，见 §4.3 `get_block_start/end`）。
- `strand`：该基因组内块的方向（`+`/`-`，§4.3 `determine_orientations`）。
- `broken_reason`：与**上一个**块的断点原因（`None` / `id_change` / `ori_change` /
  `inconsistent_order` / `indel` / `merge`，§4.6）。

## 3. 算法详解

### 3.1 级联 Bloom filter（common BF，C++）

`src/ntsynt_make_common_bf.cpp`：找"所有基因组共有 k-mer"的内存友好方法（§4.1）。

### 3.2 初始图构建与简化/过滤（`ntsynt_synteny.py`）

1. **load_minimizers**：读各基因组 minimizer TSV，建 `list_mx_info[assembly][hash] = (contig, pos)`。
   `read_minimizers`（`ntjoin_utils.py:167`）有个**内在去重**：同一基因组内若某个 minimizer
   在多个位置出现（重复/多拷贝），它会被加入 `dup_mxs` 并从 `mx_info` **整体剔除**（不只删
   一个拷贝）——即"一遇重复即整删"，这是管线自带的重复序列过滤，无需 repeat BF 就生效。
   （去重判据用 `mx` 哈希；若传入 repeat_bf，则按 `seq` 序列命中 `repeat_bf.contains(seq)` 也加入
   `dup_mxs`。另外 `NtSyntSynteny.__init__` 会把 `FILES` 按字典序**逆序**排序，保证确定性。）
   注：common BF 的 `s` 过滤在此步前只作用到 indexlr 产出。
2. **make_minimizer_graph**（ntJoin）：先 `filter_minimizers`（`ntjoin_utils.py:152`）做
   **跨基因组 set 交集**——只保留**每个基因组都出现**的 minimizer（共线前提），再按 §2.3 建图。
   common BF 的 `s` 过滤在此步前只作用到 indexlr 产出，此处交集是最终把关。
3. **run_graph_simplification**（`ntsynt_synteny.py:586`）：简化气泡。
   - `node_partially_anchored`：顶点只有**一条** max-weight（= 基因组数）的关联边。
   - 对两端都是"度 3 + 部分锚定"的边：若 source→target 恰好有**两条**简单路径
     （一条直边 + 一条 3 节点路径），删掉中间节点（`path[1]`），把直边权重提到 max。
   - 效果：消除单点噪声（重复/错误 minimizer 造成的分叉），把"真实共线"的直连边权重拉满。
   - **门控**：受 `--simplify-graph` 控制——`ntsynt_run.py` 直接调用默认**关**，但顶层/管线默认
     **开**（除非 `--no-simplify-graph`）；细化每轮 `w_rounds` 也会再简化一次
     （`ntsynt_synteny.py:504-505`）。
4. **filter_graph_global**（`ntjoin.py:78`）：删 `weight < n` 的边。注意**守卫**：
   若 `n <= min(weights)` 则**直接返回不改图**（ntJoin 单参考场景用）；ntSynt 里各基因组权重恒为 1
   （`weights_list = [1]*len(FILES)`）、`n` 默认 = 基因组数 ≥2，故实际总是过滤。
   注意**同名的 `filter_graph_global_flag_overlaps`**（`ntsynt_synteny.py:312`）是细化末轮用的变体，
   会额外记录被删边的两端 `flagged_node_pairs`（见 §3.5），勿与全局过滤混淆。

### 3.3 路径发现与共线性块提取

- **ntjoin_find_paths()**（ntJoin）：`find_paths`（`ntjoin.py:139`）对每个连通分量调
  `find_paths_process`——**不是简单遍历，而是"迭代去分支化"**：只要分量非线性（存在度 ≥3
  的节点），就 `filter_graph` 删掉**分支节点（度>2）上权重 < 递增阈值**的边（阈值从 `-n` 起
  每轮 +1），直到分量线性；再对每个线性子分量找两个度=1 端点，用 `determine_source_vertex`
  （以**最大权重基因组**上位置最小/最大者定 source/target）定向，`get_shortest_paths` 取
  简单路径，且**仅当路径覆盖该子分量的全部节点与边**（`len(path)==len(vs)` 且
  `num_edges==len(es)`）才接受。
  - **并行怪癖**：`find_paths`（`ntjoin.py:145-150`）本会在 `t > 1` 时用
    `multiprocessing.Pool` 对分量并行找路径，但 `NtSyntSynteny.__init__` 强制
    `self.args.t = 1`（`ntsynt_synteny.py:37`）——**ntSynt 实际总是单进程找路径**；
    整条管线的并行只来自 btllib 的 indexlr 线程 + snakemake 的 job 调度。pgr 启示：
    找路径天然按连通分量可并行，移植时用 rayon `par_iter` 对分量并行即可，不必学它强制串行。
- **find_paths_synteny_blocks**（`ntsynt_synteny.py:563`）：把每条路径交给
  `find_synteny_blocks`。
- **find_synteny_blocks**（`ntsynt_synteny.py:70`）：沿路径逐个 minimizer 走：
  - `continue_block(mx)`：若该 minimizer 在**每个**基因组都映射到**当前块的 contig**
    → `extend_block` 延展；否则结束当前块、开新块。
  - 块结束时 `determine_orientations()`：若 `all_oriented`（每基因组都能定 ±）→ 保留；
    否则把该块的 minimizer 从图里删除（`to_remove_nodes`），即未定向的块丢弃。
  - 关键：**路径上的相邻 minimizer 是跨基因组共线的**，块只在"跨 contig 跳变"处断开。
- **check_for_indels**（`ntsynt_synteny.py:411`）：对每个块，相邻 minimizer 对的
  `max_difference`（各基因组相邻 minimizer 位置差的 max − min）> `--bp`（indel 阈值）
  就在该处断块，并把对应的图边删掉（`remove_flagged_edges`）。
- **filter_synteny_blocks**（`ntsynt_synteny.py:431`）：删除 minimizer 数 < 阈值（4）的块。
  **魔法数怪癖**：初轮（`ntsynt_synteny.py:649`）与细化每轮（`:514`）都硬编码 `4`
  （源码 `# TODO: magic number`），且过滤条件是 `all(len >= 4)`（每基因组的块内 minimizer
  都要 ≥4，任一不足则整块丢弃）。pgr 移植应把这个阈值参数化。

**调试/健壮性细节**：
- `--interarrivals`（`ntsynt_run.py:40`）可在初轮后导出每块相邻 minimizer 的间隔距离
  （`print_interarrivals`，`ntsynt_synteny.py:577`），用于诊断 indel/断块阈值是否合理。
- `main_synteny` 在初轮无路径时**打印错误并 `sys.exit(1)`**（`ntsynt_synteny.py:654-656`，
  "no paths found. Try adjusting the specified k/w parameters"）；并在入口校验
  `w_rounds` 无重复（`:621-623`）、`--filter` 必须配 `--repeat`（`:625-626`）。
- **igraph 版本兼容**：`run_graph_simplification` 里 `get_all_simple_paths` 的
  cutoff 参数名随 igraph 版本变化（`ntsynt_synteny.py:595`，`<1.0.0` 用 `cutoff`、否则 `maxlen`）。
  pgr 用 petgraph 无此烦恼，但提示了"上游库 API 漂移"的维护成本。

### 3.4 坐标细化（`w_rounds` 递减窗口）

初始块用的是大窗口（`-w`，如 1000）minimizer，块端点粗糙。`refine_block_coordinates`
（`ntsynt_synteny.py:497`）用**递减窗口**（如 250、100）把块端点磨细：

1. **mask_assemblies_with_synteny_extents**：把每个基因组中"已被共线性块覆盖"的区间
   （长度 > `max(2*w, w+k+1)` 的块）掩膜掉（slop `-(w+k)` 后 mask_fasta）。
2. **generate_new_minimizers**：对掩膜后的序列用更小的 `w` 重新跑 indexlr（仍过滤 common），
   得到更密集的 minimizer。
3. **find_mx_in_blocks**：收集各块两端 minimizer（terminal）与内部 minimizer（internal）。
4. **filter_minimizers_synteny_blocks**：新 minimizer 里，凡是落在块区间重叠处
   （intervaltree/NCLS）或属于 internal black_list 的去掉，只留**块与块之间**的新 minimizer。
5. **build_graph**：把这些新 minimizer 以权重 1 加进现有图（`black_list=terminal_mxs`）。
   `build_graph`（ntjoin_utils.py:83-141）在**增量重建**（图已存在）时还会调
   `check_added_edges_incident_weights`（ntjoin_utils.py:70-80）加一道守卫：
   对每个新加的边，若其任一端点的**总关联权重** `sum(incident weight) >
   sum(weights.values())*2`（即超出"全基因组满权重×2"），就把这条新边删掉——
   防止细化过程中某些 minimizer 变成超度 hub 节点。注意 `sum(weights)*2` 用的是
   各基因组恒 1 的权重表，故阈值恒等于 `2×基因组数`。
6. 重新找路径 → 断 indel → 过滤 → 输出 `<prefix>.pre-collinear-merge.synteny_blocks.tsv`。
7. 最后一轮（`new_w == w_rounds[-1]`）：额外做 `filter_graph_global_flag_overlaps` +
   `refine_graph`（§3.5），然后 `merge_collinear_blocks`（§3.6）。

> **细化的简化顺序怪癖**：每轮 `build_graph` 返回的新图存到**局部变量** `graph`，而
> `run_graph_simplification(self.graph)` 作用在**旧的** `self.graph`（不含本轮新边），随后
> `self.graph = filter_graph_global(graph)` 又用局部 `graph` 覆盖——即每轮简化结果实际被丢弃，
> 只有末轮的 `filter_graph_global_flag_overlaps + refine_graph` 真正落到最终图上。
>
> **新 minimizer 合并回位置表的语义**（`update_list_mx_info`，`ntsynt_synteny.py:302-310`）：
> 细化轮只把**通过 `filter_minimizers`（跨基因组交集）保留下来**的 minimizer 合并进
> `list_mx_info`（`valid_mxs = 各基因组保留 minimizer 的并集`），被过滤掉的就不加入。
> 这保证后续 `continue_block`/定向/断块只看"仍在所有基因组出现"的 minimizer——pgr 移植时
> 对"每轮新增 minimizer 的加入条件"要同样严谨，避免把单基因组特异的点引入位置表。
>
> **`--collinear-merge` 的 `Nw` 解析发生在构造函数**（`ntsynt_synteny.py:41-46`），不在
> argparse：`^(\d+)w$` → 乘以 `w`；`^(\d+)$` → 原样整数；否则 `ValueError`。故 `--merge 3w`
> 在顶层被解释为 `3×w` bp，与底层 `--collinear-merge` 解析逻辑一致。pgr 若支持 `Nw` 写法，
> 建议同样在参数归一化阶段（`libs/` 纯函数）统一处理，而不是散落在命令文件里。

### 3.5 末端/重叠精修（最后一轮，`refine_graph`）

`filter_graph_global_flag_overlaps`（`ntsynt_synteny.py:312`）先删 `< n` 边并记录被删边的
两端 `flagged_node_pairs`；`refine_graph`（`ntsynt_synteny.py:363`）只处理两端都是**度 1**
（终端）的 flagged 对：`erode_edges` 沿末端向里"侵蚀"，直到两端 minimizer 位置在任一
基因组里不再相距 < k（`has_overlap`），把侵蚀过程中经过的关联边删掉。作用是修剪
"末端 minimizer 重叠"造成的错误延伸（重复/拷贝边界）。

### 3.6 共线合并（`merge_collinear_blocks`）

把同 contig、同向、间距合适且非 indel 的相邻块合并成一个：

| broken_reason | 条件 |
|---|---|
| `id_change` | 任一基因组 contig id 变了 |
| `ori_change` | 任一基因组方向变了 |
| `inconsistent_order` | 任一 gap 为负（顺序不一致） |
| `indel` | `max(gap) − min(gap) > --bp − k` |
| `merge` | `max(gap) >= --merge`（collinear_merge） |

否则合并（把后块 minimizer 接到前块尾）。执行两次（`ntsynt_synteny.py:526-531`），
每次合并后按 `-z`（block_size）过滤短块。

独立脚本 `bin/ntsynt_merge_collinear.py` 复用同一套 **broken_reason 分类逻辑**，可单独对已有
TSV 调用；但它是**纯坐标层**合并，与主流程有两处差异需注意：
- indel 阈值直接用 `gap_range > --indel`（默认 50000），**不是**主流程的 `--bp − k`；
- 合并时按 strand 扩展 start/end 坐标（`+` 取后块 end、`-` 取前块 start），输出 `num_minimizers` 恒为 0，
  且 `--merge` 默认 1000000。脚本先按"首个出现的 assembly"的 (contig, start) 排序再合并，
  并把每块 broken_reason 重置为 `None` 后重算（首块恒为 `None`）。

## 4. 实现细节

### 4.1 级联 Bloom filter（`ntsynt_make_common_bf.cpp`）

- **BF 尺寸**（`approximate_bf_size`）：按 Broder & Mitzenmacher 公式
  `m_bits = −hashes·n / ln(1 − fpr^(1/hashes))`，n = 第一个基因组的碱基数。
- **级联插入**：BF1 = 基因组 1 的全部 k-mer；对基因组 i>1，只有当 k-mer 在 BF_{i-1} 里
  才插入 BF_i；把 BF_{i-1} 删掉、BF_i 变新 BF。最终 BF = "所有基因组共有"的近似集合。
- 输入基因组先 `std::sort`，保证输出与文件顺序无关（可复现）。
- `--hashes > 1` 需 btllib ≥ 1.7.8（`ntSynt:126` 有版本守卫）。

### 4.2 repeat BF（实验性，`ntsynt_make_repeat_bfs.py`）

默认**不启用**（`repeat=False`）。若开：找"在任一基因组里 ≥2 次"的重复 k-mer，建 BF；
配合 `--filter Filter` 把重复 minimizer 滤掉（用于去重复序列）。README 标注 experimental。

### 4.3 块坐标与方向（`assembly_block.py` / `synteny_block.py`）

- `get_block_start() = min(pos0, pos_last)`；`get_block_end() = max(pos0, pos_last) + k`
  （0-based 半开，end 含 k 个碱基）。
- `get_block_length() = end − start`。
- `determine_orientations()`：对每个基因组，看块内 minimizer 位置序列：
  - 全递增 → `+`；全递减 → `-`；
  - 否则按 `正序比例 ≥ -m`（默认 90%）定 `+`/`-`，否则 `?`（未定向）。
- `all_oriented()`：所有基因组都非 `?`。未定向块会被剔除（§3.3）。

### 4.4 块内节点（`synteny_block.py`）

`SyntenyBlockNode = (mx, positions)`：`mx` 是块内第 i 个 minimizer 哈希，
`positions` 是各基因组上的位置列表。块的"最小化器数"取任一基因组的 minimizer 列表长度。

### 4.5 重叠检查（`check_non_overlapping`，仅 --dev）

最终输出前用 intervaltree 检查同 contig 上块是否重叠，重叠 ≥ `-z` 时打 WARNING
（不做硬失败，仅提示）。

### 4.6 broken_reason 编码

见 §3.6 表。`get_block_string(verbose=True)` 在最终 TSV 输出第 8 列。

## 5. 与 pgr 的对比（pgi = syncmer 有序排列）

| 维度 | ntSynt | pgr 现状（pgi） | 差距 |
| :--- | :--- | :--- | :--- |
| 每基因组有序 k-mer 表 | indexlr（minimizer TSV） | `pgr pgi build`（sorted syncmer `.pgi`） | **已有** |
| 跨基因组共有 k-mer 集合 | common BF（C++ 级联） | `.pgi` 是**排序** k-mer 表，可多路归并求交集 | 需新增（比 BF 更精确，无假阳性） |
| 邻接图（顶点=共享 minimizer，边=跨基因组相邻，权重=出现基因组数） | ntJoin 图 | — | 需新增 |
| 图简化/过滤（去气泡、weight≥n） | igraph Python | — | 需新增 |
| 路径 → 共线性块 | ntSynt 逻辑 | — | 需新增 |
| indel 断块 / 方向判定 / 共线合并 | 同上 | — | 需新增 |
| 坐标细化（w_rounds 递减掩膜重采样） | 同上 | — | 可选（二期） |

**核心洞察**：ntSynt 的"输入"（有序 k-mer 位置表 + 跨基因组共有过滤）正好是 `.pgi`
的排序特性 + 现成的 `PgiQuery` 合并查询能力。pgi 用 **closed syncmer**（密度 ~2/(w+1)、
覆盖有界）替代 minimizer，采样更稳；pgr 还能直接复用 `pgi align` 的归并、`PgiMmap`
按需读、rayon 并行等既有基础设施。**缺的是"多基因组图 + 路径 + 共线性块"这一层。**

## 6. 对 pgr 的启示（移植可行性）

1. **算法体量可控**：真正新增的是图构建 + 路径遍历 + 块提取/合并，约 300–500 行 Rust；
   不需要 Bloom filter（排序归并更精确）、不需要 igraph（图很简单，邻接表即可）、
   不需要 snakemake（CLI 编排）。
2. **采样器可换**：ntSynt 用 minimizer（k=24, w=1000）；pgi 用 closed syncmer（k=40,
   s=8, w=5）。**参数口径完全不同**，移植时按 pgi 语义，不照搬 ntSynt 默认值。
3. **共享集合用排序归并**：`.pgi` 按 k-mer 排序，跨基因组共有集 = 多路归并求交集
   （O(总条目数)），比 Bloom filter 精确（无假阳性），且能同时拿到每基因组的
   (contig, pos, strand)。
4. **图是稀疏无向的**：顶点 = 共有 syncmer 哈希；边 = 相邻对；权重 = 出现基因组数。
   用 `HashMap<u128, node>` + 邻接表即可，无需 igraph。
5. **验证**：用 ntSynt 的 C. elegans demo 数据（`tests/` 下）或 E. coli 多株，
   对比 block 数量/覆盖与 ntSynt/UCSC chainnet。
6. **照搬两个"去噪"细节即可提升质量**（`read_minimizers` 一遇重复即整删 + `filter_minimizers`
   跨基因组交集）：pgi 用排序表能**精确**实现这两者（无 BF 假阳性），这是相对 ntSynt 的
   精度优势，移植时应保留并作为正确性测试基准（预期块数与 ntSynt 一致或更准）。
7. **两层 CLI 是反面教训**：ntSynt 用 snakemake 把"参数推导 + 管线编排"与"算法阶段"硬拆成
   两层、且参数改名 + `-n` 同名冲突，跨层调试困难。pgr 应保持单层 `pgr <cmd>` 接口，把
   "按分歧推导默认值"做成 `libs/` 里的纯函数，供 CLI 直接调用。
8. **与 `design/pgi-ntsynt.md` 交叉印证（方法前提与边界）**：pgr 侧阶段 0 PoC（mg1655 ×
   nissle1917）实测：两基因组**共享 syncmer 位置 99.8%+ 都落在 PAF 覆盖区**内（k 降到 16 也仅
   ~0.2% 在覆盖外）。这印证了 ntSynt 方法的前提——**共享 k-mer/minimizer 只存在于低分歧、已有
   比对覆盖的区域**；共享 minimizer 图不会"发现新共线性"，而是把已捕获的保守区以 **N>2 多基因组
   块**的形式一次性给出。故 pgr 移植 ntSynt 的**价值不在覆盖而在形态**：
   - 提供 `align pgi`（pairwise → chainnet）之外的**直接多基因组块视图**；
   - 复用其"断块（`max_difference` > `--bp`）/定向 / 共线合并（broken_reason 分类）"逻辑，作为
     `pgr` 多基因组 synteny 命令的块语义；
   - 若仅要多基因组块，可**基于 pgi 共享 syncmer + 现成 PAF** 做"图外"块合并（无需重建整张
     minimizer 图），比完整移植 ntSynt 轻量得多。高分歧基因组共享 minimizer 稀疏，ntSynt 本身
     也退化，不应当作"覆盖提升"手段。
9. **并行与确定性（来自 §3.3 的 t=1 观察）**：ntSynt 强制 `t=1` 找路径、`FILES` 逆序排序
   （`ntsynt_synteny.py:38`）保证确定性。pgr 移植时：找路径按连通分量用 rayon `par_iter`
   并行（比 ntSynt 更强），但**图构建/块输出要保持确定顺序**（如按首基因组的 (contig,start)
   排序块，见 `SyntenyBlock.__lt__`，`synteny_block.py:102-109`），以便字节级回归测试稳定。
10. **输入命名耦合（易踩坑）**：`find_fa_name`（`ntsynt_synteny.py:113-119`）从 minimizer
    TSV 文件名正则提取对应的 `.fai`/`.fa`（`.k<k>.w<w>.tsv` 前缀），命名不符直接 `sys.exit(1)`。
    pgr 无需此耦合（`.pgi` 自带坐标），移植时不要照搬"按文件名推断"这种脆弱约定。

---

*参考来源: [ntSynt GitHub](https://github.com/BirolLab/ntSynt) | [ntJoin](https://github.com/bcgsc/ntJoin) | [Coombe et al. 2025, BMC Biology](https://doi.org/10.1186/s12915-025-02455-w) | [ntSynt wiki](https://github.com/BirolLab/ntSynt/wiki/Description-of-the-ntSynt-algorithm)*
