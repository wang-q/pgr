# minigraph 分析笔记

> 整理于 2026-06，源自对 minigraph-master 全部源文件的通读。目的：理解 minigraph 的参考锚定增量图构建算法，为 pgr 的图构建路线提供参考。

本文档记录 minigraph 项目的架构、核心算法与数据结构，并分析其对 `pgr` 项目的启示。 minigraph
（Heng Li 开发，约 1.3 万行 C 代码）是一个**构建参考型泛基因组图**的工具，核心是用 rGFA
格式建模泛基因组图，通过增量映射-增强的方式把多个 assembly 合并进图。

与 `pgr` 的对照参考见 [[paf-pangenome.md]]（路线决策）、[[paf-pangenome.md]]（图构建层设计）、
[[impg.md]]（隐式图路线）、[[cactus.md]]（Caf 退火-熔化）。

---
## 1. 项目概览

### 1.1 设计哲学：参考锚定的增量图

minigraph 走的是"参考锚定 + 增量增强"路线，与 impg（隐式图）和 pggb（POA 物化图）都不同：

- **参考锚定**：第一个输入 assembly 作为参考骨架，后续 assembly 通过映射定位到现有图
- **增量增强**：每次映射后，把 query 中"映射不良"的区段作为新 segment 插入图
- **rGFA 坐标系**：每个 segment 携带 `SN:Z:`/`SO:i:`/`SR:i:` 三 tag，标记其在某个 stable sequence
  （参考路径）上的偏移，提供稳定坐标

这种设计让 minigraph **极快**（k-mer 免比对，不需要 lastz/wfmash）且**线性可扩展** （每个新
assembly 只需一次映射），但代价是**依赖参考选择**——参考选不好会导致图偏斜。

### 1.2 与 pgr 起点的差异

| 维度     | minigraph                                  | pgr                                |
|----------|--------------------------------------------|------------------------------------|
| 输入     | FASTA（需自己跑比对）                      | 已有 MAF/PAF（复用 pairwise 资产） |
| 比对方式 | minimizer 种子 + chaining                  | 已有 lastz→chain→net→axt→maf       |
| 图模型   | 显式 GFA（`gfa_t`）                        | 隐式图（PAF 边集 + 区间树）        |
| MSA 方式 | reference-guided 线性插入（`mg_path2seq`） | POA（`libs/poa/`）                 |
| 坐标系   | rGFA 三 tag（SN/SO/SR）                    | PAF target 坐标（0-based forward） |
| 外部依赖 | 零（纯 C，自包含）                         | 零（纯 Rust POA）                  |

pgr 已在 [[paf-pangenome.md §1]] 论证：pgr 不需要 minigraph 的比对能力，因为 pgr 已有更成熟的 pairwise
基础设施。本笔记关注 minigraph 的**图算法层**对 pgr 的启示。

---
## 2. 模块分层与代码地图

```
入口层    main.c / options.c            命令分发、选项
IO 层     bseq.c / gfa-io.c / format.c  FASTA/FASTQ/GFA/GAF 读写
索引层    index.c / sketch.c            minimizer 索引
映射层    map-algo.c / lchain.c / gchain1.c / galign.c / miniwfa.c
          种子→线性链→图链→精细对齐
图构建层  ggen.c / ggsimple.c / gfa-aug.c   增量增强
图算法层  gfa-base.c / gfa-ed.c / gfa-bbl.c / shortk.c
          基础操作/GWFA/Bubble/K最短路径
后处理层  gcmisc.c / cal_cov.c / asm-call.c 排序/过滤/覆盖度/变异调用
```

各文件行数（实际统计）：

- `miniwfa.c`（834 行）— mini WFA 实现
- `gfa-ed.c`（617 行）— GWFA 图编辑距离
- `ggsimple.c`（570 行）— 增量图构建核心
- `gchain1.c`（535 行）— 图 chaining DP
- `gfa-base.c`（526 行）— GFA 基础操作
- `map-algo.c`（502 行）— 序列到图映射
- `lchain.c`（441 行）— 线性 chaining
- `gfa-io.c`（395 行）— GFA 读写
- `gfa-bbl.c`（372 行）— Bubble calling
- `main.c`（301 行）— 命令分发
- `format.c`（291 行）— GAF/LC 输出
- `gfa-aug.c`（260 行）— 图增强
- `shortk.c`（251 行）— K 最短路径
- 其余均 < 250 行（`index.c` 230 / `gcmisc.c` 223 / `ggen.c` 182 / `algo.c` 194 / `cal_cov.c` 139 / `asm-call.c` 147 / `options.c` 134）

**CLI 预设**（`-x`，main.c + options.c）：`lr`（默认，k=17/w=11，长读映射）、`asm`（k=19/w=10，asm-to-ref）、`se`（k=21/w=10，单端短读）、`sr`（=se 且强制 `FRAG_MODE|FRAG_MERGE`，双端 FR，短读）、`ggs`（=asm 且自动开启 `--ggen` 增量图生成，并置 `best_n=0` 关闭 secondary 输出）。`se`/`sr` 均设 `MG_M_SR|MG_M_HEAP_SORT|MG_M_2_IO_THREADS`（短读走 heap-sort 收集种子）。图生成相关默认参数见 `mg_ggopt_init`（options.c）：`min_var_len=50`、`min_map_len=100k`、`min_depth_len=20k`、`min_mapq=5`、`match_pen=10`、`ggs_shrink_pen=9`、`ggs_min_end_cnt=10`、`ggs_min_end_frac=0.1`、`ggs_max_iden=0.80`、`ggs_min_inv_iden=0.95`，且默认开启 `MG_G_NO_QOVLP`（ggsimple 不接受 query 重叠区段）。索引默认 `bucket_bits=14`；`-c`（CIGAR）在 ggsimple 模式下会被警告推荐开启（main.c#L225）。

**未实现功能守卫**（options.c 的 `mg_opt_check`，L110-118）：若只设
`MG_M_FRAG_MODE` 而未同时设 `MG_M_FRAG_MERGE`，会在启动时打印
"the fragment-without-merge mode is not implemented" 并返回失败——即
"片段模式但不合并片段"这个组合在 minigraph 中**是未实现功能**，直接拒绝。
pgr 移植 `sr` 预设时无需复刻这条路径（`sr` 恒同时置两者，见上）。

---
## 3. 核心数据结构

### 3.1 `gfa_t`（GFA 图）

定义在 `gfa-priv.h`，操作在 `gfa-base.c`：

- `gfa_seg_t` 数组：节点（segment），字段 `seq`/`len`/`rank`/`snid`/`soff`/`name`，另含 `utg`（unitig 列表）与 `aux`（任意 tag 序列，`gfa_aux_t`）
- `gfa_arc_t` 数组：有向边，字段 `v_lv`（高 32 位 = 头顶点 id，低 32 位 = arc 长度 `lv`，打包便于按头排序）/`w`（尾顶点）/`rank`/`ov`/`ow`（overlap），以及位域 `link_id:61`（一对对偶 arc 共享，指向 `link_aux[]`）/`strong`/`del`/`comp`。图拓扑完全由 `arc[]` 表示，另用 `idx[]` 做按头顶点的邻接索引（`gfa_arc_a`/`gfa_arc_n`）、`link_aux[]` 存 link 的 tag
- `gfa_sseq_t` 数组：stable sequences（参考路径），字段 `name`/`rank`/`min`/`max`

**顶点编码**：`v = seg_id << 1 | strand`，每个 segment 有正反向两个顶点。 `v^1` 取反向顶点，`v>>1`
取 segment id。这个编码贯穿整个 minigraph，简化了 strand 处理。

**关键约束**：minigraph 不支持 overlap segments（`mg_gfa_overlap` 检查），所有 arc 的 `ov`/`ow`
必须为 0。这与 GFA spec 允许 overlap 的设计不同，是 minigraph 的简化。

**图加载的规范化管道**（`gfa_finalize`，gfa-base.c#L421）：`gfa_read` 读完 S/L 行后统一走
`gfa_fix_no_seg`（清掉只在 L 行出现、无 S 行的空 segment）→ `gfa_arc_sort`（按 `v_lv` 基数排序）→
`gfa_arc_index`（构建 `idx[]` 邻接索引）→ `gfa_fix_semi_arc`（补对偶弧缺失的 overlap 长）→
`gfa_fix_symm_add`（补缺失的互补弧 `w^1→v^1`，保证图斜对称）→ `gfa_arc_len`（把 `lv` 填为
`seg.len - ov`）→ `gfa_cleanup`（物理删除 `del` 弧并重建索引）。**图拓扑不变量**是"每条弧与其
互补弧成对存在"（`comp` 位标注），`link_id` 把一对对偶弧绑定共享同一 `link_aux[]` tag。
pgr 若做 GFA 物化（`pgr paf to-gfa`），这段"加载即规范化"的管道（排序+对称补边+去重+索引）
值得移植为构建后清理步骤。

### 3.2 rGFA 三 tag

rGFA 是 GFA 1.0 的扩展，给 segment 加三个 tag：

- `SN:Z:` — stable sequence name（参考路径名）
- `SO:i:` — stable offset（在参考路径上的偏移）
- `SR:i:` — rank（0=参考路径，> 0=非参考）

这三个 tag 提供**稳定坐标系**：即使图后续被增强（插入新 segment、分割旧 segment），
参考路径上的坐标仍然可追溯。这是 minigraph 区别于普通 GFA 工具的核心。

**参考的锚定方式**：首个输入既可以是 rGFA 也可以是普通 FASTA——`gfa_read`
（gfa-io.c）解析 FASTA 的 `>` 头（`gfa_parse_fa_hdr`），把每个 contig 建成一个
segment 并设 `rank=0`、`snid=contig 名`、`soff=0`，从而把参考基因组锚定为
stable sequence（`rank==0` 即参考路径）。后续插入的 segment 才取 `rank>0`。

### 3.3 `mg_idx_t`（minimizer 索引）

`mg_idx_t` 定义在 `minigraph.h`（非 mgpriv.h），操作在 `index.c`；`mg_idx_bucket_t` 定义在 `index.c` 内：

- 分桶哈希表（`mg_idx_bucket_t`），桶数 `1<<b`（`bucket_bits`，默认 14）
- 两阶段构建：先用 `mg128_t` 数组 `a[]` 收集 (minimizer, position)，再经 `worker_post`（`kt_for` 并行）把每桶转成哈希表
- **特殊编码**（`mg_idx_a2h`）：出现 1 次的 minimizer 直接在哈希 value 存位置，并置 key 的最低比特（`kh_key|1`）；出现多次的按位置写入排序数组 `p[]`，value 存 `start<<32|n`
- `gfa_edseq_t`（`mg_idx_t` 的 `es` 字段）：每个顶点 `2i`/`2i+1` 存正向与反向互补两条序列缓存，供 GWFA 使用（`gfa_edseq_init`，gfa-ed.c）

这种"出现 1 次特殊编码"的设计在 minimap2/minigraph 中一脉相承，省内存。索引按 `k`/`w`/`bucket_bits` 构建（`mg_index_core`）；**建索引前会先把所有 segment 序列统一转大写**（`mg_index`，index.c#L215）；且 `mg_gfa_overlap` 一旦检测到任何非零 overlap 的 arc 就拒绝建索引。

**minimizer 的位打包编码**（`mg_sketch`，sketch.c#L56）：每个 minimizer 是 `mg128_t{x,y}` 一个 128 位字，
`x = kMer<<8 | kmerSpan`（低 8 位记 k-mer 的 span，可 < k，用于处理 N 断点），
`y = rid<<32 | lastPos<<1 | strand`。`mg_sketch` 是标准的**对称 (w,k)-minimizer 滑窗**：
用大小为 `w` 的环形缓冲 `buf[256]` 维护当前窗口，`seq_nt4_table` 把碱基编码为 2-bit（N 记为 4 会打断窗口），
正向/反向 k-mer 同时滚动（`kmer[0]` 左移进、`kmer[1]` 右移进并补 `3^c`），取字典序小者为该窗口的 minimizer；
对"对称 k-mer"（`kmer[0]==kmer[1]`）直接跳过（无法定链）；`hash64`（sketch.c#L28）是一个
无依赖的整数散列（5 轮乘法/XOR 混合），把 2-bit k-mer 打成 64 位再截断。窗口滑出时**同时输出所有
相同最小值的重复 minimizer**（保证链上的重复不被漏）。`mg128_t` 的 `x/y` 语义在 `mg_map_frag`
（map-algo.c#L366-L368）的 collect 阶段贯穿使用，是整个映射层的数据契约。

**索引自适应的 occurrence 阈值**（`mg_opt_update`，options.c#L120）：minimizer 的重复度阈值
`occ_max1`/`lc_max_occ` 不是写死的——建索引后调用 `mg_idx_cal_quantile`（index.c#L74）统计
每个 minimizer 出现次数的分位数（0.1 分位和 `occ_max1_frac` 分位），据此把阈值抬到"能覆盖大多数
非重复 minimizer"的水平，再 clamp 到 `occ_max1_cap`。pgr 若做 k-mer/种子索引，可借鉴这种
"建索引后按实际数据统计自动调参"的做法，避免用户手调 occurrence 阈值。

### 3.4 `mg_gchains_t`（图链集合）

图映射结果的顶层容器，三层结构：

- `gc[]`（gchain）：图链，跨多个 segment
- `lc[]`（llchain）：线性链，单 segment 内的种子链
- `a[]`（anchor）：单个种子（minimizer 命中）

每层通过 `off`/`cnt` 字段引用下一层。这种"三级索引"结构在 `gcmisc.c` 的 `mg_gchain_restore_order`/
`mg_gchain_restore_offset` 中维护。

---
## 4. 核心算法流程

### 4.1 增量图构建（`mg_ggen_aug`，ggen.c）

```
for each input assembly:
  1. mg_index(g)              对当前图建索引
  2. ggen_map()               将 assembly 映射到图
  3. mg_ggsimple()            识别映射不良区段
  4. gfa_augment()            将新序列作为 segment 插入图
```

**关键点**：每加入一个 assembly 都要**重建索引**（因为图变了）。这是 minigraph 线性但
非增量的代价——索引不能复用。pgr 的 PAF 索引是静态的（构建一次查询多次），无此问题。

**两条 ggsimple 路径**：`mg_ggen_aug`（ggen.c#L89-L100）按是否启用 `-c`（CIGAR）分派——无 `-c` 走 `mg_ggsimple`（基于锚点间隔的启发式打分，ggsimple.c#L107），有 `-c` 走 `mg_ggsimple_cigar`（ggsimple.c#L392，用 CIGAR 的逐碱基比对质量过滤插入区段，插入判别更精细）。后者先把每个 gchain 的 CIGAR 拆成逐碱基的区间（`gg_count_intv`/`gg_write_intv`/`gg_score_intv`，ggsimple.c#L330 起），再用 `mg_mss_all` 找候选、`gg_merge_seg`（ggsimple.c#L378）合并相邻弱区段。main.c#L225 明确警告 "it is recommended to add -c for graph generation"。两路径最终都调 `gfa_augment` 落图（ggsimple.c#L301/#L562）。pgr 的 PAF 天然带 CIGAR，对应的是"精细路径"。

**共享机制（补齐）**：两路径都先调 `mg_gc_index`（ggsimple.c#L11）把每个 query 的 gchain 映射到 segment/query 两个区间索引（`mg_intv_index`，algo.c），并统计锚点密度 `a_dens`；随后用 `mg_mss_all`（algo.c#L40，Ruzzo-Tompa 线性时间最大评分段）识别"映射不良"的区段，做 `--gg-min-end-cnt`/`--gg-min-end-frac` 末端裁剪；之后逐候选过滤：长度差 < `min_var_len`、含 N 碱基、query/图两侧重叠数（`mg_intv_overlap`）≠ 1 的丢弃；对长度差较小的事件用 `mg_path2seq` 取 path 序列 + `mg_wfa_cmp`（algo.c#L177，miniWFA 精确比对）做一致性校验并探测 inversion（翻转后若高一致则拆成两次插入）。

**粗框架哲学（≥100bp SV 过滤）**：minigraph 在第 3 步 `mg_ggsimple` 只把长度差 ≥ `min_var_len`
的变异插入图（[ggsimple.c#L213](../../../minigraph-master/ggsimple.c#L213)，
源码默认 50，论文 L153/L384 称 100bp）。论文 L601-609 给了四条理由：

1. 图会爆炸——"composed of millions of short segments"
2. minigraph 会失败——"Not indexing minimizers across segments, minigraph will fail to seed"
3. 小变体用标准方法更易分析——"small variants are easier to analyze with the standard methods"
4. 无算法能为数百人类基因组构建这种复杂图

**关键区分**：这是**图构建层**的过滤，不是查询层。minigraph **保留**完整比对（base-level CIGAR），
只是不把小变体变成图节点。pgr 的隐式图天然避开这个问题——query–to-maf --msa 查询层全量返回同源区段
（[[paf-pangenome.md]] §2.3），graph / to-gfa GFA 物化时才需引入同等过滤（[[paf-pangenome.md]] §4.3）。

### 4.2 序列到图映射（`mg_map_frag`，map-algo.c）

```
1. collect_minimizers    对 query 做 minimizer sketch（sketch.c）
2. collect_seed_hits     在图索引中查每个 minimizer 的位置
3. mg_lchain             线性 chaining（同 segment 内，lchain.c）
4. mg_gchain1_dp         图 chaining（跨 segment，经 arc 连接，gchain1.c）
5. miniwfa / GWFA        精细对齐（segment 内 / 跨边界，miniwfa.c / gfa-ed.c）
```

**线性 chaining**（`lchain.c`）：DP 和 RMQ 两种实现，把同 segment 内的种子连成链。 **图 chaining**
（`gchain1.c`）：把线性链当节点，用最短路径算链间图距离，DP 找最优组合。**精细对齐**：segment 内用
miniwfa（WFA），跨边界用 GWFA（图扩展 WFA）。

### 4.2b 映射后处理与 MAPQ 模型（gcmisc.c）

`mg_map_frag`（map-algo.c#L340）在得到 gchain 集后，串行执行一串后处理，这是
"原始 hit → 干净 hit 集 + 质量分数"的关键步骤，也是 pgr 的 query 过滤最值得对照的部分：

1. `mg_gchain_set_parent`（gcmisc.c#L74）：按 score 从高到低，为每个 gchain 找重叠的 primary hit，
   用**未覆盖长度占比**判定是否归为 secondary（`parent` 指向 primary），并累计 `n_sub`。
2. `mg_gchain_flt_sub`（gcmisc.c#L131）：按 `pri_ratio`/`best_n` 把弱 secondary 标记为过滤
   （`flt`），同 primary 完全同区间的直接删。
3. `mg_gchain_drop_flt`（gcmisc.c#L151）：物理压缩数组（`o2n` 旧→新索引映射），同步重写
   `id/parent`。
4. `mg_gchain_set_mapq`（gcmisc.c#L191）：**MAPQ 经验模型**——`mapq = pen_cm * 40 * (1 - subsc/score) * log(score) - 4.343*log(n_sub+1)`，
   其中 `pen_cm` 是 min(score 占比, anchor 数占比) × `uniq_ratio`（primary score 占所有
   primary score 之和的比），`subsc` 是次级 hit 的最大 score。即 MAPQ 同时惩罚"次级 hit 太强"
   （`subsc/score` 接近 1）和"重复导致 uniq_ratio 低"，cap 在 60。
5. 若开 `-c`，再 `mg_gchain_cigar`/`mg_gchain_gen_ds`（map-algo.c#L475-L478）补逐碱基 CIGAR 与差异串 `ds`。

**pgr 启示**：pgr 的 `pgr paf query` 若需输出 mapq 类似的置信度，可复用这套"score 比值 + 次级
竞争 + 重复度"的三因子模型，而不必照搬 40/log 的常数。第 1-3 步的"parent/flt/压缩"模式与
[[paf-pangenome.md]] 的传递闭包去重过滤思路同构。

### 4.3 GWFA：图扩展 WFA（gfa-ed.c）

GWFA 是 WFA 的图上扩展，wavefront 推进时跨 segment 边界。核心数据结构：

- `gwf_diag_t`：对角线（vertex, diagonal, k, traceback）
- `gwf_intv_t`：对角线区间（处理 reach-end-of-vertex 的情况）
- `gwf_trace_t`：traceback 栈

**关键算法**：

- `gwf_ed_extend`：wavefront 推进，处理四种情况（中间/vertex 末/query 末/双末）
- `gwf_dedup`：wavefront 去重（interval merge + diagonal dedup）
- `gwf_prune`：剪枝（去除远落后于最远 wavefront 的对角线）

GWFA 的"forbidden bands"机制（`gwf_mixed_dedup`）用区间合并处理 vertex 边界， 这是图扩展 WFA
区别于线性 WFA 的核心。

### 4.4 Bubble calling（gfa-bbl.c）

用 Tarjan SCC 算法识别图中的 bubble 结构：

- `gfa_scc1`：单源 SCC，返回 `gfa_sub_t`（子图）
- `gfa_bubble`：遍历所有 stable sequence 的起点，找 bubble
- `bb_n_paths`：数 bubble 内的路径数（DP）

每个 bubble 记录：

- `vs`/`ve`：起止 vertex
- `ss`/`se`：起止 stable offset
- `len_min`/`len_max`：最短/最长路径长度
- `n_paths`：路径数
- `seq_min`/`seq_max`：最短/最长路径序列
- `is_bidir`：是否涉及双链（inversion）

### 4.5 K 最短路径（`mg_shortest_k`，shortk.c）

AVL 树 + Dijkstra 的 K 最短路径实现：

- 每个顶点维护大小为 `max_k` 的 max-heap，存到达该顶点的 K 条最短路径
- 用 `target_dist` + `target_hash` 支持目标导向搜索
- 返回 `mg_pathv_t[]` 回溯数组

用于图 chaining 中计算线性链之间的图距离。`MG_MAX_SHORT_K` 是上限。

### 4.6 图增强（gfa-aug.c）

`gfa_augment`（gfa-aug.c）把插入（insertion）应用到图：

1. 分割现有 segment（如果插入点在中间）
2. 创建新 segment（插入序列）
3. 更新 arc（删除旧 arc，添加新 arc）

**新 segment 命名与 tag 建立**：所有 segment（含切分出的）统一命名为 `s1/s2/...`
（`snprintf("s%d", k+1)`）。切分片段继承父 segment 的 `snid/soff/rank`，`soff` 随切分点递增；
真正插入的新序列则在此**建立 rGFA 三 tag**——`snid = gfa_sseq_add(name)`、`soff = 插入的 query 起点`
（`coff[0]`）、`rank = max_rank+1`（>0，非参考），随后 `gfa_sseq_update` 更新 stable sequence 的
min/max。纯删除（`coff[0]==coff[1]`）不建新 segment，只在两侧加一条旁路 arc。最后 `gfa_arc_sort` +
`gfa_arc_index` + `gfa_fix_multi` 重排并去重弧。

`gfa_ins_adj`（gfa-aug.c#L213）调整插入坐标：用 `gfa_ins_shrink_semi` 做带 X-drop 的端部收缩，
把插入点往内缩到与 graph 序列一致的边界，处理相邻插入的边界情况。

---
## 5. 与 pgr 路线的对照

[[paf-pangenome.md §2]] 已明确 pgr 的核心决策，下面分析 minigraph 各部分对 pgr 的适用性。

### 5.1 pgr 不复用 minigraph 的 `gfa_t` 数据结构

[[paf-pangenome.md §2]] 已分析：在 Rust 中重建 `gfa_t` 需要：

- 节点定义（CIGAR block vs reference 坐标切分）
- 边定义（strand → 四种边方向）
- 路径定义（P-line + SN/SO/SR tag 管理）
- 节点去重（`gfa_aux_update_cv`/`gfa_sort_ref_arc`）
- 自环/重复处理（`mg_gchain_set_parent`）

工作量远超 100 行，且 query 只需坐标输出不需 GFA。**结论：query 不物化 GFA，推迟到 graph / to-gfa**。

### 5.2 pgr 不复用 minigraph 的映射算法

minigraph 的映射（minimizer → linear chain → gchain → GWFA）是为"在已有图上定位 query" 设计的。pgr
的场景是**已有 pairwise 比对（MAF/PAF）**，不需要重新做比对。

[[paf-pangenome.md §1.2]] 已论证：pgr 不需要 `--sparsify`，不需要 wfmash，不需要 minimizer chaining，
因为 MAF 里的每对已经跑过 pairwise 了。

### 5.3 minigraph 对 pgr 仍有价值的部分

#### (1) Bubble 模型作为查询层后处理过滤

minigraph 的 bubble（[gfa-bbl.c](../../../minigraph-master/gfa-bbl.c)）
用 Tarjan SCC 识别，[asm-call.c](../../../minigraph-master/asm-call.c)
基于 bubble 做变异调用。pgr 虽然不构建 GFA，但**PAF 传递闭包的连通分量**在概念上等价于 bubble。

[[paf-pangenome.md §6.4]] 已记录 Caf 的 melting 过滤维度可作为传递闭包后处理。 minigraph 的 bubble 提供
**正交视角**：

- Caf 是离线全局过滤（图构建时）
- minigraph bubble 是在线局部结构（图查询时）
- pgr 的 BFS 传递闭包结果天然就是"隐式 bubble"

可借鉴 minigraph bubble 的指标作为 pgr 传递闭包的过滤维度：

- `n_paths`（路径数）→ pgr 的 `--min-degree N`
- `len_min`/`len_max`（长度区间）→ pgr 的 `--min-chain-length N`
- `is_bidir`（双链）→ pgr 的 inversion 标注

**注意**：这些是传递闭包的**后处理过滤**，不是 BFS 本身的中断条件（查询时无法做全图 SCC）。

**`--call` 的 bubble 变异调用（asm-call.c）**：minigraph 还提供 `--call` 把 bubble 落成变异——
`mg_call_asm`（asm-call.c#L21）先 `gfa_bubble` 找出 bubble 并把每个 segment 标注 `bid`/`is_stem`/`is_src`，
再遍历 gchain：凡"两段 stem 之间夹一段非 stem segment"即识别为一个候选变异（含相邻 stem 的纯删除），
用 query/图两侧的 `mg_intv_overlap`（≠1 丢弃）判 orthology，最后按 `strand` 输出
`stable\tss\tse\t>seg...<seg\t:glen:strand:qname:qs:qe` 的 BED 行；`is_src` 用于解析反向折叠的
inversion（`bid` 相同则看谁靠 src）。这展示了"bubble 是结构，变异调用是语义"的完整链路：
**先 SCC 找拓扑结构，再结合 gchain 映射证据做变异解释**。pgr 的 PAF 隐式图没有 GFA，等价物是
"传递闭包连通分量 + 每个分量内的路径计数/长度差"，`--call` 的思路可平移为
"对每个连通分量，用覆盖该分量的 PAF 行数判定是否 orthologous"。

#### (2) `mg_path2seq` 的 reference-guided 思路

`mg_path2seq`（[ggen.c](../../../minigraph-master/ggen.c)）本质是
"在参考序列上按位置依次插入 query 序列段"，不是 all-vs-all MSA。算法循环：

```
while (1) {
  1. 找 rs ≤ r ≤ re 的得分最高 chain
  2. 有 → 写 ref 片段 + query 序列，前进 v
  3. 无 → 写剩余 ref 片段，结束
}
```

[[paf-pangenome.md §1.2]] 已指出这启示 pgr：当 cohort 有明确 reference 时，
**reference-guided 线性 MSA 比 POA 更快且无分支膨胀**。pgr 可在 `pgr paf to-maf --msa` 中根据是否有
`--reference` 参数选择后端：

- 有 reference → `fas multiz`（banded DP，reference-guided）
- 无 reference → `fas consensus`（SPOA，无参考）

#### (3) K 最短路径（`mg_shortest_k`）对图距离的启发

[shortk.c](../../../minigraph-master/shortk.c) 的 K 最短路径用于 图
chaining 中计算线性链之间的图距离。pgr 的 BFS 传递闭包目前只做"可达性"，不计算"图距离"。

如果未来需要**按同源紧密度排序**传递闭包结果，可借鉴 minigraph 的：

- `target_dist` + `target_hash` 机制（目标导向搜索）
- AVL 树 + max-heap 的 K 最短路径实现

但这是**远期需求**，query 不需要。

#### (4) 覆盖度计算模型

[cal_cov.c](../../../minigraph-master/cal_cov.c) 的 `mg_cov_asm`（cal_cov.c#L55，asm 级批量）与
`mg_cov_map`（cal_cov.c#L8，单映射）计算：

- **segment 覆盖度**：先把每个 gchain 在 segment 上的区间投影到正向链（`rev` 位翻转），
  排序后**区间合并**累加覆盖长度，除以 segment 长（cal_cov.c#L123-L133）
- **link 覆盖度**：用 `gfa_find_arc` 找 lchain 间 arc，`++cnt_link`（cal_cov.c#L106-L116）

`--cov` 入口是 `mg_ggen_cov`（ggen.c#L104）：对每个输入重复 `ggen_map` → `mg_cov_asm`，
累加后按输入数归一化（`cov_seg[j] /= n_fn`），最后用 `gfa_aux_update_cv(g, "cf", ...)`
（gfa-base.c#L493）把覆盖度写回 segment/arc 的 `cf:f:` tag——即覆盖度是**落进 GFA 的 tag**，
而非独立输出。这是"把可量化的质量指标持久化到图里"的一个工程点，pgr 若物化 GFA 可借鉴
（把覆盖度/置信度写进 S 行/L 行的 tag，而不是单独一张表）。

pgr 的 PAF 区间树已支持区间查询，可类似地计算"每个 query 区间被多少 pairwise 比对覆盖"， 作为
**传递闭包置信度**的量化指标。这与 [[paf-pangenome.md §6.4]] 的 Degree 过滤对应。

#### (4b) 粗框架过滤的两种正交维度

minigraph 的 `--min-var-len`（默认 100）是按**变异长度**过滤的粗框架（§4.1）。seqwish 提供另一种
正交维度：`--repeat-max` / `--min-repeat-dist` 按**重复拷贝数**过滤——限制同一序列在图同一位置的
拷贝数，避免高拷贝重复把图吹爆（详见 [[seqwish.md]] §4.5、§7.2）。

pgr paf graph 物化粗 GFA 时可同时启用两种过滤：

- `--min-var-len 100`（minigraph 风格）——过滤 < 100bp 的小变体
- `--repeat-max N`（seqwish 风格）——限制重复拷贝数

两者维度不同，互不替代。详见 [[paf-pangenome.md]] §5.2。

#### (5) GAF 紧凑路径编码

[format.c](../../../minigraph-master/format.c) 的 `mg_write_gaf` 实现
**GAF 路径列的紧凑编码**：

- **紧凑模式**：当整条链都落在同一条 `rank==0` 的 stable sequence（`min==0`）上时，
  把第 6-9 列合并为 `name\tmax\tst\ten`（stable 名 + 总长 + 起止偏移），并省略显式 `plen/ps/pe`
- **展开模式**：否则逐 segment 展开为 `>seg[:st-en]`/`<seg[:st-en]`（`>`/`<` 表正/反链），
  再附显式 `plen\tps\tpe`
- **顶点坐标模式**（`--vc`，`MG_M_VERTEX_COOR`）：始终展开为纯 `>seg<seg...`，不带区间

另输出 `tp:A:P/S`（primary/secondary）、`cm:i:`（锚点数）、`s1/s2:i:`（得分/次级分）、
`dv:f:`（divergence）、`NM:i:`（错配数）等 tag；`-c` 时还有 `cg:Z:`（CIGAR）与
`ds:Z:`（差异串，含 `+`/`-`/`*`/`:` 算子）。

pgr 的 `pgr paf query` 输出同源区间列表时，可借鉴这种"能合并就合并，不能就展开"的双模式输出。

---
## 6. pgr 相对 minigraph 的独有优势

### 6.1 复用 pairwise 资产

minigraph 必须自己跑比对（minimizer chaining），pgr 复用已有 MAF/PAF。 见 [[paf-pangenome.md §1.2]]。

### 6.2 Chain/Net syntenic 验证

minigraph 没有 UCSC Chain/Net 体系，pgr 可用 Chain/Net 做同源置信度标注。 见 [[paf-pangenome.md §2.1]]。
这是"复用已有 pairwise 基础设施"的深层含义：不仅复用比对数据，还复用比对数据的**质量注释**。

### 6.3 查询层挑选

minigraph 的图是预先构建的（构建时就要决定参数），pgr 的隐式图支持查询时按 `--min-identity`
等参数动态过滤。见 [[paf-pangenome.md §2.3]]。

### 6.4 MSA 质量可能更优

pgr 的 `fas_multiz/`（banded DP）对 core 区段比 minigraph 的 reference-guided 线性插入 更精确。
见 [[paf-pangenome.md §2.4]]。

---
## 7. 结论与行动建议

### 7.1 结论

minigraph 的核心价值在于**证明了一条完整的"图构建→图映射→图增强"管道可行**， 但其 `gfa_t` 数据结构和
minimizer chaining 算法对 pgr **不直接适用**——pgr 已有更成熟的 pairwise 基础设施和 PAF 隐式图。

minigraph 值得借鉴的是**算法思想**而非具体实现：

- Bubble 作为传递闭包后处理过滤的结构化指标
- reference-guided vs POA 的 MSA 后端选择策略
- K 最短路径用于按紧密度排序同源结果（远期）
- GAF 的紧凑路径输出双模式
- 覆盖度作为置信度量化指标

**物化图的两条路径**：minigraph（增量增强，输入 FASTA 自跑 minimizer chaining）与 seqwish（PAF 诱导，
输入已有 pairwise 比对）是物化 GFA 的两条不同路径。pgr paf graph 输入是 PAF，与 seqwish 同源，因此
算法骨架（spanning tree → BFS → DSU → compact → links → GFA）直接复用 seqwish（详见
[[seqwish.md]] §7.2）；minigraph 的 `--min-var-len` 粗框架过滤哲学则作为正交补充（§4.1、§5.3(4b)）。

### 7.2 行动建议

对 `pgr paf to-maf --msa` 的影响：

- **不影响**：`pgr paf to-maf --msa` 继续走 PAF → POA → MSA 路线（[[paf-pangenome.md §3]]），约 150 行新代码
- **可借鉴**：`pgr paf query` 输出格式可参考 GAF 紧凑路径编码
- **可借鉴**：`pgr paf query --transitive` 的后处理过滤可参考 minigraph bubble 指标

对 `pgr paf to-maf` / `pgr paf to-maf --msa` 的影响：

- to-maf 评估是否引入 bubble 指标作为传递闭包过滤维度
- to-maf --msa 评估是否引入 rGFA 标准（届时再考虑 Rust 版 `gfa_t`）

---
## 8. 附录：与其他 notes 文档的引用关系

```
minigraph.md (本文档) ─ 架构参考 ──┐
    │ §5.1 gfa_t 不复用 → paf-pangenome.md §2  │
    │ §5.3(1) bubble → paf-pangenome.md §4.4       │
    │ §5.3(2) reference-guided → paf-pangenome.md §1.2
    │ §6 pgr 优势 → paf-pangenome.md §1-2          │
    │                                           │
cactus.md ────────────── 架构参考 ─────────────┤
    │ §8 Caf 退火-熔化 → paf-pangenome.md §6.4     │
    │ §3 Minigraph-Cactus → paf-pangenome.md §6.5  │
    │                                           │
impg.md ──────────────── 路线参考 ─────────────┤
    │ §4 传递闭包 → paf-pangenome.md §2.4          │
    │ §9 启示 → paf-pangenome.md                   │
    │                                           │
paf-pangenome.md (路线决策) ──────────────────────┤
    │ §1 起点差异                               │
    │ §2 核心决策                               │
    │ §4 存量资产优势                           │
    │                                           │
paf-pangenome.md (图构建层设计) ───────────────┘
    │ §1 三种图构建路线
    │ §2 query 不做 GFA
    │ §3 query 最小可行实现
```

