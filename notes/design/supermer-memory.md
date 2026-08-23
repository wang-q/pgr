# supermer 建表内存压缩设计（pgr 任务 1.7）

> 目标：降 multik reads 表（supermer 两段计数）建表瞬时峰值——MG1655
> all-masters 13.1 GB 峰值的大头。对应 `~/Scripts/anchr/notes/design/
> pgr-tasks.md` §1.7。参考实现：FastK（`FASTK-master/`，分桶 + 有界内存 +
> 表归并）。

## 1. 目标与成功标准

- 硬约束：输出与现状逐字节一致（L1 golden `f566e894…`；pgr 有
  direct/supermer 对照测试基础设施）。
- 内存：MG1655 all-masters 峰值显著下降（目标 <12 GB，A/B 定）。
- CPU：wall 不劣化超 ~10%（FastK 分桶的代价是归并，需量化）。
- 决策门槛：任一不达 → 关闭 1.7（挂起），记录数据。

## 2. 现状与打点数据

`finish_records` 两处"先全量收集再合并"：
- stage-1：per_chunk span 记录 → 合并 `records`（全量）→ radix 排序 →
  分组折叠。
- stage-2：逐 group 并行展开 → per_block keys+weights（全量）→ 合并全局
  keys+weights → radix 排序 → 相邻相同键权重累加。

`PGR_SUPERMER_TIMING` + VmRSS 打点（MG1655 20× 模拟 reads，21.4 M 碱基，
`examples/supermer_stages.rs`）：

| k | rss_sort1 | rss_expand | rss_sort2 | emitted |
| --- | --- | --- | --- | --- |
| 21 | 103 MB | 115 MB | 144 MB | 2.65M |
| 41 | 99 MB | 136 MB | 202 MB | 4.31M |
| 61 | 123 MB | 196 MB | 322 MB | 6.40M |
| 81 | 183 MB | 270 MB | 468 MB | 8.15M |

结论：**内存增长主因是 stage-2**（k=81 时 expand→sort2 涨 198 MB，对应
全量 keys+weights + radix）；stage-1 records 全量相对平稳（103–183 MB）。
真实 multik（10 chain 并发、k=81–128）按比例放大即 13.1 GB 峰值的构成。

## 3. FastK 参考模型（`FASTK-master/split.c` / `count.c`）

- super-mer 按**段最小 minimizer（mc）**路由到 NPARTS 个桶（`Min_Part`
  core-prefix trie；`NPARTS = ceil(gsize / SORT_MEMORY)`，`-M` 排序内存
  上限）→ 逐桶排序 span → 展开加权 k-mer → 排序 → 相邻相同键权重累加
  （`merge_thread`）→ 表分片拼接。任意大数据都能有界内存运行。
- pgr 的 span 切分条件（`mp < mc || force`）与 FastK 一致；**`pack_run`
  本来就算出 mc**，分桶路由几乎免费。
- FastK 表查询（`Find_Kmer`）前缀索引 + 二分单命中 → 最终表必须严格有序
  无重复键；跨桶重复由"排序后相邻合并"保证（pgr 移植时用字节级测试钉死）。
- khmer 固定内存哈希（Countgraph/Counttable）不适用：diginorm 近似/顺序
  相关，与 pgr 确定性输出冲突；anchr 历史（HashMap 5.3 GB 峰值 → radix）
  也排除哈希路线。

## 4. 方案：分批打包/展开 + 块级去重 + 归并

**修正（2026-08-24）**：先试过 FastK 式 mc 分桶（span 记录带桶号），发现
pgr 的 span 是"minimizer 非递减段"，展开后桶内 k-mer 的 minimizer 各异、
跨桶重复不可避免——桶没有"桶间独立归并"的优势（FastK 的 core-prefix
trie 也无法直接移植）。改用**分批块 + 块级去重 + 归并**：等效的内存收益
（限并发 + 去重表累积），实现更简单、无桶分布不均问题。

### stage-1：分批打包 + span 去重（`pack_spans`）

- 序列按 `PACK_CHUNK`（4096 条）切块，`PACK_BATCH`（16）块一批并行打包；
- 每块排序 span 记录 + compact（权重 1）→ span 子表，块缓冲随 map 结束
  释放；
- 全部子表 `merge_tables`（拼接 + 排序 + compact）→ 全局唯一 span 表
  （模拟 4Mb 数据：全量记录 → 唯一 span，~压缩 4×）。

### stage-2：分批展开 + 键去重（`expand_spans`）

- 全局 span 表按 `EXPAND_CHUNK`（8192 span）切块，`BATCH_BLOCKS`（8）块
  一批并行展开；
- 每块 keys+weights 排序 + compact → 键子表，块缓冲释放；
- 全部子表拼接 + 排序 + compact → 全局表（跨块重复键合并计数）。

### 内存模型

- 峰值 = PACK_BATCH 块 span 缓冲 + span 子表累积 + BATCH_BLOCKS 块键
  缓冲 + 键子表累积——不再有"全量 span 记录"和"全量展开键"双份并存。
- 代价：多轮排序（每块一次 + 归并一次），CPU 增加；PACK_BATCH 是
  内存/CPU 权衡旋钮（小→省内存、大→省 CPU）。

## 5. 字节级一致性

- 计数 = 发射窗口 multiset → 排序 → 分组累加，聚合顺序无关：相同
  multiset 必然得到相同表（§1.1 设计稿 3.1）。
- 分桶不改变发射集合：每个 span 仍恰好展开一次，桶归属只影响处理顺序；
  跨桶重复键在归并时合并计数（与现状"全局相邻合并"等价）。
- 验证：现有 `matches_direct_on_random_data` 等测试 + 新增分桶版对照
  direct 的字节级测试。

## 6. A/B 结果（2026-08-24，`examples/supermer_stages.rs` 4Mb 模拟 reads）

MG1655 20× 模拟 reads 4 Mb（85.7 M 碱基，571k 条），k=81，release：

| 配置 | total | 峰值 RSS | vs 改前 |
| --- | --- | --- | --- |
| 改前（HEAD） | 0.979 s | 1841 MB | — |
| 分区分批，PACK_BATCH=16 | 1.127 s | 1465 MB | CPU +15% / 内存 -20% |
| 分区分批，PACK_BATCH=64 | 1.100 s | 1529 MB | CPU +12% / 内存 -17% |

1 Mb 数据收益较小（-8%），4 Mb 放大到 -20%——收益随数据规模增长。

**端到端预测**：supermer 约占 multik 全流程 CPU 的 19%，lib +15% →
端到端 wall 劣化约 +3%（门槛 ≤10% 达标）；reads 表内存 -20% → 13.1 GB
峰值预计降到 ~11.5-12 GB（目标 <12 GB 边缘达标）。**决定性验证需 anchr
侧 multik G37/MG1655 all-masters A/B**（pgr 工作区改动经本地 patch 生效，
重建 anchr 后跑 asm-gate + /usr/bin/time -v）。

**输出**：全库 768 测试通过（含 direct/supermer/quality 字节级对照），
fmt/clippy 干净。

## 7. 风险

- 归并 CPU 成本（额外排序/合并）可能吃掉内存收益 → wall 门槛把关。
- PACK_BATCH 调大（省 CPU）会放大打包缓冲（费内存）——权衡旋钮已实测。
- 实现复杂，golden 保护是底线。

## 8. 决策点

1. PACK_BATCH 定值：当前 16（内存优先）；端到端 wall 若超标调 64。
2. 端到端（anchr 侧）A/B 是决定性验证；不达标 → 关闭 1.7，记录数据。
