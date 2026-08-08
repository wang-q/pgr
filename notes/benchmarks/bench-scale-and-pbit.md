# 规模延迟实测（#8）与 PBit 端到端试点（#14）

> 日期：2026-08-08。对应 `design/genome-nn-query.md` §7.4 #8/#14。

## #8 E. coli NR 规模：真实 HV 检索延迟

方法：E. coli NR 抽 494 个基因组（每 4 个取 1，跳过缺失组装），
`pgr pgi build + to-hv`（syncmer 40/8/5，D=4096），L2 归一化；
`bench hv_ann_real` 测精确扫描与 hnsw_rs HNSW 的延迟与 recall_HV@10
（无 ANI 真值，只看 HV 图检索层）。构建耗时 ≈ 0.7 s/基因组（含 pgi）。

| 规模 | 精确扫描 µs/查询 | HNSW ef=10 | HNSW ef=50 |
|---|---|---|---|
| 135 | 331 | 177 µs（recall_HV 0.993） | 312 µs（0.996） |
| 494 | 1,165 | 306 µs（0.985） | 616 µs（0.997） |

线性外推：2,115 NR ≈ 精确 5 ms、HNSW 1.3–2.6 ms/查询；15,574 全 NR
≈ 精确 37 ms、HNSW 10–20 ms。内存：4096 维 i32 HV = 16 KB/基因组，
15k ≈ 240 MB 向量 + 图边。结论：**≤10k 规模精确扫描仍在 ms 级，
维持 §6.4 建议**；HNSW 在 500 规模给出 ~3.8×（ef=10）且不损失 HV
排序（recall_HV≥0.985）。

### #8b E. coli NR 全量实测（2,088 基因组，2026-08-08）

方法：NR.lst 全量 2,115（可构建 2,088，27 个缺组装目录）；pgi build +
to-hv 流式（建一个删一个，避免 /tmp tmpfs 被 41 MB/pgi 打爆）；
`bench hv_ann_real`，真值 = 全量 Mash 距离矩阵（2.18M 对，
`mash triangle` 52 s）转相似度。

| 变体 | recall_Mash@10 | recall_HV@10 | 平均查询 µs |
|---|---|---|---|
| 精确扫描 | 0.090 | 1.000 | 5,549 |
| HNSW ef=10 | 0.095 | 0.958 | 454（12.2×） |
| HNSW ef=50 | 0.092 | 0.984 | 1,027 |

- 建库耗时：2,088 基因组 ≈ 15 min（8 并行，pgi+to-hv 流式），单线程
  外推 ~2 h；向量 2,088 × 16 KB ≈ 33 MB。
- 图检索误差仍小（recall_HV 0.958–0.984），但 **HV 排序 vs Mash 真值
  recall@10 = 0.09**（随机基线 ≈0.05）——2,088 个近同 E. coli（种内
  ANI>98%）中 HV 的 top-10 与 Mash/ANI 几乎脱钩。与 135 cohort 的
  "≥99% ANI 区间 ρ≈0.38"一致，规模扩大后更明显。
- 含义：**HV 不能做种内（≥98% ANI）精细排序**（任意规模）；种间/粗
  分层才可用。这再次指向 dist mash/frac 为种内距离主选。

## #14 PBit 试点：1 参考 + 5 近缘样本 E. coli

方法：`pgr pbit create -r ref -i 5 samples`（默认 segment 4096、
k=15、min-match 18）。样本为 NR 近缘株（ANI>99%）。

> **重要更正（2026-08-08）**：首次试点报出的"边际 delta ≈30 B/株、
> ~55,000× 压缩"是**假象**——`pbit create` 按 **contig 名**匹配参考，
> 而样本是 draft 组装（NZ_* contig 名），与参考染色体名不同，导致
> **全部 contig 被跳过**（create 日志逐条 WARN "not found in reference;
> skipping"），归档只存了参考，`pbit to-fa` 重建为 0 bp。因此 naive
> create 对跨组装样本**不存任何样本内容**，压缩比无意义。

正确用法：样本须经比对（PAF，CIGAR 编码）或 contig 名与参考一致。
核心场景（align pgi → pbit）中 PAF 是**必需**输入，不是可选项。

| 指标 | 值 |
|---|---|
| naive create（无 PAF） | 归档≈参考 2bit（1.36 MB）；样本全跳过，to-fa 重建 0 bp |
| create --paf（需 cg:Z PAF） | 见下节：bug 已修 + 约束明确；真实压缩率待完整管线 |

结论（修正后）：**pbit 的压缩前提是样本已被比对到参考**（PAF/CIGAR
或同名 contig）；对 draft 组装直接 create 会静默丢数据。核心场景
（聚类选参考 → align pgi → pbit）中比对步骤不可省略，且归档前应
校验样本重建覆盖率（`pbit to-fa` + 长度对比）作为质量门。

### #14b PAF 路径深挖：一个真 bug + 三条设计约束（2026-08-08）

**Bug（已修 + 回归测试）**：`append_sample_with_paf` 在尝试 CIGAR 编码
前先用**样本 contig 名**查**参考**键表（`contig_ref_groups`），跨组装
命名（draft NZ_* vs 参考）的样本永远被跳过，PAF 白给。修复：CIGAR
优先（走 PAF 索引，不要求名字一致），LZ-diff 兜底才按名字匹配；跳过
时不再注册空 contig。新增回归测试
`test_append_sample_with_paf_cross_assembly_names`，全部 617 个测试通过。

**约束 1（集成缺口）**：pbit 需要带 `cg:Z`（且为 `=/X/I/D` 语义，M 被
拒）的 PAF；`pgr psl to-paf` 只输出 12 列、无 cg:Z——pgr 自己的
align→pbit 闭环缺一环（需实现 PSL→PAF+cg 转换，块内按序列比对生成
`=/X`）。

**约束 2（段覆盖）**：CIGAR 编码要求每个 4096 bp 查询段被**单条 PAF
记录全覆盖**；draft 组装短 contig/碎片化比对（单条记录只盖几百 bp）
全部落回 LZ（名字不匹配则跳过）→ 重建 0。

**约束 3（段内目标）**：段的 target 区间必须落在**单个**参考 4096 段内
（`t_seg_idx_start == t_seg_idx_end`）；完整 E. coli 基因组间的重排/
大 gap 使 pgi 块散布，绝大多数段的目标跨边界 → 被拒。实测两个完整
E. coli（ANI≈99.9%）pair：1159 条 PSL 链去重后仍大量重叠，贪心单调链
生成的长 CIGAR 中 D gap 达 3 Mb，重建仍为 0。

**对设计的含义**：pbit 的 CIGAR 路径面向**近同、共线性、无大重排**的
样本（完整染色体 vs 完整染色体）；真实菌株间的反转/质粒差异会显著
降低覆盖率。选参考时应优先与样本共线性的参考；归档质量门
（`pbit to-fa` 覆盖率）应进工作流。真实压缩率数字需要先把约束 1 的
转换器实现（或接入 minimap2），列为后续任务。

### #14c 为什么真实数据压缩率拿不到（2026-08-08 最终诊断）

逐层排查后根因链（均已在单元测试层面验证 CIGAR 逻辑正确）：

1. **CIGAR 逻辑本身正确**：新增两个回归测试——多段无 gap 全覆盖 CIGAR
   三段（4096+4096+3808）完整往返 ✓；跨参考段边界的删除只跳过受影响
   段、其余段正常往返 ✓（约束 3 的精确定义）。
2. **真实完整 E. coli 对的 pgi PSL 高度碎片化**：00_3076 vs 00_3230
   的染色体产生 974–1,988 条重叠 PSL 记录、1,176 个去重块、去重后
   覆盖 7.2 Mbp（染色体仅 4.6 Mbp，说明块大量互相重叠/多映射）；
   `--min-shared` 500–5000 与 `--freq` 50–1000 均不减少碎片。
3. **朴素贪心链不可用**：从碎片块重构整 contig CIGAR 时，重复/重叠块
   使查询被覆盖 2–3 次（consume 超长）或目标坐标跳变（跨段拒绝）；
   染色体记录常被转换器整体丢弃。
4. **结论**：缺的是"链级 cg:Z PAF 生产者"——需要真正的比对链化
   （minimap2 式 chaining）或重写 psl to-paf 支持链选择 + 序列感知
   `=/X`。这是 pgr 的一个明确后续功能；在此之前的真实压缩率不可得。
   已修的 `--paf` 命名 bug 与该缺口相互独立，修复仍然有效。

**状态**：#14 = 修复 bug（完成、有测试）+ 约束/缺口（已精确记录）+
真实压缩率（待链级 cg:Z 功能，已列入后续任务）。

### #14d 最终根因：链粒度 << 段粒度（2026-08-08）

最后一轮尝试（每条 PSL 链单独出 PAF 记录，'−' 链按 contig 合并整 RC
CIGAR，修掉 `re.sub` 误删数字的转换器 bug 后 0 畸形 CIGAR）重建仍 0%。
实测 pgi 链跨度分布（00_3076 vs 00_3230，1,877 条链）：

- '+' 链 872 条：**中位数 1,142 bp**、p90 3,154 bp、最大 20 kb；
  覆盖 1.28 Mbp（查询的 ~25%）。
- '−' 链合并后整 RC 记录含巨大 I/D gap（目标坐标跳变），段内目标
  几乎必然跨参考段边界。

**结论**：pbit CIGAR 路径的"单条记录全覆盖 4096 bp 段"与 pgi 的
"~1 kb 链"粒度不匹配，真实基因组上无段可编码 → 重建恒为 0。要让
pbit 在真实数据上工作，需二选一：① pgi 链化产出长链（minimap2 式
chaining，属对齐器功能升级）；② pbit 支持跨多条记录组装段（设计
改动，需评估 delta 语义）。两条均列为后续功能；本实验给出的是精确的
失败模式与量化（链长分布），而非模糊的"不兼容"。

**#14 最终状态**：bug 修复 ✅（117 测试）；naive create 数据丢失
⚠（文档已警示 + 质量门建议）；CIGAR 约束与链粒度根因 ✅（#14b–d）；
真实压缩率 ❌（待上述任一功能，已挂 todo）。

### #14e 最根本的约束：段相位对齐（2026-08-08）

把 `--segment-size` 降到 2,048/1,024/512 重试（per-chain PAF + 修好的
二进制）：512 时仅 1 个段被编码，其余仍 0。

原因（单元测试钉死）：CIGAR 路径要求查询段的**目标区间落在单个参考
段内**（`t_seg_idx_start == t_seg_idx_end`）。查询段与参考段等长时，
这等价于**相位对齐**（target_start mod segment_size == 0）。新增测试
`test_append_sample_with_paf_indel_breaks_phase`：参考 12 kb、样本在
位置 100 插入 1 bp，CIGAR 全覆盖——结果只有第 0 段（插入点之前、
相位对齐）往返，**所有下游段因相位偏移被跳过**。

推论：
- 真实基因组上任何 indel 都会永久破坏其后所有段的相位；SNP 保持相位
  但受"单记录全覆盖段"（约束 2）与链粒度限制。E. coli 完整基因组对
  的碎片化比对使可编码段 ≈ 0。
- 段大小调参无法解决（相位要求与段长无关地成立）。
- **CIGAR 编码要真正工作，需要按链/按相位区间编码，而不是固定全局
  段**——即"pbit 跨记录/跨相位组装"（已挂 todo 的设计改动）。

**#14 终态（2026-08-08）**：bug ✅（118 个 pbit 测试）；约束 1（cg:Z
生产者）、约束 2（单记录全覆盖）、约束 3（段内目标 ⇔ 相位对齐）全部
精确化并有测试；真实压缩率 ❌，需"跨相位组装"或"长链链化"设计改动。

### #14f 路线 1 落地：LZ 兜底内容匹配化 → 真实压缩率（2026-08-08）

实现 §8.5 路线 1：名字匹配失败时，用 **canonical k-mer 倒排索引**找
内容最相似的参考段，LZ 编码（`Compressor::best_ref_group`，惰性建索引，
不改归档格式、无需 PAF/同名 contig）。新增单元测试
`test_append_sample_content_match_cross_assembly`；原 CIGAR 相位/跨段
测试更新为"LZ 回退无损恢复被 CIGAR 丢弃的段"（119 个 pbit 测试过）。

真实数据结果（`pbit create` 无 PAF，4096 段，k=15，min-match 18）：

| 样本（对参考） | 重建覆盖率 | 归档总大小 | delta≈归档−参考2bit | gzip-9 样本 | delta/gzip-9 |
|---|---|---|---|---|---|
| 00_3230（完整近缘，ANI>99%） | **100.0%** | 2,373,334 B | 972,592 B | 1,821,341 B | 53% |
| 13b5（draft 近缘，ANI 98.07） | 99.99% | 2,331,219 B | 930,477 B | 1,632,049 B | 57% |
| E. albertii（分歧，ANI≈90%） | **100.0%** | 2,524,441 B | 1,123,699 B | 1,446,101 B | 78% |

参考 2bit 部分 ≈ 1,400,742 B（self-archive 实测）；归档总大小对比
gzip(ref+sample) 约省 ~29%。

结论：
1. **路线 1 完整生效**：跨组装命名 + 无 PAF 也能 ~100% 无损归档；
   LZ 内容匹配不受 CIGAR 相位约束影响。
2. 压缩率随分歧上升（53%→78% of gzip-9），近缘场景收益最明显；
   归档为结构化（2bit 参考 + 每样本 delta，可随机访问）。
3. 构建耗时 ~13 s/样本（含索引），可接受；`--min-match` 可调
   压缩/覆盖率权衡。

**#14 最终状态（更新）**：全部实现并验证——bug ✅ + 三条约束精确化 ✅
（CIGAR 路径，约束仍在但被 LZ 回退掩盖）+ **真实压缩率 ✅（路线 1
落地）**。剩余：cg:Z 生产者（CIGAR 路径的长期优化，非阻塞）。

### #14g 多样本边际成本与覆盖率（2026-08-08）

参考 = 完整 E. coli 00_3076，增量 create（每样本单独建归档测边际）：

| 新增样本 | 类型 | 边际 delta | gzip-9 样本 | delta/gzip-9 |
|---|---|---|---|---|
| 00_3230 | 完整近缘 | 972,592 B | 1,821,341 B | 53% |
| 00_3305 | 完整近缘 | 899,588 B | 1,725,017 B | 52% |
| 01G17CRGN001 | 完整近缘 | 759,907 B | 1,497,480 B | 51% |
| 13b5 | draft 近缘 | 917,517 B | 1,632,049 B | 56% |
| 10432wF10 | draft 近缘 | 815,772 B | 1,534,504 B | 53% |
| E. albertii | 分歧（ANI≈90%） | 1,174,755 B | 1,446,101 B | 81% |

6 样本归档总 6,940,900 B（含参考），全部样本重建 ≈100%（32.38 Mbp）。

**draft 对缺失归因**：13b5 vs draft 参考（10432wF10）时丢 737 bp =
3 个 contig 的边缘段（NZ_JAGMPF010000119.1/122.1/066.1，53+240+444
bp）——样本特有序列段在参考中无内容匹配，被 LZ 回退跳过；换完整参考
（00_3076）后仅丢 53 bp。结论：**完整参考对 draft 样本覆盖更好**；
`to-fa` 覆盖率应作归档质量门（设计工作流第 6 步）。

**边际成本曲线（完整）**：近缘样本 delta = gzip-9 的 50–57%，分歧
（~90% ANI）≈ 81%；结构化归档 + 随机访问 + 无损（匹配段）是相对
gzip 的核心价值。

### #14h 严格无损核对 + Raw 回退落地（2026-08-08）

背景：用户要求"放入 pbit 必须无损"。此前只按碱基**计数**覆盖率
（99.99%）验收，不够严格；逐碱基核对（contig 名 + 顺序 + 每个位置）
暴露真实缺失：

- 12 个随机样本中 7 个有缺失（6–614 bp，中位 ~300 bp）；极端样本
  （Es_coli_188，含 66 万 N）缺失 20 万 bp（N 区整段丢失）。
- 缺失根因：内容匹配路径（无 PAF/同名）对**无参考匹配的段静默跳过**
  （`best_ref_group` 返回 None 即丢弃）；PAF 路径对"无 PAF 覆盖且无
  内容匹配"的 contig 整体跳过。LZ/CIGAR 两条路径都丢。

修复（v1006，`DeltaEncoding::Raw = 2`）：

- 未匹配段改为 **Raw 存储**（flate2 压缩原文，挂 ref_group 0，解码
  不读参考段）；两处"跳过"分支全部替换为 Raw 编码，不再静默丢数据。
- `test_raw_fallback_lossless`（原 `test_skip_unknown_contig` 反转）：
  无内容匹配 contig 必须逐碱基往返一致。
- pbit create 帮助文本同步：无参考匹配的序列 Raw 存储、归档对
  ACGTN 输入严格无损；**唯一有损点是 IUPAC 简并碱基 → N**（用户
  明确允许，2026-08-08）。

验证（修复后逐碱基核对 10 个样本，含此前缺失的 188/Mod1/WW252 等）：

- 9/10 与原始 FASTA 完全一致（名称、顺序、逐碱基）；188 的差异
  1,089 处全部为简并碱基 → N（文档声明行为），无其他差异。
- 压缩率影响：Raw 段只出现在无匹配内容，占比小；端到端 279 对
  delta/gzip 结论（决策 14）不变。

## #10 SQLite 向量存储路径实测（2026-08-08）

方法：2,088 个真实 HV（i32 4096，16 KB/个）写入 SQLite BLOB
（Python 标准库 sqlite3，近似 pgr"SQLite 存 BLOB + SIMD 扫描"路径；
扫描核为 numpy matmul，与 pgr SIMD 同量级）：

| 指标 | 值 |
|---|---|
| 入库 2,088 BLOB | 0.03 s；DB 35.4 MB（向量 33 MB + 索引开销） |
| 一次性取全部 BLOB | 32.8 ms |
| i32→f32 转换（一次性） | 35.3 ms |
| 扫描+top-10（预取后，逐查询） | **2.46 ms/查询**（numpy；pgr SIMD 同量级） |

结论：
1. **SQLite 存储开销不是瓶颈**：2,088 规模 DB 仅 35 MB，查询 ms 级；
   ≤10k"SQLite + 精确扫描"方案（§6.2/§6.4）成立。若工作流为批处理
   （一次性加载向量复用），每查询仅 ~2.5 ms（核开销），远优于 HNSW
   在该规模下的复杂度和召回风险。
2. **sqlite-vec（asg017）不再评测**：其核心是 C 扩展，用户明确换用
   纯 Rust 的 **sqlite-vector-rs 0.1.0**（PGVector-like vtab + usearch
   HNSW；MIT/Apache-2.0；注意与 rqlite 那个 Elastic-2.0 的
   "sqlite-vector" 是**不同项目**，许可不冲突）。以下为实测结果。

### #10b sqlite-vector-rs vtab 实测（2026-08-08）

方法：同 cohort（2,088 个真实 HV）→ f32 L2 归一化 → 建
`CREATE VIRTUAL TABLE emb USING vector(dim=4096, type=float4, metric=l2,
m=16, ef_construction=200, ef_search=…)` → 事务包裹批量 INSERT →
`knn_match(distance, ?) LIMIT 10` 查询；召回真值 = 精确 f32 L2 top-10。
200 个查询平均。crate 0.1.0 的 library `register()` 是 `todo!()`，bench
按 loadable-extension 入口的等价逻辑手动注册（`StandardModule::<VectorTable>`
+ `register_scalar_functions`），并以 `bundled` SQLite 编译
（`sqlite3-ext-vtab` 非 static 模式走运行时 API 表，library 用法会
SIGSEGV；主机又缺 libsqlite3-dev，故用 bundled）。

| 方案 | 构建 | 查询延迟（均值） | recall_HV@10 | 备注 |
|---|---|---|---|---|
| 精确扫描（简单 f32 循环） | — | 5.50 ms | 1.000 | 对照；numpy 版 2.46 ms（上表） |
| usearch HNSW ef10 | 1.42 s | 159 µs | 0.984 | 裸 HnswIndex（无 SQLite 层） |
| usearch HNSW ef20 | 1.44 s | 198 µs | 0.996 | 同上 |
| usearch HNSW ef50 | 1.42 s | 362 µs | 0.999 | 同上 |
| usearch HNSW ef100 | 1.38 s | 535 µs | 1.000 | 同上 |
| usearch HNSW ef200 | 1.38 s | 842 µs | 1.000 | 同上 |
| **vtab ef64（SQLite）** | 1.56 s | **1.58 ms** | **1.000** | 端到端：INSERT+查询+持久化 |
| vtab 重开（shadow 加载） | 58 ms 一次性 | 1.55 ms（warm） | 1.000 | HNSW 图从 `_index` 表反序列化 |

DB 文件：69.9 MB（向量 BLOB 33 MB + HNSW 图序列化 ~36 MB）；对照 #10a
纯 BLOB 方案 35.4 MB。

结论：
1. **vtab 端到端可用**：2,088 规模构建 ~1.6 s、查询 ~1.6 ms、召回
   1.000、持久化正确（重开 58 ms 加载后查询速度不变）。
2. **SQLite 层代价 ~3.5×**：vtab 查询 1.58 ms vs 裸 usearch（ef64
   量级 ~450 µs）——KNN 每个结果都要回查 shadow 表
   （`fetch_row_by_id`，10 次 SQL + BLOB 拷贝）。仍是精确扫描
   （5.5 ms）的 3.5× 加速；HNSW 收益在该规模被 SQLite 取行开销摊薄。
3. **召回与 §6.4 hnsw_rs 结论一致**：2,088 规模 usearch HNSW
   ef≥50 recall_HV@10 ≥ 0.999，图检索误差可忽略；ef10 开始掉点
   （0.984）。
4. **工程性记录**（写进 `design/genome-nn-query.md` §6.2 表）：
   vtab 的 shadow 表 `id` 是 AUTOINCREMENT，用户提供的 id 被忽略
   （返回 id 从 1 起，需 `id-1` 映射回 0 基下标）；library 模式
   `register()` 未实现 + `sqlite3-ext-vtab` 需 static/bundled，这些是
   0.1.0 的成熟度短板，作为 dev-dep 实验可接受，作为生产依赖要再评估。

## 复现

```bash
# #8: PGR_HV_REAL_DIR=/tmp/hv_calib/hv500 ... cargo bench --bench hv_ann_real
# #14: pgr pbit create -r ref.fa -i s1.fa ... -o a.pbit
#      用 1/2/5 样本归档大小差算边际成本
```
