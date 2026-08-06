# pgr align fill / rest + pgr psl lift（pgi 锚点 + LASTZ 补全）

> **2026-08-06 设计重构（作者决策）**：原单一 `pgr align hybrid` 拆分为两个
> 独立命令 + 一个坐标回移命令：
>
> * `pgr align fill`：锚点间 **2D gap fill**（query 双侧锚点推断裁剪）；
> * `pgr align rest`：两侧各自 **trim → excise → 求补集**（一维 runlist 运算，
>   与 `rept s-kmer` 管线同源），ref holes × 整套 query holes 比对；
> * `pgr psl lift`：range 形式 PSL → 基因组坐标回移（独立命令，可单独调试）。
>
> 拆分理由：fill 与 rest 是两种互补语义（快而准 vs 彻底），用户按场景组合；
> 坐标回移是"影响很大"且易静默出错的环节（多 contig、±链、hole 边界），
> 独立成命令便于测试与复用。
>
> 本文档前半为重构后的设计（fill/rest/lift），后半保留原 hybrid 的探索历程
> （gapfill → 全补集 → 两阶段）作为决策历史。
>
> **✅ 已实现（2026-08-06）**：`pgr align fill` + `pgr align rest` 落地
> （`cmd_pgr/align/{fill,rest,common}.rs`），复用 `pgr psl lift`（含
> `chr(+):` range 名修复）；集成测试 `tests/cli_align_fill.rs`（7 例）+
> `tests/cli_align_rest.rs`（4 例），全量测试通过；真实验证见
> `scripts/verify-align-fill-rest.sh`。

## 0. 命令定位

| 命令 | 语义 | query 侧 | 适用 |
|---|---|---|---|
| `pgr align fill` | 锚点间 gap 的 2D box 填充 | 双侧锚点推断区间（裁剪） | 快路径；锚点密集、共线好 |
| `pgr align rest` | 锚点未覆盖补集填充 | 每侧独立 trim+excise+holes；ref holes × 整套 query holes | 彻底路径；含首尾/无锚点区 |
| `pgr psl lift` | range 坐标回移 | — | fill/rest 共用 + 独立调试 |

两个 align 命令可串联：`fill` 输出 PSL 与 `rest` 输出 PSL 合并（cat）后一起喂
`pgr pl chainnet --syn`，重叠冗余交给链化去重。

## 1. `pgr align fill`

复用旧 `compute_boxes`（见下文历史 §3.1 的 2D 逻辑）：

1. 锚点按 `(t_name, q_name, strand)` 分组、按 t 排序；
2. 相邻同向非重叠锚点间，t/q gap 均落在 `[--min-gap, --max-gap]` 时构建
   box，双侧外扩 `--overlap`（默认 1 kb），clamp 到 contig；
3. box 提取 t/q 子序列 → LASTZ → LAV → PSL → `psl lift` 回移 → 合并。

参数：`--overlap`(1000) / `--min-gap`(100) / `--max-gap`(无限制) /
`--avail-psl` / `--preset` / `--query-depth`(50) / `--parallel`(8)。

## 2. `pgr align rest`

**两侧对称的一维区间管线**（ref 与 query 独立，不使用 2D 坐标）：

```
该侧 PSL 坐标（ref: t 侧；query: q 侧）
  → runlist JSON（区间并集）
  → span trim    -n TRIM       （锚点覆盖两端向内缩减）
  → span excise  -n MIN_ANCHOR （去掉较小的锚点块）
  → span holes                   （该侧全基因组补集）
  = 该侧 holes
```

与 `rept s-kmer` 的 `Fill → Excise → Fill`（`libs/pl/repeat.rs`
`run_repeat_runlist_pipeline`）同源，复用 `libs::runlist::span_op`
（Trim/Excise/Holes）。实现上直接调 `span_op` 库函数（或子进程
`pgr runlist span --op`，二选一，倾向库函数减少子进程开销）。

比对与坐标：

1. query holes 逐段用 `2bit range` 提取（header = `q_name:start-end`），
   合并为**一个多序列 FASTA**（lastz 支持多 query，LAV→PSL 已验证）；
2. 每个 ref hole 提取单序列，与 query holes 多序列文件跑 LASTZ
   （job = ref holes 数；query 扫描量 = query holes 总量，非整条基因组）；
3. LAV → PSL → `psl lift` 双侧回移（ref 侧 lift_target、query 侧 lift_query）→
   与锚点合并输出。

**不回退整条 query**：实测（全补集方案）hole × 整条 query 达 116 s，
不可接受；rest 的 query 裁剪 = query holes（锚点未覆盖部分）本身，
通常为基因组的 ~10% 量级。

参数：`--trim`(500) / `--min-anchor`(500) / `--max-gap`(无限制，可选过滤
超长 hole) / `--avail-psl` / `--preset` / `--query-depth`(50) / `--parallel`(8)。

## 3. `pgr psl lift`

**已存在，直接复用（2026-08-06 确认）**：`pgr psl lift` 已实现并验证通过，
无需新增。

**接口**：

```
pgr psl lift <infile> --q-sizes q.sizes --t-sizes t.sizes -o out.psl
```

* `--q-sizes` / `--t-sizes` 分开指定（`pgr fa size` / `pgr 2bit size` 同款
  `name<TAB>len`），可按需只 lift 一侧或双侧；
* 语义：`name:start-end` 形式的 qName/tName 回移到基因组坐标（内部即
  `Psl::lift_query`/`lift_target` 循环），**block qStarts/tStarts 一并回移**
  （`'-'` 链用 `size - end` RC 帧补码）；无 range 的记录原样保留；
* 越界（end > real size）告警并跳过该记录（`--strict` 可改为硬错误），
  不静默错位；
* 测试覆盖：`tests/cli_psl.rs` 4 个 lift 用例（q/t 侧、± 链 block 补码、
  `to-rg` 反查验证），stable 1.97 下全部通过。

## 4. 坐标回移专项（实现与测试依据）

lift 的边界情况（fill/rest 都会踩到，必须逐项测试）：

1. **多 contig**：sizes 表覆盖所有 query contig；qName 从 range 名解析出
   contig 名后查表；
2. **± 链**：qStart/qEnd 恒为正链坐标，仅 `'-'` 链的 block qStarts 处于
   反向互补帧，用 `size - end` 补码（`Psl::lift_query` 已实现，参考
   `coords.rs` 的 `reverse_range` 系列语义）；
3. **hole 起止**：`2bit range` 提取为 1-based 含端点（`chr:start-end`），
   PSL 为 0-based 半开；lift 的 `start_0 = start - 1` 换算必须逐一验证；
4. **提取坐标 ≠ lastz 实际比对坐标**：lastz 输出以提取序列为参照，
   PSL 坐标全部相对 range 起点；错位会静默产生错误覆盖统计——用
   "已知同源区间回移后坐标必须落在预期区间内"的断言测试兜底。
5. **带 strand 的 range 名**：`pgr rg` 生态产出 `chr(+):start-end`（如
   `NC_000913(+):85482-111492`）；`parse_subrange` 原本把 `(+)` 并入 contig
   名导致 sizes 查表失败、lift 静默跳过——2026-08-06 已修复（剥掉 `(+)/-)`
   后缀，与 PSL strand 字段冗余），并加单测与真实命令验证。

## 5. 参数汇总

| 参数 | fill | rest | psl lift |
|---|---|---|---|
| --overlap / --trim | 1000 | 500 | — |
| --min-gap / --min-anchor | 100 | 500 | — |
| --max-gap | 无限制 | 无限制 | — |
| --sizes | — | — | 必填 |
| --avail-psl / --preset / --query-depth / --parallel | ✓ | ✓ | — |

## 6. 验证计划

1. **lift**：已存在且验证通过（`pgr psl lift`，4 个测试覆盖 ±链 block
   补码/多 contig/越界）；fill/rest 编排直接子进程复用；
2. **fill 集成测试**：复用现有 `tests/cli_align_hybrid.rs` 的 gap/负链用例
   （已迁移为 `tests/cli_align_fill.rs`，7 例）；
3. **rest 集成测试**：多 contig query、trim/excise 阈值行为、ref holes ×
   query holes 套比对的覆盖（`tests/cli_align_rest.rs`，4 例，已通过）；
4. **真实验证**：`scripts/verify-hybrid-real.sh` 拆分适配为 fill/rest 两路，
   MG1655 × Sakai 对比覆盖/耗时（`scripts/verify-align-fill-rest.sh`，
   2026-08-06 release 实测，统一 `-p 8`；PSL 并集为中间产物口径，
   **MAF（`chainnet --syn`）为最终 syntenic 结果口径**）：

   | 引擎（-p 8） | 耗时 | PSL 记录 | MAF 块 | PSL 覆盖 | MAF 覆盖 |
   |---|---|---:|---:|---:|---:|
   | pgi-only | 1.3 s | 738 | 582 | 90.74% | 89.297% |
   | fill | 10.5 s | 3296 | 576 | 91.85% | **89.991%** |
   | rest（默认 syncmer 预筛） | 6.3 s | 1764 | 553 | 91.81% | **89.985%** |
   | rest（--sampler none 全量） | 23.7 s | 2005 | 557 | 92.20% | 89.997% |
   | fill + rest 合并 | ~16.8 s | 5060 | 605 | 92.10% | **89.866%** |
   | lastz-only | 133.2 s | 21793 | 382 | 93.11% | **90.237%** |

   结论：fill/rest 在 MAF 口径下几乎等价（89.99% vs 89.99%），比 pgi 高
   ~0.69 pp、比 lastz 低 ~0.25 pp，但快 13–21×（rest 6.3 s vs lastz
   133 s）；rest 默认预筛比全量（--sampler none）快 3.7× 且 MAF 口径
   无实际损失（89.985% vs 89.997%，差 0.012 pp ≈ 560 bp）。模拟灵敏度：
   rest 255/600 ≈ lastz 256/600（cell diff 1），假阳性 0.409% < 1%。

**重要发现**：fill+rest 合并（cat）在 MAF 口径下**低于 fill 或 rest 单独**
（89.87% vs 89.99%，差 0.13 pp）——两套 PSL 的重叠块（1 kb buffer 设计）
在 566–586 kb 等区域干扰 chainnet 链化，丢失 ~14.5 kb syntenic 覆盖。
PSL 并集口径的"合并更高"（92.10%）是中间产物假象。**建议直接使用
fill 或 rest 之一**（两者 MAF 覆盖几乎相同，89.99% vs 89.99%）；
如需合并，需先协调/去重重叠块。fill/rest 的差异只在非 `--syn` 场景
（rest 多补的退化同源在 syntenic 过滤中被滤除）。

### 6.1 并行与性能分析（2026-08-06 补充）

并行提取/转换优化（提取串行 → rayon 并行）后，MG1655 × Sakai 实测
（release，16C/32T）：

| 引擎 | -p 1 | -p 8 | -p 16 | -p 32 |
|---|---:|---:|---:|---:|
| fill（144 boxes） | 55.1 s | 10.4 s | — | 6.9 s |
| rest（全量模式，预筛前：197 holes × 整套 query holes） | 132.6 s | 22.1 s | 13.8 s | 13.5 s |

**fill + rest 串联**（-p 32，两路 PSL cat 合并）：20.4 s，覆盖 **92.38%**
（fill 91.85% ∪ rest 92.20%），MAF 608 块；lastz-only 134 s/93.11% 对照，
快 ~6.6×、差 0.73 pp。**注：该 92.38% 为 PSL 并集口径的旧数据（rest 当时
为全量）；以 MAF 口径复核后合并反而低于单路（见 §6 主表与重要发现），
不建议 fill+rest 串联。**

**并行加速观察**：

- lastz 计算部分 8 核接近理论 8×（rest 197 jobs × 单 job ~0.8 s / 8 ≈ 20 s）；
  墙钟达不到 8× 是串行开销（pgi 1 s + 2bit 转换 ~1 s + 提取/命令）所致；
- rest 的深层瓶颈是**每个 job 扫描整套 query holes**（~550 kb × 197 次，
  补集语义下 query 全量重复扫描是 lastz 计算量本质）；16→32 线程无增益
  （超线程 + IO 竞争），-p 16 即接近本机上限；
- 优化清单：提取并行（fill 2.35→0.37 s、rest 3.1→0.48 s）、LAV 转换并行、
  提取与 lastz 复用同一 rayon 池。

### 6.2 采样预筛配对（2026-08-06 补充，rest 默认路径）

**动机**：rest 每个 job 全量扫描 query holes（197 × 550 kb ≈ 108 Mb），
单线程与整条 lastz 相当（132.6 s vs 132.2 s）——"片段对片段"没有更快。
正解是 **k-mer 采样预筛配对**：ref hole 与 query hole 共享采样 k-mer 才配对，
只跑配对的小片段对。

**实测矩阵**（MG1655 × Sakai，-p 8，rest 单独）：

| 采样 / 参数 | 配对 job | 耗时 | 覆盖 |
|---|---:|---:|---:|
| none（全量） | 197 | 23.7 s | 92.20% |
| syncmer s8 w5 ms50 | ~17000 | 99.6 s | 92.17% |
| syncmer s15 w5 ms1 | ~1800 | 12.5 s | 92.00% |
| **syncmer s17 w5 ms1（默认）** | 757 | 6.4 s | 91.81% |
| syncmer s17 w3 ms1 | — | 7.0 s | 91.82% |
| minimizer k17 w5 ms2 | 349 | 4.0 s | 91.63% |
| syncmer s8 w5 ms100 | — | 68.8 s | 92.15% |

**决策**：默认 `syncmer s17 / window 5 / min-shared 1`——单线程 ~33 s
（vs 整条 132 s，**~4× 加速**，满足"片段对片段更快"）。

**覆盖口径（重要，2026-08-06 用户指出后复核）**：原始 PSL 并集口径下预筛
91.81% vs 全量 92.20%（-0.39 pp），但**经 `chainnet --syn` 共线性筛查后的
MAF 口径下两者几乎无差异**（89.985% vs 89.997%，差 0.012 pp ≈ 560 bp）：
那 18 kb 差异主要是 identity 58–72% 的退化同源（IS 元件、质粒拷贝），
在 syntenic 过滤中本就被滤掉。**syntenic 场景（fill/rest 的既定用途）
预筛无实际覆盖损失**；非 `--syn`（SD/重复）场景才需要
`--smer 15`（92.00% / 12.5 s）或 `--sampler none`（92.20% / 23.7 s）。

**PSL 并集差异来源**（262 个小碎片，总 18 kb）：多为高分歧退化同源
（配错/漏配），非 lastz 假阳性尾巴；`--unmatched full` 兜底无效（未匹配
ref holes 仅 11 个且无同源），top-K 限制反而降低覆盖（真同源不一定是
共享最多的 query hole）。

---

# 附录：原 hybrid 探索历程（历史，已拆分）

> 以下为 `pgr align hybrid` 从 gapfill 到两阶段的迭代记录，作为决策历史
> 保留；新代码按上文 fill/rest/lift 设计实现。

## 0. 背景与目标（原）

> 设计笔记（已实现，2026-08-06）。背景：FASTGA 论文（*FastGA: fast genome alignment*，
> Bioinformatics Advances 5(1):vbaf238，DOI 10.1093/bioadv/vbaf238）提出混合方案
> FastGA-gapfill——以 FastGA 比对为锚点，对每对连续同向锚点之间的区间跑 LASTZ
> 填 gap，最后合并两者结果。pgr 的本地对应物：`pgr align pgi`（原生 FastGA 风格）
> 快速找片段，`pgr align lastz` 精修，两套 PSL 合并喂 `pgr pl chainnet`。
>
> **2026-08-06 定稿变更（作者明确）**：hybrid 的补区间从"锚点间 gap"扩展为
> **target 全基因组补集**——pgi 锚点在 target 上未覆盖的全部区域（含锚点间
> gap、contig 首尾、无锚点 contig）都交给 LASTZ；query 侧用整条序列而非
> 锚点配对的 q 区间。gapfill 只是补集的一个子集。实现：`compute_holes` +
> hole × 全 query 的 LASTZ job。§5.2 实测覆盖追平 lastz（93.08% vs 93.11%），
> 但耗时接近 lastz 全跑（默认参数 116 s vs 135 s）——补集语义的固有代价。

> 设计笔记（已实现，2026-08-06）。背景：FASTGA 论文（*FastGA: fast genome alignment*，
> Bioinformatics Advances 5(1):vbaf238，DOI 10.1093/bioadv/vbaf238）提出混合方案
> FastGA-gapfill——以 FastGA 比对为锚点，对每对连续同向锚点之间的区间跑 LASTZ
> 填 gap，最后合并两者结果。pgr 的本地对应物：`pgr align pgi`（原生 FastGA 风格）
> 快速找片段，`pgr align lastz` 精修，两套 PSL 合并喂 `pgr pl chainnet`。
>
> **2026-08-06 定稿变更（作者明确）**：hybrid 的补区间从"锚点间 gap"扩展为
> **target 全基因组补集**——pgi 锚点在 target 上未覆盖的全部区域（含锚点间
> gap、contig 首尾、无锚点 contig）都交给 LASTZ；query 侧用整条序列而非
> 锚点配对的 q 区间。gapfill 只是补集的一个子集。实现：`compute_holes` +
> hole × 全 query 的 LASTZ job。§5.2 实测覆盖追平 lastz（93.08% vs 93.11%），
> 但耗时接近 lastz 全跑（默认参数 116 s vs 135 s）——补集语义的固有代价。
> 日期：2026-08-06。状态：已实现为 `pgr align hybrid`（算法
> `libs/align/hybrid.rs` + CLI/编排薄壳 `cmd_pgr/align/hybrid.rs`），集成测试
> `tests/cli_align_hybrid.rs`，文档 `docs/align-hybrid.md`。
> 关联：[[pgi-align.md]]、[[sd.md]]、[[references/fastga.md]]、[[paf-pangenome.md]]。
> 命令命名（2026-08-06 与用户讨论）：不做 `pgr pl align --hybrid`，直接做
> `pgr align hybrid`——`pgr pl` 是"暂时没想好该放到哪边"的命令的临时存放处，
> hybrid 放 `pgr align` 下方便用户发现。

## 1. 背景与目标

- FASTGA 论文里 FastGA 找的是"几乎横跨整段序列、最大 gap 约 40 bp"的局部比对，
  不自己做跨大 gap 的 chaining；论文用 LASTZ 填锚点之间的 gap，把 FastGA 的速度
  和 LASTZ 的灵敏度（论文实测最高）结合起来。
- pgr 侧目标：`pgr align pgi` 做全基因组快速粗比对 → 对 pgi 没覆盖到（或锚点之间）
  的区间跑 LASTZ 补比对 → 两套 PSL 合并后进 chain/net 流程，得到更完整的比对结果。
- 意义：pgi 边界一般比真实同源区短 1–11 bp（链的种子锚定边界，见 `docs/sd.md`），
  且对接近 SD 身份下限（~90–93%）的拷贝可能漏块；LASTZ 补的正是这类缺口。
  用户期望的"全基因组补集"比论文的 gapfill 更彻底：pgi 没覆盖到的所有
  target 区域（gap、首尾、无锚点 contig）都交给 LASTZ，gapfill 只是子集。

## 2. 现有资产盘点（全部已实现）

- `pgr align pgi`：FastGA 风格原生比对（syncmer 种子 → tube 链 → wave 扩展），
  输出 PSL。核心在 `src/libs/pgi/`（build/align/wave/mmap）+ `src/libs/ds/radix_sort.rs`
  （FastGA MSDsort 移植）。
- `pgr align lastz`：LASTZ 封装（`src/libs/lastz.rs`），输出 **LAV**，7 套预设
  set01..set07（set01 Human vs Chimp 最快，set07 Human vs Opossum 最远缘/灵敏），
  要求每个 FASTA 单序列，多 contig 需先 `pgr fa split name`。
- `pgr lav to-psl`：LAV → PSL 转换。
- `pgr pl chainnet`：native psl-chain-net-axt-maf 全链路，PSL 参数接受**文件或目录**
  （多个 PSL 直接合并链化）；`scripts/verify-pangenome.sh` 已验证
  FastGA PSL → chainnet 路径，pgi 的 PSL 与之格式一致。
- `pgr fa range` / `pgr fa split name`：按坐标提取区间、按记录切分序列。

## 3. 方案设计

### 3.1 流程

1. 粗比对：`pgr align pgi target.fa query.fa -o pgi.psl`
2. 确定补比对区间（**2026-08-06 定稿：target 全基因组补集**）：
   - 锚点按 target contig 合并覆盖（并集）→ 全基因组 holes（`[0, size)` 减并集）；
     锚点间 gap、contig 首尾、无锚点 contig 全部成为 holes（实现：
     `compute_holes`，内部区间运算，等价 `psl to-rg` → `rg cover` → holes）。
   - 每个 hole 外扩 `--overlap`（默认 1 kb）缓冲，clamp 到 contig；
     `--max-gap` 可选跳过超长 hole（默认不限制）。
3. 提取：每个 hole 用 `pgr 2bit range` 提取 target 子序列（单序列 FASTA）；
   query 侧用**整条 contig**（每个 hole × 每个 query contig 一个 LASTZ job）。
4. 精修：`pgr align lastz --preset <预设>`（预设由用户选择，见 §3.2/§3.6）→ LAV
5. 转换：`pgr lav to-psl` → 坐标回移（target 侧 hole 内坐标、query 侧整条
   序列坐标）→ `lastz.psl`
6. 合并（cat）：两套 PSL 直接并列输出，**不做去重**——重叠冗余交给 chainnet
   的链化处理（见 §3.7）
7. ChainNet：
   `pgr pl chainnet [--syn] target.fa query.fa psl_all/ -o out`

### 3.2 关键决策点

- **补比对区间**：2026-08-06 定稿为**全基因组补集**（holes），非论文的
  锚点间 gap——gapfill 是 holes 的子集。代价：LASTZ 处理量接近全基因组，
  耗时从 gapfill 的"接近 pgi"退化为"接近 lastz 全跑"（§5.2 实测）；
  收益：覆盖追平 lastz。真核场景建议配 `--max-gap` 过滤新序列区。
- **合并方式**：不做去重，两套 PSL cat 并列，交给 chainnet 链化时处理重叠冗余
  （§3.7；与论文 FastGA-gapfill 直接 cat 合并一致）。
- **`--syn`**：syntenic 共线性比对加；重复/SD 分析必须不加（`pgr sd align` 明确
  规定 chain/net 精修非 `--syn`，否则重排同源丢失）。
- **预设**：做成用户选项（复用 `pgr align lastz --preset set01..set07`）。
  泛基因组主场景差异小，默认贴近的 set01/set02；远缘比较由用户自选
  set06/set07——**不默认远缘预设**（2026-08-06 与用户讨论后调整）。
- **方向一致性**：pgi 与 lastz 的 target/query 顺序必须一致，PSL 的 tName/qName
  前缀统一，否则 chain 链向混乱。

### 3.3 边界处理策略（定稿，2026-08-06 与用户讨论后）

pgi 块边界一般比真实同源区短 1–11 bp。处理方式**不是缩短 pgi 的 PSL 记录**
（那会丢失已找到的片段，且需 t/q 同步变换的新函数），而是：

1. pgi 的 PSL **原样保留**，一个区间都不动；
2. 仅当计算"补集区间"（交给 LASTZ 的范围）时，对 pgi 的 target 侧区间做
   `trim(n)`（n ≈ 25–50 bp，大于 pgi 边界误差 1–11 bp 即可）——补集因此向外
   大出一个缓冲带，真实边界落进 LASTZ 的搜索范围；
3. LASTZ 跑补集，覆盖真实边界，与 pgi 完整块边界产生**少量有意重叠**；
4. 合并时 pgi 完整块 + lastz 块直接并列输出，不做去重（chain 分支的合理归并
   交由 chainnet 链化时处理，见 §3.7）。

实现上只需要一维区间运算：锚点并集（合并重叠）→ 全基因组 holes → 每 hole
外扩 `--overlap`（原 trim 缓冲语义并入 overlap），全部在 `compute_holes`
内完成（IntSpan 同源逻辑）。**不需要新增 PSL 记录变换函数**。

### 3.4 适用场景（定稿）：仅共线性搜索

Hybrid（pgi 锚点 + LASTZ 补 gap）模式**只适合共线性（syntenic）搜索**，
即 `pgr pl chainnet --syn` 一路。原因：

- 非共线性（SD/重复）场景下 net 不做 syntenic 过滤，重叠冗余块会一起保留，
  输出碎片化、覆盖重复计数——去重只能缓解，不是正解；
- SD 场景的正确解法是给 `sd search` 补 BISER 式 `MAX_EXTEND` 边界扩展
  （见 [[sd.md]] §4.8），而不是引入混合比对；
- 本方案不计划覆盖非共线性用途。

### 3.5 论文 FastGA-gapfill 参数对照（2026-08-06 通读论文后补充）

论文 §5.2 的 FastGA-gapfill 与本方案的对应关系（详见
[[references/fastga.md]] §12.3）：

- **补 gap 的对象**：每对"顺序一致、方向一致、不重叠"、间隔 ≤1 Mb（默认）的
  锚点 → 双侧 bounding box。与我们"锚点间 gap / 未覆盖区补集"一致，且同样
  隐含共线性前提（§3.4）。
- **重叠缓冲**：论文默认 box 与锚点重叠 **1 kb** 以利 LastZ 播种；ALNfill
  `alngap -e` 默认也是 1 kb。我们的 §3.3 目前用 trim 25–50 bp——比论文小两个
  数量级。**建议实测时对比 50 bp vs 500 bp vs 1 kb**：重叠越大 LastZ 播种越稳、
  但合并时冗余/重叠越多（ALNfill 已把 `-e` 造成的重叠列为已知问题）。
- **box 去嵌套**：论文只保留最小 bounding box（无包含关系）——对应我们的
  `runlist holes` 天然产出非重叠区间，无需额外处理。
- **合并方式**：论文直接把两套输出 cat 合并（PAF），未做去重；我们同样不做
  去重，两套 PSL 直接并列交给 chainnet 链化（§3.7）。重叠越大，chainnet 阶段
  冗余越多（ALNfill 已把 `-e` 造成的重叠列为已知问题），但这些都是下游链化
  的归并范畴，不在此处粗合并。

> 论文的 FastGA-gapfill 灵敏度接近 LastZ、速度比 LastZ 快 19.3×–137.5×——
> 这是本方案可行性的直接证据，也是后续验证的对照目标。
>
> **与论文的差异（2026-08-06 定稿）**：论文只补锚点间 gap（省时），pgr 补
> 全基因组 holes（彻底）。因此论文的"速度比 LastZ 快 19–137×"不适用于 pgr
> hybrid——补集语义下 LASTZ 处理量与全基因组相当，实测 E. coli 上 hybrid
> 116 s vs lastz 135 s（默认并行），仅省 ~14%，且随并行度提升优势增大
> （job 小而多）。省时形态是 gapfill 语义的属性，覆盖完整是补集语义的属性；
> 两者不可兼得，作者明确选择后者。

### 3.6 ALNfill 实现对照（2026-08-06 源码通读后补充）

ALNfill（`alnfill-main/`，Chenxi Zhou）是论文 FastGA-gapfill 的工程化实现，
源码分析见 [[references/alnfill.md]]。对本方案有直接影响的实现细节：

- **gap 过滤**：`alngap` 只补双侧 gap 都在 [100, 1M] 的区间（`-l`/`-m`），
  且用哨兵覆盖染色体首尾端。我们的 holes 方案应加同样过滤：
  小于 100 bp 的洞让 pgi 自己处理，大于 1 Mb 的洞跳过（与 §4"超长无锚点区间
  应跳过"一致）。
- **去冗余**：`alngap` 默认对 PAF 做"双侧被覆盖 ≤50% 的贪心过滤"
  （reciprocal best，`-a` 关闭）。pgi PSL 在重复区可能有多映射块，
  算 holes 前值得先按链上锚点去冗余。
- **方向**：`alngap` 不读 PAF strand 列，混合方向的锚点对也会成 box；
  论文描述是"一致顺序方向"。我们定稿为仅共线性（§3.4），实现时可加
  方向过滤或直接交给 `chainnet --syn` 收尾。
- **LastZ 选项**：ALNfill 只用 `--format=PAF:wfmash --ambiguous=iupac`
  （lastz 默认打分）；pgr 复用 `pgr align lastz --preset`，**预设由用户选择**——
  泛基因组差异小，默认贴近的 set01/set02，远缘比较才选 set06/set07，
  不默认远缘（2026-08-06 与用户讨论后调整）。
- **坐标回移**：ALNfill 提取区间、跑 lastz、把区间坐标回移成全长坐标、输出完整
  PAF 一步完成；pgr 若做成 `pgr align hybrid`，需要 `fa range` 提取时记录
  offset，lastz 输出回移后再并 PSL。
- **内存**：ALNfill 把两个基因组整库读进内存（sdict strdup），超大基因组不现实；
  pgr 的 2bit/loc 区间提取无此问题。

### 3.7 合并策略（定稿，2026-08-06 与用户讨论后）

早先计划在 `libs/align/hybrid.rs` 里做"重叠 >50% 保长去重"的粗合并。用户指出
pgr 的 chainnet 链化流程本身就处理重叠合并，这里粗合并既粗糙又多余，改为：

- **不做去重**：`run_hybrid` 把 pgi/锚点块与 LASTZ 补块**直接并列**输出
  （anchor 在前、lastz 在后），重叠冗余交给 `pgr pl chainnet` 链化时归并。
- 这与论文 FastGA-gapfill 直接 cat 两套 PAF 一致（§3.5）。
- 相应删除 `merge_dedup`/`overlap_count` 及其单元测试；集成测试改为断言
  "region 被覆盖"而非"精确块数"。
- 锚点来源参数定名 `--avail-psl`（"已有的 PSL"，不限于 pgi——FastGA、
  minimap2 等任意比对器输出均可直接喂入），复用一个已有 PSL 时跳过内部
  `align pgi`。

## 4. 已知的坑

- LASTZ 输入单序列限制 → hole 提取 + query 整条 contig 单序列化（2bit range）。
- `pgr align lastz` 只出 LAV，需 `pgr lav to-psl`。
- pgi 边界短 1–11 bp → lastz 补块与 pgi 块边界重叠 → 重叠交给 chainnet 归并
  （不做去重，见 §3.7）。
  决策见 §3.3：pgi PSL 原样保留，仅 holes 外扩 overlap，不缩 PSL 记录本身。
- 补集含物种特异插入/着丝粒等无同源序列，LASTZ 白跑——`--max-gap` 过滤
  超长 hole（真核推荐）；E. coli 上白跑量可忽略。
- `pgr align lastz` 默认 query-depth 50 是"先到先得"式截断，补 gap 场景若覆盖深
  可能丢块，必要时调大 `--query-depth`。
- 耗时随 holes × query contigs 的 job 数增长；并行度用 `-p`（默认 8）调，
  32 核机器可显著提速。

## 5. 验证方案

### 5.1 灵敏度评估（已执行，2026-08-06，`scripts/verify-hybrid-sensitivity.sh`）

口径借鉴论文 §5.1（[[references/fastga.md]] §12.1）：模拟 A、B 两个基因组
（各 6 Mb，由 10 kb 块组成；每块 = 目标区[长度 100–5000 bp × 分歧度
1–40%] + 随机填充；块序两基因组同序打乱，无跨块共线性；分歧按 80% 替换 +
10% 插入 + 10% 缺失引入）。每 (长度, 分歧度) 组合 20 重复（共 600 目标区）。
"恢复" = 目标区被比对覆盖 ≥95%（A、B 两侧都算）。结果为每格
`hybrid/pgi/lastz`（恢复数 /20）：

| L\d | 1% | 10% | 20% | 30% | 40% |
|-----|----|----|----|----|----|
| 100 | 2/1/2 | 1/0/1 | 1/0/1 | 0/0/0 | 0/0/0 |
| 200 | 3/1/3 | 4/1/4 | 6/0/6 | 2/0/2 | 1/0/1 |
| 500 | 6/6/6 | 6/6/6 | 6/3/6 | 6/0/6 | 1/0/1 |
| 1000 | 9/7/9 | 8/7/8 | 5/4/5 | 11/5/11 | 2/0/2 |
| 2000 | 15/15/15 | 15/15/15 | 14/14/16 | 15/13/16 | 13/1/14 |
| 5000 | 20/20/20 | 20/20/20 | 20/20/20 | 20/20/20 | 19/7/20 |

合计（/600）：**pgi 186 / hybrid 256 / lastz 256**（2026-08-06 补集方案定稿后
重跑：hybrid 完全追平 lastz）。假阳性碱基比例（A 侧落在目标区之外的比对碱基）：
pgi 0.061% / hybrid 0.391% / lastz 0.491%。

结论（与论文 "FastGA-gapfill 灵敏度接近 LastZ" 的结论形态一致）：

- **hybrid 灵敏度显著高于 pgi**（+65 目标区），gap-fill 主要补的是高分歧大目标区
  （2000 bp@40%: 1→13，5000 bp@40%: 7→19）——正是 §1 里 pgi 对 SD 身份下限
  附近漏块的场景。
- **hybrid 灵敏度 = lastz**（256 = 256，逐格差 0）——补集方案下 pgi 未覆盖
  区域全部交给 LASTZ，灵敏度自然追平最灵敏的 lastz。
- **三者假阳性都极低**（<1%）；hybrid 比 pgi 略高、与 lastz 相当——来自 lastz
  块边界超出 pgi 块的真实边界扩展（§3.3 的 buffer），是预期行为，非噪声 bridging。
  实证：全部假阳性碱基 100% 落在目标区边界 500 bp 内（无一条在随机填充区深处）；
  界外尾巴 pgi 中位数 2 bp/最大 15 bp，lastz 中位数 9 bp/最大 81 bp（X-drop 在
  随机序列上很快截断），hybrid 继承 lastz 尾巴。用论文 §12.1 的"按比对判定假阳性"
  口径（>95% 比对碱基在目标区外才算 false），这些尾巴不会把任何一条 lastz
  比对判假，真正的假阳性比对接近 0——碱基口径与论文口径不可直接比。

耗时（debug 构建，--parallel 8，6 Mb）：pgi-only 9.5s；hybrid 补 gap 本身
2.7s（复用 pgi 锚点，246 个 box），整链路含 pgi ≈ 12s；lastz 15s。补 gap 的
边际开销很小，真实数据里 lastz 开销随基因组规模放大更快（论文 19–137×），
hybrid 的省时优势会更强。

### 5.2 待补充验证

**✅ 已执行（2026-08-06 两轮，`scripts/verify-hybrid-real.sh`，release 构建）**。
MG1655 × Sakai（4.64 Mb × 3 replicon 5.53 Mb），三路均走
`pgr pl chainnet --syn`；覆盖 = PSL target 侧区间并集（`psl to-rg --target-coords`
→ `rg cover` → `runlist stat`），耗时 = 比对阶段。

**补集方案（2026-08-06 定稿后的最终实现）**：

| 指标 | pgi-only | hybrid | lastz-only |
|---|---:|---:|---:|
| 耗时（比对） | 1.19 s | 116.6 s | 135.0 s |
| PSL 记录数 | 738 | 10152 | 21793 |
| raw span（块区间和） | 5.02 Mb | 6.97 Mb | 8.59 Mb |
| MAF 块数（chainnet --syn） | 582 | 565 | 382 |
| target 并集覆盖 | **90.74%** | **93.08%** | **93.11%** |

结论（与论文 FastGA-gapfill 结论形态一致）：

- **hybrid 覆盖显著高于 pgi**（+2.34 pp），**几乎追平 lastz-only**
  （93.08% vs 93.11%，差 0.03 pp）——补集方案把 pgi 未覆盖的 target 区域
  （锚点间 gap、首尾、无锚点区）全部交给 lastz，覆盖目标达成。
- **hybrid 耗时接近 lastz-only**（116.6 s vs 135.0 s，省 ~14%，默认并行
  8 vs 4）——补集语义下 LASTZ 处理量接近全基因组，不再是 gapfill 的
  "接近 pgi"形态（作者明确接受此权衡，§3.5）。job 小而多（holes ×
  query contigs），高并行下优势更大（32 核可到 ~30 s 量级）。
- **无碎片化**：hybrid MAF 565 块 < pgi 582 块，lastz 补块被 chainnet 链化
  合并进既有链而非新增碎片；raw span/并集比 1.50（pgi 1.19 / lastz 1.99），
  重叠来自 holes 的 1 kb 外扩 buffer，符合设计。
- 记录数对比：lastz-only 21793 条（大量小片段）被链化吸收为 382 块 MAF，
  证明 chainnet 对噪声块归并有效；hybrid 10152 = 738 锚 + ~9400 补块
  （补集全部 hole × 3 query contig 的 lastz 输出），覆盖 93.08%。

**gapfill 旧方案（定稿前，仅供对照）**：覆盖 91.36%（+0.62 pp vs pgi）、
耗时 5.32 s（接近 pgi）——省时但不彻底；补集方案覆盖 +1.72 pp 而耗时
+111 s。作者确认采用补集方案。

**口径说明（借鉴论文 §12.4，[[references/fastga.md]]）**：覆盖统计只按比对
start/end 计（比对内 gap 也算覆盖）；本机 lastz 直接跑明文单序列 FASTA，
未喂 soft-mask（E. coli 无重复遮蔽需求，作为三路对照口径一致即可）。

## 6. 待办

- [x] 手动脚本跑通小规模验证（临时，不入库）。已验证：'+' / '-' 链的 box 计算、
      序列提取、lastz 补 gap、LAV→PSL 坐标回移、合并（cat）全链路。
- [x] 做成 `pgr align hybrid` 子命令：算法放 `src/libs/align/hybrid.rs`，
      编排（PipelineCtx + run_cmd）内联在薄壳 `cmd_pgr/align/hybrid.rs`（参考
      `pl/chainnet.rs`；编排不复杂、无共享部分，故不另立 `libs/pl/` 文件）。
- [x] 按 AGENTS.md 要求补 `cargo fmt` / `cargo clippy` 与测试：集成测试
      `tests/cli_align_hybrid.rs`（6 例）+ `hybrid.rs` 单元测试（7 例），
      全仓 1334 测试通过，fmt/clippy clean。
- [x] 灵敏度评估（fastga.md §12.1 口径）：`scripts/verify-hybrid-sensitivity.sh`
      已跑通并入库，结果见 §5.1、结论形态与论文 FastGA-gapfill 一致。
- [ ] （可选）真实数据对比：MG1655 vs Sakai 跑 pgi-only / hybrid / lastz-only，
      对比 chainnet 覆盖率与链完整性（§5.2 验证方案）。
