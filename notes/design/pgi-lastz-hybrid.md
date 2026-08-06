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

# 附录：原 `pgr align hybrid` 探索历程（已归档）

2026-08-06 单命令 `pgr align hybrid` 拆分为 `fill` + `rest`，本附录为决策
历史摘要（详细迭代记录已由本文前半 + git 历史承载）：

1. **gapfill 版**（FastGA-gapfill 语义）：锚点间 2D box 补 gap，MG1655×Sakai
   覆盖 91.36%、5.3 s——快但不彻底（首尾/大 gap/无锚点区不补）；
2. **全补集版**：target holes × 整套 query，覆盖 93.08%（≈lastz 93.11%）
   但 116 s（hole × 整条 query 重复扫描）；
3. **两阶段版**（gapfill 2D + 动态补集）：覆盖 92.10% 但合并重叠块干扰
   chainnet 链化，174 s 反而更慢；
4. **定稿拆分**：`fill`（锚点间 gap，快）+ `rest`（一维补集 + syncmer 预筛
   配对，彻底且快），`psl lift` 独立复用；论文对照与灵敏度结论
   （rest 255/600 ≈ lastz 256/600）见本文前半 §6/§6.2。
