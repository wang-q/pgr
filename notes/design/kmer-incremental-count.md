# multik 计数复用 / 按 k 增量计数设计（pgr 任务 1.1）

> 目标：降低 `asm multik` 全流程 CPU（supermer ~19% + radix ~17%，合计
> ~36%）。对应 `~/Scripts/anchr/notes/design/pgr-tasks.md` §1.1。
> 本文是**设计 + 原型计划**，按清单要求先验证收益再决定投入，勿直接开写。

## 1. 目标与成功标准

- 硬约束：输出必须字节级一致（L1 smoke golden `f566e894…`）。
- lib 级：增量建表 vs 全量建表，wall 显著下降且输出逐字节一致。
- 端到端：收益兑现到 multik 全流程 wall（历史教训：supermer lib -32%
  未完全兑现到 anchr 全流程，不能只看 lib）。
- 任一不达标 → 关闭 1.1（挂起），记录数据，转投其他方向。

## 2. 现状梳理

### multik 建表模式（`anchr src/libs/asm/multik/schedule.rs`）

- k 序列（`auto_ks`）：短读 21/41/61/81，长读 31/61/91/121（步长 20-30）。
- 每轮全量重建 reads 表：base（k0 + 每个 later k）用
  `build_supermer[_qual]_slices`；probe（60）/repeat（`REPEAT_K`）/final
  （100）已跨轮共享（R+2 优化已完成，无重复建表）。
- unitig 表每轮 `build_direct_slices`（unitigs 每轮 append carried +
  过滤/prune，内容变化）。
- 峰值大头（MG1655 13.1 GB）：中段 k（81–128）reads 表建表瞬时 +
  10 条 chain 并发 unitig 表（每张 ~175–300 MB）。

### supermer 成本结构（占比需原型测量确认）

| 阶段 | 成本 | 跨 k 可复用？ |
| --- | --- | --- |
| stage-1 pack | N-free run 划分 O(L) + mval 滚动 O(L) + run 切分 O(窗口) + span 打包 O(span) | 部分（见 §4） |
| stage-1 sort | span records（依赖 collapse 效率） | 否（span 随 k 变） |
| stage-2 expand | O(有效窗口 × key_bytes)，随 k 线性增长 | 否（键集合全新） |
| stage-2 sort | k-mer 键 radix（~17% 全流程） | 否（跨 k 无键交集） |

## 3. 关键论据（已用测试验证）

### 3.1 字节级一致性的充分条件 = 键 multiset 一致

计数管线 = 发射窗口 multiset → radix 排序 → 分组累加，聚合顺序无关：
相同 multiset 必然得到相同表。因此增量方案**不需要模拟全量的计算顺序**，
只需证明"发射的窗口集合与全量完全一致"。

### 3.2 supermer 输出与 minimizer 长度 m 无关

m 只决定 stage-1 的 run 切分（collapse 效率）；每个有效窗口恰好被展开
一次，最终表与 m 无关。已加测试 `output_independent_of_minimizer`
（k=13/21/31/61 × 全部合法 m 输出一致）。

推论：增量方案可自由选择 m（如固定 12），不破坏字节级一致性；m 的选择
只影响性能（collapse 效率），可单独 A/B。

> 附带修复：发现 `m >= 17` 时 `pack_run` 的 u32 移位溢出（`2m` 位放不下）
> 会 panic——`build_table_slices_with_m` 的 m 校验收紧为
> `2..=min(16, k-1)`，超限返回错误而非 panic。

## 4. 增量方案候选

### 4a. minimizer 提取复用（metaMDBG 式，碱基空间版）

- 复用点（m 固定时）：N-free run 边界、mval/flp 序列。
- **存储不可行**：mval 每位置 u32（4×L）+ flip，G37/MG1655 级 reads 会
  给峰值内存加数 GB——MG1655 已 13.1 GB 贴近 L2 软上限，排除"预存全表"。
- 不存储变体：每轮重算 mval（O(L) 每位置 2 移位+掩码，本身轻量）——
  收益上限 = mval 滚动在 supermer 中的占比，预计很小。
- 结论预判：stage-2（expand+sort）是大头且不可增量，4a 的收益上限是
  supermer 的 pack 占比 × 可省比例，**大概率个位数 % 全流程**。

### 4b. unitig 表复用（unitigs 未变时）

- multik 每轮 unitigs 变化，但相邻轮可能大部分未变（需实测变化率）。
- 检测"内容未变"（长度+hash）成本 O(unitigs)；命中时省一张
  `build_direct_slices`（每张 175–300 MB × 10 chain 并发）。
- 4b 比 4a 简单（不碰 supermer 核心），优先级建议更高；先打点统计
  各轮 unitigs 变化率。

### 4c. radix 排序增量

跨 k 无键交集，排除；radix ~17% 全流程是"每轮必要的排序成本"。

## 5. 原型 A/B 设计

### 阶段 0：成本分布测量（pgr lib 级，先做）

- 输入：G37 reads 子集 / Lambda 20k；k=21/41/61/81。
- 手段：`build_impl` 临时阶段计时（stage-1 pack / stage-1 sort /
  stage-2 expand / stage-2 sort），或 `perf`。
- 产出：各阶段占比表 → 确定增量可省上限；顺带打点 multik 各轮
  unitigs 变化率（4b 可行性）。
- 决策门槛：若 stage-1 pack 可复用部分 < 20%（supermer 内），直接关闭
  4a，转向 4b 或内存方向。

### 阶段 0 实测（2026-08-24，`examples/supermer_stages.rs` + `PGR_SUPERMER_TIMING`）

**MG1655 20x 模拟 reads（142,836 条 / 21.4 M 碱基）**：

| k | pack | sort1 | expand | sort2 | total | spans | emitted |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 21 | 0.018s (13%) | 0.049s (35%) | 0.015s (11%) | 0.054s (39%) | 0.139s | 2.58M | 2.65M |
| 41 | 0.015s (10%) | 0.024s (17%) | 0.016s (11%) | 0.089s (62%) | 0.143s | 1.25M | 4.31M |
| 61 | 0.013s (7%) | 0.019s (10%) | 0.021s (11%) | 0.138s (72%) | 0.192s | 0.72M | 6.40M |
| 81 | 0.012s (4%) | 0.016s (6%) | 0.028s (10%) | 0.210s (78%) | 0.268s | 0.47M | 8.15M |

**Lambda 真实 reads（36,246 条 / 3.6 M 碱基）**：pack 占比 9–29%（数据小、
噪声大），sort1+sort2 合计仍是大头。

**结论（4a 关闭）**：唯一可跨 k 复用的 stage-1 pack 只占 supermer 的
6–30%（MG1655 下 k≥41 仅 6–10%），且其中 span 打包/run 切分依赖 k、
真正可省（mval 滚动 + N-free run 边界）更少——乐观上限（pack 全部可省）
仍低于设计稿 10% 门槛（k≥41 时 pack 6–10%）。**minimizer 提取复用
收益不达标，不做增量原型（阶段 1/2 取消）**。

**结论（4b 关闭，定性）**：multik 每轮 unitigs 都 append carried +
bridge_filter/recompact/prune，内容必然变化；"整表未变"才可复用，
命中率趋近于零；按 unitig 粒度增量维护表超出 1.1 范围。且 unitig 表
非 CPU 大头（见 `multik-complexity.md`），4b 收益有限。

**真正的成本大头**：stage-2 k-mer radix 排序（sort2，k≥61 时 70%+）。
`emitted`（唯一 span 展开的键数）随 k 增长到 8.15M，而最终表 `unique`
仅 ~990k——跨 span 共享 k-mer 的冗余 ~8× 真实存在（span 折叠只压"完全
相同的 span"，7 bp 错位的邻近 read 不折叠）。但该冗余**无低成本压缩
路径**：
- 分块排序去重：每块的排序总量不变（emitted 照排一遍），归并是额外
  成本，净收益为负；
- hash 聚合去重：内存不可接受（anchr 从 per-chunk HashMap 5.3 GB 峰值
  改到 radix 排序路径的历史教训）；
- 故 sort2 是"必要的排序成本"，radix 已线性最优。此结论推翻早先
  "分块去重压冗余"的设想。

**1.1 总结论（2026-08-24）**：增量计数（4a）与 unitig 表复用（4b）均按
数据/逻辑关闭；剩余 CPU 大头（radix 排序）无低成本优化路径。1.1 挂起，
等真实宏基因组/长读数据或新的架构思路再评估。

### 阶段 1：增量版原型（pgr lib 级）

- 4a 最小实现：固定 m=12，multik 循环外预处理 N-free run 边界（不存
  mval），每轮 build 复用 run 边界（省一次 O(L) 扫描 + run 划分逻辑）。
- A/B：k=31→41 增量 vs 全量，wall + 输出逐字节。
- 成功标准：增量 wall 显著下降（>15% lib 级）且逐字节一致。

### 阶段 2：端到端（anchr 配合）

- multik 全流程 G37/MG1655：wall/峰值/输出 golden 逐字节。
- 成功标准：端到端 wall 兑现（历史教训：不能只看 lib）。

### 决策门槛（任一不达 → 关闭 1.1）

- lib 级增量收益 < 10%；
- 端到端未兑现；
- 峰值内存恶化（固定 m 改变 collapse 效率的副作用）。

## 6. 风险与备选方向

- 增量实现复杂，golden 保护是底线；任何跳过窗口的优化都须先证明
  multiset 不变。
- 固定 m=12 对 k=21（默认 m=6）的 collapse 效率影响需 A/B 测。
- 若 4a/4b 均不达标：投入转向 **reads 表建表内存**（supermer 建表瞬时，
  MG1655 13.1 GB 峰值的大头），那是独立任务，收益量级更大（见
  `multik-complexity.md` 2026-08-23 晚节）。

## 7. 决策点

1. 阶段 0（成本分布测量）是否先做——~1 天，产出数据后 4a/4b 二选一或
   直接关闭。
2. 4a vs 4b 优先级：4b 更简单且不碰 supermer 核心，建议先打点验证。
3. A/B 不达标时按 §5 门槛关闭 1.1 并转投内存方向。
