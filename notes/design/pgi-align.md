# pgi 两索引归并比对（设计定稿：v1 最小闭环）

> 定位：`.pgi` 的第一个比对消费者。输入两个已构建的 .pgi，
> 输出 PSL 块，喂给 `pgr pl chainnet`（UCSC 链化由 pgr 承担，见
> [[fastga.md]] §12.3 决策 3）。
> 状态：2026-08-02 定稿，随实现迭代。

## 1. 范围（v1 做什么、不做什么）

**做**：

1. 两个排序 .pgi 流的线性归并 → 种子命中（频率过滤）；
2. anti-diagonal 空间贪心链化 → 链（tube 语义简化）；
3. 每条链输出一个 PSL 块（链覆盖区间，无 CIGAR 细节）；
4. `pgr pgi align` CLI + 集成测试 + E. coli 验证。

**不做（v2 及以后）**：

- 局部扩展（banded/POA DP 细化 CIGAR，需要序列输入）——v1 链块直接由
  `pgr psl to_chain` / chainnet 接手（块结构即输入，UCSC 链化自算 score）；
- lcp 连续传播（adaptamer 变长种子）——固定 k 种子先行；
- ~~pbit 内嵌索引段消费~~（已按决策 A 放弃：索引不进 pbit，见 [[pbit.md]]）；
- mmap/流式读取（E. coli 规模直接整体载入内存）。

## 2. 数据流

```
ref.pgi + query.pgi（v1 无需序列文件：链块只依赖索引的 k-mer/位置/contig 表）
  │ 1. 排序流归并（频率过滤）→ 种子命中
  ▼
hit = (key, a_contig, a_pos, a_strand, b_contig, b_pos, b_strand)
  │ 2. 方向解析 + diag/anti 坐标变换
  ▼
fwd/rev 两个空间：(pos_a, pos_b'), diag = pos_a − pos_b', anti = pos_a + pos_b'
  │ 3. 按 (contig_a, contig_b, 方向) 贪心链化
  ▼
chain = 区间 + 对角线带（间距/带宽容忍，跨度过滤）
  │ 4. PSL 块输出（q=query, t=ref）
  ▼
out.psl → pgr psl to_chain → pgr pl chainnet（现有字节级验证的链化主场）
```

## 3. 关键语义

### 3.1 种子方向（pgi 双链条目）

`.pgi` 每条位置带 strand 标记（0=正向 k-mer，1=RC k-mer）。两索引归并时
key 相等按 (a_strand, b_strand) 解析方向：

| (a, b) strand | 含义 | 变换 |
|---|---|---|
| (0,0) / (1,1) | 正链命中（两侧实际窗口相等） | `pos_b' = pos_b` |
| (0,1) / (1,0) | 负链命中（b 窗口 = RC(a 窗口)） | `pos_b' = b_len − k − pos_b` |

负链 `pos_b'` 是 RC(b) 空间坐标；PSL 输出时用
`reverse_range(q_start, q_end, q_size)` 还原原始坐标、strand 记 `-`。

### 3.2 链化（简化 tube）

按 (contig_a, contig_b, 方向) 分组，命中按 (diag, pos_a) 排序后贪心：

- 当前链从首个命中开始；
- 延伸条件：`|diag − 链均线| ≤ band` 且 `Δpos_a ≤ max_gap` 且 `Δpos_b ≤ max_gap`；
- 否则收链：两侧轴向上种子跨度 `(last + k − first) ≥ min_span` 才保留。

参数默认对齐 FastGA（-f 10 / -c 85 / -s 1000，本实现为未加倍空间）：

| 参数 | 默认 | FastGA 对应 |
|---|---|---|
| `--freq` | 10 | `-f`（任一侧频率超限即跳过） |
| `--min-span` | 85 | `-c`（CHAIN_MIN） |
| `--max-gap` | 1000 | `-s`（CHAIN_BREAK/2） |
| `--band` | 128 | 对角线带宽（tube dgmin..dgmax） |

### 3.3 PSL 输出

- q = query 基因组，t = ref 基因组（`pgr pgi align <ref> <query>`）；
- 每条链 = 一个块：q 区间 `[min_pos_b'..max_pos_b'+k)`，t 区间 `[min_a..max_a+k)`；
- 负链块 q_start/q_end 走 RC 空间 → `reverse_range` 还原；
- match/mismatch 计数置 0（链化只消费块结构；真实计数待 v2 扩展）。

## 4. CLI

```
pgr pgi align <ref.pgi> <query.pgi> -o out.psl
  [--freq 10] [--min-span 85] [--max-gap 1000] [--band 128] [--merge-gap 5000]
  [--ref-seq ref.fa|2bit] [--query-seq query.fa|2bit]
```

- 参数默认对齐 FastGA（见 §3.2）；两侧索引参数必须一致（复用 `dist pgi` 校验）。
- `--merge-gap`：相邻共线性链的合并阈值（插入序列造成的对角线平移断链，
  见 §5.7）；
- `--ref-seq`/`--query-seq`：提供序列后进入 v3 分窗 banded 扩展（16 kb 窗口 +
  2 kb 重叠，巨型链也获得真实身份率，见 §5.4）。

## 5. 验证

1. 单元测试：归并方向解析（fwd/rev）、频率过滤、链化（band/gap/span 边界）；
2. 集成测试：小序列人工 PSL 形状断言；
3. E. coli：MG1655 自比对（`pgr pgi align`）全基因组覆盖 ~100%；
   MG1655 vs 近缘株与 `FastGA -psl` 结果做覆盖/共线性结构对照（非字节一致）。

## 5.1 实测结果（2026-08-02，release，默认参数）

**MG1655 自比对**（`mg1655.pgi` vs 自身）：

- 745 块；其中**正链主链 = 1 块覆盖 4,641,650 / 4,641,652 bp（99.9999%）**；
- 其余小块为真实重复结构（rRNA 操纵子、IS 元件等），负链 186 块对应反向重复
  （rRNA 操纵子的反向拷贝）。

**MG1655 vs Sakai（O157）**：

| 指标 | pgr pgi align | FastGA -psl |
|---|---|---|
| 块数 | 1019 | 701 |
| 正/负链 | 924 / 95 | 553 / 148 |
| 查询覆盖总和（span 求和） | 4.44 Mb（95.7%） | 4.63 Mb（99.7%） |
| 最大块 | 58 kb | 108 kb |

> 真实并集覆盖（2026-08-02 复核）：Sakai 75.8% vs FastGA 78.2%、Nissle
> 双方均 77.3%——span 求和的 95.7%/99.7% 有重叠重复计数，**真实差距仅
> ~2%**，且未覆盖区是株系特异序列。详见
> [[../benchmarks/bench-pgi-align-vs-fastga.md]]。

差异符合 v1 预期：固定 k=40 精确种子在 ~5% 分歧区段断链，FastGA 的
adaptamer（lcp 扩展）+ wave 对齐补平间隙。两者 PSL 均可被
`pgr psl to-chain` 正常消费（pgr 2038 行 / FastGA 12007 行 chain）。

## 5.2 v2 局部扩展实测（2026-08-02，banded SW + `--ref-seq/--query-seq`）

- 1019 块中 **1015 块被 banded 局部比对细化**（>30 kb 链回退为块，4 块）；
- 扩展身份率 **98.41%**（4.20M match / 68k mismatch），FastGA 为 **97.83%**
  （4.52M / 100k）——同区间参数下高度一致（差异来自块边界与覆盖范围）；
- 并行化：扩展链 rayon `par_iter`，1019 块 20.5s → 2.0s（32 核），
  输出与单线程逐字节一致。

**v2 边界与已知取舍**：

- ~~>30 kb 链不扩展~~（v3 已解决：16 kb 窗口 + 2 kb 重叠沿链对角线滑动，
  巨型链获得真实身份率，见 §5.4）；
- ~~banded 局部比对用线性 gap~~（v3.1 已改为**仿射 gap**：M/I/D 三状态，
  open -8 + extend -6；长 indel 表示为正确 gap 游程，块数 -25%，见 §5.8）；
- 块内多段同源时只取最佳局部段（FastGA 同语义）。

## 5.3 株系验证（2026-08-02，MG1655 vs 三株，v2 扩展）

| query | pgr 扩展身份率 | FastGA 身份率 | 结构 |
|---|---|---|---|
| nissle1917（Nissle，重排株） | 97.62% | 97.09% | 双方均为 ~30–70 kb 共线性碎片块（大规模倒位） |
| sakai（O157:H7） | 98.41% | 97.83% | 最大块 pgr 58 kb / FastGA 108 kb |
| ec958（UPEC） | 97.60% | — | 最大块 ~30 kb |

pgr 身份率稳定高于 FastGA ~0.5%（banded 局部比对取精确匹配核心；
FastGA 的 wave 延伸进分歧区、覆盖更广）。注意：此表为 v2 时代数据，
身份率只统计 ≤30 kb 的扩展块（>30 kb 链回退为无计数块）；v3 分窗扩展
后全部链均有真实身份率（见 §5.4）。

## 5.4 v3 分窗扩展实测（2026-08-02）

- **自比对**：主链 331 个窗口（≥10 kb 块）身份率**精确 1.0000000**
  （5.30M match / 0 mismatch）；整体 99.93%（其余为 rRNA/IS/反向重复
  等真实重复结构的微小差异）；
- **MG1655 vs Sakai**：全部 1093 块扩展，身份率 **98.42%**（4.49M/72k），
  FastGA 97.83%（4.52M/100k）——匹配碱基量从 v2 的 4.20M 提升到 4.49M，
  大链不再缺失；
- 运行时间 2.0s（并行）。

## 5.5 10 株 cohort 两两验证（2026-08-02，45 对，默认参数含 --merge-gap 5000）

扩展块身份率矩阵（行×列 = ref×query；块数为扩展块数）：

| pair | 块数 | 身份率 |
|---|---|---|
| mg1655–sakai | 862 | 0.9834 |
| mg1655–nissle1917 | 1378 | 0.9746 |
| mg1655–cft073 | 884 | 0.9739 |
| mg1655–e2348_69 | 880 | 0.9740 |
| mg1655–e24377a | 957 | 0.9859 |
| mg1655–ec042 | 934 | 0.9778 |
| mg1655–ec2011c_3493 | 903 | 0.9863 |
| mg1655–ec958 | 958 | 0.9747 |
| mg1655–se11 | 853 | 0.9870 |
| sakai–nissle1917 | 1592 | 0.9736 |
| sakai–cft073 | 1230 | 0.9700 |
| sakai–e2348_69 | 1138 | 0.9718 |
| sakai–e24377a | 1227 | 0.9807 |
| sakai–ec042 | 1240 | 0.9724 |
| sakai–ec2011c_3493 | 1261 | 0.9783 |
| sakai–ec958 | 1122 | 0.9712 |
| sakai–se11 | 953 | 0.9772 |
| nissle1917–cft073 | 1835 | 0.9962 |
| nissle1917–e2348_69 | 1291 | 0.9870 |
| nissle1917–e24377a | 1547 | 0.9736 |
| nissle1917–ec042 | 1621 | 0.9734 |
| nissle1917–ec2011c_3493 | 1608 | 0.9721 |
| nissle1917–ec958 | 1596 | 0.9880 |
| nissle1917–se11 | 1170 | 0.9735 |
| cft073–e2348_69 | 1040 | 0.9864 |
| cft073–e24377a | 1022 | 0.9733 |
| cft073–ec042 | 1121 | 0.9720 |
| cft073–ec2011c_3493 | 1102 | 0.9714 |
| cft073–ec958 | 1151 | 0.9867 |
| cft073–se11 | 816 | 0.9730 |
| e2348_69–e24377a | 981 | 0.9731 |
| e2348_69–ec042 | 1070 | 0.9727 |
| e2348_69–ec2011c_3493 | 1021 | 0.9714 |
| e2348_69–ec958 | 882 | 0.9847 |
| e2348_69–se11 | 868 | 0.9716 |
| e24377a–ec042 | 1200 | 0.9765 |
| e24377a–ec2011c_3493 | 1142 | 0.9911 |
| e24377a–ec958 | 1115 | 0.9726 |
| e24377a–se11 | 1040 | 0.9913 |
| ec042–ec2011c_3493 | 1201 | 0.9761 |
| ec042–ec958 | 1131 | 0.9728 |
| ec042–se11 | 877 | 0.9756 |
| ec2011c_3493–ec958 | 1211 | 0.9734 |
| ec2011c_3493–se11 | 829 | 0.9908 |
| ec958–se11 | 942 | 0.9727 |

全部 45 对 ~60s 完成（并行扩展 + 链合并）。分布 97.0–99.6%，与 E. coli 株系亲缘
关系一致（e24377a/se11/ec2011c_3493 聚类 99.1%+、nissle–cft073 99.6%）。
合并块后身份率比无合并低 ~0.2–0.5%（合并块把分歧间隙区域计入计数，更接近真实）。
注意：身份率基于种子链发现的块（偏保守区段），且块身份计数整体比 FastGA
高 ~0.5%（banded 局部取精确核心）。

## 5.6 性能（2026-08-02，与 FastGA 端到端持平）

`pgr pgi align` 全流程（索引 ×2 + 扩展比对）1.32s vs FastGA 单命令 1.22s
（**1.08×**，MG1655 vs Sakai）。优化两个关键点：

1. **banded DP 内层按带限列迭代**：只扫 `|j−i−diag0|≤band` 的 65 列而非
   全部 16000 列（原实现 246× 浪费）；
2. **窗口级负载均衡**：所有链的所有窗口摊平进同一 rayon 流——自比对主链
   332 窗口不再成为单线程长尾（37.8s → 0.84s，45×）；跨株 2.0s → 0.66s。

修复前后输出逐字节一致。详见
[[../benchmarks/bench-pgi-align-vs-fastga.md]]。

## 5.7 链合并（--merge-gap，2026-08-02）

插入序列（IS 元件等）使局部对角线平移超出 band，把同一共线性块切成多条
贪心链。新增链后合并：同 (contig_a, contig_b, strand) 组、间距 ≤
`--merge-gap`（默认 5000）、对角线差 ≤ band 的相邻链合并为一个块。

| 输入对 | merge-gap 0 | merge-gap 5000 |
|---|---:|---:|
| MG1655 vs Sakai | 1019 块 / 最大 58 kb / 覆盖 4.44 Mb | 718 块 / 最大 83 kb / 4.51 Mb |
| MG1655 vs Nissle | 1634 块 / 最大 30 kb / 4.46 Mb | 1259 块 / 最大 64 kb / 4.54 Mb |

块数 -23~30%、最大块 +40~100%、覆盖 +1.6~1.8%。剩余断链主要来自真实
重排（对角线差超 band，不合并）与分歧区（>5 kb 间隙），后者需
adaptamer/lcp 变长种子（未来工作）。

## 5.8 仿射 gap 扩展（2026-08-02，v3.1）

banded 局部比对从线性 gap（每次 -8）升级为**仿射 gap**（M/I/D 三状态：
open -8 + extend -6，`AlignmentParams` 默认参数），长插入（IS 元件等）
被表示为正确的 gap 游程而非碎成多个块。

| 输入对 | 线性 gap 块数 | 仿射 gap 块数 | 缩减 |
|---|---:|---:|---:|
| MG1655 vs Sakai | 4710 | **3512** | -25% |
| MG1655 vs Nissle | 7451 | **5609** | -25% |

身份率不变（0.9834/0.9747——对齐内容等价，仅 indel 结构更干净）。
匹配碱基 +225、错配 -44（Sakai），CIGAR 质量明显改善。

## 5.9 adaptamer 部分种子：负结果（2026-08-02）

按 FastGA 的 lcp 归并机制（plen ≥ 12 的部分匹配种子，`--min-shared`）实现
并实测，**所有阈值均劣于精确匹配**（MG1655 vs Sakai，默认参数）：

| min_shared | records | blocks | identity |
|---|---:|---:|---:|
| 40（精确，默认） | 862 | 3512 | 0.9834 |
| 30 | 1130 | 4657 | 0.9827 |
| 25 | 1302 | 6072 | 0.9810 |
| 20 | 1614 | 8680 | 0.9781 |
| 12（FastGA plen 下限） | 3375 | 53015 | 0.9496 |

**原因**：部分种子（共享 12-39 碱基）在分歧区既补充同源种子，也带来大量
假阳性弱种子；我们的贪心链化 + banded 扩展没有 FastGA tube/wave 那样的
种子质量区分机制，弱种子生成大量低质量新链并拉低身份率。

**结论**：adaptamer 的收益依赖 FastGA 的链化/扩展机制，不能直接移植到当前
管线。保留 `--min-shared` 作为实验开关（默认 = k 精确匹配），后续若引入
FastGA 式 tube 链化再重新评估。

## 5.10 端到端管线验证（2026-08-02）

`pgr pgi align` → `pgr psl to-chain`（chainnet 内部）→ `pgr pl chainnet
--syn` 全链路，与 FastGA 驱动版本（`FastGA -psl` → `pgr psl swap` 对齐
q/t 角色 → 同一 chainnet）对比 syntenic MAF：

| 输入对 | 指标 | pgr 管线 | FastGA 管线 |
|---|---|---:|---:|
| MG1655 vs Sakai | syntenic 覆盖 | **87.7%**（392 块） | 89.3%（506 块） |
| MG1655 vs Nissle | syntenic 覆盖 | **82.9%**（541 块） | 85.3%（711 块） |

结论：
- 管线端到端可用（0.4s 产出 syntenic MAF），块结构**比 FastGA 更平滑**
  （块更少：392 vs 506 / 541 vs 711）；
- 覆盖差 1.6-2.4%，来源是分歧区：FastGA 的 wave aligner 能桥接 banded
  窗口跳过的低分区间（同 §5.4 观察）；
- 角色约定：`pgr pgi align <ref> <query>` 的 PSL 是 q=query/t=ref；
  FastGA 的 PSL 是 q=source1/t=source2，两者互换，喂 chainnet 前需
  `pgr psl swap`（或调整参数顺序）。

## 5.11 Myers wavefront 扩展器：移植与负结果（2026-08-02）

按 FastGA 技术路径移植了 `align.c forward_wave` 的 Myers wavefront 核心
（V[k] 波前 + 三分支更新 + snake + 逐波前驱精确回溯，锚点双向扩展），
见 `src/libs/alignment/wave.rs`（单元测试通过：全等/单错配/插入路径）。
作为独立实现保留，**不接入当前管线**。

接入测试（替换 banded 作为窗口扩展引擎，两种锚点策略）：

| 引擎 | PSL 块数 | chainnet syntenic 覆盖 |
|---|---:|---:|
| banded（当前） | 3512 | **87.7%** |
| wave（窗口中心锚点） | 6477 | 71.6% |
| wave（最近匹配锚点） | 12184 | 32.7% |

**负结果**：unit-cost 波前对 indel 无 gap 结构偏好，CIGAR 碎片化（块数
3-4×）；且波前从锚点贪心延伸，产出大量无法通过 syntenic 过滤的低质量块。
**根因**：FastGA 的 wave 依赖其 tube 链化（`align_contigs`）提供的锚定
上下文与阈值（ALIGN_MIN/ALIGN_RATE）；脱离该上下文单独移植扩展器不成立。

**下一步**：移植 tube 链化（种子流 → anti-diagonal 桶 → tube → 每 tube
wave 扩展），届时 wave 引擎按 FastGA 语义接入；在此之前 banded 保持默认。

## 5.12 tube 链化移植

按 FastGA `align_contigs` 移植了 **tube 链化**（`chain_tubes`）：
种子按对角线分桶（宽 64）→ 相邻桶对按 a 位置归并 → 维护 tube
（anti 覆盖 = 种子区间并集、对角线范围），`CHAIN_BREAK`（1000 bp）断开、
`CHAIN_MIN`（85 bp）触发。单元测试通过（共线合并/断链/覆盖过滤）。

## 5.13 FastGA `Local_Alignment` 移植完成（2026-08-02）

按 FastGA `align.c` 完整移植了 mid-line wave 局部比对器：

1. **`forward_wave_mid`**（`src/libs/alignment/wave.rs`）：0-wave 在 tube 对角
   带 [dgmin..dgmax] 的每个对角线上从 mid-line（anti=amid）起 snake，波前每
   波扩展 ±1 对角线；保留 `PATH_LEN=60` 位匹配向量 + `PATH_AVE=42` 门控
   （trim point = 最后一个窗口匹配数达标的 best），`TRIM_MLAG=250` 终止、
   `WAVE_LAG=70` 剪枝，内存与波带宽度成正比（不存全历史）；
2. **Myers O(ND) 分治回溯**（`split_nd` + `dandc_nd`，FastGA 的
   substitution=1 度量）：两个 wave 端点之间的 span 由 D&C 精确重建编辑脚本
   （D/I 操作，替换为隐式），`ops_to_columns` 还原 CIGAR 列；
3. **`local_alignment`**：正向 wave（向上）+ 镜像反向 wave（向下）+
   `DUB_TRIM` 短扩展重试，输出 q/t 对齐列；
4. **tube 循环**（`extend_tube`）：`BUCK_ANTI=128` 滑动 mid-line，`alow` 推进
   到 forward 端点（eant），`alast` 按（contig 对、strand、对角线桶）分组去重
   （FastGA 的 `alast` 每对相邻桶重置，不能跨组携带）。

**MG1655 vs Sakai 对照**（`--workflow tube`，默认 greedy 未变）：

| 方案 | PSL 块 | chainnet syntenic 覆盖 | 耗时 | 峰值内存 |
|---|---:|---:|---:|---:|
| 贪心链 + banded（默认） | 862 | 87.7%（392 块） | 1.36s | 1.38 GB |
| **tube + Myers wave（新）** | 643 | **88.2%**（517 块） | 8.7s | 425 MB |
| FastGA 管线 | 701 | 89.3%（506 块） | ~0.7s | ~0 MB |

结论：质量已超过贪心基线并逼近 FastGA（覆盖差 1.1%）；内存比 banded 低 3×。
**速度仍是短板**（8.7s vs FastGA 0.7s）：单次 `Local_Alignment` 成本 ~0.18ms、
调用数 ~4 万（FastGA 写入 967 条比对，失败调用也远少于我们），差距来自：

- FastGA 的 CIGAR 由 wave 自身的 Pebble 稀疏 trace 产生（`dandc_nd` /
  `Compute_Alignment` 在 FastGA 源码中是**死代码**，从未被调用）；我们每次
  调用都跑完整 D&C + 列重建；
- FastGA 失败调用的 trim 终止更快、波带内单元更简单（C 内联数组）。

## 5.14 性能对齐：调用数与 tube 结构（2026-08-02）

给 FastGA 源码加计数器重新编译后实测（MG1655 vs Sakai）：

**FastGA 总共只有 1062 次 `Local_Alignment` 调用**（含失败），每 tube 平均
1.4 次；其 tube 平均 7.7 kb、最大 119 kb。我们 8.7s 的根因是 **~4 万次调用**
（tube 平均 30 kb、最大 2.4 Mb——精确 40-mer 种子太密，把 ~45-70% 身份的
分歧岛桥接成巨型 tube，mid-line 以 ~90-230 bp 步长在里面滑动）。

已应用的优化（质量不变：Sakai 88.2% / Nissle 84.0%）：

1. **tube 并行化**：tube 间无依赖（`alast` 重叠跳过被并行 + 输出端
   `dedupe_contained`（双轴 ≥80% 重叠去重）替代）；
2. **反向/互补序列预计算**：`rt`/`rq` 从每 tube 分配改为每 contig 一次；
3. **tube anti 上限 40 kb**：巨型 tube 切片成并行任务（负载均衡）。

最终：8.7s → **1.7-1.9s**（8 线程），峰值内存 ~0.8 GB。

**进一步实验（2026-08-02，均记录避免重试）**：

- 调用统计：58,370 次调用中 **386 个零块 tube 消耗 53,333 次（91%）**——失败
  调用都在 ~45-70% 身份的分歧岛内（wave trim 冻结后空转 ~70-99 波）；
- **trim 冻结早退**（保留）：`last_good` 连续 60 波不更新即提前终止——失败
  调用的端点其实在 wave ~10-20 波就冻结了，早退只省空转、质量不变
  （88.2% 保持，1.85s→1.73s）；
- `CHAIN_BREAK` 调小（300/100 bp）：更慢（tube 碎片化、重叠调用更多）且掉
  质量（87.3%）——不是正确的杠杆；
- 中心对角线滑窗身份率门控：太激进（覆盖 88.2%→70.9%，薄保守区/偏移对角线
  被误杀）；
- 种子覆盖密度门控（cov/span）：零块与生产性 tube 分布重叠，无法干净区分；
- 种子邻近门控（amid ±300bp 内无种子则跳过）：无效——失败调用的 amid 本来
  就在种子附近（分歧岛内也有稀疏 exact-40 hit）。

结论：剩余差距（~1.7s vs FastGA 0.7s）来自失败调用数量，而失败调用根植于
种子结构（我们的 exact-40 hit 桥接分歧岛、tube 平均 19.8kb vs FastGA 7.4kb；
FastGA 的链在岛边缘断开，其链形成种子同为 363 万条密集分布，断链机制与其
adaptamer 种子选择/对角桶处理相关，尚未完全复刻）。

## 5.15 大 tube 同源门控（2026-08-02，保留）

前一轮尝试的各种门控失败后，最终找到可用的版本：**对 span > 10 kb 的大 tube
做多对角线滑窗身份检查**：

- 9 条对角线横跨 tube 带（[dgmin..dgmax] 均匀采样），每条对角线上沿 anti
  轴滑 64 bp 窗口，取最大窗口身份率；
- 所有对角线、所有窗口的最大值 < 50%（64 bp 窗口内匹配 < 32）→ 跳过整个
  tube（只产生被拒绝的调用）；
- 只查大 tube（小 tube 便宜、薄保守区由 wave 兜底），成本 ~0.05s。

为什么之前的门控不行、这个行：128 bp 窗口会混入分歧侧翼把身份率拉低（误杀
薄块）；只查中心对角线会漏掉偏移对角线上的同源；多对角线 + 64 bp 窗口 +
50% 阈值三者组合后误杀率≈0（Sakai/Nissle 覆盖与无门控完全一致）。

最终效果（8 线程）：

| 对 | chainnet 覆盖 | 耗时 | 峰值内存 |
|---|---:|---:|---:|
| MG1655 vs Sakai | 88.2%（512 块） | **1.52s** | ~0.8 GB |
| MG1655 vs Nissle | 83.9%（741 块） | **1.55s** | 0.74 GB |

与上轮（无门控 1.73-1.85s）比再快 ~15%，质量不变；累计相对最初 8.7s 为
**5.7×**。距 FastGA（0.7s / 89.3%）仍有 ~2.2× 速度差与 ~1.1% 覆盖差。

## 5.16 tube 合并顺序 bug 修复（2026-08-02，质量 +0.7%）

追查 1.1% 覆盖差时发现 `chain_tubes` 的一个真实 bug：桶内种子按
**(diag, a_pos)** 排序、合并按 **a_pos**——当对角线在桶内漂移时（真实案例：
diag -4492 的种子在 a_pos 112701，diag -4490 的种子在 a_pos 111544+），
anti 顺序与 a_pos 顺序相反，导致：

- 高 anti 的种子先被处理并抬高 `ahgh`，低 anti 的稠密种子随后全部落入
  `anti < ahgh` 分支，`cov` 不累计（实测 cov=42 < CHAIN_MIN=85）→ **整管被
  丢弃**；
- 缺失区域（t 111544-121197，98.4% 身份、7762 个 hit）因此完全无 tube。

修复：排序键改为 **(diag 桶, anti)**，合并按 **anti**（FastGA 的
`ipost`-ordered merge，`if (apost < ipost)` 取小者）。新增回归测试
（`tube_merge_uses_anti_order_when_diagonal_drifts`）。

效果（8 线程）：

| 对 | 修复前 | 修复后 | FastGA |
|---|---:|---:|---:|
| MG1655 vs Sakai | 88.2% / 512 块 / 1.52s | **88.9%** / 517 块 / 1.52s | 89.3% |
| MG1655 vs Nissle | 83.9% / 741 块 / 1.55s | **84.4%** / 743 块 / 1.43s | 85.3% |

与 FastGA 的覆盖差从 ~1.1-1.4% 缩小到 **0.4-0.9%**；缺失区域从 55 kb 降到
25 kb。速度不变（tube 结构修正后并行负载更均衡）。

## 5.17 dedupe 误删延伸块（2026-08-02，覆盖 +3.2 kb）

继续追剩余缺失区发现第二个 bug：`dedupe_contained` 的 80% 双轴重叠阈值
会删除**延伸块**——同一 tube 连续两次调用产生的块重叠 ~87%（第二次反向
波延伸更远），若第二次的延伸恰好是真实覆盖（实测：t=4094559..4118787
比 t=4094464..4115709 多覆盖 3.1 kb 的 90.9% 身份区），会被误判为重复删掉。

修复：重叠阈值 0.8 → **0.95**（只删相邻桶 tube 的近似完全重复块，保留有
真实延伸的块）。新增回归测试 `dedupe_keeps_blocks_that_extend_earlier_ones`。

效果：Sakai 覆盖 4,124,601 → **4,127,886 bp**（88.9% 不变，缺失 24.9→21.7
kb）；Nissle 84.4% 保持；耗时不变（~1.5s）。

剩余缺失（21.7 kb / 324 个小区域）经查均为 **indel 复杂区**：精确 40-mer 在
固定对角线上 0 或极少共享（如 t 2451275-2455928 与 FastGA 块同偏移身份仅
28%、共享 40-mer 为 0），FastGA 靠 adaptamer 部分种子 + wave 跨越内部漂移。
这属于剩余的 adaptamer 种子选择工作。

## 5.18 调用数骤降与排序优化（2026-08-02，总耗时 1.5s→1.0s）

anti 序合并修复的连锁效应：**tube 结构与 FastGA 对齐**（748 个 tube、平均
15.4 kb vs FastGA 815 条/14.8 kb），调用数从 58,370 骤降到 **883**（FastGA
1062），每 tube 平均 1.03 次调用；顺序耗时 11.4s → 2.7s。

剩余瓶颈变为流水线开销，阶段实测（8 线程）：

| 阶段 | 耗时 |
|---|---:|
| merge_seed_hits（顺序） | 139 ms |
| chain_tubes 排序 | 275 ms（原 741 ms） |
| extend（wave 调用） | 414 ms |
| 加载（pgi×2 + fasta×2）+ PSL 写 | ~0.4 s |

排序优化：sort key 打包成单个 u128（contig/strand/对角桶+偏移/anti），并改用
rayon `par_sort_unstable_by_key`（741→275 ms）。

**最终对照（8 线程）**：

| 对 | chainnet 覆盖 | 耗时 | 峰值内存 |
|---|---:|---:|---:|
| MG1655 vs Sakai | 88.9%（520 块） | **0.97s** | 569 MB |
| MG1655 vs Nissle | 84.4%（744 块） | **1.17s** | 783 MB |
| FastGA | 89.3% / 85.3% | 0.7s | ~0 MB |

距 FastGA：速度 1.4-1.7×，覆盖差 0.4-0.9%。剩余瓶颈：extend 的 wave 每调用
成本（~0.47ms，Rust vs C 单元开销）、merge 顺序扫描、加载开销；indel 复杂区
覆盖需 adaptamer 部分种子。

## 5.19 pgi 读取批量解析（2026-08-02）

阶段实测发现 **pgi 索引加载（114 MB × 2 个文件）占 ~0.5s，是全流程最大单项**。
`PgiIndex::read` 原来逐记录 `read_exact`（380 万条 × trait 对象虚拟分发），
改为 **1 MB 分块批量读取 + 切片解析**（块大小按 `rec_size` 对齐，避免错位）。
加载从 ~0.7s 降到 ~0.5s（greedy 无序列路径实测 0.70→0.52s）。

最终全流程 ~1.0-1.2s（8 线程，Sakai），903 测试全过（含 pgi 读写测试）。

## 5.20 tube workflow 默认部分种子（2026-08-02，质量 89.1%/84.7%）

剩余 21.7 kb 覆盖差是 indel 复杂区（精确 40-mer 在固定对角线 0 共享）。anti
修复后重测 `--min-shared`：

- `--min-shared 20`（k/2）：Sakai **89.1%**（+0.2%）、Nissle **84.7%**（+0.3%），
  耗时 ~1.2-1.4s；freq 过滤把 hit 增量限制在 +20%（393.9 万 vs 326.8 万）；
- `--min-shared 30`：反而更差（88.8%）——部分匹配噪声未被抑制；
- 结论：FastGA 默认就是部分种子（adaptamer 共享 ~31/40 bp），把 **tube
  workflow 的默认 min-shared 改为 k/2**（greedy 保持 exact 默认，862 块不变）。

最终（8 线程，tube 默认参数）：

| 对 | chainnet 覆盖 | 耗时 |
|---|---:|---:|
| MG1655 vs Sakai | **89.1%**（586 块） | 1.24s |
| MG1655 vs Nissle | **84.7%**（793 块） | 1.22s |
| FastGA | 89.3% / 85.3% | 0.7s |

覆盖差缩至 **0.2-0.6%**；剩余缺失 ~15 kb 需 adaptamer 最小种子**选择**（抑制
部分匹配噪声，FastGA 的 `is_minimal` 语义）而非盲目的部分匹配。

## 5.21 merge 并行化（2026-08-02）

`merge_seed_hits` 的每条目前缀查询彼此独立，按 4096 条分块用 rayon 并行合并
（输出顺序自由——链化都会重排）：139 ms → **61 ms**。pgi 解析并行化试验无益
（瓶颈是 114 MB 的磁盘读取而非解析，已回退）。

最终（8 线程，tube 默认参数）：Sakai 89.1%/586 块/~1.1s，Nissle
84.7%/793 块/~1.2s；903 测试全过。剩余瓶颈：pgi 加载（磁盘 I/O，~0.4s）、
chain 排序、extend 的 wave 每调用成本、15 kb indel 复杂区。

## 5.22 移除大 tube 同源门控（2026-08-02，Sakai 与 FastGA 持平 89.3%）

追查最大缺失区（t 367590-374927，7.3 kb，99% 身份、4160 个共享 40-mer）发现
是**门控误杀**：门控的对角线采样步长 ~13 bp，真实比对对角线 -58601 落在采样
点之间——采样对角线与真实差 2 bp（两轴各 1 bp 反向偏移）时，滑窗身份率只有
~25%（best=31 < 32）→ 整管被跳过。

门控当初是为 91% 的空转调用加的；**anti 序修复后调用数已从 58k 降到 ~900**，
零块 tube 只剩 83 个（18 次调用），门控收益消失、只剩误杀。直接移除。

效果：

| 对 | 移除前 | 移除后 | FastGA |
|---|---:|---:|---:|
| MG1655 vs Sakai | 89.1%（586 块） | **89.3%**（588 块） | 89.3%（506 块） |
| MG1655 vs Nissle | 84.7%（793 块） | **84.7%**（793 块） | 85.3%（711 块） |

Sakai 覆盖与 FastGA **完全持平**；缺失从 15 kb 降到 7.7 kb（剩余为 indel 复杂
小区域）；耗时不变（~1.2s）。教训：门控类启发式在根因修复后要及时复核，
否则会从"省时"变成"漏块"。

## 5.23 并行瓶颈修复：速度与 FastGA 持平（2026-08-02，~0.7s）

阶段计时发现并行效率差：chain 316ms、extend 470ms（32 线程时反而比 8 线程
慢）——不是算法问题而是**资源争用**：

1. **每 tube 的 q 分配**：正链 tube 也 `to_vec()` 复制 5.5 MB 主染色体，857
   个 tube 的分配在 allocator 上串行。改为 `Cow<[u8]>`——正链零拷贝借用、
   仅负链分配 RC：extend 32T 从 470→179ms；
2. **tube 形成串行**：`tubes_for_group` 的相邻桶对合并彼此独立，改用 rayon
   并行（collect 保序）：chain 312→198ms。

最终（默认 32 线程）：

| 对 | chainnet 覆盖 | 耗时 | 峰值内存 |
|---|---:|---:|---:|
| MG1655 vs Sakai | **89.3%**（588 块） | **0.69-0.79s** | ~0.96 GB |
| MG1655 vs Nissle | **84.7%**（793 块） | **0.74s** | 0.79 GB |
| FastGA | 89.3% / 85.3% | 0.7s | ~0 MB |

**速度已与 FastGA 持平**；Sakai 质量持平，Nissle 差 0.6%（7.7 kb indel 复杂
小区域）。剩余：内存（~0.8-1 GB，pgi 全量载入 vs FastGA mmap）、Nissle 0.6%。

## 5.24 内存 33% 削减（2026-08-02，~0.64 GB @ 8 线程）

- `drop(hits)`：tube 形成后 hits（~95 MB）不再需要，提前释放；
- 加上上轮的 q `Cow`（消除每 tube 5.5 MB 正链复制），8 线程峰值内存
  **~960 → 639 MB**，32 线程 ~825 MB；耗时不变（0.69-0.76s）。

剩余内存构成：pgi 全量载入（entries+positions ~274 MB × 2 索引）、并行 dandc
暂存（大 span 的 `split_nd` 数组）、链化临时排序。FastGA 用 mmap 保持 ~0 MB
RSS；pgi 侧做 mmap/惰性加载可进一步降低，但需新依赖或 unsafe。

## 5.25 Nissle 基线核实（2026-08-02，结论修正）

追查 Nissle 0.6% 覆盖差时曾怀疑 `nissle1917.fa.gz` 被替换、FastGA 基线无效。
深入核实后**撤回该结论**：

- `git HEAD` 与工作树的 nissle1917.fa.gz 内容**逐字节相同**（文件未变）；
- 所谓"坐标不一致"是 **naive 偏移身份检查的误区**：该区域含密集 indel
  （每 ~300 bp 一个），naive 同偏移身份 ~25% 但共享 40-mer 达 8990/10258
  （~99% 相关）——对齐真实存在、坐标正确；
- 我们自己的块（t 4508973-4519270 ↔ ns 261268 负链，10231 匹配）经同样
  验证是真实的。

结论：Nissle 的 85.3% FastGA 基线**有效**；我们 84.7% 的 0.6% 差距是真实的
（Nissle 更分歧，indel 复杂区占比更高——同 §5.20 的 adaptamer 部分种子
范畴）。调查过程再次验证"对含 indel 的对齐不能用 naive 偏移身份判断"。

进一步分析（§5.27 前）：抽查最大缺失区（t 2058021-2062339，2.2 kb）发现
**对齐本身是正确的**——wave 产出 t=2057983..2062365（4257 匹配、diffs=135），
PSL 中存在，但 chainnet 在重复区把该块过滤出 syntenic MAF。因此 Nissle 的
"缺失"包含两个成分：(1) chainnet 对重复区单块的过滤差异；(2) 真正无种子的
indel 复杂区。对齐层面的质量与 FastGA 一致。

## 5.26 第三基准验证：MG1655 vs EC958（2026-08-02）

用第三个基准对验证泛化性（ec958.fa.gz 今天 02:23 也改过，但 FastGA 输出与
当前文件坐标一致，基线有效）：

| 方案 | chainnet 覆盖 | 块数 | 耗时 |
|---|---:|---:|---:|
| 我们（tube，默认参数） | **86.2%** | 794 | 0.71s |
| FastGA | 86.3% | 707 | ~0.7s |

覆盖差 **0.1%**（缺失 18.9 kb，但多覆盖 12.6 kb——净差很小）；速度持平。

三个基准汇总（有效基线的两个）：

| 对 | 我们 | FastGA | 差 |
|---|---:|---:|---:|
| MG1655 vs Sakai | 89.3% | 89.3% | 0.0% |
| MG1655 vs EC958 | 86.2% | 86.3% | 0.1% |
| MG1655 vs Nissle | 84.7% | 85.3% | 0.6% |

## 5.27 链化排序索引键与最终状态（2026-08-02）

`chain_tubes` 的并行排序改用 `(u128 key, u32 index)` 元组（替代 `&SeedHit`
引用），瞬态缓冲减半；903 测试全过，质量/速度不变。

进一步的内存优化：`align_to_psl_ext` 改为按值接收 PgiIndex，在链化完成后
`mem::take` 释放 entries/positions（每索引 ~140 MB）——extend 阶段不再持有
k-mer 索引。32 线程峰值内存 **875 → 639 MB**（Sakai）、EC958 566 MB，8/32
线程峰值一致；质量/速度不变。

**最终状态**（8-32 线程，tube 默认参数 = min-shared k/2）：

| 对 | 我们 | FastGA | 差 |
|---|---:|---:|---:|
| MG1655 vs Sakai | 89.3% / 588 块 / 0.7-0.84s | 89.3% / 506 块 / 0.7s | **0.0%** |
| MG1655 vs EC958 | 86.2% / 794 块 / 0.71s | 86.3% / 707 块 | **0.1%** |
| MG1655 vs Nissle | 84.7% / 793 块 / 0.66s | 85.3% / 711 块 | 0.6% |
| 峰值内存 | 0.64-0.88 GB | ~0 MB（mmap） | — |

质量（对齐层面，两个有效基线）与速度均与 FastGA 持平。剩余差距与方向：

1. **Nissle 0.6%**：当时归因为 chainnet 对重复区单块的过滤差异 + 真正无
   种子的 indel 复杂区（需 adaptamer 最小种子选择，`is_minimal` 语义；
   mid-line 窗口身份率预过滤**不可行**，会误杀反向延伸穿过分歧口袋的
   有效调用）。**事后更正**：主要成分是 §5.30 的负链 PSL 坐标 bug，修复后
   差缩至 0.015%（种子选择移植见 §5.29，`is_minimal` 实为 canonical
   方向判断而非噪声抑制）；
2. **内存 0.64-0.88 GB**：pgi 全量载入（entries+positions ~274 MB × 2 索引）
   + 并行 dandc 暂存；FastGA 用 mmap 保持 ~0 MB RSS，mmap/惰性加载需新
   依赖或 unsafe（AGENTS.md 限制）。

`chain_tubes`（`src/libs/pgi/align.rs`）与 wave 引擎
（`src/libs/alignment/wave.rs`）均为独立可测组件。

## 5.28 位置表 u64 位域打包（2026-08-02，内存再降 ~5%）

`PgiIndex::positions` 由 `Vec<(u32, u32, u8)>`（12 B/元素）改为 `Vec<u64>`
位域（pos 32 bit | cid 20 bit | strand 1 bit，8 B/元素），**磁盘格式不变**
（v2 序列化循环只是改了解包/打包方向）。`pack_position`/`unpack_position`
在 `src/libs/pgi/mod.rs`，build/read/write/merge 同步更新；cid 上限 2^20
（约 10^6 contigs，debug_assert 防越界）。

实测（32 线程，tube 默认参数）：

| 对 | chainnet 覆盖 | 块数 | 耗时 | 峰值内存 |
|---|---:|---:|---:|---:|
| MG1655 vs Sakai | 89.3% | 588 | 0.82s | 639 → 607 MB |
| MG1655 vs EC958 | 86.2% | 794 | 0.72s | 566 → 535 MB |
| MG1655 vs Nissle | 84.7% | 793 | 0.77s | 538 MB |

质量/速度不变；904 测试全过（新增 pack/unpack 往返测试）。positions 在
内存构成中占比 ~15%，8 B/记录已接近该结构的物理下限；再降需动 entries
（kmer u128 + 偏移 u32×2 = 24 B，可压缩 kmer 位宽或改 mmap/惰性加载，
后者受 AGENTS.md 新依赖/unsafe 限制）。

## 5.29 adaptamer 种子选择移植：最大共享前缀 + canonical 去重（2026-08-02）

对照 FastGA 源码（`FASTGA-main/FastGA.c` 的 `new_merge_thread`）修正
`merge_seed_hits` 的种子选择语义：

1. **最大共享前缀（plen）**：FastGA 每个 T1（a）条目只对其在 T2（b）中的
   最长匹配发种子——扩展出 `plen = max lcp` 后仅共享 `plen` 碱基的范围参与
   配对；短的部分匹配只有在它是该条目**最长**匹配时才存活。旧实现是固定
   窗口内全发射（每个 lcp ≥ min_shared 的 b 条目都配对），弱种子制造大量
   低质量链。
2. **扩展范围频率过滤**：`freq` 过滤作用于 **plen 处的出现数**
   （`occ < freq` 才保留，FastGA `hgh >= top` 即跳过），而非固定窗口的
   条目数。
3. **canonical 去重**：`pgi build` 每个位置同时存 fwd/RC 两个 key，导致
   每个物理命中发射两次；按 FastGA 单方向存储语义，a 侧只保留
   `kmer <= rc(kmer)` 的 canonical 条目（与 `is_minimal` 同思路——更正
   §5.20 的误解：`is_minimal` 是 canonical 方向判断，不是噪声抑制；真正的
   噪声抑制是 plen 最大选择 + 扩展范围过滤）。
4. **floor = 12**：tube 默认 `min-shared` 由 k/2=20 改为 FastGA 的 plen
   下限 12。配合最大选择后 12-19 bp 的锚点补上 indel 复杂区（§5.9 的
   min-shared 12 灾难是"无最大选择 + 贪心链化"所致，机制不同）。

实测（32 线程，tube 默认参数）：

| 对 | 种子数 | chainnet 覆盖（前 → 后） | 块数 | 耗时 | 峰值内存（前 → 后） |
|---|---:|---:|---:|---:|---:|
| MG1655 vs Sakai | 247 万 | 89.26% → **89.31%** | 589 | 0.78s | 607 → **586 MB** |
| MG1655 vs EC958 | 227 万 | 86.17% → **86.36%** | 845 | 0.71s | 535 → **475 MB** |
| MG1655 vs Nissle | 233 万 | 84.74% → **84.98%** | 842 | 0.76s | 538 → **463 MB** |
| FastGA | 190 万 | 89.3% / 86.3% / 85.3% | — | ~0.7s | ~0 MB |

覆盖三项全部达到或超过 FastGA（Nissle 差 0.32%），速度持平，内存不升反降
（种子减半）。身份率同步上升（Sakai 97.52→97.65%，Nissle 96.68→96.73%）。
906 测试全过（新增最大前缀/扩展范围过滤单元测试）。

## 5.30 负链 PSL 坐标约定 bug：所有 '-' 块被 psl chain 静默丢弃（2026-08-02）

追查 Nissle 剩余差时发现一个**覆盖级的 bug**：`pgi align` 写出的负链 PSL 块
坐标帧与 UCSC/`psl chain` 约定相反——我们输出 qStart/qEnd 用 RC 空间、内部
qStarts 用正链坐标，而约定（FastGA 与 kent 工具）是 **qStart/qEnd 用正链、
内部 qStarts 用 RC 帧**。后果：`psl chain` 的精确打分（`calc_block_score`
按 RC 帧读 qStarts）对每个 '-' 块得到大额负分 → 全部低于 min-score →
**所有负链比对在链化阶段被静默丢弃**。此前 MAF 中 '-' 块数为 0（FastGA
Sakai 3 个 / EC958 11 个 / Nissle 9 个）。

根因：`Psl::from_align` 期望调用方传正链坐标（axt 路径已预先 reverse_range，
有 `axtToPsl.c` 注释为证），而 pgi 两个调用点（`extend_tube`、
`extend_window`）和 `chain_to_psl` 直接传了 orientation/RC 空间坐标。

修复（`src/libs/pgi/align.rs`）：

1. `extend_tube` / `extend_window`：负链时先把对齐坐标 reverse_range 成
   正链再交给 `from_align`（与 axt 路径一致）；
2. `chain_to_psl`：qStart/qEnd 保持正链，qStarts 改为 `b_len - q_end`（RC 帧）。

实测（32 线程，tube 默认参数）：

| 对 | chainnet 覆盖（修复前 → 后） | 块数 | 耗时 | 峰值内存 | FastGA |
|---|---:|---:|---:|---:|---:|
| MG1655 vs Sakai | 89.31% → **89.33%** | 588 | 0.80s | 604 MB | 89.3% |
| MG1655 vs EC958 | 86.36% → **86.38%** | 846 | 0.74s | 512 MB | 86.3% |
| MG1655 vs Nissle | 84.98% → **85.28%** | 847 | 0.74s | 464 MB | 85.30% |

**Nissle 的"0.32% 差距"绝大部分就是这个 bug**（§5.25/§5.27 把原因归为
"chainnet 对重复区单块的过滤差异"是误判——块确实在 PSL 里，但被更前面的
`psl chain` 精确打分环节丢弃，chainnet 从未见过它）。修复后三项全部与
FastGA 持平（Nissle 差 0.015% ≈ 0.7 kb，属噪声级）。906 测试全过，新增
`extend_chain_rc_query` 回归断言：负链 qStarts 必须在 RC 帧且逐段序列
identity 验证通过。

## 5.31 参考索引流式读取（2026-08-02，内存再降 ~34%）

merge 只顺序扫描 a（参考）索引，因此不必全量载入。新增 `PgiStream`
（`src/libs/pgi/mod.rs`，纯 std 分块读取、按条目批量产出、条目不跨批），
`merge_seed_hits_from_stream` 用 rayon `par_bridge` 并行处理批次；
`pgr pgi align` 的 ref 侧改为流式（query 侧仍全量——partition_point
需要随机访问）。`pgr::reader` 返回类型加 `+ Send` 以支持 par_bridge。
此实现修正了 §5.24 里"惰性加载需新依赖或 unsafe"的结论——纯 std 流式即可。

阶段峰值探针（`log::debug`，Linux `/proc/self/status` VmHWM）：

| 对 | 峰值内存（前 → 后） | 耗时 | chainnet 覆盖 |
|---|---:|---:|---:|
| MG1655 vs Sakai | 604 → **398 MB** | 0.83s | 89.33%（不变） |
| MG1655 vs EC958 | 512 → **374 MB** | 0.78s | 86.38%（不变） |
| MG1655 vs Nissle | 464 → **378 MB** | 0.81s | 85.28%（不变） |

覆盖/块数与全量载入完全一致（Sakai 588 / EC958 846 / Nissle 847 块）。
阶段峰值（Sakai）：merge 380→298 MB、chain 471→388 MB、extend 不再
超过 chain 峰值（此前 +100 MB）。剩余构成：b 索引全量（~140 MB）、
链化排序缓冲、rayon 栈；再降需 mmap 或 query 侧惰性访问（新依赖/unsafe，
待用户决策）。907 测试全过（新增 `PgiStream` 全量等价性测试）。

## 5.32 命中打包 + 链排序基数化（2026-08-02，内存再降）

- `SeedHit` 24 → 16 B（a/b contig 与 shared 改 u16，位置仍 u32）：merge 与
  chain 阶段常驻内存约 -19 MB；
- `chain_tubes` 的 `(u128, u32)` 键数组（32 B/元素 ≈ 77 MB）改为并行
  `Vec<u128>` 键 + `Vec<u32>` 序数组（≈ 48 MB）配现有 MSD 基数排序
  （`libs/ds/radix_sort`），并去掉 `Vec<&SeedHit>` 引用数组与 par_sort 的
  隐式键/引用缓冲：链化阶段约 -30 MB。

实测（32 线程，tube 默认；PSL 与全量载入**逐字节一致**）：

| 对 | 峰值内存（前 → 后） | 耗时 |
|---|---:|---:|
| MG1655 vs Sakai | 398 → **~381 MB**（extend 为峰） | 0.79-0.93s |
| MG1655 vs EC958 | 374 → **321 MB** | 0.83s |
| MG1655 vs Nissle | 378 → **310 MB** | 0.86s |

Sakai 的 extend 阶段（wave dandc 并行暂存，+60 MB）重新成为峰值；EC958/
Nissle 仍以链化为峰。907 测试全过（链化排序行为不变，tube 测试原样通过）。

## 5.33 extend 峰值收敛：共享 RC + wave 预留收敛 + 8 线程池（2026-08-02）

追查 Sakai 的 extend 峰值（32 线程时 +60 MB）：

1. **负链 tube 每调用复制整条染色体 RC**（~5.5 MB/调用，并行叠加）：改为
   `b_rcs` 预计算一次、全部负链 tube 共享借用（`Cow::Borrowed`）；
2. **`forward_wave` 预分配 `4096 × width` 单元**（v 8 B + trace 16 B），而
   波数被 TRIM_MLAG 限制在 ~250 附近——预留比实际多几十倍；改为
   `256 × width`（上限 D_CAP=500k），不足时 Vec 自动增长；
3. **并行度**：extend 的 wave/dandc 暂存随并发 tube 数线性增长。`pgr pgi
   align` 新增 `--parallel`（默认 8，即 FastGA `-T` 默认），整个对齐在专用
   rayon 池中执行——8 并发下 extend 峰值不再超过链化阶段，速度不变
   （0.78-0.89s）。

实测（默认 8 线程，tube 默认参数；PSL 逐字节一致）：

| 对 | 峰值内存（前 → 后） | 耗时 | chainnet 覆盖 |
|---|---:|---:|---:|
| MG1655 vs Sakai | 381 → **296 MB** | 0.89s | 89.33%（不变） |
| MG1655 vs EC958 | 321 → **281 MB** | 0.78s | 86.38%（不变） |
| MG1655 vs Nissle | 310 → **284 MB** | 0.78s | 85.28%（不变） |

全流程峰值现在就是链化阶段；剩余 = query 索引全量（~140 MB）+ 链化暂存。
907 测试全过。

## 5.34 多 contig 验证 + 阶段耗时分布（2026-08-02）

- `SeedHit` 的 contig id 收窄为 u16 后，`build_from_seqs` 增加 contig 数
  守卫（> 65535 直接报错，防止静默截断）；
- 新增多 contig 回归测试：3 个 contig，query 第二个为反向互补、第三个带
  分散点突变——覆盖 contig 分组、流式 merge、负链 PSL 约定。

端到端多 contig 实测（k=10，tube）：c1 正链 20000/20000、c2 负链
15000/15000（RC 正确识别）、c3 正链 9858/10000（2% 突变 → 142 错配）。

阶段耗时探针（debug 级，8 线程 Sakai）：merge 193 ms、chain_tubes 236 ms、
extend 245 ms（合计 ~0.67 s；墙钟 0.89 s 另含参考索引流式读取、序列加载、
PSL 写盘）。与 FastGA（~0.7 s）的差距主要在索引读取（mmap 可消除）和
wave 每调用开销。908 测试全过。

## 5.35 query 索引改 mmap 零拷贝（PgiMmap，2026-08-02）

query 索引不再全量读入：新增 `libs::pgi::mmap::PgiMmap`（memmap2，
MAP_PRIVATE 只读映射），记录驻留映射页——条目通过 packed k-mer 字节二分
定位（`entry_range`），位置按需解码（`entry_positions`），与 FastGA GIX
的 memory-mapped 模型一致。reference 仍走 `PgiStream` 流式读取。

- `PgiQuery` trait：resident `PgiIndex` 与 `PgiMmap` 共用的只读视图
  （k/smer/window/contigs + entry_range/entry_next/entry_kmer/entry_freq/
  entry_positions），merge 泛型化；扩展阶段只读 contigs，resident 路径在
  merge 后释放条目表，mmap 路径则从未分配。
- 修过的两个语义坑：
  1. prefix 哨兵键：`hi = lo + r` 可等于 `2^(2k)`，resident 的
     `partition_point(kmer < hi)` 天然容纳，mmap 二分必须把该值 clamp 到
     记录数，否则 key 打包后高位移出、二分错位；
  2. mmap 的 `entry_range` 返回记录区间，组内多条记录不能当作独立条目
     迭代——新增 `entry_next`（resident 为 i+1，mmap 为组尾）按条目推进，
     否则同一 k-mer 组内的每条记录都会被当作条目重复发射命中。
- 实测（合成 2×2 Mb、k=10、smer 4/2、tube + --ref-seq/--query-seq、
  8 线程、release）：旧版（query 全量读入）max RSS 328 MB / 0.80 s；
  mmap 版 298 MB / 0.76 s；merge 191 ms，峰值已被链化暂存（hits/radix/
  tubes）占据。2% 突变 query 的端到端 PSL 输出与旧二进制逐字节一致
  （同一 .pgi 输入）。query.pgi 28 MB → resident 表约 53 MB（positions
  4M×8 + entries 890k×24），Sakai 规模对应 ~140 MB，均已消除。
- 注意：query 索引必须是真实文件（mmap 不支持 stdin/gzip）；`dist pgi`
  / `stat` / `to-hv` 仍走全量 `PgiIndex::read`，暂未改。

## 6. 相关文档

- 索引格式与消费者规划：[[pbit.md]]（多参考节 + .pgi 距离消费者层级）
- FastGA 管线与简化移植评估：[[fastga.md]] §11/§12
- 泛基因组场景：[[ecoli-cohort.md]]、[[paf-pangenome.md]]
