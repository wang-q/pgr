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

**剩余差距与方向**：

1. 性能：把 D&C 回溯换成 wave 内嵌的 Pebble 稀疏 trace（trace point 间隔
   `tspace=100`，消费端补全间隙），消除每次调用的 O(span) 重建开销；
2. **调用数**：真实差距在种子结构——FastGA 的链在分歧区断开（tube≈块，
   每 tube ~1 次调用）。但简单部分种子（`--min-shared 20`，3.9s / 754 块）
   实测更慢：部分匹配使 hit 爆炸、tube 更密。FastGA 的优势在其 adaptamer
   **选择**（稀疏最小种子 + 方向不对称），不是部分匹配本身，需按
   `libfastk.c` 的 `is_minimal` 语义移植；mid-line 窗口身份率预过滤
   **不可行**（会误杀反向延伸穿过分歧口袋的有效调用，覆盖 88.2%→55%）；
3. 与 FastGA 对齐 `BUCK_ANTI`（FastGA 为加倍空间 128 = 未加倍 64）后重测。

`chain_tubes`（`src/libs/pgi/align.rs`）与 wave 引擎
（`src/libs/alignment/wave.rs`）均为独立可测组件。

## 6. 相关文档

- 索引格式与消费者规划：[[pbit.md]]（多参考节 + .pgi 距离消费者层级）
- FastGA 管线与简化移植评估：[[fastga.md]] §11/§12
- 泛基因组场景：[[ecoli-cohort.md]]、[[paf-pangenome.md]]
