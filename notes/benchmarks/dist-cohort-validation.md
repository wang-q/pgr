# 距离消费者 cohort 验证（dist pgi / dist hv / dist seq）

> 目的：验证 `.pgi` 的两个距离消费者（`dist pgi` 精确归并、`pgi to-hv` +
> `dist hv` 近似向量）在 10 株 E. coli 上的距离语义，与 `dist seq` 草图距离
> 及比对身份率（`pgr align pgi` 实测，见 [[../design/pgi-align.md]] §2.4）
> 对照。日期：2026-08-02 初测、2026-08-05 复测。

## 方法与数据

- 10 株 × 45 对：`dist pgi`（k=40、syncmer 8/5、精确集合归并）、
  `dist hv`（`pgi to-hv` dim=1024/4096/8192，k=40）、
  `dist seq`（syncmer k=8/w=5 草图）；
- 参考真值：45 对比对身份率（扩展 PSL 块）。

## 排序一致性（Spearman ρ，与 1 − 身份率）

### 初测（2026-08-02）

| 方法 | ρ | 备注 |
|---|---:|---|
| `dist seq`（k=8 草图） | **0.816** | 全基因组组成采样，最贴近真值 |
| `dist pgi`（k=40 精确） | 0.539 | 确定性但有偏（见下） |
| `dist hv`（dim 1024，稠密饱和） | **−0.05** | 饱和退化，不可用（见下） |

### 复测（2026-08-05）

| 方法 | ρ | 备注 |
|---|---:|---|
| `dist seq`（k=8 syncmer 草图，首 contig 对，jaccard） | **0.616** | 最贴近身份率真值 |
| `dist pgi`（k=40 精确，mash） | 0.590 | 确定性但有偏（见发现 1） |
| `dist hv`（稀疏 v2，dim 4096，mash） | 0.547 | 与 dist pgi 的 mash 排序一致（ρ=0.969） |

身份率真值 = pooled PSL identity（`(matches+rep)/block_len`，口径见
[[../design/pgi-align.md]] §2.2），2026-08-05 用当前二进制重算 45 对：
分布 0.9544-0.9879（均值 0.9682）。排序结论不变：`dist seq` 最贴近
身份率、`dist pgi` 次之、`dist hv` 与 `dist pgi` 高度一致。

### s=1 复测（2026-08-08，5 株 × 10 对）

`pgr pgi to-hv` 默认稀疏度 s 由 3 改为 1 后（`design/hv.md` §2.7 决策：
s 不影响精度、s=1 编码快 2.5×），用 5 株（MG1655 / Sakai / SE11 /
E2348 / CFT073，GCF 基因组）× 10 对验证：

| 对 | `dist hv` s=1 mash | `dist pgi` mash | 差 |
|---|---:|---:|---:|
| mg1655–se11 | 0.0126 | 0.0124 | +0.0002 |
| mg1655–sakai | 0.0180 | 0.0176 | +0.0004 |
| sakai–se11 | 0.0186 | 0.0186 | 0.0000 |
| mg1655–cft073 | 0.0256 | 0.0260 | −0.0004 |
| se11–cft073 | 0.0261 | 0.0270 | −0.0009 |
| sakai–cft073 | 0.0269 | 0.0278 | −0.0009 |
| e2348–cft073 | 0.0448 | 0.0438 | +0.0010 |
| sakai–e2348 | 0.0568 | 0.0557 | +0.0011 |
| mg1655–e2348 | 0.0569 | 0.0553 | +0.0016 |
| se11–e2348 | 0.0594 | 0.0569 | +0.0025 |

**Spearman ρ = 0.9879**（mash 排序），最大差异 0.0025——s=1 与精确
`dist pgi` 高度一致，排序与 s=3 时代一致（原 s=3 测 ρ=0.969，45 对）。
再次确认 s 不影响距离排序（`design/hv.md` §2.7 理论 + 实测）。

### FracMinHash 无偏性验证（2026-08-08，全集合对照）

新增 `dist seq --sampler frachash`（§5.5 方向 2）后，用**全 k-mer 集合
（不经采样）**验证各采样器的 Jaccard 估计（e2348 × cft073，k=40）：

| 方法 | Jaccard | 备注 |
|---|---:|---|
| **全 k=40 集合（真值）** | **0.4508** | 所有 canonical k-mer 精确交集 |
| **FracMinHash（k=40, scale=100）** | **0.4171** | 随机采样，**无偏**（差 0.034） |
| `dist pgi`（k=40 syncmer 8/5） | 0.0948 | **syncmer 位置偏差，严重低估 4.7×** |

**重要发现**：`dist pgi` 的 syncmer 8/5 采样在两个（近缘、含共享重复元件的）
基因组上有**位置偏差**——即使全集合重叠 45%，syncmer 位置因微小序列差异
而错开，交集比例远低于真值。FracMinHash（独立随机采样）无此问题，给出
接近全集合的 Jaccard。这也解释了此前 `dist pgi` "确定性但有偏"（§2 发现 1）
的量级。**FracMinHash 是 pgr 距离体系中唯一无偏的采样器**（排序类任务
各采样器大致可用，数值 ANI 必须用 FracMinHash）。

**偏差范围细化（2026-08-08，5 株 × 10 对 vs 全 k=40 真值）**：

| 组 | 对数 | 平均 ANI 偏差 |
|---|---:|---:|
| 不含 e2348 的对 | 6 | **~0.3 pp**（`dist pgi` 可靠） |
| 涉及 e2348 的对（EPEC，噬菌体/IS/质粒丰富） | 4 | **~3.3 pp**（数值不可信） |

`dist pgi` 不是"本身不准"：对"干净"基因组对偏差小（~0.3% ANI），
偏差集中在**重复/移动元件丰富的株**（syncmer 位置在重复区更难对齐，
交集低估被放大）。排序整体仍可用（偏差方向一致、按相似度缩放）。

**Jaccard vs Containment（2026-08-08，10 对 vs 全集合真值）**：

| 指标 | 平均相对偏差 | 排序 Spearman |
|---|---:|---:|
| Containment | **34.3%** | **0.661** |
| Jaccard（→ mash） | 39.3% | 0.515 |

containment 略准：分母是单集合（syncmer 采样率 ~57%），偏差被稀释；
jaccard 的分母 union 受 inter 低估放大。**但两者相对偏差仍 ~35%**、
排序 ρ 仅 0.5–0.7——syncmer 偏差对两个指标都有显著影响，`dist pgi`
的数值判断（无论 jaccard 还是 containment）均不可靠，无偏数值 ANI
必须用 `dist seq --sampler frachash`。

**FracMinHash 的 containment（2026-08-08，10 对 vs 全 k=40 真值）**：

| 方法 | 相对偏差 | 排序 Spearman |
|---|---:|---:|
| **FracMinHash containment** | **9.9%** | **0.976** |
| pgi containment（对照） | 34.3% | 0.661 |

FracMinHash containment = inter/card1（`dist seq --sampler frachash` 的
containment 列），**期望无偏**（共享 k-mer 判定一致，分子分母采样率 p
抵消）；残差 ~10% 是比值估计的 **Jensen 偏差**——实测 scale=10/100/1000
偏差几乎相同（一阶理论：Cov 与 E² 都 ∝ p²，p 抵消），**增大采样不能消除，
需偏差校正公式才能进一步精确**（Hera et al. 2023 方向）。排序 ρ=0.976
对粗筛/质粒检测足够；数值精确的 containment 留作校正公式后续。

> **更正（2026-08-08，同 k 对照实验）**：上述 ~10% 残差是 **k 不匹配的
> 假象**——frachash 用 k=21、真值用全 k=40 集合，两套 k-mer 的共享比例
> 本来就不同。同 k=21 下（5 株 × 10 对，全 canonical k-mer 集合真值）：
> containment 偏差仅 0.3–2.3%（正负对称，即采样方差，SE≈0.7%），
> jaccard 偏差 0.1–1.6%。读 Hera et al. 2023 原文后确认：其校正因子
> (1−(1−s)^|A|) 对 |A|≥10⁵ 的基因组 ≈ 1（4.6 Mb 时 ≈ e^(−4600)≈0），
> 仅对 |A|<~100 且 scale≥0.1 的极短序列有意义；Theorem 8 的 CI 需对 p
> 数值求解，与现有正态近似差异小。**结论：大肠杆菌/默认 scale 场景
> 无需 Hera 校正**，`dist frac` 的 containment/ANI 保持现状（数据见
> `benchmarks/bench-dist-mash-compat.md` 同 k 对照节）。

### `dist hv --sampler syncmer` 偏差对照（2026-08-08，5 株 × 10 对）

此前 syncmer 偏差数据只覆盖 `dist pgi`（k=40 smer=8/w=5）。补上 `dist hv
--sampler syncmer`（smer=8/w=55，HV 投影）的对照，让 syncmer 偏差在两个
实验入口都可体验（hv 为 FASTA 直算、pgi 需先建索引）：

| 对 | hv syncmer jac | frac jac | 真值 k=21 jac | hv con | frac con | 真值 con |
|---|---:|---:|---:|---:|---:|---:|
| cf×e2 | 0.9467 | 0.5219 | 0.6680 | 0.9673 | 0.6772 | 0.5150 |
| cf×mg | 0.9175 | 0.3398 | 0.4819 | 0.9328 | 0.4861 | 0.3393 |
| cf×sa | 0.9190 | 0.3108 | 0.4875 | 0.9594 | 0.4823 | 0.3132 |
| cf×se | 0.9292 | 0.3180 | 0.4844 | 0.9605 | 0.4795 | 0.3217 |
| e2×mg | 0.9204 | 0.3484 | 0.4976 | 0.9394 | 0.5012 | 0.3462 |
| e2×sa | 0.9198 | 0.3246 | 0.5168 | 0.9652 | 0.5049 | 0.3301 |
| e2×se | 0.9216 | 0.3309 | 0.5021 | 0.9617 | 0.5004 | 0.3294 |
| mg×sa | 0.9269 | 0.4486 | 0.6690 | 0.9894 | 0.6591 | 0.4494 |
| mg×se | 0.9452 | 0.5471 | 0.7474 | 0.9947 | 0.7347 | 0.5535 |
| sa×se | 0.9281 | 0.4259 | 0.5884 | 0.9583 | 0.5840 | 0.4327 |

**发现**：hv syncmer 的 jaccard（0.92–0.95）远高于 k=21 真值（0.48–0.75），
看似"虚高"，但这不是 HV 方法的缺陷——**k=8 全集合真值本身 jaccard≈0.999**
（8-mer 空间 65536 基本饱和，近缘基因组 8-mer 几乎全重合），hv syncmer
相对 k=8 基线**低估 ~5–8%**（syncmer 位置漂移使共享 8-mer 丢失），与
`dist pgi`（k=40 锚定）的偏差**机制相同、方向一致**（都是交集低估），
只是 k=8 基线极高掩盖了量级。

**排序一致性**（10 对距离排序 Spearman）：

| 对照 | ρ |
|---|---:|
| hv syncmer vs 真值 | 0.636 |
| frac vs 真值 | 0.782 |
| hv syncmer vs frac | 0.612 |

结论强化：syncmer 家族（`dist pgi` 与 `dist hv --sampler syncmer`）数值
均不可信（pgi 严重低估、hv 因 k=8 饱和无区分度），排序也弱于 `dist frac`
（ρ 0.64 vs 0.78）。体验入口：`dist hv --sampler syncmer`（FASTA 直算）与
`dist pgi`（索引）；数值 ANI 一律 `dist frac`。

**syncmer 家族（dist pgi 与 dist seq --sampler syncmer 同机制）**
**（2026-08-08 补充）**：两者都是 closed syncmer 采样，位置偏差同源。
k=21 时 `dist seq syncmer` 比 frachash 低 ~2%（e2-cf 0.5075 vs
0.5191、mg-se 0.5451 vs 0.5548）；pgi（smer=8 锚定 k=40）更严重
（4.7×），因"小 s-mer 窗口锚长 k-mer"放大位置漂移。**数值距离
（jaccard/containment/ANI）一律用 FracMinHash**；syncmer 家族只承担
排序/粗筛（FASTA 侧 `dist seq syncmer`）与比对锚点（pgi/align）。
**分层使用建议**：**初筛用 `dist hv`**（.hv 路径，O(D) 固定比较——
大规模近邻/粗筛/聚类，hv.md §1.1 的定位）；候选对的**中等精度距离**
用 `dist seq --sampler frachash`（无偏 + CI）或 `dist pgi`（注意重复
元件株的偏差）；最终验证用 `align pgi`。

> **注意（2026-08-08 用户澄清）**：上述"位置偏差"是**采样集合交集**的
> 固有行为，不是 pgi 的缺陷——pgi 的 syncmer 是 **align pgi（FastGA
> 兼容基因组比对）的锚点基石**：closed syncmer 由 k-mer 序列自身判定
> （同步性），共享 k-mer 在两个基因组中确定性双选 → 锚点完整保留；
> FracMinHash 共享 k-mer 的双选概率是 1/scale（hash 相同、判定完全相关，
> 非 1/scale²）——但 scale=1000 时锚点只保留 1/1000（4.6Mb 约 4600 个）
> 且无窗口保证（随机分布、长 gap 无锚），链化不可用。**不能改 pgi
> 采样器**。正确分工：pgi + syncmer = 比对/排序；无偏数值 ANI =
> `dist seq --sampler frachash`（§5.5 方向 2，含 CI）。

### FracMinHash 排序验证（2026-08-08，5 株 × 10 对 vs 全集合真值）

用**全 k=40 canonical k-mer 集合**（不经任何采样）作为真值，对比
`dist seq --sampler frachash`（k=40、scale=100）的 mash 排序：

* **Spearman ρ = 1.0000**（10 对相对顺序完全一致）；
* FracMinHash mash 平均比真值高 0.0027（Jaccard 略低估）——比值估计的
  Jensen 偏差（E[X/Y] ≠ E[X]/E[Y]），scale 越大越小，不影响排序；
* 最近对 mg1655–se11、最远对 sakai–cft073 均与真值一致；e2348–cft073
  实际是第二近（J=0.45，真值），此前 `dist pgi` 的"远"是 syncmer 偏差。

## 三个发现

### 1. `dist pgi` 是"采样集合的精确距离"，不是全组成距离

两索引归并对共享 k-mer 集合的计数是**确定性的**（零采样方差），但集合本身
只含 syncmer 位置采样的 40-mer：点突变会改变局部 s-mer 哈希、进而改变
syncmer 选择，两侧位置集漂移远快于真实组成差异。MG1655 vs Sakai（97.3%
身份率）的 k-mer containment 只有 0.54、jaccard 0.33，mash 0.0175 高估
分歧。Spearman 0.590 说明排序大致跟随但受株系重复内容干扰
（如 ec2011c_3493–se11 身份率 97.9% 却给出 mash 0.038）。

### 2. `dist hv` 对大规模集合饱和退化（dim 无关）

`hash_hv_i8` 下每个 k-mer 种子向**所有维度**累加随机 i8，260 万种子使各维
值 ~±1.3e6；共享种子主导点积，`dot ≥ card` 后被 `calc_distances` 截断，
近缘株 containment 全部饱和为 1.0、mash 压缩 ~10×（0.0002–0.0046），
与身份率 Spearman ≈ 0。dim 1024 → 4096 → 8192 不改善（期望与 dim 无关）。
**结论**：`pgi to-hv` + `dist hv` 的"4 万级粗筛"定位需要重做 HV 距离
（如阈值化双极编码 + 余弦/Hamming，或按种子数缩放），当前实现不可用于
近缘基因组排序。

> **已修复（2026-08-02，复测 2026-08-05）**：改为**稀疏投影 + 余弦**
> （`.hv` v2）：
> 每个 k-mer 只更新 `--sparse`（2026-08-08 起默认 1，原默认 3）个随机
> 维度（±1），文件头存
> `sparse` 与 `n_kmer`；比较时余弦相似度估计共享数
> `inter = cos × √(n1·n2)`。默认 `--dim 4096`。
> 结果：与 `dist pgi` 的 mash **Spearman 0.969**（当前），
> 共享 k-mer 计数平均误差 2.39%（最大 12.94%，复测与初测 2.4%/13% 一致）；
> 45 对比较 8.25s → **0.12s（复测，~71×）**（初测 14.4s → 0.29s，50×）。
> 稠密编码（`hash_hv_i8`，dim 1024–8192 均饱和）与余弦改算均无效
> （余弦范围 0.994–0.999，信息已淹没），稀疏是有效修复。

### 3. 粗筛分层（2026-08-08 更新）

* **初筛（超大规模近邻过滤）**：`dist hv`（.hv 稀疏路径，**O(D) 固定
  比较**，s=1 已验证排序 ρ=0.988 vs `dist pgi`，§1）或 FASTA 侧
  `dist seq`（k=8 syncmer，Spearman 0.616 vs 身份率，快且稳健）。
* **二次距离（候选对）**：`dist seq --sampler frachash`（无偏 + CI）或
  `dist pgi`（已建索引时，注意重复元件株的偏差，见上）。
* **最终验证**：`align pgi`（精确比对）。

> **更新（2026-08-02，复测 2026-08-05）**：`.hv`（稀疏）现在是 `dist pgi`
> 的 ~70× 快、0.97 排序一致的近似层；`dist seq` 仍是与身份率最贴近的
> 草图层（当前 0.62 vs 0.59），因为 k=40 syncmer 集合受采样位置漂移限制。

## 附带交付：`dist hv` 支持 .hv 文件

按设计文档 §7.2 补齐了消费者链：`pgr dist hv a.hv b.hv` 直接比较
`pgi to-hv` 的输出（校验 k/dim 一致，复用 `calc_distances` 核心），
含集成测试（自比较 mash=0/jaccard=1、维度不匹配报错）。

## 相关文档

- 索引格式与消费者规划：[[../design/pbit.md]]（.pgi 距离消费者层级）
- 比对身份率真值：[[../design/pgi-align.md]] §2.4
