# pgr 自研 RepeatMasker 替代的可行性评估

> 设计评估笔记，日期：2026-08-03。场景笔记见 [[../repeat-masking.md]]（ir/rept/trf + fa mask
> 的现状与 SD 关系）；FastK 分析见 [[../references/fastk.md]]。

## 1. 动机

RepeatMasker 基于现有重复库（Dfam/RepBase）与基因组做比对，用的比对器是特制版
RMBlast（blastn 变体），对大型基因组非常耗时。pgr 已实现快速比对（`pgr align pgi`、
`pgr align lastz`），且 `pgr sd search --engine pgi` 已证明"小集合 vs 基因组"的
比对模式可行。因此讨论：能否自研一个比 RepeatMasker 快得多的替代品？

结论先行：**比对部分可行且会快很多；真正的成本在注释后处理。如果只做遮蔽
（不输出 family/class 注释），实现很轻，但增量价值有限，需先验证比对敏感度再决定。**

## 2. RepeatMasker 工作流拆解

RepeatMasker ≠ "库 vs 基因组比对"。它分两大部分：

### 2.1 比对（最慢，但也是 pgr 最容易超越的部分）

*   RMBlast（特制 blastn）把库中 consensus 对基因组找所有拷贝；
*   默认分两个严格度跑（先高严格度找 anchor，再补低严格度）；
*   对 3 Gb 基因组 + 几千条 query，blastn 的种子扩展极慢。

### 2.2 注释后处理（核心价值，二十年打磨所在）

*   **two-pass 碎片整合**：转座子元件常被打断成多个 hit（分歧、内部插入/删除），
    要拼回完整元件；
*   **边界精修**：Smith-Waterman / profile 校准精确边界；
*   **family/class 注释**：从库的 FASTA header 解析（Dfam 格式如 `>XXX#LTR/ERV`）；
*   **%div/%del/%ins**：Kimura 2-parameter 距离，估计转座时间；
*   **低复杂度判定**：以低复杂度为主的 hit 标为 `Simple_repeat`/`Low_complexity`；
*   **输出**：`.out` 标准表、`.tbl` 汇总、masked 序列。

## 3. pgr 现有能力映射

| RepeatMasker 步骤 | pgr 现状 | 判断 |
| :--- | :--- | :--- |
| 库-基因组比对 | `sd search --engine pgi`（或 lastz）已验证同模式 | 高可行，预计比 RMBlast 快一个数量级 |
| 碎片整合 | `libs/chain` + `pgr pl chainnet` | 高可行（pgr 独特优势） |
| 边界精修 | `libs/alignment` / POA | 高可行 |
| family/class 注释 | 需解析库 header | 简单，按 Dfam 格式写解析 |
| K2P %div | 无 | 小函数，不难 |
| 低复杂度判定 | 部分（`pgr pl trf`） | 最大缺口之一 |

## 4. 三个方案

### 方案 A：完整 RepeatMasker 替代（含注释）

可行，但不建议一步到位：转座子生物细节（LTR 5'/3' 边界、TSD、polyA、嵌套插入、
重排元件）是泥潭，字节级复现 `.out` 会陷入无尽细节。

### 方案 B：只做遮蔽（推荐先评估）

遮蔽 = 找区间 + 覆盖。第 2 步已有（`pgr fa mask --runlist`），所以核心只剩
"区间找得准不准"。**注意：pgr 现在就有遮蔽功能**——`ir/rept/trf` 输出 runlist，
配 `fa mask` 即完成遮蔽（E. coli 全基因组 0.35s）。因此方案 B 的实际工作
是"把 FastK 近似换成真比对"，增量价值取决于敏感度。

### 方案 C：维持现状

`ir + trf + fa mask`，把精力放在补 low-complexity 过滤（见 §6）。

## 5. 遮蔽版的最小实现

不做注释时实现很轻，基础设施全在：

1.  **比对**：复用 `sd search --engine pgi` 引擎（或 lastz）跑"库 vs 基因组"。
    注意 `sd search` 的过滤器（>1 kb、>90% identity）是给 SD 调的，转座子拷贝
    分歧大（70–90%），需放宽 min-len / identity 或跳过过滤器只取 hits。
2.  **区间合并**：`spanr cover / merge / fill` 一行管道。遮蔽不在乎把碎片拼回
    "完整元件"，只在乎别漏区域。
3.  **输出**：runlist → `pgr fa mask`。

工作量比完整 RepeatMasker 小一个数量级。

## 6. 关键风险

*   **比对敏感度**：k-mer（k=17）对高分歧拷贝会漏；pgi 的 syncmer 种子对
    70% identity 的拷贝同样不轻松。这是唯一值得做方案 B 的理由，必须实测。
*   **低复杂度缺口**：RepeatMasker 默认屏蔽 low complexity（polyA、卫星、
    homopolymer）。现有 `ir` 只管库内散在重复，`trf` 覆盖串联重复，polyA 这类
    不一定被覆盖。这是遮蔽质量上更实际的差距，与用 k-mer 还是比对无关。
*   **验证基准**：E. coli 几乎无转座子，无参考价值。需用拟南芥/玉米等
    转座子丰富基因组，与 RepeatMasker 的 masked 输出对比 recall。
*   **库门槛低**：做遮蔽只需 Dfam 的 consensus FASTA，不需要 HMM 或分类信息。

## 7. 建议的验证实验（先验证再决定）

1.  取转座子丰富基因组（拟南芥或玉米）；
2.  用 Dfam consensus FASTA 经 lastz/pgi 对一遍，放宽 hit 过滤；
3.  对比三组数据：
    *   时间与 hits 数量；
    *   覆盖区间 vs 现有 `ir` 的差异；
    *   与 RepeatMasker masked 输出的 recall。
4.  若覆盖明显更好且速度可接受 → 值得做轻量遮蔽命令；否则维持方案 C。

## 8. 结论

*   完整替代（方案 A）：不推荐一步到位，后处理是泥潭。
*   遮蔽版（方案 B）：实现轻、快，但增量价值取决于敏感度，**先跑 §7 的验证**。
*   现状（方案 C）：已覆盖 90% 需求；若验证显示比对版提升有限，把精力花在
    low-complexity 过滤上更划算。
