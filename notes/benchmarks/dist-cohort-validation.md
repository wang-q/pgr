# 距离消费者 cohort 验证（dist pgi / dist hv / dist seq）

> 目的：验证 `.pgi` 的两个距离消费者（`dist pgi` 精确归并、`pgi to-hv` +
> `dist hv` 近似向量）在 10 株 E. coli 上的距离语义，与 `dist seq` 草图距离
> 及比对身份率（`pgr align pgi` 实测，见 [[../design/pgi-align.md]] §2.4）
> 对照。日期：2026-08-02。

## 方法与数据

- 10 株 × 45 对：`dist pgi`（k=40、syncmer 8/5、精确集合归并）、
  `dist hv`（`pgi to-hv` dim=1024/4096/8192，k=40）、
  `dist seq`（syncmer k=8/w=5 草图）；
- 参考真值：45 对比对身份率（扩展 PSL 块）。

## 排序一致性（Spearman ρ，与 1 − 身份率）

| 方法 | ρ | 备注 |
|---|---:|---|
| `dist seq`（k=8 草图） | **0.816** | 全基因组组成采样，最贴近真值 |
| `dist pgi`（k=40 精确） | 0.539 | 确定性但有偏（见下） |
| `dist hv`（dim 1024） | **−0.05** | 饱和退化，不可用（见下） |

## 三个发现

### 1. `dist pgi` 是"采样集合的精确距离"，不是全组成距离

两索引归并对共享 k-mer 集合的计数是**确定性的**（零采样方差），但集合本身
只含 syncmer 位置采样的 40-mer：点突变会改变局部 s-mer 哈希、进而改变
syncmer 选择，两侧位置集漂移远快于真实组成差异。MG1655 vs Sakai（98.4%
身份率）的 k-mer containment 只有 0.54、jaccard 0.33，mash 0.0175 高估
分歧。Spearman 0.539 说明排序大致跟随但受株系重复内容干扰
（如 ec2011c_3493–se11 身份率 99.2% 却给出 mash 0.038）。

### 2. `dist hv` 对大规模集合饱和退化（dim 无关）

`hash_hv_i8` 下每个 k-mer 种子向**所有维度**累加随机 i8，260 万种子使各维
值 ~±1.3e6；共享种子主导点积，`dot ≥ card` 后被 `calc_distances` 截断，
近缘株 containment 全部饱和为 1.0、mash 压缩 ~10×（0.0002–0.0046），
与身份率 Spearman ≈ 0。dim 1024 → 4096 → 8192 不改善（期望与 dim 无关）。
**结论**：`pgi to-hv` + `dist hv` 的"4 万级粗筛"定位需要重做 HV 距离
（如阈值化双极编码 + 余弦/Hamming，或按种子数缩放），当前实现不可用于
近缘基因组排序。

> **已修复（2026-08-02）**：改为**稀疏投影 + 余弦**（`.hv` v2）：
> 每个 k-mer 只更新 `--sparse`（默认 3）个随机维度（±1），文件头存
> `sparse` 与 `n_kmer`；比较时余弦相似度估计共享数
> `inter = cos × √(n1·n2)`。默认 `--dim 4096`。
> 结果：与 `dist pgi` 的 mash **Spearman 0.9694**，共享 k-mer 计数
> 平均误差 2.4%（最大 13%）；45 对比较 14.4s → **0.29s（50×）**。
> 稠密编码（`hash_hv_i8`，dim 1024–8192 均饱和）与余弦改算均无效
> （余弦范围 0.994–0.999，信息已淹没），稀疏是有效修复。

### 3. `dist seq` 草图是当前最优的粗筛层

k=8 syncmer 集合（4^8 = 65k 可能值，接近全组成）与身份率 Spearman 0.816，
快且稳健——4 万级近邻过滤继续用 `dist seq`/`dist hv`（FASTA 侧，k=8），
`dist pgi` 作为"已建索引时的确定性精确层"，`dist hv`（.hv）待修复后再用。

> **更新（2026-08-02）**：`.hv`（稀疏）现在是 `dist pgi` 的 50× 快、
> 0.97 排序一致的近似层；`dist seq` 仍是与身份率最贴近的草图层
> （0.82 vs 0.51），因为 k=40 syncmer 集合受采样位置漂移限制。

## 附带交付：`dist hv` 支持 .hv 文件

按设计文档 §7.2 补齐了消费者链：`pgr dist hv a.hv b.hv` 直接比较
`pgi to-hv` 的输出（校验 k/dim 一致，复用 `calc_distances` 核心），
含集成测试（自比较 mash=0/jaccard=1、维度不匹配报错）。

## 相关文档

- 索引格式与消费者规划：[[../design/pbit.md]]（.pgi 距离消费者层级）
- 比对身份率真值：[[../design/pgi-align.md]] §2.4
