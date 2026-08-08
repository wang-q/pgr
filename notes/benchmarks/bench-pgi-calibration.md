# pgi 距离/投影标定（#16/#17）

> 日期：2026-08-08。子集 10 基因组（5 E. coli + E. albertii/fergusonii/
> marmotae + Yersinia enterocolitica/pseudotuberculosis），45 对中有
> ANI 的 30 对。对应 `design/genome-nn-query.md` §7.4 #16/#17。

## 方法

- `pgr pgi build`（默认 k=40、syncmer s=8/w=5，0.6 s/基因组，41 MB/pgi）。
- `pgr dist pgi a.pgi b.pgi`（k-mer 合并确定性 Jaccard/Mash 距离）。
- `pgr pgi to-hv`（D=4096、sparse s=1）→ `pgr dist hv a.hv b.hv`。
- 对照：FASTA 直算 `pgr dist hv --sampler syncmer -k 8 -w 5`。

## #16 pgi 距离 vs skani ANI

| 分层 | n | Spearman |
|---|---|---|
| 全部 | 30 | −0.918 |
| ANI ≥95% | 10 | −0.705 |
| ANI 90–95% | 15 | −0.551 |
| ANI <90% | 5 | −0.800 |

结论：pgi（k=40 syncmer）距离整体与 ANI 高度相关，但近缘段（≥95%）
排序弱——与 HV/mini 的规律一致（syncmer 采样同样有近缘分辨率问题）。

## #17 pgi→HV 投影的保距性

- pgi-hv vs pgi 直接距离：Spearman **0.966**（Pearson 0.66，非线性，
  斜率 ~3.6）——HV 投影忠实保留 pgi k-mer 集的距离排序。
- pgi-hv vs FASTA syncmer HV（-k 8 -w 5）：Spearman 0.785，且距离
  尺度差 ~20 倍（同一对 pgi 0.07–0.09 vs fasyn ~0.0035）——两条 HV
  路径参数语义不同（pgi 的 syncmer span 与 dist hv 的 s-mer/window
  定义不一致），**不可直接当同一管线比较**，需显式对齐参数。

## 结论

1. `pgr pgi to-hv` 可作为 pgi 索引 → HV 嵌入的保距投影（ρ≈0.97），
   HV 侧检索可复用于 pgi 距离语义。
2. pgi 距离与 ANI 的相关性与 syncmer 类采样一致（近缘弱），pgi 侧
   证据与 dist 家族结论互相印证。
3. 文档应注明：`.hv` 文件自带采样参数头，`dist hv` 跨管线比较前须
   核对参数一致（pgi to-hv 默认 syncmer 40/8/5，FASTA 直算默认
   minimizer 21/5）。
