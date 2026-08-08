# bac120 标记基因路由准确率（#19）

> 日期：2026-08-08。对应 `design/genome-nn-query.md` §7.4 #19。
> 目的：§6.5 的"生物学先验"落地方案——用保守标记基因做廉价查询路由，
> 测它在真实数据上的准确率。

## 方法

- 数据：135 cohort 基因组；`Domain/bac120.fa.gz`（369,814 个标记蛋白）
  与 `Domain/seq_asm_f3.tsv`（蛋白→组装→标记映射，旧名体系，按 GCF
  accession 映射到 cohort）。
- 路由键：**8 个全覆盖 bac120 标记**（PF01025.14、PF00466.15、
  TIGR00019、PF02576.12、TIGR00061、TIGR00054、TIGR00043、TIGR00029；
  135/135 覆盖）拼接成每基因组蛋白集。
- 距离：`pgr dist frac --protein --merge --scale 100`（默认 scale 1000
  对蛋白太稀，sketch 仅 1 元素——本身是个发现）。
- 评估：留一法——查询的最近 aa 邻居所属物种 = 路由结果；准确率 =
  路由物种 == 查询自身物种。对照：HV（pgi-to-hv 点积）路由、ANI
  top-1 物种（标签噪声上限）。

## 结果

| 路由方式 | 留一法准确率 |
|---|---|
| 8×bac120 标记蛋白（aa frac） | **0.756** |
| HV（pgi-to-hv 点积） | 0.822 |
| ANI top-1 物种（金标准上限） | 0.800 |

- 全局精确（HV 点积排序）recall@10 vs ANI = 0.711（n=135）。
- 路由后"只搜路由物种 clade"的 recall 无法统计（n=0）：路由目标物种
  clade 大多 <10 成员——再次印证 #7 的"小 clade 路由失效"。

## 结论

1. **8 个保守标记基因的 aa 距离即可路由 75.6% 的查询到正确物种**，
   接近 ANI 金标准上限（80%），只比全基因组 HV 路由低 7 pp——作为
   §6.5 的"生物学先验"，成本极低、效果可用。
2. 物种标签与 ANI top-1 只有 80% 一致（近缘株标注噪声/种内近邻跨
   物种），路由准确率的天花板就在 ~80%，标记路由已接近。
3. 蛋白 frac 默认 scale=1000 对短蛋白（~300 aa）过稀（sketch 1 元素，
   距离全 0/1），蛋白场景应降 scale（如 100）——建议 `dist frac
   --protein` 文档提示。
4. 路由 recall 仍需大 clade 场景验证（E. coli 全量 NR 是天然试验场）。

## 复现

```bash
# 1) cohort_marker.tsv：seq_asm_f3.tsv 按 accession 映射到 cohort
# 2) pgr fa some bac120.fa.gz <prots> → 按基因组拆分
# 3) pgr dist frac genomes.lst --list-files --merge --protein --scale 100
# 4) 留一法分析：/tmp/hv_calib/analyze_marker.py
```
