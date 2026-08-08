# 真实 cohort 上 HV 最近邻检索：精确 / HNSW / 物种路由（#6/#7）

> 日期：2026-08-08。135 个真实基因组（同 `bench-hv-ani-calibration.md`），
> HV 向量来自 `pgr pgi to-hv`（k=40 syncmer、D=4096、sparse=1），
> L2 归一化后以欧氏排序 ≡ cosine。真值：skani ANI top-10（生物学）与
> HV 精确 top-10（图检索误差）。对应 `design/genome-nn-query.md`
> §7.4 #6/#7。复现：`cargo bench --bench hv_ann_real`（env 指向
> /tmp/hv_calib 的数据）。

## 结果

| 变体 | recall_ANI@10 | recall_HV@10 | 平均查询 µs |
|---|---|---|---|
| 精确扫描 | 0.663 | 1.000 | 331 |
| 全局 HNSW ef=10 | 0.664 | 0.993 | 177 |
| 全局 HNSW ef=50 | 0.664 | 0.996 | 312 |
| 全局 HNSW ef=100 | 0.664 | 0.996 | 363 |
| 路由 R=1（物种，ef 无关） | 0.542 | 0.699 | 88 |
| 路由 R=2（物种） | 0.641 | 0.899 | 132 |

## 结论

1. **图检索层误差可忽略**：全局 HNSW 的 recall_HV@10 ≥ 0.993（ef≥10），
   recall_ANI 与精确扫描完全相同（0.664 vs 0.663）——HV→ANI 的 0.66
   差距全部来自**距离层**（HV 排序 vs ANI 排序），与
   `bench-hv-ani-calibration.md` 一致。P1 ② 收尾：检索层证据齐了。
2. **物种硬路由在本 cohort 上反而有害**：R=1 时 recall_HV 跌到 0.699，
   R=2 回升到 0.899——原因是很多物种 clade 成员 <10（E. whittamii 3、
   多数 Yersinia 种 1–4），真 top-10 必然跨 clade，硬路由漏掉真近邻。
   这补上了 §6.5 合成实验没覆盖的失败模式：**路由键的 clade 必须足够
   大（≥K 甚至数倍 K），否则路由变成截断**；小 cohort 或小 clade 下
   全局检索（或精确扫描）更稳。
3. **规模提示**：135 个基因组时精确扫描 331 µs/查询（SIMD 4096 维），
   HNSW 约 2× 更快（177 µs）——规模账维持"≤10k 精确够用"结论。
4. 局限：ANI-truth 仅统计有 ≥10 个已知 ANI 邻居的查询（远缘查询被
   排除）；路由实验的物种键来自标签，真实部署若用聚类键需保证
   clade 规模。

## 后续

- 用 E. coli NR 子集（≥1,000 成员的单物种）重测路由：此时物种内
  clade 足够大，路由应恢复 §6.5 合成实验的正收益——验证"clade 规模
  是路由生效前提"。
