# 真实数据聚类验证：Necom 聚类 vs 物种标签、距离噪声稳定性（#11/#12）

> 日期：2026-08-08。cohort 同 `bench-hv-ani-calibration.md`（135 基因组，
> 16 个物种目录/14 个有效物种标签）。对应 `design/genome-nn-query.md`
> §7.4 #11/#12。

## #11 Necom 聚类 vs 物种标签

方法：mash（NWR .msh）与 HV（k21/w5/D4096）距离矩阵 → `necom mat
to-phylip` → `necom clust hier --method ward` → `necom cut simple -k K`
→ `necom eval partition` vs 物种标签（cluster format）。

| 距离 | K | ARI | 同质度 homogeneity | RI |
|---|---|---|---|---|
| mash | 10 | **0.739** | 0.785 | 0.931 |
| mash | 16 | 0.681 | 0.878 | 0.935 |
| mash | 25 | 0.588 | 0.902 | 0.924 |
| mash | 40 | 0.333 | 0.924 | 0.900 |
| HV | 10 | 0.573 | 0.711 | 0.899 |
| HV | 16 | 0.646 | 0.835 | 0.924 |
| HV | 25 | 0.525 | 0.909 | 0.916 |
| HV | 40 | 0.256 | 0.928 | 0.894 |

结论：物种级划分（K≈10–16）两种距离都能较好恢复（ARI 0.57–0.74），
Mash 各 K 下都优于 HV（+0.07–0.17）；K 过大（40）时 ARI 骤降（过度
细分，同质度上升但一致率崩）。注意 E. coli 与 Escherichia sp. 混聚在
一起（sp. 标签实际多为 E. coli），物种标签本身有噪声。

## #12 距离噪声下的聚类稳定性

方法：HV 距离矩阵加乘性高斯噪声（σ = 5/10/20/40%），每水平 3 个种子
重聚类（ward, K=16），与原聚类算 ARI。

| 噪声 σ | ARI（3 次重复均值 ± SD） |
|---|---|
| 5% | 0.918 ± 0.014 |
| 10% | 0.821 ± 0.141 |
| 20% | 0.726 ± 0.110 |
| 40% | 0.363 ± 0.101 |

结论：物种级聚类对 ≤20% 距离噪声稳健（ARI ≥0.73）；40% 噪声破坏聚类。
结合 `bench-hv-ani-calibration.md`（HV 近缘区间排序噪声大），物种级
划分不受影响，但**种内亚结构（如 E. coli 内部的分群）对距离噪声更
敏感**，后续应专门用 E. coli 子集验证。

## 复现

```bash
# 距离矩阵 → phylip → ward 聚类 → cut → eval：
# necom mat to-phylip mash.pair -o mash.phylip
# necom clust hier --method ward mash.phylip -o tree.nwk
# necom cut simple -k 16 tree.nwk -o clust.tsv
# necom eval partition clust.tsv --other species.pair \
#     --input-format cluster --other-format pair
# 稳定性脚本：/tmp/hv_calib/run_stability.py
```
