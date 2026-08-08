# 真实数据第二批验证：frac CI 校准、groups 分组、ANI 阈值、树一致性

> 日期：2026-08-08。数据与 cohort 同 `bench-hv-ani-calibration.md`
> （135 基因组，skani ANI 全矩阵）。对应 `design/genome-nn-query.md`
> §7.4 的 #5/#13/#15/#18。

## #5 frac 的 ANI 置信区间 vs skani ANI（覆盖率）

方法：`pgr dist frac --merge --ci` 全 cohort 输出每对 ANI 95% CI
（0–1 标度），与 skani ANI（0–100）对比，CI 乘 100 后计算覆盖率。

| ANI 区间 | n | CI 覆盖率 | CI 宽度中位数 | \|CI 中心 − skani ANI\| 中位数 |
|---|---|---|---|---|
| ≥99% | 122 | 15.6% | 0.10 | 0.21 |
| 95–99% | 1,339 | 21.4% | 0.17 | 0.23 |
| 90–95% | 2,427 | 1.4% | 0.51 | 0.92 |
| 85–90% | 721 | 0.0% | 0.67 | 1.27 |
| <85% | 1,328 | 11.8% | 2.42 | 1.91 |
| 全部 | 5,937 | **8.4%** | 0.53 | 0.96 |

结论：frac 的 CI 覆盖的是 **FracMinHash 估计自身的采样误差**（Hera et
al. 2023 公式），不能当作"覆盖金标准 ANI"的区间——frac 与 skani 的
系统偏差（近缘 ~0.2、远缘 ~1.9 ANI 单位）远大于 CI 宽度。文档/帮助文本
应注明 CI 语义，避免用户误读。

## #13 NWR groups.tsv（minhash 树 height 0.4 分组）与物种一致性

groups.tsv 覆盖 680 个代表基因组，**只有 13 个组**（组大小中位数 34，
最大 169）——在 Enterobacterales + Pasteurellales 全数据集上 height
0.4 切出的是**科/目级大组**，不是物种级 clade。

- 物种纯度：中位数 0.029（几乎每组都混合几十个物种）；属纯度中位数
  0.472（31% 的组单属）。
- mash.dist 组内距离中位数 0.197 vs 组间 0.284——区分度弱。

结论：groups.tsv 不能直接当物种级路由键；#7 的真实 clade 实验要用
更低的切割高度（或物种标签）在 cohort 自身距离上重新聚。

## #15 ANI 物种阈值（物种标签为真值的实证）

skani ANI 按物种标签分层（排除 "sp."）：

| 关系 | n | ANI q05 / q50 / q95 |
|---|---|---|
| 同种 | 1,034 | 96.1 / 98.2 / 99.1 |
| 同属异种 | 3,827 | 83.4 / 91.0 / 97.8 |
| 跨属 | 1,076 | 80.6 / 81.1 / 82.6 |

阈值判定（同种 vs 异种）：

| ANI 阈值 | 异种误判为同种（<阈值占比） | 同种误伤（<阈值占比） |
|---|---:|---:|
| 90 | 25.4% | 0.0% |
| 92 | 68.1% | 0.2% |
| 95 | 88.6% | 0.8% |
| 97 | 93.2% | 19.9% |

结论：~95% 阈值与经典物种边界一致（95% 阈值下同种误伤仅 0.8%、
异种漏判 11.4%）；90–92% 更保守。同种内 ANI 下限可达 91.7（个别
异常株/标注噪声），建议 CLI/文档以 95% 为默认"同种"判据、90% 为
"近似同种"提示。

## #18 minhash 距离树 vs bac120 标记基因树（物种级 cophenetic）

注意：两棵树 tip 命名体系不同（旧版 `Atl_hermannii_...` vs 新版
`At_herm_...`），需按 GCF accession 映射；bac120 树另有物种级
`Species___N` tip（`[S=...]` 注释使 necom 无法解析 condensed 版，
用 order 版 + accession 映射解决）。

- 共有物种：137（minhash 675 / bac120 146），物种对 9,316。
- cophenetic Spearman = **0.568**，Pearson = 0.919（深度分裂主导）。
- 分距离段：minhash 距离 0–0.05（近缘物种对）ρ=0.31；0.05–0.15
  ρ=0.39；0.15–0.4 ρ=0.71；>0.4 ρ=0.005。

结论：minhash 距离树与标记基因树在中远缘一致尚可，**近缘物种对的
排序一致性弱（ρ≈0.3–0.4）**——树层面再次印证 k-mer 距离在近缘下的
分辨率缺陷（与 `bench-hv-ani-calibration.md` 的距离层结论一致）。
bac120（保守标记基因）树更适合做物种级参考拓扑。

## 复现

```bash
# #5:  pgr dist frac cohort.fa.lst --list-files --merge --ci -o frac.ci.tsv
# #13: groups.tsv + genome.taxon.tsv + mash.dist.tsv（纯 python 分析）
# #15: ani.full.tsv + genome.taxon.tsv（纯 python 分析）
# #18: necom nwk label/distance + accession 映射（analyze_trees4.py）
```
