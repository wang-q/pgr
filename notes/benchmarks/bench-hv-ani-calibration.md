# HV 距离 vs ANI 金标准标定（真实 Enterobacterales 基因组）

> 目的：给 `design/genome-nn-query.md` §7 的 P1 补第一块硬证据——HV
> 距离与真实 ANI 的相关性、分辨率区间与召回。日期：2026-08-08。
> 结论先行：HV 距离只在 ANI 90–98% 区间中等可靠（Spearman 0.5–0.6），
> **≥99% 与 <85% 失效**；同条件下 Mash 与 ANI 几乎完全相关（ρ≈0.97–0.99）。
> 物种内（≥98% ANI）选参考 / 聚类不应以 HV 距离为主。

## 数据与 cohort（严格挂靠 NWR 指导文件）

- 数据源：`~/data/Escherichia/`（Enterobacterales + Pasteurellales，
  150k 组装，132,572 个通过 QC，见 `~/Scripts/genomes/groups/Escherichia.md`）。
- cohort：135 个基因组，全部来自各物种 **NR.lst**（非冗余代表）且都在
  `summary/pass.lst`（QC 通过）内；物种标签取自 `summary/genome.taxon.tsv`。
  亲缘分层：E. coli NR 40（近缘，ANI≥98%）+ 其他 Escherichia 种 60
  （E. albertii/fergusonii/marmotae/ruysiae/whittamii/sp、Pseudescherichia，
  中等，ANI 88–97%）+ Yersinia 36（远缘，ANI<88%）。全两两 = 9,045 对。
- 样本清单与映射：`/tmp/hv_calib/cohort.meta.tsv`（name/species/path，
  临时目录，可重建）。

## 方法

- **HV 距离**：`pgr dist hv --list-files --parallel 8`，默认 DNA minimizer
  k=21/w=5；D=4096 与 D=16384 各跑一遍（输出第 7 列 = Mash 式距离
  d = −(1/k)ln(2J/(1+J))，1−d ≈ ANI 估计值）。
- **Mash 距离**：直接用 NWR 已算好的每基因组 `.msh` sketch
  （`MinHash/<species>/msh/<name>.msh`），`mash triangle -E -p 8`（与
  NWR `dist.sh` 同款参数），9045 对全覆盖。
- **ANI 金标准**：skani 0.1.0 `dist --ql/--rl -t 8 --min-af 0`（取满覆盖；
  极远缘无比对命中者视为未知，5,937/9,045 对有 ANI）。skani 是全基因组
  ANI（近似 BLAST-ANI），本文以它为真值；GSearch 用 BLAST-ANI/FastANI，
  二者同级别。
- **指标**：Spearman（排序）/Pearson/RMSE（HV 估计 ANI vs skani ANI，
  RMSE 按 0–1 标度）；recall@10 = 以 ANI 为真值取 top-10，与 HV/Mash
  距离 top-10 的交集比例（自比对剔除，按查询新颖度分层）。分析脚本：
  `/tmp/hv_calib/analyze_ani.py`、`analyze_dim.py`（临时）。

## 结果

### HV(1−d) vs skani ANI（D=4096）

| 分层 | n | Spearman | Pearson | RMSE(0–1) |
|---|---|---|---|---|
| 全部 | 5,937 | 0.882 | 0.481 | 0.191 |
| 同种内 | 1,115 | 0.610 | 0.124 | 0.168 |
| 种间 | 4,822 | 0.861 | 0.504 | 0.195 |
| ANI ≥99% | 122 | **0.383** | 0.238 | 0.045 |
| ANI 95–99% | 1,339 | 0.608 | 0.115 | 0.154 |
| ANI 90–95% | 2,427 | 0.496 | 0.184 | 0.072 |
| ANI 85–90% | 721 | 0.462 | 0.378 | 0.137 |
| ANI <85% | 1,328 | **0.054** | 0.068 | 0.344 |

### Mash vs skani ANI（参考，同数据）

| 分层 | n | Spearman | Pearson |
|---|---|---|---|
| 全部 | 5,937 | −0.990 | −0.982 |
| 同种内 | 1,115 | **−0.974** | −0.987 |
| 种间 | 4,822 | −0.983 | −0.977 |

### HV 距离分位数（分辨率直观检查，D=4096）

| ANI 区间 | n | hv_dist q05 / q50 / q95 |
|---|---|---|
| ≥99% | 122 | 0.016 / 0.043 / 0.065 |
| 95–99% | 1,339 | 0.034 / 0.058 / 0.099 |
| 90–95% | 2,427 | 0.087 / 0.109 / 0.139 |
| 85–90% | 721 | 0.108 / 0.126 / 0.206 |
| <85% | 1,328 | 0.152 / 0.202 / 1.000 |

相邻 ANI 区间的 hv_dist 分布大量重叠（95–99% 与 90–95% 的 q05–q95
区间几乎连续）——排序分辨率差的直接体现。

### recall@10（真值 = skani ANI top-10，135 个查询）

| 方法 | 总体 | 新颖度 ≥98% | 新颖度 90–98% |
|---|---|---|---|
| HV（D=4096） | 0.622 | 0.612 (n=124) | 0.727 (n=11) |
| Mash | 0.762 | 0.762 (n=124) | 0.764 (n=11) |

### D=4096 vs D=16384（分辨率是否随维度改善）

| ANI 区间 | Spearman D=4096 | Spearman D=16384 |
|---|---|---|
| ≥99% | 0.383 | 0.387（**无改善**） |
| 95–99% | 0.608 | 0.608（无改善） |
| 90–95% | 0.496 | 0.577 |
| 85–90% | 0.462 | 0.558 |
| <85% | 0.054 | 0.322 |
| recall@10 总体 | 0.622 | 0.629 |

## 结论

1. **HV 距离的可靠区间是 ANI 90–98%（及部分 85–90%），且仅中等
   （Spearman 0.5–0.6）**；≥99% 的近缘株排序几乎失效（ρ≈0.38，且
   D=16384 无改善——不是维度饱和度，是方法固有噪声）；<85% 远缘完全
   失效（ρ≈0.05–0.32）。
2. **Mash 是 ANI 的可靠代理**（同种内 ρ=−0.97），recall@10 比 HV 高
   14 pp——种内近缘排序任务 Mash 明显更优。
3. **对设计的直接影响**：物种内（≥98% ANI）聚类 / 选参考应以
   `dist mash` / `dist frac` 为主；HV 适合做**嵌入 / 粗筛 / 查询路由**
   （85–98% 带），不适合做 ANI 级精排，更不能替代 skani/fastANI 标定。
4. 增大 D（4096→16384）只改善中远缘，代价是 4× 内存/时间；对
   近缘分辨率无帮助，参数上不必为近缘场景升级 D。
5. 待补：完整度鲁棒性、sampler/k 扫描、HNSW 检索在真实 HV 上的
   ANI-truth 召回（§7.2 ② 的图检索部分）。

## 复现

```bash
mkdir -p /tmp/hv_calib
# 1. 从 NR.lst/pass.lst 抽样 135 基因组 -> cohort.meta.tsv / cohort.fa.lst
# 2. HV:  pgr dist hv cohort.fa.lst --list-files --parallel 8 -o hv.tsv
# 3. Mash: mash triangle -E -p 8 -l cohort.msh.lst > mash.tsv
# 4. ANI:  skani dist --ql cohort.fa.lst --rl cohort.fa.lst -t 8 --min-af 0 -o ani.full.tsv
# 5. 分析: python3 analyze_ani.py / analyze_dim.py
```

## 补充：四种距离统一对标 ANI（#1，2026-08-08）

同 cohort、同 ANI 真值，追加 `pgr dist frac --merge`（k=21, scale=1000）
与 `pgr dist mini --merge`（k=21/w=5, rapid）。Spearman 绝对值：

| 方法 | 全部 | 同种内 | ≥99% | 95–99% | 90–95% | 85–90% | <85% |
|---|---|---|---|---|---|---|---|
| HV | 0.882 | 0.610 | 0.383 | 0.608 | 0.496 | 0.462 | 0.054 |
| Mash | 0.990 | 0.974 | 0.805 | 0.972 | 0.961 | 0.924 | 0.676 |
| **frac** | **0.991** | **0.973** | 0.796 | **0.970** | **0.961** | 0.919 | 0.673 |
| mini | 0.917 | 0.612 | 0.401 | 0.607 | 0.632 | 0.583 | 0.600 |

recall@10（真值 = skani ANI top-10）：

| 方法 | 总体 | ≥98% | 90–98% |
|---|---|---|---|
| HV | 0.621 | 0.612 | 0.727 |
| Mash | 0.762 | 0.762 | 0.764 |
| frac | 0.757 | 0.759 | 0.736 |
| mini | 0.629 | 0.613 | 0.809 |

**结论补充**：① `dist frac` 与 Mash 同为 ANI 的可靠代理（ρ 0.97–0.99，
recall 与 Mash 相当），"frac 用于 ANI 估计"的既有建议得到实证支持；
② **minimizer 采样（mini）与 HV 有相同的近缘分辨率缺陷**（同种内与
≥99% 区间 ρ≈0.4–0.6）——问题是采样层（minimizer）而非 HV 编码本身；
③ mini 在 <85% 区间（ρ 0.60）明显好于 HV（0.05），远缘下 minimizer
草图仍保有信息，HV 的 4096 维编码是远缘失效的主因。
