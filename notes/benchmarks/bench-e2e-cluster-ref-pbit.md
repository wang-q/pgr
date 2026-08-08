# 端到端验证：聚类 → 选参考 → 比对 → PBit 归档（2026-08-08）

对应 `design/genome-nn-query.md` §7.2 ③④⑤（P2）与核心工作流
"物种内聚类选参考 + PBit 归档"。目的：量化**参考选择策略**对
归档压缩率的影响，为实施方案提供决策依据。

## 数据与子集

- 数据源：`~/data/Escherichia/`（Enterobacterales NR，E. coli 全库
  51,318 组装；cohort = 2,088 个 E. coli NR 基因组，
  `/tmp/hv_calib/meta2115.tsv` + `mash2115.tsv` 全对距离现成）。
- Pilot：随机 30 个；全量：**farthest-point 采样 100 个**（基于
  mash2115.tsv 距离，seed 42，保证覆盖物种内多样性）。
- 所有样本 FASTA 路径校验存在（`ASSEMBLY/Escherichia_coli/<name>/`）。

## 方法

1. `pgr dist mash --merge --list-files` 全对距离（30 样本 435 对 /
   100 样本 4,950 对；0.01–0.04 mash 距离 = ANI ~96–99%）。
2. `necom mat to-phylip` + `necom clust upgma` + `necom cut simple`：
   30 样本 k=4（23/3/3/1）；100 样本 k=7（44/35/9/8/2/1/1）。
3. 每簇按三种策略各选 1 个参考：
   - **center**：簇内到其他成员 mash 距离和最小；
   - **longest**：FASTA 解压总长最大；
   - **random**：固定 seed 随机。
4. 簇内每个非参考样本 × 参考：
   `pgr align pgi`（PSL）→ `pgr psl to-paf` → `pgr pbit create`
   （纯 LZ 与 `--paf` CIGAR 两条路径）。
5. delta 压缩率 =（pbit 归档 − 参考 self-archive）/ gzip-9 样本大小；
   `pgr pbit to-fa` 抽查覆盖率。

## 结果

### Pilot（30 样本，26 查询 × 3 策略 = 78 对）

| 策略 | LZ delta/gzip（mean） | 备注 |
|---|---|---|
| center | 0.499 | 中心参考最差 |
| longest | **0.466** | 簇内最长参考最优 |
| random | 0.490 | 介于两者 |

### 全量（100 样本，98 查询 × 3 策略 = 279 对）

| 策略 | LZ delta/gzip mean | median | min | max |
|---|---|---|---|---|
| center | 0.5535 | 0.5505 | 0.4486 | 0.6327 |
| longest | **0.5198** | 0.5172 | 0.4334 | 0.6378 |
| random | 0.5213 | 0.5262 | 0.4302 | 0.6276 |

逐查询配对（n=98）：**longest < center 71/98**、random < center 61/98；
longest vs random 无稳定差异（39/49）。按簇看，longest 在 4/5 个
有样本的簇最优或接近最优；random 在 clade 1 偶然更优（0.496 vs
longest 0.532），说明**单参考的随机性影响可达 3–4 pp**。

### 参考特征（解释 longest 优势）

| 策略 | 参考 contig 数（中位） | 参考总长（平均） |
|---|---|---|
| center | 52 | 4.81 Mb |
| longest | **28** | **5.38 Mb** |
| random | 74 | 4.99 Mb |

center 参考（簇内距离和最小）倾向选"典型 draft"，内容覆盖量反而
最少；压缩率主要由**参考的内容覆盖量（总长 × 完整度）**决定，而非
与样本的平均相似度。

### 覆盖与 CIGAR 路径

- to-fa 覆盖率抽查（8/8）：100 版 ≥0.9984（KTE66 draft 极端 0.16%；
  其余 ≤0.02%），完整参考下多为 100%——LZ 内容匹配基本无损。
- CIGAR（`--paf`）与纯 LZ 归档大小差平均仅 32 B（max 143 B）：
  已知段相位约束使 CIGAR 路径基本回退 LZ（#14e），端到端无差异。

### Complete vs draft 参考（配对对照，2026-08-08 补充）

用户提出生产规则："选作参考的应该必须是 Complete"。本 cohort
（farthest-point 100 个 E. coli）中组装级别分布：Contig 48 /
Scaffold 27 / Complete Genome 22 / Chromosome 3（仅 25% 达标），
7 个簇中有 **2 簇完全没有 Complete/Chromosome 成员、1 簇仅 1 个**。

配对实验（控制"簇内最长"变量，只变组装级别）：clade 0/1 各选
complete-longest（48/GF60，Complete Genome，5–9 contigs 含质粒）与
draft-longest（49832/531，Contig，~200 contigs）两个参考，对同一批
75 个查询做 pbit（LZ）：

| 参考 | delta/gzip（mean） | complete 更优的对数 |
|---|---|---|
| complete-longest | 0.5249 | 19/75 |
| draft-longest | **0.5049** | —（draft 更优 56/75） |

即**压缩率维度 Complete 参考不占优**（draft-longest 反而 -2 pp）。
机制假说：① 查询 75% 是 draft，其 contig 颗粒与 draft 参考的 contig
天然对齐（LZ 内容匹配的 `best_ref_group` 命中率高）；② draft 参考含
质粒/未装配片段，与 draft 查询共享更多内容；③ Complete 参考的 4 kb
参考段是染色体上切出的，与查询 contig 边界错位。覆盖率两种参考都
≈100%（抽查），差异在编码效率而非覆盖。

**方案含义**：Complete 参考的正当理由不在压缩率，而在比对质量、
坐标可解释性与下游一致性；生产规则应定为"参考优先 Complete
（Chromosome 次之），若簇内无 Complete 则用簇内最长 draft（压缩率
不差）或并入相邻簇"，而不是"必须 Complete"。

## 结论与决策建议

1. **参考选择影响压缩率 ~3.4 pp（longest vs center）**；在 E. coli
   物种内（ANI 96–99%）最长/高完整度参考是稳定优选的简单启发式。
2. **不要用"距离中心"当参考**：中心参考是典型 draft，内容覆盖少，
   压缩率反而最差。压缩率维度上"最长 draft" ≥ "Complete"（配对
   实验，见上节）——**完整度不是压缩率的充分条件**；但生产规则仍
   建议 Complete 优先（比对/可解释性/一致性理由，见上节方案含义）。
3. 单参考随机性影响 3–4 pp：小簇（<10 成员）时参考选择比簇划分更
   影响结果；多参考（每簇 2–3 个）可摊平随机性（待验证）。
4. 端到端流程（dist → necom → pgi → pbit）全部就绪、耗时可控：
   100 样本 × 3 参考策略 ≈ 12 分钟（8 线程），pgi align ~1.5 s/对、
   pbit ~3.5 s/样本。

## 复现

```bash
# 子集选择（farthest-point，seed 42）
cd /tmp/e2e   # 脚本与中间产物
python3 - <<'PY'   # 见 bench 文档注释；cohort100.tsv/list100.txt 已生成
PY
pgr dist mash --merge --list-files -p 8 list100.txt -o dist100.pair
necom mat to-phylip dist100.pair.clean -o dist100.phylip
necom clust upgma dist100.phylip -o tree100.nwk
necom cut simple tree100.nwk -k 7 -o clust100_k7.tsv
python3 run_e2e_100.py   # 279 任务：align + pbit（LZ/CIGAR）
```
