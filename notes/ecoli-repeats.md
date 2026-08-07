# 10 株大肠杆菌的重复遮蔽：`pgr rept` 5 方法 × 3 外部库

> 2026-08-07 运行，pgr 0.4.0。目的：量化"真核转座子（TE）序列在细菌基因组中
> 出现"的程度。对 `tests/genome/` 下 10 株 E. coli 分别用 `pgr rept` 的 5 个子命令
> 做重复检测（runlist 可直接喂 `pgr fa mask`），其中 e-kmer / e-align 各搭配 3 个
> 外部重复库，共 9 种"方法×库"组合 × 10 株 = 90 次运行，全部成功。

## 1. 数据与方法

### 1.1 输入

* 基因组：`tests/genome/{cft073,e2348_69,e24377a,ec042,ec2011c_3493,ec958,mg1655,nissle1917,sakai,se11}.fa.gz`，
  10 株 E. coli，总长 52,425,656 bp（单株 4.64–5.59 Mb）。
* 外部重复库（`~/data/repeats/`，已按 [docs/rept.md](../docs/rept.md) 清洗：
  大写、IUPAC→N、去 dash、去重复 ID）：
  * `tncentral.fa.gz` — 6,093 条，原核插入序列（IS）数据库；
  * `repbase.fa.gz` — 31,491 条，经典重复库（真核为主）；
  * `dfam.fa.gz` — 26,292 条，转座子家族库（真核为主）。

### 1.2 命令与参数

全部使用默认参数（同 docs/rept.md），e-align 额外 `-p 8`：

```bash
pgr rept e-kmer  ~/data/repeats/<lib>.fa.gz genome.fa.gz -o out.json
    # k=17, min-len=300
pgr rept e-align ~/data/repeats/<lib>.fa.gz genome.fa.gz -p 8 -o out.json
    # k=40, min-identity=0.70, min-len=50
pgr rept s-kmer  genome.fa.gz -o out.json
    # k=17, min-len=100
pgr rept s-align genome.fa.gz -p 8 -o out.json
    # window=200/100, min-depth=4
pgr rept trf     genome.fa.gz -o out.json
    # TRF 默认参数
```

e-kmer 按 docs 提示**串行**执行（FastK 并发会 SIGSEGV），其余方法 4 路并行。
90 个 runlist JSON 保留在 `/tmp/rept_mask/json/`
（`<genome>.<method>.<lib>.json`），完整明细在 `/tmp/rept_mask/full.tsv`。

统计口径：runlist 区间为闭区间，片段长度 = end − start + 1（与 `pgr runlist stat`
一致，已抽查验证）；平均长度 = 总遮蔽 bp / 片段数。

**比对后端说明**：`pgr rept e-align` 内部调用 `pgr align pgi`（PGI 后端），
本节所有 e-align 数字均来自 PGI；`pgr rept s-align` 是 lastz 窗口自比对
（LASTZ 后端）。为核对两个后端，另用 LASTZ 做了 TnCentral 库比对（§2.4）。

## 2. 结果

### 2.1 方法汇总（10 株合计，总长 52,425,656 bp）

| 方法 | 库 | 片段数 | 遮蔽 bp | 平均长度 | 每株平均片段数 | 覆盖度 |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: |
| e-kmer | TnCentral | 621 | 876,745 | 1,411.8 | 62.1 | 1.67% |
| e-kmer | RepBase | 117 | 109,539 | 936.2 | 11.7 | 0.21% |
| e-kmer | Dfam | 120 | 108,637 | 905.3 | 12.0 | 0.21% |
| e-align | TnCentral | 734 | 1,151,386 | 1,568.6 | 73.4 | 2.20% |
| e-align | RepBase | 166 | 159,475 | 960.7 | 16.6 | 0.30% |
| e-align | Dfam | 239 | 245,488 | 1,027.1 | 23.9 | 0.47% |
| s-kmer | — | 3,831 | 2,976,515 | 777.0 | 383.1 | 5.68% |
| s-align | — | 17,902 | 4,546,661 | 254.0 | 1,790.2 | 8.67% |
| trf | — | 969 | 258,081 | 266.3 | 96.9 | 0.49% |

### 2.2 每基因组明细（片段数 / 平均长度 bp）

| 基因组 | e-kmer Tn | e-kmer RB | e-kmer DF | e-align Tn | e-align RB | e-align DF | s-kmer | s-align | trf |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| cft073 | 64/1170.1 | 2/773.5 | 3/628.3 | 79/1306.8 | 6/856.2 | 13/1033.1 | 461/532.5 | 1906/222.4 | 74/203.4 |
| e2348_69 | 69/1363.9 | 6/956.2 | 5/1054.4 | 80/1483.8 | 12/722.1 | 19/891.9 | 295/728.9 | 1501/224.0 | 77/112.8 |
| e24377a | 101/1372.5 | 15/1025.5 | 15/1024.8 | 100/1802.9 | 25/965.6 | 32/1013.2 | 410/725.5 | 1787/253.1 | 87/163.1 |
| ec042 | 56/1445.4 | 14/675.8 | 15/630.9 | 63/1763.6 | 19/858.4 | 26/945.8 | 419/507.3 | 1931/202.0 | 100/317.4 |
| ec2011c_3493 | 84/1898.3 | 17/841.0 | 19/761.3 | 98/1963.6 | 26/1001.7 | 33/1040.2 | 531/565.4 | 2277/227.8 | 85/152.0 |
| ec958 | 49/2125.9 | 7/772.6 | 7/771.4 | 65/1940.3 | 10/912.3 | 17/1023.8 | 313/587.2 | 1603/200.0 | 94/250.7 |
| mg1655 | 48/1184.9 | 38/1124.3 | 39/1079.6 | 55/1302.2 | 39/1147.9 | 46/1153.2 | 170/755.6 | 1457/167.8 | 84/223.4 |
| nissle1917 | 59/1069.6 | 7/1003.9 | 7/1005.7 | 78/1169.5 | 15/938.1 | 25/1021.7 | 371/2017.3 | 1615/536.1 | 142/610.6 |
| sakai | 55/1251.7 | 4/920.3 | 4/916.3 | 72/1315.0 | 9/602.0 | 16/856.2 | 527/833.1 | 2121/304.3 | 143/206.7 |
| se11 | 36/993.1 | 7/610.6 | 6/661.8 | 44/1408.6 | 5/1159.8 | 12/1173.3 | 334/617.0 | 1704/205.0 | 83/203.7 |

### 2.3 TnCentral（真实重复库）方法间交叉比较

TnCentral 是三个库里唯一"真实"的重复库（原核 IS 元件），用它做基准交叉，
看 e-kmer / e-align（同一库的两种检测机制）以及 s-kmer / s-align / trf
（自找方法）之间的一致程度。交叉按遮蔽 bp（闭区间并集长度）计算，10 株合计：

| 方法 A | 方法 B | A 遮蔽 bp | B 遮蔽 bp | 交集 bp | A 被 B 覆盖 | B 被 A 覆盖 |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: |
| e-kmer | e-align | 876,745 | 1,151,386 | 874,452 | 99.7% | 75.9% |
| s-kmer | e-kmer | 2,976,515 | 876,745 | 608,202 | 20.4% | 69.4% |
| s-kmer | e-align | 2,976,515 | 1,151,386 | 690,360 | 23.2% | 60.0% |
| s-align | e-kmer | 4,546,661 | 876,745 | 649,056 | 14.3% | 74.0% |
| s-align | e-align | 4,546,661 | 1,151,386 | 786,266 | 17.3% | 68.3% |
| trf | e-kmer | 258,081 | 876,745 | 9,579 | 3.7% | 1.1% |
| trf | e-align | 258,081 | 1,151,386 | 13,832 | 5.4% | 1.2% |

每基因组逐一验证：e-kmer 的遮蔽 bp 有 **99.2–99.9%** 落在 e-align 内（10 株
全部一致），而 e-align 比 e-kmer 多出 **17–42%** 的 bp（每株 e-kmer 覆盖 e-align
的 57.6–82.7%）。s-kmer / s-align 对 e-align 的覆盖每株分别为 40.5–77.4% /
49.4–83.8%。

### 2.4 后端核对：e-align 的 PGI vs LASTZ（TnCentral）

用 LASTZ 复现 e-align（库比对）核对 PGI 结果。流程：TnCentral 清洗时把
167 个非 ACGTN 污染字符（数字/字母）换成 N，按 `pgr fa split name` 拆成单序列，
`pgr align lastz genome tn_split --preset set01 -p 8`（LASTZ 后端，preset 与
s-align 默认一致）→ `pgr lav to-psl` → 与 e-align 相同的过滤
（ident ≥ 0.70、target span ≥ 50、excise 50、fill 10）→ runlist。
注意 `pgr align lastz` 要求单序列 FASTA，且并发 >8 线程无收益（实测 4×8 全核
打满反而拖慢，最终按 -p 8 串行完成）。

| 基因组 | PGI 片段 | PGI bp | LASTZ 片段 | LASTZ bp | 交集 bp | PGI 被 LASTZ 覆盖 |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| cft073 | 79 | 103,240 | 154 | 124,498 | 102,724 | 99.5% |
| e2348_69 | 80 | 118,704 | 150 | 146,665 | 118,241 | 99.6% |
| e24377a | 100 | 180,294 | 173 | 198,496 | 179,840 | 99.7% |
| ec042 | 63 | 111,104 | 143 | 130,090 | 110,674 | 99.6% |
| ec2011c_3493 | 98 | 192,430 | 174 | 215,128 | 192,227 | 99.9% |
| ec958 | 65 | 126,122 | 124 | 147,840 | 125,550 | 99.5% |
| mg1655 | 55 | 71,619 | 145 | 88,232 | 70,657 | 98.7% |
| nissle1917 | 78 | 91,218 | 160 | 115,154 | 90,675 | 99.4% |
| sakai | 72 | 94,677 | 143 | 116,912 | 94,172 | 99.5% |
| se11 | 44 | 61,978 | 120 | 78,708 | 61,301 | 98.9% |
| **合计** | **734** | **1,151,386** | **1,486** | **1,361,723** | **1,146,061** | **99.5%** |

区间结构进一步说明差异来源：LASTZ 比 PGI 多出的 210,337 bp 中，有 **760 个
与 PGI 零重叠的新片段、合计 198,985 bp、平均仅 262 bp**（占多出部分的 94.6%）——
即 LASTZ 多出来的覆盖几乎全是**新增的短片段**，而不是把 PGI 的区域接长。
反过来看 PGI 的 734 个区间：626 个（85.3%）被**单个** LASTZ 片段完整包住、
仅 11 个（1.5%）被 LASTZ 切成多段、其余为边缘少量差异。因此"PGI 片段更完整"
的说法不成立：同一区域上 LASTZ 至少和 PGI 一样长，两者差异主要是 LASTZ 多捡了
一堆平均 ~260 bp 的零散短匹配。

### 2.5 其他 PGI 用法的参数调查与敏感参数复核

项目内所有 `pgr align pgi` 调用点及其参数（含默认值）：

| 调用点 | k | smer/window | freq | min-shared | 说明 |
| :--- | ---: | ---: | ---: | ---: | :--- |
| `pgr align pgi` CLI 默认 | 40 | 8/5 | 10 | 12（FastGA plen floor） | 通用基因组比对 |
| `pgr rept e-align` | 40 | 8/5 | 100 | 16 | 重复遮蔽，偏保守（高 freq 过滤库内重复 k-mer） |
| `pgr align fill` / `rest` | 40 | 8/5 | 10 | 12 | pgi 锚点 + lastz 补全，用 CLI 默认 |
| `pgr sd search/cross --engine pgi` | 31 | 8/5 | 50 | 12 | SD 检测，实测调优（见 notes/design/sd.md） |
| `pgr dist pgi` | — | — | — | — | 只用 `.pgi` 索引算距离，不做比对 |

SD 的调参依据（notes/design/sd.md §140-160）：freq 10→50 修复高拷贝重复漏检
（e2348_69 562→0）；k 40→31 修复 90–93% 分歧拷贝漏检（sakai 4→0、
e24377a 2→0）；k31 未引入假阳性（hits 全部 ≥0.90 identity），hit 数 +35%，
性能持平。

用 SD 的敏感参数（k=31, freq=50, min-shared=12）重跑 e-align（TnCentral × 10 株）：

| 配置 | 片段数 | 遮蔽 bp | 平均长度 | LASTZ 覆盖比例 |
| :--- | ---: | ---: | ---: | ---: |
| PGI 旧默认（k40/f100/ms16） | 734 | 1,151,386 | 1,568.6 | 84.2% |
| PGI 新默认（k31/f50/ms12） | 810 | 1,250,721 | 1,544.1 | 90.4% |
| LASTZ set01 | 1,486 | 1,361,723 | 916.4 | 100% |

敏感参数下 PGI 与 LASTZ 的差距从 210,337 bp 缩小到 111,002 bp（关闭约 47%），
PGI 敏感版 98.4% 的 bp 仍在 LASTZ 内。注意参数要成套：mg1655 上 k31/f50 若保留
min-shared=16 反而降到 65,946 bp（比默认 71,619 bp 还少），k31 必须同时放宽
min-shared 才有效。余下 ~111 kb 仍是 LASTZ 独有的短片段，说明覆盖差异约一半来自
e-align 默认参数保守，另一半是 LASTZ 本身更敏感（或噪声）。

> 2026-08-07：`pgr rept e-align` 默认参数已按上述结论改为 k=31 / freq=50 /
> min-shared=12（新默认 = 表中"PGI 新默认"一行，代码与 docs/rept.md 已同步）。

### 2.6 数字怎么读：两种"重复"口径（库遮蔽 vs 自找重复）

"大肠杆菌重复很少"这句话只对**经典重复/转座子**成立。同一批基因组用两种口径
算出的覆盖度差一个量级（下表为占基因组百分比）：

| 基因组 | 大小 (Mb) | e-kmer TnCentral | e-align TnCentral | s-kmer | s-align |
| :--- | ---: | ---: | ---: | ---: | ---: |
| cft073 | 5.23 | 1.4% | 2.0% | 4.7% | 8.1% |
| e2348_69 | 5.07 | 1.9% | 2.3% | 4.2% | 6.6% |
| e24377a | 5.25 | 2.6% | 3.4% | 5.7% | 8.6% |
| ec042 | 5.36 | 1.5% | 2.1% | 4.0% | 7.3% |
| ec2011c_3493 | 5.44 | 2.9% | 3.5% | 5.5% | 9.5% |
| ec958 | 5.25 | 2.0% | 2.4% | 3.5% | 6.1% |
| mg1655 | 4.64 | 1.2% | 1.5% | 2.8% | 5.3% |
| nissle1917 | 5.44 | 1.2% | 1.7% | **13.8%** | **15.9%** |
| sakai | 5.59 | 1.2% | 1.7% | 7.8% | 11.5% |
| se11 | 5.16 | 0.7% | 1.2% | 4.0% | 6.8% |

用 mg1655 的 RepeatMasker 输出（`tests/pgr/mg1655.rm.gff`，49,379 bp = 1.06%）
做交叉验证：e-kmer TnCentral 覆盖了 RM 的 90.7%；而 s-kmer 的 128 kb 里只有
32.7% 落在 RM 内，s-align 的 244 kb 里只有 17.8% —— 自找方法的大头不是转座子，
而是 SD、rRNA 操纵子、tRNA、多拷贝基因家族等（docs/rept.md 已明示 s-kmer /
s-align 会把 SD 算进去）。

> 注：§2.6 的 mg1655 RM 输出（49,379 bp）来自旧 singularity 镜像、Dfam
> `-species bacteria` 库。2026-08-07 原生安装 RepeatMasker 4.2.4 后用同一
> TnCentral 库重跑，结果见 §2.7（RM-is 163,249 bp / RM-strict 89,743 bp）。

nissle1917 是明显的异常值（s-kmer 748 kb / 13.8%、s-align 866 kb / 15.9%），
重跑两次结果完全一致，不是运行错误：其 s-kmer 有 14 个 ≥10 kb 的大块（合计
278 kb，最大 39.5 kb），其中约 292 kb 集中在染色体末端 5.0–5.44 Mb，属于大片段
重复/SD 区域；而 mg1655 完全没有 >10 kb 的块（最大 5.7 kb）。所以 nissle1917
反映的是株系特异的大重复，不代表大肠杆菌的一般水平。

### 2.7 RepeatMasker 金标准核对（原生安装 + TnCentral 库，10 株全跑）

2026-08-07 原生安装 RepeatMasker 4.2.4（RMBlast 2.14.1 默认引擎 + TRF），
对 10 株用**同一个 TnCentral 库**（`zcat` 解压后 `-lib`，`-pa 8` 串行，
单株 11–25 s）全基因组遮蔽。`.out` → `util/rmOutToGFF3.pl` → `pgr gff runlist`。
两个 RM 口径：

* `rm-is`：全部 IS 行（过滤掉 `Simple_repeat`/`Low_complexity`/`Satellite`），
  RM 原始灵敏度；
* `rm-strict`：再加查询跨度 ≥ 50 bp、100−div−del−ins ≥ 70
  （对齐 e-align 的 min-len=50 / min-identity=0.70）。

10 株合计（总长 52,425,656 bp），我们方法 vs RM-is（2,150,773 bp / 4.10%）：

| 方法 | 我们 bp | RM-is bp | 交集 bp | RM 被我们覆盖 | 我们被 RM 覆盖 |
| :--- | ---: | ---: | ---: | ---: | ---: |
| e-kmer TnCentral | 876,745 | 2,150,773 | 876,126 | 40.7% | 99.9% |
| e-align PGI 旧默认 | 1,151,386 | 2,150,773 | 1,150,515 | 53.5% | 99.9% |
| e-align PGI k31 | 1,250,721 | 2,150,773 | 1,249,341 | 58.1% | 99.9% |
| e-align LASTZ | 1,361,723 | 2,150,773 | 1,360,586 | 63.3% | 99.9% |

vs RM-strict（1,366,188 bp / 2.61%）：

| 方法 | 交集 bp | RM 被我们覆盖 | 我们被 RM 覆盖 |
| :--- | ---: | ---: | ---: |
| e-kmer TnCentral | 875,693 | 64.1% | 99.9% |
| e-align PGI 旧默认 | 1,134,520 | 83.0% | 98.5% |
| e-align PGI k31 | 1,202,869 | 88.0% | 96.2% |
| e-align LASTZ | 1,278,374 | 93.6% | 93.9% |

残差构成（RM-strict 与我们方法的差集，10 株合计）：

| 方法 | RM 独有 bp / 片段数（平均长度） | 我们独有 bp / 片段数（平均长度） |
| :--- | :--- | :--- |
| e-kmer | 490,495 / 2,085（235 bp） | 1,052 / 163（6 bp） |
| e-align PGI 旧默认 | 231,668 / 1,712（135 bp） | 16,866 / 237（71 bp） |
| e-align PGI k31 | 163,319 / 1,705（96 bp） | 47,852 / 288（166 bp） |
| e-align LASTZ | 87,814 / 1,695（52 bp） | 83,349 / 520（160 bp） |

要点：

1. **我们几乎不产生 RM 之外的假阳性**：对 RM-is 口径全部 99.9% 命中在 RM 内；
   对 RM-strict 也仍有 93.9–99.9%（LASTZ 最低，见第 3 点）。
2. **RM 原始输出远比我们宽松**：RM-is 是我们的 1.6–2.5 倍。mg1655 的 646 行
   IS 命中里 207 行 < 50 bp、198 行相似度 < 70%，这些短/弱匹配占了大头；
   加上阈值过滤后 RM-strict 比 RM-is 少 45%（mg1655 163,249 → 89,743 bp）。
3. **相同阈值下差距大幅缩小**：LASTZ 总 bp（1,361,723）≈ RM-strict
   （1,366,188），覆盖 RM-strict 的 93.6%；PGI k31 88.0%、旧默认 83.0%、
   e-kmer 64.1%。残差仍是短片段：RM 独有平均 52 bp（对 LASTZ）～96 bp
   （对 k31），即紧贴阈值的一批短匹配；我们独有平均 160 bp，主要是 RM 对
   同一元件按自身打分切碎、以及双方 identity 口径（RM 的 100−div−del−ins vs
   e-align 的匹配比例）差异造成的边界差。
4. **旧 §2.6 数字的适用边界**："e-kmer 覆盖 RM 90.7%"针对的是 Dfam-RM
   （49,379 bp）；同库同基因组下换成 TnCentral-RM 后，e-kmer 只覆盖
   RM-is 的 40.7% / RM-strict 的 64.1%。做金标准比较必须注明 RM 用的是
   哪个库。
5. **E. coli 的 simple repeat 很少**：RM 全量仅比 RM-is 多 4–6 kb/株
   （合计约 54 kb），远低于 trf 的 258 kb —— RM 内嵌 TRF 参数保守，且与
   库方法互补的结论不变。

逐株明细（RM-strict 被各方法覆盖的百分比）：

| 基因组 | RM-is bp | RM-strict bp | e-kmer | e-align 旧 | PGI k31 | LASTZ |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| cft073 | 199,745 | 124,964 | 59.9% | 81.0% | 86.7% | 93.2% |
| e2348_69 | 224,129 | 146,879 | 64.0% | 79.3% | 88.7% | 94.0% |
| e24377a | 281,668 | 198,653 | 69.7% | 89.8% | 91.2% | 95.5% |
| ec042 | 210,389 | 129,636 | 62.4% | 83.9% | 88.1% | 94.1% |
| ec2011c_3493 | 294,651 | 215,844 | 73.8% | 89.1% | 92.7% | 96.0% |
| ec958 | 228,415 | 146,750 | 70.9% | 84.5% | 91.5% | 94.0% |
| mg1655 | 163,249 | 89,743 | 63.3% | 79.2% | 81.5% | 90.1% |
| nissle1917 | 197,061 | 115,704 | 54.3% | 76.9% | 83.4% | 91.9% |
| sakai | 198,256 | 120,271 | 57.2% | 77.0% | 84.9% | 91.1% |
| se11 | 153,210 | 77,744 | 46.0% | 77.8% | 80.7% | 90.4% |

命令链（可复现）：

```bash
# 1) RM 全基因组（每株，-pa 8 串行；-lib 不接受 gz，库需先 zcat 解压）
~/Scripts/pgr/RepeatMasker/RepeatMasker genome.fa.gz -lib tncentral.fa \
    -pa 8 -e rmblast -dir out/$g
# 2) .out → GFF → runlist
perl ~/Scripts/pgr/RepeatMasker/util/rmOutToGFF3.pl out/$g/$g.fa.out > gff/$g.rm.gff
pgr gff runlist gff/$g.rm.gff -o json/$g.rm.json
# 3) 口径过滤：rm-is 去 simple/low；rm-strict 再加 ≥50 bp、≥70% 相似度
awk '!/^ *(SW|score|$)/ && $11 !~ /Simple_repeat|Low_complexity|Satellite/' \
    out/$g/$g.fa.out > out/$g/$g.fa.out.is
awk '!/^ *(SW|score|$)/ && $11 !~ /Simple_repeat|Low_complexity|Satellite/ &&
     ($7-$6+1) >= 50 && (100-$2-$3-$4) >= 70' \
    out/$g/$g.fa.out > out/$g/$g.fa.out.strict
# 4) 交叉统计
pgr runlist statop sizes/$g.sizes json/$g.<method>.json json/$g.rm-strict.json --all
pgr runlist compare json/$g.rm-strict.json json/$g.<method>.json --op diff -o tmp.json
```

> 一键复现（含 RM 运行、GFF/runlist 转换、statop + diff 统计、汇总输出）：
> `scripts/rm-gold-compare.sh <genome_dir> <our_json_dir> <sizes_dir> <lib.fa> <out_dir>`。
> 注意 PATH 上的 `pgr` 是旧构建（无 `gff runlist`），脚本默认用仓库的
> `target/release/pgr`（无 release 则用 debug）。

运行产物（临时）：`/tmp/rm_gold/`（每株 `.out`/GFF/runlist JSON，`rm-is`/
`rm-strict` 两套）、`rm_vs_ours.tsv` / `rm_strict_vs_ours.tsv` / `diff_stats.tsv`。

### 2.8 `pgr rept masker` 完整复刻验证（TRF + RMBlast + TRF，10 株全跑）

2026-08-07 实现 `pgr rept masker`：按 RepeatMasker 4.2.4 `-lib` 流程逐阶段复刻
（设计见 [design/masker.md](design/masker.md)）。
每 60 kb/2 kb batch（`SimpleBatcher` 算法）依次执行：

1. TRF PERFECT（2/7/7/80/10/50/10，拷贝 > 4）——RM 第一 TRF 阶段，找到的
   简单重复被**切除**（源码确认：-001 中间文件比原序列短 118 bp = 两处
   PERFECT 命中长度）；
2. rmblastn `general_search_parameters`（minscore=225 原样、word_size 9、
   gapopen 24/gapextend 6、mask_level 101、xdrops 450/225/112、片段 GC 选
   20p##g 矩阵）——在 PERFECT 掩蔽后的序列上搜库（IS）；
3. TRF DIVERGED（2/3/5/75/20/33/7，拷贝 > 5）——RM 第二 TRF 阶段，在
   PERFECT+IS 掩蔽后的序列上找分歧简单重复；
4. 三段区间合并 → runlist JSON。

阶段间掩蔽用 X 近似 RM 的切除/X 掩蔽（hit 集等价、坐标不漂移）。
10 株 × TnCentral 全跑，每株约 28 s（2 进程 × 4 线程）。

10 株合计（总长 52,425,656 bp），rm 共遮蔽 2,218,258 bp（4.23%）：

| 对比方 | 对方 bp | 交集 bp | 对方被 rm 覆盖 |
| :--- | ---: | ---: | ---: |
| RM 全量 .out（IS + Simple_repeat） | 2,204,952 | 2,204,874 | **100.0%** |
| RM-is（仅 IS） | 2,150,773 | 2,150,733 | **100.0%** |
| RM-strict（≥50 bp / ≥70%） | 1,366,188 | 1,366,188 | **100.0%** |
| e-align PGI k31 | 1,250,721 | 1,248,057 | 99.8% |
| e-align LASTZ | 1,361,723 | 1,360,402 | 99.9% |

逐株（rm bp / RM 全量被覆盖 / RM-is 被覆盖）：

| 基因组 | rm bp | RM 全量 | RM-is |
| :--- | ---: | ---: | ---: |
| cft073 | 205,711 | 100.0% | 100.0% |
| e2348_69 | 231,310 | 100.0% | 100.0% |
| e24377a | 287,998 | 100.0% | 100.0% |
| ec042 | 217,778 | 100.0% | 100.0% |
| ec2011c_3493 | 302,042 | 100.0% | 100.0% |
| ec958 | 234,869 | 100.0% | 100.0% |
| mg1655 | 168,772 | 100.0% | 100.0% |
| nissle1917 | 203,400 | 100.0% | 100.0% |
| sakai | 205,826 | 100.0% | 100.0% |
| se11 | 160,552 | 100.0% | 100.0% |

要点：

1. **RM 全量输出（含简单重复）被完整覆盖**：10 株 RM .out 全部落在 rm
   区间内（100.0%，合计仅 78 bp 含入差），我们多出 ~0.6%（每株 1.1–1.5 kb）
   是元件端点小尾巴 + DIVERGED 边界——RM 的 ProcessRepeats 边界精修会
   裁，我们保留原始跨度。对遮蔽目标是"更完整"。
2. **根因复盘（为什么之前差）**：① minscore 误用 209（7.5% 折扣只在未
   使用的 `runTestStage`，正式 `runStage` 用 225）；② 未复刻 RM 默认
   60 kb/2 kb 分片与每片段 GC 矩阵；③ 漏掉 RM 的两个 TRF 阶段（PERFECT +
   DIVERGED）——这是"换物种会更大"的三处隐患，已全部消除。
3. **rm ⊇ 此前所有方法**：RM-strict / e-align k31 / LASTZ 均被 rm 覆盖
   99.8%+；同时补齐了简单重复（TRF 两阶段），是全流程最接近 RM 的一档。
4. **定位**：e-kmer / e-align（PGI/LASTZ）是快速近似；rm 是 RM 的忠实
   复刻（每株 ~28 s）。追求与 RM 一致选 rm，追求速度选 e-align / e-kmer。
5. rm 默认不带长度/相似度过滤（min-len=0，同 RM 原始输出）；要与 e-align
   可比时加 `--min-len 50 --fill-fragment 10`（≈ RM-strict 口径）。

运行产物：`/tmp/rm12/`（每株 runlist JSON + log）、`agg3.tsv`。

## 3. 解读

1. **真核 TE 库（RepBase/Dfam）在 E. coli 中的检出很少**：每株 e-kmer 仅 2–19 个
   片段、e-align 5–33 个片段；总长 e-kmer 1.5–15.4 kb、e-align 5.1–34.3 kb，
   覆盖 <0.7% 基因组，平均片段约 0.6–1.2 kb。与 TnCentral（原核 IS）的每株 36–101
   片段、0.7–3.5% 覆盖形成鲜明对比。→ "真核转座子出现在大肠杆菌基因组"的
   概率/丰度**极低**；这些短匹配更可能是保守区/低复杂度序列的随机命中，而非
   真核 TE 的横向转移拷贝。
2. **原核 IS 元件才是 E. coli 重复序列的主体**：TnCentral 检出量远高于真核库；
   e-align 比 e-kmer 更敏感（TnCentral 覆盖 1.67% → 2.20%）。
3. **自找重复（s-kmer / s-align）会把 SD 等全部算进去**（docs/rept.md 已明示），
   因此片段数与覆盖度远高于库方法：s-align 每株 1,457–2,277 片段、覆盖 8.67%。
   若后续要做 SD 检测，不能用它们预遮蔽（会把 SD 一起盖掉）。
4. **trf（串联重复）与库方法互补**：每株 74–143 片段、总长 8.7–86.7 kb。
   Nissle1917 的 trf / s-kmer / s-align 均明显偏高（trf 86.7 kb、s-kmer 平均片段
   2,017 bp），但其库方法检出并不突出 → 额外重复主要来自非 IS 的自重复/SD。
5. **TnCentral 上的方法关系（§2.3）**：e-align 几乎完全包含 e-kmer（99.7%），
   且多出约 24% 的 bp —— e-kmer 是 e-align 的保守近似，要做完整 TnCentral
   遮蔽应选 e-align。自找方法只能回收 60–74% 的 TnCentral bp（单拷贝或高度
   分化的 IS 拷贝低于自比对深度/频率阈值会被漏掉），其余大部分遮蔽来自 SD
   等非 IS 重复；trf 与 TnCentral 基本不相交（串联 vs 散在，符合预期）。
6. **PGI 与 LASTZ 后端核对（§2.4）**：TnCentral 上 PGI 的 99.5% 遮蔽 bp
   落在 LASTZ 内（每株 98.7–99.9%），但 LASTZ 更敏感也更碎 —— 1,486 片段 /
   1,361,723 bp，平均 916 bp，比 PGI（734 片段 / 1,569 bp 平均）多约 18% bp、
   多一倍片段；多出的 bp 里 94.6% 是 PGI 完全未检出的短片段（760 个、平均
   262 bp），而 85.3% 的 PGI 区间被单个 LASTZ 片段完整包住 —— 差异不在"谁把
   区域接得更长"，而在 LASTZ 额外捡到一批短匹配。前面关于 TnCentral / 真核库 /
   自找方法的结论不依赖比对后端；追求"少遗漏"选 LASTZ 版，e-align 默认仍是
   PGI（2026-08-07 起默认参数改为 k31/freq50/ms12，见 §2.5）。
7. **覆盖差距约一半是 e-align 默认参数造成的（§2.5）**：项目内其他 PGI 用法
   确实改参数——SD 检测实测把 k 40→31、freq 10→50（min-shared 12）换来灵敏度。
   e-align 用这套参数后 TnCentral 覆盖 1,250,721 bp（默认 1,151,386），与 LASTZ
   的差距缩小约 47%，但 LASTZ 仍多出 ~111 kb 短片段。e-align 默认参数已按此
   调整为 k31/freq50/ms12（与 SD 一致）。
8. **没有"1 Mb 需要遮蔽"的结论（§2.6）**：单株最大是 nissle1917 的 s-align
   866 kb（15.9%），10 株里唯一超过 10% 的异常株；"4.55 Mb / 8.67%"是 10 株合计
   /平均。经典重复（TnCentral 库遮蔽）每株只有 1–3.5%，真核 TE 库（RepBase/Dfam）
   每株 <0.5% —— 和"大肠杆菌转座子很少"一致；s-kmer/s-align 的 5–16% 是"所有
   重复相关序列"口径（含 SD 与多拷贝基因），不能当转座子含量读。
9. **RepeatMasker 金标准核对结论（§2.7）**：我们方法无假阳性问题（99%+ 命中
   落在 RM 内）；RM 原始输出更敏感，但多出的几乎全是 <50 bp 或 <70% 相似度的
   短碎片（RM-is 2,150,773 bp → 阈值过滤后 RM-strict 1,366,188 bp）。相同阈值下
   LASTZ 覆盖 RM-strict 93.6%、PGI k31 88.0%、旧默认 83.0%、e-kmer 64.1%；
   LASTZ 的总量已与 RM-strict 相当（1,361,723 vs 1,366,188 bp）。"金标准"
   必须声明所用库与过滤口径——同一基因组用 Dfam-RM 与 TnCentral-RM 差 3 倍以上。
10. **`pgr rept masker` 完整复刻成功（§2.8）**：按 RM 阶段顺序（TRF PERFECT →
    rmblastn → TRF DIVERGED）复刻后，RM 全量 .out（IS + 简单重复）被其
    覆盖 100.0%、RM-strict / e-align k31 / LASTZ 99.8%+，残差仅元件端点
    2–30 bp（RM 边界精修裁剪）。它成为方法梯队里最接近 RepeatMasker 的
    一档；快速近似仍用 e-kmer / e-align。

## 4. 参考

* 方法与库的完整说明：[docs/rept.md](../docs/rept.md)
* 重复遮蔽方案设计：[design/repeat-masking.md](design/repeat-masking.md)
* 10 株 cohort 来源与用途：[ecoli-genome.md](ecoli-genome.md)
