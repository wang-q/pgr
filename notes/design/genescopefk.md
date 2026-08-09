# genescopefk.R（GenomeScope 2.0）对齐记录

> 2026-08 整理。目标：让 `pgr kmer gsize --model` 与 anchr 用的
> `genescopefk.R`（GenomeScope 2.0）处理一致，输出可被
> `anchr templates/2_fastk.tera.sh`（`grep '^kmercov' model.txt` +
> `sed '1,6 d' summary.txt`）直接消费。参考源码：仓库内
> `minpack.lm-master/`（minpack.lm 1.2-4 源码）与用户下载的 R 4.4.2
> 源码（`src/nmath/`）。关联：[fastk.md](../references/fastk.md)、
> [merqury-fk.md](../references/merqury-fk.md)。

## 1. 根本 bug：lmdif 移植的列主序索引反转

`src/libs/lm.rs` 里 `qrsolv`/`lmpar` 有 **4 处** `r[i*n+j]` 与
`r[j*n+i]` 写反（R 矩阵按列主序存，`r[j*n+i]` = 第 j 列第 i 行 =
Fortran 的 `r(i,j)`）：

| 位置 | 错误写法 | 正确写法（= Fortran `r(i,j)` 语义） |
|---|---|---|
| qrsolv 初始拷贝 | `r[i*n+j] = r[j*n+i]`（下→上） | `r[j*n+i] = r[i*n+j]`（上→下） |
| qrsolv Givens i 循环 | `r[i*n+k]`（写进上三角） | `r[k*n+i]`（写进下三角） |
| qrsolv 回代 | `r[i*n+j]`（读上三角） | `r[j*n+i]`（读下三角） |
| lmpar 牛顿修正 | `r[i*n+j]`（读上三角） | `r[j*n+i]`（读下三角） |

后果：R 矩阵的上三角被逐步清零/污染，LM 轨迹完全错误——pgr 收敛到
length≈852（RSS 19.0M），而 R 收敛到 988（RSS 15.7M）。这就是最初
"Model Fit 62.13% vs 64.37%"、错误率 49% 等所有数字对不上的根源。

注意：早期简单的单测（J=I 的二次型）测不出这个 bug——R=I 时上/下三角
都是 0，索引反了也无差异。真正暴露问题的是对照 Fortran 源码的逐值差分。

## 2. 修正清单

### 2.1 `qrsolv` / `lmpar` 索引（见上）

修完后用 minpack.lm Fortran 的精确输入做单测，`lmpar` 输出 par 与
`x` 逐位一致（见 `src/libs/lm.rs` 的 `lmpar_matches_minpack_lm_reference`
与 `qrsolv_matches_minpack_lm_reference`）。

### 2.2 `qrfac`：对齐 minpack.lm 的实际实现

**注意：minpack.lm 的 qrfac.f 与经典 netlib 版本不同**：

- 主元前**没有**循环顶部的列范数更新（经典版有）；
- 列范数重算条件是指数 **2**（`p05*(rdiag/wa)**2 <= epsmch`），不是经典的 3；
- 列范数工作数组是独立的 `wa`，`acnorm` 保持原始列序（lmdif 后面要用它做
  diag 缩放和 gnorm）。原来的实现把 `acnorm` 当工作数组用（随主元交换、
  被更新），导致 diag/gnorm 拿到错误的列范数。

### 2.3 `enorm`：缩放的 minpack 版本

minpack 的 enorm 用 rdwarf/rgiant 分段的三和累积防溢出/下溢，不是朴素
`sqrt(sum(x²))`。已按 `enorm.f` 移植。

### 2.4 qtf 循环

`fjac(j,j)==0` 时 Fortran 仍执行 `fjac(j,j)=wa1(j); qtf(j)=wa4(j)`，
原实现 `continue` 跳过了这两行。

### 2.5 错误率 / score / summary / write_outputs（早前已修，见 git 历史）

- 错误率：R 公式 `1-(1-total_error_kmers/total_kmers)^(1/k)`（旧实现把
  "低覆盖 k-mer 比例"当错误率，得到 49%）；
- `score_model`：去掉直方图尾位（R 的 `kmer_hist_orig`）；
- `err_cut`：R 语义 `tail(which(x <= kcovfloor), 1)`；
- `summary.txt`：5 行头 + 空行、min 列 NA、千分位、signif 百分比、
  `Read Error Rate` 行；p=2 的 Homozygous/Heterozygous 列用 r1 的 2-SE 区间；
- 失败时只写头、不写 model.txt。

## 3. R dnbinom 移植：`src/libs/kmer/nbinom.rs`

R 的 `dnbinom(x, size, mu)` 走 `dnbinom_mu`，内部用 Loader 的无消减算法
（bd0/stirlerr），比朴素 `lgamma` 公式在病态问题上稳得多。拟合问题的
Hessian 条件数高达 6e15，残差 1e-12 级的舍入差异会被放大成完全不同的
局部最优——所以必须逐位复刻 R 的数值。

移植内容（全部来自 R 4.4.2 `src/nmath/`）：

- `stirlerr`（stirlerr.c：半整数表 + 分段级数）；
- `bd0`（bd0.c：Taylor 级数）；
- `logcf` / `log1pmx` / `lgamma1p`（pgamma.c）；
- `dbinom_raw`（dbinom.c：`lc - 0.5*lf`，`lf = ln(2π)+ln x+ln1p(-x/n)`）；
- `ebd0` + `BD0_SCALE` 表（bd0.c：128 项 × 4 个 f32 的对数表，Rust 不认
  十六进制浮点，需转十进制）；
- `dpois_raw`（dpois.c：size=Inf 时的 Poisson 极限）；
- `dnbinom_mu` 的完整分支（x==0、x<1e-10*size 的 MM 公式、主路径）。

**坑**：提取 `BD0_SCALE` 时正则要带上符号——表里大量 `-0x...` 项，漏掉
负号会让 `dpois_raw` 差 1e-7 量级（ebd0 的 yh+yl 对不上真值）。

验证：`dnbinom`/`dpois_raw` 对照 R 逐值（8.8e-17 内），单测见
`nbinom.rs`。

## 4. 验证方法（差分对照）

1. **gfortran 编译 minpack.lm-master 真实源码**（lmdif.f/qrfac.f/lmpar.f/
   qrsolv.f/fdjac2.f/enorm.f/dpmpar.f + 移植的 bd0 dnbinom），跑同一份
   hist.tsv：第一步接受 (0.0044, 54.34, 0, 942.4)，与 R 的轨迹一致，
   收敛到 length≈1019.9——证明算法+残差到位后轨迹可复现；
2. **逐值差分**：对比 qrfac 后的 rdiag/acnorm/qtf/delta、lmpar 的
   par/x/sdiag，定位到 qrsolv 的索引反转；
3. **R 对照**：`dnbinom` 逐 x 对比、`dpois_raw` 对照 R 输出；
4. **单元测试**：`lmpar_matches_minpack_lm_reference`、
   `qrsolv_matches_minpack_lm_reference`（用 Fortran 的精确输入锁定参考值）、
   `dnbinom_matches_r_values`、`dpois_raw_matches_r`。

## 5. 最终对齐状态（hist.tsv，k=17，p=1，单起点）

| 项 | pgr（修复后） | R 参考 | gfortran 参照 |
|---|---|---|---|
| d | 0 | 0 | 0 |
| kmercov | 55.76 | 55.73 | 55.78 |
| bias | 0 | 0 | 0 |
| length 参数 | 1017 | 988 | 1019.9 |
| Genome Haploid Length | 899 bp | 899 bp | — |
| Model Fit | 63.98% | 64.37% | — |
| Read Error Rate | 0.0331% | 0.0348% | — |

summary.txt 除 Model Fit / Error Rate 两行外与 R 逐字节一致；模型结构
（d=0、bias=0）完全一致。剩余差异是同一盆地内 Rust 与 f2c C 的语言级
浮点路径差异（1017 vs 988 vs 1019.9），无法也不必要逐位消除。

**p=2 端到端对照**（同 hist.tsv，2026-08-10）：pgr d=0/r1=0.001197/
kcov=27.91/bias=0/length=1041.8 vs R d=0/r1=0.001252/kcov=27.86/
bias=0/length=1004——结构与 p=1 相同的对齐程度（同盆地，summary
Haploid 903 vs 904、Repeat/Unique 一致，Model Fit 64.71% vs 65.25%）。

**anchr 2_fastk 消费验证**（2026-08-10）：`grep '^kmercov' model.txt |
tr -s ' ' '\t' | cut -f 2` → COV=55.8；`summary.txt` 经 `sed '1,6 d'` +
表格式处理的输出与 R 过同一管线**逐行结构一致**（仅 Model Fit /
Error Rate 数值不同）。

## 6. 注意事项 / 教训

- **rsstrace 陷阱**：minpack.lm 的 `rsstrace[0]` 会被首轮 fdjac2 的扰动
  评估污染（fdjac2 用同一个 niter 槽位写 RSS），所以 R 打印的 "It. 0 RSS"
  不是起点的真实 RSS——别拿它当残差对照基准；
- **Fortran 打印**：固定格式（fixed-form）写语句超过 72 列会被截断，调试
  打印要拆短或用 list-directed；
- **合成数据测试**：`command_kmer_gsize_model_fit` 的 genome_size 区间
  (300, 3000) 是按坏 LM 校准的，修复后真实拟合（同盆地）给 4647，区间已
  放宽到 (500, 10000)；
- **单起点 vs 多起点**：早期为了绕开坏 LM 加了多初值（d×length 网格），
  LM 修好后恢复 R 的单起点配置（d=0.10, length=est/p），与参考的处理一致。
