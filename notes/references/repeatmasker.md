# RepeatMasker 安装与自定义库（TnCentral）使用记录

## 结论

RepeatMasker 本身**不绑定任何数据库**。开箱即可用自定义 FASTA 库：

```bash
RepeatMasker genome.fa -lib my_library.fa
```

TnCentral 的 FASTA（`~/data/repeats/tncentral.fa.gz`）完全走这条路，**不需要也不能"装进"程序**。
只有按物种自动选库的 `-species` 模式才需要 FamDB（Dfam/RepBase 的 H5 格式），TnCentral 不是这种格式。

## 本地安装状态（2026-08-07）

- 源码：`/home/wangq/Scripts/pgr/RepeatMasker/`，版本 4.2.4（官方最新）。
- 搜索引擎：RMBlast 2.14.1（CBP 编译版，`~/.cbp/bin`，要求 glibc ≤2.16，
  CentOS 7 可跑），经合并目录
  `/home/wangq/Scripts/pgr/RepeatMasker/rmblast-cbp-bin` 配置为默认引擎。
- TRF：`/home/wangq/.cbp/bin/trf`（系统已有）。
- FamDB：未配置（因此只有 `-lib` 模式可用；`-species` 不可用）。
- 已通过 `perl ./configure` 完成配置（RMBLAST_DIR 指向 `rmblast-cbp-bin`），
  `./RepeatMasker -h` 可运行。RMBlast 的 tar 包仍在 `~/Downloads/`，/tmp 里
  的解压副本可删。

> **为什么是合并目录**：configure 校验 RMBLAST_DIR 需同时存在
> rmblastn / dustmasker / makeblastdb / blastdbcmd / blastdb_aliastool /
> blastn 六个可执行文件（`RepeatMaskerConfig::validateParam` 逐个检查
> `-x`），而 CBP 只装了前两个。合并目录里 rmblastn+makeblastdb 软链 CBP 版
> （真正使用的引擎），其余四个软链官方包、仅为 configure 校验占位——
> RepeatMasker 4.2.4 运行时只调 rmblastn 和 makeblastdb（源码核实：
> dustmasker/blastdbcmd/blastdb_aliastool/blastn 除校验名单外无任何引用）。
> **注意**：四个占位软链指向官方预编译包（glibc ≥2.29），在 CentOS 7 上
> 本身跑不起来，只是不会被调用；若想彻底干净，可把这四个换成同年代老构建
> （如 blast+ 2.2.28）或 bioconda 包里的对应二进制。搬到 CentOS 7 时软链
> 目标路径需保持一致或改硬拷贝。

### 重新 configure（若以后目录再移动）

```bash
cd /home/wangq/Scripts/pgr/RepeatMasker
perl ./configure -perlbin "$(which perl)" \
  -trf_prgm /home/wangq/.cbp/bin/trf \
  -rmblast_dir /home/wangq/Scripts/pgr/RepeatMasker/rmblast-cbp-bin \
  -default_search_engine rmblast
```

configure 期间回答 "Configure FamDB now?" 为 `n`（当前不需要物种库）。

## 冒烟测试（TnCentral 库）

MG1655 前 100 kb（`tests/genome/mg1655.fa.gz`），`-pa 8`：

```bash
zcat ~/data/repeats/tncentral.fa.gz > /tmp/rmtest/tncentral.fa
RepeatMasker mg1655_chunk.fa -lib /tmp/rmtest/tncentral.fa \
  -pa 8 -e rmblast -dir /tmp/rmtest/out2
```

结果：3.12% 被遮蔽；检出 IS621、IS186B、IS1A、ISPpu12、ISEc39、
ISSoEn2、Tn7243 等真实 IS 序列，与 MG1655 已知内容吻合。

## 两个实测坑

1. **`-lib` 不接受 gzip**：直接把 `tncentral.fa.gz` 传给 `-lib` 会失败，
   makeblastdb 报 "Input doesn't start with a defline"。
   必须先 `zcat` 解压成 `.fa`。
2. **TnCentral 源库有 24/6093 条记录格式瑕疵**：部分序列行开头黏着
   accession 前缀（如 `In1223` 的序列以 `NX784502...` 开头，少数以
   `_PAJ...` 开头）。RepeatMasker 能容忍（按非法字符处理），
   但正式做金标准比对前建议把这些前缀清掉。
3. **RM 会把 .fa.gz 输入解压到文件旁**：RepeatMasker 对 `.fa.gz` 输入
   自动 `gunzip -c file.gz > file`（写到输入同目录，RepeatMasker:748-755）。
   **不要在仓库内直接对 .gz 跑 RM**，否则会在源码目录留下一份解压副本
  （2026-08-07 实测：`tests/genome/mg1655.fa` 因此反复出现）。做法：
  先 zcat 解压到 /tmp 再喂 RM，或用软链；若文件旁已存在同名解压版，
  RM 会直接报错退出。

## RepeatMasker 源码算法流水线

> 以下内容经核对 `RepeatMasker` 主脚本（7426 行 Perl）与 `NCBIBlastSearchEngine.pm`/
> `TRF.pm` 源码整理，用于 `pgr rept masker` 复刻对照（配合 [design/masker.md](../design/masker.md)）。

### 整体流程

```
输入 FASTA（.gz 先 gunzip 到输入同目录，见「两个实测坑」第 3 条）
  │  按 batch 分片：SimpleBatcher，fragmentSize=60000、overlapLen=2000（RepeatMasker:629,638）
  │  多进程并行（fork，-pa 控制），每 batch 独立跑「搜索阶段」
  ▼
每 batch 的搜索阶段（runSearchStages，rmblast 引擎）：
  1. TRF 阶段：识别 Simple Repeats（PERFECT 参数，见下）
  2. rmblast 阶段：对 custom/species 库做高复杂度搜索（runStage）
  3. 若干后续 stage（SINE/retro/cut 等）——注意：这是 **物种库（-species）模式**才有的
     多阶段（按类分库逐步搜索）；`-lib` 自定义库模式在 PERFECT/DIVERGED 两次 TRF 后
     只有 `general_search_parameters`（stage 001）一个 rmblast stage（`RepeatMasker:3612-3656`）
  → 产出该 batch 的批注 *.cat
  ▼
合并所有 batch 的 *.cat → $file.cat（带头部，见下）
  ▼
ProcessRepeats（外部脚本）：读 .cat，做 final 处理（合并重叠、算 divergence/
  Kimura、subelement 划分、IS 剪切），输出 .out / .masked / .tbl / .gff
```

- **batch 分片**：`fragmentSize=60000`、`overlapLen=2000`（`RepeatMasker:629,638`，`-frag` 可改）。
  batch 间 2000 bp 重叠，结束后 `adjustFragmentPositions` 修正批注坐标。
- **GC 相关矩阵**：搜索矩阵名实为 `{div}p{GC}g.matrix`（如 `20p43g.matrix`），`div` 来自
  `searchParams->{matrix}` 的 `\d+` 前缀、`GC` 由 `chooseMatrices` 按实测 GC 分档
  （35/37/39/…/51/53g，`RepeatMasker:4229-4266`）。GC 背景默认 43%，仅当
  `-gccalc` 或"单序列 && 长度>2000 bp"时用 batch 实测平均 GC（`RepeatMasker:2855-2881`，
  上限 `max_matrix_gc`）。默认 `general_search_parameters`（-lib 走这个 recipe）为
  `minscore=225、minmatch=[8,9,11,13]、matrix=20p##g.matrix、gap_init=-30、bandwidth=14、
  masklevel=90、filterContained=0`（`RepeatMasker:2122-2138`），可供 pgr 复刻时对齐默认阈值。

### TRF 两套参数（TRF.pm + runTRFStage）

RepeatMasker 把"简单/低复杂度重复"交给 TRF 分两阶段：

| 阶段(searchParams) | 用途 | match/mismatch/delta | pm/pi | minscore | maxperiod | minCopyNumber |
|------|------|----------------------|-------|----------|-----------|---------------|
| `PERFECT` | 年轻串联重复（young） | 2 / 7 / 7 | 80 / 10 | 50 | 10 | 4 |
| `DIVERGED` | 古老串联重复（old） | 2 / 3 / 5 | 75 / 20 | 33 | 7 | 5 |

`runTRFStage` 用 `$searchParams eq "PERFECT"` 区分两套（`RepeatMasker:2665-2698`）；
注意第二阶段实际传入的字符串是 **`"DIVERGED"`**（`RepeatMasker:3356,4086`），并非表里的
"OLD"——那是用途描述，源码只在 `eq "PERFECT"` 处分支，其余一律走"古老串联重复"分支。
两阶段编号/chooseClass：PERFECT → stage 251、chooseClass `"simple"`（`RepeatMasker:3602-3603`），
DIVERGED → stage 252、chooseClass `"masking"`（`RepeatMasker:3362-3363`）。
`lambda`（0.41/0.32）与 `mu` 数组用于把 TRF 原始分折算成与 blast 可比的 bitScore
（`RepeatMasker:2679-2697`）。`pgr rept masker` 的 TRF 两阶段（young/old）即对应这两套参数。

**TRF 结果的后处理链**（`runTRFStage`，`RepeatMasker:2729-2785`）：先过
`copyNumber > minCopyNumber` 阈值，再用 crossmatch 的 `simple1.matrix` 对每个 hit 重打分
（gapOpen=-30、gapExt=-15、xDrop=500），simple1.matrix 分阈值 20，再经 lambda/mu 折算
bitScore；仅 PERFECT 阶段额外做 `maskLevelFilter(1)`（因它要 excision，不允许重叠，
`RepeatMasker:2780-2785`）。这是 pgr 复刻 TRF 时可直接对齐的阈值链。

### rmblast 搜索与 outfmt

- **outfmt**：`-outfmt="6 score perc_sub perc_query_gap perc_db_gap qseqid qstart
  qend qlen sstrand sseqid sstart send slen kdiv cpg_kdiv transi transv cpg_sites"`
  （`NCBIBlastSearchEngine.pm:611-613`，默认无 alignment 时为 18 列；带 `-a` 生成 alignment
  时 611 行在末尾追加 `qseq sseq` 成 20 列）。2.13+ 的 kdiv/cpg 列在此确实被 RM 使用——
  但仅用于计算 divergence，`pgr rept` 的 `parse_tab_row` 只取 qseqid/qstart/qend 三列做
  区间，无需这些列。
- **过滤链**（`runStage`，`RepeatMasker:2918-2952`）：
  1. `search`：调 rmblastn + makeblastdb（对自定义库 `processCustomLib` 先建库）；
  2. `filterContainedResults`（可选，滤掉被更长 hit 包含的）；
  3. `maskLevelFilter(value => $masklevel)`：按 masklevel 阈值过滤；
  4. `filterResults(...)`：按 `$filterType`（"simple"/"masking"/...）与 divergence 做最终过滤。
  5. 各 stage 的搜索参数来自 `getSearchRecipes`（minscore/minmatch/bandwidth 等）。
- **`-lib` 自定义库路径**：`processCustomLib` 先 `makeblastdb` 建库（`RepeatMasker:6551`），
  因此 `-lib` 必须传未压缩的 `.fa`（见「两个实测坑」第 1 条）。

### .cat 批注格式（供 ProcessRepeats）

`$file.cat` 头部为若干 `##` 行（`## RepeatMasker version …`、`## RM Library: …`、
`## Total Sequences/Length/NonMask(排除>20bp N/X 段)/NonSub(排除所有非 ACGT 碱基)`，
batch 有重叠时还有 `## Batch Overlap Boundaries`，`RepeatMasker:1228-1250`）+ 结束标记
`## RAW Annotations:`，之后按 batch 追加原始批注（每行一个 hit）。`ProcessRepeats` 把 raw
批注做合并/打分/分级（`parseCATFile`→`processSequence`，片段合并 `cycleReJoin`，各类别专一
处理器 `preProcessLINE/preProcessDNATransp/preProcessLTR/scoreSINEPair/scoreLTRPair` 等），
最终 `divergence` 在 `filterResults` 处折算为 `100 * pctDiverge / (100 - pctInsert)`
（`RepeatMasker:4588-4595`），由 `-div` 阈值过滤。

### 对 pgr rept 的启示

- `pgr rept masker` 复刻的是 **TRF 两阶段 + rmblast 区间化** 的最小闭环，不必复刻 batch
  分片/ProcessRepeats/GC 矩阵那套工程外壳；但"TRF young/old 两套参数 + minCopyNumber 阈值"
  与"rmblast 区间去重 + divergence 过滤链"的顺序值得对齐，保证与 RM 金标准可比。
- `.cat` 的原始批注格式（按 hit 逐行）与 pgr 的 PAF/区间中间表示可以互相转换，便于对拍。
- **rmblast 版本特性开关**（`NCBIBlastSearchEngine.pm:setPathToEngine`）：2.13+ 启用
  `hasTabFormat`（18 列 outfmt）+ `hasQueryThreading`，2.14.1+ 追加 `hasDBSoftMasking`
  （允许库侧软屏蔽）。`pgr rept masker` 若直接消费 rmblastn，应至少要求 2.13+ 的 tab
  输出；2.14.1+ 的 DB soft masking 对处理低复杂度库有参考价值。
- **ProcessRepeats 的"分裂片段合并 + 家族亲缘评分"是 pgr 若要对拍 RM 的 .out/.tbl/.gff
  才需要的重头戏**：`cycleReJoin` 把被 batch/剪切拆碎的同一元素碎片按 (Score,PctSub,
  PctDel,PctIns) 四元组重连，`preProcessLINE/preProcessDNATransp/preProcessLTR` 与
  `scoreSINEPair/scoreLTRPair/scoreLINEPair` 做各类别的双端元件拼接判定，
  `joinDNATransposonFragments` 处理 DNA 转座子断裂。仅做"区间去重"无法复现 RM 的 final
  注释语义；若 pgr 只做 masker 层面（区间输出）则不需要，但做 `.out` 金标准对拍时需要。
- **makeblastdb 建库在运行时只发生在 processCustomLib（-lib）与 createLib（库合成）两处**，
  NCBI 引擎实际只调 rmblastn + makeblastdb（`RepeatMasker:6549,7332`），其余四个二进制
  （dustmasker/blastdbcmd/blastdb_aliastool/blastn）仅存于 configure 校验名单
  （`RepeatMaskerConfig.pm:119-124` 的 `expected_binaries`），pgr 复刻引擎路径时无需理会。

## CentOS 7（glibc 2.17）部署兼容性（2026-08-07）

**问题**：NCBI 官方 `rmblast-2.14.1+-x64-linux-GLIBC_2.31.tar.gz`
预编译包要求的最高 glibc 符号为 **GLIBC_2.29**（本机 readelf 实测），
CentOS 7 只有 glibc 2.17，直接跑不起来。

**我们代码的真实依赖**（比之前判断宽松）：
- `parse_tab_row` 只消费 rmblastn tab 输出的 qseqid/qstart/qend 三列；
  18 列 outfmt 里 2.13+ 才新增的 kdiv/cpg 等列**没有被使用**。
- RepeatMasker 4.2.4 自身对 rmblastn <2.13 也只是退回 legacy 解析
  （NCBIBlastSearchEngine.pm setPathToEngine 里的特性开关），并非硬性拒绝。

**可选方案**（按推荐顺序）：
1. **CBP 安装的 rmblast 2.14.1（`~/.cbp/bin/rmblastn` + `makeblastdb`）**：
   2026-08-07 实测，版本 2.14.1、要求的最高 glibc 符号仅 **GLIBC_2.16**
   （官方预编译包是 2.29），CentOS 7（2.17）可直接运行。用它跑 60 kb MG1655
   片段 × TnCentral：63 条原始命中 / 25 条去重区间，与官方 2.14.1 结果
   完全一致；18 列 outfmt 与 v5 库格式均为正式行为，代码零改动。
2. **bioconda rmblast 2.14.1**：linux-64 conda 包按 glibc 2.17 兼容构建，
   CentOS 7 可直接运行，版本与本地验证完全一致、零结果差异。
   服务器装 micromamba（单静态二进制，无需 root）：
   `micromamba create -p ~/rmblast-env -c conda-forge -c bioconda rmblast=2.14.1`，
   然后 `pgr rept masker ... --rmblast-dir ~/rmblast-env/bin`。
3. **CentOS 7 源码编译 2.14.1**：需要 devtoolset-8+（C++14）与 boost、
   zlib/bzip2 开发包，编译耗时长，只作无 conda 时的备选。

**pgr 本体**：glibc 兼容性由项目发布流程保证（`.github/workflows/publish.yml`
用 `cargo zigbuild` 交叉到 glibc 2.17），此处不赘述。TRF 4.09 为静态二进制。

## 相关笔记

* [design/masker.md](../design/masker.md)：`pgr rept masker`
  实现设计（参数表、TRF 两阶段、验证结果）
* [design/repeat-masking.md](../design/repeat-masking.md)：重复标记总体方案与
  RepeatMasker 源码梳理（附录 A）
* [ecoli-repeats.md](../ecoli-repeats.md) §2.7/§2.8：RepeatMasker 金标准核对与
  masker 复刻对拍

## 参考

- 官方安装页：https://www.repeatmasker.org/RepeatMasker/（依赖、configure 流程、`-lib` 声明）
- RMBlast 下载：https://www.repeatmasker.org/rmblast/
- GitHub README："You can use it immediately with a custom library (`-lib mylib.fa`)"
- GitHub issue #289：额外库（Dfam partition 等）均通过 `-lib` 传入
