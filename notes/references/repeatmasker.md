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
  │  按 batch 分片：SimpleBatcher，fragmentSize=60000、overlapLen=2000（RepeatMasker:628,637）
  │  多进程并行（fork，-pa 控制），每 batch 独立跑「搜索阶段」
  ▼
每 batch 的搜索阶段（runSearchStages，rmblast 引擎）按固定顺序执行：
  1. TRF PERFECT 阶段（stage 251，chooseClass=simple）：识别**年轻** Simple Repeats
     （`-nolow` 或 `-alu` 时跳过；`RepeatMasker:3594-3605`）
  2. rmblast 高复杂度搜索阶段（runStage）：
     - **`-lib` 自定义库模式**：只有一个 stage —— `general_search_parameters`
       （stage 001，chooseClass=masking，`RepeatMasker:3619-3656`），随后若
       `seqDB->getSubtLength() < 15` 则提前 `last`。
     - **`-species` 物种库模式**：按类分库的**多阶段**逐步搜索（sinecutlib/shortcutlib/
       cutlib/shortlib/longlib/retrolib/at.lib 等，stage 35x/40x/45x/50x/60x/70x）。
  3. TRF DIVERGED 阶段（stage 252，chooseClass=masking）：识别**古老** Simple Repeats
     （`-nolow` 时跳过；`RepeatMasker:4082-4097`）——**位置在所有 rmblast 阶段之后**。
  → 产出该 batch 的批注 *.cat
  ▼
合并所有 batch 的 *.cat → $file.cat（带头部，见下）
  ▼
ProcessRepeats（外部脚本）：读 .cat，做 final 处理（合并重叠、算 divergence/
  Kimura、subelement 划分、IS 剪切），输出 .out / .masked / .tbl / .gff
```

- **batch 分片**：`fragmentSize=60000`、`overlapLen=2000`（`RepeatMasker:628,637`，`-frag` 可改，
  且 `-frag` 不得 < 2×overlapLen，`RepeatMasker:642-651`）。`SimpleBatcher::_packBatches` 的
  打包规则（`SimpleBatcher.pm:592-800`）：短序列累积直到 `batchLenCtr + seqLength > fragmentLen*1.25`
  才封存本 batch（故 batch 可略超 60000）；单条序列长度 > fragmentLen 时被切分，切分数 `divisor`
  为满足 `(len+(divisor-1)*overlapLen)/divisor <= fragmentLen` 的最小整数，每片实际长度
  `size=int((len+(divisor-1)*overlapLen)/divisor)`，相邻片重叠 2000 bp；含切片的 batch
  `completeSeqs=0`（`isBatchFragmented` 为真），其批注坐标后续由 `adjustFragmentPositions`
  修正回原序列（`RepeatMasker:1399`），并用 `overlapMiddle=overlapLen/2=1000` 做边界换算。
- **fork 并行调度**（`RepeatMasker:834-1193`）：主进程维护 `%batchStatus`/`%children`，
  `numberChildren = min(-pa, batchCount)`，`JOBLOOP` 里 `wait()` 回收子进程。子进程
  产出一个 `_batch-N.masked`（`SimpleBatcher::writeBatchFile`）+ 一个 `_batch-N.cat`；父进程
  用 `nextBatchToConcatenate` 按 batch 序号**严格有序**地追加进 `$file.cat.all`（并行完成
  顺序不定，未就绪就 `last` 等下一轮，`RepeatMasker:993-1016`）。batch 失败自动重试
  `retryLimit=2` 次，连续 `badForkMax=20` 次坏 fork 才退出（`RepeatMasker:837-839,921-932`）；
  单序列 batch 完成后若 `isBatchFragmented` 由子进程就地跑 `adjustFragmentPositions`
  （`RepeatMasker:1153-1160`）。
- **.cat 输出与合并**：`totseqlen > 10Mb` 时整个 `.cat` 用 gzip 输出为 `$file.cat.gz`
  （`RepeatMasker:831,1222-1226`），否则为 `$file.cat`；头部 `##` 行与 `## RAW Annotations:`
  结束标记由父进程写（`RepeatMasker:1216-1251`），各 batch 的 raw 批注追加其后。
- **GC 相关矩阵**：搜索矩阵名实为 `{div}p{GC}g.matrix`（如 `20p43g.matrix`），`div` 来自
  `searchParams->{matrix}` 的 `\d+` 前缀、`GC` 由 `chooseMatrices` 按实测 GC 分档
  （35/37/39/…/51/53g，`RepeatMasker:4229-4266`）。GC 背景默认 43%，仅当
  `-gccalc` 或"单序列 && 长度>2000 bp"时用 batch 实测平均 GC（`RepeatMasker:2858-2881`，
  上限 `max_matrix_gc`）。默认 `general_search_parameters`（-lib 走这个 recipe）为
  `minscore=225、minmatch=[8,9,11,13]、matrix=20p##g.matrix、gap_init=-30、
  ins_gap_ext=-6、del_gap_ext=-5、bandwidth=14、masklevel=90、filterContained=0、
  chooseClass=masking、excise=0`（`RepeatMasker:2122-2138`），可供 pgr 复刻时对齐默认阈值。

### TRF 两套参数（TRF.pm + runTRFStage）

RepeatMasker 把"简单/低复杂度重复"交给 TRF 分两阶段：

| 阶段(searchParams) | 用途 | match/mismatch/delta | pm/pi | minscore | maxperiod | minCopyNumber |
|------|------|----------------------|-------|----------|-----------|---------------|
| `PERFECT` | 年轻串联重复（young） | 2 / 7 / 7 | 80 / 10 | 50 | 10 | 4 |
| `DIVERGED` | 古老串联重复（old） | 2 / 3 / 5 | 75 / 20 | 33 | 7 | 5 |

`runTRFStage` 用 `$searchParams eq "PERFECT"` 区分两套（`RepeatMasker:2665-2698`）；
注意第二阶段实际传入的字符串是 **`"DIVERGED"`**（`RepeatMasker:4086`，rmblast 路径；
HMMER 路径在 3356），并非表里的"OLD"——那是用途描述，源码只在 `eq "PERFECT"` 处分支，
其余一律走"古老串联重复"分支。
两阶段编号/chooseClass：PERFECT → stage 251、chooseClass `"simple"`（`RepeatMasker:3602-3603`），
DIVERGED → stage 252、chooseClass `"masking"`（rmblast 路径 `RepeatMasker:4092-4093`；
HMMER 路径 `RepeatMasker:3362-3363`）。
`lambda`（0.41/0.32）与 `mu` 数组用于把 TRF 原始分折算成与 blast 可比的 bitScore
（`RepeatMasker:2679-2697`）。`pgr rept masker` 的 TRF 两阶段（young/old）即对应这两套参数。

**TRF 结果的后处理链**（`runTRFStage`，`RepeatMasker:2729-2785`）：先过
`copyNumber > minCopyNumber` 阈值，再用 crossmatch 的 `simple1.matrix` 对每个 hit 重打分
（gapOpen=-30、gapExt=-15、xDrop=500），simple1.matrix 分阈值 20，再经 lambda/mu 折算
bitScore；仅 PERFECT 阶段额外做 `maskLevelFilter(1)`（因它要 excision，不允许重叠，
`RepeatMasker:2780-2785`）。这是 pgr 复刻 TRF 时可直接对齐的阈值链。

### E. coli 插入元件（IS）专项处理（locateISElements）

每个 batch 在正式 `runSearchStages` 之前，若存在 `generalLibDir/is.lib` 且未 `-no_is`，
先跑 `locateISElements`（`RepeatMasker:1508,1087-1104`）。E. coli 专用，pgr 处理细菌基因组
时可参考：

- **搜索参数**：`minscore=17、minmatch=15、matrix=identity.matrix`、无 gap/bandwidth/masklevel
  （`RepeatMasker:1535-1541`），用 `is.lib`（各类 IS 元件序列）对每个 batch 搜一次，输出
  `*.iscat`（`RepeatMasker:1545`）。
- **完整元件判定**：`beginis==1 && leftis==0`（比对覆盖到 IS 的 5' 端到 3' 端）才算完整
  （`RepeatMasker:1668`）；跨测序 gap 的两段 IS 若 `begin - lastend <= 2` 且方向/名一致可合并
  （`RepeatMasker:1645-1660`）。
- **TSD（靶位点重复）检测**：按元件类型查表 `%dupLengths`——IS1→9/8/10 bp、IS2→5、IS3→3、
  IS5→4、IS10→9、IS30→2、IS150→3/4、IS186→10/11、Tn1000→5（`RepeatMasker:1673-1703`），
  取左右侧翼 `dupLength` 长的序列验证是否成对；`-is_clip` 把 IS+TSD 整段从序列中剪掉
  （输出 `$file.withoutIS` 供后续 maskSource 使用，`RepeatMasker:1284-1286`）。

### rmblast 搜索与 outfmt

- **outfmt**：`-outfmt="6 score perc_sub perc_query_gap perc_db_gap qseqid qstart
  qend qlen sstrand sseqid sstart send slen kdiv cpg_kdiv transi transv cpg_sites"`
  （`NCBIBlastSearchEngine.pm:611-613`）。
  **修正**：现版本（11/20/12 起）`search()` 恒调用 `setGenerateAlignments(1)`
  （`RepeatMasker:2019`），故实际 outfmt **恒为 20 列**（611 行，末尾追加 `qseq sseq`）；
  18 列版本（613 行）仅当引擎不生成 alignment 时才走，属**死分支**。`-a` 如今只控制
  ProcessRepeats 是否产出最终 `*.align` 文件，不再决定 rmblast outfmt。
  20 列解析（`parseTabOutput`，`NCBIBlastSearchEngine.pm:980-1012`）：`kdiv→pctRawKimuraDiverge`、
  `cpg_kdiv→pctKimuraDiverge`（CpG 校正）、`cpg_sites→cpGSites` 直接装入 SearchResult；
  主 `pctDiverge` 用的是 `perc_sub`（flds[1]），`qseq/sseq`（flds[18/19]）即为比对序列。
  `pgr rept` 的 `parse_tab_row` 只取 qseqid/qstart/qend 三列做区间，无需这些列。
- **search() 失败重试**（`RepeatMasker:2037-2087`）：rmblastn 返回错误时循环降参重试。
  **注意：各分支是互斥的 if/elsif，不是顺序降档**——`bandwidth > 14 → 14`；
  `bandwidth == 4 → 1`（注释为 "Extreme measures for very long simple satellites"，仅 TRF
  长简单卫星的二/三次重检走这里）；`minmatch < 10 → minmatch++`；均不满足则打印引擎参数后
  `exit(-1)`（HMMER 引擎任何错误直接退出）。
- **搜索参数映射**（`NCBIBlastSearchEngine::getParameters`，`NCBIBlastSearchEngine.pm:444-668`）：
  恒定追加 `-num_alignments 9999999`（`NCBIBlastSearchEngine.pm:456`，让 rmblastn 返回所有
  可能的比对而非默认 top hits）；
  `minscore/bandwidth → -xdrop_ungap/-xdrop_gap/-xdrop_gap_final/-min_raw_gapped_score`（带宽
  "+"走 MaskerAid 旧换算、"-"按 gap 罚分推 band、"0"按 minScore 倍数：`xdrop_ungap=minScore*2、
  xdrop_gap_final=minScore*4、xdrop_gap=int(minScore/2)`）；
  `gap_init/ins_gap_ext → -gapopen abs(gap_init − ins_gap_ext)`（如 -lib 的 -30−(-6)=24，
  **并非直接取 gap_init**）与 `-gapextend abs(ins_gap_ext)`（=6）、`minmatch → -word_size`、
  `masklevel → -mask_level`（>0 才加）、非 basic 评分模式追加 `-complexity_adjust`、
  `-dust no`、`-num_threads 4`（有 cores 则用之）。
- **过滤链**（`runStage`，`RepeatMasker:2918-2952`）：
  1. `search`：调 rmblastn + makeblastdb（对自定义库 `processCustomLib` 先建库）；
  2. `filterContainedResults`（可选，`filterContained=1` 时滤掉被更长 hit 包含的）；
  3. `preMaskLevelFilter`（仅 stage 401/501/502/452）；
  4. `maskLevelFilter(value => $masklevel)`：按 masklevel 阈值过滤；
  5. `filterResults(...)`：按 `$filterType`（"simple"/"masking"/...）与 divergence 做最终过滤。
  6. 各 stage 的搜索参数来自 `getSearchRecipes`（minscore/minmatch/bandwidth 等）。
- **`-lib` 自定义库路径**：`processCustomLib` 先 `makeblastdb` 建库（`RepeatMasker:6549-6551`），
  因此 `-lib` 必须传未压缩的 `.fa`（见「两个实测坑」第 1 条）。

### .cat 批注格式（供 ProcessRepeats）

`$file.cat` 头部为若干 `##` 行（`## RepeatMasker version …`、`## run with <engine> version …`、
`## RM Library: …`、`## Total Sequences/Length/NonMask(排除>20bp N/X 段)/NonSub(排除所有
非 ACGT 碱基)`，batch 有重叠时还有 `## Batch Overlap Boundaries`（逐序列列出边界坐标），
`RepeatMasker:1216-1250`；总长 >10 Mb 时整个 .cat 用 gzip 输出为 `$file.cat.gz`）+ 结束标记
`## RAW Annotations:`，之后按 batch 追加原始批注（每行一个 hit，`SearchResult::AlignWithQuerySeq`
格式）。原始批注的 ID 由 `postProcessSearch` 统一编为 `m_b{batch}s{stage}i{index}`（mask）
或 `c_b…`（cut）（`RepeatMasker:5559-5567`），并把 mask 的碱基写回 seqDB 为 `X`、cut 的写回
`x`（`RepeatMasker:5493-5533`），同时把 divergence 折算为 `100*pctDiverge/(100-pctInsert)`
（`RepeatMasker:5541-5547`）。`ProcessRepeats` 把 raw 批注做合并/打分/分级
（`parseCATFile`→`processSequence`，片段合并 `cycleReJoin`，各类别专一
处理器 `preProcessLINE/preProcessDNATransp/preProcessLTR/scoreSINEPair/scoreLTRPair` 等），
最终 `divergence` 在 `filterResults` 处再次折算为 `100 * pctDiverge / (100 - pctInsert)`
（`RepeatMasker:4588-4595`），由 `-div` 阈值过滤。

### ProcessRepeats 处理流程与输出格式

`ProcessRepeats`（364KB Perl）读 .cat，先 `parseCATFile` 把 raw 批注读成
`PRSearchResult` 集合并按序列切分，再对每条序列跑 `processSequence`
（`ProcessRepeats:557`）。其内部是多轮 cycle：

- **cycle 1**：`cycleReJoin`（`ProcessRepeats:6383`）把被 batch 分片/剪切拆碎的同一元素
  碎片重连（依据 Score/PctSub/PctDel/PctIns 四元组判断是否同一元件的两段）。
- **cycle 2+**：去边缘效应批注（batch overlap 中点的重复批注）、去 masklevel 违例、
  给移位卫星改名、构建 DNA 转座子等价结构；随后按类别调用专一处理器
  `preProcessLINE`/`preProcessDNATransp`/`preProcessLTR` 与配对评分
  `scoreSINEPair`/`scoreLTRPair`/`scoreLINEPair`（LINE 两段/5' 3' 端、LTR 两 LTR 拼接、
  SINE 双端判定），`joinDNATransposonFragments` 处理 DNA 转座子断裂；refinement 片段用
  `replaceRMFragmentChainWithRefinement` 替换。
- **分度值**：`.out` 中 `div.` 为主 `pctDiverge`（`perc_sub`），经 `100*pctDiv/(100-pctIns)`
  折算；另有 Kimura（CpG 校正）值存于 SearchResult 但不直接打印到 `.out`。
- **Kimura 公式**（`SearchResult::calcKimuraDivergence`，`SearchResult.pm:1607-1683`）：
  `kimura = 100 * | −0.5 * ln( (1−2p−q) * √(1−2q) ) |`，`p`=转换率、`q`=颠换率（以 well
  characterized 碱基计）。CpG 修正：CpG 位单个转换计 1/10 转换、两个转换计 1 个颠换
  （`divCpGMod=1`）。该函数仅用于 legacy 解析路径（<2.13）；2.13+ tab 路径直接用 rmblast
  的 `kdiv/cpg_kdiv` 列。

**`.out` 格式**（`generateOutput`，`ProcessRepeats:5269-5756`）：两行表头——
`SW  perc  perc  perc  query  position in query  matching  repeat  position in repeat` 与
`score  div.  del.  ins.  sequence  begin  end  (left)  repeat  class/family  begin  end  (left)  ID`
（HMMER 引擎第一列用 `bit` 而非 `SW`；`-no_id` 不输出末列 ID）。`generateOutput` 同时统计
aggregateStats（`ProcessRepeats:5292-5317`）：按类别正则把批注归入 SINE/LINE(细分
LINEI/LINERTE/LINECR1/LINER2)/LTR(细分 LTRERV/LTRBEL/LTRGYP)/DNA(细分 DNATC1/DNAHAT/DNAPIG/
DNAH/DNAP)/PLE/RC/RNA/SATEL/SIMPLE/LOWCOMP/OTHER，覆盖长度对重叠批注去重（被包含者计入条数
不计碱基）。

**`.tbl` 格式**（`generateTableOutput`，`ProcessRepeats:5756-6348`）：首部汇总 sequences 数、
total length、GC level、bases masked(百分比)；正文按 SINE/LINE/LTR/DNA/PLE/Retroposon/RC/RNA/
Satellite/Simple/Low_complexity 分类列出 `number of elements / length occupied / percentage of
sequence`。`-excln` 时以 `NonMask`（排除 ≥20 bp 的 N/X 段）为分母算百分比。`-maskSource` 缺省时
百分比无法计算（打警告略过）。其他输出：`-a`→`.align`、`-xm`→`.out.xm`、`-ace`→`.out.ace`、
`-poly`→`.polyout`、`-gff`→`.out.gff`（gff3）、`-html`→`.out.html`。

**`.out.gff` GFF3 字段**（`generateOutput`，`ProcessRepeats:5533-5536,5667-5694`）：首行
`##gff-version 3`；每条序列首个批注前打 `##sequence-region <query> 1 <seqLen>`；
feature 行 9 列依次为——
`<query>  RepeatMasker  dispersed_repeat  <qstart>  <qend>  <PctSubst>  <+/->  .  ID=<id>;Target "Motif:<hitName>" <subjStart> <subjEnd>`。
注意 **score 列用的是 `PctSubst`（div. 百分比）而非 SW score**；`ID` 为 RM 内部分配的
`printid`（`-no_id` 时打印"empty"占位），strand 由 orientation（"C"→`-`，其余→`+`）决定。

**`.masked` 屏蔽输出**（`ProcessRepeats::maskSequence`，`ProcessRepeats:9459-9552`）：默认把
批注区间替换为 `N`；`-x` 替换为 `X`；`-xsmall` 替换为小写（`lc`）。屏蔽前按 `getQueryName`
归并批注，对重叠区间做**包含裁剪**——被前一条完全包含者跳过、部分重叠者从 `prevEnd+1` 起算
（`ProcessRepeats:9483-9490`），与 `pgr fa mask` 的区间合并语义一致。序列输出按 50 bp 折行。

**FastaDB 长度统计口径**（`FastaDB.pm`，`SeqDBI.pm` 同）：`getSeqLength`=全长；
`getSubtLength`=去掉所有非 ACGT 模糊碱基（`XNRYMK`，`FastaDB.pm:1035`）后的可替换长度；
`getXNLength`=再去掉 **≥20 bp 连续 X/N 段**（`([X,N]{20,})`，`FastaDB.pm:1032`）；`getGCLength`=
G/C 计数。三者在读取/重建 FASTA 索引时一并算好缓存（`FastaDB.pm:1007-1009`）。`.cat` 头部
`Total NonMask`/`Total NonSub` 即分别对应 `getXNLength`/`getSubtLength`。`maxIDLength`（RM 用 50）
超长会 `croak`（`FastaDB.pm:1052-1058`）。pgr `fa mask`/`rept` 若要复现 RM 的 `NonMask`
百分比分母，需按"≥20bp 连续 X/N 段"这一口径实现。

### 对 pgr rept 的启示

- `pgr rept masker` 复刻的是 **TRF 两阶段 + rmblast 区间化** 的最小闭环，不必复刻 batch
  分片/ProcessRepeats/GC 矩阵那套工程外壳；但"TRF young/old 两套参数 + minCopyNumber 阈值"
  与"rmblast 区间去重 + divergence 过滤链"的顺序值得对齐，保证与 RM 金标准可比。
- `.cat` 的原始批注格式（按 hit 逐行）与 pgr 的 PAF/区间中间表示可以互相转换，便于对拍。
- **rmblast 版本特性开关**（`NCBIBlastSearchEngine.pm:setPathToEngine`）：2.13+ 启用
  `hasTabFormat`（18/20 列 tab outfmt）+ `hasQueryThreading`，2.14.1+ 追加 `hasDBSoftMasking`
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

### 其它可执行程序与 util 脚本

源码目录无 `lib/` 子目录（依赖模块直接平铺在根目录，如 `SearchResult.pm`、
`SearchEngineI.pm`、`FastaDB.pm`、`SimpleBatcher.pm`、`TRF.pm`、`Matrix.pm`、
`NCBIBlastSearchEngine.pm` 等）。RepeatMasker 自带几个独立可执行程序：

- **DupMasker**：用已知片段重复（segmental duplication）库注释序列（`util/` 下也有配套）。
- **RepeatProteinMask**：用 RepeatPeps 蛋白库（从 FamDB 生成）做蛋白水平屏蔽
  （转座子蛋白比对，引擎为 `rmblastx`/blastp，`RepeatProteinMask:328-408`）。
- **DateRepeats**：按插入时间给重复元素定年（读取 `.align`/`.cat`，用 Kimura divergence 推断）。

`util/` 辅助脚本：`RM2Bed.py`/`rmToUCSCTables.pl`（转 BED/UCSC 表）、`bigRmskAlignBed.as`、
`rmOutToGFF3.pl`/`rmToTrackHub.pl`/`dupliconToSVG.pl`（转 GFF3/track hub/SVG）、
`buildSummary.pl`（按类别汇总）、`calcDivergenceFromAlign.pl`（重算 divergence 生成 `.divsum`，
供 `createRepeatLandscape.pl` 画 repeat landscape）、`maskFile.pl`（按 .out 屏蔽序列）、
`combineRMFiles.pl`/`renumberRMFiles.pl`、`rmOut2Fasta.pl`、`trfMask`、
`wublastToCrossmatch.pl`、`getRepeatMaskerBatch.pl`、`buildRMLibFromEMBL.pl`。
`Libraries/RepeatAnnotationData.pm`（25 MB）是 ProcessRepeats 用的重复注解数据库
（subelement 划分、families 等，`ProcessRepeats:1119,1264` 引用）。
`Matrices/` 下按引擎分 `crossmatch`/`ncbi`/`wublast` 三套打分矩阵（`ncbi/nt/*.matrix`）。

## CentOS 7（glibc 2.17）部署兼容性（2026-08-07）

**问题**：NCBI 官方 `rmblast-2.14.1+-x64-linux-GLIBC_2.31.tar.gz`
预编译包要求的最高 glibc 符号为 **GLIBC_2.29**（本机 readelf 实测），
CentOS 7 只有 glibc 2.17，直接跑不起来。

**我们代码的真实依赖**（比之前判断宽松）：
- `parse_tab_row` 只消费 rmblastn tab 输出的 qseqid/qstart/qend 三列；
  20 列 outfmt 里 2.13+ 才新增的 kdiv/cpg/qseq/sseq 等列**没有被使用**。
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
