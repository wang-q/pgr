# App-Egaz 流程梳理与 pgr 对照（2026-08-06）

> 对课题组旧流程 `~/Scripts/App-Egaz`（Perl，`egaz` = **E**asy **G**enome
> **Aligner**，`$VERSION=0.2.9`）的系统性梳理：命令清单、核心流程、以及与 pgr 现有能力的逐项
> 对照，评估复现可行性。核心组合是 **LASTZ + UCSC kent 工具链**——这正是
> pgr `align lastz` + `pl chainnet` 已经原生复现的部分。
>
> 本版（2026-08-11）对照源码逐行核对并补充各命令实现细节、`parameters.yml` 预设、
> `lpcnam` 编排逻辑、`template` 模式，修正了此前与源码不符的描述。

## 1. 概述

- 仓库：`https://github.com/wang-q/App-Egaz`，Perl（App::Egaz），2018 起维护，
  容器化分发（Singularity/Docker，`wangq/egaz:master`）；
- 定位：**易用的基因组比对流程编排器**——用 LASTZ 做两两/自比对，用 UCSC
  kent 工具链（axtChain → … → axtToMaf）做 chain/net 精修，产出 axt/maf/block FA，
  配合 faops/fasr/spanr 做序列准备与区间分析；
- 命令层是 App::Cmd 框架（`App::Cmd::Setup -app`），入口 `script/egaz` 即 `App::Egaz->run`；
  各命令 `lib/App/Egaz/Command/*.pm` 用 MCE 并行，模板用 Template-Toolkit；
- **macOS 注意**：Egaz.pm 的 POD 明确警告 `egaz lpcnam` 中 `axtChain` 等在 macOS 下
  工作不可靠，此时应改用 `egaz lastz` 自带链化（`C=2`）而非 kent chain/net；
- 依赖的外部工具：lastz、kent-tools（lpcnam 用到 14 个二进制）、faops、
  faToTwoBit、RepeatMasker、spanr、fasr、fasops、samtools、snp-sites、bcftools、
  Mash、FastTree、RAxML(raxmlHPC)、multiz、wfmash/pafplot、blast 系（blastn/makeblastdb）、
  mummer/sparsemem、mafft、gzip/pigz、tsv-utils、linkr/rgr、nw_utils、circos 等；
- 两个主流程文档：`doc/Scer.md`（S288c vs 5 株多基因组两两比对）与
  `doc/Scer-self.md`（S288c 自比对）。

## 2. 命令清单（18 个）

| 类别 | 命令 | 功能 |
|---|---|---|
| 序列准备 | `prepseq` | faops filter(-N, 简化名, 可选 -a min) → split-name/about → 可选 RepeatMasker → chr.sizes → faToTwoBit → chr.fasta.fai（faops filter -U + samtools faidx）；`--gi` 用 perl 正则去 GI 号 |
| | `partition` | 按大小分块（默认 `--chunk 10010000 --overlap 10000`，输出 `infile[start,end]` 1-based 坐标） |
| | `maskfasta` | soft/hard masking（输入 fasta + runlist.yml，`--hard` 变 N） |
| | `repeatmasker` | RepeatMasker 包装（--species Fungi 等，默认 `-xsmall` 软遮蔽） |
| 比对 | `lastz` | lastz 包装：`parameters.yml` 预设 set01..set10（共 10 套，均默认 `C=2` 即 lastz 内置链化；doc 流程显式 `-C 0` 覆盖以关闭内置链化改由 kent chain/net 处理）、`--isself`、`--paired`（按名字相似度挑最相近 chr 对）、`--tp/--qp`（分块后触发 normalize）、LAV 命名 `[t]vs[q].N.lav` |
| | `blastn` / `blastmatch` / `blastlink` / `exactmatch` | blast 路线（另类同源发现） |
| 格式转换 | `lav2axt` / `lav2psl` / `normalize` / `formats` | LAV → axt/psl、LAV 归一化（移植自 kentUtils blastz-normalizeLav）、格式说明 |
| 流程 | `lpcnam` | **lav-psl-chain-net-axt-maf**（UCSC 14 个 kent 命令链，`--syn` 产 synNet.maf，`--lineargap/--minscore`，默认 loose/1000） |
| | `multiz` | 多序列比对（以 target 为中心、按树 ladder 逐对 profile-profile 渐进合并，`M=10` 最小输出宽度） |
| | `fas2vcf` | block FA → VCF（fasops split → snp-sites → bcftools concat） |
| 其他 | `template` | 多基因组模板（multi/self/prep 三模式，产 1_pair/2_mash/3_multi/4_vcf/9_pack_up 等 .sh） |
| | `raxml` | 树构建包装（fasops → raxmlHPC GTRGAMMA） |

### 2.1 命令实现细节（逐命令核对源码）

- **`prepseq`**：输入可为单文件或目录。
  - 单文件模式：`faops filter -N -s [-a min]`（IUPAC→N、简化序列名、可选 `-a min` 去短）
    → 可选 `--gi` 用 perl 正则 `s/gi\|\d+\|\w+\|//` 去 GI 号 → 有 `--about` 走
    `faops split-about stdin <about> outdir`（按约大小分块），否则 `faops split-name stdin outdir`
    （按名字拆分）；gzip 输入 OK。
  - 目录模式：`--outdir` 强制为输入目录，`--min/--about` 被忽略，仅处理 `*.fa`
    （要求每个文件单序列且文件名=序列名）；在临时目录里逐文件 `faops filter -N -s`
    后拷回。
  - 可选 `--repeatmasker "..."` → 调 `egaz repeatmasker $outdir/*.fa -o $outdir`。
  - 产出四样：`chr.sizes`（`faops size`）、`chr.2bit`（`faToTwoBit`）、
    `chr.fasta`+`chr.fasta.fai`（`cat *.fa | faops filter -U` 后 `samtools faidx`）。
- **`partition`**：读 fasta（App::Fasops），校验单序列；若 `size > chunk+overlap`
  用 `overlap_ranges(1, size, chunk, overlap)` 生成 1-based `[start,end]` 区间并
  `touch` 创建空占位文件 `infile[start,end]`（lastz 据此取子区间），否则仅 `infile[1,size]`。
  默认 `chunk=10_010_000, overlap=10_000`。
- **`maskfasta`**：读 fasta + runlist.yml（YAML→AlignDB::IntSpan），对 runlist 覆盖区间
  soft（小写）或 `--hard`（N）；`--len` 行宽默认 80。
- **`repeatmasker`**：在 tempdir 中对每个输入跑
  `RepeatMasker <in> -dir . [-species] [--opt] -xsmall --parallel N`（默认 `-xsmall` 软遮蔽）；
  把 `<name>.masked` 拷为 `<basename>.fa`、`.out` 拷为 `<basename>.rm.out`；
  `--gff` 用 `rmOutToGFF3.pl`（自动定位 libexec/util）产出 `.rm.gff`。
- **`blastn`**：`makeblastdb -dbtype nucl` → `blastn -task megablast
  -max_target_seqs 20 -culling_limit 20 -dust no -soft_masking false
  -evalue(0.01) -word_size(40) -outfmt "7 qseqid sseqid qstart qend sstart send qlen slen nident"`，
  并按 `blastn -h` 探测版本自动补 `-max_hsps 10`（2.6.0+）或 `-max_hsps_per_subject 10`。
- **`blastmatch`**：解析 blastn 9 列报表，保留 query 同源覆盖 `nident/qlen ≥ --coverage(0.9)`
  的命中合并为 runlist；`--perchr` 逐染色体输出一条（可重叠）；邻近（vicinity=10）区间
  用 IntSpan set 包含关系裁剪嵌套。
- **`blastlink`**：解析报表，挑 query/hit 长度与 nident 相对 max 均 ≥0.9 的 pair →
  `query\thit\tstrand` 链接，去重。
- **`exactmatch`**：mummer/sparsemem（`-maxmatch -F -l length -b -n [-k 4]`，默认长度 100）
  找全基因组精确匹配，仅保留 query 起点=1 且长度=整条 query 的"整条命中"；`--discard`
  丢弃拷贝数超限者。
- **`lav2axt`**：解析 LAV 各 stanza，按 h/a/b/e/l 生成 axt（反向链对 query 做 revcom），
  缓存 fasta，header 只留首 token。
- **`lav2psl`**：解析 LAV → psl 21 列（计算 match/mismatch、双侧插入、strand、
  sizes/qStarts/tStarts），反向链按 `qSeqStops` 换算 q 坐标。
- **`normalize`**：移植自 kentUtils `blastz-normalizeLav`；s-stanza 的 start/stop 延伸到
  `--tlen/--qlen`，a/b/e/l/m 坐标按 `t_from` 偏移、`qlen`（含反向链 `qlen-to+x`）重映射，
  输出 1-based。
- **`multiz`**：扫描 `.maf`/`.maf.gz`，`gzip -dcf | perl` 取每文件物种（≤2），聚合成
  species→chr→files 映射并选 target（同 chr 出现最多的"潜在 target"）；物种数 <3 时退化为
  直接拷贝 pairwise maf；有 `--tree` 用 Bio::Phylo 按 ladder（从 target 逐层向根、按
  patristic 距离排序）定 stitch 顺序，否则按命令行/出现次数顺序；对每个 chr 逐对
  `multiz M=10 maf1 maf2 1 out1 out2 > stepN.maf` 渐进合并（profile-profile），末步得
  `chr.maf` → gzip；产物 `info.yml`、`steps.csv`（`--keeptmp` 保留中间 .step/.out 文件）。
- **`fas2vcf`**：`fasops subset --list`（可选）→ `fasops split --simple` 到 tempdir →
  每个 .fas 用 `snp-sites -v` 转 vcf，并把假染色体名（chrUn/1）改为真实 chr、位置偏移回
  基因组坐标 → `bcftools concat`。
- **`raxml`**：`fasops names` + `gzip -dcf * | fasops concat --relaxed` 生成 phylip →
  自动探测 `raxmlHPC[-PTHREADS-*]` 变体，`-T N -f a -m GTRGAMMA -p seed -x seed -N bootstrap(100)`
  （可选 `-o outgroup`）→ 输出 `RAxML_bipartitions`。

## 3. 核心流程（doc/Scer.md 提取）

### 3.1 两两比对（LASTZ + UCSC）

```
prepseq：faops filter -N -s（去 N、简化名）→ split-name/about
  → 可选 RepeatMasker（--species Fungi）→ chr.sizes → chr.2bit → chr.fasta.fai
lastz --set set01 -C 0：全基因组对 → [t]vs[q].N.lav
lpcnam：LAV → PSL → axtChain(loose/1000) → chainAntiRepeat → chainMergeSort
  → chainPreNet → chainNet → netSyntenic → netChainSubset → chainStitchId
  → netSplit → netToAxt → axtSort → axtToMaf（+ netFilter/chainSplit）
  → axtNet/*.net.axt.gz 或 mafSynNet/*.synNet.maf.gz（--syn）
后处理：fasr axt2fas/maf2fas → fasr filter --ge 1000 → fasr check 校验
  → spanr cover/stat（覆盖/N50/depth）→ dotplot（wfmash + pafplot）
```

实测结果（S288c vs RM11_1a，Scer.md 表）：lav2axt 直连覆盖 96.3%、
lpcnam（chain/net）95.9%、lpcnam --syn 94.9%、partition 95.9%。

### 3.2 自比对（Scer-self.md）

`lastz --isself` + `lpcnam`（同文件对）→ axt → block FA → spanr/rgr 分析；
另含 minimap2 自同源图、wfmash -X、blast 找旁系同源（axt→`egaz blastn`→
`blastmatch`/`blastlink` + linkr/rgr/fasr 合并链接）等路线。

### 3.3 多基因组（template）

`egaz template S288c RM11_1a YJM789 Spar ... --multi` 生成
`1_pair.sh`（全两两 lastz+lpcnam，`--set set01 -C 0`，`--partition` 时加 `--tp --qp`）→
`2_mash.sh`（Mash 距离）→ `3_multi.sh`（multiz，引导树 `--tree>--order>--mash`）→
`4_vcf.sh`（fas2vcf，**需 `--vcf` 才生成**）→ `9_pack_up.sh`（打包 tar.gz）。
`--multi` 默认用 RAxML 建树，`--fasttree` 改用 FastTree；`--aligndb` 额外生成
`6_chr_length.sh`/`7_multi_aligndb.sh`。`--self` 模式生成 `1_self/2_mash/3_proc[/4_circos(需
--circos)]/9_pack_up`；`--prep` 模式生成 `0_prep.sh`（对 NCBI 下载目录逐文件 prepseq）。

### 3.4 lpcnam 编排细节（逐命令核对）

输入 `<path/lav>` 可为单个 `.lav`、`lav.tar.gz` 或目录（目录时 `--outdir` 即该目录）；
`--tname/--qname` 默认取 `<path/target>/<path/query>` basename，作为 MAF 的 `-tPrefix/-qPrefix`。
`gzip` 在检测到 pigz 时替换为 `pigz -p N`。流程：

1. **lavToPsl**（egaz 自带，非 kent）：逐 LAV → psl（并行）。
2. **axtChain**：`-minScore=--minscore -linearGap=--lineargap -psl in.psl t.chr.2bit
   q.chr.2bit stdout | chainAntiRepeat t.2bit q.2bit stdin out.chain`（并行）。
3. **chainMergeSort / chainPreNet**：因单进程会超系统 maxfile 限制，先按 **100 个 .chain 一批**
   合并成 `all.N.chain.tmp`，再全部合并成 `all.chain`；随后
   `chainPreNet all.chain t.chr.sizes q.chr.sizes all.pre.chain`。
4. **chain-net**：`chainNet -minSpace=1 all.pre.chain t.sizes q.sizes stdout query.chainnet |
   netSyntenic stdin noClass.net` → `netChainSubset noClass.net all.chain stdout |
   chainStitchId stdin over.chain` → `netSplit noClass.net net/`。
5. **netToAxt**：逐 `net/*.net`：`netToAxt net t.2bit q.2bit stdout | axtSort stdin stdout |
   gzip > axtNet/*.axt.gz`（并行）。
6. **打包清理**：`lav.tar.gz`、`net.tar.gz`、`psl.tar.gz`、`chain.tar.gz`，删除中间 .lav/.psl/
   .chain/.tmp 文件。
7. **非 --syn**：`axtToMaf -tPrefix/-qPrefix *.axt(.gz) t.sizes q.sizes stdout | gzip →
   mafNet/*.maf.gz`。
8. **--syn**：`netFilter -syn noClass.net | netSplit stdin synNet/`；`chainSplit chain/
   all.chain.gz`；逐 synNet 的 net：`netToAxt | axtSort | axtToMaf | gzip →
   mafSynNet/*.synNet.maf.gz`；清理 `synNet/`、`chain/`。

## 4. egaz 命令 → pgr 对照

**工具级映射（三个旧工具被 pgr 命令族整体吸收）**：

| 旧工具 | pgr 对应 | 说明 |
|---|---|---|
| **faops**（序列操作） | **`pgr fa`** | filter/size/n50/masked/split/window 等子命令 |
| **spanr**（runlist 区间操作） | **`pgr runlist`**（+ `pgr rg`） | merge/stat/statop 在 runlist；cover/coverage 迁至 rg（同源） |
| **fasr**（block FA 分析） | **`pgr fas`** | axt2fas/maf2fas → `fas`；filter/check/subset 同族 |

以下为命令级对照：

| egaz | pgr 对应 | 状态 |
|---|---|---|
| `lastz`（set01..10 预设、均 C=2；doc 流程 `-C 0`；`--isself`/`--paired`/`--tp`/`--qp`） | `align lastz`（仅预设 set01..07 且硬编码 C=0；`--self`；Cactus 风格） | ◐ 核心覆盖，参数有差异（见 §5.1） |
| `lpcnam`（UCSC 14 个 kent 命令链） | `pl chainnet`（psl→chain→net→axt→maf，`--syn`） | ✅ 字节级一致（verify-ucsc-pipeline.sh，pgr 共复现 16 个含 faToTwoBit/lavToPsl） |
| `lav2psl` | `lav to-psl` | ✅ |
| `lav2axt` / `normalize` | 无直接命令；流程可 lav→psl→chainnet→axt 替代 | ◐ 可组装 |
| `multiz` | `fas multiz` | ✅ |
| `fas2vcf` | `fas to-vcf` | ✅ |
| `maskfasta` | `fa mask`（soft/hard） | ✅ |
| `repeatmasker`（包装） | `rept`（e-kmer/e-align/s-kmer/trf） | ✅ 原生替代（含遮蔽验证） |
| `partition` | `fa window`（`--window/--step` 对应 chunk/overlap） | ◐ 参数映射 |
| `prepseq`（faops filter/size/split + faToTwoBit） | `fa filter`/`size`/`split`/`window`、`fa to-2bit`、`2bit size` | ✅ 可组装 |
| `fasr axt2fas` / `maf2fas` | `axt to-fas` / `maf to-fas` | ✅ |
| `fasr filter` / `check` / `subset` | `fas filter` / `check` / `subset` | ✅ |
| `spanr gff/cover/merge/stat/coverage` | `gff runlist/rg`、`runlist merge/stat`、`rg cover/coverage` | ✅ |
| `faops n50` / `masked` | `fa n50` / `fa masked` | ✅ |
| `template`（pair/mash/multi/vcf） | 无直接命令；`verify-pangenome.sh` 已演示 10 基因组 pair→chainnet→paf→query 流程；Mash 为外部 | ◐ 脚本组装 |
| `blastn`/`blastmatch`/`blastlink`/`exactmatch` | 无（blast 是 NCBI 外部工具，pgr 不做） | ✗ 不覆盖 |
| `raxml` | 无（phylogeny 已迁 necom） | ✗ 不覆盖 |
| `wfmash`/`pafplot`（dotplot） | `plot dot` | ✅ |

## 5. 复现评估

**核心流程（LASTZ + UCSC chain/net）pgr 已完整覆盖且字节级一致**——`lpcnam`
的 14 个 kent 命令链就是 `pgr pl chainnet`（`verify-ucsc-pipeline.sh` 固化）。
因此 App-Egaz 的主线（prepseq → lastz → lpcnam → fas 分析）在 pgr 上可以
端到端复现，且去掉 kent-tools/faops/spanr/fasr 依赖。

需要脚本组装的部分（pgr 无单命令编排）：

- `prepseq` 序列准备流水（fa filter → split/window → to-2bit → mask），
  可写成 `scripts/` 的 bash 流程或 pgr 命令链；
- `template` 多基因组全两两 + multiz 模板，可用
  `scripts/verify-pangenome.sh` 的模式扩展（该脚本已演示 10 基因组
  FastGA/chainnet → PAF → query → graph 的等价流程）；
- `partition` 大基因组分块（`fa window` 参数映射）；
- `lav normalize` 若确需，可评估移植（当前 lastz 输出无需归一化也能走
  `lav to-psl`）。

### 5.1 `lastz` 包装差异与可借鉴点

egaz 的 `lastz` 与 pgr `align lastz` 都是 lastz→LAV 的薄包装，但参数面并不完全对齐：

- **预设集**：egaz `parameters.yml` 有 set01..set10（10 套，全部 `C=2`，即 lastz 内置链化）；
  pgr 预设仅移植 set01..set07 且硬编码 `C=0`（关闭内置链化，等价 egaz doc 流程的 `-C 0`
  覆盖，配合 kent chain/net）。set08..set10 是远缘物种的 `Q=distant` 变体，pgr 未移植，可补。
- **`--paired`**：egaz 用 `String::Similarity` 按文件名相似度为每个 target 挑最相近的一个
  query，做一对一（近缘菌株 chr 一一对应）比对；pgr 无此选项，只能全两两笛卡尔积。
- **`--tp/--qp` 分块 + 自动 normalize**：egaz 在 target/query 被 `partition` 分块后，自动对每个
  LAV 调 `normalize`（以 `chr.sizes` 的 tlen/qlen 重映射 a/b/e/l 坐标，移植自 kentUtils
  blastz-normalizeLav），再喂 `lav2psl`。pgr `align lastz` 不触发此步；分块比对需自行
  `fa window` 分块并按块重映射坐标。
- **多序列文件**：egaz 静默取 target 文件第一条序列（其余丢弃）；pgr 严格要求单序列文件
  （LAV 与 lastz `[multiple]` 不兼容），多 contig 直接报错——更安全。
- **Cactus 风格**：pgr 额外带 `--querydepth=keep,nowarn:N`、`[nameparse=darkspace]`、
  `--markend`、`--ambiguous=iupac`，面向 Cactus RepeatMasking 工作流，是 egaz 没有的。

借鉴点：若要复现 Scer-self / 旁系同源场景，`--paired`（相似度挑对）可作为 `align lastz`
的可选模式；分块 + normalize 目前由 `fa window` 分块 + 坐标重映射手工替代。

不覆盖：blast 系（NCBI 外部工具）、raxml（phylogeny 已迁 necom）。
若需复现，建议以 `doc/Scer.md` 的 S288c vs RM11_1a 为样板，用 pgr 命令链
重跑一遍并对比文中的覆盖/N50 表（0.9592 vs 0.9632 等），作为端到端验证。
