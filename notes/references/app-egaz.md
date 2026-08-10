# App-Egaz 流程梳理与 pgr 对照（2026-08-06）

> 对课题组旧流程 `~/Scripts/App-Egaz`（Perl，`egaz` = **E**asy **G**enome
> **Aligner**）的系统性梳理：命令清单、核心流程、以及与 pgr 现有能力的逐项
> 对照，评估复现可行性。核心组合是 **LASTZ + UCSC kent 工具链**——这正是
> pgr `align lastz` + `pl chainnet` 已经原生复现的部分。

## 1. 概述

- 仓库：`https://github.com/wang-q/App-Egaz`，Perl（App::Egaz），2018 起维护，
  容器化分发（Singularity/Docker，`wangq/egaz:master`）；
- 定位：**易用的基因组比对流程编排器**——用 LASTZ 做两两/自比对，用 UCSC
  kent 工具链（axtChain → … → axtToMaf）做 chain/net 精修，产出 axt/maf/block FA，
  配合 faops/fasr/spanr 做序列准备与区间分析；
- 依赖的外部工具：lastz、kent-tools（16 个二进制）、faops、faToTwoBit、
  RepeatMasker、spanr、fasr、Mash、FastTree、wfmash/pafplot、blast 系等；
- 两个主流程文档：`doc/Scer.md`（S288c vs 5 株多基因组两两比对）与
  `doc/Scer-self.md`（S288c 自比对）。

## 2. 命令清单（18 个）

| 类别 | 命令 | 功能 |
|---|---|---|
| 序列准备 | `prepseq` | faops filter(-N, 简化名, 可选 -a min) → split-name/about → 可选 RepeatMasker → chr.sizes → faToTwoBit → chr.fasta.fai（faops filter -U + samtools faidx）；`--gi` 用 perl 正则去 GI 号 |
| | `partition` | 按大小分块（默认 `--chunk 10010000 --overlap 10000`，输出 `infile[start,end]` 1-based 坐标） |
| | `maskfasta` | soft/hard masking（输入 fasta + runlist.yml，`--hard` 变 N） |
| | `repeatmasker` | RepeatMasker 包装（--species Fungi 等） |
| 比对 | `lastz` | lastz 包装：`parameters.yml` 预设 set01..set10（共 10 套，均默认 `C=2` 即 lastz 内置链化；doc 流程显式 `-C 0` 覆盖以关闭内置链化改由 kent chain/net 处理）、`--isself`、`--paired`（按名字相似度挑最相近 chr 对）、`--tp/--qp`（分块后触发 normalize）、LAV 命名 `[t]vs[q].N.lav` |
| | `blastn` / `blastmatch` / `blastlink` / `exactmatch` | blast 路线（另类同源发现） |
| 格式转换 | `lav2axt` / `lav2psl` / `normalize` / `formats` | LAV → axt/psl、LAV 归一化（移植自 kentUtils blastz-normalizeLav）、格式说明 |
| 流程 | `lpcnam` | **lav-psl-chain-net-axt-maf**（UCSC 14 个 kent 命令链，`--syn` 产 synNet.maf，`--lineargap/--minscore`，默认 loose/1000） |
| | `multiz` | 多序列比对（以 target 为中心、按树 ladder 逐对 profile-profile 渐进合并，`M=10` 最小输出宽度） |
| | `fas2vcf` | block FA → VCF |
| 其他 | `template` | 多基因组模板（1_pair → 2_mash → 3_multi → 4_vcf → 9_pack_up） |
| | `raxml` | 树构建包装 |

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
另含 minimap2 自同源图、wfmash -X、blast 找旁系同源等路线。

### 3.3 多基因组（template）

`egaz template S288c RM11_1a YJM789 Spar ... --multi` 生成
`1_pair.sh`（全两两 lastz+lpcnam）→ `2_mash.sh`（Mash 距离）→
`3_multi.sh`（multiz）→ `4_vcf.sh`（fas2vcf）→ `9_pack_up.sh`。

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
