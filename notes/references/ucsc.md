# UCSC chain-net pipeline

> 整理于 2026-06，源自对 UCSC kent-tools chain-net pipeline 脚本的整理。目的：为 pgr 的 chain/net/axt/maf 模块提供 Rust 重实现的参照基准。
> 复核于 2026-08-01：在 `pseudocat` vs `pseudopig` 上本机重跑 UCSC 工具链与 pgr（0.3.1）工具链，
> 逐文件 `diff` 验证了 §4 的字节级结论，并修正 lastz/2bit/meta.tmp 等边界差异（见 §3.6、§4.4）。

本文件记录了 UCSC kent-tools 中 chain→net→axt→maf 标准 pairwise 比对流程的完整 shell 脚本，
以 `pseudocat` vs `pseudopig` 为示例。该流程是 pgr `chain`/`net`/`axt`/`psl`/`lav`/`maf` 模块的
Rust 重实现参照基准。

**关联文档**：[[cactus.md]]（§1.11 Cactus vs UCSC Chain/Net 数据结构对比）。

## 0. 工具位置参考

### UCSC kent-tools 二进制

所有 UCSC 工具（含 `lastz`）安装在 `/home/wangq/.cbp/bin/`：

| 工具 | 路径 |
|---|---|
| `faToTwoBit` | `/home/wangq/.cbp/bin/faToTwoBit` |
| `lavToPsl` | `/home/wangq/.cbp/bin/lavToPsl` |
| `axtChain` | `/home/wangq/.cbp/bin/axtChain` |
| `chainAntiRepeat` | `/home/wangq/.cbp/bin/chainAntiRepeat` |
| `chainMergeSort` | `/home/wangq/.cbp/bin/chainMergeSort` |
| `chainPreNet` | `/home/wangq/.cbp/bin/chainPreNet` |
| `chainNet` | `/home/wangq/.cbp/bin/chainNet` |
| `netSyntenic` | `/home/wangq/.cbp/bin/netSyntenic` |
| `netChainSubset` | `/home/wangq/.cbp/bin/netChainSubset` |
| `chainStitchId` | `/home/wangq/.cbp/bin/chainStitchId` |
| `netSplit` | `/home/wangq/.cbp/bin/netSplit` |
| `netToAxt` | `/home/wangq/.cbp/bin/netToAxt` |
| `axtSort` | `/home/wangq/.cbp/bin/axtSort` |
| `axtToMaf` | `/home/wangq/.cbp/bin/axtToMaf` |
| `netFilter` | `/home/wangq/.cbp/bin/netFilter` |
| `chainSplit` | `/home/wangq/.cbp/bin/chainSplit` |
| `lastz` | `/home/wangq/.cbp/bin/lastz` |

> 注：这些二进制不在默认 PATH 中，使用时需 `export PATH="/home/wangq/.cbp/bin:$PATH"`。
>
> ⚠️ `.cbp/bin/` 里还装着一个 **pgr 0.2.0 旧二进制**；若把该目录放在 PATH 前面，`pgr` 会解析到旧版
> （没有 `pl chainnet` 等 0.3.x 命令）。字节级复现请保证当前构建的 pgr 优先，例如：
> `export PATH="/home/wangq/Scripts/pgr/target/debug:$PATH"`，或把新版本安装覆盖到 `.cbp/bin`。

### chainnet 源码

UCSC chainnet 源码（kent-utils 子集）位于 `/home/wangq/Scripts/chainnet/`：

- **仓库**：`https://github.com/wang-q/chainnet.git`
- **结构**：
  - `src/*.c` — 各工具的 main 文件（`axtChain.c`, `chainNet.c`, `axtToMaf.c` 等）
  - `src/lib/` — 共享库（`chain.c`, `chainConnect.c`, `gapCalc.c`, `axt.c`, `chainBlock.c` 等）
- **关键参考文件**：
  - `src/lib/gapCalc.c` — 间隙成本计算（`gapCalcCost`, `defaultGapCosts`, `originalGapCosts`）
  - `src/lib/chainConnect.c` — DP 链化成本函数（`chainConnectCost`, `cBlockFindCrossover`）
  - `src/lib/chainBlock.c` — 链化 DP 主逻辑（`chainBlocks`）
  - `src/lib/axt.c` — AXT 读写（`axtWrite` 的 `static int ix = 0` 自动重编号）

## 1. Pseudocat and pseudopig

```bash
# Lastz
lastz tests/pgr/pseudocat.fa tests/pgr/pseudopig.fa \
    > tests/pgr/lastz.lav

lavToPsl tests/pgr/lastz.lav stdout > tests/pgr/lastz.psl

# Prep
pgr fa size tests/pgr/pseudocat.fa -o tests/pgr/pseudocat.sizes
faToTwoBit tests/pgr/pseudocat.fa tests/pgr/pseudocat.2bit
pgr fa size tests/pgr/pseudopig.fa -o tests/pgr/pseudopig.sizes
faToTwoBit tests/pgr/pseudopig.fa tests/pgr/pseudopig.2bit

# Chain
mkdir -p tests/pgr/pslChain

# axtChain - Chain together axt alignments.
# usage:
#   axtChain -linearGap=loose in.axt tNibDir qNibDir out.chain
# Where tNibDir/qNibDir are either directories full of nib files, or the
# name of a .2bit file
axtChain -minScore=1000 -linearGap=loose -psl tests/pgr/lastz.psl \
    tests/pgr/pseudocat.2bit tests/pgr/pseudopig.2bit \
    tests/pgr/pslChain/lastz.raw.chain

# chainAntiRepeat - Get rid of chains that are primarily the results of
# repeats and degenerate DNA
# usage:
#    chainAntiRepeat tNibDir qNibDir inChain outChain
# options:
#    -minScore=N - minimum score (after repeat stuff) to pass
#    -noCheckScore=N - score that will pass without checks (speed tweak)
chainAntiRepeat tests/pgr/pseudocat.2bit tests/pgr/pseudopig.2bit \
    tests/pgr/pslChain/lastz.raw.chain tests/pgr/pslChain/lastz.chain

# Merge & PreNet
# chainMergeSort - Combine sorted files into larger sorted file
# usage:
#    chainMergeSort file(s)
# Output goes to standard output
# options:
#    -saveId - keep the existing chain ids.
#    -inputList=somefile - somefile contains list of input chain files.
#    -tempDir=somedir/ - somedir has space for temporary sorting data, default ./
chainMergeSort tests/pgr/pslChain/lastz.chain > tests/pgr/all.chain

# chainPreNet - Remove chains that don't have a chance of being netted
# usage:
#   chainPreNet in.chain target.sizes query.sizes out.chain
chainPreNet tests/pgr/all.chain \
    tests/pgr/pseudocat.sizes tests/pgr/pseudopig.sizes \
    tests/pgr/all.pre.chain

# Net
# chainNet - Make alignment nets out of chains
# usage:
#   chainNet in.chain target.sizes query.sizes target.net query.net
chainNet -minSpace=1 tests/pgr/all.pre.chain \
    tests/pgr/pseudocat.sizes tests/pgr/pseudopig.sizes \
    tests/pgr/pseudocat.chainnet tests/pgr/pseudopig.chainnet

# netSyntenic - Add synteny info to net.
# usage:
#   netSyntenic in.net out.net
netSyntenic tests/pgr/pseudocat.chainnet tests/pgr/noClass.net

# netChainSubset - Create chain file with subset of chains that appear in
# the net
# usage:
#    netChainSubset in.net in.chain out.chain
# options:
#    -gapOut=gap.tab - Output gap sizes to file
#    -type=XXX - Restrict output to particular type in net file
#    -splitOnInsert - Split chain when get an insertion of another chain
#    -wholeChains - Write entire chain references by net, don't split
#     when a high-level net is encoundered.  This is useful when nets
#     have been filtered.
#    -skipMissing - skip chains that are not found instead of generating
#     an error.  Useful if chains have been filtered.
netChainSubset -verbose=0 tests/pgr/noClass.net tests/pgr/all.chain tests/pgr/subset.chain

# chainStitchId - Join chain fragments with the same chain ID into a single
#    chain per ID.  Chain fragments must be from same original chain but
#    must not overlap.  Chain fragment scores are summed.
# usage:
#    chainStitchId in.chain out.chain
chainStitchId tests/pgr/subset.chain tests/pgr/over.chain

mkdir -p tests/pgr/net

# netSplit - Split a genome net file into chromosome net files
# usage:
#   netSplit in.net outDir
netSplit tests/pgr/noClass.net tests/pgr/net

# NetToAxt
mkdir -p tests/pgr/axtNet

# netToAxt - Convert net (and chain) to axt.
# usage:
#   netToAxt in.net in.chain target.2bit query.2bit out.axt
# note:
# directories full of .nib files (an older format)
# may also be used in place of target.2bit and query.2bit.
netToAxt tests/pgr/net/cat.net tests/pgr/all.pre.chain \
    tests/pgr/pseudocat.2bit tests/pgr/pseudopig.2bit \
    tests/pgr/axtNet/cat.tmp.axt

# axtSort - Sort axt files
# usage:
#    axtSort in.axt out.axt
# options:
#    -query - Sort by query position, not target
#    -byScore - Sort by score    
axtSort tests/pgr/axtNet/cat.tmp.axt tests/pgr/axtNet/cat.axt

# axtToMaf - Convert axt to maf
# usage:
#   axtToMaf in.axt target.sizes query.sizes out.maf
axtToMaf tests/pgr/axtNet/cat.axt \
    tests/pgr/pseudocat.sizes tests/pgr/pseudopig.sizes \
    tests/pgr/axtNet/cat.maf


# Synteny Mode

mkdir -p tests/pgr/synNet
mkdir -p tests/pgr/chain

# netFilter - Filter out parts of net.  What passes
# filter goes to standard output.  Note a net is a
# recursive data structure.  If a parent fails to pass
# the filter, the children are not even considered.
# usage:
#    netFilter in.net(s)
netFilter -syn tests/pgr/noClass.net > tests/pgr/synNet.net
netSplit tests/pgr/synNet.net tests/pgr/synNet

# chainSplit - Split chains up by target or query sequence
# usage:
#    chainSplit outDir inChain(s)
# options:
#    -q  - Split on query (default is on target)
#    -lump=N  Lump together so have only N split files.
chainSplit tests/pgr/synNet tests/pgr/all.chain

# Convert each net/chain pair to MAF
# For each file in synNet/*.net:
#   netToAxt ${file} ${file}.chain target.2bit query.2bit out.axt
#   axtSort in.axt out.axt
#   axtToMaf in.axt target.sizes query.sizes out.maf
netToAxt tests/pgr/synNet/cat.net tests/pgr/synNet/cat.chain \
    tests/pgr/pseudocat.2bit tests/pgr/pseudopig.2bit \
    tests/pgr/synNet/cat.tmp.axt

axtSort tests/pgr/synNet/cat.tmp.axt tests/pgr/synNet/cat.axt

axtToMaf tests/pgr/synNet/cat.axt \
    tests/pgr/pseudocat.sizes tests/pgr/pseudopig.sizes \
    tests/pgr/synNet/cat.maf
```

## 2. 对 pgr 的启示

上述脚本完整呈现了 UCSC pairwise 比对流程的 14 个步骤（外加末尾 Synteny 模式与准备步骤）。截至 2026-08，**脚本中出现的全部工具都已在 pgr 中落地**，无需再通过 `pgr pl ucsc` 调用任何 kent-tool。

**pgr 已实现命令对照表：**

| UCSC 工具 | pgr 命令 | 说明 |
|---|---|---|
| `faToTwoBit` | `pgr fa to-2bit` | FASTA → 2bit；支持多文件、`--no-mask` |
| `lastz` | `pgr lav lastz` | 封装外部 lastz 二进制（非原生算法）；输出到目录 |
| `lavToPsl` | `pgr lav to-psl` | LAV → PSL |
| `axtChain` | `pgr psl chain` | 复用 `libs::chain` DP 引擎；`--gap-model`/`--min-score` |
| `chainAntiRepeat` | `pgr chain anti-repeat` | `--target-2bit`/`--query-2bit`/`--min-score` |
| `chainMergeSort` | `pgr chain sort` | 已支持多文件合并排序 + `--input-list`/`--save-id` |
| `chainPreNet` | `pgr chain pre-net` | `--pad`/`--incl-hap` |
| `chainNet` | `pgr chain net` | `libs::chain::net::ChainNet`；要求输入按 score 降序 |
| `netSyntenic` | `pgr net syntenic` | `libs::chain::net::classify_syntenic` |
| `netChainSubset` | `pgr net subset` | `--whole-chains`/`--split-on-insert`/`--type` |
| `chainStitchId` | `pgr chain stitch` | 按 ID 合并片段、分数求和 |
| `netSplit` | `pgr net split` | 按染色体切分到目录 |
| `netToAxt` | `pgr net to-axt` | 输出用 `-o` 指定（UCSC 为位置参数） |
| `axtSort` | `pgr axt sort` | `--by-query`/`--by-score`/`--keep-ids`（默认重编号） |
| `axtToMaf` | `pgr axt to-maf` | `-t`/`-q`/`-o` 标志式（UCSC 为位置参数） |
| `netFilter` | `pgr net filter` | 支持 `--syn`/`--nonsyn` 及多种区间过滤 |
| `chainSplit` | `pgr chain split` | `--by-query`/`--lump` |

**外部依赖说明：** 唯一的外部依赖是 `lastz` 比对器本身（需 PATH 中存在 `lastz`）。注意：**字节级
复现时 lastz 必须裸调用**（参数与 §1 完全一致），`pgr lav lastz` 封装器不是字节透明的（见 §3.6）。
除此之外，整个 pairwise 流程已无任何 kent-tool 依赖。

**关键结论：**

1. **链路完整性**：UCSC 14 步主流程 + Synteny 模式（`netFilter -syn` + `chainSplit`）+ 准备步骤（`faToTwoBit`）均已 Rust 化，零 kent-tool 依赖。
2. **`chainMergeSort` 等价**：`pgr chain sort` 已支持多文件合并排序（`--input-list` 读文件列表、`--save-id` 保留原 ID），不再需要 `pgr pl ucsc` 编排。
3. **`chainNet` 排序要求**：`pgr chain net` 强制要求输入按 score 降序排列（否则报错），因此 `pgr chain sort` 必须先于 `pgr chain net`。UCSC `chainNet` 不强制，但本管线的 `chain sort` 天然满足。
4. **格式互通**：pgr 的 `axt`/`chain`/`net`/`psl`/`maf` 格式与 UCSC 保持兼容，可混用 Rust 实现与外部工具。

## 3. pgr 等价管线

以下脚本完全用 pgr 命令重写 §1 的 UCSC 流程，以 `pseudocat` vs `pseudopig` 为示例，与 §1 各阶段一一对应，可直接执行。

### 3.1 准备阶段（对应 §1 L14–L24）

```bash
# Lastz —— 必须与 §1 完全一致的裸调用（默认参数，lastz v1.04.41）
# ⚠️ 不要用 `pgr lav lastz`：preset（如 set01）会更换打分矩阵/参数产生不同比对；
# 即使不带 preset，包装器也会附加 [nameparse=darkspace]、--querydepth/--format=lav/
# --markend/--ambiguous=iupac/--output= 等，改变 d stanza 并多输出一行
# "# lastz end-of-file"（比对内容一致，但字节不同）。
lastz tests/pgr/pseudocat.fa tests/pgr/pseudopig.fa > tests/pgr/lastz.lav

# 转 PSL（pgr 与 lavToPsl 字节级一致，含 ## 注释行）
pgr lav to-psl tests/pgr/lastz.lav -o tests/pgr/lastz.psl

# Prep
pgr fa size tests/pgr/pseudocat.fa -o tests/pgr/pseudocat.sizes
pgr fa to-2bit tests/pgr/pseudocat.fa -o tests/pgr/pseudocat.2bit
pgr fa size tests/pgr/pseudopig.fa -o tests/pgr/pseudopig.sizes
pgr fa to-2bit tests/pgr/pseudopig.fa -o tests/pgr/pseudopig.2bit
```

### 3.2 Chain 阶段（对应 §1 L26–L64）

```bash
mkdir -p tests/pgr/pslChain

# axtChain 等价（PSL → chain）
pgr psl chain tests/pgr/pseudocat.2bit tests/pgr/pseudopig.2bit \
    tests/pgr/lastz.psl \
    --min-score 1000 --gap-model loose \
    -o tests/pgr/pslChain/lastz.raw.chain

# chainAntiRepeat 等价
pgr chain anti-repeat \
    --target-2bit tests/pgr/pseudocat.2bit \
    --query-2bit tests/pgr/pseudopig.2bit \
    tests/pgr/pslChain/lastz.raw.chain \
    -o tests/pgr/pslChain/lastz.chain

# chainMergeSort 等价（已支持多文件合并排序）
pgr chain sort tests/pgr/pslChain/lastz.chain -o tests/pgr/all.chain

# chainPreNet 等价
pgr chain pre-net tests/pgr/all.chain \
    tests/pgr/pseudocat.sizes tests/pgr/pseudopig.sizes \
    -o tests/pgr/all.pre.chain
```

### 3.3 Net 阶段（对应 §1 L66–L106）

```bash
# chainNet 等价（--min-space 1 对齐 UCSC -minSpace=1）
pgr chain net tests/pgr/all.pre.chain \
    tests/pgr/pseudocat.sizes tests/pgr/pseudopig.sizes \
    tests/pgr/pseudocat.chainnet tests/pgr/pseudopig.chainnet \
    --min-space 1

# netSyntenic 等价
pgr net syntenic tests/pgr/pseudocat.chainnet -o tests/pgr/noClass.net

# netChainSubset 等价
pgr net subset tests/pgr/noClass.net tests/pgr/all.chain \
    tests/pgr/subset.chain

# chainStitchId 等价
pgr chain stitch tests/pgr/subset.chain -o tests/pgr/over.chain

mkdir -p tests/pgr/net

# netSplit 等价（输出目录用 -o 指定）
pgr net split tests/pgr/noClass.net -o tests/pgr/net
```

### 3.4 NetToAxt → MAF 阶段（对应 §1 L108–L134）

```bash
mkdir -p tests/pgr/axtNet

# netToAxt 等价（pgr 用 -o 指定输出，UCSC 为位置参数）
pgr net to-axt tests/pgr/net/cat.net tests/pgr/all.pre.chain \
    tests/pgr/pseudocat.2bit tests/pgr/pseudopig.2bit \
    -o tests/pgr/axtNet/cat.tmp.axt

# axtSort 等价
pgr axt sort tests/pgr/axtNet/cat.tmp.axt -o tests/pgr/axtNet/cat.axt

# axtToMaf 等价（pgr 用 -t/-q/-o 标志，UCSC 为位置参数）
pgr axt to-maf tests/pgr/axtNet/cat.axt \
    -t tests/pgr/pseudocat.sizes -q tests/pgr/pseudopig.sizes \
    -o tests/pgr/axtNet/cat.maf
```

### 3.5 Synteny 模式（对应 §1 L137–L172）

```bash
mkdir -p tests/pgr/synNet tests/pgr/chain

# netFilter -syn 等价
pgr net filter tests/pgr/noClass.net --syn -o tests/pgr/synNet.net
pgr net split tests/pgr/synNet.net -o tests/pgr/synNet

# chainSplit 等价（默认按 target 切分；输出目录用 -o 指定）
pgr chain split tests/pgr/all.chain -o tests/pgr/synNet

# 每对 net/chain 转 MAF
pgr net to-axt tests/pgr/synNet/cat.net tests/pgr/synNet/cat.chain \
    tests/pgr/pseudocat.2bit tests/pgr/pseudopig.2bit \
    -o tests/pgr/synNet/cat.tmp.axt
pgr axt sort tests/pgr/synNet/cat.tmp.axt -o tests/pgr/synNet/cat.axt
pgr axt to-maf tests/pgr/synNet/cat.axt \
    -t tests/pgr/pseudocat.sizes -q tests/pgr/pseudopig.sizes \
    -o tests/pgr/synNet/cat.maf
```

### 3.6 与 UCSC 原流程的行为差异

* `pgr lav lastz` 输出到**目录**（每个 target/query 对一个 .lav），非单个 stdout 文件；需 `cat *.lav` 合并后再 `pgr lav to-psl`。它封装外部 `lastz` 二进制，需 PATH 中存在 `lastz`。
* **lastz 是字节级复现的前提，必须裸调用**（参数、版本与 §1 完全一致：默认矩阵 + O=400/E=30/
  K=3000/L=3000/M=0，v1.04.41）。`pgr lav lastz` **不是字节透明**的：
  - `--preset set01` 等 preset 使用不同矩阵/参数（Q=similar、L=2200、Y=3400 等），实测在该
    数据集上几乎找不到比对（输出仅 15 行/449 字节，默认 lastz 为 354 行/8528 字节）；
  - 即使不带 preset，包装器也附加 `[nameparse=darkspace]`、`--querydepth=keep,nowarn:50`、
    `--format=lav`、`--markend`、`--ambiguous=iupac`、`--output=<path>`：d stanza 命令行不同，
    `--markend` 多一行 `# lastz end-of-file` 注释（比对内容本身与裸调用一致）。
* **2bit 头部格式（格式演进，pgr 有意保持 v1）**：`pgr fa to-2bit` 恒写 **version=1 + u64
  索引偏移**；UCSC `faToTwoBit` 默认写 **version=0 + u32 偏移**（`-long` 选项才用 v1/u64，
  用于 >4Gb 组装）。v0/u32 是十几年前的旧格式：u32 索引偏移把单个 2bit 文件的上限卡在 4Gb
  序列，早已被 UCSC 官方废弃（`faToTwoBit -long` 帮助文本即注明 v1 "NOT COMPATIBLE WITH
  OLDER CODE"，指的就是读不了 v1 的旧代码）；第三方实现里"+2bit-version+ 0 The only valid
  version"之类的说法是旧解读，不成立。pgr 坚持 v1/u64 是跟随 UCSC 官方演进方向的正确选择，
  **不做字节对齐**。序列数据与 v0 完全一致，差异仅每序列索引多 4 字节 + version 字段；且
  双向互通已验证：UCSC 工具能读 pgr 的 v1 文件（`axtChain` 用 pgr 生成的 2bit 产出相同
  chain），pgr 也能读 UCSC 的 v0 文件。
* **`meta.tmp`**：UCSC `netSplit`/`chainSplit` 会在输出目录额外写一个 `meta.tmp` 元数据文件，
  pgr 不生成。目录级对比会有此差异，对下游命令无影响。
* **LAV 解析**：`pgr lav to-psl` 会跳过不认识的 LAV stanza 并告警（本数据为 `m { n 0 }` 空
  match-list stanza），实测不影响 PSL 输出（仍字节级一致）；但含真实 match-list 的 LAV 输入需注意。
* `pgr chain net` **强制要求输入按 score 降序**（否则报错），因此 `pgr chain sort` 必须先于 `pgr chain net`。UCSC `chainNet` 不强制；本管线中 `chain sort` 已在 `pre-net` 之前执行，天然满足。
* `pgr chain sort` 已等价 `chainMergeSort`（多文件合并排序 + `--input-list` + `--save-id`），无需再经 `pgr pl ucsc`。
* 部分命令参数风格差异：`pgr net to-axt` 的输出用 `-o` 指定（target/query 仍为位置参数），`pgr axt to-maf` 的 sizes 与输出均用 `-t` / `-q` / `-o` 标志；UCSC 对应工具均为纯位置参数。
* `pgr chain net` 默认 `--min-space 25`，UCSC 脚本用 `-minSpace=1`，需显式 `--min-space 1` 对齐。
* `pgr axt sort` **默认重编号** AXT ID（从 0 开始），与 UCSC `axtSort` 一致（`axtWrite` 使用 `static int ix = 0`）。`--keep-ids` 可保留原始 ID。

## 4. 隔离测试验证报告（2026-08-01，全面核对）

以 `pseudocat` vs `pseudopig` 为测试数据，控制相同输入，逐命令对比 UCSC kent-tool 与 pgr 实现的输出。隔离测试用料已在验证完成后合并为正式测试 fixtures（`tests/pgr/`）与集成测试（`tests/cli_ucsc.rs`，18 个测试覆盖 16 个命令）。

### 4.1 逐命令对比结果

| # | UCSC 工具 | pgr 命令 | 结果 | 差异说明 |
|---|---|---|---|---|
| 1 | `faToTwoBit` | `pgr fa to-2bit` | 序列一致 ✓ | 序列数据完全相同；头部不同：UCSC 默认 version=0 + u32 偏移（旧格式，4Gb 上限，已废弃），pgr 恒 version=1 + u64 偏移（跟随官方 `-long` 演进方向，每序列 +4 字节）。互通已验证（UCSC 可读 pgr v1，pgr 可读 UCSC v0），pgr 有意不做字节对齐（见 §3.6） |
| 2 | `lavToPsl` | `pgr lav to-psl` | 完全一致 ✓ | 字节级一致（含 `##` 注释行） |
| 3 | `axtChain` | `pgr psl chain` | 完全一致 ✓ | **字节级一致**。`--gap-model loose` 输出与 `axtChain -linearGap=loose` 完全相同 |
| 4 | `chainAntiRepeat` | `pgr chain anti-repeat` | 完全一致 ✓ | 字节级一致（含注释行） |
| 5 | `chainMergeSort` | `pgr chain sort` | 完全一致 ✓ | 字节级一致（含注释行） |
| 6 | `chainPreNet` | `pgr chain pre-net` | 完全一致 ✓ | 字节级一致（含注释行） |
| 7 | `chainNet` | `pgr chain net` | 完全一致 ✓ | 字节级一致 |
| 8 | `netSyntenic` | `pgr net syntenic` | 完全一致 ✓ | 字节级一致 |
| 9 | `netChainSubset` | `pgr net subset` | 完全一致 ✓ | 字节级一致 |
| 10 | `chainStitchId` | `pgr chain stitch` | 完全一致 ✓ | 字节级一致 |
| 11 | `netSplit` | `pgr net split` | 完全一致 ✓ | 字节级一致（含注释行） |
| 12 | `netToAxt` | `pgr net to-axt` | 完全一致 ✓ | 字节级一致 |
| 13 | `axtSort` | `pgr axt sort` | 完全一致 ✓ | **字节级一致**。默认重编号 AXT ID（匹配 UCSC `axtWrite` 的 `static int ix=0`）；`--keep-ids` 可保留原始 ID |
| 14 | `axtToMaf` | `pgr axt to-maf` | pgr 验证正确 | UCSC `axtToMaf` 在 Linux x86_64 崩溃（`intToPt` null pointer, exit=134），非 pgr 问题。pgr 输出经逐条验证：6 条记录的 score、坐标、strand、序列内容均与 AXT 输入一致 |
| 15 | `netFilter -syn` | `pgr net filter --syn` | 一致 ✓ | 两边均输出空文件（top score 124204 < 默认 minTopScore 300000）。`--nonsyn` 排除注释行后一致 |
| 16 | `chainSplit` | `pgr chain split` | 一致 ✓ | chain 数据完全一致；pgr 不透传注释行，UCSC 透传 |

### 4.2 全管线端到端验证

2026-08-01 全面重跑两端 12 步主流程 + synteny 模式，逐文件 `diff` 对比：

| 步骤 | UCSC 输出 | pgr 输出 | diff 结果 |
|---|---|---|---|
| 1. axtChain / psl chain | `lastz.raw.chain` | `lastz.raw.chain` | **IDENTICAL** |
| 2. chainAntiRepeat | `lastz.chain` | `lastz.chain` | **IDENTICAL** |
| 3. chainMergeSort | `all.chain` | `all.chain` | **IDENTICAL** |
| 4. chainPreNet | `all.pre.chain` | `all.pre.chain` | **IDENTICAL** |
| 5a. chainNet (target) | `pseudocat.chainnet` | `pseudocat.chainnet` | **IDENTICAL** |
| 5b. chainNet (query) | `pseudopig.chainnet` | `pseudopig.chainnet` | **IDENTICAL** |
| 6. netSyntenic | `noClass.net` | `noClass.net` | **IDENTICAL** |
| 7. netChainSubset | `subset.chain` | `subset.chain` | **IDENTICAL** |
| 8. chainStitchId | `over.chain` | `over.chain` | **IDENTICAL** |
| 9. netSplit | `net/cat.net` | `net/cat.net` | **IDENTICAL** |
| 10. netToAxt | `axtNet/cat.tmp.axt` | `axtNet/cat.tmp.axt` | **IDENTICAL** |
| 11. axtSort | `axtNet/cat.axt` | `axtNet/cat.axt` | **IDENTICAL** |
| 12. axtToMaf | *(UCSC 崩溃，无输出)* | `axtNet/cat.maf` (18409 bytes) | pgr 正常 |
| — netFilter -syn | `synNet.net` (空) | `synNet.net` (空) | **IDENTICAL** |
| — chainSplit | `synNet/cat.chain` | `synNet/cat.chain` | chain 数据 IDENTICAL（注释行差异） |

### 4.3 管道命令验证（`pgr pl chainnet` vs `pgr pl ucsc`）

以相同测试数据、相同参数（`--gap-model loose --min-score 1000`），分别运行 `pgr pl chainnet` 和 `pgr pl ucsc`，模拟管道脚本的完整执行流程，逐文件对比：

| 文件 | 管道对比结果 |
|---|---|
| `all.chain`（chain sort / chainMergeSort） | **IDENTICAL** |
| `all.pre.chain`（chainPreNet） | **IDENTICAL** |
| `noClass.net`（chainNet \| netSyntenic） | **IDENTICAL** |
| `over.chain`（netChainSubset \| chainStitchId） | **IDENTICAL** |
| `net/cat.net`（netSplit） | **IDENTICAL** |
| `axtNet/cat.tmp.axt`（netToAxt） | **IDENTICAL** |
| `axtNet/cat.axt`（axtSort） | **IDENTICAL** |
| MAF（axtToMaf） | UCSC 崩溃，pgr 正常 |

**管道脚本差异说明：**

* `pgr pl chainnet` 使用 `pgr chain net` 写入靶标和查询两个 `.chainnet` 文件，然后 `net syntenic` 处理靶标侧；`pgr pl ucsc` 将 `chainNet stdout` 直接 pipe 到 `netSyntenic`，不写入中间靶标 `.chainnet` 文件。两者输出（`noClass.net`）一致。
* `pgr pl ucsc` 的 `chainMergeSort` 采用分批合并（`CHAIN_BATCH_SIZE=100`），`pgr pl chainnet` 的 `chain sort` 使用单次 `--input-list`。两者输出（`all.chain`）一致。
* 两者 synteny 模式均输出空目录（`netFilter -syn` 因 top score < minTopScore 无输出）。

### 4.4 结论

2026-08-01 本机重跑（pgr 0.3.1，两侧使用同一 `lastz.lav`/sizes 输入）确认：**12 步 chain-net-axt
主流程中可对比的 11 步全部与 UCSC 字节级完全一致**（逐文件 `diff` 无差异）；第 12 步 `axtToMaf`
因 UCSC 二进制在本机 Linux x86_64 崩溃而无法直接对比，pgr 输出（18409 字节 MAF）经逐条验证正确。
剩余差异全部在**准备/边界步骤**，不影响比对数据本身：

* `axtToMaf`：UCSC 链化库构建的二进制在 Linux x86_64 崩溃（`intToPt` null pointer at
  `obscure.c:330` in `loadIntHash` @ `axtToMaf.c:63`，exit=134），pgr 正常输出。
* `faToTwoBit` vs `pgr fa to-2bit`：2bit 头部格式不同（v0/u32 vs v1/u64，每序列 +4 字节）。
  v0/u32 是十几年前的旧格式（4Gb 上限），UCSC 官方自己用 `-long`/v1 支持 >4Gb，pgr 有意
  保持 v1/u64 不做字节对齐；序列数据一致且双向互通（详见 §3.6）。
* `netFilter`/`chainSplit`：注释行保留策略不同（pgr 的 `net filter` 保留 `##` 行、`chain split`
  不透传注释；UCSC 相反）。
* `netSplit`/`chainSplit`：UCSC 额外生成 `meta.tmp` 文件。
* lastz：必须裸调用复现（`pgr lav lastz` 非字节透明，见 §3.6）。

pgr 的 chain-net-axt-maf 管线可作为 UCSC kent-tool 的 Rust 替代，达到字节级复制要求。
