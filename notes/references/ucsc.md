# UCSC chain-net pipeline

> 整理于 2026-06，源自对 UCSC kent-tools chain-net pipeline 脚本的整理。目的：为 pgr 的 chain/net/axt/maf 模块提供 Rust 重实现的参照基准。

本文件记录了 UCSC kent-tools 中 chain→net→axt→maf 标准 pairwise 比对流程的完整 shell 脚本，
以 `pseudocat` vs `pseudopig` 为示例。该流程是 pgr `chain`/`net`/`axt`/`psl`/`lav`/`maf` 模块的
Rust 重实现参照基准。

**关联文档**：[[cactus.md]]（§1.11 Cactus vs UCSC Chain/Net 数据结构对比）。

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
| `axtSort` | `pgr axt sort` | `--by-query`/`--by-score`/`--renumber` |
| `axtToMaf` | `pgr axt to-maf` | `-t`/`-q`/`-o` 标志式（UCSC 为位置参数） |
| `netFilter` | `pgr net filter` | 支持 `--syn`/`--nonsyn` 及多种区间过滤 |
| `chainSplit` | `pgr chain split` | `--by-query`/`--lump` |

**外部依赖说明：** 唯一的外部依赖是 `lastz` 比对器本身（由 `pgr lav lastz` 封装调用，需 PATH 中存在 `lastz`）。这属于比对器而非 kent-tool，符合预期。除此之外，整个 pairwise 流程已无任何 kent-tool 依赖。

**关键结论：**

1. **链路完整性**：UCSC 14 步主流程 + Synteny 模式（`netFilter -syn` + `chainSplit`）+ 准备步骤（`faToTwoBit`）均已 Rust 化，零 kent-tool 依赖。
2. **`chainMergeSort` 等价**：`pgr chain sort` 已支持多文件合并排序（`--input-list` 读文件列表、`--save-id` 保留原 ID），不再需要 `pgr pl ucsc` 编排。
3. **`chainNet` 排序要求**：`pgr chain net` 强制要求输入按 score 降序排列（否则报错），因此 `pgr chain sort` 必须先于 `pgr chain net`。UCSC `chainNet` 不强制，但本管线的 `chain sort` 天然满足。
4. **格式互通**：pgr 的 `axt`/`chain`/`net`/`psl`/`maf` 格式与 UCSC 保持兼容，可混用 Rust 实现与外部工具。

## 3. pgr 等价管线

以下脚本完全用 pgr 命令重写 §1 的 UCSC 流程，以 `pseudocat` vs `pseudopig` 为示例，与 §1 各阶段一一对应，可直接执行。

### 3.1 准备阶段（对应 §1 L14–L24）

```bash
# Lastz（封装外部 lastz；输出到目录 lastz_out/，每个 target/query 对一个 .lav）
pgr lav lastz tests/pgr/pseudocat.fa tests/pgr/pseudopig.fa \
    --preset set01 -o tests/pgr/lastz_out

# 合并目录内所有 .lav 为单个文件，再转 PSL
cat tests/pgr/lastz_out/*.lav > tests/pgr/lastz.lav
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
* `pgr chain net` **强制要求输入按 score 降序**（否则报错），因此 `pgr chain sort` 必须先于 `pgr chain net`。UCSC `chainNet` 不强制；本管线中 `chain sort` 已在 `pre-net` 之前执行，天然满足。
* `pgr chain sort` 已等价 `chainMergeSort`（多文件合并排序 + `--input-list` + `--save-id`），无需再经 `pgr pl ucsc`。
* 部分命令参数风格差异：`pgr net to-axt` 的输出用 `-o` 指定（target/query 仍为位置参数），`pgr axt to-maf` 的 sizes 与输出均用 `-t` / `-q` / `-o` 标志；UCSC 对应工具均为纯位置参数。
* `pgr chain net` 默认 `--min-space 25`，UCSC 脚本用 `-minSpace=1`，需显式 `--min-space 1` 对齐。

## 4. 隔离测试验证报告（2026-08-01，GapCalc 修复后状态）

以 `pseudocat` vs `pseudopig` 为测试数据，控制相同输入，逐命令对比 UCSC kent-tool 与 pgr 实现的输出。隔离测试用料已在验证完成后合并为正式测试 fixtures（`tests/pgr/`）与集成测试（`tests/cli_ucsc.rs`，18 个测试覆盖 16 个命令）。

### 4.1 逐命令对比结果

| # | UCSC 工具 | pgr 命令 | 结果 | 差异说明 |
|---|---|---|---|---|
| 1 | `faToTwoBit` | `pgr fa to-2bit` | 序列一致 ✓ | 2bit 文件大小差 4 字节：UCSC version=0，pgr version=1；序列数据完全相同 |
| 2 | `lavToPsl` | `pgr lav to-psl` | 完全一致 ✓ | 修复后 pgr 保留 `##` 注释行，字节级一致 |
| 3 | `axtChain` | `pgr psl chain` | 完全一致 ✓ | **字节级一致**。修复 GapCalc 表互换 + `both=dq+dt` + ID 全局重编号 + axtChain 头注释后，`pgr psl chain --gap-model loose` 输出与 `axtChain -linearGap=loose` 完全相同 |
| 4 | `chainAntiRepeat` | `pgr chain anti-repeat` | 完全一致 ✓ | 修复后注释行透传，排除注释后字节级一致 |
| 5 | `chainMergeSort` | `pgr chain sort` | 完全一致 ✓ | 修复后注释行透传 |
| 6 | `chainPreNet` | `pgr chain pre-net` | 完全一致 ✓ | 修复后注释行透传 |
| 7 | `chainNet` | `pgr chain net` | 完全一致 ✓ | 字节级一致 |
| 8 | `netSyntenic` | `pgr net syntenic` | 完全一致 ✓ | 字节级一致 |
| 9 | `netChainSubset` | `pgr net subset` | 完全一致 ✓ | 修复后 pgr 按 UCSC `chainFastSubsetOnT` 逻辑重算子链 score（t 坐标跨度比），字节级一致 |
| 10 | `chainStitchId` | `pgr chain stitch` | 完全一致 ✓ | 字节级一致 |
| 11 | `netSplit` | `pgr net split` | 完全一致 ✓ | 修复后注释行排序一致 |
| 12 | `netToAxt` | `pgr net to-axt` | 完全一致 ✓ | 修复后输出顺序匹配 UCSC net 树 pre-order 遍历，`cat.tmp.axt` 字节级一致 |
| 13 | `axtSort` | `pgr axt sort` | 完全一致 ✓ | 字节级一致 |
| 14 | `axtToMaf` | `pgr axt to-maf` | pgr 验证正确 | UCSC `axtToMaf` 在 macOS arm64 崩溃（Trace/BPT trap, exit=133）。pgr 输出经逐条验证：6 条记录的 score、坐标、strand、序列内容均与 AXT 输入一致 |
| 15 | `netFilter -syn` | `pgr net filter --syn` | 一致 ✓ | 两边均输出空文件（top score 124204 < 默认 minTopScore 300000）。`--nonsyn` 排除注释行后一致 |
| 16 | `chainSplit` | `pgr chain split` | 一致 ✓ | chain 数据完全一致；UCSC 按 ID 排序，pgr 保留输入顺序，排序后一致 |

### 4.2 GapCalc 修复详情（2026-08-01）

此前 `axtChain` vs `pgr psl chain` 存在链数差异（UCSC 5 条 vs pgr 7 条），根因是 `GapCalc` 实现中有两处 bug：

**Bug 1：`loose` 与 `medium` 的间隙成本表互换**

UCSC `gapCalc.c` 中：
- `defaultGapCosts`（`-linearGap=loose` 别名）：qGap `[325, 360, 400, 450, 600, 1100, 3600, 7600, 15600, 31600, 56600]`
- `originalGapCosts`（`-linearGap=medium` 别名）：qGap `[350, 425, 450, 600, 900, 2900, 22900, 57900, 117900, 217900, 317900]`

pgr 的 `GapCalc::loose()` 和 `GapCalc::medium()` 恰好将两张表互换，导致 `--gap-model loose` 实际使用的是 UCSC `medium`（更高成本，不桥接大 gap），`--gap-model medium` 实际使用 UCSC `loose`。修复后两表与 UCSC 语义对齐。

**Bug 2：`calc()` 同时间隙用 `max(dq, dt)` 而非 `dq + dt`**

UCSC `gapCalcCost` 在 `dq > 0 && dt > 0` 时使用 `both = dq + dt`（两轴间隙长度之和）查 `bothGap` 表。pgr 错误地使用 `max(dq, dt)`，导致同时间隙的成本被低估。修复后改用 `dq + dt`，并加上 UCSC 的 `BIGNUM` 外推保护（外推为负值时返回 `0x3fffffff`）。

**附带修复：chain ID 全局重编号**

UCSC `axtChain` 在全局按 score 降序排序所有 chain 后，通过 `chainWriteHead → chainIdNext` 按排序顺序分配 ID 1, 2, 3, ...。pgr 此前按分组处理顺序分配 ID，导致 ID 与 UCSC 不一致。修复后在 `chain_psl` 中全局排序后统一重编号。

**附带修复：axtChain 头注释**

UCSC `axtChain` 通过 `axtScoreSchemeDnaWrite` 在文件开头写入 `##matrix=axtChain` 和 `##gapPenalties=axtChain` 两行元数据。pgr 此前不输出这两行。修复后 `SubMatrix::axt_chain_header()` 生成对应格式，在 PSL 注释行之前写入。

### 4.3 剩余差异

**`netFilter` 注释行保留——pgr 更保守**

pgr `net filter` 保留输入 net 中的 `##matrix`/`##gapPenalties` 注释行，UCSC `netFilter` 会剥离。排除注释行后数据完全一致。pgr 的行为是设计选择（注释透传不丢失元数据）。

### 4.4 结论

经 GapCalc 两处 bug 修复后，`pgr psl chain --gap-model loose` 的输出与 `axtChain -linearGap=loose` **字节级完全一致**（`diff` 无差异）。16 个命令的隔离测试表明 pgr 重实现与 UCSC kent-tool 的输出在核心数据层面完全一致。剩余 1 处差异为预期行为（`netFilter` 注释保留是 pgr 设计选择）。pgr 的 chain-net-axt-maf 管线可作为 UCSC kent-tool 的 Rust 替代。
