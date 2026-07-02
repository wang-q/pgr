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

上述脚本完整呈现了 UCSC pairwise 比对流程的 14 个步骤。pgr 对这些步骤的覆盖情况如下：

**已有 Rust 原生实现的步骤：**

| UCSC 工具            | pgr 对应命令              | 说明                                     |
|----------------------|---------------------------|------------------------------------------|
| `lastz`              | `pgr lav lastz`           | Rust 封装，内置 UCSC 风格 preset          |
| `lavToPsl`           | `pgr lav to-psl`          | LAV → PSL 转换                           |
| `axtChain`           | `pgr psl chain`           | Rust 重实现，复用 `libs::chain` DP 引擎   |
| `chainAntiRepeat`    | `pgr chain anti-repeat`   | Rust 重实现                              |
| `chainPreNet`        | `pgr chain pre-net`       | Rust 重实现                              |
| `netSyntenic`        | `pgr net syntenic`        | Rust 重实现                              |
| `chainStitchId`      | `pgr chain stitch`        | Rust 重实现                              |
| `netSplit`           | `pgr net split`           | Rust 重实现                              |
| `netToAxt`           | `pgr net to-axt`          | Rust 重实现                              |
| `axtSort`            | `pgr axt sort`            | Rust 重实现                              |

**仍依赖外部 kent-tools 的步骤（通过 `pgr pl ucsc` 编排）：**

| UCSC 工具            | 原因                                             |
|----------------------|--------------------------------------------------|
| `chainMergeSort`     | 多文件合并排序，`pgr chain sort` 仅处理单文件     |
| `chainNet`           | net 构建算法复杂，尚未 Rust 化                    |
| `netChainSubset`     | net→chain 反向投影，逻辑较复杂                    |
| `axtToMaf`           | 可用 `pgr axt to-fas` + `pgr fas` 管线替代       |

**关键结论：**

1. **链路完整性**：pgr 已覆盖 UCSC 流程 14 步中的 10 步，剩余 4 步通过 `pgr pl ucsc` 调用外部工具。
2. **纯 Rust 优势**：已 Rust 化的步骤实现了零 panic、友好错误处理，且无需安装 kent-tools。
3. **格式互通**：pgr 的 `axt`/`chain`/`net`/`psl`/`maf` 格式与 UCSC 保持兼容，可混用 Rust 实现与外部工具。
4. **Synteny 模式**：脚本末尾的 `netFilter -syn` + `chainSplit` 分染色体流程，在 pgr 中可通过 `pgr net filter --syn` + `pgr chain split` 组合实现。
