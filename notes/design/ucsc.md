# UCSC chain-net pipeline

## Pseudocat and pseudopig

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
