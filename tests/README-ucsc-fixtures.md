# UCSC kent fixtures

Test data imported byte-for-byte from the UCSC Genome Browser `kent` source
tree (master zip, 2026-08). Each section lists the `tests/` location and the
original kent paths. Files already present before this import (2bit, lav,
psl->chain, axtChain pipeline) share the same origin and are not repeated here.

## `tests/axt`

- `input/blockBug.axt`, `input/blockBug.sizes`, `expected/blockBug.psl`
  - `src/hg/mouseStuff/axtToPsl/tests/{input,expected}/blockBug.*`
- `input/zeroScore.axt`, `input/hg15.sizes`, `input/rn3.sizes`,
  `expected/zeroScore.maf`
  - `src/hg/mouseStuff/axtToMaf/tests/{input,expected}/zeroScore.*`
- `input/dropOverlap1.axt`, `expected/dropOverlap1.axt`
  - `src/hg/mouseStuff/axtDropOverlap/{input,expected}/test1.axt`
- `input/hgCountAlign.axt` — `src/hg/makeDb/hgCountAlign/test.axt`
- `input/pslPretty_{empty,justComments,bad,dnax}.axt`
  - `src/hg/pslPretty/test/expected/{empty,justComments,bad,dnax}.axt`

## `tests/chain`

- `input/hg19ToHg38.{awful,noChange,someChange}.chain`,
  `expected/*.bridged.chain`
  - `src/hg/mouseStuff/chainBridge/tests/{input,expected}/`
- `input/chr22.chain`, `expected/chr22.out`
  - `src/hg/mouseStuff/chainInfo/{input,expected}/chr22.*`
- `input/{mm6ToHg17,negQ}.chain` — `src/utils/pslMap/tests/input/`
- `input/testDog.chain` — `src/hg/mouseStuff/netToAxt/tests/input/`
- `input/testIn.chain` — `src/hg/liftUp/test/input/`

## `tests/net`

- `input/testDog.net`, `expected/testDog.{noSplit,split}.axt`
  - `src/hg/mouseStuff/netToAxt/tests/{input,expected}/`
- `input/testIn.net`, `expected/testOut{Q,T}.net`
  - `src/hg/liftUp/test/{input,expected}/`

## `tests/maf`

- `input/split_test1.{maf,correct.0.maf,correct.1.maf}`
  - `src/hg/ratStuff/mafSplit/test/test1.*`
- `input/coverage_test1.{maf,txt}` — `src/hg/mouseStuff/mafCoverage/`
- `input/chrM.maf` — `src/hg/mouseStuff/mafRanges/tests/input/`
- `input/hgLoad_testDot.maf`, `input/hgLoad_testPipe.maf`
  - `src/hg/makeDb/hgLoadMaf/tests/input/test{Dot,Pipe}.maf`
- `input/addIRows_{input,expected}.maf`
  - `src/hg/ratStuff/mafAddIRows/tests/{input,expected}/`

## `tests/psl/input`

- `pslScore_test.psl`, `NP_149062.psl` — `src/utils/pslScore/tests/`
- `{annotProt,gapBothMRna,gencode.src,gencodeGenome,kgMRna,mm6ToHg17,
  mrnaRefSeq,mrnaRefSeqX,mrnaToMm6,negQ.refSeq,parent,protGencode,refSeqGen,
  retro,spAnnotIso,spGencode}.psl` — `src/utils/pslMap/tests/input/`
- `pslCTCBad.psl`, `pslCTCGood.psl` — `src/hg/checkTableCoords/tests/input/`
- `blat_refRna.psl` — `src/blat/test/basic/refRna.psl`

## `tests/2bit`

- `input/{s1,s2,s1.s2}.2bit`, `input/{s1,s2,s1.s2}.fa.gz`,
  `expected/{s1,s2,s1.s2,s1.s2_2bit}.txt`
  - `src/utils/faCount/tests/{input,expected}/` (same files as faSize tests)
- `input/{creaGeno,testPcr}.2bit` — `src/gfServer/tests/input/`

## `tests/fasta`

- `input/{ucsc_basic.fa,ucsc_uniqIc.fa,ucsc_basicIds.lst}`, `expected/*.fa`
  - `src/utils/faFilter/tests/{input,expected}/`

## `tests/gff/input`

- `*.gff3` — `src/hg/utils/gff3ToGenePred/tests/input/`
- `{acembly,ceSangerGene,frameBug,ncbi,regress,spaceInName,tigr}.gff`,
  `{flybase,nscan,regress,twinscan,vegaGene,vegaPseudo}.gtf`
  - `src/hg/lib/tests/input/genePred/`

## `tests/fastq`

- `goodHg19.fastq` — `src/hg/encode3/validateFiles/tests/`
- `{good,bad}.fastq` — `src/hg/encode/encodeValidate/test/input/badDdf/`
- `FoxP2_SL167.fastq` — `src/hg/encode/encodeValidate/test/input/chipseq/`

## Known UCSC quirks (verified while importing)

- `axtToMaf`'s own test passes `rn3.sizes hg15.sizes` (swapped relative to the
  documented `tSizes qSizes` order), so its `expected/zeroScore.maf` srcSizes do
  not match the axt target/query. See `command_axt_to_maf_ucsc_zero_score_real_sizes`.
- `faFilter` gold outputs drop header descriptions; pgr keeps them, so the fa
  filter tests compare record names only.
