#!/usr/bin/env bash

# verify-ucsc-pipeline.sh
# Byte-for-byte regression check: pgr native chainnet vs UCSC kent-tools pipeline.
#
# Uses the saved E. coli LAV (tests/genome/mg1655-sakai.lastz.lav) as input so
# lastz (~2 min) does not have to be re-run.  Requires:
#   - pgr binary (defaults to target/debug/pgr; use PGR=... to override)
#   - kent-tools in PATH (axtChain, chainAntiRepeat, ..., axtToMaf)
#
# Usage:
#   export PATH="/home/wangq/.cbp/bin:$PATH"
#   scripts/verify-ucsc-pipeline.sh
#   GAP_MODEL=medium MIN_SCORE=5000 scripts/verify-ucsc-pipeline.sh

set -euo pipefail

cd "$(dirname "$0")/.."
PGR="${PGR:-$PWD/target/debug/pgr}"
GAP_MODEL="${GAP_MODEL:-loose}"
MIN_SCORE="${MIN_SCORE:-1000}"
LAV="$PWD/tests/genome/mg1655-sakai.lastz.lav"
FA_T="$PWD/tests/genome/mg1655.fa.gz"
FA_Q="$PWD/tests/genome/sakai.fa.gz"

if [ ! -x "$PGR" ]; then
    echo "Error: pgr binary not found at $PGR (run 'cargo build' first)." >&2
    exit 1
fi
for cmd in axtChain chainAntiRepeat chainMergeSort chainPreNet chainNet \
    netSyntenic netChainSubset chainStitchId netSplit netToAxt axtSort axtToMaf lavToPsl; do
    command -v "$cmd" >/dev/null || { echo "Error: kent tool '$cmd' not in PATH." >&2; exit 1; }
done

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

echo "==> Prepare (lavToPsl, sizes, 2bit)"
lavToPsl "$LAV" "$WORK/in.psl"
"$PGR" fa size "$FA_T" -o "$WORK/target.chr.sizes"
"$PGR" fa size "$FA_Q" -o "$WORK/query.chr.sizes"
"$PGR" fa to-2bit "$FA_T" -o "$WORK/target.chr.2bit"
"$PGR" fa to-2bit "$FA_Q" -o "$WORK/query.chr.2bit"

echo "==> UCSC kent-tools pipeline"
mkdir -p "$WORK/kent"
axtChain -minScore="$MIN_SCORE" -linearGap="$GAP_MODEL" -psl "$WORK/in.psl" \
    "$WORK/target.chr.2bit" "$WORK/query.chr.2bit" "$WORK/kent/01.chain"
chainAntiRepeat "$WORK/target.chr.2bit" "$WORK/query.chr.2bit" \
    "$WORK/kent/01.chain" "$WORK/kent/02.ar.chain"
chainMergeSort "$WORK/kent/02.ar.chain" > "$WORK/kent/03.all.chain"
chainPreNet "$WORK/kent/03.all.chain" "$WORK/target.chr.sizes" "$WORK/query.chr.sizes" \
    "$WORK/kent/04.pre.chain"
chainNet "$WORK/kent/04.pre.chain" "$WORK/target.chr.sizes" "$WORK/query.chr.sizes" \
    "$WORK/kent/05.target.net" "$WORK/kent/05.query.net"
netSyntenic "$WORK/kent/05.target.net" "$WORK/kent/06.syn.net"
netChainSubset "$WORK/kent/06.syn.net" "$WORK/kent/04.pre.chain" "$WORK/kent/07.subset.chain"
chainStitchId "$WORK/kent/07.subset.chain" "$WORK/kent/08.over.chain"
netSplit "$WORK/kent/06.syn.net" "$WORK/kent/net"
netToAxt "$WORK/kent/net/NC_000913.net" "$WORK/kent/04.pre.chain" \
    "$WORK/target.chr.2bit" "$WORK/query.chr.2bit" "$WORK/kent/09.axt"
axtSort "$WORK/kent/09.axt" "$WORK/kent/10.sorted.axt"
axtToMaf -tPrefix=mg1655. -qPrefix=sakai. "$WORK/kent/10.sorted.axt" \
    "$WORK/target.chr.sizes" "$WORK/query.chr.sizes" "$WORK/kent/11.maf"

echo "==> pgr native pipeline (chainnet)"
mkdir -p "$WORK/pgr" "$WORK/pgr/pslChain"
"$PGR" psl chain "$WORK/target.chr.2bit" "$WORK/query.chr.2bit" "$WORK/in.psl" \
    --min-score "$MIN_SCORE" --gap-model "$GAP_MODEL" -o "$WORK/pgr/pslChain/01.tmp"
"$PGR" chain anti-repeat --target-2bit "$WORK/target.chr.2bit" \
    --query-2bit "$WORK/query.chr.2bit" "$WORK/pgr/pslChain/01.tmp" \
    -o "$WORK/pgr/pslChain/02.ar.chain"
"$PGR" chain sort "$WORK/pgr/pslChain/02.ar.chain" -o "$WORK/pgr/03.all.chain"
"$PGR" chain pre-net "$WORK/pgr/03.all.chain" \
    "$WORK/target.chr.sizes" "$WORK/query.chr.sizes" -o "$WORK/pgr/04.pre.chain"
"$PGR" chain net "$WORK/pgr/04.pre.chain" \
    "$WORK/target.chr.sizes" "$WORK/query.chr.sizes" \
    "$WORK/pgr/05.target.net" "$WORK/pgr/05.query.net"
"$PGR" net syntenic "$WORK/pgr/05.target.net" -o "$WORK/pgr/06.syn.net"
"$PGR" net subset "$WORK/pgr/06.syn.net" "$WORK/pgr/04.pre.chain" "$WORK/pgr/07.subset.chain"
"$PGR" chain stitch "$WORK/pgr/07.subset.chain" -o "$WORK/pgr/08.over.chain"
"$PGR" net split "$WORK/pgr/06.syn.net" -o "$WORK/pgr/net"
"$PGR" net to-axt "$WORK/pgr/net/NC_000913.net" "$WORK/pgr/04.pre.chain" \
    "$WORK/target.chr.2bit" "$WORK/query.chr.2bit" -o stdout |
    "$PGR" axt sort stdin -o "$WORK/pgr/10.sorted.axt"
"$PGR" axt to-maf "$WORK/pgr/10.sorted.axt" \
    --t-prefix mg1655. --q-prefix sakai. \
    -t "$WORK/target.chr.sizes" -q "$WORK/query.chr.sizes" -o "$WORK/pgr/11.maf"

echo "==> Compare intermediate files and final MAF"
FAIL=0
compare() {
    if cmp -s "$1" "$2"; then
        echo "  OK   $(basename "$1")"
    else
        echo "  FAIL $(basename "$1")"
        FAIL=1
    fi
}
compare "$WORK/kent/01.chain"          "$WORK/pgr/pslChain/01.tmp"
compare "$WORK/kent/02.ar.chain"       "$WORK/pgr/pslChain/02.ar.chain"
compare "$WORK/kent/03.all.chain"      "$WORK/pgr/03.all.chain"
compare "$WORK/kent/04.pre.chain"      "$WORK/pgr/04.pre.chain"
compare "$WORK/kent/05.target.net"     "$WORK/pgr/05.target.net"
compare "$WORK/kent/05.query.net"      "$WORK/pgr/05.query.net"
compare "$WORK/kent/06.syn.net"        "$WORK/pgr/06.syn.net"
compare "$WORK/kent/07.subset.chain"   "$WORK/pgr/07.subset.chain"
compare "$WORK/kent/08.over.chain"     "$WORK/pgr/08.over.chain"
compare "$WORK/kent/net/NC_000913.net" "$WORK/pgr/net/NC_000913.net"
compare "$WORK/kent/10.sorted.axt"     "$WORK/pgr/10.sorted.axt"
compare "$WORK/kent/11.maf"            "$WORK/pgr/11.maf"

echo "==> Syntenic (--syn) mode"
mkdir -p "$WORK/syn"
netFilter -syn "$WORK/kent/06.syn.net" > "$WORK/syn/syn.net"
netSplit "$WORK/syn/syn.net" "$WORK/syn/net"
chainSplit "$WORK/syn/chains" "$WORK/kent/04.pre.chain"
netToAxt "$WORK/syn/net/NC_000913.net" "$WORK/syn/chains/NC_000913.chain" \
    "$WORK/target.chr.2bit" "$WORK/query.chr.2bit" "$WORK/syn/out.axt"
axtSort "$WORK/syn/out.axt" "$WORK/syn/sorted.axt"
axtToMaf -tPrefix=mg1655. -qPrefix=sakai. "$WORK/syn/sorted.axt" \
    "$WORK/target.chr.sizes" "$WORK/query.chr.sizes" "$WORK/syn/out.maf"

"$PGR" pl chainnet --syn --gap-model "$GAP_MODEL" --min-score "$MIN_SCORE" \
    "$FA_T" "$FA_Q" "$WORK/in.psl" -o "$WORK/syn_pgr" >/dev/null
compare "$WORK/syn/out.maf" "$WORK/syn_pgr/NC_000913.maf"

if [ "$FAIL" -eq 0 ]; then
    echo "PASS: pgr chainnet (gap-model=$GAP_MODEL, min-score=$MIN_SCORE, normal + --syn) is byte-for-byte identical."
else
    echo "FAIL: differences found (see above)." >&2
    exit 1
fi
