#!/usr/bin/env bash
# Generate FastGA golden .1aln files for pgr 1aln migration tests.
#
# Uses the FastGA toolchain (FAtoGDB -> FastGA) on the E. coli MG1655 / Sakai
# genomes already under tests/genome/. The resulting .1aln is the golden input
# for reading, and the reference for round-trip (write) verification.
#
# FastGA writes its scratch files (.1gdb/.gix/.bps next to the input, plus
# .ktab/.post in the cwd), so the genomes are copied into a throwaway temp dir
# and FastGA runs there. Everything except the .1aln is cleaned on exit.
#
# Usage: scripts/gen-1aln-golden.sh [outdir]
#   outdir    directory to keep the .1aln (default: tests/genome/)

set -euo pipefail

cd "$(dirname "$0")/.."          # repo root
fastga="$PWD/FASTGA-main"
outdir="${1:-$PWD/tests/genome}"
out="$outdir/mg1655-sakai.1aln"

# FastGA invokes FAtoGDB via system(), so the toolchain dir must be on PATH.
export PATH="$fastga:$PATH"

mg1655="$PWD/tests/genome/mg1655.fa.gz"
sakai="$PWD/tests/genome/sakai.fa.gz"

scratch="$(mktemp -d)"
trap 'rm -rf "$scratch"' EXIT

# Copy genomes into scratch so FastGA's intermediates (next to input + cwd)
# all land there and are cleaned up, not in tests/genome.
cp "$mg1655" "$scratch/"
cp "$sakai" "$scratch/"

echo "== FastGA: mg1655 vs sakai =="
(
    cd "$scratch"
    "$fastga/FastGA" -k -T8 -1:"$scratch/mg1655-sakai.1aln" \
        "$scratch/mg1655.fa.gz" "$scratch/sakai.fa.gz"
)

mkdir -p "$outdir"
mv "$scratch/mg1655-sakai.1aln" "$out"

echo "== verify: ALNtoPAF reads it back =="
"$fastga/ALNtoPAF" "$out" 2>/dev/null | head -3
ls -la "$out"