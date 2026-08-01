#!/usr/bin/env bash

# verify-pangenome.sh
# Three-E. coli pangenome smoke test: FastGA pairwise PSL -> pgr pl chainnet
# (chain -> net -> axt -> maf, the mandatory syntenic route) -> maf to-paf ->
# multi-file index -> transitive query -> coarse graph -> stat.  Requires:
#   - pgr binary (defaults to target/debug/pgr; use PGR=... to override)
#   - FastGA in PATH (pangenome-route upstream aligner, see notes/ecoli-cohort.md)
#
# Usage:
#   scripts/verify-pangenome.sh

set -euo pipefail

cd "$(dirname "$0")/.."
PGR="${PGR:-$PWD/target/debug/pgr}"
FA_M="$PWD/tests/genome/mg1655.fa.gz"
FA_K="$PWD/tests/genome/sakai.fa.gz"
FA_S="$PWD/tests/genome/se11.fa.gz"

if [ ! -x "$PGR" ]; then
    echo "Error: pgr binary not found at $PGR (run 'cargo build' first)." >&2
    exit 1
fi
command -v FastGA >/dev/null || { echo "Error: FastGA not in PATH." >&2; exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

echo "==> 1. FastGA pairwise PSL (query-target)"
FastGA -v -psl "$FA_K" "$FA_M" > "$WORK/sakai-mg1655.psl" 2>/dev/null
FastGA -v -psl "$FA_S" "$FA_M" > "$WORK/se11-mg1655.psl" 2>/dev/null
FastGA -v -psl "$FA_S" "$FA_K" > "$WORK/se11-sakai.psl" 2>/dev/null

echo "==> 2. pgr pl chainnet --syn (mandatory syntenic route)"
"$PGR" pl chainnet --syn "$FA_M" "$FA_K" "$WORK/sakai-mg1655.psl" -o "$WORK/cn_sm" >/dev/null 2>&1
"$PGR" pl chainnet --syn "$FA_M" "$FA_S" "$WORK/se11-mg1655.psl" -o "$WORK/cn_sm11" >/dev/null 2>&1
"$PGR" pl chainnet --syn "$FA_K" "$FA_S" "$WORK/se11-sakai.psl" -o "$WORK/cn_ss" >/dev/null 2>&1

echo "==> 3. MAF -> PAF (per-pair dirs avoid same-name collisions) and merge"
mkdir -p "$WORK/paf_sm" "$WORK/paf_sm11" "$WORK/paf_ss"
for m in "$WORK"/cn_sm/*.maf; do
    "$PGR" maf to-paf "$m" -o "$WORK/paf_sm/$(basename "$m" .maf).paf" >/dev/null 2>&1
done
for m in "$WORK"/cn_sm11/*.maf; do
    "$PGR" maf to-paf "$m" -o "$WORK/paf_sm11/$(basename "$m" .maf).paf" >/dev/null 2>&1
done
for m in "$WORK"/cn_ss/*.maf; do
    "$PGR" maf to-paf "$m" -o "$WORK/paf_ss/$(basename "$m" .maf).paf" >/dev/null 2>&1
done
cat "$WORK"/paf_sm/*.paf "$WORK"/paf_sm11/*.paf "$WORK"/paf_ss/*.paf > "$WORK/all.paf"

echo "==> 4. Multi-file PAF index"
"$PGR" paf index "$WORK"/paf_sm/*.paf "$WORK"/paf_sm11/*.paf "$WORK"/paf_ss/*.paf \
    -o "$WORK/pangenome.paf.idx" >/dev/null 2>&1

echo "==> 5. Transitive query (three-way synteny from MG1655)"
"$PGR" paf query "$WORK/pangenome.paf.idx" mg1655.NC_000913:100000-110000 \
    --transitive -o "$WORK/tri.paf" >/dev/null 2>&1
grep -q "se11.NC_011415" "$WORK/tri.paf"  # SE11 chromosome reached via BFS
grep -q "sakai.NC_002695" "$WORK/tri.paf" # Sakai chromosome reached via BFS

echo "==> 6. Coarse graph + stats"
"$PGR" paf graph "$WORK/all.paf" -o "$WORK/pangenome.gfa" >/dev/null 2>&1
"$PGR" paf stat "$WORK/all.paf" -o "$WORK/pangenome.stat" >/dev/null 2>&1
grep -q "^segments" "$WORK/pangenome.stat"
grep -q "^paths" "$WORK/pangenome.stat"

echo "PASS: 3-genome pangenome pipeline (FastGA -> chainnet --syn -> maf-to-paf -> index -> query -> graph -> stat)."
