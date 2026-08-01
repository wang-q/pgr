#!/usr/bin/env bash

# verify-pangenome.sh
# Ten-E. coli pangenome test: FastGA pairwise PSL -> pgr pl chainnet
# (chain -> net -> axt -> maf, the mandatory syntenic route) -> maf to-paf ->
# multi-file index -> transitive query -> coarse graph -> stat -> downstream
# (to-maf --msa / to-vcf / to-fas --msa / to-gfa) on one 10 kb MG1655 region.
# Requires:
#   - pgr binary (prefers target/release/pgr for the POA-heavy downstream
#     steps, falls back to target/debug/pgr; use PGR=... to override)
#   - FastGA in PATH (pangenome-route upstream aligner, see notes/ecoli-cohort.md)
#
# Usage:
#   scripts/verify-pangenome.sh

set -euo pipefail

cd "$(dirname "$0")/.."
PGR="${PGR:-}"
if [ -z "$PGR" ] && [ -x "$PWD/target/release/pgr" ]; then
    PGR="$PWD/target/release/pgr"
else
    PGR="${PGR:-$PWD/target/debug/pgr}"
fi

GENOMES=(
    "$PWD/tests/genome/mg1655.fa.gz"
    "$PWD/tests/genome/sakai.fa.gz"
    "$PWD/tests/genome/se11.fa.gz"
    "$PWD/tests/genome/cft073.fa.gz"
    "$PWD/tests/genome/e2348_69.fa.gz"
    "$PWD/tests/genome/ec042.fa.gz"
    "$PWD/tests/genome/ec2011c_3493.fa.gz"
    "$PWD/tests/genome/e24377a.fa.gz"
    "$PWD/tests/genome/ec958.fa.gz"
    "$PWD/tests/genome/nissle1917.fa.gz"
)

if [ ! -x "$PGR" ]; then
    echo "Error: pgr binary not found at $PGR (run 'cargo build' first)." >&2
    exit 1
fi
command -v FastGA >/dev/null || { echo "Error: FastGA not in PATH." >&2; exit 1; }
for f in "${GENOMES[@]}"; do
    [ -f "$f" ] || { echo "Error: missing genome $f" >&2; exit 1; }
done

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

echo "==> 1-3. FastGA pairwise PSL -> pgr pl chainnet --syn -> MAF -> PAF"
shopt -s nullglob
mkdir -p "$WORK/paf"
n=${#GENOMES[@]}
p=0
for ((i = 0; i < n; i++)); do
    for ((j = i + 1; j < n; j++)); do
        a="${GENOMES[$i]}"
        b="${GENOMES[$j]}"
        pa=$(basename "$a" .fa.gz)
        pb=$(basename "$b" .fa.gz)
        p=$((p + 1))
        echo "==> [pair $p/$((n * (n - 1) / 2))] $pa x $pb"
        # Orientation (verified against the 3-genome run): FastGA(b, a) then
        # chainnet(a, b) yields PAF query=b / target=a, so MG1655 (always "a"
        # here, index 0) is the target of every pair.  The seed region is then
        # found in `trees` regardless of strand ('+' or '-'); BFS from the
        # query side relies on '+' mirrors only and would miss '-' records.
        FastGA -v -psl "$b" "$a" > "$WORK/$pa-$pb.psl" 2>/dev/null
        "$PGR" pl chainnet --syn "$a" "$b" "$WORK/$pa-$pb.psl" -o "$WORK/cn_$p" >/dev/null 2>&1
        for m in "$WORK"/cn_$p/*.maf; do
            # Pair-prefixed names: the same replicon appears as target in
            # multiple pairs (e.g. NC_004431 for cft073 x mg1655/sakai/se11),
            # so bare basenames would collide and overwrite each other.
            "$PGR" maf to-paf "$m" -o "$WORK/paf/$pa-$pb-$(basename "$m" .maf).paf" >/dev/null 2>&1
        done
    done
done
shopt -u nullglob
cat "$WORK"/paf/*.paf > "$WORK/all.paf"

echo "==> 4. Multi-file PAF index"
"$PGR" paf index "$WORK"/paf/*.paf -o "$WORK/pangenome.paf.idx" >/dev/null 2>&1

echo "==> 5. Transitive query (10-way synteny from MG1655)"
"$PGR" paf query "$WORK/pangenome.paf.idx" mg1655.NC_000913:100000-110000 \
    --transitive -o "$WORK/tri.paf" >/dev/null 2>&1
for f in "${GENOMES[@]}"; do
    strain=$(basename "$f" .fa.gz)
    grep -qF "$strain." "$WORK/tri.paf" || { echo "Error: $strain not reached by transitive query" >&2; exit 1; }
done

echo "==> 6. Coarse graph + stats"
"$PGR" paf graph "$WORK/all.paf" -o "$WORK/pangenome.gfa" >/dev/null 2>&1
"$PGR" paf stat "$WORK/all.paf" -o "$WORK/pangenome.stat" >/dev/null 2>&1
grep -q "^segments" "$WORK/pangenome.stat"
grep -q "^paths" "$WORK/pangenome.stat"

echo "==> 7. Downstream: prefixed BGZF FASTA + fasta-tsv for MSA emitters"
mkdir -p "$WORK/fas"
for f in "${GENOMES[@]}"; do
    s=$(basename "$f" .fa.gz)
    gzip -dc "$f" \
        | awk -v p="$s" '/^>/{sub(/^>/,"",$1); print ">"p"."$1; next} {print}' \
        | "$PGR" fa gz stdin -o "$WORK/fas/$s.fa.gz" >/dev/null 2>&1
done
cut -f1,6 "$WORK/all.paf" | tr '\t' '\n' | sort -u \
    | while read -r n; do printf '%s\t%s/fas/%s.fa.gz\n' "$n" "$WORK" "${n%%.*}"; done \
    > "$WORK/seqs.tsv"

echo "==> 8. to-maf --msa (10-way POA MSA)"
"$PGR" paf to-maf -f "$WORK/seqs.tsv" "$WORK/pangenome.paf.idx" \
    mg1655.NC_000913:100000-110000 --transitive --msa -o "$WORK/region.maf" >/dev/null 2>&1
[ "$(grep -c '^a' "$WORK/region.maf")" -eq 1 ] # one block per queried region
[ "$(awk '/^s/{print $2}' "$WORK/region.maf" | sort -u | wc -l)" -eq 10 ] # one s-line per genome
[ "$(awk '/^s/{c++} END{print c}' "$WORK/region.maf")" -eq 10 ]           # no duplicate names in block
for f in "${GENOMES[@]}"; do
    strain=$(basename "$f" .fa.gz)
    grep -qF "$strain." "$WORK/region.maf" \
        || { echo "Error: $strain missing from 10-way MSA" >&2; exit 1; }
done

echo "==> 9. to-vcf (10-way variants)"
"$PGR" paf to-vcf -f "$WORK/seqs.tsv" "$WORK/pangenome.paf.idx" \
    mg1655.NC_000913:100000-110000 --transitive -o "$WORK/region.vcf" >/dev/null 2>&1
[ "$(grep -vc '^#' "$WORK/region.vcf")" -gt 0 ]
[ "$(grep -v '^##' "$WORK/region.vcf" | head -n1 | awk -F'\t' '{print NF-9}')" -eq 10 ] # unique sample columns

echo "==> 10. to-fas --msa (block FA) + to-gfa (local graph)"
"$PGR" paf to-fas -f "$WORK/seqs.tsv" "$WORK/pangenome.paf.idx" \
    mg1655.NC_000913:100000-110000 --transitive --msa -o "$WORK/region.fas" >/dev/null 2>&1
[ "$(grep -c '^>' "$WORK/region.fas")" -eq 10 ]
"$PGR" paf to-gfa -f "$WORK/seqs.tsv" "$WORK/pangenome.paf.idx" \
    mg1655.NC_000913:100000-110000 --transitive -o "$WORK/region.gfa" >/dev/null 2>&1
grep -q '^S' "$WORK/region.gfa"
grep -q '^P' "$WORK/region.gfa"

echo "PASS: 10-genome pangenome pipeline (FastGA -> chainnet --syn -> maf-to-paf -> index -> query -> graph -> stat -> msa/vcf/fas/gfa)."
