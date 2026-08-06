#!/bin/bash

# rm-gold-compare.sh
# RepeatMasker gold-standard comparison against pgr e-kmer / e-align runlists.
#
# Workflow:
# 1. Run RepeatMasker (-lib <library>, rmblast engine, -pa N) on every genome.
# 2. Convert .out to GFF3 (full, IS-only, strict >=50 bp & >=70% identity)
#    and then to runlist JSON via `pgr gff runlist`.
# 3. Cross our runlists with the RM runlists via `pgr runlist statop` and
#    compute RM-only / our-only diff stats via `pgr runlist compare`.
# 4. Write per-genome TSVs and print an aggregate summary.
#
# Notes:
# * The library must be plain FASTA (RepeatMasker -lib does not accept .gz).
# * `rm-is` filters out Simple_repeat/Low_complexity/Satellite rows;
#   `rm-strict` additionally requires span >= 50 bp and identity >= 70%
#   (100 - div - del - ins), matching e-align's min-len / min-identity.
#
# Usage:
#   scripts/rm-gold-compare.sh <genome_dir> <our_json_dir> <sizes_dir> \
#       <library.fa> <out_dir> [rm_dir] [threads] [pgr]

set -euo pipefail

if [ "$#" -lt 5 ]; then
    echo "Usage: $0 <genome_dir> <our_json_dir> <sizes_dir> <library.fa> <out_dir> [rm_dir] [threads] [pgr]" >&2
    exit 1
fi

GENOME_DIR=$1
OUR_JSON_DIR=$2
SIZES_DIR=$3
LIB=$4
OUT=$5
RM_DIR=${6:-/home/wangq/Scripts/pgr/RepeatMasker}
THREADS=${7:-8}
if [ -n "${8:-}" ]; then
    PGR=$8
elif [ -x "$(dirname "$0")/../target/release/pgr" ]; then
    PGR="$(dirname "$0")/../target/release/pgr"
else
    PGR="$(dirname "$0")/../target/debug/pgr"
fi

RM="$RM_DIR/RepeatMasker"
RM2GFF="$RM_DIR/util/rmOutToGFF3.pl"

METHODS="e-kmer.tncentral e-align.tncentral e-align-pgi-k31 e-align-lastz.tncentral"

# Check dependencies.
for dep in "$RM" "$RM2GFF"; do
    [ -f "$dep" ] || { echo "Error: '$dep' not found" >&2; exit 1; }
done
if ! command -v "$PGR" >/dev/null 2>&1 && [ ! -x "$PGR" ]; then
    echo "Error: '$PGR' not found" >&2
    exit 1
fi

mkdir -p "$OUT/gff" "$OUT/json"
printf 'genome\tmethod\tour_bp\trm_is_bp\tinter_bp\trm_cov_by_ours\tours_cov_by_rm\n' > "$OUT/rm_vs_ours.tsv"
printf 'genome\tmethod\tour_bp\trm_strict_bp\tinter_bp\trm_cov_by_ours\tours_cov_by_rm\n' > "$OUT/rm_strict_vs_ours.tsv"
printf 'genome\tmethod\trm_only_bp\trm_only_n\tavg\tour_only_bp\tour_only_n\tavg\n' > "$OUT/diff_stats.tsv"

for gfa in "$GENOME_DIR"/*.fa.gz; do
    g=$(basename "$gfa" .fa.gz)
    mkdir -p "$OUT/$g"

    echo "==> RepeatMasker: $g"
    "$RM" -lib "$LIB" -pa "$THREADS" -e rmblast -dir "$OUT/$g" "$gfa" > "$OUT/$g.log" 2>&1
    out="$OUT/$g/$g.fa.out"

    # Full output, IS-only, and strict (>=50 bp, >=70% identity).
    perl "$RM2GFF" "$out" > "$OUT/gff/$g.rm.gff"
    awk '!/^ *(SW|score|$)/ && $11 !~ /Simple_repeat|Low_complexity|Satellite/' "$out" > "$out.is"
    perl "$RM2GFF" "$out.is" > "$OUT/gff/$g.rm-is.gff"
    awk '!/^ *(SW|score|$)/ && $11 !~ /Simple_repeat|Low_complexity|Satellite/ &&
         ($7-$6+1) >= 50 && (100-$2-$3-$4) >= 70' "$out" > "$out.strict"
    perl "$RM2GFF" "$out.strict" > "$OUT/gff/$g.rm-strict.gff"

    for kind in rm rm-is rm-strict; do
        "$PGR" gff runlist "$OUT/gff/$g.$kind.gff" -o "$OUT/json/$g.$kind.json"
    done

    for m in $METHODS; do
        ours="$OUR_JSON_DIR/$g.$m.json"
        [ -f "$ours" ] || { echo "skip missing $ours" >&2; continue; }

        for kind in is strict; do
            line=$("$PGR" runlist statop "$SIZES_DIR/$g.sizes" "$ours" "$OUT/json/$g.rm-$kind.json" --all | tail -1)
            our_bp=$(echo "$line" | cut -f2)
            rm_bp=$(echo "$line" | cut -f3)
            inter=$(echo "$line" | cut -f4)
            c2=$(echo "$line" | cut -f6)
            oc=$(awk -v a="$our_bp" -v i="$inter" 'BEGIN{printf "%.4f", i/a}')
            if [ "$kind" = is ]; then
                printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$g" "$m" "$our_bp" "$rm_bp" "$inter" "$c2" "$oc" >> "$OUT/rm_vs_ours.tsv"
            else
                printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$g" "$m" "$our_bp" "$rm_bp" "$inter" "$c2" "$oc" >> "$OUT/rm_strict_vs_ours.tsv"
            fi
        done

        # Diff stats vs rm-strict: RM-only and our-only fragments.
        "$PGR" runlist compare "$OUT/json/$g.rm-strict.json" "$ours" --op diff -o "$OUT/tmp_diff.json"
        rl=$("$PGR" runlist stat "$SIZES_DIR/$g.sizes" "$OUT/tmp_diff.json" --all | tail -1)
        rm_only_bp=$(echo "$rl" | cut -f2)
        rm_only_n=$("$PGR" runlist convert "$OUT/tmp_diff.json" | wc -l)

        "$PGR" runlist compare "$ours" "$OUT/json/$g.rm-strict.json" --op diff -o "$OUT/tmp_diff2.json"
        ol=$("$PGR" runlist stat "$SIZES_DIR/$g.sizes" "$OUT/tmp_diff2.json" --all | tail -1)
        our_only_bp=$(echo "$ol" | cut -f2)
        our_only_n=$("$PGR" runlist convert "$OUT/tmp_diff2.json" | wc -l)

        ravg=$(awk -v b="$rm_only_bp" -v n="$rm_only_n" 'BEGIN{if(n>0) printf "%.0f", b/n; else print 0}')
        oavg=$(awk -v b="$our_only_bp" -v n="$our_only_n" 'BEGIN{if(n>0) printf "%.0f", b/n; else print 0}')
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$g" "$m" "$rm_only_bp" "$rm_only_n" "$ravg" "$our_only_bp" "$our_only_n" "$oavg" >> "$OUT/diff_stats.tsv"
    done
done

echo "==> Aggregate (10 genomes)"
awk -F'\t' 'NR>1{
    our[$2]+=$3; rm[$2]+=$4; inter[$2]+=$5
} END{
    for (m in our)
        printf "%-24s our=%d\trm=%d\tinter=%d\trm_cov_by_ours=%.1f%%\tours_cov_by_rm=%.1f%%\n",
            m, our[m], rm[m], inter[m], 100*inter[m]/rm[m], 100*inter[m]/our[m]
}' "$OUT/rm_strict_vs_ours.tsv"
