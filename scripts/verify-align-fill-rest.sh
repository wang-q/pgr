#!/usr/bin/env bash

# verify-align-fill-rest.sh
# Real-data comparison on MG1655 x Sakai: pgi-only / fill / rest / lastz-only.
# For each engine: alignment -> PSL stats + `pgr pl chainnet --syn` ->
# target-side union coverage (design/pgi-lastz-hybrid.md §6).
#
# Requires:
#   - pgr binary (prefers target/release/pgr, falls back to target/debug/pgr;
#     use PGR=... to override)
#   - lastz in PATH
#
# Usage:
#   scripts/verify-align-fill-rest.sh

set -euo pipefail

cd "$(dirname "$0")/.."
PGR="${PGR:-}"
if [ -z "$PGR" ] && [ -x "$PWD/target/release/pgr" ]; then
    PGR="$PWD/target/release/pgr"
else
    PGR="${PGR:-$PWD/target/debug/pgr}"
fi
[ -x "$PGR" ] || { echo "Error: pgr binary not found at $PGR (run 'cargo build' first)." >&2; exit 1; }
command -v lastz >/dev/null || { echo "Error: lastz not in PATH." >&2; exit 1; }

TARGET="$PWD/tests/genome/mg1655.fa.gz"
QUERY="$PWD/tests/genome/sakai.fa.gz"
SIZES="$PWD/tests/genome/mg1655.chr.sizes"
[ -f "$TARGET" ] && [ -f "$QUERY" ] && [ -f "$SIZES" ]

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# --- 0. Prepare plain-text single-sequence inputs (lastz needs these) ---
gzip -dc "$TARGET" > "$WORK/mg1655.fa"
mkdir -p "$WORK/sakai_split" "$WORK/lastz_raw" "$WORK/lastz_psl"
"$PGR" fa split name "$QUERY" -o "$WORK/sakai_split" >/dev/null 2>&1

# coverage of a PSL on the target genome (union, bp + fraction)
coverage() {
    local psl=$1
    "$PGR" psl to-rg --target-coords "$psl" | "$PGR" rg cover stdin -o "$WORK/union.json"
    "$PGR" runlist stat "$SIZES" "$WORK/union.json" | tail -1
}

# target-side union coverage (%) from a `chainnet --syn` MAF (the final
# syntenic-filtered result). The MAF src names are `<strain>.<contig>`
# (e.g. mg1655.NC_000913); pass the strain prefix ("mg1655.").
maf_coverage() {
    local cn_dir=$1 prefix=$2
    local maf
    maf=$(ls "$cn_dir"/*.maf 2>/dev/null | head -1)
    [ -n "$maf" ] || { echo "no-maf"; return; }
    local genome_len
    genome_len=$(awk -v p="$prefix" '$2 ~ "^"p {print $6; exit}' "$maf")
    grep '^s ' "$maf" \
        | awk -v p="$prefix" '$2 ~ "^"p {print $3, $3 + $4}' \
        | sort -n \
        | awk -v L="$genome_len" '
            NR == 1 { s = $1; e = $2; next }
            $1 <= e { if ($2 > e) e = $2; next }
            { cov += e - s; s = $1; e = $2 }
            END { cov += e - s; printf "%.4f (%d)", cov / L * 100, cov }
          '
}

TARGET_PREFIX="$(basename "$TARGET" .fa.gz)."

echo "==> Engine 1/5: pgi-only"
t0=$(date +%s%N)
"$PGR" align pgi "$TARGET" "$QUERY" -o "$WORK/pgi.psl" >/dev/null 2>&1
t1=$(date +%s%N)
"$PGR" pl chainnet --syn "$TARGET" "$QUERY" "$WORK/pgi.psl" -o "$WORK/cn_pgi" >/dev/null 2>&1
echo "runtime (align): $(( (t1 - t0) / 1000000 )) ms"
echo "records: $(wc -l < "$WORK/pgi.psl")"
echo "raw span bp: $(awk '{s+=$17-$16} END{print s}' "$WORK/pgi.psl")"
echo "maf blocks: $(grep -c '^a' "$WORK/cn_pgi"/*.maf)"
echo "maf cov: $(maf_coverage "$WORK/cn_pgi" "$TARGET_PREFIX")"
coverage "$WORK/pgi.psl"

echo "==> Engine 2/5: fill (2D anchor-gap fill)"
t0=$(date +%s%N)
"$PGR" align fill "$TARGET" "$QUERY" -o "$WORK/fill.psl" >/dev/null 2>&1
t1=$(date +%s%N)
"$PGR" pl chainnet --syn "$TARGET" "$QUERY" "$WORK/fill.psl" -o "$WORK/cn_fill" >/dev/null 2>&1
echo "runtime (align): $(( (t1 - t0) / 1000000 )) ms"
echo "records: $(wc -l < "$WORK/fill.psl")"
echo "raw span bp: $(awk '{s+=$17-$16} END{print s}' "$WORK/fill.psl")"
echo "maf blocks: $(grep -c '^a' "$WORK/cn_fill"/*.maf)"
echo "maf cov: $(maf_coverage "$WORK/cn_fill" "$TARGET_PREFIX")"
coverage "$WORK/fill.psl"

echo "==> Engine 3/5: rest (whole-genome complement fill)"
t0=$(date +%s%N)
"$PGR" align rest "$TARGET" "$QUERY" -o "$WORK/rest.psl" >/dev/null 2>&1
t1=$(date +%s%N)
"$PGR" pl chainnet --syn "$TARGET" "$QUERY" "$WORK/rest.psl" -o "$WORK/cn_rest" >/dev/null 2>&1
echo "runtime (align): $(( (t1 - t0) / 1000000 )) ms"
echo "records: $(wc -l < "$WORK/rest.psl")"
echo "raw span bp: $(awk '{s+=$17-$16} END{print s}' "$WORK/rest.psl")"
echo "maf blocks: $(grep -c '^a' "$WORK/cn_rest"/*.maf)"
echo "maf cov: $(maf_coverage "$WORK/cn_rest" "$TARGET_PREFIX")"
coverage "$WORK/rest.psl"

echo "==> Engine 4/5: fill + rest combined (PSL cat)"
cat "$WORK/fill.psl" "$WORK/rest.psl" > "$WORK/fill_rest.psl"
echo "records: $(wc -l < "$WORK/fill_rest.psl")"
"$PGR" pl chainnet --syn "$TARGET" "$QUERY" "$WORK/fill_rest.psl" -o "$WORK/cn_fr" >/dev/null 2>&1
echo "maf blocks: $(grep -c '^a' "$WORK/cn_fr"/*.maf)"
echo "maf cov: $(maf_coverage "$WORK/cn_fr" "$TARGET_PREFIX")"
coverage "$WORK/fill_rest.psl"

echo "==> Engine 5/5: lastz-only"
t0=$(date +%s%N)
for c in "$WORK"/sakai_split/*.fa; do
    name=$(basename "$c" .fa)
    "$PGR" align lastz "$WORK/mg1655.fa" "$c" --parallel 8 -o "$WORK/lastz_raw/$name" >/dev/null 2>&1
done
t1=$(date +%s%N)
for lav in "$WORK"/lastz_raw/*/*.lav; do
    name=$(basename "$lav" .lav | sed 's/\[mg1655\]vs\[//; s/\]$//')
    "$PGR" lav to-psl "$lav" -o "$WORK/lastz_psl/$name.psl" >/dev/null 2>&1
done
cat "$WORK"/lastz_psl/*.psl > "$WORK/lastz.psl"
"$PGR" pl chainnet --syn "$TARGET" "$QUERY" "$WORK/lastz.psl" -o "$WORK/cn_lastz" >/dev/null 2>&1
echo "runtime (align+convert): $(( (t1 - t0) / 1000000 )) ms"
echo "records: $(wc -l < "$WORK/lastz.psl")"
echo "raw span bp: $(awk '{s+=$17-$16} END{print s}' "$WORK/lastz.psl")"
echo "maf blocks: $(grep -c '^a' "$WORK/cn_lastz"/*.maf)"
echo "maf cov: $(maf_coverage "$WORK/cn_lastz" "$TARGET_PREFIX")"
coverage "$WORK/lastz.psl"
