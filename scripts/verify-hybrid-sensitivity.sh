#!/usr/bin/env bash

# verify-hybrid-sensitivity.sh
# Sensitivity evaluation for `pgr align hybrid` following fastga.md §12.1
# (adapted from the FastGA paper §5.1): simulate two genomes of 10 kb blocks,
# each block = a similar region (length L at divergence d) + random filler,
# block order shuffled identically so no long-range colinear homology exists.
# A target region is "recovered" iff it is covered >= 95% on BOTH genomes.
#
# Produces a per-(L, d) recovery matrix for pgi-only / hybrid / lastz-only,
# plus the false-positive aligned-base fraction on the A genome, then asserts
# the design's qualitative claims:
#   * hybrid sensitivity >= pgi sensitivity (strictly better on >= 1 cell)
#   * hybrid sensitivity ~= lastz sensitivity (within a small cell tolerance)
#   * all three false-positive fractions < 1%
#
# Requires:
#   - pgr binary (prefers target/release/pgr, falls back to target/debug/pgr;
#     use PGR=... to override)
#   - lastz in PATH
#   - python3 (analysis + data generation)
#
# Usage:
#   scripts/verify-hybrid-sensitivity.sh

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
command -v python3 >/dev/null || { echo "Error: python3 not in PATH." >&2; exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# --- 1. Generate simulated genomes (scaled-down §12.1: 6 Mb, 20 repeats) ---
python3 - "$WORK" <<'PY'
import random, sys, os
OUT = sys.argv[1]
LENGTHS  = [100, 200, 500, 1000, 2000, 5000]
DIVERGEN = [1, 10, 20, 30, 40]
REPEATS  = 20
BLOCK    = 10000
def rand_seq(n, rng):
    return ''.join(rng.choice('ACGT') for _ in range(n))
def diverge(seq, rate, rng):
    out = []
    for b in seq:
        r = rng.random()
        if r < 0.80 * rate:      out.append(rng.choice('ACGT')) # substitution
        elif r < 0.90 * rate:    continue                       # deletion
        elif r < rate:           out.append(rng.choice('ACGT')); out.append(b) # insertion
        else:                    out.append(b)
    return ''.join(out)
rng = random.Random(1)
blocks = []
for d in DIVERGEN:
    for L in LENGTHS:
        for rep in range(REPEATS):
            t_sim = rand_seq(L, rng)
            q_sim = diverge(t_sim, d / 100.0, rng)
            blocks.append((t_sim + rand_seq(BLOCK - L, rng),
                           q_sim + rand_seq(BLOCK - L, rng),
                           f"L{L}_d{d}_{rep}", L, d))
rng.shuffle(blocks)   # same order drives both genomes -> homolog at same start
t_seq = ''.join(b[0] for b in blocks)
q_seq = ''.join(b[1] for b in blocks)
targets = []
pos = 0
for (t_block, q_block, name, L, d) in blocks:
    targets.append((name, pos, pos + L, L, d)); pos += BLOCK
os.makedirs(OUT, exist_ok=True)
with open(f"{OUT}/A.fa", "w") as f: f.write(">A\n" + t_seq + "\n")
with open(f"{OUT}/B.fa", "w") as f: f.write(">B\n" + q_seq + "\n")
with open(f"{OUT}/targets.tsv", "w") as f:
    for (name, s, e, L, d) in targets: f.write(f"{name}\t{s}\t{e}\t{L}\t{d}\n")
PY

# --- 2. Three engines ---
echo "==> pgi-only"
"$PGR" align pgi "$WORK/A.fa" "$WORK/B.fa" -o "$WORK/pgi.psl" --parallel 8 >/dev/null 2>&1
echo "==> hybrid (pgi anchors + LASTZ gap fill)"
"$PGR" align hybrid "$WORK/A.fa" "$WORK/B.fa" --pgi-psl "$WORK/pgi.psl" -o "$WORK/hybrid.psl" --parallel 8 >/dev/null 2>&1
echo "==> lastz-only"
"$PGR" align lastz "$WORK/A.fa" "$WORK/B.fa" -o "$WORK/lastz_out" --parallel 8 >/dev/null 2>&1
"$PGR" lav to-psl "$WORK"/lastz_out/*.lav -o "$WORK/lastz.psl" >/dev/null 2>&1

# --- 3. Analyze + assert ---
python3 - "$WORK" <<'PY'
import sys, os
BASE = sys.argv[1]
REPEATS, LENGTHS, DIVERGEN = 20, [100,200,500,1000,2000,5000], [1,10,20,30,40]
targets = {}
with open(f"{BASE}/targets.tsv") as f:
    for line in f:
        name, s, e, L, d = line.rstrip("\n").split("\t")
        targets[name] = (int(s), int(e), int(L), int(d))

def load_psl(path):
    recs = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            fld = line.split(); qn, tn = fld[9], fld[13]
            qs, qe, ts, te = int(fld[11]), int(fld[12]), int(fld[15]), int(fld[16])
            if tn == "A": recs.append(("A", ts, te))
            if qn == "B": recs.append(("B", qs, qe))
    return recs

def cov_frac(intervals, s, e):
    if e <= s: return 0.0
    cov = [0]*(e-s)
    for _, st, en in intervals:
        for i in range(max(st,s), min(en,e)):
            if i >= s and i < e: cov[i-s] = 1
    return sum(cov)/(e-s)

def evaluate(path):
    recs = load_psl(path)
    a = [r for r in recs if r[0]=="A"]; b = [r for r in recs if r[0]=="B"]
    res = {(L,d):0 for L in LENGTHS for d in DIVERGEN}
    for name,(s,e,L,d) in targets.items():
        if cov_frac(a,s,e) >= 0.95 and cov_frac(b,s,e) >= 0.95: res[(L,d)] += 1
    return res

def spec(path):
    bits = bytearray(6000000)
    for (s, e, L, d) in targets.values(): bits[s:e] = b"\x01" * (e - s)
    tot = fp = 0
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            fld = line.split()
            if fld[13] != "A": continue
            ts, te = int(fld[15]), int(fld[16])
            tot += max(0, te-ts); fp += sum(1 for i in range(ts,te) if bits[i]==0)
    return (fp/tot*100) if tot else 0.0

eng = {k: evaluate(f"{BASE}/{k}.psl") for k in ("pgi","hybrid","lastz")}
fs = {k: spec(f"{BASE}/{k}.psl") for k in ("pgi","hybrid","lastz")}

print("recovered/20 per (L, d); cells shown as hybrid/pgi/lastz:")
print("L\\d    " + " ".join(f"{d:>12}%" for d in DIVERGEN))
for L in LENGTHS:
    print(f"{L:<5} " + " ".join(
        f"{eng['hybrid'][(L,d)]:>3}/{eng['pgi'][(L,d)]:>2}/{eng['lastz'][(L,d)]:>2}"
        for d in DIVERGEN))
print("false-positive aligned-base %% on A: pgi=%.3f hybrid=%.3f lastz=%.3f"
      % (fs['pgi'], fs['hybrid'], fs['lastz']))

tp = sum(eng['hybrid'].values()); pg = sum(eng['pgi'].values()); lz = sum(eng['lastz'].values())
assert tp >= pg, f"hybrid ({tp}) not >= pgi ({pg})"
diff = sum(abs(eng['hybrid'][k]-eng['lastz'][k]) for k in eng['hybrid'].keys())
assert tp > pg, "hybrid must strictly beat pgi on sensitivity"
assert diff <= 5, f"hybrid vs lastz cell diff {diff} too large"
assert fs['pgi'] < 1.0 and fs['hybrid'] < 1.0 and fs['lastz'] < 1.0, "false-positive too high"
print(f"PASS: hybrid={tp}/600 >= pgi={pg}/600, ~lastz={lz}/600 (cell diff {diff}); FP all < 1%")
PY