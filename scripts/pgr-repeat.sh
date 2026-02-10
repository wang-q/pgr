#!/bin/bash

# pgr-cactus-mask.sh
# An example pipeline demonstrating the pgr-based repeat masking workflow (Cactus-style).
#
# Workflow:
# 1. Windowing: Split genome into overlapping fragments
# 2. Alignment: Self-align fragments using lastz
# 3. Conversion: LAV -> PSL -> Lift to original coordinates
# 4. Coverage: Calculate depth from alignments
# 5. Output: JSON file with high-depth regions

set -e

# Check dependencies
for cmd in pgr spanr; do
    if ! command -v $cmd &> /dev/null; then
        echo "Error: '$cmd' is not installed or not in PATH."
        exit 1
    fi
done

# Usage
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <input.fa> [work_dir] [threads]"
    exit 1
fi

INPUT_FA=$1
WORK_DIR=${2:-"masking_work_dir"}
THREADS=${3:-4}

OUTPUT_JSON="$WORK_DIR/mask_regions.json"

# Create output directory
mkdir -p "$WORK_DIR"
echo "Working directory: $WORK_DIR"

# Cleanup on exit
# trap "rm -rf $WORK_DIR" EXIT

echo "==> Step 1: Windowing (Fragmentation)"
# Split input into overlapping windows (200bp window, 100bp step = 2x coverage)
# --chunk 10000 ensures low memory usage for large genomes and enables parallel alignment
# Each chunk covers approx 1Mb genomic region (10,000 records * 100bp step)
# Fragments are output to a dedicated directory for clean query input
mkdir -p "$WORK_DIR/fragments"
pgr fa window "$INPUT_FA" -l 200 -s 100 --chunk 10000 -o "$WORK_DIR/fragments/fragments.fa"

echo "==> Step 2: Alignment (Genome vs Fragments)"
# Split genome into chromosomes (supports .gz automatically)
# This allows fine-grained parallelization and avoids single large file issues
echo "    Splitting genome into sequences..."
mkdir -p "$WORK_DIR/genome"
pgr fa split name "$INPUT_FA" -o "$WORK_DIR/genome"

# Run lastz alignment with preset parameters (Human vs Chimp set01 is a good default)
# Target is the directory of chromosomes, Query is the directory of fragments
pgr lav lastz "$WORK_DIR/genome" "$WORK_DIR/fragments" \
    --preset set01 \
    --parallel "$THREADS" \
    -o "$WORK_DIR/lastz_out"

echo "==> Step 3: Format Conversion & Coordinate Lifting"
# Convert LAV to PSL individually to avoid large file issues
: > "$WORK_DIR/fragments.psl"
for lav_file in "$WORK_DIR/lastz_out"/*.lav; do
    pgr lav to-psl "$lav_file" >> "$WORK_DIR/fragments.psl"
done

# Generate chrom.sizes for the input file (needed for lifting)
pgr fa size "$INPUT_FA" > "$WORK_DIR/chrom.sizes"

# Lift coordinates from fragments back to original genome
# --q-sizes is required to handle negative strand coordinates correctly
pgr psl lift "$WORK_DIR/fragments.psl" --q-sizes "$WORK_DIR/chrom.sizes" -o "$WORK_DIR/lifted.psl"

echo "==> Step 4: Range Extraction"
# Convert PSL alignments to depth coverage ranges (.rg)
# This extracts the query coordinates from the alignments
pgr psl to-range "$WORK_DIR/lifted.psl" > "$WORK_DIR/coverage.rg"

echo "==> Step 5: Depth Calculation"
# Use spanr to identify regions with depth >= 4
# Note: Since we used 2x coverage fragmentation (50% overlap), the baseline depth is 2.
# Single copy = depth 2. Two copies = depth 4 (2 from self + 2 from paralog).
# A threshold of 4 strictly captures regions with >= 2 copies.
spanr coverage "$WORK_DIR/coverage.rg" -m 4 > "$OUTPUT_JSON"

echo "Success! Mask regions saved to $OUTPUT_JSON"
