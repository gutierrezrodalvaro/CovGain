#!/bin/bash
set -euo pipefail

# This script sorts and index the 5x deduplicated bam files using samtools

# Paths
PROJECT_ROOT="${EMSeq_cov:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"

SAMTOOLS="/opt/software/samtools-1.11/samtools"
INPUT_DIR="$PROJECT_ROOT/10x/deduplicated_10x"
OUTPUT_DIR="$PROJECT_ROOT/10x/dedup_sorted_indexed_10x"

mkdir -p "$OUTPUT_DIR"

# Threads
THREADS=12

# Loop through all BAM files
for bam in "$INPUT_DIR"/*.bam; do
    filename=$(basename "$bam" .bam)

    echo "Processing $filename with $THREADS threads..."

    # Sort BAM
    "$SAMTOOLS" sort -@ "$THREADS" \
        -o "$OUTPUT_DIR/${filename}.sorted.bam" \
        "$bam"

    # Index BAM
    "$SAMTOOLS" index -@ "$THREADS" "$OUTPUT_DIR/${filename}.sorted.bam"

    echo "Finished $filename"
done

echo "All BAM files sorted and indexed successfully using $THREADS threads."

