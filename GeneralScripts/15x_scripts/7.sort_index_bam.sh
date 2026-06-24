#!/bin/bash
set -euo pipefail

# This script sorts and index the deduplicated bam files using samtools

# Paths
SAMTOOLS="/opt/software/samtools-1.11/samtools"
INPUT_DIR="../../15x/deduplicated"
OUTPUT_DIR="../../15x/dedup_sorted_indexed_15x"

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

