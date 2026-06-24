#!/bin/bash
set -euo pipefail

# AGR
# 11/02/2026

# This script align the trimed 5x dlabrax seqs
# The first time you use the reference genome you need to bisulphite convert it

# /opt/software/Bismark-0.24.2/bismark_genome_preparation --verbose --path_to_aligner /opt/software/bowtie2-2.4.3-linux-x86_64 --parallel 4 --genomic_composition ../../genomic_ref_data

# 1) Set directories

# Project root detection
PROJECT_ROOT="${EMSeq_cov:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"

SEABASS_GENOME="$PROJECT_ROOT/genomic_ref_data"
input_dir="$PROJECT_ROOT/5x/trimmed_5x"
ALIGNMENT="$PROJECT_ROOT/5x/aligned_5x"

# Create ALIGNMENT directory if it doesn't exist
mkdir -p "$ALIGNMENT"

# 2) Loop through files in TRIMMED folder
for r1_file in "$input_dir"/*_R1_val_1.fq.gz; do
    # Extract the base name of the file (everything before _val_1.fq.gz)
    base_name=$(basename "$r1_file" _R1_val_1.fq.gz)

    # Construct the name of the corresponding read 2 file
    r2_file="${input_dir}/${base_name}_R2_val_2.fq.gz"

    # Check if the read 2 file exists
    if [ -f "$r2_file" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing $base_name"

        # Run Bismark alignment
        /opt/software/Bismark-0.24.2/bismark -q \
            --multicore 8 \
            -N 1 \
            --genome "$SEABASS_GENOME" \
            --path_to_bowtie2 /opt/software/bowtie2-2.4.3-linux-x86_64 \
            --samtools_path /opt/software/samtools-1.11 \
            -1 "$r1_file" \
            -2 "$r2_file" \
            --score_min L,0,-0.4 \
            --rg_tag \
            -o "$ALIGNMENT"
    else
        echo "Warning: Corresponding read 2 file not found for $base_name"
    fi
done
