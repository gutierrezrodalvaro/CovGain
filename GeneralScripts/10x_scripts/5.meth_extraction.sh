#!/bin/bash
set -e
set -u
set -o pipefail

# AGR
# 20/02/2025

# Methylation extraction

# 1) Set directories
SEABASS_GENOME="../../genomic_ref_data"
input_dir="../../10x/deduplicated_10x"
EXTRACT_PATH="../../10x/meth_extracted_10x"

mkdir -p "$EXTRACT_PATH"

# Run bismark_methylation_extractor with the --bam option and specify the output directory
        /opt/software/Bismark-0.24.2/bismark_methylation_extractor "$input_dir"/*.deduplicated.bam \
                              --paired-end --comprehensive \
                              --cytosine_report \
                              --cutoff 1 --multicore 12 \
			      --samtools_path /opt/software/samtools-1.11 \
                              --genome_folder "$SEABASS_GENOME" \
                              --output "$EXTRACT_PATH"

mkdir -p "$EXTRACT_PATH/CpG_reports"
mv "$EXTRACT_PATH"/*CpG_report.txt "$EXTRACT_PATH"/CpG_reports
