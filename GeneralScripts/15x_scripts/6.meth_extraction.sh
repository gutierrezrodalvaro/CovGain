#!/bin/bash
set -e
set -u
set -o pipefail

# AGR
# 04/12/2025

# Methylation extraction

# 1) Set directories
SEABASS_GENOME="../../genomic_ref_data"
input_dir="/data2/gonzalo2/alvaro_data2/EMSeq_cov/15x/deduplicated"
EXTRACT_PATH="../../15x/meth_extracted_15x"

mkdir -p "$EXTRACT_PATH"

# Run bismark_methylation_extractor with the --bam option and specify the output directory
/opt/software/Bismark-0.24.2/bismark_methylation_extractor "$input_dir"/*.bam \
	--paired-end --comprehensive \
        --cytosine_report \
        --cutoff 1 --multicore 4 \
	--samtools_path /opt/software/samtools-1.11 \
        --genome_folder "$SEABASS_GENOME" \
        --output "$EXTRACT_PATH"

