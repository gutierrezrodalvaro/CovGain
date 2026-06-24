#!/bin/bash
set -e
set -u
set -o pipefail

# AGR
# 24/11/2025

# This script align the trimed 15x dlabrax seqs
# This script is adapted from the ones of Nuria Sanchez Baizan (Ifremer)

# The first time you use the reference genome you need to bisulphite convert it

/opt/software/Bismark-0.24.2/bismark_genome_preparation --verbose --path_to_aligner /opt/software/bowtie2-2.4.3-linux-x86_64 --parallel 4 --genomic_composition ../../genomic_ref_data

# 1) Set directories

SEABASS_GENOME="../../genomic_ref_data"
input_dir="../../15x/trimmed_15x"
ALIGNMENT="../../15x/aligned_15x"

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
        echo "Processing $base_name"

        # Run Bismark alignment
        /opt/software/Bismark-0.24.2/bismark -q \
            --multicore 12 \
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
        echo "Warning: Corresponding read 2 file not found for $r1_file"
    fi
done

