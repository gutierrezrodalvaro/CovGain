#!/bin/bash
# Remove PCR duplicates and optical duplicates that can introduce bias in downstream analyses.
# Reduce false positive variant calls and improve accuracy of quantitative analyses.
# EM-seq libraries can incorporate UMIs to more accurately identify and remove PCR duplicates

# Script adapted from Nuria SB scripts 
# Alvaro
# 20/02/2026

# Project root detection     
PROJECT_ROOT="${EMSeq_cov:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"

src="$PROJECT_ROOT/10x/aligned_10x"
output="$PROJECT_ROOT/10x/deduplicated_10x"

mkdir -p "$output"

#Run deduplication for all reads
for bam_file in "$src"/*.bam; do
    # Extract the base name of the BAM file (without the directory and extension)
    base_name=$(basename "$bam_file" .bam)
    
    echo "Processing: $base_name"

    # Run deduplicate_bismark with the --barcode option and specify the output directory
    /opt/software/Bismark-0.24.2/deduplicate_bismark --samtools_path /opt/software/samtools-1.11 --output_dir "$output" --paired --bam "$bam_file"
    
done

echo "Deduplication of the 10x bam files done!"
