#!/bin/bash
set -euo pipefail

# AGR
# 07/02/2026

# Trimming and filtering using Trimgalore! of the raw sequencing data of the 2x dlabrax
# This script is adapted from the ones of Nuria Sanchez Baizan (Ifremer)

src1="../../2x/raw_fastq_2x_subsampled"
outputdir1="../../2x/trimmed_2x"

# Create the outputdir1 if it does not exists
mkdir -p "$outputdir1"

for r1 in "$src1"/*_R1.fastq.gz; do
	r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"
    	
	echo "Processing pair:"
	echo "  R1 = $r1"
	echo "  R2 = $r2"
    
	/opt/software/TrimGalore-0.6.10/trim_galore -q 30 --illumina --phred33 -j 4 --paired \
		--clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10 \
		--path_to_cutadapt /home/alvaro/anaconda3/bin/cutadapt \
		"$r1" "$r2" \
		-o "$outputdir1"
done
    
# I do not run a FastQC after trimming as this is a subset and the total quality was checked for the 15x dataset

# Unification of the trimming reports

report_dir="$outputdir1"
output_file=""$report_dir"/trimming_summary_2x.txt"

# Clear output file (create if it does not exists, clear if it exists)
: > "$output_file"

# Use find to locate all trimming report files in the current directory
find "$report_dir" -maxdepth 1 -name "*trimming_report.txt" | while read -r report; do
    # Extract the filename
    filename=$(basename "$report")

    # Extract Quality-trimmed information
    quality_trimmed=$(grep "Quality-trimmed:" "$report" | awk '{print $2, $3, $4}')

    # Extract Total written information
    total_written=$(grep "Total written (filtered):" "$report" | awk '{print $4, $5, $6}')

    # Write the extracted information to the output file
    echo "File: $filename" >> "$output_file"
    echo "Quality-trimmed: $quality_trimmed" >> "$output_file"
    echo "Total written (filtered): $total_written" >> "$output_file"
    echo "" >> "$output_file"
done
echo "Summary has been written to $output_file"


