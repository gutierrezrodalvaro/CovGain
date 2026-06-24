#!/bin/bash
set -e
set -u
set -o pipefail

# AGR
# 19/11/2025

# Trimming and filtering using Trimgalore! of the raw sequencing data of the 15x dlabrax
# This script is adapted from the ones of Nuria Sanchez Baizan (Ifremer)

src1="../../15x/raw_15x"
outputdir1="../../15x/trimmed_15x"

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
    
# Fastqc of the trimmed sequencing data 

src2="../../15x/trimmed_15x"
outputdir2="../../15x/fastQC_reports/trimmed_15x_reports"

/opt/software/FastQC/fastqc -t 4 -o "$outputdir2" \
	"$src2"/*fq.gz


