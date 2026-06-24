#!/bin/bash
set -e
set -u
set -o pipefail

# AGR
# 19/11/2025

# Fastqc of the raw sequencing data files 15x from Nuria SB

src="../../15x/raw_15x"
outputdir="../../15x/fastQC_reports/raw_15x_reports/"

mkdir -p "$outputdir"

# Do the fastQC for each fastq file in src with 4 threads
/opt/software/FastQC/fastqc -t 4 -o "$outputdir" \
        "$src"/*fastq.gz

# Do the multiQC report combining all the fastqc reports in the directory
/home/alvaro/anaconda3/bin/multiqc "$outputdir" \
    -o "$outputdir" \
    -i "raw15x_multi-QC"

# Print in the console the number of reads per file
cd "$src"
for i in *.fastq.gz; do echo -n "$i "; echo $(zcat $i | wc -l)/4 | bc; done
