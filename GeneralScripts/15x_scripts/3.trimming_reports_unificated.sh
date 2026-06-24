#!/bin/bash
set -e
set -u
set -o pipefail

# AGR
# 24/11/2025

# This script extract the percentage trimmed and generates a new summary reports of all the sequences, 15x dlabrax
# This script is adapted from the ones of Nuria Sanchez Baizan (Ifremer)

report_dir="../../15x/trimmed_15x"
output_file="../../15x/trimmed_15x/trimming_summary_15x.txt"

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
