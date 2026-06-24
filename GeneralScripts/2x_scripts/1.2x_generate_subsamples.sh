#!/bin/bash

# Alvaro Gutierrez Rodriguez
# 07/02/2026

# 2x downsampling for paired-end FASTQs

# Directories
INPUT_DIR="/data2/gonzalo2/alvaro_data2/EMSeq_cov/15x/raw_15x"
OUTPUT_DIR="/data2/gonzalo2/alvaro_data2/EMSeq_cov/2x/raw_fastq_2x_subsampled"

# Fraction for 2x from 15x (exact calculation)
FRACTION=$(echo "scale=6; 2/15" | bc)
SEED=42  # reproducible

mkdir -p "$OUTPUT_DIR"

# Use a single Docker container for all samples (more efficient)
docker run --rm \
  -u "$UID" \
  -v "$INPUT_DIR":/input \
  -v "$OUTPUT_DIR":/output \
  staphb/seqtk \
  bash -c '
    echo "Starting 2x downsampling from 15x..."
    
    for R1 in /input/*_R1.fastq.gz; do
      BASENAME=$(basename "$R1" "_R1.fastq.gz")
         
      R2="/input/${BASENAME}_R2.fastq.gz"
      
      # Validate that R2 file exists
      if [[ ! -f "$R2" ]]; then
        echo "ERROR: R2 file not found for $BASENAME. Skipping..." >&2
        continue
      fi
      
      echo "Downsampling $BASENAME to 2x coverage"
      
      zcat "$R1" | seqtk sample -s'"$SEED"' - '"$FRACTION"' | gzip > "/output/${BASENAME}_R1.fastq.gz"
      zcat "$R2" | seqtk sample -s'"$SEED"' - '"$FRACTION"' | gzip > "/output/${BASENAME}_R2.fastq.gz"
    done
    
    echo "2x downsampling complete in container!"
  '

echo "2x downsampling complete! Output in $OUTPUT_DIR"
