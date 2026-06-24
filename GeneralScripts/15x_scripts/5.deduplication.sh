# Remove PCR duplicates and optical duplicates that can introduce bias in downstream analyses.
# Reduce false positive variant calls and improve accuracy of quantitative analyses.
# EM-seq libraries can incorporate UMIs to more accurately identify and remove PCR duplicates

# Script adapted from Nuria SB scripts 
# Alvaro
# 03/12/2025

#Check if your files have UMIsInspect BAM file for UMI presence:
# samtools view your_file.bam | head -n 10
#@HISEQ:87:00000000_AATT
#In this example, AATT is the UMI.

src="/data2/gonzalo2/alvaro_data2/EMSeq_cov/15x/aligned_15x"
output="../../15x/deduplicated"

mkdir -p "$output"

#Run deduplication for all reads
for bam_file in "$src"/*.bam; do
    # Extract the base name of the BAM file (without the directory and extension)
    base_name=$(basename "$bam_file" .bam)

    # Run deduplicate_bismark with the --barcode option and specify the output directory
    /opt/software/Bismark-0.24.2/deduplicate_bismark --samtools_path /opt/software/samtools-1.11 --output_dir "$output" --paired --bam "$bam_file"
done
