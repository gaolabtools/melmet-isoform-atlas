#! /bin/bash

# ================================
# Generate BAM files as input for IsoQuant
# Software:SAMtools v1.18

# Purpose:
# Add cell barcode tags,merge BAM files
# keep standard chromosomes, and remove reads from selected regions (NEAT1 and MALAT1)
  
# ================================

# Sample name, run for all samples
SAMPLE="sample_name"

# Directory containing BAM files (output files from scNanoGPS)
BAM_DIR="/home/isoquant/${SAMPLE}/bam"

# Output directory
OUTPUT_DIR="/home/isoquant/${SAMPLE}/bam"

# Step 1: Add CB tag to each BAM file
for bam_file in "$BAM_DIR"/*.bam; do
  filename=$(basename "$bam_file" .bam | cut -d '.' -f 1)
  samtools view -h "$bam_file" | awk -v tag="CB:Z:$filename" '{print $0"\t"tag}' | samtools view -Sb - > "$OUTPUT_DIR/$filename.tagged.bam"
done


# Step 2: Merge tagged BAM files
samtools merge "${OUTPUT_DIR}/${SAMPLE}.bam" "${OUTPUT_DIR}"/*.tagged.bam
samtools index "${OUTPUT_DIR}/${SAMPLE}.bam"


# Step 3: Keep only standard chromosomes
samtools view -b "${OUTPUT_DIR}/${SAMPLE}.bam" chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > "${OUTPUT_DIR}/${SAMPLE}_chr.bam"
samtools index "${OUTPUT_DIR}/${SAMPLE}_chr.bam"


# Step 4: Remove reads overlapping selected regions (NEAT1 and MALAT1)
samtools view "${OUTPUT_DIR}/${SAMPLE}_chr.bam" chr11:65422774-65445540 | cut -f1 > "${OUTPUT_DIR}/${SAMPLE}_reads.txt"
samtools view "${OUTPUT_DIR}/${SAMPLE}_chr.bam" chr11:65497640-65508073 | cut -f1 >> "${OUTPUT_DIR}/${SAMPLE}_reads.txt"

sort "${OUTPUT_DIR}/${SAMPLE}_reads.txt" | uniq > "${OUTPUT_DIR}/unique_${SAMPLE}_reads.txt"

samtools view -h "${OUTPUT_DIR}/${SAMPLE}_chr.bam" | grep -v -F -f "${OUTPUT_DIR}/unique_${SAMPLE}_reads.txt" | samtools view -b > "${OUTPUT_DIR}/${SAMPLE}_filtered.bam"
samtools index "${OUTPUT_DIR}/${SAMPLE}_filtered.bam"


echo "Processing for ${SAMPLE} complete."





