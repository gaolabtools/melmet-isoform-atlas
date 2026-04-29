#! /bin/bash

# ================================
# Bulk Nanopore RNA-seq data processing
# Software: NanoQC v0.9.4, NanoPlot v1.41.6, NanoFilt v2.8.0, minimap2 v2.26
#
# Purpose:
# Generate BAM files as input for IsoQuant
# ================================

# Sample name
SAMPLE="MEL01ECM"

# Raw FASTQ directory
FASTQ_DIR="/data/nanopore/MBM_bulk/fastq_pass/${SAMPLE}"

# Reference genome and annotation
REF_GENOME="/home/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
GTF="/home/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf"

# Merge FASTQ files
cat "$FASTQ_DIR"/*.fastq.gz > "pass_${SAMPLE}.fastq.gz"

# Quality control 
conda activate nanopack
nanoQC "pass_${SAMPLE}.fastq.gz" -o nanoQC

conda activate nanoplot_env       
NanoPlot -t 10 --fastq "pass_${SAMPLE}.fastq.gz" --maxlength 40000 --plots hex dot  -o nanoplot

# Filtering and Trimming
conda activate nanofilt
gunzip -c "pass_${SAMPLE}.fastq.gz" | NanoFilt -l 100 --headcrop 10 | gzip > "nanofilt_${SAMPLE}.fastq.gz"

#Align reads to the reference genome using minimap2
conda activate minimap2 
minimap2 -ax splice -t 10 "$REF_GENOME" "nanofilt_${SAMPLE}.fastq.gz" | samtools sort -o "sorted_${SAMPLE}.bam"
samtools index "sorted_${SAMPLE}.bam"

