#!/bin/bash

# ================================
# Isoform-level quantification
# Software: IsoQuant v3.3.1

# Purpose:
# Detect novel transcript isoforms and quantify isoform expression in individual cells/spots 
# ================================


# Activate isoquant environment
conda activate /home/software/miniconda3/envs/isoquant

# Run isoquant pipeline
# Note: For bulk RNA-seq data, remove the --read_group option
isoquant.py \
--data_type nanopore \
--bam_list bam_file_list.txt \
--read_group tag:CB \
--genedb /home/genome/refdata-gex-GRCh38-2024-A/genes/genes.gtf \
--complete_genedb \
--reference /home/genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
--output allsamples \
--prefix allsamples \
--threads 20 \
--clean_start
