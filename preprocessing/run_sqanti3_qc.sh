#!/bin/bash

# ================================
# Isoform structural annotation and QC
# Software: SQANTI3 v5.2
#
# Purpose:
# Classify transcript structure (FSM,ISM,NIC,and NNC)
# QC for downstream analysis
# ================================


conda activate SQANTI3.env

# Input isoform annotation gtf file 
# Generated from IsoQuant output after removing novel genes
GTF="allsamples_extended_annotation_filtered_novelgene.gtf"

# Reference annotation and genome
REF_GTF="/home/genome/refdata-gex-GRCh38-2024-A/genes/genes.gtf"
REF_GENOME="/home/genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa"

#QC resources
CAGE="/home/tools/SQANTI3-5.2/data/ref_TSS_annotation/refTSS_v4.1_human_coordinate.hg38.bed"
POLYA="/home/tools/SQANTI3-5.2/data/polyA_motifs/mouse_and_human.polyA_motif.txt"

# Run SQANTI3 QC
python /home/tools/SQANTI3-5.2/sqanti3_qc.py $GTF \
$REF_GTF \
$REF_GENOME \
-o allsamples \
-d allsamples \
--cpus 10 \
--CAGE_peak $CAGE \
--polyA_motif_list $POLYA

