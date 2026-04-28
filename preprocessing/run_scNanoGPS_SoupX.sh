#! /bin/bash

# ================================
# snNanoRNA-seq data preprocessing for SoupX input generation
# Software: scNanoGPS v0.14
#
# Purpose:
# This script generates gene expression profiles required as input for SoupX,
# a tool used to estimate and remove ambient RNA contamination in single cell data.
#
# ================================

# Path to scNanoGPS installation
P_DIR="/home/tools/scNanoGPS_v0.14/"

# Reference genome (GRCh38)
REF_GENOME="/home/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"

# Minimap2 genome index
IND_GENOME="/home/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.mmi"

# Gene annotation (GTF)
GENOME_ANNOTATION="/home/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf"

# Number of cores
ncores=20

#forces inclusion of up to 20,000 barcodes, ensuring sufficient droplets are retained for ambient RNA estimation
python3 $P_DIR/assigner.py -t $ncores --forced_no 20000 &> run_assigner.log.txt
python3 $P_DIR/curator.py -t $ncores --ref_genome $REF_GENOME --idx_genome $IND_GENOME --skip_curation 1 &> run_curator.log.txt
python3 $P_DIR/reporter_expression.py --gtf $GENOME_ANNOTATION -t $ncores --min_gene_no 1 &> run_reporter_expression.log.txt
python3 $P_DIR/reporter_summary.py --ref_genome $REF_GENOME --gtf $GENOME_ANNOTATION --qualimap_param "--java-mem-size=300G" &> run_reporter_summary.log.txt

