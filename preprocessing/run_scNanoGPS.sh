#! /bin/bash

# ================================
# sc/snNanoRNA-seq data processing
# Software: scNanoGPS v0.14
#
# Purpose:
# 1) Profile Nanopore read length distribution (data plotted in Extended Data Fig. 1a)
# 2) Identify cell barcodes and UMIs
# 3) Curate reads and generate consensus reads (used as input for IsoQuant)
# 4) Quantify gene expression in individual cells
# 5) Detect and annotate SNVs in individual cells
# 6) Generate summary QC reports
# ================================


# Path to scNanoGPS installation directory
P_DIR="/home/tools/scNanoGPS_v0.14/"

# Directory containing raw Nanopore FASTQ files
FASTQ="/data/nanopore/fastq_pass/"

# Human reference genome FASTA file, GRCh38
REF_GENOME="/home/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"

# Minimap2 index file for the GRCh38 reference genome
IND_GENOME="/home/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.mmi"

# Gene annotation file in GTF format
GENOME_ANNOTATION="/home/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf"

# BED file containing rRNA and hemoglobin regions to exclude
rRNA_HB_BED="/home/genome/gencode/rRNA_HB.gencode_v44.gene.bed"

# ANNOVAR installation directory and annotation database
ANNOVAR="/home/tools/annovar"
ANNOVAR_DB="/home/tools/annovar/hg38db/"

# Genome version used for ANNOVAR annotation
ANNOVAR_GV="hg38"

# ANNOVAR annotation databases used for SNV annotation
ANNOVAR_PROTOCOL="refGene,cytoBand,gnomad30_genome,avsnp150,dbnsfp42c,cosmic96_coding,cosmic96_noncoding"

# ANNOVAR operation types corresponding to each annotation database
ANNOVAR_OP="gx,r,f,f,f,f,f"

# Gene cross-reference file for additional gene annotation
ANNOVAR_XREF="/home/tools/annovar/hg38db/omim/gene_xref.txt"




# Number of cores used for parallel processing
ncores=20


# Step 1: Generate read length profile from raw Nanopore FASTQ files
python3 $P_DIR/other_utils/read_length_profiler.py -i $FASTQ &> run_read_length_profiler.log.txt &

# Step 2: Scan raw FASTQ files to identify putative cell barcodes and UMIs
python3 $P_DIR/scanner.py -i $FASTQ -t $ncores &> run_scanner.log.txt

# Step 3: Assign reads to individual cells/nuclei based on detected barcodes
python3 $P_DIR/assigner.py -t $ncores &> run_assigner.log.txt

# Step 4: Curate and error-correct reads, then generate consensus reads
python3 $P_DIR/curator.py -t $ncores --ref_genome $REF_GENOME --idx_genome $IND_GENOME --exc_bed $rRNA_HB_BED &> run_curator.log.txt

# Step 5: Generate gene-level expression matrix
python3 $P_DIR/reporter_expression.py --gtf $GENOME_ANNOTATION -t $ncores &> run_reporter_expression.log.txt

# Step 6: Detect SNVs from consensus reads and annotate variants using ANNOVAR
python3 $P_DIR/reporter_SNV.py --ref_genome $REF_GENOME -t $ncores --annovar $ANNOVAR --annovar_db $ANNOVAR_DB --annovar_gv $ANNOVAR_GV --annovar_protocol $ANNOVAR_PROTOCOL --annovar_operation $ANNOVAR_OP --annovar_xref $ANNOVAR_XREF &> run_reporter_SNV.log.txt

# Step 7: Convert ANNOVAR VCF output to a tab-delimited table
python3 $P_DIR/other_utils/parse_annovar_column.py -i scNanoGPS_res/annovar.hg38_multianno.vcf > scNanoGPS_res/annovar.hg38_multianno.tsv

# Step 8: Generate final summary report and QC metrics
python3 $P_DIR/reporter_summary.py --ref_genome $REF_GENOME --gtf $GENOME_ANNOTATION --qualimap_param "--java-mem-size=300G" &> run_reporter_summary.log.txt

