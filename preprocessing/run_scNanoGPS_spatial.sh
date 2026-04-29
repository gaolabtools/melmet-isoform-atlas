#! /bin/bash

# ================================
# Nanopore long-read spatial RNA-seq data processing 
# Software: scNanoGPS v2.0
#
# Purpose:
# 1) Curate reads and generate consensus reads (used as input for IsoQuant)
# 4) Quantify gene expression in individual spots
# ================================

P_DIR="/home/tools/scNanoGPS_v2"
FASTQ="/data/nanopore/sample/fastq_pass"
REF_GENOME="/home/genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
IND_GENOME="/home/genome/refdata-gex-GRCh38-2024-A/fasta/refdata-gex-GRCh38-2024-A.mmi"
GENOME_ANNOTATION="/home/genome/refdata-gex-GRCh38-2024-A/genes/genes.gtf"
WBC="/home/tools/scNanoGPS_v2/10x_barcodes/visium-v1.txt"
EXC_BED="/home/genome/gencode/rRNA_HB.gencode_v44.gene.bed"

ISOQUANT="/home/tools/IsoQuant/isoquant.py"
ANNOVAR="/home/tools/annovar"
ANNOVAR_DB="/home/tools/annovar/hg38db"
ANNOVAR_GV="hg38"
ANNOVAR_PROTOCOL="refGene,cytoBand,gnomad30_genome,avsnp150,dbnsfp42c,cosmic96_coding,cosmic96_noncoding"
ANNOVAR_OP="gx,r,f,f,f,f,f"
ANNOVAR_XREF="/home/tools/annovar/hg38db/omim/gene_xref.txt"
PT_SEQ="TTTTTTTTTTTT"

ncores=10
nanobc_no=6000

python3 $P_DIR/other_utils/read_length_profiler.py -i $FASTQ &> logs/run_read_length_profiler.log.txt &
python3 $P_DIR/scanner.py -t $ncores -i $FASTQ --pT $PT_SEQ &> logs/run_scanner.log.txt
python3 $P_DIR/assigner.py -t $ncores --whitelist $WBC --forced_no $nanobc_no &> logs/run_assigner.log.txt
python3 $P_DIR/spatial_bc_converter.py -t $ncores --whitelist $WBC &> logs/run_spatial_bc_converter.log.txt
python3 $P_DIR/curator.py -t $ncores --ref_genome $REF_GENOME --idx_genome $IND_GENOME --exc_bed $EXC_BED &> logs/run_curator.log.txt
python3 $P_DIR/reporter_expression.py -t $ncores --gtf $GENOME_ANNOTATION --min_gene_no 1 &> logs/run_reporter_expression.log.txt
python3 $P_DIR/reporter_isoform.py -t $ncores --ref_genome $REF_GENOME --gtf $GENOME_ANNOTATION --isoquant $ISOQUANT &> logs/run_reporter_isoform.log.txt
python3 $P_DIR/reporter_SNV.py -t $ncores --ref_genome $REF_GENOME --annovar $ANNOVAR --annovar_db $ANNOVAR_DB --annovar_gv $ANNOVAR_GV --annovar_protocol $ANNOVAR_PROTOCOL --annovar_operation $ANNOVAR_OP --annovar_xref $ANNOVAR_XREF &> logs/run_reporter_SNV.log.txt
python3 $P_DIR/other_utils/parse_annovar_column.py -i scNanoGPS_res/annovar.hg38_multianno.vcf > scNanoGPS_res/annovar.hg38_multianno.tsv
python3 $P_DIR/reporter_summary.py --ref_genome $REF_GENOME --gtf $GENOME_ANNOTATION --mrg_bam scNanoGPS_res/IsoQuant_res/merged.curated.minimap2.bam --qualimap_param "--java-mem-size=300G" &> logs/run_reporter_summary.log.txt
