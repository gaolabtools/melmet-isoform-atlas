# ============================================================
# Adapted IsoformSwitchAnalyzeR workflow for single-cell isoform data
# ============================================================
#
# This script adapts IsoformSwitchAnalyzeR, originally designed for bulk RNA-seq with biological replicates, to single-cell isoform data.
#
# Pseudo-replicate expression matrices are generated only to construct a compatible IsoformSwitchAnalyzeR object and enable downstream
# transcript annotation, ORF integration, sequence extraction, and visualization.
#
# Final isoform fractions (IF), differential isoform usage (dIF), and statistical significance (FDR)  are replaced with values
# derived from the single-cell DIU analysis.
# As a result, all reported switching metrics reflect the single-cell-based analysis.

# ============================================================
#Example: MEL13 tumor cells ECM vs BM
# ============================================================

library(IsoformSwitchAnalyzeR)
library(rtracklayer)
library(DESeq2)
library(tidyverse)
library(readr)
library(dplyr)

# ----------------------------
# Load DIU result
# ----------------------------
######result of DIU analysis
DIU <- read.table("/home/BrainMET/5.DIU_ECM_vs_BM/individual samples/1.MEL13_tumor_isousage.txt",sep="\t",header=T) 
 

# ----------------------------
# Prepare isoform expression matrix
# ----------------------------
iso_expression = DIU %>%
  dplyr::select(isoform, ECM, BM) %>%
  dplyr::rename(isoform_id = "isoform")

set.seed(123)

######Pseudo-replicate expression matrices are generated only to construct a compatible IsoformSwitchAnalyzeR object
iso_expression <- iso_expression %>%
  dplyr::mutate(
    ECM_v1 = ECM,
    BM_v1  = BM,
    ECM_v2 = ECM + runif(n(), min = 0, max = 1),
    BM_v2  = BM  + runif(n(), min = 0, max = 1),
    ECM_v3 = ECM + runif(n(), min = 0, max = 1),
    BM_v3  = BM  + runif(n(), min = 0, max = 1)
  ) %>%
  dplyr::select(isoform_id, ECM_v1, BM_v1, ECM_v2, BM_v2, ECM_v3, BM_v3)


# ----------------------------
# Prepare design matrix
# ----------------------------
myDesign = tribble(
  ~sampleID, ~condition, ~donor,
  "ECM_v1", "ECM", "v1",
  "BM_v1", "BM", "v1",
  "ECM_v2", "ECM", "v2",
  "BM_v2", "BM", "v2",
  "ECM_v3", "ECM", "v3",
  "BM_v3", "BM", "v3",
) %>%
  dplyr::mutate(
    dplyr::across(c(condition, donor), as_factor)
  )


# ----------------------------
# Filter GTF annotation to DIU genes/isoforms
# ----------------------------
genes <- DIU$gene
######gtf file from isoquant output
gtf<-readGFF("/home/BrainMET/2.isoform_filtering/runindividual/MEL13/MEL13.transcript_models.gtf")
gtf<-gtf[gtf$gene_name %in% genes | gtf$transcript_id %in% iso_expression$isoform_id,]
flag <- na.omit(gtf[,c("gene_id","gene_name")])
gtf$gene_name <- flag[match(gtf$gene_id,flag$gene_id),]$gene_name
export(gtf, "filtered_genes.gtf")


# ----------------------------
# Create SwitchAnalyzeR object
# ----------------------------
######Create pre-filtered switchAnalyzeRlist
isoSwitchList <- importRdata(
  isoformCountMatrix   = iso_expression,
  designMatrix         = myDesign,
  isoformExonAnnoation = "filtered_genes.gtf",
  isoformNtFasta       = '/home/BrainMET/2.isoform_filtering/runindividual/MEL13/qc/MEL13_corrected.fasta', #squanti output
  addAnnotatedORFs     = FALSE,
  fixStringTieAnnotationProblem = FALSE # otherwise will mess up gene_ids
)
summary(isoSwitchList)


# ----------------------------
# Add ORF annotation from SQANTI3 output
# ----------------------------
gff_data<-import("/home/BrainMET/2.isoform_filtering/runindividual/MEL13/qc/MEL13_corrected.gtf.cds.gff")
gff_data$type <- as.character(gff_data$type)
export(gff_data, "/home/BrainMET/2.isoform_filtering/runindividual/MEL13/qc/MEL13_corrected.gtf.cds.gtf")

isoSwitchList <- addORFfromGTF(
  switchAnalyzeRlist     = isoSwitchList,
  pathToGTF              = '/home/BrainMET/2.isoform_filtering/runindividual/MEL13/qc/MEL13_corrected.gtf.cds.gtf'
)


# ----------------------------
# Pre-filter switch list
# ----------------------------
isoSwitchList <- preFilter(
  switchAnalyzeRlist         = isoSwitchList,
  geneExpressionCutoff       = 0,     # default
  isoformExpressionCutoff    = 0,     # default
  IFcutoff                   = 0,  # default
  removeSingleIsoformGenes   = TRUE,  # default
  reduceToSwitchingGenes     = FALSE, # default (we didn't run DEXSeq yet)
  keepIsoformInAllConditions = TRUE   # we only have 2 conditions so doesn't matter
)


# ----------------------------
# Run DEXSeq-based isoform switch testing
# ----------------------------
isoSwitchList_part1 <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist         = isoSwitchList,
  reduceToSwitchingGenes     = FALSE
)


# ----------------------------
# Replace switch statistics with original DIU results
# ----------------------------
matched_idx <- match(isoSwitchList_part1$isoformFeatures$isoform_id,DIU$isoform)

isoSwitchList_part1$isoformFeatures$IF1<-DIU[matched_idx,]$ECM_pct
isoSwitchList_part1$isoformFeatures$IF2<-DIU[matched_idx,]$BM_pct 
isoSwitchList_part1$isoformFeatures$iso_value_1<-DIU[matched_idx,]$ECM
isoSwitchList_part1$isoformFeatures$iso_value_2<-DIU[matched_idx,]$BM 
isoSwitchList_part1$isoformFeatures$dIF<-DIU[matched_idx,]$delta
isoSwitchList_part1$isoformFeatures$gene_switch_q_value<-DIU[matched_idx,]$FDR


# ----------------------------
# Extract amino acid sequences
# ----------------------------
isoSwitchList_part1$aaSequence = NULL
isoSwitchList_part1 <- extractSequence(
  switchAnalyzeRlist = isoSwitchList_part1,
  pathToOutput       = "/home/BrainMET/7.isoformSwitchAnalyzeR/MEL13/",
  extractNTseq       = TRUE,
  extractAAseq       = TRUE,
  removeShortAAseq   = FALSE,
  removeLongAAseq    = FALSE,
  onlySwitchingGenes = FALSE,
  dIFcutoff = 0.1,
)

saveRDS(isoSwitchList_part1, file = "isoSwitchList_part1.rds")
summary(isoSwitchList_part1)


# ----------------------------
# Run IsoformSwitchAnalyzeR part 2
# ----------------------------
isoSwitchList_part2 <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist        = isoSwitchList_part1, 
  n                         = 10, # number of PDF plots to generate
  removeNoncodinORFs        = FALSE,
  pathToCPC2resultFile      = "MEL13_cpc2.txt",
  pathToPFAMresultFile      = "MEL13_pfam.txt",
  pathToNetSurfP2resultFile  = "MEL13_NetsurfP.csv",
  pathToSignalPresultFile   = "MEL13_signalP.txt",
  pathToDeepLoc2resultFile  = "MEL13_DeepLoc.csv",
  pathToDeepTMHMMresultFile = "MEL13_TMRs.gff3",
  pathToOutput              = "output",
  outputPlots               = T,
  consequencesToAnalyze = c(
    'intron_retention',
    'coding_potential',
    'ORF_seq_similarity',
    'NMD_status',
    'domains_identified',
    'domain_isotype',
    'signal_peptide_identified'
  ),
)


# ----------------------------
# Analyze predicted consequences of isoform switches
# ----------------------------
isoSwitchList_part2 <- analyzeSwitchConsequences(isoSwitchList_part2, 
                                                 onlySigIsoforms = T, 
                                                 dIFcutoff = 0.1, 
                                                 consequencesToAnalyze = c(
                                                   'tss',
                                                   'tts',
                                                   'last_exon',
                                                   'isoform_length',
                                                   'exon_number',
                                                   'intron_structure',
                                                   'ORF_length', 
                                                   '5_utr_seq_similarity',
                                                   '5_utr_length', 
                                                   '3_utr_seq_similarity', 
                                                   '3_utr_length',
                                                   'coding_potential',
                                                   'ORF_seq_similarity',
                                                   'NMD_status',
                                                   'domains_identified',
                                                   'signal_peptide_identified',
                                                   'sub_cell_location'))

saveRDS(isoSwitchList_part2, file = "isoSwitchList_part2.rds")

pdf("switch_consequences.pdf",width=11,height=7)
switch_consequences <- extractConsequenceEnrichment(
  isoSwitchList_part2,dIFcutoff = 0.1,countGenes = F,
  returnResult = T # if TRUE returns a data.frame with the summary statistics
)
dev.off()

# ----------------------------
# Example transcript structure plot 
# ----------------------------
pdf("STX3.pdf",width=7,height=4)
switchPlot(
  isoSwitchList_part2,dIFcutoff = .1,logYaxis = TRUE,
  gene='STX3',rescaleTranscripts=FALSE,reverseMinus=FALSE)
dev.off()


  
  
  
