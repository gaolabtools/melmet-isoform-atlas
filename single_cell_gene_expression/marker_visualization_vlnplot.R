# ================================
# Marker gene visualization for integrated snRNA-seq data
# Software: Seurat v4.4.0
#
# Purpose:
# Identify marker genes for major cell types across in-house long-read data and public short-read data
# Visualize selected marker genes 
# ================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(purrr)
library(grid)

source("StackedVlnPlot.R")


setwd("/home/BrainMET/1.integrated_all/sn/")
all<-readRDS("sn_rpca_dataset.RDS")

# ----------------
# Create combined cell type-dataset annotation
# ----------------

all$type<-paste0(all$celltype_main_final,"_",all$dataset)
all$type<-factor(all$type,levels=c("tumor cells_inhouse","tumor cells_jana","epithelial cells_inhouse","epithelial cells_jana","endothelial cells_inhouse","endothelial cells_jana","fibroblasts_inhouse","fibroblasts_jana","myocytes_inhouse","myeloid cells_inhouse",
                                   "myeloid cells_jana","T/NK cells_inhouse","T/NK cells_jana","plasma cells_inhouse","plasma cells_jana","B cells_inhouse","B cells_jana","mast cells_jana","oligodendrocytes_inhouse","oligodendrocytes_jana",
                                   "astrocytes_inhouse","astrocytes_jana","neurons_inhouse","neurons_jana"))


Idents(all)<-all$type

# ----------------
# Identify marker genes
# ----------------
markers<-FindAllMarkers(all,min.pct = 0.10,logfc.threshold = 1.5)

# ----------------
# Define marker genes used in Figure 1 and Extended Figure 1 
# ----------------
features_shared<-c("PAX3","TYR","DCT","PMEL","MUC5B","KRT19","DSP","KRT10","VWF","FLT1","PLVAP","CLDN5","COL1A1","COL3A1","COL6A3","DCN","IGFN1","MLIP","ITIH4","RYR1",
            "CD163","MRC1","CD86","C1QB","CD247","CD2","CD3G","IL7R","MZB1","IGHGP","IGHG1","IGLC1","MS4A1","BANK1","BLK","FCRL1","MS4A2","CTSG","IL18R1","HPGD",
            "TF","PCSK6","CNDP1","MOG","GPC5","SLC1A2","GFAP","ADGRV1","RIMS2","SYT1","OPCML","CDH18")


features_unique<-c("MUC5B","DUOX2","FCGBP","KRT19","DSP","KRT10","PKP1","KRT14","IGFN1","MLIP","ITIH4","RYR1","MS4A2","CTSG","IL18R1","HPGD")


# ----------------
# Plot stacked violin plots
# ----------------

color_inhouse<-"#74267D"
color_jana<-"#DEBA47"

group_colors <- c(
  color_inhouse, color_jana,
  color_inhouse, color_jana,
  color_inhouse, color_jana,
  color_inhouse, color_jana,
  color_inhouse,
  color_inhouse, color_jana,
  color_inhouse, color_jana,
  color_inhouse, color_jana,
  color_inhouse, color_jana,
  color_jana,
  color_inhouse, color_jana,
  color_inhouse, color_jana,
  color_inhouse, color_jana
)



pdf("marker_shared_vlnplot.pdf",width=8,height = 11)
StackedVlnPlot(all,features=features_shared,pt.size=0,cols=group_colors)
dev.off()

pdf("marker_unique_vlnplot.pdf",width=8,height = 5)
StackedVlnPlot(all,features=features_unique,pt.size=0,cols=group_colors)
dev.off()

