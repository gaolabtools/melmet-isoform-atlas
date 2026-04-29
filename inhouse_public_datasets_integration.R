# ================================
# Integration of in-house and public snRNA-seq datasets
# Software: Seurat v4.4.0
#
# Purpose:
# Integrate datasets using reciprocal PCA (RPCA)
# Perform dimensionality reduction, clustering, and visualization
# ================================


library(Seurat)
library(dplyr)

setwd("/home/BrainMET/1.integrated_all/sn/")
inhouse<-readRDS("../inhouse/inhouse.RDS") 
jana<-readRDS("../public_data/cell/filtered.RDS")
inhouse$dataset<-"inhouse"
jana$dataset<-"jana"


# ----------------
# Filter public single-nucleus dataset
# ----------------

jana <- jana %>%
  subset(
    cell_type_main != "Low-quality cells" &
      !cell_type_int %in% c("Doublets", "Contamination", "Low-quality cells") &
      sequencing == "Single nuclei"
  )


# ----------------
# Retain shared expressed genes
# ----------------
jana_genes<-row.names(jana)[rowSums(jana@assays$RNA@counts)!=0]
inhouse_genes<-row.names(inhouse)[rowSums(inhouse@assays$RNA@counts)!=0]
universe <- intersect(jana_genes, inhouse_genes)
inhouse<-inhouse[universe,]
jana<-jana[universe,]


# ----------------
# Merge datasets
# ----------------
all<-merge(jana,c(inhouse))


# ----------------
# RPCA-based integration
# ----------------

all.list <- SplitObject(all, split.by = "dataset")

all.list <- lapply(X = all.list, FUN = function(x) {
  x <- NormalizeData(x) 
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = all.list)

all.list <- lapply(X = all.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

all.anchors <- FindIntegrationAnchors(object.list = all.list, anchor.features = features, reduction = "rpca")
all.combined <- IntegrateData(anchorset = all.anchors)

DefaultAssay(all.combined) <- "integrated"


# ----------------
# Dimensionality reduction and clustering
# ----------------

all.combined <- ScaleData(all.combined, verbose = FALSE)
all.combined <- RunPCA(all.combined, npcs = 50, verbose = FALSE)
all.combined <- RunUMAP(all.combined, reduction = "pca", dims = 1:50)
all.combined <- FindNeighbors(all.combined, reduction = "pca", dims = 1:50)
all.combined <- FindClusters(all.combined, resolution = 0.5)

DefaultAssay(all.combined)<-"RNA"
saveRDS(all.combined,"sn_rpca_dim50.RDS")


# ----------------
# Harmonize cell-type annotations
# ----------------
all<-readRDS("sn_rpca_dim50.RDS")

all$celltype_main_final <- case_when(
  all$celltype_main_final == "gallbladder epithelial" ~ "epithelial cells",
  all$dataset == "jana" & all$cell_type_int == "Mast cells" ~ "mast cells",
  all$dataset == "jana" & all$cell_type_int == "Plasma cells" ~ "plasma cells",
  all$dataset == "jana" & all$cell_type_int == "B cells" ~ "B cells",
  all$dataset == "jana" & all$cell_type_fine == "Astrocytes" ~ "astrocytes",
  all$dataset == "jana" & all$cell_type_fine == "Oligodendrocytes" ~ "oligodendrocytes",
  all$dataset == "jana" & all$cell_type_fine == "Neurons" ~ "neurons",
  all$dataset == "jana" & all$cell_type_main == "T/NK cells" ~ "T/NK cells",
  all$dataset == "jana" & all$cell_type_main == "Myeloid cells" ~ "myeloid cells",
  all$dataset == "jana" & all$cell_type_main == "Endothelial cells" ~ "endothelial cells",
  all$dataset == "jana" & all$cell_type_main == "Stromal cells" ~ "fibroblasts",
  all$dataset == "jana" & all$cell_type_main == "Tumor cells" ~ "tumor cells",
  all$dataset == "jana" & all$cell_type_main == "Epithelial cells" ~ "epithelial cells",
  TRUE ~ as.character(all$celltype_main_final)
)

all$celltype_main_final<-factor(all$celltype_main_final,levels=c("tumor cells","epithelial cells","endothelial cells","fibroblasts","myocytes","myeloid cells","T/NK cells","plasma cells","B cells","mast cells","oligodendrocytes","astrocytes","neurons"))


# ----------------
# Plot UMAP by cell type
# ----------------


celltype_colors<-c("#3E7470","#7BC4C5" ,"#8C529D","#D0AFC4","#AB3282","#AB884B" ,"#DEBF7E", "#3F4395","#8FA2D3","#2171b5","#dc7629","#f4d2aB","#B53E2b")


pdf("all_celltype.pdf",width=10,height=10)
DimPlot(all, reduction = "umap", label = TRUE,group.by="celltype_main_final",raster=FALSE,cols = celltype_colors)
DimPlot(all, reduction = "umap", label = TRUE,group.by="celltype_main_final",raster=FALSE,cols = celltype_colors)+theme_test()+ theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),legend.position = "none")
DimPlot(all, reduction = "umap", label = FALSE,group.by="celltype_main_final",raster=FALSE,cols = celltype_colors)+theme_test()+ theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),legend.position = "none")
dev.off()


# ----------------
# Plot UMAP by dataset
# ----------------

all$dataset<-factor(all$dataset,levels=c("jana","inhouse"))


pdf("output_dataset.pdf",width=10,height=10)
DimPlot(all, reduction = "umap", label = TRUE,group.by="dataset",raster=FALSE,cols=c("#DEBA47","#74267D"))
DimPlot(all, reduction = "umap", label = TRUE,group.by="dataset",raster=FALSE,cols=c("#DEBA47","#74267D"))+theme_test()+ theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),legend.position = "none")
DimPlot(all, reduction = "umap", label = FALSE,group.by="dataset",raster=FALSE,cols=c("#DEBA47","#74267D"))+theme_test()+ theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),legend.position = "none")
dev.off()




