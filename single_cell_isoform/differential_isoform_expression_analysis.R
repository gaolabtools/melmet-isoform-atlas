# ============================================================
# Differential isoform expression analysis
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)


# ----------------------------
# Load Seurat object
# ----------------------------
inhouse<-readRDS("/home/BrainMET/1.integrated_all/inhouse/inhouse.RDS")

celltype_colors <-c("#3E7470","#7BC4C5" ,"#8C529D","#D0AFC4","#AB3282","#AB884B" ,"#DEBF7E", "#3F4395","#8FA2D3","#DC762D","#F4D2A8","#B53E2B")



# ----------------------------
# Load isoform count matrix
# ----------------------------
iso<-readRDS("/home/BrainMET/2.isoform_filtering/allsamples/isoform_counts.RDS")

# Add isoform count matrix as a new assay
inhouse[["iso"]] <- CreateAssayObject(counts = iso)
DefaultAssay(inhouse)<-"iso"
inhouse <- NormalizeData(inhouse,normalization.method = "LogNormalize",scale.factor = 10000)
inhouse <- ScaleData(object = inhouse)
inhouse <- FindVariableFeatures(inhouse, selection.method = "vst", nfeatures = 10000)
inhouse <- RunPCA(inhouse,features = VariableFeatures(inhouse)) 
pcSelect=50
inhouse <- FindNeighbors(inhouse,dims = 1:pcSelect)   
inhouse <- FindClusters(inhouse, resolution = 0.5)
inhouse <- RunUMAP(inhouse,dims = 1:pcSelect)  

# ----------------------------
# Identify differentially expressed isoforms
# ----------------------------
Idents(inhouse)<-inhouse$celltype_main_final 

markers<-FindAllMarkers(inhouse,min.pct=0.1,logfc.threshold = 0.25,only.pos = TRUE)

top20<-markers %>% group_by(cluster) %>% top_n(20,avg_log2FC)

# Scale selected marker isoforms
inhouse<-ScaleData(inhouse,features = top20$gene)


# ----------------------------
# order cells and generate heatmap
# ----------------------------
inhouse$group <-paste0(inhouse$site,"_",inhouse$celltype_main_final)
group_levels <- c("ECM_tumor cells","BM_tumor cells","ECM_gallbladder epithelial","ECM_endothelial cells","BM_endothelial cells","ECM_fibroblasts","BM_fibroblasts","ECM_myocytes","ECM_myeloid cells",
                    "BM_myeloid cells","ECM_T/NK cells","BM_T/NK cells","ECM_plasma cells","BM_plasma cells","ECM_B cells","BM_B cells","BM_oligodendrocytes","BM_astrocytes","BM_neurons")
inhouse$group<-factor(inhouse$group,levels =group_levels)
                        
                        
ordered_cells <- unlist(
  lapply(group_levels, function(x) {
    Cells(subset(inhouse, subset = group == x))
  })
)

pdf("heatmap_top20_ordered.pdf",width=70,height = 20)
DoHeatmap(inhouse,features = top20$gene,group.colors = celltype_colors,cells = ordered_cells,draw.lines = TRUE,label = F)+scale_fill_viridis()
dev.off()


