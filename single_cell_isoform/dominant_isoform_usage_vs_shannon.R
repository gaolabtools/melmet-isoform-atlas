############################################################
# Relationship between dominant isoform usage and isoform diversity (Shannon index)

############################################################


suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(purrr)
  library(ggtrendline)
})


#############################
#parameters
#############################
# Paths
seurat_rds <- "/home/BrainMET/1.integrated_all/inhouse/inhouse.RDS"
iso_rds    <- "/home/BrainMET/2.isoform_filtering/allsamples/isoform_counts.RDS"

# Cell types excluded due to fewer than 100 cells
exclude_celltypes <- c("neurons", "myocytes")

#keep genes detected  in at least this % of cells in each group
min_pct_detected <- 5

# Marker for major cell types:
deg <- "DEGs_major_celltype.txt"  
marker_padj <- 0.05
marker_logfc <- 1.5

#############################
# Load data
#############################
seu <- readRDS(seurat_rds)
iso <- readRDS(iso_rds)

#############################
# Build annotation table
#############################
anno <- seu@meta.data %>%
  rownames_to_column("cellbarcode") %>%
  transmute(
    cellbarcode,
    group = paste0(celltype_main_final, "_", site))

collapse_map <- c(
  "B cells_BM"       = "B cells",
  "B cells_ECM"  = "B cells",
  "plasma cells_BM"           = "plasma cells",
  "plasma cells_ECM"      = "plasma cells",
  "astrocytes_BM"               = "astrocytes",
  "gallbladder epithelial_ECM" = "epithelial",
  "neurons_BM"                  = "neurons",
  "oligodendrocytes_BM"         = "oligodendrocytes",
  "myocytes_ECM"   = "myocytes"
)  

anno$group <- ifelse(anno$group %in% names(collapse_map),
                     collapse_map[anno$group],
                     anno$group)

anno <- anno %>%
  filter(!group %in% exclude_celltypes)

# align isoform matrix
iso  <- iso[, anno$cellbarcode, drop = FALSE]


group_order<-c("tumor cells_ECM","tumor cells_BM","epithelial","endothelial cells_ECM","endothelial cells_BM",
               "fibroblasts_ECM","fibroblasts_BM","myeloid cells_ECM","myeloid cells_BM",
               "T/NK cells_ECM","T/NK cells_BM","plasma cells","B cells",
               "oligodendrocytes","astrocytes")


#############################
# Aggregate isoform counts by group
#############################
# Result: counts_cs has columns:
#   gene, isoform, <group1>, <group2>, ...

counts_list <- list()

for (i in as.character(unique(anno$group))) {
  barcodes_of_type <- filter(anno, group == i)$cellbarcode
  isoforms_of_type <- iso %>% dplyr::select(barcodes_of_type)
  counts<-rowSums(isoforms_of_type)
  counts_list[[i]] <- counts
}

counts_cs <- do.call(cbind, counts_list) %>% as.data.frame() %>% tibble::rownames_to_column(var = "row_name") %>%
  tidyr::separate(row_name, into=c("gene","isoform"), sep="_")

#############################
# Extract top-2 expressed isoforms
#############################

#Extract expression of dominant isoform  per gene per group
dominant_tx = counts_cs %>% group_by(gene) %>% summarise(across(where(is.numeric), max)) 
#Extract expression of second dominant isoform per gene per group
second_dominant_tx <- counts_cs %>%
  group_by(gene) %>%
  summarize(across(where(is.numeric), ~nth(sort(., decreasing = TRUE), 2)))
second_dominant_tx [is.na(second_dominant_tx)] = 0

two_dominant_tx<-rbind(dominant_tx,second_dominant_tx)

#############################
# Calculate Shannon index 
#############################

#shannon index
calc_shannon <- function(x) {
  p <- x / sum(x)
  p <- p[p > 0]  
  shannon_index <- -sum(p * log(p))
  return(shannon_index)
}

shannon_index <-two_dominant_tx %>%
  pivot_longer(-gene, names_to = "group", values_to = "counts") %>%
  group_by(gene, group) %>%
  summarise(Shannon_Index = calc_shannon(counts), .groups = 'drop') %>%
  pivot_wider(names_from = group, values_from = Shannon_Index)


#############################
# most dominant isoform uasge 
#############################
#calculate isoform usage
long_counts <- counts_cs  %>% pivot_longer(-c("gene", "isoform"),names_to = "celltype", values_to = "counts")
gene_counts = long_counts  %>% group_by(gene, celltype ) %>% summarise(gene_counts=sum(counts))
long_counts  <- long_counts  %>% left_join(gene_counts, by=c("gene", "celltype")) %>% mutate(iso_usage=counts/gene_counts)
long_counts$iso_usage[is.na(long_counts$iso_usage)] = 0
iso_usage = long_counts %>% dplyr::select(gene,isoform, celltype , iso_usage) %>% pivot_wider(names_from = celltype, values_from = iso_usage)

#most dominant iso usage
dominant_tx = iso_usage %>% group_by(gene) %>% summarise(across(where(is.numeric), max)) 
dominant_tx = dominant_tx[,c("gene",celltype_order)]



#############################
# Detection filter (gene expressed in >= min_pct_detected cells)
#############################
pct_cell<-iso %>% tibble::rownames_to_column(var = "row_name") %>%
  tidyr::separate(row_name, into=c("gene","isoform"), sep="_") %>%
  dplyr::group_by(gene) %>% dplyr::summarise(across(where(is.numeric), sum))  %>%
  pivot_longer(cols = -gene, names_to = "cellbarcode", values_to = "counts") %>%
  left_join(anno, by = "cellbarcode")   

pct_cell <- pct_cell %>%
  group_by(gene, group) %>%
  summarise(percent_nonzero = mean(counts > 0) * 100, .groups = 'drop')

pct_cell<-pct_cell %>% filter(!group %in% exclude_celltypes)

gene_keep<-pct_cell %>% dplyr::filter(percent_nonzero >= min_pct_detected)

plot_shannon_index<-shannon_index %>% pivot_longer(-gene, names_to = "group", values_to = "shannon") %>% 
  dplyr::filter(!is.na(shannon)) 
plot_dominant_tx<-dominant_tx %>% pivot_longer(-gene, names_to = "group", values_to = "usage") %>% 
  dplyr::filter(usage != 0) 


dominant_tx_long <- dominant_tx %>%
  pivot_longer(-gene, names_to = "group", values_to = "usage")

plot_shannon_dominant <- plot_shannon_index %>%
  left_join(dominant_tx_long, by = c("gene", "group"))  


filtered_plot_shannon_dominant<-inner_join(plot_shannon_dominant, gene_keep, by = c("gene", "group"))


#############################
# plotting
#############################

pdf("correlation_dominant_usage_shannon.pdf",width=6,height = 8)
ggtrendline(filtered_plot_shannon_dominant$usage,filtered_plot_shannon_dominant$shannon,model="power3P",linecolor = "red",linetype = 1,linewidth = 1,
            CI.level=0.99,CI.fill = "#5ab4ac",CI.alpha = 0.2,CI.color = NA,CI.lty = 2,CI.lwd = 1.5,eq.x=0.7,rrp.x=0.7)+
  geom_point(aes(filtered_plot_shannon_dominant$usage,filtered_plot_shannon_dominant$shannon),size=0.1,color="grey",shape=1)+
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(title = "Scatter Plot of Dominant Usage vs Shannon Index",
       x = "Most dominant isoform usage",
       y = "Shannon index") +
  theme_bw()
dev.off()











