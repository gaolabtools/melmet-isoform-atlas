############################################################
# Isoform diversity (Shannon index) across cell types in ECM and in BM 

#Purpose:
#    Quantify isoform diversity per gene per group using Shannon index
#    Visualize distributions of isoform diversity across different groups
#    Perform enrichment analysis using Fisher’s exact test

# Input requirements:
#   A) Seurat object with metadata:
#         - group; e.g., "tumor cells_BM" and "tumor cells_ECM" 
#
#   B) Isoform count matrix:
#        - R object "iso" (matrix/data.frame)
#        - rownames(iso) are "gene_isoform",i.e.,CCL4-202
#        - colnames(iso) are cell barcodes matching Seurat cell name
#
# Output:
#   - Density plots of per-gene Shannon index per group
#   - Tables of Shannon index per gene per group
#.  - Odds ratio enrichment plots

############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(purrr)
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

group_colors<-c("#3E7470","#3E7470","#7BC4C5" ,"#8C529D","#8C529D","#D0AFC4","#D0AFC4","#AB884B" ,"#AB884B" ,"#DEBF7E","#DEBF7E", "#3F4395","#8FA2D3","#dc7629","#f4d2aB","#B53E2B")


#############################
# Step 1: Aggregate isoform counts by group
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
# Step 2: Extract top-2 expressed isoforms
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
# Step 3: Calculate Shannon index per gene per cell type
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
# Step 4: Detection filter (gene expressed in >= min_pct_detected cells)
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
filtered_plot_shannon_index<-inner_join(plot_shannon_index, gene_keep, by = c("gene", "group"))

#dominance classification
#clasify genes into "with_dominant" and "without_dominant" based on shannon_index

filtered_plot_shannon_index <- filtered_plot_shannon_index %>%
  mutate(dominance = ifelse(shannon <= 0.4, "with_dominant", "without_dominant"))

#############################
# Step 5: Marker gene filter (keep only cell type marker genes)
#############################
degs <- read.table(deg, sep = "\t", header = TRUE)
# Expected columns :
#   - gene
#   - celltype
#   - p_val_adj
#   - avg_log2FC
markers <- degs %>%
  filter(p_val_adj < marker_padj, avg_log2FC > marker_logfc) %>%
  select(gene, celltype) %>%
  distinct()

markers<-markers %>%
  filter(!celltype %in% exclude_celltypes)

filtered_plot_shannon_index$celltype<-gsub("_.*", "", filtered_plot_shannon_index$group)
markers_filtered_plot_shannon_index<-inner_join(filtered_plot_shannon_index, markers, by = c("gene", "celltype"))


#############################
# Step 6: Plot isoform diversity distributions
#############################
p1 <- filtered_plot_shannon_index  %>%  mutate(group = factor(group, levels = group_order)) %>% 
  ggplot( aes(x=shannon, fill=group)) +
  geom_density( alpha=0.9, position = 'identity') +
  scale_fill_manual(values=group_colors) +
  labs(fill="",x ="shannon index")+
  facet_wrap(~group) +
  theme_test()
print(p1)

p2 <- markers_filtered_plot_shannon_index  %>% mutate(group = factor(group, levels = group_order)) %>% 
  ggplot( aes(x=shannon, fill=group)) +
  geom_density( alpha=0.9, position = 'identity') +
  scale_fill_manual(values=group_colors) +
  labs(fill="",x ="shannon index")+
  facet_wrap(~group) +
  theme_test()
print(p2)


#############################
# Save outputs 
#############################
write.table(
  filtered_plot_shannon_index,
  file = "shannon_index_group.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  markers_filtered_plot_shannon_index,
  file = "shannon_index_group_marker.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ggsave(
  filename = "shannon_index_density.pdf",
  plot = p1,
  width = 12,
  height = 8
)

ggsave(
  filename = "shannon_index_density_marker.pdf",
  plot = p2,
  width = 12,
  height = 8
)

#############################
# fisher test and odd ratios
#############################

dominance_percentage_allgenes <- filtered_plot_shannon_index %>%
  group_by(group) %>%
  summarize(
    with_dominant = sum(dominance == "with_dominant") / n() * 100,
    without_dominant = sum(dominance == "without_dominant") / n() * 100
  )


dominance_percentage_markers <- markers_filtered_plot_shannon_index %>%
  group_by(group) %>%
  summarize(
    with_dominant = sum(dominance == "with_dominant") / n() * 100,
    without_dominant = sum(dominance == "without_dominant") / n() * 100
  )


results <- lapply(1:nrow(dominance_percentage_allgenes), function(i) {
  control_row <- dominance_percentage_allgenes[i, ]
  case_row <- dominance_percentage_markers[i, ]
  
  # Construct the 2x2 contingency table
  contingency_table <- matrix(c(control_row$with_dominant, control_row$without_dominant,
                                case_row$with_dominant, case_row$without_dominant), 
                              nrow = 2, byrow = TRUE)
  
  # Perform Fisher's exact test
  fisher_test <- fisher.test(contingency_table)
  
  odds_ratio_results <- oddsratio(contingency_table, method = "wald")
  odds_ratio <- odds_ratio_results$measure[2, "estimate"]
  lower_ci <- odds_ratio_results$measure[2, "lower"]
  upper_ci <- odds_ratio_results$measure[2, "upper"]
  
  # Return the results
  list(group = control_row$group, p.value = fisher_test$p.value, odds_ratio = odds_ratio,
       lower_ci = lower_ci, upper_ci = upper_ci)
})

# Print the results
results_df <- do.call(rbind, lapply(results, as.data.frame))
results_df$color <- ifelse(results_df$p.value < 0.05, "#F8B496", "grey")
results_df$group<-factor(results_df$group, levels = rev(group_order))


p3 <- ggplot(results_df, aes(x = group, y = odds_ratio, color = color)) +
  geom_segment(aes(x = group, xend = group, y = lower_ci, yend = upper_ci, color = color), linewidth = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  geom_point(size = 1.5) +
  scale_color_identity() +
  coord_flip() +
  labs(title = "Odds Ratios with 95% CI", x = "Group", y = "Odds Ratio") +
  theme_test()

ggsave(
  filename = "diversity_odds_ratio.pdf",
  plot = p3,
  width = 4,
  height = 3
)










