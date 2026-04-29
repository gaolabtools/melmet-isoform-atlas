#Differential Isoform Usage (DIU) script
#Purpose:
  #   Compare isoform usage between two conditions (e.g., BM vs ECM, between different cell types, between different cell states)
#Input requirements:
  #   1) Seurat object with metadata columns:
  #        - group/condition; e.g., "BM" and "ECM" or cell type/cell state labels 
  #
  #   2) Isoform count matrix:
  #        - R object "iso" (matrix/data.frame)
  #        - rownames(iso) are "gene_isoform",i.e.,CCL4-202
  #        - colnames(iso) are cell barcodes matching Seurat cell names
  #
# Output:
  #   A table of isoform usage, delta usage, chi-square p-values, BH-FDR per gene.
  #

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(Matrix)
})

#############################
# parameters
#############################

# Paths 
seurat_rds <- "/path/to/inhouse.RDS"
iso_rds    <- "/path/to/isoform_counts.RDS"

# Subset parameters 
patient    <- "MEL13"
celltype   <- "tumor cells" 

# Group labels 
group_col  <- "site"
group1     <- "ECM"
group2     <- "BM"

# Filtering thresholds
min_pct_detected <- 5    # minimum % of cells expressing gene in EACH group
min_gene_umi      <- 10  # minimum gene-level UMI in EACH group (aggregated)
fdr_cutoff        <- 0.05
delta_cutoff      <- 0.10


# Output
out_tsv <- "MEL13_tumor cells_isoform_usage.tsv"


#############################
# Load data
#############################

seu <- readRDS(seurat_rds)
iso <- readRDS(iso_rds)

# Subset Seurat object by cell type
seu_sub <- subset(seu, patient == patient_id & celltype_main_final == celltype)

# Ensure both groups present
grp_tab <- table(seu_sub@meta.data[[group_col]])
if (!(group1 %in% names(grp_tab) && group2 %in% names(grp_tab))) {
  stop("Selected subset does not contain both groups. Check patient/celltype/group labels.")
}

message("Cells in subset: ", ncol(seu_sub))
message("Group counts:\n", paste(capture.output(print(grp_tab)), collapse = "\n"))

#############################
# Align isoform matrix columns to selected cells
#############################
cells_use <- rownames(seu_sub@meta.data)

# Build annotation dataframe
anno <- seu_sub@meta.data[cells_use, , drop = FALSE] %>%
  rownames_to_column("cellbarcode") %>%
  transmute(cellbarcode, group = .data[[group_col]])

# Subset isoform matrix to selected cells
iso_sub <- iso[, anno$cellbarcode, drop = FALSE]

## Drop isoforms with all zeros (faster downstream)
keep_rows <- rowSums(iso_sub) > 0
iso_sub <- iso_sub[keep_rows, , drop = FALSE]

message("Isoforms after nonzero filter: ", nrow(iso_sub))


#############################
# Step 1: Filter genes by detection fraction in BOTH groups
#############################
# Convert to long gene counts per cell by summing isoforms for each gene in each cell
gene_by_cell <- iso_sub %>% as.data.frame() %>% tibble::rownames_to_column(var = "row_name") %>%
  tidyr::separate(row_name, into=c("gene","isoform"), sep="_") %>%
  dplyr::group_by(gene) %>% dplyr::summarise(across(where(is.numeric), sum))  %>%
  pivot_longer(cols = -gene, names_to = "cellbarcode", values_to = "counts") %>%
  left_join(anno, by = "cellbarcode")   

# only test genes that are detected in a minimum fraction of min.pct cells in both two groups
pct_cell <- gene_by_cell %>%
  group_by(gene, group) %>%
  summarise(percent_nonzero = mean(counts > 0) * 100, .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = percent_nonzero)

genes_keep <- pct_cell %>%
  filter(.data[[group1]] >= min_pct_detected, .data[[group2]] >= min_pct_detected) %>%
  pull(gene)

message("Genes passing detection filter: ", length(genes_keep))


#############################
# Step 2: Aggregate isoform counts per group (pseudo-bulk per group)
#############################
groups <- c(group1, group2)

counts_list <- list()

for (i in groups) {
  barcodes_of_type <- filter(anno, group == i)$cellbarcode
  isoforms_of_type <- iso_sub %>% dplyr::select(barcodes_of_type)
  counts<-rowSums(isoforms_of_type)
  counts_list[[i]] <- counts
}

counts_cs <- do.call(cbind, counts_list) %>% as.data.frame() %>% tibble::rownames_to_column(var = "row_name") %>%
  tidyr::separate(row_name, into=c("gene","isoform"), sep="_") %>% dplyr::filter(gene %in% genes_keep)


#############################
# Step 3: Filter genes by minimum total UMI per group
#############################

geneUMI<-counts_cs %>%
  group_by(gene) %>%
  dplyr::summarise(across(where(is.numeric), sum)) %>%
  filter(across(all_of(groups), ~ . >= min_gene_umi))

counts_cs<- counts_cs %>% dplyr::filter(gene %in% geneUMI$gene)

geneUMI <- counts_cs %>%
  group_by(gene) %>%
  summarise(across(all_of(groups), sum), .groups = "drop") %>%
  filter(across(all_of(groups), ~ . >= min_gene_umi))

counts_cs <- counts_cs %>%
  filter(gene %in% geneUMI$gene)

message("Genes passing min UMI filter: ", length(unique(counts_cs$gene)))

#############################
# Step 4: Keep genes with multiple isoform
#############################

genes_with_multiple_isoforms <-  counts_cs %>% dplyr::select(gene, isoform) %>% group_by(gene) %>% 
  summarise(numIso=n_distinct(isoform ))  %>% filter(numIso>1) %>% dplyr::select(gene) %>% pull()

message("Genes with multiple isoforms: ", length(genes_with_multiple_isoforms))

counts_cs<- counts_cs %>% filter(gene %in% genes_with_multiple_isoforms) 


#############################
# Step 5: Compute isoform usage and delta
#############################
long_counts <- counts_cs  %>% pivot_longer(-c("gene", "isoform"),names_to = "group", values_to = "counts")
gene_counts <- long_counts  %>% group_by(gene, group ) %>% summarise(gene_counts=sum(counts))
iso_usage  <- long_counts  %>% left_join(gene_counts, by=c("gene", "group")) %>% mutate(iso_usage=counts/gene_counts) %>% 
  dplyr::select(gene,isoform, group , iso_usage) %>% pivot_wider(names_from = group, values_from = iso_usage)

colnames(iso_usage)[match(group1, colnames(iso_usage))] <- paste0(group1, "_pct")
colnames(iso_usage)[match(group2, colnames(iso_usage))] <- paste0(group2, "_pct")

# Merge back counts + usage
res <- counts_cs %>%
  left_join(iso_usage, by = c("gene", "isoform"))

# Delta
res <- res %>%
  mutate(delta = .data[[paste0(group2, "_pct")]] - .data[[paste0(group1, "_pct")]]) %>%
  arrange(desc(delta))  %>% as.data.frame()


#############################
# Step 6: Gene-level chi-square test on isoform composition
#############################
# For each gene, build an isoform-by-group contingency table:
#   rows = isoforms, cols = groups, entries = aggregated counts.
# Then run chisq.test to get p-value, BH adjust across genes.

chisq_by_gene <- res %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    pvalue = chisq.test(
      matrix(c(.data[[group1]], .data[[group2]]), ncol = 2, byrow = FALSE)
    )$p.value
  ) %>%
  dplyr::mutate(FDR = p.adjust(pvalue, method = "BH"))

res <- res %>%
  left_join(chisq_by_gene, by = "gene")

head(res)


#############################
# Step 7: Summarize significant DIU
#############################

message("Significant genes (FDR <= ", fdr_cutoff, ", |delta| >= ", delta_cutoff, "): ",
        n_distinct(sig$gene))
message("Significant isoforms: ", n_distinct(sig$isoform))

#############################
# Write output
#############################

write.table(res, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote: ", out_tsv)
