
# ============================================================
# Summarize isoform switch functional consequences across melanoma samples
# ============================================================

library(IsoformSwitchAnalyzeR)
library(tidyverse)
library(ggplot2)


# ----------------------------
# Define paths and samples
# ----------------------------

setwd("/home/BrainMET/7.isoformSwitchAnalyzeR/")

sample_ids <- c(
  "MEL01", "MEL02", "MEL03", "MEL04",
  "MEL07", "MEL11", "MEL12", "MEL13"
)

# ----------------------------
# Load IsoformSwitchAnalyzeR objects
# ----------------------------

iso_switch_lists <- lapply(sample_ids, function(sample_id) {
  readRDS(file.path(sample_id, "isoSwitchList_part2.rds"))
})

names(iso_switch_lists) <- sample_ids


# ----------------------------
# Extract isoform switch functional consequences
# ----------------------------

extract_switch_consequences <- function(iso_switch_list, sample_id) {
  iso_switch_list$switchConsequence %>%
    dplyr::mutate(
      sample = sample_id,
      event_id = paste(
        gene_name,
        isoformUpregulated,
        isoformDownregulated,
        sep = "_"
      )
    ) %>%
    dplyr::select(
      sample,
      event_id,
      switchConsequence
    ) %>%
    dplyr::filter(!is.na(switchConsequence)) %>%
    dplyr::distinct()
}

switch_consequence_df <- purrr::map2_dfr(
  iso_switch_lists,
  sample_ids,
  extract_switch_consequences
)


# ----------------------------
# Count switch consequences
# ----------------------------

excluded_consequences <- c(
  "Signal peptide gain",
  "Signal peptide loss",
  "SubCell location gain",
  "SubCell location loss",
  "SubCell location switch"
)

consequence_counts <- switch_consequence_df %>%
  dplyr::filter(!switchConsequence %in% excluded_consequences) %>%
  dplyr::count(switchConsequence, name = "count")


# ----------------------------
# Define paired switch consequences
# ----------------------------

consequence_pairs <- tibble::tribble(
  ~gain,                     ~loss,
  "3UTR is longer",           "3UTR is shorter",
  "5UTR is longer",           "5UTR is shorter",
  "Complete ORF gain",        "Complete ORF loss",
  "Domain gain",              "Domain loss",
  "Exon gain",                "Exon loss",
  "Last exon more downstream", "Last exon more upstream",
  "Length gain",              "Length loss",
  "NMD insensitive",          "NMD sensitive",
  "ORF is longer",            "ORF is shorter",
  "Transcript is coding",     "Transcript is Noncoding",
  "Tss more downstream",      "Tss more upstream",
  "Tts more downstream",      "Tts more upstream"
)


# ----------------------------
# Test directional bias within each consequence pair
# ----------------------------

test_consequence_pair <- function(gain, loss, counts_df) {
  
  count_gain <- counts_df$count[counts_df$switchConsequence == gain]
  count_loss <- counts_df$count[counts_df$switchConsequence == loss]
  
  count_gain <- ifelse(length(count_gain) == 0, 0, count_gain)
  count_loss <- ifelse(length(count_loss) == 0, 0, count_loss)
  
  total_count <- count_gain + count_loss
  
  if (total_count == 0) {
    return(NULL)
  }
  
  test_result <- prop.test(
    x = count_gain,
    n = total_count,
    p = 0.5
  )
  
  tibble(
    consequence = gain,
    paired_consequence = loss,
    proportion = as.numeric(test_result$estimate),
    ci_low = test_result$conf.int[1],
    ci_high = test_result$conf.int[2],
    p_value = test_result$p.value,
    switch_n = total_count
  )
}

consequence_test_results <- purrr::pmap_dfr(
  consequence_pairs,
  test_consequence_pair,
  counts_df = consequence_counts
) %>%
  dplyr::mutate(
    fdr = p.adjust(p_value, method = "fdr"),
    significant = fdr < 0.05,
    consequence = factor(consequence, levels = rev(consequence))
  )


pdf("Switch_Consequences.pdf",width=11,height=7)
ggplot(consequence_test_results, aes(x = proportion, y = consequence)) +
  geom_point(aes(size = switch_n, color = significant))  +  # Size and color based on p-value
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2) +  # Horizontal error bars
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), 
                     name = "FDR < 0.05") +  # Color scale
  scale_size_continuous(name = "switches") +  # Size legend
  scale_x_continuous(limits = c(0, 1)) + 
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "gray") +  # D
  labs(
    x = "Proportion (with 95% Confidence Interval)",
    y = "Comparison (Consequence Pairs)",
    title = "Proportions of Switch Consequences with 95% Confidence Intervals"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )
dev.off()

