############################################################
# Figure 5c: Stage-specific transcript heatmap — End (S3)
#
# Comparison:
#   Inoculated versus Spontaneous fermentations at S3
#
# Biological replicates:
#   Two independent tanks per fermentation strategy
#
# DESeq2 design:
#   ~ Vinification
#
# Differential-expression criteria:
#   Benjamini–Hochberg adjusted p-value < 0.05
#   Absolute log2 fold-change > 1
#
# Heatmap gene selection:
#   All tested genes were ranked according to their
#   Benjamini–Hochberg adjusted p-values, and the 20
#   top-ranked genes were selected for visualization.
#   Heatmap selection was therefore based on ranking by
#   adjusted p-value and did not require all displayed genes
#   to satisfy the differential-expression thresholds above.
#
# Heatmap values:
#   DESeq2-normalized counts scaled independently for each
#   gene using row-wise range standardization (0–1)
#
# Input files:
#   Final_names_end.txt
#   TheDesign2ViniEnd.txt
#
# Output:
#   Figure5c_S3_GeneExpressionHeatmap.pdf
############################################################

library(DESeq2)
library(dplyr)
library(pheatmap)
library(vegan)

# ----------------------------------------------------------
# 1. Import gene-count data
# ----------------------------------------------------------

cts <- read.table(
  file = "Final_names_end.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE,
  row.names = 1,
  check.names = FALSE
)

# ----------------------------------------------------------
# 2. Import sample metadata
# ----------------------------------------------------------

Information <- read.table(
  file = "TheDesign2ViniEnd.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  row.names = 1,
  check.names = FALSE
)

# Remove accidental leading or trailing spaces
Information$Vinification <- trimws(
  Information$Vinification
)

Information$Stage <- trimws(
  Information$Stage
)

# Confirm that metadata sample IDs match count-table columns
if (!all(colnames(cts) %in% rownames(Information))) {
  stop(
    "Sample names in Final_names_end.txt and ",
    "TheDesign2ViniEnd.txt do not match."
  )
}

# Reorder metadata to match count-table columns
Information <- Information[
  colnames(cts),
  ,
  drop = FALSE
]

# ----------------------------------------------------------
# 3. Define the fermentation factor
# ----------------------------------------------------------

Information$Vinification <- factor(
  Information$Vinification,
  levels = c(
    "spontaneous",
    "inoculated"
  )
)

if (anyNA(Information$Vinification)) {
  stop(
    "Unexpected values were detected in the ",
    "Vinification column."
  )
}

# ----------------------------------------------------------
# 4. Construct the DESeq2 dataset
# ----------------------------------------------------------

dds_2VinEnd_Func <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(cts)),
  colData = Information,
  design = ~ Vinification
)

# ----------------------------------------------------------
# 5. Filter genes by relative contribution
# ----------------------------------------------------------

# Convert counts to sample-wise relative abundance (%)
cnts_perc <- 100 * decostand(
  counts(dds_2VinEnd_Func),
  MARGIN = 2,
  method = "total"
)

# Retain genes whose summed relative contribution across
# the four S3 samples is at least 0.1 percentage points
keep_genes <- rowSums(
  cnts_perc,
  na.rm = TRUE
) >= 0.1

dds_2VinEnd_Func <- dds_2VinEnd_Func[
  keep_genes,
]

message(
  sum(keep_genes),
  " genes retained from ",
  length(keep_genes),
  " annotated genes for the S3 analysis."
)

# ----------------------------------------------------------
# 6. Run stage-specific DESeq2 analysis
# ----------------------------------------------------------

dds_2VinEnd_Func <- DESeq(
  dds_2VinEnd_Func
)

saveRDS(
  dds_2VinEnd_Func,
  file = "dds_2VinEnd_Func.rds"
)

# Explicit contrast:
# log2(Inoculated / Spontaneous)
deseq_results <- results(
  dds_2VinEnd_Func,
  contrast = c(
    "Vinification",
    "inoculated",
    "spontaneous"
  ),
  alpha = 0.05
)

deseq_results_df <- as.data.frame(
  deseq_results
) %>%
  tibble::rownames_to_column(
    "Gene"
  )

# Export all DESeq2 results
write.table(
  deseq_results_df,
  file = "S3_DESeq2_results_all.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------------------------------------
# 7. Identify genes meeting the DE significance criteria
# ----------------------------------------------------------

# These thresholds define differential expression but are
# NOT used to restrict the genes entering the heatmap.
significant_genes <- deseq_results_df %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    !is.na(log2FoldChange),
    abs(log2FoldChange) > 1
  ) %>%
  arrange(
    padj
  )

write.table(
  significant_genes,
  file = "S3_DESeq2_results_significant.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(
  nrow(significant_genes),
  " genes met padj < 0.05 and |log2FC| > 1 at S3."
)

# ----------------------------------------------------------
# 8. Select the 20 top-ranked genes for the heatmap
# ----------------------------------------------------------

# Rank all tested genes by Benjamini–Hochberg adjusted
# p-value. Genes with NA padj values are placed last.
#
# IMPORTANT:
# Heatmap selection is based on ranking by padj and does not
# require the selected genes to pass padj < 0.05 and
# |log2FC| > 1.

heatmap_ranked_genes <- deseq_results_df %>%
  arrange(
    is.na(padj),
    padj
  )

top_genes <- heatmap_ranked_genes %>%
  slice_head(
    n = min(20, nrow(heatmap_ranked_genes))
  ) %>%
  pull(
    Gene
  )

# Export the genes selected for Figure 5c
heatmap_top20_results <- heatmap_ranked_genes %>%
  slice_head(
    n = min(20, nrow(heatmap_ranked_genes))
  )

write.table(
  heatmap_top20_results,
  file = "S3_DESeq2_top20_by_adjusted_pvalue.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------------------------------------
# 9. Extract DESeq2-normalized counts
# ----------------------------------------------------------

normalized_counts <- counts(
  dds_2VinEnd_Func,
  normalized = TRUE
)

selected_normalized_counts <- normalized_counts[
  top_genes,
  ,
  drop = FALSE
]

# Export all normalized counts
write.table(
  normalized_counts,
  file = "S3_DESeq2_normalized_counts.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# Export normalized counts for the 20 heatmap genes
write.table(
  selected_normalized_counts,
  file = "S3_Figure5c_top20_normalized_counts.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# ----------------------------------------------------------
# 10. Apply row-wise range standardization
# ----------------------------------------------------------

heatmap_matrix <- decostand(
  selected_normalized_counts,
  MARGIN = 1,
  method = "range"
)

# Remove genes with undefined scaled values, if any
heatmap_matrix <- heatmap_matrix[
  complete.cases(heatmap_matrix),
  ,
  drop = FALSE
]

# Export the actual 0–1 values displayed in the heatmap
write.table(
  heatmap_matrix,
  file = "S3_Figure5c_top20_range_scaled_values.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# ----------------------------------------------------------
# 11. Prepare sample annotation
# ----------------------------------------------------------

annotation_col <- data.frame(
  Vinification = colData(
    dds_2VinEnd_Func
  )$Vinification
)

rownames(annotation_col) <- colnames(
  heatmap_matrix
)

# ----------------------------------------------------------
# 12. Generate the heatmap
# ----------------------------------------------------------

pdf(
  file = "Figure5c_S3_GeneExpressionHeatmap.pdf",
  width = 10,
  height = 10
)

pheatmap(
  heatmap_matrix,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 10,
  border_color = "grey30",
  cellwidth = 18,
  cellheight = 15,
  legend_breaks = c(0, 1),
  legend_labels = c(
    "Minimum per gene",
    "Maximum per gene"
  ),
  main = "End"
)

dev.off()
```
