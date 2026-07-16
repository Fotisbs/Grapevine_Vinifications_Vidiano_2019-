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
# Gene-selection criteria:
#   Benjamini–Hochberg adjusted p-value < 0.05
#   Absolute log2 fold-change > 1
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
Information$Vinification <- trimws(Information$Vinification)
Information$Stage <- trimws(Information$Stage)

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
# 7. Select significant genes
# ----------------------------------------------------------

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

if (nrow(significant_genes) == 0) {
  stop(
    "No genes met the criteria padj < 0.05 ",
    "and |log2FC| > 1 at S3."
  )
}

# Select up to 20 most significant genes
top_genes <- significant_genes %>%
  slice_head(
    n = min(20, nrow(significant_genes))
  ) %>%
  pull(
    Gene
  )

# ----------------------------------------------------------
# 8. Extract DESeq2-normalized counts
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

write.table(
  normalized_counts,
  file = "S3_DESeq2_normalized_counts.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# ----------------------------------------------------------
# 9. Apply row-wise range standardization
# ----------------------------------------------------------

heatmap_matrix <- decostand(
  selected_normalized_counts,
  MARGIN = 1,
  method = "range"
)

heatmap_matrix <- heatmap_matrix[
  complete.cases(heatmap_matrix),
  ,
  drop = FALSE
]

# ----------------------------------------------------------
# 10. Prepare sample annotation
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
# 11. Generate the heatmap
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