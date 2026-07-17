############################################################
# Figure 4: Volcano plot of differential gene expression
#
# Global DESeq2 model:
#   ~ Vinification + Stage
#
# Contrast:
#   Inoculated versus Spontaneous
#
# Positive log2FC:
#   Higher expression in inoculated fermentations
#
# Negative log2FC:
#   Higher expression in spontaneous fermentations
#
# Volcano-plot thresholds:
#   nominal p-value < 0.05
#   absolute log2 fold-change > 1
############################################################

library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(vegan)

# ----------------------------------------------------------
# 1. Import gene-count table
# ----------------------------------------------------------

cts <- read.table(
  file = "FinalNames.txt",
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
  file = "TheDesign.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  row.names = 1,
  check.names = FALSE
)

# Match metadata order to count-table columns
if (!all(colnames(cts) %in% rownames(Information))) {
  stop(
    "Sample names in FinalNames.txt and TheDesign.txt do not match."
  )
}

Information <- Information[
  colnames(cts),
  ,
  drop = FALSE
]

# ----------------------------------------------------------
# 3. Define factor levels
# ----------------------------------------------------------

Information$Vinification <- factor(
  Information$Vinification,
  levels = c(
    "spontaneous",
    "inoculated"
  )
)

Information$Stage <- factor(
  Information$Stage,
  levels = c(
    "Start",
    "Middle",
    "End"
  )
)

if (anyNA(Information$Vinification)) {
  stop(
    "Unexpected values were found in the Vinification column."
  )
}

if (anyNA(Information$Stage)) {
  stop(
    "Unexpected values were found in the Stage column."
  )
}

# ----------------------------------------------------------
# 4. Construct DESeq2 object
# ----------------------------------------------------------

dds_ProcFunc <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(cts)),
  colData = Information,
  design = ~ Vinification + Stage
)

# ----------------------------------------------------------
# 5. Filter genes
# ----------------------------------------------------------

# Convert counts to sample-wise relative abundance (%)
cnts_perc <- 100 * decostand(
  counts(dds_ProcFunc),
  MARGIN = 2,
  method = "total"
)

# Retain genes whose summed relative contribution across all
# samples is at least 0.1 percentage points
keep_genes <- rowSums(
  cnts_perc,
  na.rm = TRUE
) >= 0.1

dds_ProcFunc <- dds_ProcFunc[
  keep_genes,
]

message(
  sum(keep_genes),
  " genes retained from ",
  length(keep_genes),
  " annotated genes."
)

# ----------------------------------------------------------
# 6. Run DESeq2
# ----------------------------------------------------------

dds_ProcFunc <- DESeq(
  dds_ProcFunc
)

saveRDS(
  dds_ProcFunc,
  file = "dds_ProcFunc_global_model.rds"
)

# The contrast is:
# log2(inoculated / spontaneous)
#
# If "inoculated" corresponds to the inoculated fermentation,
# positive values indicate higher expression in inoculated
# fermentations.

res <- results(
  dds_ProcFunc,
  contrast = c(
    "Vinification",
    "inoculated",
    "spontaneous"
  )
)

# Convert results to a data frame
res_df <- as.data.frame(
  res
) %>%
  tibble::rownames_to_column(
    "Gene"
  ) %>%
  filter(
    !is.na(log2FoldChange),
    !is.na(pvalue)
  )

# ----------------------------------------------------------
# 7. Define volcano-plot categories
# ----------------------------------------------------------

res_df <- res_df %>%
  mutate(
    Expression_pattern = case_when(
      pvalue < 0.05 &
        log2FoldChange > 1 ~
        "Higher in Inoculated",
      
      pvalue < 0.05 &
        log2FoldChange < -1 ~
        "Higher in Spontaneous",
      
      TRUE ~
        "Not significant"
    ),
    Expression_pattern = factor(
      Expression_pattern,
      levels = c(
        "Higher in Spontaneous",
        "Not significant",
        "Higher in Inoculated"
      )
    )
  )

# ----------------------------------------------------------
# 8. Select ten genes in each direction
# ----------------------------------------------------------

top_inoculated <- res_df %>%
  filter(
    Expression_pattern ==
      "Higher in Inoculated"
  ) %>%
  arrange(
    pvalue
  ) %>%
  slice_head(
    n = 10
  )

top_spontaneous <- res_df %>%
  filter(
    Expression_pattern ==
      "Higher in Spontaneous"
  ) %>%
  arrange(
    pvalue
  ) %>%
  slice_head(
    n = 10
  )

top_genes <- bind_rows(
  top_spontaneous,
  top_inoculated
)

# ----------------------------------------------------------
# 9. Generate volcano plot
# ----------------------------------------------------------

volcano_plot <- ggplot(
  res_df,
  aes(
    x = log2FoldChange,
    y = -log10(pvalue),
    color = Expression_pattern
  )
) +
  geom_point(
    alpha = 0.7,
    size = 2
  ) +
  scale_color_manual(
    values = c(
      "Higher in Spontaneous" = "blue",
      "Not significant" = "grey",
      "Higher in Inoculated" = "red"
    )
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "black"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "black"
  ) +
  geom_text_repel(
    data = top_genes,
    aes(
      label = Gene
    ),
    size = 3,
    max.overlaps = 20
  ) +
  theme_minimal(
    base_size = 14
  ) +
  labs(
    title = "Volcano Plot",
    x = expression(
      log[2](
        "Inoculated / Spontaneous"
      )
    ),
    y = expression(
      -log[10](
        "nominal p-value"
      )
    ),
    color = ""
  )

print(
  volcano_plot
)

ggsave(
  filename = "Figure4_VolcanoPlot.pdf",
  plot = volcano_plot,
  width = 9,
  height = 7
)

# ----------------------------------------------------------
# 10. Export the complete DESeq2 results
# ----------------------------------------------------------

write.table(
  res_df,
  file = "Figure4_DESeq2_results.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
