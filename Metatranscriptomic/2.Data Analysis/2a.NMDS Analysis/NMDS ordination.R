############################################################
## Nature-style NMDS of functional transcriptomic profiles
##
## Input files:
##   FinalNames.txt
##     Rows    = functional annotations
##     Columns = samples
##
##   ExperimentalDesign.txt
##     Rows = samples
##
## Filtering:
##   Counts are converted to sample-wise relative abundance (%).
##   Functional annotations with a summed relative abundance
##   of at least 0.1% across all samples are retained.
##
## Distance:
##   Bray-Curtis
############################################################


## =========================================================
## 1. Load packages
## =========================================================

library(phyloseq)
library(vegan)
library(RColorBrewer)


## =========================================================
## 2. Import transcriptomic count data
## =========================================================

## Rows = 80,046 functional annotations
## Columns = 12 samples

cts <- read.table(
  file = "FinalNames.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE,
  row.names = 1,
  check.names = FALSE
)

## Check raw count-table dimensions

dim(cts)

## Expected:
## 80046 rows × 12 columns


## =========================================================
## 3. Import experimental metadata
## =========================================================

TheDesign <- read.table(
  file = "ExperimentalDesign.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE,
  row.names = 1,
  check.names = FALSE
)

## Remove accidental leading or trailing spaces

TheDesign$Vinification <- trimws(
  TheDesign$Vinification
)

TheDesign$Stage <- trimws(
  TheDesign$Stage
)


## =========================================================
## 4. Check and match sample names
## =========================================================

## Check whether every sample in the count table occurs
## in the experimental metadata

if (!all(colnames(cts) %in% rownames(TheDesign))) {

  missing_samples <- setdiff(
    colnames(cts),
    rownames(TheDesign)
  )

  stop(
    paste(
      "The following count-table samples are missing",
      "from ExperimentalDesign.txt:",
      paste(missing_samples, collapse = ", ")
    )
  )
}

## Match metadata order to the count-table sample order

TheDesign <- TheDesign[
  colnames(cts),
  ,
  drop = FALSE
]

## Confirm identical sample order

stopifnot(
  identical(
    colnames(cts),
    rownames(TheDesign)
  )
)


## =========================================================
## 5. Build the phyloseq object
## =========================================================

## Functional annotations are rows.
## Samples are columns.

ps_Transcript <- phyloseq(
  otu_table(
    as.matrix(cts),
    taxa_are_rows = TRUE
  ),
  sample_data(TheDesign)
)

ps_Transcript


## =========================================================
## 6. Convert counts to relative abundance (%)
## =========================================================

ps_Transcript_100 <- transform_sample_counts(
  ps_Transcript,
  function(x) {

    if (sum(x) == 0) {
      return(x)
    }

    100 * x / sum(x)
  }
)

## Check that each sample sums to approximately 100%

sample_sums(ps_Transcript_100)


## =========================================================
## 7. Filter low-abundance functional annotations
## =========================================================

## Calculate the summed relative abundance of each annotation
## across all 12 samples

taxa_total <- taxa_sums(
  ps_Transcript_100
)

## Retain annotations with a summed relative abundance
## of at least 0.1% across all samples

ps_Transcript_100_filtered <- prune_taxa(
  taxa_total >= 0.1,
  ps_Transcript_100
)

## Remove any zero-abundance annotations, if present

ps_Transcript_100_filtered <- prune_taxa(
  taxa_sums(ps_Transcript_100_filtered) > 0,
  ps_Transcript_100_filtered
)


## =========================================================
## 8. Diagnostics
## =========================================================

ps_Transcript_100_filtered

cat(
  "\nNumber of samples:",
  nsamples(ps_Transcript_100_filtered),
  "\n"
)

cat(
  "Number of retained annotations:",
  ntaxa(ps_Transcript_100_filtered),
  "\n"
)

cat(
  "Annotations are stored in rows:",
  taxa_are_rows(ps_Transcript_100_filtered),
  "\n"
)

cat(
  "Sample variables:",
  paste(
    sample_variables(ps_Transcript_100_filtered),
    collapse = ", "
  ),
  "\n"
)

## Expected:
##
## Number of samples: 12
## Number of retained annotations: 1095
## Annotations are stored in rows: TRUE


## =========================================================
## 9. Calculate Bray-Curtis NMDS
## =========================================================

set.seed(123)

ord_nmds_transcript <- ordinate(
  physeq = ps_Transcript_100_filtered,
  method = "NMDS",
  distance = "bray",
  k = 2,
  trymax = 999,
  autotransform = FALSE,
  trace = FALSE
)

cat(
  "\nNMDS stress:",
  ord_nmds_transcript$stress,
  "\n"
)

## Expected stress:
## approximately 0.075


## =========================================================
## 10. Extract NMDS sample coordinates
## =========================================================

nmds_sites <- scores(
  ord_nmds_transcript,
  display = "sites"
)

nmds_sites <- as.matrix(
  nmds_sites
)

## Check the coordinates

nmds_sites


## =========================================================
## 11. Extract and match metadata
## =========================================================

metadata <- data.frame(
  sample_data(ps_Transcript_100_filtered),
  check.names = FALSE
)

metadata <- metadata[
  rownames(nmds_sites),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(nmds_sites),
    rownames(metadata)
  )
)


## =========================================================
## 12. Set factor levels
## =========================================================

## Set vinification order.
## Modify the levels below only if your metadata uses
## different names.

metadata$Vinification <- factor(
  metadata$Vinification
)

## Check stage names before assigning their order

unique(metadata$Stage)

## Set the biologically meaningful stage order.
## Replace these names if ExperimentalDesign.txt uses
## alternative labels.

metadata$Stage <- factor(
  metadata$Stage,
  levels = c(
    "Start",
    "Middle",
    "End"
  )
)

## Stop if the specified stage names did not match the data

if (any(is.na(metadata$Stage))) {

  stop(
    paste(
      "Some Stage values did not match the specified levels.",
      "Check unique(TheDesign$Stage) and update the factor levels."
    )
  )
}


## =========================================================
## 13. Nature-style NMDS plotting function
## =========================================================

plot_nmds_nature <- function(
    nmds_object,
    site_scores,
    grouping_variable,
    output_file,
    panel_label,
    legend_title = NULL,
    palette_name = "RdBu",
    width = 7,
    height = 6
) {

  ## Convert grouping variable to factor

  grouping_variable <- factor(
    grouping_variable
  )

  ## Check correspondence between samples and groups

  if (length(grouping_variable) != nrow(site_scores)) {

    stop(
      "The number of grouping values does not match the number of samples."
    )
  }

  ## Check for missing grouping values

  if (any(is.na(grouping_variable))) {

    stop(
      paste(
        "The grouping variable contains missing values.",
        "Check the metadata factor levels."
      )
    )
  }

  ## Number of groups

  number_groups <- nlevels(
    grouping_variable
  )

  ## RColorBrewer palettes generally require at least
  ## three colours

  palette_size <- max(
    3,
    number_groups
  )

  ## Check that the requested number of colours is supported

  maximum_colours <- RColorBrewer::brewer.pal.info[
    palette_name,
    "maxcolors"
  ]

  if (palette_size > maximum_colours) {

    stop(
      paste0(
        "Palette ",
        palette_name,
        " supports a maximum of ",
        maximum_colours,
        " colours."
      )
    )
  }

  ## Prepare group colours

  plot_colours <- RColorBrewer::brewer.pal(
    n = palette_size,
    name = palette_name
  )[seq_len(number_groups)]

  ## Convert factor groups to numeric indices

  group_numbers <- as.numeric(
    grouping_variable
  )

  ## Calculate plotting ranges

  x_range <- range(
    site_scores[, 1],
    na.rm = TRUE
  )

  y_range <- range(
    site_scores[, 2],
    na.rm = TRUE
  )

  ## Add 15% plotting space around sample coordinates

  x_padding <- diff(x_range) * 0.15
  y_padding <- diff(y_range) * 0.15

  if (x_padding == 0) {
    x_padding <- 0.1
  }

  if (y_padding == 0) {
    y_padding <- 0.1
  }

  ## Open PDF graphics device

  cairo_pdf(
    filename = output_file,
    width = width,
    height = height
  )

  ## Ensure the device closes even if plotting fails

  on.exit(
    dev.off(),
    add = TRUE
  )

  ## Set graphical parameters

  par(
    mar = c(5, 5, 2, 2),
    mgp = c(2.8, 0.8, 0),
    las = 1,
    family = "sans"
  )

  ## Draw an empty NMDS plotting area

  plot(
    site_scores,
    type = "n",
    frame.plot = FALSE,
    xlab = "NMDS1",
    ylab = "NMDS2",
    xlim = c(
      x_range[1] - x_padding,
      x_range[2] + x_padding
    ),
    ylim = c(
      y_range[1] - y_padding,
      y_range[2] + y_padding
    ),
    cex.lab = 1.2,
    cex.axis = 1
  )

  ## Add convex-hull ellipses before points so that
  ## they remain in the background

  vegan::ordiellipse(
    site_scores,
    groups = grouping_variable,
    kind = "ehull",
    draw = "lines",
    lty = 2,
    lwd = 1
  )

  ## Add sample points

  points(
    site_scores,
    bg = plot_colours[group_numbers],
    col = "black",
    pch = 21,
    cex = 1.6,
    lwd = 0.8
  )

  ## Add NMDS stress to lower-left corner

  par(adj = 0)

  title(
    sub = paste0(
      "stress ",
      sprintf(
        "%.2f",
        nmds_object$stress
      )
    ),
    cex.sub = 1.1,
    line = 3
  )

  ## Add panel label to lower-right corner

  par(adj = 1)

  title(
    sub = panel_label,
    cex.sub = 1.1,
    line = 3
  )

  par(adj = 0.5)

  ## Add legend

  graphics::legend(
    "topright",
    bty = "n",
    title = legend_title,
    legend = levels(grouping_variable),
    pch = 21,
    pt.bg = plot_colours,
    col = "black",
    pt.cex = 1.5,
    cex = 1
  )
}


## =========================================================
## 14. Plot NMDS according to vinification strategy
## =========================================================

plot_nmds_nature(
  nmds_object = ord_nmds_transcript,
  site_scores = nmds_sites,
  grouping_variable = metadata$Vinification,
  output_file = "NMDS_Transcript_Vinification_NatureStyle.pdf",
  panel_label = "Transcriptomic profiles",
  legend_title = "Vinification",
  palette_name = "RdBu",
  width = 7,
  height = 6
)


## =========================================================
## 15. Plot NMDS according to fermentation stage
## =========================================================

plot_nmds_nature(
  nmds_object = ord_nmds_transcript,
  site_scores = nmds_sites,
  grouping_variable = metadata$Stage,
  output_file = "NMDS_Transcript_Stage_NatureStyle.pdf",
  panel_label = "Transcriptomic profiles",
  legend_title = "Fermentation stage",
  palette_name = "RdBu",
  width = 7,
  height = 6
)


## =========================================================
## 16. Optional: export NMDS coordinates
## =========================================================

nmds_coordinates <- data.frame(
  Sample = rownames(nmds_sites),
  NMDS1 = nmds_sites[, 1],
  NMDS2 = nmds_sites[, 2],
  Vinification = metadata$Vinification,
  Stage = metadata$Stage,
  row.names = NULL
)

write.table(
  nmds_coordinates,
  file = "NMDS_Transcript_Coordinates.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


## =========================================================
## 17. Optional: save the filtered phyloseq object
## =========================================================

saveRDS(
  ps_Transcript_100_filtered,
  file = "ps_Transcript_100_filtered.rds"
)

saveRDS(
  ord_nmds_transcript,
  file = "ord_nmds_transcript.rds"
)


## =========================================================
## 18. Session information
## =========================================================

sessionInfo()
