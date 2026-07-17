############################################################
## NMDS plots for transcriptomic profiles
##
## Input object:
## ps_Transcript_100_filtered
##
## Data:
## 12 samples × 1,095 retained functional annotations
##
## Distance:
## Bray-Curtis
############################################################

library(phyloseq)
library(vegan)
library(RColorBrewer)

## =========================================================
## 1. Build and Check the phyloseq object
## =========================================================

#Build the object with transcriptomic data (80046 Genes, counts)

cts <- read.table (file = "FinalNames.txt", header = TRUE, sep = "\t",quote = "", stringsAsFactors = FALSE, row.names = 1)
View(data.frame(cts))

TheDesign <- read.table (file = "ExperimentalDesign.txt", header = TRUE, sep = "\t",quote = "", stringsAsFactors = FALSE, row.names = 1)
View(data.frame(TheDesign))

ps_Transcript

ps_Transcript <- phyloseq(otu_table(cts, taxa_are_rows=FALSE), sample_data(TheDesign))

##Transform RA 100%
ps_Transcript_100 <- transform_sample_counts(ps_Transcript, function(OTU) 100*OTU/sum(OTU))

taxa_total <- taxa_sums(ps_Transcript_100)

#Filter Genes and keep those contributing at least 0.1%
ps_Transcript_100_filtered <- prune_taxa(taxa_total >= 0.1, ps_Transcript_100)

ps_Transcript_100_filtered



ps_Transcript_100_filtered

nsamples(ps_Transcript_100_filtered)  # Expected: 12
ntaxa(ps_Transcript_100_filtered)     # Expected: 1095
taxa_are_rows(ps_Transcript_100_filtered)  # Expected: TRUE

sample_variables(ps_Transcript_100_filtered)

## =========================================================
## 2. Calculate NMDS
## =========================================================

set.seed(123)

ord_nmds_transcript <- ordinate(
  ps_Transcript_100_filtered,
  method = "NMDS",
  distance = "bray",
  k = 2,
  trymax = 999,
  autotransform = FALSE
)

ord_nmds_transcript$stress

## Sample coordinates
nmds_sites <- scores(
  ord_nmds_transcript,
  display = "sites"
)

## Ensure sample order matches metadata
metadata <- data.frame(
  sample_data(ps_Transcript_100_filtered)
)

metadata <- metadata[
  rownames(nmds_sites),
  ,
  drop = FALSE
]

## =========================================================
## 3. Reusable plotting function
## =========================================================

plot_nmds_nature <- function(
    nmds_object,
    nmds_sites,
    grouping_variable,
    output_file,
    panel_label,
    legend_title = NULL,
    width = 7,
    height = 6
) {
  
  ## Convert grouping variable to factor
  grouping_variable <- factor(grouping_variable)
  
  ## Number of groups
  number_groups <- nlevels(grouping_variable)
  
  ## Prepare palette
  ## brewer.pal requires at least three colours
  palette_size <- max(4, number_groups)
  
  plot_colours <- RColorBrewer::brewer.pal(
    n = palette_size,
    name = "RdBu"
  )[seq_len(number_groups)]
  
  ## Numeric group codes for point colours
  group_numbers <- as.numeric(grouping_variable)
  
  ## Add some space around the plotted points
  x_range <- range(nmds_sites[, 1], na.rm = TRUE)
  y_range <- range(nmds_sites[, 2], na.rm = TRUE)
  
  x_padding <- diff(x_range) * 0.15
  y_padding <- diff(y_range) * 0.15
  
  ## Protect against a zero coordinate range
  if (x_padding == 0) {
    x_padding <- 0.1
  }
  
  if (y_padding == 0) {
    y_padding <- 0.1
  }
  
  ## Open PDF device
  cairo_pdf(
    filename = output_file,
    width = width,
    height = height
  )
  
  ## Plot margins and text
  par(
    mar = c(5, 5, 2, 2),
    mgp = c(2.8, 0.8, 0),
    las = 1,
    family = "sans"
  )
  
  ## Empty plot
  plot(
    nmds_sites,
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
  
  ## Add group ellipses in the background
  vegan::ordiellipse(
    nmds_sites,
    groups = grouping_variable,
    kind = "ehull",
    draw = "lines",
    lty = 2,
    lwd = 1
  )
  
  ## Add samples
  points(
    nmds_sites,
    bg = plot_colours[group_numbers],
    col = "black",
    pch = 21,
    cex = 1.6,
    lwd = 0.8
  )
  
  ## Add stress in lower-left corner
  par(adj = 0)
  
  title(
    sub = paste0(
      "stress ",
      sprintf("%.2f", nmds_object$stress)
    ),
    cex.sub = 1.1,
    line = 3
  )
  
  ## Add panel label in lower-right corner
  par(adj = 1)
  
  title(
    sub = panel_label,
    cex.sub = 1.1,
    line = 3
  )
  
  par(adj = 0.5)
  
  ## Legend
  if (is.null(legend_title)) {
    
    graphics::legend(
      "topright",
      bty = "n",
      legend = levels(grouping_variable),
      pch = 21,
      pt.bg = plot_colours,
      pt.cex = 1.5,
      cex = 1
    )
    
  } else {
    
    graphics::legend(
      "topright",
      bty = "n",
      title = legend_title,
      legend = levels(grouping_variable),
      pch = 21,
      pt.bg = plot_colours,
      pt.cex = 1.5,
      cex = 1
    )
  }
  
  dev.off()
}

## =========================================================
## 4. NMDS according to vinification strategy
## =========================================================

plot_nmds_nature(
  nmds_object = ord_nmds_transcript,
  nmds_sites = nmds_sites,
  grouping_variable = metadata$Vinification,
  output_file = "NMDS_Transcript_Vinification_NatureStyle.pdf",
  panel_label = "Transcriptomic profiles",
  legend_title = "Vinification"
)

## =========================================================
## 5. NMDS according to fermentation stage
## =========================================================

## Set the desired stage order
metadata$Stage <- factor(
  metadata$Stage,
  levels = c(
    "Start",
    "Middle",
    "End"
  )
)

plot_nmds_nature(
  nmds_object = ord_nmds_transcript,
  nmds_sites = nmds_sites,
  grouping_variable = metadata$Stage,
  output_file = "NMDS_Transcript_Stage_NatureStyle.pdf",
  panel_label = "Transcriptomic profiles",
  legend_title = "Fermentation stage"
)

