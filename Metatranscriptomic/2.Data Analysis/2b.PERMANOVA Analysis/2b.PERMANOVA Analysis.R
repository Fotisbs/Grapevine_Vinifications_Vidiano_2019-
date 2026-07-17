############################################################
## PERMANOVA of functional transcriptomic profiles
##
## Input files:
##   FinalNames.txt
##     Rows    = functional annotations
##     Columns = samples
##
##   ExperimentalDesign.txt
##     Rows = samples
##
## Data preparation:
##   Counts are converted to relative abundance (%).
##   Functional annotations with a summed relative abundance
##   of at least 0.1% across all samples are retained.
##
## Dissimilarity:
##   Bray-Curtis
##
## Model:
##   Vinification + Stage
############################################################


## =========================================================
## 1. Load packages
## =========================================================

library(phyloseq)
library(vegan)


## =========================================================
## 2. Import transcriptomic count data
## =========================================================

## Rows = functional annotations
## Columns = samples

cts <- read.table(
  file = "FinalNames.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE,
  row.names = 1,
  check.names = FALSE
)

## Check dimensions

dim(cts)

## Expected:
## 80,046 annotations × 12 samples


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

## Verify that all count-table samples are present
## in the metadata

if (!all(colnames(cts) %in% rownames(TheDesign))) {

  missing_samples <- setdiff(
    colnames(cts),
    rownames(TheDesign)
  )

  stop(
    paste(
      "The following samples are missing from",
      "ExperimentalDesign.txt:",
      paste(missing_samples, collapse = ", ")
    )
  )
}

## Reorder metadata to match the count table

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
## 5. Set explanatory variables as factors
## =========================================================

TheDesign$Vinification <- factor(
  TheDesign$Vinification
)

## Check the exact stage names

unique(TheDesign$Stage)

## Set the biological stage order

TheDesign$Stage <- factor(
  TheDesign$Stage,
  levels = c(
    "Start",
    "Middle",
    "End"
  )
)

## Stop if stage labels do not match those in the metadata

if (any(is.na(TheDesign$Stage))) {

  stop(
    paste(
      "Some Stage values do not match Start, Middle and End.",
      "Check unique(TheDesign$Stage) and modify the factor levels."
    )
  )
}


## =========================================================
## 6. Build the phyloseq object
## =========================================================

## Functional annotations are rows and samples are columns

ps_Transcript <- phyloseq(
  otu_table(
    as.matrix(cts),
    taxa_are_rows = TRUE
  ),
  sample_data(TheDesign)
)

ps_Transcript


## =========================================================
## 7. Convert counts to relative abundance (%)
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

## Confirm that each sample sums to approximately 100%

sample_sums(ps_Transcript_100)


## =========================================================
## 8. Filter low-abundance annotations
## =========================================================

## Calculate the summed relative abundance of each annotation
## across all samples

taxa_total <- taxa_sums(
  ps_Transcript_100
)

## Retain annotations with a summed relative abundance
## of at least 0.1% across all samples

ps_Transcript_100_filtered <- prune_taxa(
  taxa_total >= 0.1,
  ps_Transcript_100
)

## Remove zero-abundance annotations, if any remain

ps_Transcript_100_filtered <- prune_taxa(
  taxa_sums(ps_Transcript_100_filtered) > 0,
  ps_Transcript_100_filtered
)


## =========================================================
## 9. Diagnostics
## =========================================================

cat(
  "Number of samples:",
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

## Expected:
##
## Number of samples: 12
## Number of retained annotations: 1095
## Annotations are stored in rows: TRUE


## =========================================================
## 10. Prepare abundance and metadata tables
## =========================================================

## Extract the filtered relative-abundance matrix

transcript_abundance <- as(
  otu_table(ps_Transcript_100_filtered),
  "matrix"
)

## PERMANOVA requires:
## rows    = samples
## columns = annotations
##
## Since taxa_are_rows = TRUE, transpose the matrix.

if (taxa_are_rows(ps_Transcript_100_filtered)) {

  transcript_abundance <- t(
    transcript_abundance
  )
}

## Extract metadata

transcript_metadata <- data.frame(
  sample_data(ps_Transcript_100_filtered),
  check.names = FALSE
)

## Match metadata order to abundance-table rows

transcript_metadata <- transcript_metadata[
  rownames(transcript_abundance),
  ,
  drop = FALSE
]

## Confirm that samples are aligned

stopifnot(
  identical(
    rownames(transcript_abundance),
    rownames(transcript_metadata)
  )
)

## Final dimensions should be:
## 12 samples × 1095 annotations

dim(transcript_abundance)


## =========================================================
## 11. Bray-Curtis PERMANOVA
## =========================================================

set.seed(123)

permanova_transcript <- adonis2(
  transcript_abundance ~ Vinification + Stage,
  data = transcript_metadata,
  method = "bray",
  permutations = 999,
  by = "margin"
)

## Display results

permanova_transcript


## =========================================================
## 12. Export PERMANOVA results
## =========================================================

permanova_output <- data.frame(
  Term = rownames(permanova_transcript),
  permanova_transcript,
  row.names = NULL,
  check.names = FALSE
)

write.table(
  permanova_output,
  file = "PERMANOVA_Transcript_Vinification_Stage.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


## =========================================================
## 13. Optional: save analysis objects
## =========================================================

saveRDS(
  ps_Transcript_100_filtered,
  file = "ps_Transcript_100_filtered.rds"
)

saveRDS(
  permanova_transcript,
  file = "PERMANOVA_Transcript_results.rds"
)


## =========================================================
## 14. Session information
## =========================================================

sessionInfo()
