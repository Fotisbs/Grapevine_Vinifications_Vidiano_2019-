##PERMANOVA##
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

mypermanova_ps_Transcript_100_filtered <- adonis2(
  ps_Transcript_100_filtered@otu_table ~ Vinification + Stage,
  data = data.frame(ps_Transcript_100_filtered@sam_data),
  method = "bray",
  by = "margin"
)

write.table(data.frame(mypermanova_ps_Transcript_100_filtered), file="mypermanova_ps_Transcript_100_filtered.txt", quote = F,col.names = NA, sep="\t")
