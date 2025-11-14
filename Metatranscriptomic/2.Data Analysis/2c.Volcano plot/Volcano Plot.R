###VOLCANO PLOT#
library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(gplots)
library(apeglm)

#Set the working direct

#Load the count data file from the previous step (annotation) called FunctionalAnnotation.txt)
cts <- read.table (file = "FunctionalAnnotation.txt", header = TRUE, sep = "\t",quote = "", stringsAsFactors = FALSE, row.names = 1)

#Load the sample information (I created a txt with the experimental design informations of the samples and named it "ExperimentalDesign.txt")
Information <- read.table (file = "ExperimentalDesign.txt", header = TRUE, sep = "\t")

#Set factor levels 
Information$Vinification <- factor(Information$Vinification)
Information$Stage <- factor(Information$Stage)

#Create a deseq object and import the count data and sample information
dds_ProcFunc<- DESeqDataSetFromMatrix(countData = cts, colData = Information, design = ~ Vinification + Stage)

#Set the reference for the Vinification factor
dds_ProcFunc$Vinification <- factor(dds_ProcFunc$Vinification, levels = c("spontaneous", "inoculated"))
dds_ProcFunc$Stage <- factor(dds_ProcFunc$Stage, levels = c("Start", "Middle", "End"))

#Filter the genes
# keep <- rowSums(counts(dds_ProcFunc)) >= 5

# transform the counts into percentages 
library(vegan)
cnts_perc <- 100*decostand(counts(dds_ProcFunc), MARGIN = 2,method = "total")
cnts_perc_sel <- cnts_perc[which(rowSums(cnts_perc) >= 0.1),]
keep2 <- row.names(cnts_perc_sel)
# dds_ProcFunc <- dds_ProcFunc[keep,]
dds_ProcFunc <- dds_ProcFunc[keep2,]


library(ggplot2)
library(ggrepel)
library(dplyr)
#Perform the statistical test(s) to identify differential expressed genes
dds_ProcFunc <- DESeq(dds_ProcFunc)
saveRDS(dds_ProcFunc, file = "dds_ProcFunc.RDS")


res <- results(dds_ProcFunc, contrast = c("Vinification", "commercial", "inoculated"))

# Order by adjusted p-value
res <- res[order(res$pvalue), ]

# Convert results to a data.frame
res_df <- as.data.frame(res)

# Define significance thresholds
res_df$threshold <- "NotSig"
res_df$threshold[res_df$pvalue < 0.05 & res_df$log2FoldChange > 1]  <- "Up"
res_df$threshold[res_df$pvalue < 0.05 & res_df$log2FoldChange < -1] <- "Down"
res_df$threshold <- factor(res_df$threshold, levels = c("Down", "NotSig", "Up"))

# Select top 10 up and top 10 down genes by smallest padj
top_up <- res_df %>% filter(threshold == "Up") %>% arrange(pvalue) %>% head(10)
top_down <- res_df %>% filter(threshold == "Down") %>% arrange(pvalue) %>% head(10)
top_genes <- rbind(top_up, top_down)

# Volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = threshold)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Down" = "blue", "NotSig" = "grey", "Up" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes, aes(label = rownames(top_genes)),
                  size = 3, max.overlaps = 20) +
  theme_minimal(base_size = 14) +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 p-value",
       color = "Regulation")




