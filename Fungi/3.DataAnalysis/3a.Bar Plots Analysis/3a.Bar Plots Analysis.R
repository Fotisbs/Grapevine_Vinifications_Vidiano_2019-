##Start the analysis by Merging the replicates of each sample per stage of Spontaneous and Inoculated Vinifications##
MUST_Fungi_Vidiano_2019_Merged <- merge_samples(fungi_vinification_Annotated, "stage")

MUST_Fungi_Vidiano_2019_Merged <- prune_taxa(taxa_sums(MUST_Fungi_Vidiano_2019_Merged)>0,MUST_Fungi_Vidiano_2019_Merged)

View(data.frame(sample_data(MUST_Fungi_Vidiano_2019_Merged)))
View(data.frame(otu_table(MUST_Fungi_Vidiano_2019_Merged)))
View(data.frame(tax_table(MUST_Fungi_Vidiano_2019_Merged)))

##Rename sample data file after merging
write.table(data.frame(sample_data(MUST_Fungi_Vidiano_2019_Merged)), file="MUST_Fungi_Vidiano_2019_Merged.txt", quote = F,col.names = NA, sep="\t")

SampleDataNew67 <- read.table("MUST_Fungi_Vidiano_2019_Merged.txt", header=T,sep = "\t",row.names = 1)

sample_data(MUST_Fungi_Vidiano_2019_Merged) <- SampleDataNew67

MUST_Fungi_Vidiano_2019_Merged

View(data.frame(sample_data(MUST_Fungi_Vidiano_2019_Merged)))

##Transform RA 100%
MUST_Fungi_Vidiano_2019_Merged #(Raw Data Merged)

MUST_Fungi_Vidiano_2019_Merged_100 <- transform_sample_counts(MUST_Fungi_Vidiano_2019_Merged, function(OTU) 100*OTU/sum(OTU))

MUST_Fungi_Vidiano_2019_Merged_100 #(Transformation 100%)

#Most abundant taxa (Top 20 at species level) 
MUST_Fungi_Vidiano_2019_Merged
MUST_Fungi_Vidiano_2019_Merged_100

##tax glom at species level first#
rank_names(MUST_Fungi_Vidiano_2019_Merged_100)

MUST_Fungi_Vidiano_2019_Merged_100

MUST_Fungi_Species <- tax_glom(MUST_Fungi_Vidiano_2019_Merged_100, taxrank = "Species")

MUST_Fungi_Species

MUST_Fungi_Species <- prune_taxa(taxa_sums(MUST_Fungi_Species)>0,MUST_Fungi_Species)

##top20
myTaxa20_MUST_Fungi_Species <- names(sort(taxa_sums(MUST_Fungi_Species), decreasing = TRUE)[1:20])  

Top20_MUST_Fungi_Species <- prune_taxa(myTaxa20_MUST_Fungi_Species, MUST_Fungi_Species)

taxa_names(Top20_MUST_Fungi_Species)

Top20_MUST_Fungi_Species <- prune_taxa(taxa_sums(Top20_MUST_Fungi_Species)>0,Top20_MUST_Fungi_Species)

mytax20_MUST_Fungi_Species <- data.frame(tax_table(Top20_MUST_Fungi_Species), stringsAsFactors = F)

# For ITS - Remove letter from taxonomy
for (i in c(1:nrow(mytax20_MUST_Fungi_Species))) {
  for(j in c(1:ncol(mytax20_MUST_Fungi_Species))) {
    mytax20_MUST_Fungi_Species[i,j] <- gsub("[a-z]__","",mytax20_MUST_Fungi_Species[i,j])
  }
}

mytxplot20_MUST_Fungi_Species <- data.frame(OTU = row.names(mytax20_MUST_Fungi_Species), 
                                            txplt = paste(row.names(mytax20_MUST_Fungi_Species), " ", mytax20_MUST_Fungi_Species$Genus,  ":", mytax20_MUST_Fungi_Species$Species,  sep = ""))

row.names(mytxplot20_MUST_Fungi_Species) <- mytxplot20_MUST_Fungi_Species$OTU

taxa_names(Top20_MUST_Fungi_Species) <- mytxplot20_MUST_Fungi_Species[taxa_names(Top20_MUST_Fungi_Species),"txplt"]

mycols <- c("burlywood2","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","dodgerblue2","#E31A1C")

pdf(file = "Top20_MUST_Fungi_Species.pdf", width = 14, height = 7)

plot_bar(Top20_MUST_Fungi_Species, fill="OTU", title = "") + facet_grid(cols = vars(vinification)) + geom_col() + scale_fill_manual(values = mycols)

dev.off()