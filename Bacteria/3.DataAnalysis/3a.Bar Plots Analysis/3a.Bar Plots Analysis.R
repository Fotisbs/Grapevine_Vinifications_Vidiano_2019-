##Start the analysis by Merging the replicates of each sample per stage of Spontaneous and Inoculated Vinifications##
MUST_Bacteria_Vidiano_2019_Merged <- merge_samples(bacteria_vinification_Annotated, "stage")

MUST_Bacteria_Vidiano_2019_Merged <- prune_taxa(taxa_sums(MUST_Bacteria_Vidiano_2019_Merged)>0,MUST_Bacteria_Vidiano_2019_Merged)

View(data.frame(sample_data(MUST_Bacteria_Vidiano_2019_Merged)))
View(data.frame(otu_table(MUST_Bacteria_Vidiano_2019_Merged)))
View(data.frame(tax_table(MUST_Bacteria_Vidiano_2019_Merged)))

##Rename sample data file after merging
write.table(data.frame(sample_data(MUST_Bacteria_Vidiano_2019_Merged)), file="MUST_Bacteria_Vidiano_2019_Merged.txt", quote = F,col.names = NA, sep="\t")

SampleDataNew67 <- read.table("MUST_Bacteria_Vidiano_2019_Merged.txt", header=T,sep = "\t",row.names = 1)

sample_data(MUST_Bacteria_Vidiano_2019_Merged) <- SampleDataNew67

MUST_Bacteria_Vidiano_2019_Merged

View(data.frame(sample_data(MUST_Bacteria_Vidiano_2019_Merged)))

##Transform RA 100%
MUST_Bacteria_Vidiano_2019_Merged #(Raw Data Merged)

MUST_Bacteria_Vidiano_2019_Merged_100 <- transform_sample_counts(MUST_Bacteria_Vidiano_2019_Merged, function(OTU) 100*OTU/sum(OTU))

MUST_Bacteria_Vidiano_2019_Merged_100 #(Transformation 100%)

#Most abundant taxa (Top 20 at genus level) 
MUST_Bacteria_Vidiano_2019_Merged
MUST_Bacteria_Vidiano_2019_Merged_100

##tax glom at genus level first#
rank_names(MUST_Bacteria_Vidiano_2019_Merged_100)

MUST_Bacteria_Vidiano_2019_Merged_100

MUST_Bacteria_Genus <- tax_glom(MUST_Bacteria_Vidiano_2019_Merged_100, taxrank = "Genus")

MUST_Bacteria_Genus

MUST_Bacteria_Genus <- prune_taxa(taxa_sums(MUST_Bacteria_Genus)>0,MUST_Bacteria_Genus)

##top20
myTaxa20_MUST_Bacteria_Genus <- names(sort(taxa_sums(MUST_Bacteria_Genus), decreasing = TRUE)[1:20])  

Top20_MUST_Bacteria_Genus <- prune_taxa(myTaxa20_MUST_Bacteria_Genus, MUST_Bacteria_Genus)

taxa_names(Top20_MUST_Bacteria_Genus)

Top20_MUST_Bacteria_Genus <- prune_taxa(taxa_sums(Top20_MUST_Bacteria_Genus)>0,Top20_MUST_Bacteria_Genus)

mytax20_MUST_Bacteria_Genus <- data.frame(tax_table(Top20_MUST_Bacteria_Genus), stringsAsFactors = F)

mytxplot20_MUST_Bacteria_Genus <- data.frame(OTU = row.names(mytax20_MUST_Bacteria_Genus), 
                                             txplt = paste(row.names(mytax20_MUST_Bacteria_Genus), " ", mytax20_MUST_Bacteria_Genus$Genus,  ":", mytax20_MUST_Bacteria_Genus$Genus,  sep = ""))

row.names(mytxplot20_MUST_Bacteria_Genus) <- mytxplot20_MUST_Bacteria_Genus$OTU

taxa_names(Top20_MUST_Bacteria_Genus) <- mytxplot20_MUST_Bacteria_Genus[taxa_names(Top20_MUST_Bacteria_Genus),"txplt"]

mycols <- c("dodgerblue2","#E31A1C","green4","#6A3D9A","#FF7F00","black", "gold1","skyblue2","#FB9A99","palegreen2","#CAB2D6","#FDBF6F","gray70", "khaki2","maroon","orchid1","deeppink1","blue1","steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")

pdf(file = "Top20_MUST_Bacteria_Genus.pdf", width = 14, height = 7)

plot_bar(Top20_MUST_Bacteria_Genus, fill="OTU", title = "") + facet_grid(cols = vars(vinification)) + geom_col() + scale_fill_manual(values = mycols)


dev.off()
