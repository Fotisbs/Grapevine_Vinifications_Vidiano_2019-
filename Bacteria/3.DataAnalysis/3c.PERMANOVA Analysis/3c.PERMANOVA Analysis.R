##PERMANOVA##
##PERMANOVA Variables (perm. 999)##
bacteria_vinification_Annotated <- readRDS(".................")

MUST_Bacteria_Vidiano_2019 <- Bacteria_vinification_Annotated

##Transform RA 100%
MUST_Bacteria_Vidiano_2019_100 <- transform_sample_counts(MUST_Bacteria_Vidiano_2019, function(OTU) 100*OTU/sum(OTU))

MUST_Bacteria_Vidiano_2019_100 #(Transformation 100%)

mypermanova_MUST_Bacteria_Vidiano_2019_100 <- adonis2(MUST_Bacteria_Vidiano_2019_100@otu_table ~ vinification + stage2, method = "bray", data = data.frame(MUST_Bacteria_Vidiano_2019_100@sam_data))

write.table(data.frame(mypermanova_MUST_Bacteria_Vidiano_2019_100), file="mypermanova_MUST_Bacteria_Vidiano_2019_100.txt", quote = F,col.names = NA, sep="\t")
