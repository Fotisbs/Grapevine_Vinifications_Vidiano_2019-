##Differential Abundance Heatmaps##
##Vinifications
bacteria_vinification_Annotated <- readRDS(".................")

MUST_Bacteria_Vidiano_2019 <- bacteria_vinification_Annotated

MUST_Bacteria_Genus_2 <- tax_glom(MUST_Bacteria_Vidiano_2019, taxrank = "Genus")

MUST_Bacteria_Genus_2

psdt <- MUST_Bacteria_Genus_2

##transform phyloseq object raw counts to relative abundance##
psdt_ra <- transform_sample_counts(psdt, function(x) x / sum(x))

##Define the number of OTUs to be used in the analysis, 1060 ASVs in our case   
mynumOTUs <- 1060

## Prepare working object and select experimental variable of interest##
ps_htmp <- prune_taxa(names(taxa_sums(psdt_ra)[order(taxa_sums(psdt_ra), decreasing = T)][1:mynumOTUs]),psdt_ra)

## Prepare experimental variables table##
mydesign <- data.frame(sample_data(ps_htmp))

## Prepare the relative abundance table for the top ASVs
dtcounts <- decostand(data.frame(otu_table(ps_htmp))[row.names(mydesign),], method = "range", MARGIN = 2)

## ASVs 
myotus <- colnames(dtcounts)

## Run Kruskal-Wallis, non-parametrix multiple comparison test##
mystatsout <- list()

for(myotu in myotus){
  if(sum(dtcounts[,myotu]) == 0){
    mystatsout[[myotu]]$krusk$statistics$p.chisq <- 1
  } else {
    mykrusk <- kruskal(dtcounts[,myotu], mydesign[rownames(dtcounts),]$vinification, group = T)
    mystatsout[[myotu]][["krusk"]] <- mykrusk
  }
}

## Extract p-values from the test results without p.adj##
mystatsoutkruskpvals <- unlist(lapply(myotus, function(x) mystatsout[[x]]$krusk$statistics$p.chisq))
names(mystatsoutkruskpvals) <- myotus

## Adjust p-values with "False Discovery Rate" method
mystatsoutkruskpvals.adj <- p.adjust(mystatsoutkruskpvals, method = "fdr")

# Check if there are any ASVs that significantly differ with p.adjust
for (i in c(1:length(mystatsoutkruskpvals.adj))) {
  if (mystatsoutkruskpvals.adj[i] < 0.05) {
    print(names(mystatsoutkruskpvals.adj[i]))
  }
}
###without p.adjust
for (i in c(1:length(mystatsoutkruskpvals))) {
  if (mystatsoutkruskpvals[i] < 0.05) {
    print(names(mystatsoutkruskpvals[i]))
  }
}

##p.adjust was preferred to continue if was applicable to the data## 
##Tables preparation of the most abundant micororganisms## 
# prepare the ASV names
mytax <- data.frame(tax_table(ps_htmp), stringsAsFactors = F)

# For ITS - Remove letter from taxonomy
for (i in c(1:nrow(mytax))) {
  for(j in c(1:ncol(mytax))) {
    mytax[i,j] <- gsub("[a-z]__","",mytax[i,j])
  }
}


mytxplot <- data.frame(OTU = row.names(mytax), 
                       txplt = paste(row.names(mytax), " ", mytax$Genus,  " ", mytax$Genus,  sep = ""))

##Prepare stars, select among (mystatsoutkruskpvals,mystatsoutkruskpvals.adj)##
##without p-adjustment)##
mystarspheat <- array()
for(myotu_sig in mytxplot$OTU){
  mystarspheat[myotu_sig] <- stars.pval(get(paste("mystatsoutkruskpvals", sep = ""))[myotu_sig])
}

## OR with p-adjustment##
mystarspheat <- array()
for(myotu_sig in mytxplot$OTU){
  mystarspheat[myotu_sig] <- stars.pval(get(paste("mystatsoutkruskpvals.adj", sep = ""))[myotu_sig])
}

## Remove samples with missing values (NA)
mystarspheat1 <- mystarspheat[complete.cases(mystarspheat)]
mystarspheat1_tbl <- data.frame(star=mystarspheat1, row.names = names(mystarspheat1))
### Prepare final data and names and stars table
##transpose data##
dtcounts <- t(dtcounts) ##antimetathesi##


counts_nms <- merge(mytxplot,dtcounts, by.x = "OTU", by.y = "row.names", all = T)


counts_nms_stars <- merge(mystarspheat1_tbl,counts_nms, by.x = "row.names",by.y = "OTU", all = T)


counts_nms_stars$nms_star <- paste(counts_nms_stars$txplt,counts_nms_stars$star)


counts_fin <- counts_nms_stars[,4:(ncol(counts_nms_stars)-1)]
row.names(counts_fin) <- counts_nms_stars$nms_star

# Prepeare the mean RA % for each taxon in all samples##
myrelabund <- data.frame("mean RA %" = 100*rowMeans(counts_fin), check.names = F)

# Merge the taxon abundance per sample (counts_fin) with the mean RA
counts_fin_relabund <- merge(counts_fin,myrelabund, by = "row.names", all.x = T)
# Remove not-needed columns
counts_fin_relabund1 <- counts_fin_relabund[,-c(1,ncol(counts_fin_relabund))]

rownames(counts_fin_relabund1) <- counts_fin_relabund$Row.names

##Aggregate table according to Vinifications##
##Aggregated table
counts_fin_t <- t(counts_fin)

counts_fin_agg <- aggregate(counts_fin_t, by = list(mydesign$vinification), mean)
rownames(counts_fin_agg) <- counts_fin_agg$Group.1


counts_fin_agg_t <- t(counts_fin_agg[,-grep("Group.1",colnames(counts_fin_agg))])


counts_fin_agg_relabund <- merge(counts_fin_agg_t,myrelabund, by = "row.names", all.x = T)

# Prepare the table with only the treatments per column
counts_fin_agg_relabund1 <- counts_fin_agg_relabund[,grep(paste(levels(mydesign$vinification),collapse = "|"),colnames(counts_fin_agg_relabund))]
# Remove not-needed columns
counts_fin_agg_relabund1 <- counts_fin_agg_relabund[,-c(1,ncol(counts_fin_agg_relabund))]

row.names(counts_fin_agg_relabund1) <- paste(counts_fin_agg_relabund$Row.names)

##Plot Aggregated table##
mat_fr_plot_agg <- counts_fin_agg_relabund1[grep("\\*|\\.",row.names(counts_fin_agg_relabund1)),]

##you can write the table##
write.table(mat_fr_plot_agg, file = "Heatmap_Bacteria_Vidiano_Vinifications.txt")

##Prepare function variables####Aggregated Table##
htmp_mat <- mat_fr_plot_agg  

mat_fr_plot_agg_top20 <- mat_fr_plot_agg[order(rowSums(mat_fr_plot_agg), decreasing = T),][1:20,]

htmp_mat <- mat_fr_plot_agg_top20

htmp_mat <- decostand(mat_fr_plot_agg_top20, "range", MARGIN = 1)

###HeatmapPlot##
htmp_plot <- function(htmp_mat, mycolors, mycuttreecols, mycuttreerows, mean_RA) { 
  pheatmap(htmp_mat, color = mycolors, cluster_rows = T,
           fontsize = 15,
           cutree_cols = mycuttreecols, cutree_rows = mycuttreerows,
           annotation_row = mean_RA, 
           border_color = "grey30",
           cellwidth = 15, cellheight = 15)
}
## Heatmap colors##
mycolors <- colorRampPalette(c("grey95","lightsteelblue4","darkred","darkorange2"))(n = 119)
#mycuttreecols 
mycutree_col <- cutree(hclust(vegan::vegdist(htmp_mat)), h = 0.8)
mycuttreecols <- max(mycutree_col)
#mycuttreerows
mycutree_row <- cutree(hclust(vegan::vegdist(htmp_mat)), h = 0.8)
mycuttreerows<- max(mycutree_row)
#mean_RA
mean_RA <- myrelabund

cairo_pdf(paste("Heatmap_Bacteria_Vidiano_Vinifications.pdf", sep = ""), width =17, height =35)

htmp_plot(htmp_mat = htmp_mat,
          mycolors = mycolors,
          mycuttreecols = mycuttreecols,
          mycuttreerows = mycuttreerows,
          mean_RA = mean_RA)

dev.off()

##Differential Abundance Heatmaps##
##Stages
psdt <- MUST_Bacteria_Genus_2

##transform phyloseq object raw counts to relative abundance##
psdt_ra <- transform_sample_counts(psdt, function(x) x / sum(x))

##Define the number of OTUs to be used in the analysis, 1060 ASVs in our case## 
mynumOTUs <- 1060

## Prepare working object and select experimental variable of interest##
ps_htmp <- prune_taxa(names(taxa_sums(psdt_ra)[order(taxa_sums(psdt_ra), decreasing = T)][1:mynumOTUs]),psdt_ra)

## Prepare experimental variables table##
mydesign <- data.frame(sample_data(ps_htmp))

## Prepare the relative abundance table for the top 200 ASVs
dtcounts <- decostand(data.frame(otu_table(ps_htmp))[row.names(mydesign),], method = "range", MARGIN = 2)

## ASVs 
myotus <- colnames(dtcounts)

## Run Kruskal-Wallis, non-parametrix multiple comparison test##
mystatsout <- list()

for(myotu in myotus){
  if(sum(dtcounts[,myotu]) == 0){
    mystatsout[[myotu]]$krusk$statistics$p.chisq <- 1
  } else {
    mykrusk <- kruskal(dtcounts[,myotu], mydesign[rownames(dtcounts),]$stage2, group = T)
    mystatsout[[myotu]][["krusk"]] <- mykrusk
  }
}

## Extract p-values from the test results without p.adj##
mystatsoutkruskpvals <- unlist(lapply(myotus, function(x) mystatsout[[x]]$krusk$statistics$p.chisq))
names(mystatsoutkruskpvals) <- myotus

## Adjust p-values with "False Discovery Rate" method
mystatsoutkruskpvals.adj <- p.adjust(mystatsoutkruskpvals, method = "fdr")

# Check if there are any ASVs that significantly differ with p.adjust
for (i in c(1:length(mystatsoutkruskpvals.adj))) {
  if (mystatsoutkruskpvals.adj[i] < 0.05) {
    print(names(mystatsoutkruskpvals.adj[i]))
  }
}
###without p.adjust
for (i in c(1:length(mystatsoutkruskpvals))) {
  if (mystatsoutkruskpvals[i] < 0.05) {
    print(names(mystatsoutkruskpvals[i]))
  }
}
##p.adjust was preferred to continue if was applicable to the data## 
##Tables preparation of the most abundant micororganisms## 
# prepare the ASV names
mytax <- data.frame(tax_table(ps_htmp), stringsAsFactors = F)

# For ITS - Remove letter from taxonomy
for (i in c(1:nrow(mytax))) {
  for(j in c(1:ncol(mytax))) {
    mytax[i,j] <- gsub("[a-z]__","",mytax[i,j])
  }
}


mytxplot <- data.frame(OTU = row.names(mytax), 
                       txplt = paste(row.names(mytax), " ", mytax$Genus,  " ", mytax$Genus,  sep = ""))

##Prepare stars, select among (mystatsoutkruskpvals,mystatsoutkruskpvals.adj)##
##without p-adjustment)##
mystarspheat <- array()
for(myotu_sig in mytxplot$OTU){
  mystarspheat[myotu_sig] <- stars.pval(get(paste("mystatsoutkruskpvals", sep = ""))[myotu_sig])
}

## OR with p-adjustment##
mystarspheat <- array()
for(myotu_sig in mytxplot$OTU){
  mystarspheat[myotu_sig] <- stars.pval(get(paste("mystatsoutkruskpvals.adj", sep = ""))[myotu_sig])
}

## Remove samples with missing values (NA)
mystarspheat1 <- mystarspheat[complete.cases(mystarspheat)]
mystarspheat1_tbl <- data.frame(star=mystarspheat1, row.names = names(mystarspheat1))
### Prepare final data and names and stars table
##transpose data##
dtcounts <- t(dtcounts) ##antimetathesi##


counts_nms <- merge(mytxplot,dtcounts, by.x = "OTU", by.y = "row.names", all = T)


counts_nms_stars <- merge(mystarspheat1_tbl,counts_nms, by.x = "row.names",by.y = "OTU", all = T)


counts_nms_stars$nms_star <- paste(counts_nms_stars$txplt,counts_nms_stars$star)


counts_fin <- counts_nms_stars[,4:(ncol(counts_nms_stars)-1)]
row.names(counts_fin) <- counts_nms_stars$nms_star

# Prepeare the mean RA % for each taxon in all samples##
myrelabund <- data.frame("mean RA %" = 100*rowMeans(counts_fin), check.names = F)

# Merge the taxon abundance per sample (counts_fin) with the mean RA
counts_fin_relabund <- merge(counts_fin,myrelabund, by = "row.names", all.x = T)
# Remove not-needed columns
counts_fin_relabund1 <- counts_fin_relabund[,-c(1,ncol(counts_fin_relabund))]

rownames(counts_fin_relabund1) <- counts_fin_relabund$Row.names

##Aggregate table according to GTDs_condition##
##Aggregated table
counts_fin_t <- t(counts_fin)

counts_fin_agg <- aggregate(counts_fin_t, by = list(mydesign$stage2), mean)
rownames(counts_fin_agg) <- counts_fin_agg$Group.1


counts_fin_agg_t <- t(counts_fin_agg[,-grep("Group.1",colnames(counts_fin_agg))])


counts_fin_agg_relabund <- merge(counts_fin_agg_t,myrelabund, by = "row.names", all.x = T)

# Prepare the table with only the treatments per column
counts_fin_agg_relabund1 <- counts_fin_agg_relabund[,grep(paste(levels(mydesign$stage2),collapse = "|"),colnames(counts_fin_agg_relabund))]
# Remove not-needed columns
counts_fin_agg_relabund1 <- counts_fin_agg_relabund[,-c(1,ncol(counts_fin_agg_relabund))]

row.names(counts_fin_agg_relabund1) <- paste(counts_fin_agg_relabund$Row.names)

##Plot Aggregated table## 
mat_fr_plot_agg <- counts_fin_agg_relabund1[grep("\\*|\\.",row.names(counts_fin_agg_relabund1)),]

##you can write the table##
write.table(mat_fr_plot_agg, file = "Heatmap_Bacteria_Vidiano_Stages.txt")

##Prepare function variables####Aggregated Table##
htmp_mat <- mat_fr_plot_agg  

mat_fr_plot_agg_top20 <- mat_fr_plot_agg[order(rowSums(mat_fr_plot_agg), decreasing = T),][1:20,]

htmp_mat <- mat_fr_plot_agg_top20

htmp_mat <- decostand(mat_fr_plot_agg_top20, "range", MARGIN = 1)


###HeatmapPlot##
htmp_plot <- function(htmp_mat, mycolors, mycuttreecols, mycuttreerows, mean_RA) { 
  pheatmap(htmp_mat, color = mycolors, cluster_rows = T,
           fontsize = 15,
           cutree_cols = mycuttreecols, cutree_rows = mycuttreerows,
           annotation_row = mean_RA, 
           border_color = "grey30",
           cellwidth = 15, cellheight = 15)
}
## Heatmap colors##
mycolors <- colorRampPalette(c("grey95","lightsteelblue4","darkred","darkorange2"))(n = 119)
#mycuttreecols 
mycutree_col <- cutree(hclust(vegan::vegdist(htmp_mat)), h = 0.8)
mycuttreecols <- max(mycutree_col)
#mycuttreerows
mycutree_row <- cutree(hclust(vegan::vegdist(htmp_mat)), h = 0.8)
mycuttreerows<- max(mycutree_row)
#mean_RA
mean_RA <- myrelabund

cairo_pdf(paste("Heatmap_Bacteria_Vidiano_Stages.pdf", sep = ""), width =17, height =35)

htmp_plot(htmp_mat = htmp_mat,
          mycolors = mycolors,
          mycuttreecols = mycuttreecols,
          mycuttreerows = mycuttreerows,
          mean_RA = mean_RA)

dev.off()