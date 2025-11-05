##Transform RA 100%
MUST_Fungi_Vidiano_2019_Merged #(Raw Data Merged)

MUST_Fungi_Vidiano_2019_Merged_100 <- transform_sample_counts(MUST_Fungi_Vidiano_2019_Merged, function(OTU) 100*OTU/sum(OTU))

MUST_Fungi_Vidiano_2019_Merged_100 #(Transformation 100%)

##NMDS##
##ordination matrix
ord.nmds.bray1 <- ordinate(MUST_Fungi_Vidiano_2019_Merged_100, method="NMDS", distance="bray")

#Vinification
plot_ordination(MUST_Fungi_Vidiano_2019_Merged_100, ord.nmds.bray1, color="vinification", shape ="", label = "", title=paste("NMDS (stress ",round(ord.nmds.bray1$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("#984EA3","#FF7F00","#33A02C"))

##Stage
plot_ordination(MUST_Fungi_Vidiano_2019_Merged_100, ord.nmds.bray1, color="stage", shape ="", label = "", title=paste("NMDS (stress ",round(ord.nmds.bray1$stress, 2),")", sep = "")) + geom_point(size = 4) + scale_color_manual(values=c("#984EA3","#FF7F00","#33A02C","#80b1d3","Blue"))

##Manuscript Format Alernative Style##
##Vinification##
MUST_Fungi_Vidiano_2019_Merged_100

ord.nmds.bray1 <- ordinate(MUST_Fungi_Vidiano_2019_Merged_100, method="NMDS", distance="bray")

mytax_tbl <- data.frame(tax_table(MUST_Fungi_Vidiano_2019_Merged_100))

mytax_tbl$forplt <- paste(mytax_tbl$Genus, mytax_tbl$Species, sep = ":")

myplotcols <- RColorBrewer::brewer.pal(n = 4, name = 'RdBu')

mynmdssit <- ord.nmds.bray1$points

mynmdsspe <- ord.nmds.bray1$species

my_sel_var <- factor(MUST_Fungi_Vidiano_2019_Merged_100@sam_data$vinification)

myterr_sel <- as.numeric(my_sel_var)

cairo_pdf("Fungi_Vinification_NMDS.pdf", height = 6, width = 7)

plot(mynmdssit, frame = F, cex = 0, pch = 21, xlim = c(min(1.3*mynmdssit[,1]),max(1.3*mynmdssit[,1])))

vegan::ordiellipse(mynmdssit, groups = my_sel_var, kind = "ehull", lty = 2, lwd=1)

points(mynmdssit, bg = myplotcols[myterr_sel], pch = 21, cex = 1.5)

par(adj = 0)
title(sub = paste("stress ", round(ord.nmds.bray1$stress,2), sep = ""), cex.sub = 1.2)
par(adj = 1)
title(sub = "")
par(adj = .5)

library(TeachingDemos)

graphics::legend("topright",bty = "n", legend = levels(my_sel_var), pch = 21, pt.bg = myplotcols[1:length(levels(my_sel_var))], pt.cex = 1.5)

dev.off()

##Stages##
ord.nmds.bray1 <- ordinate(MUST_Fungi_Vidiano_2019_Merged_100, method="NMDS", distance="bray")

mytax_tbl <- data.frame(tax_table(MUST_Fungi_Vidiano_2019_Merged_100))

mytax_tbl$forplt <- paste(mytax_tbl$Genus, mytax_tbl$Species, sep = ":")

myplotcols <- RColorBrewer::brewer.pal(n = 6, name = 'RdBu')

mynmdssit <- ord.nmds.bray1$points

mynmdsspe <- ord.nmds.bray1$species

my_sel_var <- factor(MUST_Fungi_Vidiano_2019_Merged_100@sam_data$stage)

myterr_sel <- as.numeric(my_sel_var)

cairo_pdf("Fungi_Stages_NMDS.pdf", height = 6, width = 7)

plot(mynmdssit, frame = F, cex = 0, pch = 21, xlim = c(min(1.3*mynmdssit[,1]),max(1.3*mynmdssit[,1])))

vegan::ordiellipse(mynmdssit, groups = my_sel_var, kind = "ehull", lty = 2, lwd=1)

points(mynmdssit, bg = myplotcols[myterr_sel], pch = 21, cex = 1.5)

par(adj = 0)
title(sub = paste("stress ", round(ord.nmds.bray1$stress,2), sep = ""), cex.sub = 1.2)
par(adj = 1)
title(sub = "")
par(adj = .5)

library(TeachingDemos)

graphics::legend("topright",bty = "n", legend = levels(my_sel_var), pch = 21, pt.bg = myplotcols[1:length(levels(my_sel_var))], pt.cex = 1.5)

dev.off()