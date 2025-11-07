##Rarefaction Curves##
fungi_vinification_Annotated <- readRDS(".................")

MUST_Fungi_Vidiano_2019 <- fungi_vinification_Annotated

a <- data.frame(t(otu_table(MUST_Fungi_Vidiano_2019)))

rarecurve(a, step=50, cex=0.5, label = F, xlim=c(0,400), ylim =c(0,20))