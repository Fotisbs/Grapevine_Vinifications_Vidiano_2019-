##Rarefaction Curves##
bacteria_vinification_Annotated <- readRDS(".................")

MUST_Bacteria_Vidiano_2019 <- bacteria_vinification_Annotated

a <- data.frame(t(otu_table(MUST_Bacteria_Vidiano_2019)))


rarecurve(a, step=50, cex=0.5, label = F, xlim=c(0,650), ylim =c(0,20))
