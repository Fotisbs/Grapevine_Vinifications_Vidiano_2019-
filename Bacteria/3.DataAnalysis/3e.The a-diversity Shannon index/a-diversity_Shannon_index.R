##ALPHADIVERSITY Shannon Index with anova/kruskal
bacteria_vinification_Annotated <- readRDS(".................")

MUST_Bacteria_Vidiano_2019 <- bacteria_vinification_Annotated

MUST_Bacteria_Vidiano_2019

View(data.frame(sample_data(MUST_Bacteria_Vidiano_2019)))

MUST_Bacteria_Vidiano_2019

adiv <- plot_richness(MUST_Bacteria_Vidiano_2019, measures=c("Shannon"))

adiv <- plot_richness(MUST_Bacteria_Vidiano_2019, x="stage2", measures=c("InvSimpson", "Shannon"), color="vinification", shape ="vinification")


## calculate Good's coverage estimate as well
install.packages("entropart")
library(entropart)

good <- MetaCommunity(t(MUST_Bacteria_Vidiano_2019@otu_table))$SampleCoverage.communities

good_tbl <- data.frame(samples = names(good), variable = rep("coverage",length(good)), value=good)

alpha_long <- rbind(adiv$data[,c("samples", "variable", "value")], good_tbl)

# convert long to wide
library(reshape2)
alpha_wide <- dcast(alpha_long, samples ~ variable, value.var="value")

row.names(alpha_wide) <- alpha_wide$samples

alpha_wide <- alpha_wide[,c("InvSimpson", "Shannon")]

colnames(alpha_wide) <- c("InvSimpson", "Shannon")

alpha_wide_fact <- merge(alpha_wide, data.frame(sample_data(MUST_Bacteria_Vidiano_2019)), by = "row.names")

row.names(alpha_wide_fact) <- alpha_wide_fact$Row.names
alpha_wide_fact <- alpha_wide_fact[-which(colnames(alpha_wide_fact)%in%"Row.names")]
alpha_wide_fact$vinification <- factor(alpha_wide_fact$vinification)

library("agricolae")
#### perform anova or equivalent for the alpha diversity indices ----

##### from
mytestvars <- colnames(alpha_wide)[2:1]


library("agricolae")

mystatsout <- list()

mytestfacts <- colnames(alpha_wide_fact[, c(6,9)])
mytestfacts <- colnames(alpha_wide_fact[,6:9])


for(mytestfact in mytestfacts){
  cairo_pdf(paste("alpha_div_plot_",mytestfact,".pdf", sep = ""), height = 3, width = 9, onefile = T)
  par(mfrow = c(2,3))
  for(mytestvar in mytestvars){
    
    
    
    # create the aov matrix
    myaovmatrix <- alpha_wide_fact
    
    # run a shapiro wilk test to select parametric or non parametric analysis
    shap_out <- shapiro.test(myaovmatrix[,mytestvar])
    mystatsout[[paste(mytestvar, sep = " // ")]][["shap"]] <- shap_out
    
    # run the parametric or non-parametric analysis according to the shapiro.test results
    if(shap_out$p.value < 0.05){
      # non-parametric
      mykrusk <- kruskal(myaovmatrix[,mytestvar], myaovmatrix[,mytestfact], group = T, p.adj = "BH")
      
      mystatsout[[paste(mytestvar, sep = " // ")]][["krusk"]] <- mykrusk
      # prepare the barplot
      myaovmatrix[,mytestfact] <- factor(myaovmatrix[,mytestfact])
      mytestvarord <- levels(myaovmatrix[,mytestfact])
      par(mar = c(2,10,4,2))
      barerrplt <- bar.err(mykrusk$means[mytestvarord[length(mytestvarord):1],], variation="SD", xlim=c(0,1.2*max(mykrusk$means[mytestvarord[length(mytestvarord):1],1] + mykrusk$means[mytestvarord[length(mytestvarord):1],3])),horiz=TRUE, bar=T, col="grey60", space=1, main= paste(mytestvar), las=1)
      
      par(xpd = T)
      if(mykrusk$statistics$p.chisq <= 0.05){
        text(x = mykrusk$means[mytestvarord[length(mytestvarord):1],1] + mykrusk$means[mytestvarord[length(mytestvarord):1],3], barerrplt$x,labels = mykrusk$groups[mytestvarord[length(mytestvarord):1],2], pos = 4, font = 3)
        par(xpd = F)
      }
    } else { 
      # perform the parametric
      
      # select the alphadiv matrix and design rows
      myform <- as.formula(paste("`",mytestvar,"` ~ ",mytestfact, sep = ""))
      mymod <- aov(myform, data = myaovmatrix)
      mysumaov <- summary(mymod)
      
      mystatsout[[paste(mytestvar, sep = " // ")]][["ANOVA"]] <- mysumaov
      
      # order the matrices etc
      mytestvarord <- levels(factor(myaovmatrix[,mytestfact]))
      
      # run the Tukey test
      myHSDtest <- HSD.test(mymod, mytestfact, group=T)
      
      mystatsout[[paste(mytestvar, sep = " // ")]][["HSD test"]] <- myHSDtest
      
      # prepare the barplot
      par(mar = c(2,10,4,2))
      barerrplt <- bar.err(myHSDtest$means[mytestvarord[length(mytestvarord):1],], variation="SD", xlim=c(0,1.2*max(myHSDtest$means[mytestvarord[length(mytestvarord):1],1] + myHSDtest$means[mytestvarord[length(mytestvarord):1],2])),horiz=TRUE, bar=T, col="grey60", space=1, main= paste(mytestvar), las=1)
      
      par(xpd = T)
      if(mysumaov[[1]]$`Pr(>F)`[1] <= 0.05){
        text(x = myHSDtest$means[mytestvarord[length(mytestvarord):1],1] + myHSDtest$means[mytestvarord[length(mytestvarord):1],2], barerrplt$x,labels = myHSDtest$groups[mytestvarord[length(mytestvarord):1],2], pos = 4)
        par(xpd = F)
      }
    }
    
    
  }
  dev.off()
}

capture.output(mystatsout,file = paste("alpha_div_stats_",mytestfact,".txt", sep = ""))
##### till here

##Start Extra Figures Alpha Diversity, lines with errors bars##
##Bacteria
MUST_Bacteria_Vidiano_2019.NMS <- MUST_Bacteria_Vidiano_2019
sample_names(MUST_Bacteria_Vidiano_2019.NMS) <- paste("samp",sample_names(MUST_Bacteria_Vidiano_2019), sep = "")

myrichness <- estimate_richness(MUST_Bacteria_Vidiano_2019.NMS, split = TRUE, measures = c("Shannon", "InvSimpson", "Observed")) # calculate the Inverse 

mysamdat <- data.frame(sample_data(MUST_Bacteria_Vidiano_2019.NMS)[,c("vinification","stage","stage2")]) # obtain the sample info data.frame

myrichness_fin <- merge(myrichness, mysamdat, by = "row.names") # merge the two tables
row.names(myrichness_fin) <- myrichness_fin$Row.names # add the row.names after merging
myrichness_fin <- myrichness_fin[,-which(colnames(myrichness_fin)%in%"Row.names")] # remove the remnant row names column

# prepare the plot Shannon
plt <- ggplot(myrichness_fin, aes(stage2, Shannon, group = vinification, colour = vinification)) + facet_wrap(~vinification, scales = "free_x", nrow = 2) + ylim(0,3)
plt <- plt + stat_summary(fun = "mean", geom="line", size = 1, position=position_dodge(0))  
plt <- plt + stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) + theme_light()
# prepare the plot Observed
plt <- ggplot(myrichness_fin, aes(stage2, Observed, group = vinification, colour = vinification)) + facet_wrap(~vinification, scales = "free_x", nrow = 2) + ylim(0,25)
plt <- plt + stat_summary(fun = "mean", geom="line", size = 1, position=position_dodge(0))  
plt <- plt + stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) + theme_light()
# prepare the plot InvSimpson
plt <- ggplot(myrichness_fin, aes(stage2, InvSimpson, group = vinification, colour = vinification)) + facet_wrap(~vinification, scales = "free_x", nrow = 2) + ylim(0,10)
plt <- plt + stat_summary(fun = "mean", geom="line", size = 1, position=position_dodge(0))  
plt <- plt + stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) + theme_light()