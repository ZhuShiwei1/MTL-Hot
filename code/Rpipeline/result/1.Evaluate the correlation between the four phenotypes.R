########################################################
#' evaluated the correlations between the four phenotypes  
#' 
##################################################
setwd("../../")

load("data/TCGA_data/pan.cancer.phenotype.RData")

###  
#' spearman correlation 
#' 
################################
library(Hmisc)
cor.phenotypes <- rcorr(pan.cancer.phenotype, type="spearman")
# APM.score Tcell.score IFNgamma.score PDL1.exp
# APM.score           1.00        0.34           0.74     0.53
# Tcell.score         0.34        1.00           0.61     0.49
# IFNgamma.score      0.74        0.61           1.00     0.60
# PDL1.exp            0.53        0.49           0.60     1.00
# All associations were significant



###  
#' Visualize the correlation between the four phenotypes
#' 
################################
library(corrplot)
pdf("result/data_prepare/phenotypes_corr.pdf", width = 5, height = 5)
corrplot(cor.phenotypes$r,type="upper",tl.col ="black",
         tl.srt = 45,addCoef.col = "black", tl.cex=0.8)
dev.off()


###  
#' 
#' Heat map visualization of the phenotypic score after normalization
################################

library(pheatmap)
#breaks
bk <- c(seq(-4,-0.1,by=0.01),seq(0,4,by=0.01))
plot.data <- t(pan.cancer.phenotype)
pheatmap(plot.data, scale="row",
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         breaks=bk, cluster_rows=F,
         treeheight_row=0,show_rownames=T,show_colnames=F,
         width=8, height=4,
         filename="result/data_prepare/phenotypes_heatmap1.pdf")


###  
#' 
#' The boxplots visualize the distribution of scores in histological immune subtypes
################################
# Braun.NatMed.2021
load("data/validation_data/Braun.NatMed.2021/processed_data/cancer.phenotype.RData")
clincial_data <- read.csv("data/validation_data/Braun.NatMed.2021/processed_data/sub_clinical_data.csv", header = T)
df <- as.data.frame(cancer.phenotype)
df$ImmunoPhenotype <- clincial_data$ImmunoPhenotype

df <- df[-which(df$ImmunoPhenotype=="No_IF"),]  # 103   6
pdf('result/data_prepare/Braun_boxplot.pdf', width=10, height=4)
grid.newpage()  
pushViewport(viewport(layout = grid.layout(1,4)))
for(i in 1:4){
  vp1 <- viewport(layout.pos.row = 1, layout.pos.col = i, 
                  width = unit(0,"npc"),
                  height = unit(0,"npc"))
  e <- ggplot(df, aes(x = factor(ImmunoPhenotype, levels=c("Infiltrated", "Excluded", "Desert")), 
                      y = df[,i], fill = ImmunoPhenotype))+  
    geom_boxplot()
  p <- e + geom_dotplot(binaxis = "y", stackdir = "center",position = position_dodge(1),dotsize = 0.5)+
    geom_signif(comparisons = list(c("Infiltrated", "Excluded"),
                                   c("Infiltrated", "Desert")),
                test=wilcox.test, step_increase=0.1)+
    scale_fill_manual(values= c("Infiltrated"="#d53e4f", "Excluded"= "#66c2a5", "Desert"="#4575b4"))+
    labs(x = "", y = colnames(df)[i])+
    ggtitle(colnames(df)[i])+
    guides(fill="none")+
    theme_bw()
  print(p, vp=vp1)
}
dev.off()


# Mariathasan.Nature.2018
load("data/validation_data/Mariathasan.Nature.2018/processed_data/cancer.phenotype.RData")
clincial_data <- read.csv("data/validation_data/Mariathasan.Nature.2018/processed_data/sub_clinical_data.csv", header = T)

df <- as.data.frame(cancer.phenotype)
df$ImmunoPhenotype <- clincial_data$Immune.phenotype
df <- df[-which(is.na(df$ImmunoPhenotype)),]  # 103   6
pdf('result/data_prepare/Mariathasan_boxplot.pdf', width=10, height=4)
grid.newpage()  
pushViewport(viewport(layout = grid.layout(1,4)))
for(i in 1:4){
  vp1 <- viewport(layout.pos.row = 1, layout.pos.col = i, 
                  width = unit(0,"npc"),
                  height = unit(0,"npc"))
  e <- ggplot(df, aes(x = factor(ImmunoPhenotype, levels=c("inflamed", "excluded", "desert")), 
                      y = df[,i], fill = ImmunoPhenotype))+  
    geom_boxplot()
  p <- e + geom_dotplot(binaxis = "y", stackdir = "center",position = position_dodge(1),dotsize = 0.5)+
    geom_signif(comparisons = list(c("inflamed", "excluded"),
                                   c("inflamed", "desert")), 
                test=wilcox.test, step_increase=0.1)+
    scale_fill_manual(values= c("inflamed"="#d53e4f", "excluded"= "#66c2a5", "desert"="#4575b4"))+
    labs(x = "", y = colnames(df)[i])+
    ggtitle(colnames(df)[i])+
    guides(fill="none")+
    theme_bw()
  print(p, vp=vp1)
}
dev.off()

