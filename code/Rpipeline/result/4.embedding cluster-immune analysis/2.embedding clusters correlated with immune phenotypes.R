setwd("../../")

library(ggplot2)
library(dplyr)
library(psych)
library(Rtsne)

load("data/TCGA_data/TCGA_clinical_data.RData")
cluster.info <- read.csv("result/last_shared_layer_clustering/TCGA/0/y_pred_2.csv", row.names = 1)
TCGA_clinical_data$embedding_cluster <- cluster.info[,1]
TCGA_clinical_data$embedding_cluster[which(TCGA_clinical_data$embedding_cluster==0)] <- "Cluster_0"
TCGA_clinical_data$embedding_cluster[which(TCGA_clinical_data$embedding_cluster==1)] <- "Cluster_1"
TCGA_clinical_data$embedding_cluster <- factor(TCGA_clinical_data$embedding_cluster)
#Cluster_0 Cluster_1 
#7032      1723 


load("data/TCGA_data/pan.cancer.phenotype.RData")
TCGA_clinical_data <- cbind(TCGA_clinical_data,pan.cancer.phenotype )
res_path <- "result/embedding_clusters_immune_characterizing/imm_phenotype_diff/"

output_dir <- file.path(res_path)
if (!dir.exists(output_dir))  dir.create(output_dir)

###  
#' 
#' t-SNE visualization
################################
last_layer = read.table("result/model_training/MTL/MTL_Hot/final_embeddings/0.txt", header = F,sep=" ")
data <- last_layer
set.seed(321)
tsne_out = Rtsne(
  data,
  dims = 2,
  pca = F, 
  max_iter = 1000,
  theta = 0.2, 
  perplexity = 30, 
  verbose = F
)
tsne_result = as.data.frame(tsne_out$Y)
colnames(tsne_result) = c("tSNE1","tSNE2")

pdf(paste0(res_path,'t_SNE.pdf'))
ggplot(tsne_result,aes(tSNE1,tSNE2)) +
  geom_point(aes(colour= TCGA_clinical_data$embedding_cluster)) +
  scale_color_manual(values=c("#4575b4","#d53e4f"))+
  ggtitle("Embedding clusters") +
  theme_bw()+ theme(panel.grid=element_blank(), 
                    plot.title = element_text(hjust = 0.5), 
                    legend.position="top")
dev.off()



###  
#' 
#' Differences in phenotypes' scores between these two clusters
################################
# 秩和检验组间差异
source("code/RPipeline/function/WilcoxTest.R")
imm.features <- c("APM.score", "Tcell.score","IFNgamma.score", "PDL1.exp",
                  "Leukocyte.Fraction", "CYT_score",
                  "Stromal.Fraction", "TIL.Regional.Fraction",
                  "Nonsilent.Mutation.Rate", "Neoantigens",
                  "Homologous.Recombination.Defects",
                  "BCR.Shannon", "BCR.Richness", "TCR.Shannon","TCR.Richness")
imm.wilcox.test.res<- do.call(rbind,lapply(imm.features, function(feature){
  score <- TCGA_clinical_data[which(!is.na(TCGA_clinical_data[,feature])), feature]
  embedding_cluster <- TCGA_clinical_data[which(!is.na(TCGA_clinical_data[,feature])), "embedding_cluster"]
  x <- as.numeric(score[which(embedding_cluster=="Cluster_0")])
  y <- as.numeric(score[which(embedding_cluster=="Cluster_1")])
  res <- WilcoxTest(x, y, paired = FALSE,types=c("Cluster_0","Cluster_1"), alternative = "two.sided")
  res
}))
rownames(imm.wilcox.test.res) <- imm.features 
save(imm.wilcox.test.res, file=paste0(res_path, "imm.wilcox.test.res.RData"))
imm.wilcox.test.res$p.value
#<2.220446e-16 <2.220446e-16 <2.220446e-16 <2.220446e-16 <2.220446e-16 <2.220446e-16 <2.220446e-16 <2.220446e-16
#3.593971e-10

pdf(paste0(res_path, "imm_features_compare_boxplot.pdf"),height=4, width=4)
for(feature in imm.features){
  plot.data <- TCGA_clinical_data[,c("embedding_cluster",feature)]
  p <- ggplot(plot.data, aes(embedding_cluster,plot.data[,feature]))+
    geom_violin(aes(fill=embedding_cluster))+  
    scale_fill_manual(values = c("#377eb8", "#e41a1c"))+
    geom_boxplot(width=0.1)+
    theme_classic(base_size = 20)+
    theme(axis.text = element_text(color = 'black'),
          legend.position = 'none')+
    labs(x="", y=feature, title="")
  print(p)
}
dev.off()


pdf(paste0(res_path, "imm_features_T_BCR_TMB_Neoantigens_compare_boxplot.pdf"), height=4, width=4)
for(feature in imm.features[c(9,10,13,15)]){
  plot.data <- TCGA_clinical_data[,c("embedding_cluster",feature)]
  Cluster_0_values <- na.omit(plot.data[which(plot.data$embedding_cluster=="Cluster_0"), feature])
  Cluster_0_upper <-  quantile(Cluster_0_values, 0.75)+ 1.5*(quantile(Cluster_0_values, 0.75) - quantile(Cluster_0_values, 0.25))
  Cluster_0_lower <-  quantile(Cluster_0_values, 0.25)- 1.5*(quantile(Cluster_0_values, 0.75) - quantile(Cluster_0_values, 0.25))
  Cluster_1_values <- na.omit(plot.data[which(plot.data$embedding_cluster=="Cluster_1"), feature])
  Cluster_1_upper <-  quantile(Cluster_1_values, 0.75)+ 1.5*(quantile(Cluster_1_values, 0.75) - quantile(Cluster_1_values, 0.25))
  Cluster_1_lower <-  quantile(Cluster_1_values, 0.25)- 1.5*(quantile(Cluster_1_values, 0.75) - quantile(Cluster_1_values, 0.25))
  upper_cutoff <- max(Cluster_0_upper, Cluster_1_upper)
  lower_cutoff <- min(Cluster_0_lower, Cluster_1_lower)
  plot.data[which(plot.data[,feature]> upper_cutoff),feature] <- NA
  plot.data[which(plot.data[,feature]< lower_cutoff),feature] <- NA
  plot.data <- plot.data[complete.cases(plot.data),]
  p <- ggplot(plot.data, aes(embedding_cluster,plot.data[,feature]))+
    geom_violin(aes(fill=embedding_cluster))+  
    scale_fill_manual(values = c("#377eb8", "#e41a1c"))+
    geom_boxplot(width=0.1, outlier.size=0.2)+
    theme_classic(base_size = 20)+
    theme(axis.text = element_text(color = 'black'),
          legend.position = 'none')+
    labs(x="", y=feature, title="")
  print(p)
}
dev.off()



###  
#' 
#' Differences in lymphocyte fractions measured by H&E between these two clusters
################################
score <- TCGA_clinical_data[which(!is.na(TCGA_clinical_data[,"TIL.Regional.Fraction"])), "TIL.Regional.Fraction"]
embedding_cluster <- TCGA_clinical_data[which(!is.na(TCGA_clinical_data[,"TIL.Regional.Fraction"])), "embedding_cluster"]
x <- as.numeric(score[which(embedding_cluster=="Cluster_0")])
y <- as.numeric(score[which(embedding_cluster=="Cluster_1")])
res <- WilcoxTest(x, y, paired = FALSE,types=c("HOT","COLD"), alternative = "two.sided")
# p.value z.statistic
#1 <2.220446e-16   -20.32815

pdf(paste0(res_path, "TCGA_H&E_TIL_compare_boxplot.pdf"), width=4, height=4)
ggplot(TCGA_clinical_data,aes(embedding_cluster,TCGA_clinical_data[,"TIL.Regional.Fraction"]))+
    geom_boxplot(aes(fill=embedding_cluster))+  
    #geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5)+
    scale_fill_manual(values = c("#377eb8", "#e41a1c"))+
    #geom_boxplot(width=0.1)+
    theme_classic(base_size = 20)+
    theme(axis.text = element_text(color = 'black'),
          legend.position = 'none')+
    labs(x="", y="TIL region percentages(H&E)", title="")
dev.off()

