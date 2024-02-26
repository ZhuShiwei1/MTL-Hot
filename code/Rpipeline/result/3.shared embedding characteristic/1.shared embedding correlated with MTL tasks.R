
setwd("../../")
library(Rtsne)
library(ggplot2)
library(scales)
library(plot3D)
library(RColorBrewer)
res_path <- "result/last_shared_layer_related_to_immune_phenotypes/"
output_dir <- file.path(res_path)
if (!dir.exists(output_dir))  dir.create(output_dir)

###  
#' data preparation

#' 
################################
last_layer = read.table("result/model_training/MTL/MTL_Hot/final_embeddings/0.txt", header = F,sep=" ")

load("data/TCGA_data/pan.cancer.phenotype.RData")


load("data/TCGA_data/TCGA_clinical_data.RData")
##########################visualization#####################################################
# The phenotypic label values were converted into percentiles

pan.cancer.phenotype_percentage <- apply(pan.cancer.phenotype, 2, function(x){
  y = rep(0,length(x))
  for(i in 1:length(x)) {
    y[i] = round(length(which(x<=x[i]))/length(x)*100,4)}
  y
})
mean_percentage <- apply(pan.cancer.phenotype_percentage, 1, mean)
pan.cancer.phenotype_percentage <- cbind(pan.cancer.phenotype_percentage, mean_percentage)

###  
#' t-SNE visualization
#' 
################################
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


#######
# 
#' Visualize the association of the last shared layer with the four phenotypes on the model
#######
pdf(paste0(res_path,'MTL_t_SNE.pdf'))
for( i in 1:5){
  Phenotype_percentile <- pan.cancer.phenotype_percentage[,i]
  p=ggplot(tsne_result,aes(tSNE1,tSNE2)) +
    geom_point(aes(colour= Phenotype_percentile)) +
    scale_color_gradientn(colors=c("#08306b", "#2171b5","white", "#d7301f", "#7f0000"))+ 
    ggtitle(colnames(pan.cancer.phenotype_percentage)[i]) +
    theme_bw()+ theme(panel.grid=element_blank(), 
                      plot.title = element_text(hjust = 0.5), 
                      legend.position="top")
  print(p)
}
dev.off()

#######
# 
#' Visualize the association of the last shared layer with other immune-related phenotypes
#######
Immune_features <- TCGA_clinical_data[, c("Nonsilent.Mutation.Rate", "Neoantigens",
                                          "Leukocyte.Fraction","CYT_score")]
Immune_features_percentage <- apply(Immune_features, 2, function(x){
  y = rep(0,length(x))
  for(i in 1:length(x)){
    if(is.na(x[i])){
      y[i] <- NA
    }else{
      y[i] = round(length(which(x[which(!is.na(x))]<=x[i]))/length(which(!is.na(x)))*100,4)}
  }
  y
})
pdf(paste0(res_path,'/MTL_Immune_features_t_SNE1.pdf'))
for( i in 1:6){
  Phenotype_percentile <- Immune_features_percentage[,i]
  p=ggplot(tsne_result,aes(tSNE1,tSNE2)) +
    geom_point(aes(colour= Phenotype_percentile)) +
    scale_color_gradientn(colors=c("#4575b4", "#91bfdb","white", "#fc8d59", "#d73027"))+ 
    ggtitle(colnames(Immune_features_percentage)[i]) +
    theme_bw()+ theme(panel.grid=element_blank(), 
                      plot.title = element_text(hjust = 0.5), 
                      legend.position="top")
  print(p)
}
dev.off()

#######
# 
#' Visualize the last shared layer's association with an immune subtyping defined by others
#######
pdf(paste0(res_path,'/MTL_MFP_t_SNE.pdf'))
MFP=factor(TCGA_clinical_data$MFP, levels=c("IE", "IE/F", "F", "D"))
p=ggplot(tsne_result,aes(tSNE1,tSNE2)) +
  geom_point(aes(colour= MFP)) +
  scale_color_manual(values=c("#d53e4f", "#fdae61", "#66c2a5", "#4575b4"))+ 
  ggtitle("MFP") +
  theme_bw()+ theme(panel.grid=element_blank(), 
                    plot.title = element_text(hjust = 0.5), 
                    legend.position="top")
print(p)
dev.off()

pdf(paste0(res_path,'/MTL_Immune_Subtype_t_SNE.pdf'))
TCGA_Immune_Subtype=factor(TCGA_clinical_data$Immune.Subtype, 
                           levels=c("C1", "C2", "C3", "C4", "C5", "C6"))
p=ggplot(tsne_result,aes(tSNE1,tSNE2)) +
  geom_point(aes(colour= TCGA_Immune_Subtype)) +
  scale_color_manual(values=c("#E83A17", "#F0EA39", "#60B52F", "#6CC6D6", "#285DAB", "#AF589E"))+ 
  ggtitle("TCGA_Immune_Subtype") +
  theme_bw()+ theme(panel.grid=element_blank(), 
                    plot.title = element_text(hjust = 0.5), 
                    legend.position="top")
print(p)
dev.off()


#######
# 
#' Visualize the association of the last shared layer with the cancer type
#######
library(ggsci)
pdf(paste0(res_path,'/MTL_cancer_type_t_SNE.pdf'))
p=ggplot(tsne_result,aes(tSNE1,tSNE2)) +
    geom_point(aes(colour= TCGA_clinical_data$cancer_type)) +
    scale_color_manual(values=pal_ucscgb("default", alpha = 0.6)(25))+ 
    ggtitle("Cancer type") +
    theme_bw()+ theme(panel.grid=element_blank(), 
                      plot.title = element_text(hjust = 0.5), 
                      legend.position="top")
print(p)
dev.off()



###  
#' 
#' The significant associations of the last shared layer with tumor immune-related features were calculated
################################
data <- last_layer

clinical_factors <- c("Leukocyte.Fraction","CYT_score", 
                      "Nonsilent.Mutation.Rate", 
                      "Neoantigens")
clinical_info <- TCGA_clinical_data[, clinical_factors]
rownames(clinical_info) <- TCGA_clinical_data$samples


spearman_corr_res <- lapply(1:ncol(last_layer), function(i){

  inflamed_corr <- apply(pan.cancer.phenotype, 2, function(score){
    res=cor.test(data[which(!is.na(score)),i], score[which(!is.na(score))], method="spearman",exact=FALSE)
    c(res$p.value, res$estimate)
  })

  clinical_corr <- apply(clinical_info, 2, function(score){
    res=cor.test(data[which(!is.na(score)),i], score[which(!is.na(score))], method="spearman",exact=FALSE)
    c(res$p.value, res$estimate)
  })
  
  data.frame(pvalues=c(inflamed_corr[1,], clinical_corr[1,]),
             corr = c(inflamed_corr[2,], clinical_corr[2,]))
})
spearman_corr_pvalue = do.call(rbind, lapply(spearman_corr_res, function(x)x$pvalues))
colnames(spearman_corr_pvalue) <- c("APM.score", "Tcell.score", "IFNgamma.score",  "PDL1.exp",
                                    "Leukocyte.Fraction","CYT_score", 
                                    "TMB", "Neoantigens")
spearman_corr_rho = do.call(rbind, lapply(spearman_corr_res, function(x)x$corr))
colnames(spearman_corr_rho) <- c("APM.score", "Tcell.score", "IFNgamma.score",  "PDL1.exp",
                                 "Leukocyte.Fraction","CYT_score", 
                                 "TMB", "Neoantigens")


spearman_corr_adj_pvalue <- t(apply(spearman_corr_pvalue, 1, function(x)p.adjust(x, "BH")))
apply(spearman_corr_adj_pvalue, 2, function(x)length(which(x<0.05)))
# APM.score         Tcell.score      IFNgamma.score            PDL1.exp  Leukocyte.Fraction 
# 98                  97                  97                  99                  99 
# CYT_score                   TMB         Neoantigens 
# 97                            94                  93 
save(spearman_corr_rho, spearman_corr_pvalue, spearman_corr_adj_pvalue, 
     file=paste0(res_path,"/MTL_nodes_output_corr_imm.RData" ))



library(pheatmap)
plot.data <- t(spearman_corr_rho)

mark_matrix <- matrix(t(spearman_corr_adj_pvalue) <= 0.05, 8)
mark_matrix <- ifelse(t(spearman_corr_adj_pvalue)> 0.05, "ns", "")
  
bk <- c(seq(-1,-0.001,by=0.001),seq(0,1,by=0.001))
pheatmap(plot.data,cluster_row = FALSE,cluster_cols = TRUE, 
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         display_numbers = mark_matrix,
         breaks=bk,
         show_colnames = FALSE,scale="none",treeheight_col=0,
         filename = paste0(res_path,'MTL_imm_corr_heatmap1.pdf'),
         width=15,
         height=2.5
         ) 
