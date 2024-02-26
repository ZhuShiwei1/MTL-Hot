
setwd("../../")
library(limma)
library(dplyr)


load("data/TCGA_data/log.exp.matrix.RData") # 16884  8755

load("data/TCGA_data/TCGA_clinical_data.RData")
cluster.info <- read.csv("result/last_shared_layer_clustering/TCGA/0/y_pred_2.csv", row.names = 1)
load("data/TCGA_data/cancer_samples.RData")
cancer_samples$cancer_type[which(cancer_samples$cancer_type %in% c("COAD", "READ"))] <- "COREAD"

res_path <- "result/embedding_clusters_biology_survival_diff/diff_expression_analysis/"
output_dir <- file.path(res_path)
if (!dir.exists(output_dir))  dir.create(output_dir)


###  
#' 
#' Identification of differentially expressed genes using limma
################################
groups <- model.matrix(~factor(cluster.info[,1])+0) 
colnames(groups) <- c("Cluster_0", "Cluster_1")
df.fit <- lmFit(log.exp.matrix, groups)  # lmFit fits a linear model using weighted least squares for each gene
df.matrix <- makeContrasts(Cluster_1 - Cluster_0 , levels = groups)
fit <- contrasts.fit(df.fit, df.matrix) # Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models
fit <- eBayes(fit, trend=TRUE) #Empirical Bayes smoothing of standard errors ， (shrinks standard errors that are much larger or smaller than those from other genes towards the average standard error) 
Output <- topTable(fit,n = Inf, adjust = "fdr")
Output = subset(Output, select=c("adj.P.Val","P.Value","logFC"))
colnames(Output)=c("FDR","P.Value","logFC")
#Raw p-value (based on t) from test that logFC differs from 0
# adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value
# logFC: log2 fold change of Cluster_1/Cluster_0


logFC = 1
fdr <- 0.05
Output[fdr > Output[,"FDR"]  &  Output[,"logFC"] >= logFC, ncol(Output)+1] = "Up−regulated"
Output[fdr > Output[,"FDR"]  & -logFC >= Output[,"logFC"], ncol(Output)] = "Down−regulated"
Output[Output[,"FDR"] >= fdr | logFC > abs(Output[,"logFC"]) , ncol(Output)] = "Unchanged"
colnames(Output)[ncol(Output)] = "Regulate"

diff_exp_genes_df <- Output
diff_exp_genes_df$logP <- -log10(diff_exp_genes_df$FDR)
diff_exp_genes_df$logP[which(diff_exp_genes_df$logP>=350)] <- 350
save(diff_exp_genes_df, file=paste0(res_path,"diff_exp_genes_df.RData"))
table(diff_exp_genes_df$Regulate)
#Down−regulated      Unchanged   Up−regulated 
#159          15712           1013
table(diff_exp_genes_df$FDR<0.05)
#FALSE  TRUE 
#2582 14302

###  
#' 
#' volcano plots visualize the differentially expressed genes
################################
library(dplyr) 
library(ggrepel) 
library(ggplot2)
plot.data <- diff_exp_genes_df
plot.data$Gene <- rownames(plot.data)

diff_gene <- rbind(subset(plot.data, Regulate %in% c("Up−regulated")),
                   subset(plot.data, Regulate %in% c("Down−regulated")))
text_data <- rbind(plot.data[order(plot.data$logFC)[1:10],], plot.data[order(plot.data$logFC, decreasing=T)[1:10],])

pdf(paste0(res_path, "diff_exp_genes_volcPlot1.pdf"), height=6, width=8)
ggplot(plot.data, aes(logFC, logP)) +
  geom_point(size = 1, aes(color = Regulate)) +
  #  The symbol-marked dots indicate the five genes with the largest negative or positive standardized mean difference.
  #geom_label_repel(data = diff_gene, aes(logFC,logP, label = Gene), size = 2, box.padding=0.2)+
  geom_text_repel(data=text_data, aes(label=Gene), box.padding=0.2, size = 2)+
  geom_hline(aes(yintercept=(-log10(0.05))), colour="#404040", linetype="dashed")+
  xlab(expression("log"[2]*" fold change"))+
  ylab(expression("-log"[10]*" fdr")) + 
  scale_color_manual(values = c("steelblue","grey","red")) +
  theme_minimal()
dev.off()


###  
#' 
#' The differentially expressed gene rank list was calculated for specific cancer types
################################
cancers <- c("HNSC", "LUAD", "CESC", "KIRC", "SKCM", "STAD", "LUSC", "BLCA", "BRCA" ,
             "ESCA","OV" ,"THCA", "COREAD", "UCEC", "LIHC")

cancers_diff_exp_genes_df <- lapply(cancers, function(cancer){
  cancer.idx <- which(cancer_samples$cancer_type==cancer)
  per.exp.matrix <- log.exp.matrix[,cancer.idx]
  groups <- model.matrix(~factor(cluster.info[cancer.idx,1])+0) 
  colnames(groups) <- c("Cluster_0", "Cluster_1")
  df.fit <- lmFit(per.exp.matrix, groups)  # lmFit fits a linear model using weighted least squares for each gene
  df.matrix <- makeContrasts(Cluster_1 - Cluster_0 , levels = groups)
  fit <- contrasts.fit(df.fit, df.matrix) # Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models
  fit <- eBayes(fit, trend=TRUE) #Empirical Bayes smoothing of standard errors ， (shrinks standard errors that are much larger or smaller than those from other genes towards the average standard error) 
  Output <- topTable(fit,n = Inf, adjust = "fdr")
  Output = subset(Output, select=c("adj.P.Val","P.Value","logFC"))
  colnames(Output)=c("FDR","P.Value","logFC")
  Output
})
names(cancers_diff_exp_genes_df) <- cancers
save(cancers_diff_exp_genes_df, file=paste0(res_path,"cancers_diff_exp_genes_df.RData"))
