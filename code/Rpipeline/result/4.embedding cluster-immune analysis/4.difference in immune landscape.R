
setwd("../../")
res_path <- "result/embedding_clusters_immune_characterizing/immune_landscape_diff/"
output_dir <- file.path(res_path)
if (!dir.exists(output_dir))  dir.create(output_dir)



load("data/TCGA_data/log.exp.matrix.RData")
cluster.info <- read.csv("result/last_shared_layer_clustering/TCGA/0/y_pred_2.csv", row.names = 1)
load("data/TCGA_data/cancer_samples.RData")
cancer_samples$cancer_type[which(cancer_samples$cancer_type %in% c("COAD", "READ"))] <- "COREAD"

# load expression data and clustering result of expO project
load("data/validation_data/ExpO_project/processed_data/exp_data.RData")
val_exp_data <- t(exp_data)  
val.cluster.info <- read.csv("result/last_shared_layer_clustering/ExpO_project/0/y_pred_2.csv", row.names = 1)
#0    1 
#1863  295
val.cancer.samples <- read.csv("data/validation_data/ExpO_project/processed_data/cancer.samples.csv")



###  
#' 
#' The difference of the expression level of the IMs 
################################
# load IM genes
load("data/immune_genes/IMs.info.RData")
IMs.exp <- log.exp.matrix[which(rownames(log.exp.matrix) %in% IMs.info[, "HGNC Symbol"]), ]  # 67 8755
val.IMs.exp <- val_exp_data[, which(colnames(val_exp_data) %in% IMs.info[, "HGNC Symbol"])]  # 2158   72
val.IMs.exp  <- log2(val.IMs.exp+1)
#######
# 
#' Calculate differences between clusters
#######
source("code/RPipeline/function/WilcoxTest.R")
# TCGA
data <- as.data.frame(t(IMs.exp))
data$embedding_cluster <- ifelse(cluster.info[,1]==0, "Cluster_0", "Cluster_1")
IMs.wilcox.test.res <- do.call(rbind,lapply(rownames(IMs.exp), function(feature){
  score <- data[which(!is.na(data[,feature])), feature]
  embedding_cluster <- data[which(!is.na(data[,feature])), "embedding_cluster"]
  x <- as.numeric(score[which(embedding_cluster=="Cluster_0")])
  y <- as.numeric(score[which(embedding_cluster=="Cluster_1")])
  res <- WilcoxTest(x, y, paired = FALSE,types=c("Cluster_0","Cluster_1"), alternative = "two.sided")
  log2FC <- log2(mean(y, rm.na=T)/mean(x, rm.na=T))
  data.frame(z.statistic = res$z.statistic, p.value = res$p.value, log2FC=log2FC)
}))
rownames(IMs.wilcox.test.res) <- rownames(IMs.exp) 
IMs.wilcox.test.res$adj.pvalue <- p.adjust(IMs.wilcox.test.res$p.value, method="BH")
save(IMs.wilcox.test.res, file=paste0(res_path, "IMs.wilcox.test.res.RData"))

# expO
data <- as.data.frame(val.IMs.exp)
data$embedding_cluster <- ifelse(val.cluster.info[,1]==0, "Cluster_0", "Cluster_1")
val.IMs.wilcox.test.res<- do.call(rbind,lapply(colnames(val.IMs.exp), function(feature){
  score <- data[which(!is.na(data[,feature])), feature]
  embedding_cluster <- data[which(!is.na(data[,feature])), "embedding_cluster"]
  x <- as.numeric(score[which(embedding_cluster=="Cluster_0")])
  y <- as.numeric(score[which(embedding_cluster=="Cluster_1")])
  res <- WilcoxTest(x, y, paired = FALSE,types=c("Cluster_0","Cluster_1"), alternative = "two.sided")
  log2FC <- log2(mean(y, rm.na=T)/mean(x, rm.na=T))
  data.frame(z.statistic = res$z.statistic, p.value = res$p.value, log2FC=log2FC)
}))
rownames(val.IMs.wilcox.test.res) <- colnames(val.IMs.exp)
val.IMs.wilcox.test.res$adj.pvalue <- p.adjust(val.IMs.wilcox.test.res$p.value, method="BH")
save(val.IMs.wilcox.test.res, file=paste0(res_path, "val.IMs.wilcox.test.res.RData"))


#######
# 
#' visualization
#######
library(ComplexHeatmap)
# TCGA
Immune_Checkpoint <- c("Stimulatory", "Inhibitory", "MHC_Class_I_or_II")
per.types.IMs <- lapply(Immune_Checkpoint, function(x)IMs.info[IMs.info[["HGNC Symbol"]] %in% rownames(IMs.wilcox.test.res) & IMs.info[["Immune Checkpoint"]] == x,"HGNC Symbol"])
IMs.order <- unlist(per.types.IMs)
IMs.exp <- IMs.exp[IMs.order, ]
plot.data <- t(scale(t(IMs.exp)))
plot.data[which(plot.data>5)] <- 5
plot.data[which(plot.data<(-5))] <- -5
pdf(paste0(res_path, "IMs_exp_heatmap.pdf"), height=7, width=12)
top_annotation = HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = c("#2166ac","#b2182b")),
                                                        labels = c("Cluster_0", "Cluster_1"),
                                                        labels_gp = gpar(col = "black")))
Heatmap(plot.data,
        col = colorRampPalette(c("navy","white","firebrick3"))(100),
        show_row_names = T,show_column_dend = F,show_row_dend = F,
        top_annotation = top_annotation,
        column_split = factor(cluster.info[,1]),
        row_split = c(rep("Stimulatory", 34), rep("Inhibitory", 20), rep("MHC_Class_I_or_II", 13)), 
        column_title = NULL,
        show_column_names = F)
dev.off()

# expO
per.types.IMs <- lapply(Immune_Checkpoint, function(x)IMs.info[IMs.info[["HGNC Symbol"]] %in% rownames(val.IMs.wilcox.test.res) & IMs.info[["Immune Checkpoint"]] == x,"HGNC Symbol"])
IMs.order <- unlist(per.types.IMs)
val.IMs.exp <- val.IMs.exp[,IMs.order ]
plot.data <- t(scale(val.IMs.exp))
plot.data[which(plot.data>5)] <- 5
plot.data[which(plot.data<(-5))] <- -5
pdf(paste0(res_path, "val.IMs_exp_heatmap.pdf"), height=7, width=12)
top_annotation = HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = c("#2166ac","#b2182b")),
                                                        labels = c("Cluster_0", "Cluster_1"),
                                                        labels_gp = gpar(col = "black")))
Heatmap(plot.data,
        col = colorRampPalette(c("navy","white","firebrick3"))(100),
        show_row_names = T,show_column_dend = F,show_row_dend = F,
        top_annotation = top_annotation,
        column_split = factor(val.cluster.info[,1]),
        row_split = c(rep("Stimulatory", 37), rep("Inhibitory", 23), rep("MHC_Class_I_or_II", 12)), 
        column_title = NULL,
        show_column_names = F)
dev.off()

#######
# 
#' the difference of the expression level of the IMs in specific cancer types
#######
cancers <- c("HNSC", "LUAD", "CESC", "KIRC", "SKCM", "STAD", "LUSC", "BLCA", "BRCA" ,
             "ESCA","OV" ,"THCA", "COREAD", "UCEC", "LIHC")

all_data <- as.data.frame(t(IMs.exp))
all_data$embedding_cluster <- ifelse(cluster.info[,1]==0, "Cluster_0", "Cluster_1")
cancers.IMs.wilcox.test.res <- lapply(cancers, function(cancer){
  data <- all_data[which(cancer_samples[,1]==cancer),]
  IMs.wilcox.test.res <- do.call(rbind,lapply(rownames(IMs.exp), function(feature){
    score <- data[which(!is.na(data[,feature])), feature]
    embedding_cluster <- data[which(!is.na(data[,feature])), "embedding_cluster"]
    x <- as.numeric(score[which(embedding_cluster=="Cluster_0")])
    y <- as.numeric(score[which(embedding_cluster=="Cluster_1")])
    res <- WilcoxTest(x, y, paired = FALSE,types=c("Cluster_0","Cluster_1"), alternative = "two.sided")
    log2FC <- log2(mean(y, rm.na=T)/mean(x, rm.na=T))
    data.frame(z.statistic = res$z.statistic, p.value = res$p.value, log2FC=log2FC)
  }))
  rownames(IMs.wilcox.test.res) <- rownames(IMs.exp) # 这几个特征在组间均有显著差异
  IMs.wilcox.test.res$adj.pvalue <- p.adjust(IMs.wilcox.test.res$p.value, method="BH")
  IMs.wilcox.test.res
})
names(cancers.IMs.wilcox.test.res) <- cancers
save(cancers.IMs.wilcox.test.res, file=paste0(res_path, "cancers.IMs.wilcox.test.res.RData"))

library(pheatmap)
plot.data <- do.call(rbind,lapply(cancers.IMs.wilcox.test.res, function(x)x[,"log2FC"]))
colnames(plot.data) <- rownames(cancers.IMs.wilcox.test.res[[1]])
Immune_Checkpoint <- c("Stimulatory", "Inhibitory", "MHC_Class_I_or_II")
per.types.IMs <- lapply(Immune_Checkpoint, function(x)IMs.info[IMs.info[["HGNC Symbol"]] %in% colnames(plot.data) & IMs.info[["Immune Checkpoint"]] == x,"HGNC Symbol"])
IMs.order <- unlist(per.types.IMs)
plot.data <- plot.data[,IMs.order]


get.sigInfo <- function(p.value){ return(ifelse(p.value>0.05, "ns" ,""))}
sig.info <-  matrix(rep("", nrow(plot.data)*ncol(plot.data)), nrow(plot.data))
for(i in 1: nrow(plot.data)){
  sig.info[i, ] <-  get.sigInfo(cancers.IMs.wilcox.test.res[[i]][IMs.order,"adj.pvalue"])
}

rownames(IMs.info) <- IMs.info[["HGNC Symbol"]]
annotation_col = data.frame(Immune_Checkpoint = factor(IMs.info[IMs.order, "Immune Checkpoint"]))
rownames(annotation_col) = unlist(IMs.order)
ann_colors = list( Immune_Checkpoint = c(Stimulatory = "#b2182b", Inhibitory = "#4d4d4d", MHC_Class_I_or_II = "#e0e0e0"))

bk <-  c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01))
colors = c(colorRampPalette(colors = c("#0072B5FF","white"))(length(bk)/2),colorRampPalette(colors = c("white","#E18727FF"))(length(bk)/2))
pheatmap::pheatmap(plot.data, scale = 'none', 
         cluster_row = FALSE, cluster_cols = FALSE, 
         treeheight_row = 0, treeheight_col = 0,
         gaps_col= cumsum(unlist(lapply(per.types.IMs,length)))[-3], 
         display_numbers =sig.info, fontsize_number=6, 
         breaks=bk,
         color = colors, 
         annotation_names_col =FALSE,
         border_color = "grey", fontsize_row = 5, fontsize_col = 6, legend = TRUE,
         legend_break=seq(-1, 1, by=0.5), fontsize=1,
         height = 4, width = 12, 
         filename = paste0(res_path, "cancers_IMs_FC_heatmap.pdf"))


###  
#' 
#' the difference of immune related signature score between these two clusters
################################
# load immune related signatures
library(openxlsx)
TME_Fges <- read.xlsx("data/TCGA_data/raw_data/1-s2.0-S1535610821002221-mmc2.xlsx", sheet="SolidBlood Cancer (GMT)", startRow =2)

load("data/TCGA_data/pan.cancer.exp.data.RData")

#######
# 
#'  calculates the individual phenotypic scores for each sample
#######
library(limma)
library(org.Hs.eg.db)
source("code/Rpipeline/function/ssGSEA.R")

# 1) update symbol
Fges <- unique(TME_Fges[,2])
TME_Fges_list <- lapply(Fges, function(signature){
  gene_list <- na.omit(as.character(TME_Fges[which(TME_Fges[,2]==signature),3:ncol(TME_Fges)]))
  update.symbols <- alias2Symbol(gene_list,species = "Hs",expand.symbols = F)
  update.symbols
})
names(TME_Fges_list) <- Fges


# 2) calculate ssGSEA score
# TCGA
exp.matrix<- pan.cancer.exp.data 
score <- ssGSEA_Score(expr=as.matrix(exp.matrix), sig.list = TME_Fges_list) # 29 8755
save(score, file=paste0(res_path, "TME_Fges_score.RData"))

# ExpO
val_score <- ssGSEA_Score(expr=as.matrix(t(val_exp_data)), sig.list = TME_Fges_list) # 29 8755
save(val_score, file=paste0(res_path, "val_TME_Fges_score.RData"))



#######
# 
#' the difference between these two clusters
#######

source("code/RPipeline/function/WilcoxTest.R")
data <- as.data.frame(t(score))
data$embedding_cluster <- ifelse(cluster.info[,1]==0, "Cluster_0", "Cluster_1")

imm.wilcox.test.res<- do.call(rbind,lapply(rownames(score), function(feature){
  score <- data[which(!is.na(data[,feature])), feature]
  embedding_cluster <- data[which(!is.na(data[,feature])), "embedding_cluster"]
  x <- as.numeric(score[which(embedding_cluster=="Cluster_0")])
  y <- as.numeric(score[which(embedding_cluster=="Cluster_1")])
  res <- WilcoxTest(x, y, paired = FALSE,types=c("Cluster_0","Cluster_1"), alternative = "two.sided")
  res
}))
rownames(imm.wilcox.test.res) <- rownames(score) 
imm.wilcox.test.res$adj.pvalue <- p.adjust(imm.wilcox.test.res$p.value, method="BH")
save(imm.wilcox.test.res, file=paste0(res_path, "TME_Fges_score.wilcox.test.res.RData"))

data <- as.data.frame(t(val_score))
data$embedding_cluster <- ifelse(val.cluster.info[,1]==0, "Cluster_0", "Cluster_1")
imm.wilcox.test.res<- do.call(rbind,lapply(rownames(val_score), function(feature){
  score <- data[which(!is.na(data[,feature])), feature]
  embedding_cluster <- data[which(!is.na(data[,feature])), "embedding_cluster"]
  x <- as.numeric(score[which(embedding_cluster=="Cluster_0")])
  y <- as.numeric(score[which(embedding_cluster=="Cluster_1")])
  res <- WilcoxTest(x, y, paired = FALSE,types=c("Cluster_0","Cluster_1"), alternative = "two.sided")
  res
}))
rownames(imm.wilcox.test.res) <- rownames(val_score) 
imm.wilcox.test.res$adj.pvalue <- p.adjust(imm.wilcox.test.res$p.value, method="BH")
save(imm.wilcox.test.res, file=paste0(res_path, "val_TME_Fges_score.wilcox.test.res.RData"))


#######
# 
#' visualization
#######
plot.data <- t(scale(t(score)))
plot.data[which(plot.data>5)] <- 5
plot.data[which(plot.data<(-5))] <- -5
library(ComplexHeatmap)
pdf(paste0(res_path, "TME_Fges_score_heatmap.pdf"), height=6, width=12)
top_annotation = HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = c("#2166ac","#b2182b")),
                                                        labels = c("Cluster_0", "Cluster_1"),
                                                        labels_gp = gpar(col = "black")))
Heatmap(plot.data,
        col = colorRampPalette(c("navy","white","firebrick3"))(100),
        show_row_names = T,show_column_dend = F,show_row_dend = F,
        top_annotation = top_annotation,
        column_split = factor(cluster.info[,1]),
        row_split = c(rep(1, 11), rep(2, 10), rep(3, 6), rep(4,2)), 
        column_title = NULL,
        show_column_names = F)
dev.off()

plot.data <- t(scale(t(val_score)))
plot.data[which(plot.data>5)] <- 5
plot.data[which(plot.data<(-5))] <- -5
library(ComplexHeatmap)
pdf(paste0(res_path, "val_TME_Fges_score_heatmap.pdf"), height=6, width=12)
top_annotation = HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = c("#2166ac", "#b2182b")),
                                                        labels = c("Cluster_0","Cluster_1"),
                                                        labels_gp = gpar(col = "black")))
Heatmap(plot.data,
        col = colorRampPalette(c("navy","white","firebrick3"))(100),
        show_row_names = T,show_column_dend = F,show_row_dend = F,
        top_annotation = top_annotation,
        column_split = factor(val.cluster.info[,1]),
        row_split = c(rep(1, 11), rep(2, 10), rep(3, 6), rep(4,2)), 
        column_title = NULL,
        show_column_names = F)
dev.off()

#######
# 
#' The difference in immune landscape across different cancer types
#######
cancers <- c("HNSC", "LUAD", "CESC", "KIRC", "SKCM", "STAD", "LUSC", "BLCA", "BRCA" ,
             "ESCA","OV" ,"THCA", "COREAD", "UCEC", "LIHC")

all_data <- as.data.frame(t(score))
all_data$embedding_cluster <- ifelse(cluster.info[,1]==0, "Cluster_0", "Cluster_1")
cancers.TME_Fges.wilcox.test.res <- lapply(cancers, function(cancer){
  data <- all_data[which(cancer_samples[,1]==cancer),]
  wilcox.test.res <- do.call(rbind,lapply(rownames(score), function(feature){
    #score <- data[which(!is.na(data[,feature])), feature]
    embedding_cluster <- data[, "embedding_cluster"]
    x <- as.numeric(data[which(embedding_cluster=="Cluster_1") , feature])
    y <- as.numeric(data[which(embedding_cluster=="Cluster_0"), feature])
    res <- WilcoxTest(x, y, paired = FALSE,types=c("Cluster_0","Cluster_1"), alternative = "two.sided")
    data.frame(z.statistic = res$z.statistic, p.value = res$p.value)
  }))
  rownames(wilcox.test.res) <- rownames(score) 
  wilcox.test.res$adj.pvalue <- p.adjust(wilcox.test.res$p.value, method="BH")
  wilcox.test.res$log10.adj.p <- -log10(wilcox.test.res$adj.pvalue)
  wilcox.test.res$log10.adj.p[which(wilcox.test.res$z.statistic<0)] <- -(wilcox.test.res$log10.adj.p[which(wilcox.test.res$z.statistic<0)])
  wilcox.test.res$log10.adj.p[which(wilcox.test.res$log10.adj.p>10)] <- 10
  wilcox.test.res
})
names(cancers.TME_Fges.wilcox.test.res) <- cancers
save(cancers.TME_Fges.wilcox.test.res, file=paste0(res_path, "cancers.TME_Fges.wilcox.test.res.RData"))


library(pheatmap)
plot.data <- do.call(rbind,lapply(cancers.TME_Fges.wilcox.test.res, function(x)x[,"log10.adj.p"]))
colnames(plot.data) <- rownames(score)

sig.info <-  ifelse(abs(plot.data)<abs(log10(0.05)), "ns", "")

bk <-  c(seq(-10,-0.1,by=0.01),seq(0,10,by=0.01))
colors = c(colorRampPalette(colors = c("#0072B5FF","white"))(length(bk)/2),colorRampPalette(colors = c("white","#E18727FF"))(length(bk)/2))
pheatmap::pheatmap(plot.data, scale = 'none', 
         cluster_row = FALSE, cluster_cols = FALSE, 
         treeheight_row = 0, treeheight_col = 0,
         gaps_col= c(11, 21, 27), 
         display_numbers =sig.info,fontsize_number=6, 
         breaks=bk,
         color = colors, 
         annotation_names_col =FALSE,
         border_color = "grey", fontsize_row = 5, fontsize_col = 6, legend = TRUE,
         legend_break=seq(-10, 10, by=5), fontsize=1,
         height = 4, width = 8, 
         filename = paste0(res_path, "cnacers_TME_Fges_FC_heatmap.pdf"))
