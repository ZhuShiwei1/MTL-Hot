setwd("../../")



load("data/TCGA_data/TCGA_clinical_data.RData")
cluster.info <- read.csv("result/last_shared_layer_clustering/TCGA/0/y_pred_2.csv", row.names = 1)
TCGA_clinical_data$embedding_cluster <- cluster.info[,1]
TCGA_clinical_data$embedding_cluster[which(TCGA_clinical_data$embedding_cluster==0)] <- "Cluster_0"
TCGA_clinical_data$embedding_cluster[which(TCGA_clinical_data$embedding_cluster==1)] <- "Cluster_1"
TCGA_clinical_data$embedding_cluster <- factor(TCGA_clinical_data$embedding_cluster)
#Cluster_0 Cluster_1 
#7032      1723 


res_path <- "result/embedding_clusters_immune_characterizing/cancer_specific_molecular_diff/"
output_dir <- file.path(res_path)
if (!dir.exists(output_dir))  dir.create(output_dir)

###  
#' 
#' Hypergeometric tests determined the enrichment of each cluster in each cancer type
################################
TCGA_clinical_data$cancer_type[which(TCGA_clinical_data$cancer_type %in% c("COAD", "READ"))] <- "COREAD"
hyper_p <- lapply(unique(TCGA_clinical_data$cancer_type), function(cancer){
  samples_num <- nrow(TCGA_clinical_data)
  Cluster_0_total <- length(which(TCGA_clinical_data$embedding_cluster=="Cluster_0"))
  cancer_num <- length(which(TCGA_clinical_data$cancer_type==cancer))
  Cluster_0_in_cancer <- length(which(TCGA_clinical_data$cancer_type==cancer & TCGA_clinical_data$embedding_cluster=="Cluster_0"))
  Cluster_1_total <- samples_num- Cluster_0_total
  Cluster_1_in_cancer <- cancer_num-Cluster_0_in_cancer
  Cluster_0_p <- phyper(Cluster_0_in_cancer-1,m=Cluster_0_total,n=Cluster_1_total,k=cancer_num, lower.tail=F)
  
  Cluster_1_p <- phyper(Cluster_1_in_cancer-1, Cluster_1_total, Cluster_0_total, cancer_num, lower.tail=F)

  data.frame(Cluster_0=Cluster_0_p,Cluster_1=Cluster_1_p)
})
names(hyper_p) <- unique(TCGA_clinical_data$cancer_type)
hyper_p <- t(do.call(rbind, hyper_p))


adj.p <- p.adjust(c(hyper_p[1,],hyper_p[2,]) , method = "BH")
hyper_p_adj <- hyper_p
hyper_p_adj[1,] <- -log10(adj.p[1:(length(adj.p)/2)])
hyper_p_adj[2,] <- -log10(adj.p[(length(adj.p)/2+1):(length(adj.p))])

# The number of samples of two clusters in each cancer type
clusters_in_cancers <- do.call(cbind,lapply(unique(TCGA_clinical_data$cancer_type), function(cancer){
  table(subset(TCGA_clinical_data, cancer_type==cancer)$embedding_cluster)
}))
colnames(clusters_in_cancers) <- unique(TCGA_clinical_data$cancer_type)

order_data <- hyper_p_adj
order_data[1,] <- -order_data[1,]
plot_order <- order(colSums(order_data), decreasing = T)

# visualization
library(pheatmap) 
bk <- seq(0,10,by=0.001)  
pheatmap(hyper_p_adj[,plot_order], 
         cluster_row = FALSE,sclae ="none",cluster_cols=FALSE,
         color = colorRampPalette(colors = c("#91bfdb","#ffffbf","#fc8d59"))(length(bk)),
         breaks=bk,
         display_numbers = clusters_in_cancers[,plot_order],
         filename = paste0(res_path, "hyper_cancers_heatmap.pdf"),
         height=2,
         width=12)

###  
#' 
#'Association of different clusters within different cancer types with viruses or MSI
################################
pdf(paste0(res_path,"specific_cancer_MSI.pdf"), width=4, height=4)
colors <- c("#e34a33", "#fee8c8")
plot.data <- TCGA_clinical_data[which(TCGA_clinical_data$cancer_type == "COREAD"),c("MSI", "embedding_cluster")]
plot.data$MSI <- factor(plot.data$MSI, levels = c("MSI", "MSS"))
plot.data$embedding_cluster <- factor(plot.data$embedding_cluster)
plot.data <- plot.data[complete.cases(plot.data),]
ggplot(plot.data,aes(embedding_cluster,fill=MSI))+
  geom_bar(stat="count",position="fill")+ scale_fill_manual(values= colors)+
  theme(axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+ 
  labs(x="", y="Fraction of cases(%)", title="COREAD")
chisq.test(plot.data$MSI, plot.data$embedding_cluster) #  p-value =1.631e-10
table(plot.data$embedding_cluster) # 317        35 

plot.data <- TCGA_clinical_data[which(TCGA_clinical_data$cancer_type=="STAD"),c("MSI", "embedding_cluster")]
plot.data$MSI <- factor(plot.data$MSI, levels = c("MSI", "MSS"))
plot.data$embedding_cluster <- factor(plot.data$embedding_cluster)
plot.data <- plot.data[complete.cases(plot.data),]
ggplot(plot.data,aes(embedding_cluster,fill=MSI))+
  geom_bar(stat="count",position="fill")+ scale_fill_manual(values= colors)+
  theme(axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+ 
  labs(x="", y="Fraction of cases(%)", title="STAD")
chisq.test(plot.data$MSI, plot.data$embedding_cluster) #  p-value =1.378e-06
table(plot.data$embedding_cluster) # 285       123
dev.off()


pdf(paste0(res_path,"specific_cancer_virus.pdf"), width=4, height=4)
colors <- c("#8856a7", "#e0ecf4")
plot.data <- TCGA_clinical_data[which(TCGA_clinical_data$cancer_type =="HNSC"), c("HPV", "embedding_cluster")]
plot.data$HPV <- factor(plot.data$HPV, levels = c("Positive", "Negative"))
plot.data$embedding_cluster <- factor(plot.data$embedding_cluster)
plot.data <- plot.data[complete.cases(plot.data),]
ggplot(plot.data,aes(embedding_cluster,fill=HPV))+
  geom_bar(stat="count",position="fill")+ scale_fill_manual(values= colors)+
  theme(axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+ 
  labs(x="", y="Fraction of cases(%)", title="HNSC(HPV)")
chisq.test(plot.data$HPV, plot.data$embedding_cluster) #  p-value =0.0009385
table(plot.data$embedding_cluster) # 206       178 

plot.data <- TCGA_clinical_data[which(TCGA_clinical_data$cancer_type=="STAD"),c("EBV", "embedding_cluster")]
plot.data$EBV <- factor(plot.data$EBV, levels = c("Positive", "Negative"))
plot.data$embedding_cluster <- factor(plot.data$embedding_cluster)
plot.data <- plot.data[complete.cases(plot.data),]
ggplot(plot.data,aes(embedding_cluster,fill=EBV))+
  geom_bar(stat="count",position="fill")+ scale_fill_manual(values= colors)+
  theme(axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+ 
  labs(x="", y="Fraction of cases(%)", title="STAD(EBV)")
chisq.test(plot.data$EBV, plot.data$embedding_cluster) #  p-value = 6.622e-07
table(plot.data$embedding_cluster) # 165        84 

dev.off()

