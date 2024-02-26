
setwd("../../")

res_path <- "result/last_shared_layer_related_to_immune_phenotypes/normed_KMeans_medoids/"
output_dir <- file.path(res_path)
if (!dir.exists(output_dir))  dir.create(output_dir)


load("data/TCGA_data/pan.cancer.phenotype.RData")
num_centres = 50
last_layer_consensus_nodes <- read.table(paste0("result/model_training/MTL/MTL_Hot/final_embeddings/Cluster_",num_centres,"_embedding.txt"), header=F, sep=" ")
rownames(last_layer_consensus_nodes) <- rownames(pan.cancer.phenotype)
colnames(last_layer_consensus_nodes) <- 1:num_centres

last_layer_activations <- lapply(0:99,function(rep){
  read.table(paste0('result/model_training/MTL/MTL_Hot/final_embeddings/',rep, '.txt'), header = F,sep=" ")
})
names(last_layer_activations) <- 1:100
save(last_layer_activations, file="result/model_training/MTL/MTL_Hot/final_embeddings/last_layer_activations.RData")


consensus_nodes_info <- read.table(paste0("result/model_training/MTL/MTL_Hot/final_embeddings/Cluster_",num_centres,"_medoids_info.csv"), header=T, sep=",")
consensus_nodes_output <- apply(consensus_nodes_info, 1, function(x){
  run <- x[1]
  node_idx <- x[2]
  last_layer_activations[[run+1]][,node_idx+1]
})
rownames(consensus_nodes_output) <- rownames(pan.cancer.phenotype)
save(consensus_nodes_output, file=paste0("result/model_training/MTL/MTL_Hot/final_embeddings/consensus_nodes_output_",num_centres,".RData"))


load("data/TCGA_data/TCGA_clinical_data.RData")
###  
#' 
#' The significant associations of the last shared layer with tumor immune-related features were calculated 
################################
data <- consensus_nodes_output
last_layer <- consensus_nodes_output
clinical_factors <- c("Leukocyte.Fraction","CYT_score", "Nonsilent.Mutation.Rate", 
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
                                    "Leukocyte.Fraction","CYT_score", "TMB",
                                    "Neoantigens")
spearman_corr_rho = do.call(rbind, lapply(spearman_corr_res, function(x)x$corr))
colnames(spearman_corr_rho) <- c("APM.score", "Tcell.score", "IFNgamma.score",  "PDL1.exp",
                                 "Leukocyte.Fraction","CYT_score", "TMB",
                                 "Neoantigens")


spearman_corr_adj_pvalue <- t(apply(spearman_corr_pvalue, 1, function(x)p.adjust(x, "BH")))
apply(spearman_corr_adj_pvalue, 2, function(x)length(which(x<0.05)))
#APM.score        Tcell.score     IFNgamma.score           PDL1.exp Leukocyte.Fraction          CYT_score 
#50                 48                 49                 50                 49                 48 
#TMB        Neoantigens 
#46                 47
save(spearman_corr_rho, spearman_corr_pvalue, spearman_corr_adj_pvalue, 
     file=paste0(res_path, "nodes_output_",num_centres,"_corr_imm.RData" ))



library(pheatmap)
plot.data <- t(spearman_corr_rho)
mark_matrix <- matrix(ifelse(t(spearman_corr_adj_pvalue) <= 0.05, " ","ns"), 8)

bk <- c(seq(-1,-0.001,by=0.001),seq(0,1,by=0.001))
pheatmap(plot.data,cluster_row = FALSE,cluster_cols = TRUE, 
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
         display_numbers = mark_matrix,
         breaks=bk,
         show_colnames = FALSE,scale="none",treeheight_col=0,
         filename =paste0(res_path,"imm_corr_heatmap_",num_centres,".pdf"),
         width=15,
         height=2.5
) 

plot.data <- t(-log10(spearman_corr_adj_pvalue))
plot.data[which(plot.data=="Inf")] <- 300
plot.data[which(plot.data>300)] <- 300
mark_matrix <- matrix(ifelse(plot.data < -log10(0.05), "ns"," "), 8)

pheatmap(plot.data,
         cluster_row = FALSE,cluster_cols = FALSE, 
         color = colorRampPalette(c("white", "#d73027"))(100),
         na_col = "white",
         display_numbers = mark_matrix,
         show_colnames = FALSE,scale="none",treeheight_col=0 ,
         filename = paste0(res_path,"imm_corr_pvalue_heatmap_",num_centres,".pdf"),
         width=15,
         height=2.5)
