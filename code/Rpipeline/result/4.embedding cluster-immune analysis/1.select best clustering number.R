
setwd("../../")
library(fpc)
library(ggplot2)


last_layer_activations <- lapply(0:99,function(rep){
  read.table(paste0('result/model_training/MTL/MTL_Hot/final_embeddings/',rep, '.txt'), header = F,sep=" ")
})
names(last_layer_activations) <- 1:100
save(last_layer_activations, file="result/model_training/MTL/MTL_Hot/final_embeddings/last_layer_activations.RData")


last_layer_clusters <- lapply(0:99, function(rep){
  cluster_res <- list()
  for(cluster_num in 2:5){
    per.res <- read.csv(paste0("result/last_shared_layer_clustering/TCGA/", rep, "/y_pred_", cluster_num, ".csv"), row.names=1)
    cluster_res <- c(cluster_res, list(per.res[,1]))
  }
  cluster_res
})
names(last_layer_clusters) <- 1:100
save(last_layer_clusters, file="result/last_shared_layer_clustering/TCGA/last_layer_clusters.RData")


###  
#' 
#'  evaluate the clustering quality under different clustering numbers
################################
cluster_measures <- lapply(1:100, function(rep){
  last_layer <- last_layer_activations[[rep]] # 8755  100
  cluster_res <- last_layer_clusters[[rep]]
  distance_matrix <- dist(last_layer)
  measures <- do.call(rbind,lapply(cluster_res, function(per_res){
    stats <- cluster.stats(distance_matrix, per_res)
    data.frame(sli = stats$avg.silwidth, ch = stats$ch)
  }))
  measures
})
names(cluster_measures) <- 1:100
save(cluster_measures, file="result/last_shared_layer_clustering/TCGA/cluster_measures.RData")

best_cluster_num <- lapply(cluster_measures, function(measures){
  calinski.best <- c(2:6)[which.max(measures$ch)]
  sli.best <- c(2:6)[which.max(measures$sli)]
  list(best_cluster_num = calinski.best, sli.best=sli.best)
})
names(best_cluster_num) <- 1:100

table(unlist(lapply(best_cluster_num, function(x)x$best_cluster_num))) 

table(unlist(lapply(best_cluster_num, function(x)x$sli.best))) 

save(best_cluster_num, file="result/last_shared_layer_clustering/TCGA/best_cluster_num.RData")

# plot 
plot.data <- data.frame(Cluster_number = rep(2:5, times=100), 
                        times = rep(1:100, each=4), 
                        sli = unlist(lapply(cluster_measures, function(x)x$sli)),
                        ch = unlist(lapply(cluster_measures, function(x)x$ch)))
pdf("result/last_shared_layer_clustering/TCGA/cluster_measures.pdf", height=4, width=6)
ggplot(data = plot.data, mapping = aes(x = Cluster_number, y = sli,group=times)) + geom_line(colour = "grey")+
  xlab("Cluster numbers") + ylab("Silhouette index") + theme_bw()
ggplot(data = plot.data, mapping = aes(x = Cluster_number, y = ch,group=times)) + geom_line(colour = "grey")+
  xlab("Cluster numbers") + ylab("Calinski-Harabaz") + theme_bw()
dev.off()




###  
#' 
#' The classification of samples that are simultaneously clustered into the same class in 100 training sessions
################################
cluster_res_mat <- do.call(rbind, lapply(last_layer_clusters, function(x)x[[1]]))

cluster0_idx <- which(cluster_res_mat[1,]==0)
cluster1_idx <- which(cluster_res_mat[1,]==1)
cluster_res_mat <- apply(cluster_res_mat, 1, function(x){

  per_cluster0_idx <- which(x==0)
  per_cluster1_idx <- which(x==1)
  cluster00_overlap_num <- length(intersect(cluster0_idx, per_cluster0_idx))# 
  cluster10_overlap_num <- length(intersect(cluster1_idx, per_cluster0_idx))# 

  if(cluster00_overlap_num > cluster10_overlap_num){ 
    x[per_cluster0_idx] <- 0  
    x[per_cluster1_idx] <- 1
  }else{
    x[per_cluster1_idx] <- 0 
    x[per_cluster0_idx] <- 1
  }
  x
}) # 8755  100

diff.clu.count <- apply(cluster_res_mat, 1, function(x){length(unique(x))})
table(diff.clu.count)
#1    2 
#7896  859(90.18846%)
table(cluster_res_mat[which(diff.clu.count==1),1])  
#0    1 
#6702 1194 


plot.data <- data.frame(times = rep(1:100, each=2),
                        clusters = factor(rep(0:1, times=100)),
                        nums = as.vector(apply(cluster_res_mat, 2, table)))
pdf("result/last_shared_layer_clustering/TCGA/clusters_100run_areaPlot.pdf", height=5, width=12)
ggplot(plot.data, aes(x = times)) +  
  geom_area(aes(y = nums, fill = clusters), alpha = 0.4, position = 'fill') +
  scale_fill_manual(values = c("#2166ac", "#b2182b")) + theme_bw()
dev.off()

pdf("result/last_shared_layer_clustering/TCGA/clusters_100run_areaPlot_value.pdf", height=5, width=12)
ggplot(plot.data, aes(x = times)) +  
  geom_area(aes(y = nums, fill = clusters), alpha = 0.5, position = 'stack') +
  scale_fill_manual(values = c("#2166ac", "#b2182b")) + theme_bw()+
  geom_hline(aes(yintercept=1194), colour="#b2182b", linetype="dashed")+
  geom_hline(aes(yintercept=nrow(cluster_res_mat)-6702), colour="#2166ac", linetype="dashed")

dev.off()




###  
#' 
#' The stability of the model clustering results was evaluated in a set of pan-cancer validation sets-EXPO project
################################

last_layer_activations <- lapply(0:99,function(rep){
  read.table(paste0('result/model_validation/ExpO_project/final_embeddings/',rep, '/val_rescale.txt'), header = F,sep=" ")
})
names(last_layer_activations) <- 1:100
save(last_layer_activations, file="result/model_validation/ExpO_project/final_embeddings/last_layer_activations.RData")

last_layer_clusters <- lapply(0:99, function(rep){
  cluster_res <- list()
  for(cluster_num in 2:6){
    per.res <- read.csv(paste0("result/last_shared_layer_clustering/ExpO_project/", rep, "/y_pred_", cluster_num, ".csv"), row.names=1)
    cluster_res <- c(cluster_res, list(per.res[,1]))
  }
  cluster_res
})
names(last_layer_clusters) <- 1:100
save(last_layer_clusters, file="result/last_shared_layer_clustering/ExpO_project/last_layer_clusters.RData")


###  
#' 
#' evaluate the clustering quality under different clustering numbers
################################
cluster_measures <- lapply(1:100, function(rep){
  last_layer <- last_layer_activations[[rep]] # 8755  100
  cluster_res <- last_layer_clusters[[rep]]
  distance_matrix <- dist(last_layer)
  measures <- do.call(rbind,lapply(cluster_res, function(per_res){
    stats <- cluster.stats(distance_matrix, per_res)
   
    data.frame(sli = stats$avg.silwidth, ch = stats$ch)
  }))
  measures
})
names(cluster_measures) <- 1:100
save(cluster_measures, file="result/last_shared_layer_clustering/ExpO_project/cluster_measures.RData")


best_cluster_num <- lapply(cluster_measures, function(measures){
  calinski.best <- c(2:6)[which.max(measures$ch)]
  sli.best <- c(2:6)[which.max(measures$sli)]

  list(best_cluster_num = calinski.best, sli.best=sli.best)
})
names(best_cluster_num) <- 1:100
table(unlist(lapply(best_cluster_num, function(x)x$best_cluster_num))) 
#2  3 
#66 34
table(unlist(lapply(best_cluster_num, function(x)x$sli.best))) 
#2  3 
#74 26
save(best_cluster_num, file="result/last_shared_layer_clustering/ExpO_project/best_cluster_num.RData")


## plot
plot.data <- data.frame(Cluster_number = rep(2:6, times=100), 
                        times = rep(1:100, each=5), 
                        sli = unlist(lapply(cluster_measures, function(x)x$sli)),
                        ch = unlist(lapply(cluster_measures, function(x)x$ch)))
pdf("result/last_shared_layer_clustering/ExpO_project/cluster_measures.pdf", height=4, width=6)
ggplot(data = plot.data, mapping = aes(x = Cluster_number, y = sli,group=times)) + geom_line(colour = "grey")+
  xlab("Cluster numbers") + ylab("Silhouette index") + theme_bw()
ggplot(data = plot.data, mapping = aes(x = Cluster_number, y = ch,group=times)) + geom_line(colour = "grey")+
  xlab("Cluster numbers") + ylab("Calinski-Harabaz") + theme_bw()
dev.off()


cluster_res_mat <- do.call(rbind, lapply(last_layer_clusters, function(x)x[[1]]))
#100 2158

#### 1) Because the class name is different each time the cluster is defined, the class name is unified
cluster0_idx <- which(cluster_res_mat[1,]==0)
cluster1_idx <- which(cluster_res_mat[1,]==1)
cluster_res_mat <- apply(cluster_res_mat, 1, function(x){
  per_cluster0_idx <- which(x==0)
  per_cluster1_idx <- which(x==1)
  cluster00_overlap_num <- length(intersect(cluster0_idx, per_cluster0_idx))
  cluster10_overlap_num <- length(intersect(cluster1_idx, per_cluster0_idx))
  
  if(cluster00_overlap_num > cluster10_overlap_num){ 
    x[per_cluster0_idx] <- 0  
    x[per_cluster1_idx] <- 1
  }else{
    x[per_cluster1_idx] <- 0 
    x[per_cluster0_idx] <- 1
  }
  x
}) #  2158  100



plot.data <- data.frame(times = rep(1:100, each=2),
                        clusters = factor(rep(0:1, times=100)),
                        nums = as.vector(apply(cluster_res_mat, 2, table)))
pdf("result/last_shared_layer_clustering/ExpO_project/clusters_100run_areaPlot.pdf", height=5, width=12)
ggplot(plot.data, aes(x = times)) +  
  geom_area(aes(y = nums, fill = clusters), alpha = 0.4, position = 'fill') +
  scale_fill_manual(values = c("#2166ac", "#b2182b")) + theme_bw()
dev.off()

pdf("result/last_shared_layer_clustering/ExpO_project/clusters_100run_areaPlot_value.pdf", height=5, width=12)
ggplot(plot.data, aes(x = times)) +  
  geom_area(aes(y = nums, fill = clusters), alpha = 0.5, position = 'stack') +
  scale_fill_manual(values = c("#2166ac", "#b2182b")) + theme_bw()
dev.off()





