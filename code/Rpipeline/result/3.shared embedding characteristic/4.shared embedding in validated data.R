
setwd("../../")

res_path <- "result/last_shared_layer_related_to_immune_phenotypes/validate_data/"
output_dir <- file.path(res_path)
if (!dir.exists(output_dir))  dir.create(output_dir)


###  
#' Expression Project for Oncology (expO)
################################
val_res_path = "result/model_validation/ExpO_project/"

# Visualize how the phenotype scores predicted by the model map on the embedding layer
val_predict_score <- read.csv(paste0(val_res_path, 'rescaled_final_model/val_preds.csv'), header = T, row.names=1)
last_layer_activations <- read.table(paste0(val_res_path, 'final_embeddings/0/val_rescale.txt'), header = F,sep=" ")

val_predict_score_percentage <- apply(val_predict_score, 2, function(x){
  y = rep(0,length(x))
  for(i in 1:length(x)) {
    y[i] = round(length(which(x<=x[i]))/length(x)*100,4)}
  y
})
mean_percentage <- apply(val_predict_score_percentage, 1, mean)
val_predict_score_percentage <- cbind(val_predict_score_percentage, mean_percentage)

library(Rtsne)
data <- last_layer_activations
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

library(ggplot2)
pdf(paste0(res_path, 'ExpO_project_t-SNE.pdf'))
for( i in 1:5){
  Phenotype_percentile <- val_predict_score_percentage[,i]
  p=ggplot(tsne_result,aes(tSNE1,tSNE2)) +
    geom_point(aes(colour= Phenotype_percentile), size=3) +
    scale_color_gradientn(colors=c("#08306b", "#2171b5","white", "#d7301f", "#7f0000"))+ 
    ggtitle(colnames(val_predict_score_percentage)[i]) +
    theme_bw()+ theme(panel.grid=element_blank(), 
                      plot.title = element_text(hjust = 0.5), 
                      legend.position="top")
  print(p)
}
dev.off()





###  
#' Braun.NatMed.2021
################################
val_res_path = "result/model_validation/Braun.NatMed.2021/"


#######
# 
#' Visualize how the phenotype scores predicted by the model map on the embedding layer
#######
val_predict_score <- read.csv(paste0(val_res_path, 'rescaled_final_model/val_preds.csv'), header = T, row.names=1)
last_layer_activations <- read.table(paste0(val_res_path, 'final_embeddings/0/val_rescale.txt'), header = F,sep=" ")
val_clinical_data <- read.csv("data/validation_data/Braun.NatMed.2021/processed_data/sub_clinical_data.csv")

val_predict_score_percentage <- apply(val_predict_score, 2, function(x){
  y = rep(0,length(x))
  for(i in 1:length(x)) {
    y[i] = round(length(which(x<=x[i]))/length(x)*100,4)}
  y
})
mean_percentage <- apply(val_predict_score_percentage, 1, mean)
val_predict_score_percentage <- cbind(val_predict_score_percentage, mean_percentage)

library(Rtsne)
data <- last_layer_activations
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

library(ggplot2)
pdf(paste0(res_path, 'Braun.NatMed.2021_t-SNE.pdf'))
for( i in 1:5){
  Phenotype_percentile <- val_predict_score_percentage[,i]
  p=ggplot(tsne_result,aes(tSNE1,tSNE2)) +
    geom_point(aes(colour= Phenotype_percentile), size=3) +
    scale_color_gradientn(colors=c("#08306b", "#2171b5","white", "#d7301f", "#7f0000"))+ 
    ggtitle(colnames(val_predict_score_percentage)[i]) +
    theme_bw()+ theme(panel.grid=element_blank(), 
                      plot.title = element_text(hjust = 0.5), 
                      legend.position="top")
  print(p)
}
dev.off()


pdf(paste0(res_path, 'Braun.NatMed.2021_tsne_ImmunePhenotypes.pdf'))
val_clinical_data$ImmunoPhenotype <- factor(val_clinical_data$ImmunoPhenotype, levels = c("No_IF", "Infiltrated","Excluded", "Desert"))
p=ggplot(tsne_result,aes(tSNE1,tSNE2)) +
  geom_point(aes(colour= val_clinical_data$ImmunoPhenotype), size=3) +
  scale_color_manual(values =c("#999999", "#d53e4f","#66c2a5","#4575b4"))+ 
  ggtitle("ImmunoPhenotype") +
  theme_bw()+ theme(panel.grid=element_blank(), 
                    plot.title = element_text(hjust = 0.5), 
                    legend.position="top")
print(p)
dev.off()


pdf(paste0(res_path, 'Braun.NatMed.2021_tsne_benefit.pdf'))
Benefit = val_clinical_data$Benefit
Benefit[which(val_clinical_data$Arm=="EVEROLIMUS")] <- "No_IF"
Benefit <- factor(Benefit, levels = c("No_IF", "CB", "ICB", "NCB"))
p=ggplot(tsne_result,aes(tSNE1,tSNE2)) +
  geom_point(aes(colour= Benefit), size=3) +
  scale_color_manual(values =c("#999999", "#d53e4f","#4575b4", "#66c2a5"))+ 
  ggtitle("Benefit") +
  theme_bw()+ theme(panel.grid=element_blank(), 
                    plot.title = element_text(hjust = 0.5), 
                    legend.position="top")
print(p)
dev.off()



###  
#' Mariathasan.Nature.2018
################################
val_res_path = "result/model_validation/Mariathasan.Nature.2018/"

#######
# 
#' Visualize how the phenotype scores predicted by the model map on the embedding layer
#######

val_predict_score <- read.csv(paste0(val_res_path, 'rescaled_final_model/val_preds.csv'), header = T, row.names=1)

last_layer_activations <- read.table(paste0(val_res_path, 'final_embeddings/0/val_rescale.txt'), header = F,sep=" ")

val_clinical_data <- read.csv("data/validation_data/Mariathasan.Nature.2018/processed_data/sub_clinical_data.csv")

val_predict_score_percentage <- apply(val_predict_score, 2, function(x){
  y = rep(0,length(x))
  for(i in 1:length(x)) {
    y[i] = round(length(which(x<=x[i]))/length(x)*100,4)}
  y
})
mean_percentage <- apply(val_predict_score_percentage, 1, mean)
val_predict_score_percentage <- cbind(val_predict_score_percentage, mean_percentage)

library(Rtsne)
data <- last_layer_activations
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
library(ggplot2)
pdf(paste0(res_path,'Mariathasan.Nature.2018_t-SNE.pdf'))
for( i in 1:5){
  Phenotype_percentile <- val_predict_score_percentage[,i]
  p=ggplot(tsne_result,aes(tSNE1,tSNE2)) +
    geom_point(aes(colour= Phenotype_percentile), size=3) +
    scale_color_gradientn(colors=c("#08306b", "#2171b5","white", "#d7301f", "#7f0000"))+ 
    ggtitle(colnames(val_predict_score_percentage)[i]) +
    theme_bw()+ theme(panel.grid=element_blank(), 
                      plot.title = element_text(hjust = 0.5), 
                      legend.position="top")
  print(p)
}
dev.off()


pdf(paste0(res_path,'Mariathasan.Nature.2018_t-SNE_ImmunePhenotypes.pdf'))

val_clinical_data$Immune.phenotype[which(is.na(val_clinical_data$Immune.phenotype))] <- "No_IF"
val_clinical_data$Immune.phenotype <- factor(val_clinical_data$Immune.phenotype, levels = c("No_IF", "inflamed","excluded", "desert"))
p=ggplot(tsne_result,aes(tSNE1,tSNE2)) +
  geom_point(aes(colour= val_clinical_data$Immune.phenotype), size=3) +
  scale_color_manual(values =c("#999999", "#d53e4f", "#66c2a5","#4575b4"))+ 
  ggtitle("ImmunoPhenotype") +
  theme_bw()+ theme(panel.grid=element_blank(), 
                    plot.title = element_text(hjust = 0.5), 
                    legend.position="top")
print(p)
dev.off()

pdf(paste0(res_path,'Mariathasan.Nature.2018_t-SNE_IC.pdf'))
val_clinical_data$IC.Level <- factor(val_clinical_data$IC.Level, levels = c("IC0", "IC1","IC2+"))
p=ggplot(tsne_result,aes(tSNE1,tSNE2)) +
  geom_point(aes(colour= val_clinical_data$IC.Level), size=3) +
  scale_color_manual(values =c("#A7A7A7", "#AEDBD7", "#344B9E"))+ 
  ggtitle("IC") +
  theme_bw()+ theme(panel.grid=element_blank(), 
                    plot.title = element_text(hjust = 0.5), 
                    legend.position="top")
print(p)
dev.off()
pdf(paste0(res_path,'Mariathasan.Nature.2018_t-SNE_TC.pdf'))
val_clinical_data$TC.Level <- factor(val_clinical_data$TC.Level, levels = c("TC0", "TC1","TC2+"))
p=ggplot(tsne_result,aes(tSNE1,tSNE2)) +
  geom_point(aes(colour= val_clinical_data$TC.Level), size=3) +
  scale_color_manual(values =c("#A7A7A7", "#AEDBD7", "#344B9E"))+ 
  ggtitle("TC") +
  theme_bw()+ theme(panel.grid=element_blank(), 
                    plot.title = element_text(hjust = 0.5), 
                    legend.position="top")
print(p)
dev.off()



pdf(paste0(res_path,'Mariathasan.Nature.2018_t-SNE_response.pdf'))

val_clinical_data$Best.Confirmed.Overall.Response <- factor(val_clinical_data$Best.Confirmed.Overall.Response, levels = c("CR", "PD"))
p=ggplot(tsne_result,aes(tSNE1,tSNE2)) +
  geom_point(aes(colour= val_clinical_data$Best.Confirmed.Overall.Response), size=3) +
  scale_color_manual(values =c("#5ab4ac","#d8b365"))+ 
  ggtitle("Response") +
  theme_bw()+ theme(panel.grid=element_blank(), 
                    plot.title = element_text(hjust = 0.5), 
                    legend.position="top")
print(p)
dev.off()
