
setwd("../../")

res_path <- "result/embedding_clusters_biology_survival_diff/diff_in_survival/"


###  
#' Molecular stratification of metastatic melanoma using gene expression profiling: Prediction of survival outcome and benefit from molecular targeted therapy.
#' 215 metastatic melanoma， microarray profiles(quantile-normalization)
#' GSE65904
################################

clinical_data <- read.csv("data/validation_data/Cirenajwis.Oncotarget.2015/processed_data/sub_clinical_data.csv",
                          header=T, check.names = F)

cluster.info <- read.csv("result/last_shared_layer_clustering/Cirenajwis.Oncotarget.2015/y_pred_2.csv", row.names = 1)
#0   1 
# 42 172
clinical_data$embedding_cluster <- ifelse(cluster.info[,1]==0, "Cluster_1", "Cluster_0")

library(survival)
library(survminer)

#######1) DSS

outcome <- subset(clinical_data, select = c("gender", "age","disease_specific_survival", "disease_specific_survival_in_days", "embedding_cluster"))
colnames(outcome) <- c("gender", "age","status","time", "embedding_cluster")
samples <- rownames(clinical_data)
if(any(is.na(outcome)) | length(which(outcome$time==0))>0){
  na_row <- which(!complete.cases(outcome))
  time0_row <- which(outcome$time==0)
  outcome <- outcome[-c(na_row, time0_row),]
  samples <- clinical_data$samples[-c(na_row, time0_row)]
}
table(outcome$embedding_cluster)
#Cluster_0 Cluster_1 
#170        39
cox_new <- coxph(Surv(time, status) ~ embedding_cluster, data =outcome)
summary(cox_new)
# HR = 0.4270(0.2382-0.7655),P=0.00428

pdf(paste0(res_path,"GSE65904_dss_cox_KM_forestplot.pdf"))

ggsurvplot(survfit(Surv(time, status) ~ embedding_cluster, data =outcome),
           data=outcome,legend = c(0.8,0.8),legend.title="",
           palette=c("#377eb8", "#e41a1c"),
           pval = TRUE, 
           xlab = "Follow up time(d)",  
           ylab = "Disease Specific Survival Proportion(%)",
           title = paste0("GSE65904, N=", nrow(outcome)))
cox.res <- coxph(Surv(time, status) ~ embedding_cluster+gender+age, data =outcome)
summary(cox.res)
plot.data <- summary(cox.res)$coef[,c(2, 5)]
plot.data <- as.data.frame(cbind(plot.data, summary(cox.res)$conf.int[, c(3,4)]))
plot.data$Varnames <- factor(rownames(plot.data))
colnames(plot.data) <- c("HR", "p.value", "Lower", "Upper","Varnames")
p <- ggplot(plot.data, aes(HR, Varnames))
p + geom_point(size=3) +
  geom_errorbarh(aes(xmax =Upper, xmin = Lower), height = 0.1) +
  scale_x_continuous(limits= c(0, 2), breaks= seq(0, 2, 0.5)) +
  geom_vline(aes(xintercept = 1)) +
  xlab('Hazard Ratio') + ylab(' ') + theme_bw()
dev.off()


####### distant_metastasis_free_survival

outcome <- subset(clinical_data, select = c("gender", "age","distant_metastasis_free_survival", "distant_metastasis_free_survival_in_days", "embedding_cluster"))
colnames(outcome) <- c("gender", "age","status","time", "embedding_cluster")
samples <- rownames(clinical_data)
if(any(is.na(outcome)) | length(which(outcome$time==0))>0){
  na_row <- which(!complete.cases(outcome))
  time0_row <- which(outcome$time==0)
  outcome <- outcome[-c(na_row, time0_row),]
  samples <- clinical_data$samples[-c(na_row, time0_row)]
}
table(outcome$embedding_cluster)
#Cluster_0 Cluster_1 
#115        34
cox_new <- coxph(Surv(time, status) ~ embedding_cluster, data =outcome)
summary(cox_new)
# HR = 0.4841(0.2673-0.8768),P=0.0167

pdf(paste0(res_path,"GSE65904_DMFS_cox_KM_forestplot.pdf"))

ggsurvplot(survfit(Surv(time, status) ~ embedding_cluster, data =outcome),
           data=outcome,legend = c(0.8,0.8),legend.title="",
           palette=c("#377eb8", "#e41a1c"),
           pval = TRUE, 
           xlab = "Follow up time(d)",  
           ylab = "Distant Metastasis Free Survival Proportion(%)",
           title = paste0("GSE65904, N=", nrow(outcome)))
cox.res <- coxph(Surv(time, status) ~ embedding_cluster+gender+age, data =outcome)
summary(cox.res)
plot.data <- summary(cox.res)$coef[,c(2, 5)]
plot.data <- as.data.frame(cbind(plot.data, summary(cox.res)$conf.int[, c(3,4)]))
plot.data$Varnames <- factor(rownames(plot.data))
colnames(plot.data) <- c("HR", "p.value", "Lower", "Upper","Varnames")
p <- ggplot(plot.data, aes(HR, Varnames)) 
p + geom_point(size=3) +
  geom_errorbarh(aes(xmax =Upper, xmin = Lower), height = 0.1) +
  scale_x_continuous(limits= c(0, 3.2), breaks= seq(0, 2, 0.5)) +
  geom_vline(aes(xintercept = 1)) +
  xlab('Hazard Ratio') + ylab(' ') + theme_bw()
dev.off()
