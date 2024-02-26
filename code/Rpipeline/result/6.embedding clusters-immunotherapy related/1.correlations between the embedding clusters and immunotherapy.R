
setwd("../../")

res_path <- "result/embedding_clusters_immunotherapy_response/"


###  
#' Braun.NatMed.2021
################################
per_res_path <- paste0(res_path, "Braun.NatMed.2021/")
output_dir <- file.path(per_res_path)
if (!dir.exists(output_dir))  dir.create(output_dir)

clinical.data <- read.csv("data/validation_data/Braun.NatMed.2021/processed_data/sub_clinical_data.csv",
                          header=T)

cluster.info <- read.csv("result/last_shared_layer_clustering/Braun.NatMed.2021/y_pred_2.csv", row.names = 1)
#0   1 
#264  44

clinical.data$embedding_cluster <- ifelse(cluster.info[,1]==1, "Cluster_1", "Cluster_0")

sub.clinical.data <- subset(clinical.data, Arm=="NIVOLUMAB")
table(sub.clinical.data$embedding_cluster)
#0   1 
#153  25
library(survival)
library(survminer)
#######1) data preparation
outcome <- subset(clinical.data, Arm=="NIVOLUMAB", select = c("OS_CNSR", "OS", "embedding_cluster"))
colnames(outcome) <- c("status","time", "embedding_cluster")
samples <- rownames(clinical.data)
if(any(is.na(outcome)) | length(which(outcome$time==0))>0){
  outcome <- outcome[-c(which(!complete.cases(outcome)), which(outcome$time==0)),]
}
cox_new <- coxph(Surv(time, status) ~ embedding_cluster, data =outcome)
summary(cox_new)
# HR = 1.6493(1.019-2.669),P=0.0417

pdf(paste0(per_res_path, "KMcurve.pdf"))
ggsurvplot(survfit(Surv(time, status) ~ embedding_cluster, data =outcome),
           data=outcome,legend = c(0.8,0.8),legend.title="",
           palette=c("#377eb8", "#e41a1c"),
           pval = TRUE, 
           xlab = "Time(months)", 
           ylab = "Survival Proportion", 
           title = paste0("Anti-PD-1 therapy, ccRCC, N=", nrow(outcome)))
dev.off()

# PFS
outcome <- subset(clinical.data, Arm=="NIVOLUMAB", select = c("PFS_CNSR", "PFS", "embedding_cluster"))
colnames(outcome) <- c("status","time", "embedding_cluster")
samples <- rownames(clinical.data)
if(any(is.na(outcome)) | length(which(outcome$time==0))>0){
  outcome <- outcome[-c(which(!complete.cases(outcome)), which(outcome$time==0)),]
}
cox_new <- coxph(Surv(time, status) ~ embedding_cluster, data =outcome)
summary(cox_new)
# HR = 1.6493,P=0.0417

pdf(paste0(per_res_path, "PFS_KMcurve.pdf"))
ggsurvplot(survfit(Surv(time, status) ~ embedding_cluster, data =outcome),
           data=outcome,legend = c(0.8,0.8),legend.title="",
           palette=c("#377eb8", "#e41a1c"),
           pval = TRUE, 
           xlab = "Time(months)",  
           ylab = "PFS Proportion", 
           title = paste0("Anti-PD-1 therapy, ccRCC, N=", nrow(outcome)))
dev.off()


table(sub.clinical.data[,c("embedding_cluster","Benefit")])
#embedding_cluster CB ICB NCB
#Cluster_0 47  50  56
#Cluster_1  8   6  11


# Multivariate cox regression analysis
outcome <- subset(sub.clinical.data, select = c("TMB_Counts","Age","OS_CNSR", "OS", "embedding_cluster"))
colnames(outcome) <- c("TMB", "Age","status","time", "embedding_cluster")
if(any(is.na(outcome)) | length(which(outcome$time==0))>0){
  outcome <- outcome[-c(which(!complete.cases(outcome)), which(outcome$time==0)),]
}
cox.res <- coxph(Surv(time, status) ~ embedding_cluster+TMB+Age, data =outcome)
summary(cox.res)
# Forest map visualization
pdf(paste0(per_res_path, "cox_forestplot.pdf"), height=4, width=4)
plot.data <- summary(cox.res)$coef[,c(2, 5)]
plot.data <- as.data.frame(cbind(plot.data, summary(cox.res)$conf.int[, c(3,4)]))
plot.data$Varnames <- factor(rownames(plot.data))
colnames(plot.data) <- c("HR", "p.value", "Lower", "Upper","Varnames")
p <- ggplot(plot.data, aes(HR, Varnames)) 
  geom_errorbarh(aes(xmax =Upper, xmin = Lower), height = 0.1) +
  scale_x_continuous(limits= c(0, 3.5), breaks= seq(0, 3.5, 0.5)) +
  geom_vline(aes(xintercept = 1)) +
  xlab('Hazard Ratio') + ylab(' ') + theme_bw()
dev.off()


library(grid)
library(gridExtra)
library(plyr)
colors <- c("#1a9641","#a6d96a","#d7191c")
plot.data <- clinical.data[,c("Benefit", "embedding_cluster")]
plot.data$Benefit <- factor(plot.data$Benefit, levels = c("CB","ICB","NCB"))
plot.data$embedding_cluster <- factor(plot.data$embedding_cluster)
pdf(paste0(per_res_path,"correlated_with_response.pdf"), width=4, height=4)
ggplot(plot.data,aes(embedding_cluster,fill=Benefit))+
  geom_bar(stat="count",position="fill")+ scale_fill_manual(values= colors)+
  theme(axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        legend.title=element_blank(),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+ 
  labs(x="", y="Fraction of cases(%)", title="BLCA")
dev.off()

chisq.test(plot.data$Benefit, plot.data$embedding_cluster) #  p-value =  0.5987



###  
#' Mariathasan.Nature.2018
################################
per_res_path <- paste0(res_path, "Mariathasan.Nature.2018/")

output_dir <- file.path(per_res_path)
if (!dir.exists(output_dir))  dir.create(output_dir)

clinical.data <- read.csv("data/validation_data/Mariathasan.Nature.2018/processed_data/sub_clinical_data.csv",
                          header=T)

cluster.info <- read.csv("result/last_shared_layer_clustering/Mariathasan.Nature.2018/y_pred_2.csv", row.names = 1)
clinical.data$embedding_cluster <- ifelse(cluster.info[,1]==0, "Cluster_1", "Cluster_0")

library(survival)
library(survminer)
#######1) data preparation
outcome <- subset(clinical.data, select = c("censOS", "os", "embedding_cluster"))
colnames(outcome) <- c("status","time", "embedding_cluster")
if(any(is.na(outcome)) | length(which(outcome$time==0))>0){
  outcome <- outcome[-c(which(!complete.cases(outcome)), which(outcome$time==0)),]
}
cox_new <- coxph(Surv(time, status) ~ embedding_cluster, data =outcome)
summary(cox_new)
# HR = 1.9123(1.19-3.073),P=0.00739
# 0.4610(0.2942    0.7225) p=0.000731
# KM curves
pdf(paste0(per_res_path, "KMcurve.pdf"))
ggsurvplot(survfit(Surv(time, status) ~ embedding_cluster, data =outcome),
           data=outcome,legend = c(0.8,0.8),legend.title="",
           palette=c("#377eb8", "#e41a1c"),
           pval = TRUE, 
           xlab = "Time(months)", 
           ylab = "Survival Proportion", 
           title = paste0("Anti-PD-L1 therapy, BLCA, N=", nrow(outcome)))
dev.off()
table(clinical.data[,c("embedding_cluster", "Best.Confirmed.Overall.Response")])
#embedding_cluster  CR  NE  PD  PR  SD
#Cluster_0  15  44 150  39  52
#Cluster_1  10   6  17   4  11
# Multivariate cox regression analysis
outcome <- subset(clinical.data, select = c("FMOne.mutation.burden.per.MB","Tobacco.Use.History", "censOS", "os", "embedding_cluster"))
colnames(outcome) <- c("TMB", "Smoke","status","time", "embedding_cluster")
if(any(is.na(outcome)) | length(which(outcome$time==0))>0){
  outcome <- outcome[-c(which(!complete.cases(outcome)), which(outcome$time==0)),]
}
cox.res <- coxph(Surv(time, status) ~ embedding_cluster+TMB+Smoke, data =outcome)
summary(cox.res)
#0.40935(0.2442 -   0.6862) p=0.000701

# TMB stratification
outcome$TMB <- ifelse(outcome$TMB>median(outcome$TMB), "TMB_H", "TMB_L")
outcome$strat_TMB <- paste(outcome$TMB, outcome$embedding_cluster, sep="_")
outcome$strat_TMB <- factor(outcome$strat_TMB, levels=c("TMB_H_Cluster_1", "TMB_L_Cluster_1", "TMB_H_Cluster_0", "TMB_L_Cluster_0"))
# KM curves
pdf(paste0(per_res_path, "TMB_stratify_KMcurve.pdf"))
ggsurvplot(survfit(Surv(time, status) ~ TMB, data =outcome),
           data=outcome,legend = c(0.8,0.8),legend.title="",
           palette=c( "#d7191c","#2c7bb6"),
           pval = TRUE, 
           xlab = "Time(Months)",  
           ylab = "Survival Proportion", 
           title = paste0("Anti-PD-L1 therapy, BLCA, N=", nrow(outcome)))

ggsurvplot(survfit(Surv(time, status) ~ strat_TMB, data =outcome),
           data=outcome,legend = c(0.8,0.8),legend.title="",
           palette=c( "#d7191c","#fdae61", "#abd9e9","#2c7bb6"),
           pval = TRUE, 
           xlab = "Time(Months)",  
           ylab = "Survival Proportion", 
           title = paste0("Anti-PD-L1 therapy, BLCA, N=", nrow(outcome)))
dev.off()

# Forest map visualization
pdf(paste0(per_res_path, "cox_forestplot.pdf"), height=4, width=4)
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



library(grid)
library(gridExtra)
library(plyr)
colors <- c("#1a9641","#a6d96a","#fdae61","#d7191c")

plot.data <- clinical.data[,c("Best.Confirmed.Overall.Response", "embedding_cluster")]
plot.data <- plot.data[which(plot.data$Best.Confirmed.Overall.Response %in% c("PR", "CR", "PD", "SD")),]
plot.data$Best.Confirmed.Overall.Response <- factor(plot.data$Best.Confirmed.Overall.Response, levels = c("CR","PR","SD", "PD"))
plot.data$embedding_cluster <- factor(plot.data$embedding_cluster)
pdf(paste0(per_res_path,"correlated_with_response.pdf"), width=4, height=4)
ggplot(plot.data,aes(embedding_cluster,fill=Best.Confirmed.Overall.Response))+
  geom_bar(stat="count",position="fill")+ scale_fill_manual(values= colors)+
  theme(axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        legend.title=element_blank(),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+ 
  labs(x="", y="Fraction of cases(%)", title="BLCA")
dev.off()

chisq.test(plot.data$Best.Confirmed.Overall.Response, plot.data$embedding_cluster) #  p-value = 0.0005914
chisq.test(plot.data$Best.Confirmed.Overall.Response, plot.data$embedding_cluster,  simulate.p.value=T) #  p-value = 0.002499



###  
#' Cui.NPJGenomMed.2016
################################
per_res_path <- paste0(res_path, "Cui.NPJGenomMed.2016/")
output_dir <- file.path(per_res_path)
if (!dir.exists(output_dir))  dir.create(output_dir)

clinical.data <- read.csv("data/validation_data/Cui.NPJGenomMed.2016/processed_data/sub_clinical_data.csv",
                          header=T)

PDL1_clinical.data <- subset(clinical.data, data_set != "VanAllen")
data_set <- unique(PDL1_clinical.data$data_set)



cluster.info <- lapply(data_set, function(x){
  per_info <- read.csv(paste0("result/last_shared_layer_clustering/Cui.NPJGenomMed.2016_",x,"/y_pred_2.csv"), row.names = 1)
  max_idx <- c(0,1)[which.max(table(per_info[,1]))]  
  per_info[,1] <- ifelse(per_info[,1] == max_idx, 0, 1)
  per_info
  })
names(cluster.info) <- data_set
lapply(cluster.info, table)
PDL1_clinical.data$per_embedding_cluster <- ifelse(do.call(rbind,cluster.info)[,1]==1, "Cluster_1", "Cluster_0")

combind.cluster.info <- read.csv(paste0("result/last_shared_layer_clustering/combind_skcm_datasets/y_pred_2.csv"), row.names = 1)
PDL1_clinical.data$embedding_cluster <- ifelse(combind.cluster.info[,1]==1, "Cluster_1", "Cluster_0")
#Cluster_0 Cluster_1 
#195        33 


library(survival)
library(survminer)
pdf(paste0(per_res_path, "KMcurve.pdf"))
ggsurvplot(survfit(Surv(OS, status) ~ embedding_cluster, data =PDL1_clinical.data),
           data=PDL1_clinical.data,legend = c(0.8,0.8),legend.title="",
           palette=c("#377eb8", "#e41a1c"),
           pval = TRUE, 
           xlab = "Time(months)",  
           ylab = "Survival Proportion", 
           title = paste0("Anti-PD1 therapy, Melanoma, N=", nrow(PDL1_clinical.data)))
for(cohort in data_set){
  plot.data <- subset(PDL1_clinical.data, data_set==cohort)
  P <- ggsurvplot(survfit(Surv(OS, status) ~ per_embedding_cluster, data =plot.data),
             data=plot.data,legend = c(0.8,0.8),legend.title="",
             palette=c("#377eb8", "#e41a1c"),
             pval = TRUE, 
             xlab = "Time(months)", 
             ylab = "Survival Proportion",
             title = paste0(cohort, ", N=", nrow(plot.data)))
  print(P)
} 

dev.off()
table(PDL1_clinical.data[,c("response", "embedding_cluster")])
#response Cluster_0 Cluster_1
#-1        30         4
#0        101         8
#1         64        21
# Multivariate cox regression analysis
outcome <- subset(PDL1_clinical.data, select = c("data_set","status", "OS", "embedding_cluster"))
colnames(outcome) <- c("cohort", "status","time", "embedding_cluster")
if(any(is.na(outcome)) | length(which(outcome$time==0))>0){
  outcome <- outcome[-c(which(!complete.cases(outcome)), which(outcome$time==0)),]
}
cox.res <- coxph(Surv(time, status) ~ cohort+ embedding_cluster, data =outcome)
summary(cox.res)
# 0.480(0.2477-0.9324) p value:0.0302

# Forest map visualization
pdf(paste0(per_res_path, "cox_forestplot.pdf"), height=4, width=4)
plot.data <- summary(cox.res)$coef[,c(2, 5)]
plot.data <- as.data.frame(cbind(plot.data, summary(cox.res)$conf.int[, c(3,4)]))
plot.data$Varnames <- factor(rownames(plot.data))
colnames(plot.data) <- c("HR", "p.value", "Lower", "Upper","Varnames")
p <- ggplot(plot.data, aes(HR, Varnames)) 
p + geom_point(size=3) +
  geom_errorbarh(aes(xmax =Upper, xmin = Lower), height = 0.1) +
  scale_x_continuous(limits= c(0, 2.2), breaks= seq(0, 2, 0.5)) +
  geom_vline(aes(xintercept = 1)) +
  xlab('Hazard Ratio') + ylab(' ') + theme_bw()
dev.off()


library(grid)
library(gridExtra)
library(plyr)
colors <- c("#1a9641","#d7191c")

plot.data <- PDL1_clinical.data[,c("response", "embedding_cluster")]
plot.data$response <- ifelse(plot.data$response==1,  "CR/PR", "SD/PD")
plot.data$embedding_cluster <- factor(plot.data$embedding_cluster)
pdf(paste0(per_res_path,"correlated_with_response.pdf"), width=4, height=4)
ggplot(plot.data,aes(embedding_cluster,fill=response))+
  geom_bar(stat="count",position="fill")+ scale_fill_manual(values= colors)+
  theme(axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        legend.title=element_blank(),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+ 
  labs(x="", y="Fraction of cases(%)", title="Melanoma")
dev.off()

chisq.test(plot.data$response, plot.data$embedding_cluster) #  p-value =  0.001418
#table(plot.data)
#response Cluster_0 Cluster_1
#CR/PR        64        21
#SD/PD       131        12
