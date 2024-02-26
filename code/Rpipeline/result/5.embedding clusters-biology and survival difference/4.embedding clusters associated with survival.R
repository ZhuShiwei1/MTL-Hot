
setwd("../../")


load("data/TCGA_data/TCGA_clinical_data.RData")
cluster.info <- read.csv("result/last_shared_layer_clustering/TCGA/0/y_pred_2.csv", row.names = 1)
TCGA_clinical_data$embedding_cluster <- cluster.info[,1]
TCGA_clinical_data$embedding_cluster[which(TCGA_clinical_data$embedding_cluster==0)] <- "Cluster_0"
TCGA_clinical_data$embedding_cluster[which(TCGA_clinical_data$embedding_cluster==1)] <- "Cluster_1"
TCGA_clinical_data$embedding_cluster <- factor(TCGA_clinical_data$embedding_cluster)
#Cluster_0 Cluster_1 
#7032      1723 
TCGA_clinical_data$cancer_type[which(TCGA_clinical_data$cancer_type %in% c("COAD", "READ"))] <- "COREAD"
res_path <- "result/embedding_clusters_biology_survival_diff/diff_in_survival/"
output_dir <- file.path(res_path)
if (!dir.exists(output_dir))  dir.create(output_dir)




###  
#' 
#' Compare the survival differences between these two clusters
################################
library(survival)
library(survminer)
library(ggpubr)
library(forestplot)
### OS

fit.data <- TCGA_clinical_data[complete.cases(TCGA_clinical_data[,c("OS", "OS.time","embedding_cluster", "age_at_initial_pathologic_diagnosis","gender")]),]
cox.res <- coxph(Surv(OS.time,OS) ~ embedding_cluster+age_at_initial_pathologic_diagnosis+gender+cancer_type, data = fit.data)
summary(cox.res) 
#HR = exp(coef) = 0.805317(0.72894-0.8933),P=3.53e-05
ggforest(cox.res, data = fit.data)

pdf(paste0(res_path,"imm_corr_OS_forestplot.pdf"), height=5, width=5)
# visualization
HR <- signif(exp(coef(cox.res))[1], 3)
HR.confint.lower <- signif(exp(confint(cox.res))[1,1], 3)
HR.confint.upper <- signif(exp(confint(cox.res))[1,2], 3)
test_data <- data.frame(coef = c( NA, 1, HR),
                        low = c(NA, 1, HR.confint.lower),
                        high = c(NA, 1, HR.confint.upper))
tabletext <- cbind(c("\n",paste0("Cluster_0(", table(fit.data$embedding_cluster)[1],")"), 
                     paste0("Cluster_1(", table(fit.data$embedding_cluster)[2],")")),
                   c("HR (95% CI for HR)", "reference", paste0(HR,"(", HR.confint.lower, "-", HR.confint.upper, ")")),
                   c("P Value", NA, signif(summary(cox.res)$coef[1,5], digits=3)))
forestplot(labeltext = tabletext,graph.pos=2,
           mean = test_data$coef,lower = test_data$low,upper = test_data$high,
           zero = 1, colgap=unit(5,"mm"),
           title="Hazard Ratio Plot", clip = c(0.5, 1.5), 
           txt_gp=fpTxtGp(label=gpar(cex=0.5),
                          ticks=gpar(cex=0.5),
                          xlab=gpar(cex = 0.5),
                          title=gpar(cex = 0.5)),
           xticks = seq(0.6, 1.2, 0.2),
           col = fpColors(box = "#BC3C28", lines="#BC3C28",zero = "gray50"),
           lineheight = "auto",
           xlab = "HR")
dev.off()


pdf(paste0(res_path,"imm_corr_OS_KMcurve.pdf"))
# KM curves for all samples
nsamples <- nrow(fit.data)
fit <- survfit(Surv(OS.time,OS) ~ embedding_cluster, data = fit.data) 
p <- ggsurvplot(fit, data = fit.data, legend = c(0.8,0.8),legend.title="",
                palette=c("#377eb8", "#e41a1c"),
                title= paste("embedding_cluster, N=", nsamples, sep=""),
                xlab = "Follow up time(d)", ylab = "Survival probability", pval = TRUE)
print(p)
# KM curves for each cancers
cluster_count <- as.data.frame(table(subset(TCGA_clinical_data, select=c(cancer_type, embedding_cluster))))
cancers <- unique(cluster_count[(cluster_count[,3]>5 & cluster_count[,2]=="Cluster_1"),1])
for(cancer in cancers){
  data <- subset(TCGA_clinical_data, cancer_type==cancer, select=c("OS", "OS.time", "embedding_cluster"))
  colnames(data)[1:2] <- c("status", "time")
  
  fit.data <- data[complete.cases(data[,c("status", "time","embedding_cluster")]),]
  fit <- survfit(Surv(time,status) ~ embedding_cluster, data = fit.data)
  p <- ggsurvplot(fit, data = fit.data, legend = c(0.8,0.8),legend.title="",
                  palette=c("#377eb8", "#e41a1c"),
                  title= paste(cancer,": embedding_cluster, N=", nrow(fit.data), sep=""),
                  xlab = "Follow up time(d)", ylab = "Survival probability", pval = TRUE)
  print(p)
}
dev.off()






###PFI
fit.data <- TCGA_clinical_data[complete.cases(TCGA_clinical_data[,c("PFI", "PFI.time","embedding_cluster","age_at_initial_pathologic_diagnosis")]),]

cox.res <- coxph(Surv(PFI.time,PFI) ~ embedding_cluster+age_at_initial_pathologic_diagnosis+cancer_type, data = fit.data)# 创建生存对象 
summary(cox.res) 
#HR = exp(coef) = 0.837997(0.7594-0.9248),P=0.000439

pdf(paste0(res_path,"imm_corr_PFI_forestplot.pdf"), height=5, width=5)

HR <- signif(exp(coef(cox.res))[1], 3)
HR.confint.lower <- signif(exp(confint(cox.res))[1,1], 3)
HR.confint.upper <- signif(exp(confint(cox.res))[1,2], 3)
test_data <- data.frame(coef = c( NA, 1, HR),
                        low = c(NA, 1, HR.confint.lower),
                        high = c(NA, 1, HR.confint.upper))
tabletext <- cbind(c("\n",paste0("Cluster_0(", table(fit.data$embedding_cluster)[1],")"), 
                     paste0("Cluster_1(", table(fit.data$embedding_cluster)[2],")")),
                   c("HR (95% CI for HR)", "reference", paste0(HR,"(", HR.confint.lower, "-", HR.confint.upper, ")")),
                   c("P Value", NA, signif(summary(cox.res)$coef[1,5], digits=3)))
forestplot(labeltext = tabletext,graph.pos=2,
           mean = test_data$coef,lower = test_data$low,upper = test_data$high,
           zero = 1, colgap=unit(5,"mm"),
           title="Hazard Ratio Plot", clip = c(0.5, 1.5), 
           txt_gp=fpTxtGp(label=gpar(cex=0.5),
                          ticks=gpar(cex=0.5),
                          xlab=gpar(cex = 0.5),
                          title=gpar(cex = 0.5)),
           xticks = seq(0.6, 1.2, 0.2),
           col = fpColors(box = "#BC3C28", lines="#BC3C28",zero = "gray50"),
           lineheight = "auto",
           xlab = "HR")
dev.off()


pdf(paste0(res_path,"imm_corr_PFI_KMcurve_temp.pdf"))

nsamples <- nrow(fit.data)
fit <- survfit(Surv(PFI.time,PFI) ~ embedding_cluster, data = fit.data)
p <- ggsurvplot(fit, data = fit.data, legend = c(0.8,0.8),legend.title="",
                palette=c("#377eb8", "#e41a1c"),
                title= paste("embedding_cluster, N=", nsamples, sep=""),
                xlab = "Follow up time(d)", ylab = "PFI probability", pval = TRUE)
print(p)

cluster_count <- as.data.frame(table(subset(TCGA_clinical_data, select=c(cancer_type, embedding_cluster))))
cancers <- unique(cluster_count[(cluster_count[,3]>5 & cluster_count[,2]=="Cluster_1"),1])
for(cancer in cancers){
  data <- subset(TCGA_clinical_data, cancer_type==cancer, select=c("PFI", "PFI.time", "embedding_cluster"))
  colnames(data)[1:2] <- c("status", "time")
 
  fit.data <- data[complete.cases(data[,c("status", "time","embedding_cluster")]),]
  fit <- survfit(Surv(time,status) ~ embedding_cluster, data = fit.data)
  p <- ggsurvplot(fit, data = fit.data, legend = c(0.8,0.8),legend.title="",
                  palette=c("#377eb8", "#e41a1c"),
                  title= paste(cancer,": embedding_cluster, N=", nrow(fit.data), sep=""),
                  xlab = "Follow up time(d)", ylab = "PFI probability", pval = TRUE)
  print(p)
}
dev.off()



###  
#' 
#' Independence of prognostic prediction in specific cancer types
################################
cluster_count <- as.data.frame(table(subset(TCGA_clinical_data, select=c(cancer_type, embedding_cluster))))
cancers <- unique(cluster_count[(cluster_count[,3]>5 & cluster_count[,2]=="Cluster_1"),1])
colnames(TCGA_clinical_data)[3] <- "Age"
colnames(TCGA_clinical_data)[49] <- "Clusters"
# BLCA, HNSC, SKCM
cancer <- "SKCM"
sub_clinical_data <- subset(TCGA_clinical_data, cancer_type==cancer)
fit.data <- sub_clinical_data[complete.cases(sub_clinical_data[,c("gender","Age" , "OS", "OS.time","Clusters")]),]
cox.res <- coxph(Surv(OS.time,OS) ~ Clusters+gender+Age, data = fit.data) 
summary(cox.res) 
ggforest(cox.res, data = fit.data)

pdf(paste0(res_path, cancer, "_cox_forestplot.pdf"), height=4, width=4)
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


cancer <- "BLCA"
sub_clinical_data <- subset(TCGA_clinical_data, cancer_type==cancer)
fit.data <- sub_clinical_data[complete.cases(sub_clinical_data[,c("gender","Smoking","Age","P_stage","OS", "OS.time","Clusters")]),]
cox.res <- coxph(Surv(OS.time,OS) ~ Clusters+gender+Age+Smoking+P_stage, data = fit.data) 
summary(cox.res)
ggforest(cox.res, data = fit.data)
pdf(paste0(res_path, cancer, "_cox_forestplot.pdf"), height=4, width=4)
plot.data <- summary(cox.res)$coef[,c(2, 5)]
plot.data <- as.data.frame(cbind(plot.data, summary(cox.res)$conf.int[, c(3,4)]))
plot.data$Varnames <- factor(rownames(plot.data))
colnames(plot.data) <- c("HR", "p.value", "Lower", "Upper","Varnames")
p <- ggplot(plot.data, aes(HR, Varnames)) 
p + geom_point(size=3) +
  geom_errorbarh(aes(xmax =Upper, xmin = Lower), height = 0.1) +
  scale_x_continuous(limits= c(0, 3.2), breaks= seq(0, 3.2, 0.5)) +
  geom_vline(aes(xintercept = 1)) +
  xlab('Hazard Ratio') + ylab(' ') + theme_bw()
dev.off()


cancer <- "HNSC"
sub_clinical_data <- subset(TCGA_clinical_data, cancer_type==cancer)
fit.data <- sub_clinical_data[complete.cases(sub_clinical_data[,c("gender","Smoking","Age","HPV","P_stage","OS", "OS.time","Clusters")]),]
cox.res <- coxph(Surv(OS.time,OS) ~ Clusters+gender+Age+Smoking+HPV+P_stage, data = fit.data)
summary(cox.res) 
ggforest(cox.res, data = fit.data)
pdf(paste0(res_path, cancer, "_cox_forestplot.pdf"), height=4, width=4)
plot.data <- summary(cox.res)$coef[,c(2, 5)]
plot.data <- as.data.frame(cbind(plot.data, summary(cox.res)$conf.int[, c(3,4)]))
plot.data$Varnames <- factor(rownames(plot.data))
colnames(plot.data) <- c("HR", "p.value", "Lower", "Upper","Varnames")
p <- ggplot(plot.data, aes(HR, Varnames)) 
p + geom_point(size=3) +
  geom_errorbarh(aes(xmax =Upper, xmin = Lower), height = 0.1) +
  scale_x_continuous(limits= c(0, max(plot.data$Upper)), breaks= seq(0, max(plot.data$Upper), 0.5)) +
  geom_vline(aes(xintercept = 1)) +
  xlab('Hazard Ratio') + ylab(' ') + theme_bw()
dev.off()

#### PFI
cancer <- "SKCM"
sub_clinical_data <- subset(TCGA_clinical_data, cancer_type==cancer)
fit.data <- sub_clinical_data[complete.cases(sub_clinical_data[,c("gender","Age" , "PFI", "PFI.time","Clusters")]),]
cox.res <- coxph(Surv(PFI.time,PFI) ~ Clusters+gender+Age, data = fit.data) 
summary(cox.res) 
ggforest(cox.res, data = fit.data)

pdf(paste0(res_path, cancer, "_cox_PFI_forestplot.pdf"), height=4, width=4)
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


cancer <- "BLCA"
sub_clinical_data <- subset(TCGA_clinical_data, cancer_type==cancer)
fit.data <- sub_clinical_data[complete.cases(sub_clinical_data[,c("gender","Smoking","Age","P_stage","PFI", "PFI.time","Clusters")]),]
cox.res <- coxph(Surv(PFI.time,PFI) ~ Clusters+gender+Age+Smoking+P_stage, data = fit.data)
summary(cox.res)
ggforest(cox.res, data = fit.data)
pdf(paste0(res_path, cancer, "_cox_PFI_forestplot.pdf"), height=4, width=4)
plot.data <- summary(cox.res)$coef[,c(2, 5)]
plot.data <- as.data.frame(cbind(plot.data, summary(cox.res)$conf.int[, c(3,4)]))
plot.data$Varnames <- factor(rownames(plot.data))
colnames(plot.data) <- c("HR", "p.value", "Lower", "Upper","Varnames")
p <- ggplot(plot.data, aes(HR, Varnames)) 
p + geom_point(size=3) +
  geom_errorbarh(aes(xmax =Upper, xmin = Lower), height = 0.1) +
  scale_x_continuous(limits= c(0, 3.2), breaks= seq(0, 3.2, 0.5)) +
  geom_vline(aes(xintercept = 1)) +
  xlab('Hazard Ratio') + ylab(' ') + theme_bw()
dev.off()


cancer <- "HNSC"
sub_clinical_data <- subset(TCGA_clinical_data, cancer_type==cancer)
fit.data <- sub_clinical_data[complete.cases(sub_clinical_data[,c("gender","Smoking","Age","HPV","P_stage","PFI", "PFI.time","Clusters")]),]
cox.res <- coxph(Surv(PFI.time,PFI) ~ Clusters+gender+Age+Smoking+HPV+P_stage, data = fit.data) 
summary(cox.res)
ggforest(cox.res, data = fit.data)
pdf(paste0(res_path, cancer, "_cox_PFI_forestplot.pdf"), height=4, width=4)
plot.data <- summary(cox.res)$coef[,c(2, 5)]
plot.data <- as.data.frame(cbind(plot.data, summary(cox.res)$conf.int[, c(3,4)]))
plot.data$Varnames <- factor(rownames(plot.data))
colnames(plot.data) <- c("HR", "p.value", "Lower", "Upper","Varnames")
p <- ggplot(plot.data, aes(HR, Varnames)) 
p + geom_point(size=3) +
  geom_errorbarh(aes(xmax =Upper, xmin = Lower), height = 0.1) +
  scale_x_continuous(limits= c(0, max(plot.data$Upper)), breaks= seq(0, max(plot.data$Upper), 0.5)) +
  geom_vline(aes(xintercept = 1)) +
  xlab('Hazard Ratio') + ylab(' ') + theme_bw()
dev.off()

