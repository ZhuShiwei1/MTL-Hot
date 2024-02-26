
setwd("../../")
function.files <- list.files(path="code/Rpipeline/function/", recursive=TRUE, pattern="*.R")
for(i in function.files) source(paste("code/Rpipeline/function/", i, sep=""))
res_path = "result/embedding_clusters_immunotherapy_response/Mariathasan.Nature.2018/"

###  
#' TGF-β attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells
#' Mariathasan.Nature.2018
################################

load("data/validation_data/Mariathasan.Nature.2018/processed_data/cancer.phenotype.RData")



clinical.data <- read.csv("data/validation_data/Mariathasan.Nature.2018/processed_data/sub_clinical_data.csv",
                          header=T)


cluster.info <- read.csv("result/last_shared_layer_clustering/Mariathasan.Nature.2018/y_pred_2.csv", row.names = 1)
clinical.data$embedding_cluster <- ifelse(cluster.info[,1]==0, "Cluster_1", "Cluster_0")
clinical.data <- cbind(clinical.data, cancer.phenotype)

for(pheno in colnames(cancer.phenotype)){
  clinical.data[,pheno] <- ifelse(clinical.data[,pheno]>median(clinical.data[,pheno]), "High", "Low")
}

clinical.data$good_markers_counts <- apply(clinical.data[,colnames(cancer.phenotype)], 1, function(x){
  if(all(x=="High")){
    res="group A"
  }else if(length(which(x=="High"))%in% c(2,3)){
    res="group B"
  }else if(length(which(x=="High"))%in% c(0,1)){
    res="group C"  
  }
  res
})
table(clinical.data[which(clinical.data$embedding_cluster=="Cluster_1"),"good_markers_counts" ])
#group A group B 
#45       3 


library(survival)
library(survminer)
pdf(paste0(res_path, "phenotype_score_KMcurve1.pdf"))
outcome <- subset(clinical.data, select = c("censOS", "os", "embedding_cluster"))
colnames(outcome) <- c("status","time", "embedding_cluster")
if(any(is.na(outcome)) | length(which(outcome$time==0))>0){
  outcome <- outcome[-c(which(!complete.cases(outcome)), which(outcome$time==0)),]
}

for(pheno in colnames(cancer.phenotype)){
  pheno_score <- clinical.data[,pheno]
  outcome1 <- cbind(outcome, pheno_score)
  outcome1$pheno_score <- factor(outcome1$pheno_score, levels=c("Low", "High"))
  p <- ggsurvplot(survfit(Surv(time, status) ~ pheno_score, data =outcome1),
             data=outcome1,legend = c(0.8,0.8),legend.title="",
             palette=c("#377eb8", "#e41a1c"),
             pval = TRUE, 
             xlab = "Time(months)",  
             ylab = "Survival Proportion", 
             title = pheno)
  print(p)
}
pheno_score <- clinical.data[,"good_markers_counts"]
outcome1 <- cbind(outcome, pheno_score)
p <- ggsurvplot(survfit(Surv(time, status) ~ pheno_score, data =outcome1),
                data=outcome1,legend = c(0.8,0.8),legend.title="",
                pval = TRUE,
                xlab = "Time(months)",  
                ylab = "Survival Proportion", 
                title = "")
print(p)
dev.off()


#######
# 
#'  Univariate cox regression results
#######
pdf(paste0(res_path, "phenotype_score_frestplot.pdf"),  height=1, width=3)
outcome <- subset(clinical.data, select = c("censOS", "os", "embedding_cluster"))
colnames(outcome) <- c("status","time", "embedding_cluster")
if(any(is.na(outcome)) | length(which(outcome$time==0))>0){
  outcome <- outcome[-c(which(!complete.cases(outcome)), which(outcome$time==0)),]
}
#  Univariate cox regression results for embedding clusters
cox.res <- coxph(Surv(time, status) ~ embedding_cluster, data =outcome)
summary(cox.res)
# 0.5229(0.3254-0.8403) p=0.00739 

plot.data <- summary(cox.res)$coef[,c(2, 5)]
plot.data <- as.data.frame(c(plot.data, summary(cox.res)$conf.int[, c(3,4)]))
plot.data <- as.data.frame(t(plot.data))
colnames(plot.data) <- c("HR", "p.value", "Lower", "Upper")
plot.data$Varnames <- "embedding_clusterCluster_1"
p <- ggplot(plot.data, aes(HR, Varnames)) 
p <- p + geom_point(size=3) +
  geom_errorbarh(aes(xmax =Upper, xmin = Lower), height = 0.1) +
  scale_x_continuous(limits= c(0, 2.5), breaks= seq(0, 2.5, 0.5)) +
  geom_vline(aes(xintercept = 1)) +
  theme(axis.text.y = element_blank(),
        panel.background = element_rect(fill = "white"), 
        panel.border = element_rect(color = "black",fill=NA))+
  xlab('Hazard Ratio') + ylab(' ')
print(p)

#  Univariate cox regression results for each phenotype
for(pheno in colnames(cancer.phenotype)){
  pheno_score <- clinical.data[,pheno]
  outcome1 <- cbind(outcome, pheno_score)
  outcome1$pheno_score <- factor(outcome1$pheno_score, levels=c("Low", "High"))
  cox.res <- coxph(Surv(time, status) ~ pheno_score, data =outcome1)
  print(summary(cox.res))
  plot.data <- summary(cox.res)$coef[,c(2, 5)]
  plot.data <- as.data.frame(c(plot.data, summary(cox.res)$conf.int[, c(3,4)]))
  plot.data <- as.data.frame(t(plot.data))
  colnames(plot.data) <- c("HR", "p.value", "Lower", "Upper")
  plot.data$Varnames <- pheno
  p <- ggplot(plot.data, aes(HR, Varnames))
  p <- p + geom_point(size=3) +
    geom_errorbarh(aes(xmax =Upper, xmin = Lower), height = 0.1) +
    scale_x_continuous(limits= c(0, 2.5), breaks= seq(0, 2.5, 0.5)) +
    geom_vline(aes(xintercept = 1)) +
    theme(axis.text.y = element_blank(),
          panel.background = element_rect(fill = "white"), 
          panel.border = element_rect(color = "black",fill=NA))+
    xlab('Hazard Ratio') + ylab(' ') 
  print(p)
}
dev.off()

# Univariate cox regression results for good markers
pdf(paste0(res_path, "good_markers_counts_frestplot1.pdf"),  height=2, width=3)
outcome <- subset(clinical.data, select = c("censOS", "os", "good_markers_counts"))
colnames(outcome) <- c("status","time", "good_markers_counts")
#outcome$good_markers_counts <- factor(outcome$good_markers_counts, levels=c("group E","group D", "group C", "group B", "group A"))
outcome$good_markers_counts <- factor(outcome$good_markers_counts, levels=c("group C", "group B", "group A"))

if(any(is.na(outcome)) | length(which(outcome$time==0))>0){
  outcome <- outcome[-c(which(!complete.cases(outcome)), which(outcome$time==0)),]
}
cox.res <- coxph(Surv(time, status) ~ good_markers_counts, data =outcome)
summary(cox.res)
plot.data <- summary(cox.res)$coef[,c(2, 5)]
plot.data <- as.data.frame(cbind(plot.data, summary(cox.res)$conf.int[, c(3,4)]))
colnames(plot.data) <- c("HR", "p.value", "Lower", "Upper")
plot.data$Varnames <- rownames(plot.data)
p <- ggplot(plot.data, aes(HR, Varnames))
p <- p + geom_point(size=3) +
  geom_errorbarh(aes(xmax =Upper, xmin = Lower), height = 0.1) +
  scale_x_continuous(limits= c(0, 2), breaks= seq(0, 2, 0.5)) +
  geom_vline(aes(xintercept = 1)) +
  theme(axis.text.y = element_blank(),
        panel.background = element_rect(fill = "white"), 
        panel.border = element_rect(color = "black",fill=NA))+
  xlab('Hazard Ratio') + ylab(' ')
print(p)

dev.off()



library(grid)
library(gridExtra)
library(plyr)
colors <- c("#1a9641","#a6d96a","#fdae61","#d7191c")

plot.data <- clinical.data[,c("Best.Confirmed.Overall.Response", "good_markers_counts", colnames(cancer.phenotype))]
plot.data <- plot.data[which(plot.data$Best.Confirmed.Overall.Response %in% c("PR", "CR", "PD", "SD")),]
plot.data$Best.Confirmed.Overall.Response <- factor(plot.data$Best.Confirmed.Overall.Response, levels = c("CR","PR","SD", "PD"))
plot.data$good_markers_counts <- factor(plot.data$good_markers_counts)
pdf(paste0(res_path,"good_markers_correlated_with_response.pdf"), width=4, height=4)
p <- ggplot(plot.data,aes(good_markers_counts,fill=Best.Confirmed.Overall.Response))+
  geom_bar(stat="count",position="fill")+ scale_fill_manual(values= colors)+
  theme(axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        legend.title=element_blank(),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+ 
  labs(x="", y="Fraction of cases(%)", title="BLCA")
print(p)
for(pheno in colnames(cancer.phenotype)){
  p <- ggplot(plot.data,aes(plot.data[,pheno],fill=Best.Confirmed.Overall.Response))+
    geom_bar(stat="count",position="fill")+ scale_fill_manual(values= colors)+
    theme(axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
          legend.title=element_blank(),
          panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank())+ 
    labs(x="", y="Fraction of cases(%)", title=pheno)
  print(p)
}
dev.off()

table(plot.data)
#Best.Confirmed.Overall.Response group A group B group C
#CR      11       9       5
#PR      13      10      20
#SD      22      17      24
#PD      40      45      82

chisq.test(plot.data$Best.Confirmed.Overall.Response, plot.data$good_markers_counts) #  p-value = 0.1056 #  0.1247
chisq.test(plot.data$Best.Confirmed.Overall.Response, plot.data$good_markers_counts, simulate.p.value=T) #  p-value =  0.09895  # 0.1364
chisq.test(plot.data$Best.Confirmed.Overall.Response, plot.data$APM.score) #  p-value = 0.1316
chisq.test(plot.data$Best.Confirmed.Overall.Response, plot.data$Tcell.score) # p-value = 0.3107
chisq.test(plot.data$Best.Confirmed.Overall.Response, plot.data$IFNgamma.score) # p-value = 0.02429
chisq.test(plot.data$Best.Confirmed.Overall.Response, plot.data$PDL1.exp) # p-value = 0.233
