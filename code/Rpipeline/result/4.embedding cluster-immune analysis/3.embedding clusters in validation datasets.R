setwd("../../")
library(ggplot2)

res_path <- "result/embedding_clusters_immune_characterizing/val_imm_phenotype_diff/"
output_dir <- file.path(res_path)
if (!dir.exists(output_dir))  dir.create(output_dir)


###  
#' Braun.NatMed.2021
################################
clinical.data <- read.csv("data/validation_data/Braun.NatMed.2021/processed_data/sub_clinical_data.csv",
                          header=T)

cluster.info <- read.csv("result/last_shared_layer_clustering/Braun.NatMed.2021/y_pred_2.csv", row.names = 1)
#0   1 
#264  44
clinical.data$embedding_cluster <- cluster.info[,1]

pdf(paste0(res_path,"Braun.NatMed.2021_embedding_clustering.pdf"), width=4, height=4)
plot.data <- clinical.data[-which(clinical.data$ImmunoPhenotype=="No_IF"),c("ImmunoPhenotype", "embedding_cluster")]
plot.data$ImmunoPhenotype <- factor(plot.data$ImmunoPhenotype, levels = c("Infiltrated", "Excluded", "Desert"))
plot.data$embedding_cluster <- factor(plot.data$embedding_cluster)
colors <- c("#d53e4f","#66c2a5","#4575b4")
p1 <- ggplot(plot.data,aes(embedding_cluster,fill=ImmunoPhenotype))+
  geom_bar(stat="count",position="fill")+ scale_fill_manual(values= colors)+
  theme(axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+ 
  labs(x="", y="Fraction of cases(%)", title="Immune.Subtype")
print(p1)
dev.off()
chisq.test(plot.data$ImmunoPhenotype, plot.data$embedding_cluster,simulate.p.value = FALSE ) #  p-value =0.6916 
chisq.test(plot.data$ImmunoPhenotype, plot.data$embedding_cluster,simulate.p.value = TRUE ) #  p-value =0.8116 
table(plot.data)
#ImmunoPhenotype  0  1
#Infiltrated 68 11
#Excluded     3  0
#Desert      19  2

###  
#' Mariathasan.Nature.2018
################################
clinical.data <- read.csv("data/validation_data/Mariathasan.Nature.2018/processed_data/sub_clinical_data.csv",
                          header=T)

cluster.info <- read.csv("result/last_shared_layer_clustering/Mariathasan.Nature.2018/y_pred_2.csv", row.names = 1)
clinical.data$embedding_cluster <- cluster.info[,1]
clinical.data$embedding_cluster[which(clinical.data$embedding_cluster==0)] <- "cluster_1"
clinical.data$embedding_cluster[which(clinical.data$embedding_cluster==1)] <- "cluster_0"

pdf(paste0(res_path,"Mariathasan.Nature.2018_embedding_clustering.pdf"), width=4, height=4)
plot.data <- clinical.data[,c("Immune.phenotype","IC.Level", "TC.Level","embedding_cluster")]
plot.data$Immune.phenotype <- factor(plot.data$Immune.phenotype, levels = c("inflamed", "excluded", "desert"))
plot.data$embedding_cluster <- factor(plot.data$embedding_cluster)
colors <- c("#d53e4f","#66c2a5","#4575b4")
p1 <- ggplot(subset(plot.data, Immune.phenotype != "NA"),aes(embedding_cluster,fill=Immune.phenotype))+
  geom_bar(stat="count",position="fill")+ scale_fill_manual(values= colors)+
  theme(axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+
  labs(x="", y="Fraction of cases(%)", title="Immune.Subtype")
print(p1)

colors <- c("#A7A7A7", "#AEDBD7", "#344B9E")
p2 <- ggplot(subset(plot.data, IC.Level != "NA"),aes(embedding_cluster,fill=IC.Level))+
  geom_bar(stat="count",position="fill")+ scale_fill_manual(values= colors)+
  theme(axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+
  labs(x="", y="Fraction of cases(%)", title="Immune.Subtype")
print(p2)
p3 <- ggplot(subset(plot.data, TC.Level != "NA"),aes(embedding_cluster,fill=TC.Level))+
  geom_bar(stat="count",position="fill")+ scale_fill_manual(values= colors)+
  theme(axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+
  labs(x="", y="Fraction of cases(%)", title="Immune.Subtype")
print(p3)
dev.off()

chisq.test(plot.data$Immune.phenotype, plot.data$embedding_cluster) #  p-value = 7.129e-13
chisq.test(plot.data$IC.Level, plot.data$embedding_cluster) #  p-value =  1.675e-10
chisq.test(plot.data$TC.Level, plot.data$embedding_cluster,simulate.p.value = FALSE) # p-value = 3.433e-06 
chisq.test(plot.data$TC.Level, plot.data$embedding_cluster,simulate.p.value = TRUE) #  p-value = 0.0004998

table(subset(plot.data, Immune.phenotype != "NA", select=c(Immune.phenotype, embedding_cluster)))
#Immune.phenotype cluster_0 cluster_1
#inflamed        45        29
#excluded       125         9
#desert          75         1
table(subset(plot.data, IC.Level != "NA", select=c(IC.Level, embedding_cluster)))
#IC.Level cluster_0 cluster_1
#IC0         95         2
#IC1        123         9
#IC2+        82        36
table(subset(plot.data, TC.Level != "NA", select=c(TC.Level, embedding_cluster)))
####
#TC.Level cluster_0 cluster_1
#TC0        248        27
#TC1         20         2
#TC2+        32        18
