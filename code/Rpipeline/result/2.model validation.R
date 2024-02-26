########################################################
#' 
#' 模型验证
##################################################
setwd("../../")


###  
#' Braun.NatMed.2021
################################
predict.score <- read.csv("result/model_validation/Braun.NatMed.2021/rescaled_final_model/val_preds.csv",
                          header=T, row.names=1)
clinical.data <- read.csv("data/validation_data/Braun.NatMed.2021/processed_data/sub_clinical_data.csv",
                          header=T)

######### Differences in predicted phenotypic scores were directly compared between predicted immune subtypes##################
#1) The phenotypic label values were converted into percentiles
predict_phenotype_scores_percentage <- apply(predict.score, 2, function(x){
  y = rep(0,length(x))
  for(i in 1:length(x)) {
    y[i] = round(length(which(x<=x[i]))/length(x)*100,4)}
  y
})
mean_percentage <- apply(predict_phenotype_scores_percentage, 1, mean)
predict_phenotype_scores_percentage <- cbind(predict_phenotype_scores_percentage, mean_percentage)

#2) difference in phenotypes'scores
imm_subtype <- clinical.data[,"ImmunoPhenotype"]
#Desert    Excluded Infiltrated       No_IF 
#21           3          79         205 
df <- as.data.frame(predict_phenotype_scores_percentage)
df$imm_subtype <- imm_subtype
df <- df[-which(imm_subtype=="No_IF"),]  # 103   6

p_I_D <- rep(0,5)
p_I_E <- rep(0,5)
for(i in 1:5){
  p_I_D[i] = wilcox.test(df[which(df$imm_subtype=="Infiltrated"), i], df[which(df$imm_subtype=="Desert"), i], alternative ="two.sided")$p.value
  p_I_E[i] = wilcox.test(df[which(df$imm_subtype=="Infiltrated"), i], df[which(df$imm_subtype=="Excluded"), i], alternative ="two.sided")$p.value
  }
p_I_D
# 0.049607094 0.044893632 0.007490784 0.570715524 0.018219994
p_I_E
# 0.44387534 0.10839971 0.02622255 0.18228909 0.05719136

# 3) visualization
res_path = "result/model_validation/Braun.NatMed.2021/"
library(ggplot2)
pdf(paste(res_path, 'ImmunoPhenotype_boxplot_val.pdf', sep=""), width=5, height=5)
for(i in 1:5){
  e <- ggplot(df, aes(x = factor(imm_subtype, levels=c("Infiltrated", "Excluded", "Desert")), y = df[,i], fill = imm_subtype))+  geom_boxplot()
  p <- e + geom_dotplot(binaxis = "y", stackdir = "center",position = position_dodge(1),dotsize = 0.5)+
    scale_fill_manual(values =  c("#4575b4","#66c2a5",  "#d53e4f")) +
    labs(x = "", y = colnames(df)[i])+
    ggtitle(colnames(df)[i])+
    theme_bw()
  print(p)
}
dev.off()

df <- as.data.frame(predict.score)
df$imm_subtype <- imm_subtype
df <- df[-which(imm_subtype=="No_IF"),]  # 103   6
pdf(paste(res_path, 'ImmunoPhenotype_predictScore_boxplot_val.pdf', sep=""), width=5, height=5)
for(i in 1:4){
  e <- ggplot(df, aes(x = factor(imm_subtype, levels=c("Infiltrated", "Excluded", "Desert")), y = df[,i], fill = imm_subtype))+  geom_boxplot()
  p <- e + geom_dotplot(binaxis = "y", stackdir = "center",position = position_dodge(1),dotsize = 0.5)+
    scale_fill_manual(values =  c("#4575b4","#66c2a5",  "#d53e4f")) +
    labs(x = "", y = colnames(df)[i])+
    ggtitle(colnames(df)[i])+
    theme_bw()
  print(p)
}
dev.off()

###  
#' Mariathasan.Nature.2018
################################
predict.score <- read.csv("result/model_validation/Mariathasan.Nature.2018/rescaled_final_model/val_preds.csv",
                          header=T, row.names=1)

clinical.data <- read.csv("data/validation_data/Mariathasan.Nature.2018/processed_data/sub_clinical_data.csv",
                          header=T)

######### Differences in predicted phenotypic scores were directly compared between predicted immune subtypes##################
#1) The phenotypic label values were converted into percentiles
predict_phenotype_scores_percentage <- apply(predict.score, 2, function(x){
  y = rep(0,length(x))
  for(i in 1:length(x)) {
    y[i] = round(length(which(x<=x[i]))/length(x)*100,4)}
  y
})
mean_percentage <- apply(predict_phenotype_scores_percentage, 1, mean)
predict_phenotype_scores_percentage <- cbind(predict_phenotype_scores_percentage, mean_percentage)

#2) difference in phenotypes'scores
imm_subtype <- clinical.data[,"Immune.phenotype"]
#desert excluded inflamed 
#76      134       74 
df <- as.data.frame(predict_phenotype_scores_percentage)
df$imm_subtype <- imm_subtype
df <- df[-which(is.na(imm_subtype)),]  # 284   6

p_I_D <- rep(0,5)
p_I_E <- rep(0,5)
for(i in 1:5){
  p_I_D[i] = wilcox.test(df[which(df$imm_subtype=="inflamed"), i], df[which(df$imm_subtype=="desert"), i], alternative ="two.sided")$p.value
  p_I_E[i] = wilcox.test(df[which(df$imm_subtype=="inflamed"), i], df[which(df$imm_subtype=="excluded"), i], alternative ="two.sided")$p.value
}
p_I_D
#4.130719e-19 1.338331e-11 3.440413e-21 7.297675e-18 1.096921e-20
p_I_E
# 3.372960e-17 7.661503e-05 1.649166e-16 8.296102e-14 6.585898e-16
# 3) visualization
res_path = "result/model_validation/Mariathasan.Nature.2018/"

library(ggplot2)
pdf(paste(res_path, 'ImmunoPhenotype_boxplot_val.pdf', sep=""), width=5, height=5)
for(i in 1:5){
  e <- ggplot(df, aes(x = factor(imm_subtype, levels=c("inflamed", "excluded", "desert")), y = df[,i], fill = imm_subtype))+  geom_boxplot()
  p <- e + geom_dotplot(binaxis = "y", stackdir = "center",position = position_dodge(1),dotsize = 0.5)+
    scale_fill_manual(values =  c("#4575b4","#66c2a5",  "#d53e4f")) +
    labs(x = "", y = colnames(df)[i])+
    ggtitle(colnames(df)[i])+
    theme_bw()
  print(p)
}
dev.off()

#3) Visualize different PDL1 levels
IC <- clinical.data[,"IC.Level"]
#IC0 IC1 IC2+
#97  132  118 
df <- as.data.frame(predict_phenotype_scores_percentage)
df$IC <- IC
df <- df[-which(is.na(IC)),]

# 计算显著性
p01 <- rep(0,5)
p12 <- rep(0,5)
p02 <- rep(0,5)
for(i in 1:5){
  p01[i] = wilcox.test(df[which(df$IC=="IC0"), i], df[which(df$IC=="IC1"), i], alternative ="two.sided")$p.value
  p12[i] = wilcox.test(df[which(df$IC=="IC1"), i], df[which(df$IC=="IC2+"), i], alternative ="two.sided")$p.value
  p02[i] = wilcox.test(df[which(df$IC=="IC0"), i], df[which(df$IC=="IC2+"), i], alternative ="two.sided")$p.value
}
p01
#7.422852e-08 1.157490e-05 5.700053e-13 1.993698e-09 1.649281e-11
p12
#  1.713506e-10 9.820543e-06 9.936788e-13 4.796590e-07 2.869549e-11
p02
#3.917907e-20 1.200005e-13 2.504570e-24 3.038498e-18 3.423624e-23


pdf(paste(res_path, 'IC_boxplot_val.pdf', sep=""), width=5, height=5)
for(i in 1:5){
  e <- ggplot(df, aes(x = factor(IC, levels=c("IC0", "IC1", "IC2+")), y = df[,i], fill = IC))+  
    geom_boxplot()
  p <- e + geom_dotplot(binaxis = "y", stackdir = "center",position = position_dodge(1),dotsize = 0.5)+
    scale_fill_manual(values =  c("#A7A7A7", "#AEDBD7", "#344B9E")) +
    labs(x = "", y = colnames(df)[i])+
    ggtitle(colnames(df)[i])+
    theme_bw()
  print(p)
}
dev.off()


TC <- clinical.data[,"TC.Level"]
#TC0  TC1 TC2+
#153   12   26
df <- as.data.frame(predict_phenotype_scores_percentage)
df$TC <- TC
df <- df[-which(is.na(TC)),]


p01 <- rep(0,5)
p12 <- rep(0,5)
p02 <- rep(0,5)
for(i in 1:5){
  p01[i] = wilcox.test(df[which(df$TC=="TC0"), i], df[which(df$TC=="TC1"), i], alternative ="two.sided")$p.value
  p12[i] = wilcox.test(df[which(df$TC=="TC1"), i], df[which(df$TC=="TC2+"), i], alternative ="two.sided")$p.value
  p02[i] = wilcox.test(df[which(df$TC=="TC0"), i], df[which(df$TC=="TC2+"), i], alternative ="two.sided")$p.value
  }
p01
#0.36859985 0.40685460 0.05097347 0.05913273 0.10245079
p12
#2.913901e-03 8.150875e-02 1.473622e-02 1.387518e-05 1.942254e-03
p02
#2.511401e-10 6.797212e-05 2.704623e-10 4.882089e-17 3.111755e-12

pdf(paste(res_path, 'TC_boxplot_val.pdf', sep=""), width=5, height=5)
for(i in 1:5){
  e <- ggplot(df, aes(x = factor(TC, levels=c("TC0", "TC1", "TC2+")), y = df[,i], fill = TC))+  
    geom_boxplot()
  p <- e + geom_dotplot(binaxis = "y", stackdir = "center",position = position_dodge(1),dotsize = 0.5)+
    scale_fill_manual(values =  c("#A7A7A7", "#AEDBD7", "#344B9E")) +
    labs(x = "", y = colnames(df)[i])+
    ggtitle(colnames(df)[i])+
    theme_bw()
  print(p)
}
dev.off()
