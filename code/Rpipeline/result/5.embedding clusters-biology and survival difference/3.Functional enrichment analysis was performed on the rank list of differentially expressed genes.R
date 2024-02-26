
setwd("../../")
source("code/Rpipeline/function/get_pathway_enrich.R")
library(enrichplot)


res_path <- "result/embedding_clusters_biology_survival_diff/RankedGeneList_enriched_pathways/"
output_dir <- file.path(res_path)
if (!dir.exists(output_dir))  dir.create(output_dir)

load("result/embedding_clusters_biology_survival_diff/diff_expression_analysis/diff_exp_genes_df.RData")
logFC <- diff_exp_genes_df$logFC
names(logFC) <- rownames(diff_exp_genes_df)
genelist <- logFC[order(logFC, decreasing = T)]

# load KEGG classification(Note: 2023-9-15 Data download athttps://www.genome.jp/brite/hsa00001.keg）
KEGG_classification <- read.table("data/pathway/KEGG-pathway-classification.txt",
                                  header =T, sep="\t", as.is=TRUE,
                                  colClasses= "character")
# import Reactome pathways (（2023-09-01）https://reactome.org/download/current/)
homo_Reactome_pathways <- read.table("data/pathway/homo_Reactome_pathways.txt",
                                     header =T, sep="\t",stringsAsFactors=FALSE,quote = "")


###  
#' 
#' GOs of enrichment of rank gene list
################################
set.seed(1234)
gsetSource <- "GO"
GO_gse_res <- get_pathway_enrich(genelist=genelist, gsetSource=gsetSource, keyType='SYMBOL',
                                 method="gse", pvalueCutoff=0.01) 
length(GO_gse_res$ID) #   1271


gsetSource <- "KEGG"
KEGG_gse_res <- get_pathway_enrich(genelist=genelist, gsetSource=gsetSource, keyType='SYMBOL',
                                   method="gse", pvalueCutoff=0.01) 
length(KEGG_gse_res$ID) #  92

disease_pathway <- paste0("hsa",KEGG_classification$Pathway.ID[which(KEGG_classification$Pathway.Class.1=="Human Diseases")])
select_KEGG_gse_res <- KEGG_gse_res@result[-which(KEGG_gse_res@result$ID %in% disease_pathway),] # 46


gsetSource <- "Reactome"
Reactome_gse_res <- get_pathway_enrich(genelist=genelist, gsetSource=gsetSource, keyType='SYMBOL',
                                       method="gse", pvalueCutoff=0.01) 
length(Reactome_gse_res$ID) #  212

pancancer_gse_res <- list(GO=GO_gse_res@result, KEGG=select_KEGG_gse_res, Reactome=Reactome_gse_res@result)
save(pancancer_gse_res, file=paste0(res_path,"pancancer_gse_res.RData"))

outputs <- rbind(rbind(pancancer_gse_res$KEGG, pancancer_gse_res$Reactome), pancancer_gse_res$GO[,-1])
outputs$Pathway_source <- c(rep("KEGG", nrow(pancancer_gse_res$KEGG)),
                            rep("Reactome", nrow(pancancer_gse_res$Reactome)),
                            rep("GO", nrow(pancancer_gse_res$GO)))
outputs <- outputs[,c(12, 1:11)]
write.csv(outputs, file=paste0(res_path, "pancancer_gse_res.csv"), quote =F, row.names=F)
lapply(pancancer_gse_res, function(df)table(df$NES>0))
#$GO
#FALSE  TRUE 
#117  1154 
#$KEGG
#FALSE  TRUE 
#10    36 
#$Reactome
#FALSE  TRUE 
#42   170 
lapply(pancancer_gse_res, function(df)table(abs(df$NES)>1))
#$GO
#TRUE 
#1271 
#$KEGG
#TRUE 
#46 
#$Reactome
#TRUE 
#212 


###  
#' 
#' KEGG and Reactome pathways of enrichment of rank gene list
################################

pdf(paste0(res_path, "go_gse_top20_barplot.pdf"), width=6, height=4)
bar <- get_pathway_barplot(gse_res=subset(pancancer_gse_res$GO, NES<0), gsetSource="GO", 
                           top=20)
print(bar)
bar <- get_pathway_barplot(gse_res=subset(pancancer_gse_res$GO, NES>0), gsetSource="GO", 
                           top=20)
print(bar)
dev.off()


# Barplots visualizes enrichment pathway results
pdf(paste0(res_path, "KEGG_gse_top20_barplot.pdf"), width=6, height=4)
bar <- get_pathway_barplot(gse_res=subset(pancancer_gse_res$KEGG, NES<0), gsetSource="KEGG", 
                           extra_data=KEGG_classification, top=20)
print(bar)
bar <- get_pathway_barplot(gse_res=subset(pancancer_gse_res$KEGG, NES>0), gsetSource="KEGG", 
                           extra_data=KEGG_classification, top=20)
print(bar)
dev.off()
pdf(paste0(res_path, "KEGG_gse_barplot.pdf"))
# 10
bar <- get_pathway_barplot(gse_res=subset(pancancer_gse_res$KEGG, NES<0), gsetSource="KEGG", 
                           extra_data=KEGG_classification)
print(bar)
# 36
bar <- get_pathway_barplot(gse_res=subset(pancancer_gse_res$KEGG, NES>0), gsetSource="KEGG", 
                           extra_data=KEGG_classification)
print(bar)
dev.off()

### Reactome
pdf(paste0(res_path, "Reactome_gse_top20_barplot.pdf"), width=6, height=4)
bar <- get_pathway_barplot(gse_res=subset(pancancer_gse_res$Reactome, NES<0), gsetSource="Reactome", 
                           extra_data=homo_Reactome_pathways, top=20)
print(bar)
bar <- get_pathway_barplot(gse_res=subset(pancancer_gse_res$Reactome, NES>0), gsetSource="Reactome", 
                           extra_data=homo_Reactome_pathways, top=20)
print(bar)
dev.off()

pdf(paste0(res_path, "Reactome_gse_barplot.pdf"))
bar <- get_pathway_barplot(gse_res=subset(pancancer_gse_res$Reactome, NES<0), gsetSource="Reactome", 
                           extra_data=homo_Reactome_pathways)
print(bar)
bar <- get_pathway_barplot(gse_res=subset(pancancer_gse_res$Reactome, NES>0), gsetSource="Reactome", 
                           extra_data=homo_Reactome_pathways)
print(bar)
dev.off()



###  
#' 
#' Identify enrichment of up-down-regulation pathways for specific cancer types
################################
cancers <- c("HNSC", "LUAD", "CESC", "KIRC", "SKCM", "STAD", "LUSC", "BLCA", "BRCA" ,
             "ESCA","OV" ,"THCA", "COREAD", "UCEC", "LIHC")
load("result/embedding_clusters_biology_survival_diff/diff_expression_analysis/cancers_diff_exp_genes_df.RData")

source("code/Rpipeline/function/get_pathway_enrich.R")
set.seed(1234)

cancers_gse_res <- lapply(cancers, function(cancer){
  
  logFC <- cancers_diff_exp_genes_df[[cancer]]$logFC
  names(logFC) <- rownames(cancers_diff_exp_genes_df[[cancer]])
  genelist <- logFC[order(logFC, decreasing = T)]
  
  gse_res <- lapply(c("GO", "KEGG", "Reactome"), function(gsetSource){
    gse_res <- get_pathway_enrich(genelist=genelist, gsetSource=gsetSource, keyType='SYMBOL',
                                  method="gse", pvalueCutoff=1) 
    gse_res@result
  })
  names(gse_res) <- c("GO", "KEGG", "Reactome")
  gse_res
})
names(cancers_gse_res) <- cancers
save(cancers_gse_res, file=paste0(res_path, "cancers_gse_res.RData"))


#######
# 
#' Only the enrichment results of KEGG on the differentially expressed gene rank list of each cancer type were visualized
#######
load(paste0(res_path, "cancers_gse_res.RData"))
load(paste0(res_path,"pancancer_gse_res.RData"))
cancers_gse_res <- c(list(pancancer_gse_res), cancers_gse_res)
names(cancers_gse_res)[1] = "PANCANCER"
disease_pathway <- paste0("hsa",KEGG_classification$Pathway.ID[which(KEGG_classification$Pathway.Class.1=="Human Diseases")])

sig_pathways <- unique(unlist(lapply(cancers_gse_res, function(x){
  KEGG_res <-x$KEGG 
  KEGG_res <- KEGG_res[-(which(KEGG_res$ID %in% disease_pathway)),]
  KEGG_res[which(KEGG_res$p.adjust<0.01), "ID"]
})))  # 101

sig_pathway_info <- KEGG_classification[which(paste0("hsa",KEGG_classification$Pathway.ID) %in% sig_pathways),]  #96  4

id.order <- sig_pathway_info$Pathway.Name


plot.data <- do.call(cbind, lapply(cancers_gse_res, function(x){
  KEGG_res <-x$KEGG
  rownames(KEGG_res)<- KEGG_res$Description
  sig_KEGG_pathway_p <- KEGG_res[id.order, "p.adjust"]
  log10p <- (-log10(sig_KEGG_pathway_p))* sign(KEGG_res[id.order, "NES"])
  log10p
}))
colnames(plot.data) <- names(cancers_gse_res)
rownames(plot.data) <- id.order

# visualization
plot.data[which(is.na(plot.data))] <- 0


sig.info <- apply(plot.data, 2, function(x) ifelse(abs(x)>abs(log10(0.01)), "*", ""))


library(pheatmap)
bk <-  c(seq(-5,-0.1,by=0.01),seq(0,5,by=0.01))

annotation_row  = data.frame(pathway_Class=sig_pathway_info$Pathway.Class2)
annotation_row$pathway_Class[which(sig_pathway_info$Pathway.Class.1=="Metabolism")] <- "Metabolism"
row.names(annotation_row) <- sig_pathway_info$Pathway.Name

colors = c(colorRampPalette(colors = c("#0072B5FF","white"))(length(bk)/2),colorRampPalette(colors = c("white","#E18727FF"))(length(bk)/2))
pheatmap(plot.data, scale = 'none', 
         cluster_row = FALSE, cluster_cols = FALSE, 
         treeheight_row = 0, treeheight_col = 0,
         display_numbers =sig.info, fontsize_number=6, 
         breaks=bk,
         annotation_row = annotation_row,
         color = colors, 
         annotation_names_col =FALSE,
         border_color = NA, fontsize_row = 3, fontsize_col = 6, legend = TRUE,
         legend_break=seq(-5, 5, by=5), fontsize=3,
         height = 8, width = 10, 
         filename = paste0(res_path, "cancers_KEGG_pathway_GSEA_heatmap.pdf"))

