
setwd("../../")
source("code/Rpipeline/function/get_pathway_enrich.R")
library(enrichplot)
set.seed(1234)


res_path <- "result/embedding_clusters_biology_survival_diff/DEGs_enriched_pathways/"
output_dir <- file.path(res_path)
if (!dir.exists(output_dir))  dir.create(output_dir)
load("result/embedding_clusters_biology_survival_diff/diff_expression_analysis/diff_exp_genes_df.RData")
diff.exp.geneSet <- rownames(diff_exp_genes_df)[which(diff_exp_genes_df$Regulate %in% c("Up−regulated", "Down−regulated"))]


# load KEGG classification(Note: 2023-9-15 Data download athttps://www.genome.jp/brite/hsa00001.keg）
KEGG_classification <- read.table("data/pathway/KEGG-pathway-classification.txt",
                                  header =T, sep="\t", as.is=TRUE,
                                  colClasses= "character")
# import Reactome pathways (（2023-09-01）https://reactome.org/download/current/)
homo_Reactome_pathways <- read.table("data/pathway/homo_Reactome_pathways.txt",
                                     header =T, sep="\t",stringsAsFactors=FALSE,quote = "")


###  
#' 
#' enrichment analysis

################################
### GO
gsetSource <- "GO"
GO_enrich_res <- get_pathway_enrich(genelist=diff.exp.geneSet, gsetSource=gsetSource, keyType='SYMBOL',
                                    method="enrich", pvalueCutoff=0.01) 
length(GO_enrich_res$ID) #  1037

pdf(paste0(res_path, "go_enrich_geneset_top20_dotplot.pdf"))
dotplot(GO_enrich_res, showCategory=20, color = 'p.adjust', font.size=8)
dev.off()

### KEGG
gsetSource <- "KEGG"
KEGG_enrich_res <- get_pathway_enrich(genelist=diff.exp.geneSet, gsetSource=gsetSource, keyType='SYMBOL',
                                      method="enrich", pvalueCutoff=0.01) 
length(KEGG_enrich_res$ID) # 59
# The function uses the KEGG database source：https://github.com/SciLicium/kegg-db （The data collection time is：2023-03-21）
# Disease-related pathways are not considered
disease_pathway <- paste0("hsa",KEGG_classification$Pathway.ID[which(KEGG_classification$Pathway.Class.1=="Human Diseases")])
categorys <- KEGG_enrich_res$Description[-which(KEGG_enrich_res$ID %in% disease_pathway)]
# 27
pdf(paste0(res_path, "KEGG_enrich_geneset_all_dotplot.pdf"))
enrichplot::dotplot(KEGG_enrich_res, showCategory=categorys,  color = 'p.adjust', font.size=8)
dev.off()
KEGG_enrich_res <- subset(as.data.frame(KEGG_enrich_res), Description %in% categorys)

### Reactome
gsetSource <- "Reactome"
Reactome_enrich_res <- get_pathway_enrich(genelist=diff.exp.geneSet, gsetSource=gsetSource, keyType='SYMBOL',
                                          method="enrich", pvalueCutoff=0.01) 
length(Reactome_enrich_res$ID) # 63

pdf(paste0(res_path, "Reactome_enrich_geneset_top20_dotplot.pdf"))
dotplot(Reactome_enrich_res, showCategory=20, color = 'p.adjust', font.size=8)
dev.off()


pancancer_enrich_res <- list(GO=GO_enrich_res, KEGG=KEGG_enrich_res, Reactome=Reactome_enrich_res)
save(pancancer_enrich_res, file=paste0(res_path,"pancancer_enrich_res.RData"))
outputs <- rbind(rbind(as.data.frame(pancancer_enrich_res$KEGG), 
                       as.data.frame(pancancer_enrich_res$Reactome)), 
                 as.data.frame(pancancer_enrich_res$GO)[,-1])
outputs$Pathway_source <- c(rep("KEGG", nrow(pancancer_enrich_res$KEGG)),
                            rep("Reactome", nrow(pancancer_enrich_res$Reactome)),
                            rep("GO", nrow(pancancer_enrich_res$GO)))
outputs <- outputs[,c(10, 1:9)]
write.csv(outputs, file=paste0(res_path, "pancancer_enrich_res.csv"), quote =F, row.names=F)




###  
#' 
#' Barplots visualizes enrichment results
################################

pdf(paste0(res_path, "go_enrich_top20_barplot.pdf"), width=6, height=4)
bar <- get_pathway_barplot(gse_res=GO_enrich_res, gsetSource="GO", 
                           top=20)
print(bar)
dev.off()


pdf(paste0(res_path, "KEGG_enrich_all_barplot.pdf"))
bar <- get_pathway_barplot(gse_res=KEGG_enrich_res, gsetSource="KEGG", 
                           extra_data=KEGG_classification)
print(bar)
dev.off()


### Reactome
pdf(paste0(res_path, "Reactome_enrich_top20_barplot.pdf"), width=6, height=4)
bar <- get_pathway_barplot(gse_res=Reactome_enrich_res, gsetSource="Reactome", 
                           extra_data=homo_Reactome_pathways, top=20)
print(bar)
dev.off()

pdf(paste0(res_path, "Reactome_enrich_all_barplot.pdf"))
bar <- get_pathway_barplot(gse_res=Reactome_enrich_res, gsetSource="Reactome", 
                           extra_data=homo_Reactome_pathways)
print(bar)
dev.off()
