
get_pathway_enrich <- function(genelist, keyType='SYMBOL', gsetSource, 
                               method="enrich",pvalueCutoff=0.05,qvalueCutoff=0.05){
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(fgsea)
  # GO
  if(gsetSource=="GO" & keyType %in% c('SYMBOL', 'ENTREZID')){ 
    if(method=="enrich" & is.null(names(genelist))){
      gse_res <- enrichGO(genelist, ont = "ALL", OrgDb=org.Hs.eg.db, keyType = keyType,  
                          pvalueCutoff = pvalueCutoff, pAdjustMethod = "BH", qvalueCutoff=qvalueCutoff)
    }else if(method=="gse" & (!is.null(names(genelist)))){ 
      genelist <- sort(genelist,decreasing=T)
      gse_res <- gseGO(genelist, ont = "ALL", OrgDb=org.Hs.eg.db, keyType = keyType, eps=0,
                       pvalueCutoff = pvalueCutoff, pAdjustMethod = "BH")
    }else{
      print("Error input for genelist. geneset for enrich, ranked genelist for gse")
      gse_res=NULL
    }
  }
  
  # KEGG
  else if(gsetSource %in% c("KEGG", "Reactome") & keyType %in% c('SYMBOL', 'ENTREZID')){

    if(keyType=='SYMBOL' & is.null(names(genelist))){
      gene.df <- bitr(genelist,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                      OrgDb = org.Hs.eg.db)
      genelist <- gene.df$ENTREZID
    }else if(keyType=='SYMBOL' & (!is.null(names(genelist)))){
      gene.tx <- bitr(names(genelist),fromType="SYMBOL",toType=c("ENTREZID"),
                      OrgDb = org.Hs.eg.db)
      if(length(which(duplicated(gene.tx[,1])))>0) gene.tx <- gene.tx[-which(duplicated(gene.tx[,1])),]
      colnames(gene.tx)[1] <- "gene_name"
      mydata = data.frame(gene_name=names(genelist), values = genelist)
      gene.tx <- merge(gene.tx,mydata,by="gene_name")
      genelist <- gene.tx$values #numeric vector
      names(genelist) <- as.character(gene.tx$ENTREZID) #named vector
      genelist <- sort(genelist,decreasing=T) #decreasing order
    }
    
    
    if(gsetSource == 'KEGG'){
      ## KEGG
      library(KEGG.db)
      if(method=="enrich" & is.null(names(genelist))){
       
        gse_res <- enrichKEGG(genelist, organism='hsa',keyType="kegg",
     
                              pvalueCutoff=pvalueCutoff,pAdjustMethod='BH',
                              qvalueCutoff=qvalueCutoff,use_internal_data=T)
      }else if(method=="gse" & (!is.null(names(genelist)))){

        gse_res <- gseKEGG(genelist, organism='hsa',keyType="kegg",

                           pvalueCutoff=pvalueCutoff,pAdjustMethod='BH', use_internal_data=T)
      }else{
        print("Error input for genelist. geneset for enrich, ranked genelist for gse")
        gse_res=NULL
      }
    }else{
      library(ReactomePA)
      # Reactome
      if(method=="enrich" & is.null(names(genelist))){

        gse_res = enrichPathway(gene=genelist,organism= "human",
               
                                pvalueCutoff=pvalueCutoff, pAdjustMethod= "BH",
                                qvalueCutoff=qvalueCutoff, readable=T)
        
      }else if(method=="gse" & (!is.null(names(genelist)))){

        gse_res = gsePathway(geneList=genelist,organism= "human",eps=0,

                             pvalueCutoff=pvalueCutoff, pAdjustMethod= "BH")
      }else{
        print("Error input for genelist. geneset for enrich, ranked genelist for gse")
        gse_res=NULL
      }
    }
  }else{
    print("gsetSource must be GO, KEGG or Reactome, and keyType must be SYMBOL, ENTREZID")
    gse_res=NULL
  }
  
  return(gse_res)
  
}





get_pathway_barplot <- function(gse_res, gsetSource, extra_data=NULL, var="",top=NULL){
  library(ggplot2)
  library(stringr)
  library(ggsci)
  

  gse_res.df <- as.data.frame(gse_res)
  gse_res.df$p.adjust <- -log10(gse_res.df$p.adjust)
  if(!is.null(top)){
    gse_res.df <- gse_res.df[order(gse_res.df$p.adjust, decreasing = T),]
    if(dim(gse_res.df)[1]<top) top=dim(gse_res.df)[1]
    gse_res.df <- gse_res.df[1:top,]
  }

  if(gsetSource=="GO"){

    names(gse_res.df)[1] <- "Father_type"

    plot.data <- matrix(rep(0, dim(gse_res.df)[2]), nrow=1)
    colnames(plot.data) <- colnames(gse_res.df)
    for(types in c("BP", "MF", "CC")){
      per.df <- subset(gse_res.df, Father_type==types)
      per.df <- per.df[order(per.df$p.adjust, decreasing = F),] # 按p值大小排序
      plot.data <- rbind(plot.data, per.df)
    }
    plot.data <- as.data.frame(plot.data[-1,])
    plot.data$Description <- factor(plot.data$Description, levels=plot.data$Description)
    plot.data$Father_type <- factor(plot.data$Father_type, levels=c("BP", "MF", "CC"))
    colors = pal_d3("category20", alpha=0.7)(3)
    
  }else if(gsetSource=="KEGG"){

    Father_type <- unlist(lapply(gse_res.df$Description, function(x){
      father <- extra_data[which(extra_data$Pathway.Name==x),4]
      if(length(father)==0) father <- ""
      father[1]
    }))
    gse_res.df$Father_type <- Father_type
    

    plot.data <- matrix(rep(0, dim(gse_res.df)[2]), nrow=1)
    colnames(plot.data) <- colnames(gse_res.df)
    for(types in unique(extra_data[,4])){
      per.df <- subset(gse_res.df, Father_type==types)
      per.df <- per.df[order(per.df$p.adjust, decreasing = F),] 
      plot.data <- rbind(plot.data, per.df)
    }
    plot.data <- as.data.frame(plot.data[-1,])
    plot.data$Description <- factor(plot.data$Description, levels=plot.data$Description)
    plot.data$Father_type <- factor(plot.data$Father_type, levels=unique(plot.data$Father_type))
    colors = pal_d3("category20", alpha=0.7)(length(unique(gse_res.df$Father_type)))
    
  }else if(gsetSource=="Reactome"){

    Father_type <- unlist(lapply(gse_res.df$Description, function(x){
      father <- extra_data[which(extra_data[,2]==x),4]
      if(length(father)==0) father <- ""
      unique(father)
    }))
    gse_res.df$Father_type <- Father_type

    plot.data <- matrix(rep(0, dim(gse_res.df)[2]), nrow=1)
    colnames(plot.data) <- colnames(gse_res.df)
    for(types in unique(gse_res.df$Father_type)){
      per.df <- subset(gse_res.df, Father_type==types)
      per.df <- per.df[order(per.df$p.adjust, decreasing = F),] 
      plot.data <- rbind(plot.data, per.df)
    }
    plot.data <- as.data.frame(plot.data[-1,])
    plot.data$Description <- factor(plot.data$Description, levels=plot.data$Description)
    plot.data$Father_type <- factor(plot.data$Father_type, levels=unique(plot.data$Father_type))
    colors = pal_d3("category20", alpha=0.7)(length(unique(gse_res.df$Father_type)))
  }else{
    stop("ERROR:gsetSource must be GO, KEGG or Reactome.")
  }

  bar_plot <- ggplot(data = plot.data, 
                   aes(x = Description, y = p.adjust,fill = Father_type))+
    geom_bar(stat = "identity")+
    coord_flip()+theme_bw()+ 
    scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ 
    xlab(gsetSource)+ylab("-log10(adj.P)")+ ggtitle(var)+ 
    scale_fill_manual(values=colors)+
    theme(axis.title = element_text(size = 13), 
          axis.text = element_text(size = 11), 
          plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), 
          legend.title = element_text(size = 8), 
          legend.text = element_text(size = 5),
          axis.text.y = element_text(size = 5))

  
  return(bar_plot)
  
}

