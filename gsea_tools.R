#gsea function
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
run_gsea_kegg = function(correlated_gene_vector){
  #input: correlated_gene_vector
  #input_type: named numerical vector with name as genes
  ids<-bitr(names(correlated_gene_vector), fromType = "SYMBOL", toType = "ENTREZID", 
            OrgDb=org.Hs.eg.db)
  #remove duplicated gene symbol
  dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
  correlated_gene_vector = correlated_gene_vector[names(correlated_gene_vector)
                                                  %in% dedup_ids$SYMBOL]
  #rename vector with ENTREZID
  names(correlated_gene_vector) = dedup_ids$ENTREZID
  #sort correlated gene by descending order
  correlated_gene_vector = sort(correlated_gene_vector, decreasing = T)
  
  #run GSEA KEGG terms
  GSEA_kegg <- gseKEGG(geneList = correlated_gene_vector,
                       organism     = "hsa",
                       minGSSize    = 10,
                       maxGSSize    = 1000,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "none",
                       keyType       = "ncbi-geneid")
  
  return(GSEA_kegg)
}

#dotplot function
plot_dot = function(gsea_obj, plot_title){
  return(dotplot(gsea_obj, showCategory = 10, 
          title = plot_title , split=".sign")+ ggplot2::facet_grid(.~.sign)+
          theme(plot.title = element_text(color = "blue", size = 20, face = "bold", hjust = 0.5))
  )
}