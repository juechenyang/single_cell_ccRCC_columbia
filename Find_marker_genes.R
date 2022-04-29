source("call_libraries.R")
source("Initialization.R")
#plan("multicore", workers = 6)
integrated_cancer_cells = readRDS("integrated_cancer_cells.rds")
Idents(integrated_cancer_cells) = as.character(integrated_cancer_cells$integrated_snn_res.0.2)
inCluster6 = integrated_cancer_cells$integrated_snn_res.0.2==6
DefaultAssay(integrated_cancer_cells) = "RNA"
integrated_cancer_cells = ScaleData(integrated_cancer_cells)
mat = data.frame(integrated_cancer_cells[["RNA"]]@scale.data, check.names = F)
t_value = apply(mat, MARGIN = 1, FUN = function(x){
  cluster6 = x[inCluster6]
  other = x[!inCluster6]
  t_result = t.test(cluster6, other)
  return(t_result$statistic)
})

#gsea function
library(clusterProfiler)
library(org.Hs.eg.db)
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

gsea_obj = run_gsea_kegg(t_value)
dotplot(gsea_obj, showCategory = 10, 
        title = "Enriched Pathways" , split=".sign")+ facet_grid(.~.sign)







all_markers = FindAllMarkers(integrated_cancer_cells, only.pos = TRUE, 
                             min.pct = 0.25, logfc.threshold = 0.25)
cluster6_roc_markers = FindMarkers(integrated_cancer_cells, 
                                   ident.1 = 6, test.use = "t", logfc.threshold = 0,
                                   min.pct = 0.25)
top_n = 50
top_n_markers = all_markers %>% group_by(cluster) %>% top_n(n = top_n, wt = avg_log2FC)

# output top n markers
write.csv(top_n_markers, "top_n_markers.csv", row.names = F)
write.csv(all_markers, "all_markers.csv", row.names = F)