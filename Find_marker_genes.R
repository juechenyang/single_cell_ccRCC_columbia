source("call_libraries.R")
source("Initialization.R")
source("gsea_tools.R")
#plan("multicore", workers = 6)
integrated_cancer_cells = readRDS("integrated_cancer_cells.rds")
DefaultAssay(integrated_cancer_cells) = "RNA"
integrated_cancer_cells = NormalizeData(integrated_cancer_cells)
integrated_cancer_cells$stage = 
  ifelse(integrated_cancer_cells$patient=="Patient5", "Metastatic", integrated_cancer_cells$stage)
#Idents(integrated_cancer_cells) = as.character(integrated_cancer_cells$integrated_snn_res.0.2)
plot_list = list()
for(x in 0:9){
  inCluster = integrated_cancer_cells$integrated_snn_res.0.2==x
  mat = data.frame(integrated_cancer_cells[["RNA"]]@data, check.names = F)
  t_value = apply(mat, MARGIN = 1, FUN = function(x){
    cluster = x[inCluster]
    other = x[!inCluster]
    t_result = t.test(cluster, other)
    return(t_result$statistic)
  })
  gsea_obj = run_gsea_kegg(t_value)
  #png("gsea_cluster6.png", width = 10, height = 10, units = "in", res = 300)
  plot_title = paste0("cluster ", x, "enrichment_dotplot")
  each_cluster_plot = dotplot(gsea_obj, showCategory = 10, 
                              title = plot_title , split=".sign")+ facet_grid(.~.sign)
  plot_list = append(plot_list, list(each_cluster_plot))
}

png("gsea_plots.png", units = "in", res = 400, 
    width = 32, height = 24)
ggarrange(
  plotlist = plot_list,
  ncol = 4,
  nrow = 3
)
dev.off()


####################################################################################
###Compare cluster 4 between stage 1 and stage 3
####################################################################################
for(x in c("pT1b","pT3a","Metastatic")){
  integrated_cancer_cell_stage = subset(integrated_cancer_cell, subset = stage == x)
  
}

cluster_number = 2
stage1_cluster <- (integrated_cancer_cells$integrated_snn_res.0.2==cluster_number & integrated_cancer_cells$stage == "pT1b")
stage3_cluster <- (integrated_cancer_cells$integrated_snn_res.0.2==cluster_number & integrated_cancer_cells$stage == "pT3a")
mat = data.frame(integrated_cancer_cells[["RNA"]]@data, check.names = F)
t_value = apply(mat, MARGIN = 1, FUN = function(x){
  stage1 = x[stage1_cluster]
  stage3 = x[stage3_cluster]
  t_result = t.test(stage3, stage1)
  return(t_result$statistic)
})
gsea_obj = run_gsea_kegg(t_value)
png("gsea_cluster2_stage3_vs_stage1.png", width = 10, height = 10, units = "in", res = 300)
dotplot(gsea_obj, showCategory = 10, split=".sign")+ 
        facet_grid(.~.sign)+
        ggtitle("cluster2 stage3 vs stage1")
dev.off()




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