source("call_libraries.R")
source("Initialization.R")

integrated_cancer_cells = readRDS("integrated_cancer_cells.rds")
Idents(integrated_cancer_cells) = as.character(integrated_cancer_cells$integrated_snn_res.0.2)
all_markers = FindAllMarkers(integrated_cancer_cells, only.pos = TRUE, 
                             min.pct = 0.25, logfc.threshold = 0.25)
top_n = 50
top_n_markers = all_markers %>% group_by(cluster) %>% top_n(n = top_n, wt = avg_log2FC)

# output top n markers
write.csv(top_n_markers, "top_n_markers.csv", row.names = F)
write.csv(all_markers, "all_markers.csv", row.names = F)