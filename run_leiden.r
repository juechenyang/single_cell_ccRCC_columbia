source("call_libraries.R")
plan("multicore", workers = 12)
integrated_cancer_cell = readRDS("integrated_cancer_cells.rds")
integrated_cancer_cell = FindClusters(integrated_cancer_cell, 
                                      resolution = seq(0.05,0.3, 0.05),
                                      algorithm = 4,
				      method = "igraph")
saveRDS(integrated_cancer_cell, "integrated_cancer_cells_SCT_leiden.rds")
