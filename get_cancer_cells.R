source("call_libraries.R")
plan("multicore", workers = 12)
tissue = "Tumor"
cd45_status = "CD45-"
integrated = readRDS("cd45neg_tumor_integrated_SCT_nonRPCA_azimuth.rds")
Idents(integrated) = integrated$integrated_snn_res.0.05
#defining cancer cells by cluster 0 + CA9 expression
DefaultAssay(integrated) = "RNA"
integrated_cancer_cell = subset(x = integrated, idents = 0)
#focus on potential cancer and integrated them to remove batch effect
integrated_cancer_cell = DietSeurat(integrated_cancer_cell, assays = "RNA")
integrated_cancer_cell = Integrate_seurat_obj_by_key(integrated_cancer_cell, key = "patient",
                                                    method_SCT = T, method_rpca = T)
# run subclustering of these cancer cells and save to rds
DefaultAssay(integrated_cancer_cell) = "integrated"
integrated_cancer_cell = run_pca_and_umap(integrated_cancer_cell)
#integrated_cancer_cell = find_neighbours_and_clusters(integrated_cancer_cell, res = seq(0.05,0.3, 0.05))
integrated_cancer_cell <- FindNeighbors(integrated_cancer_cell, dims = 1:30, verbose = F)
integrated_cancer_cell = FindClusters(integrated_cancer_cell, 
                                      resolution = seq(0.05,0.3, 0.05),
                                      algorithm = 4)
saveRDS(integrated_cancer_cell, "integrated_cancer_cells_SCT_leiden.rds")