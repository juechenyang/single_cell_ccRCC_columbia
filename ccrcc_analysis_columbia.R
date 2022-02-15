# prepare data
source("ccrcc_analysis_columbia_prepare_data.R")
source("single_cell_azimuth_g23_enrichment.R")
# remove redundant data
rm(list=setdiff(ls(), c("seurat_list_cd45pos","seurat_list_cd45neg", "plot_dir", "rds_dir", "test.combined.sct")))
rm(list=setdiff(ls(), c("seurat_list_cd45pos","seurat_list_cd45neg", "plot_dir", 
                        "rds_dir", "seurat_list_cd45neg_tumor", "seurat_list_cd45neg_normal")))
rm(list=setdiff(ls(), c("big_list")))

# Integrate all CD45- samples
set.seed(1234)
options(future.globals.maxSize = 15000 * 1024^2)
#test_list = seurat_list_cd45pos
test_list = SplitObject(seurat_list_cd45neg[[2]], split.by = "tissue")
regular = T
if(regular){
  SCT_test_list <- lapply(X = test_list, FUN = SCTransform)
}else{
  SCT_test_list <- lapply(X = test_list, FUN = SCTransform, method = "glmGamPoi")
}
rm(test_list)
features <- SelectIntegrationFeatures(object.list = SCT_test_list, nfeatures = 8000)
PrepSCT_test_list <- PrepSCTIntegration(object.list = SCT_test_list, anchor.features = features)
rm(SCT_test_list)
PCA_test_list <- lapply(X = PrepSCT_test_list, FUN = RunPCA, features = features)
rm(PrepSCT_test_list)
if(regular){
  test.anchors <- FindIntegrationAnchors(object.list = PCA_test_list, reference = 3,
                                         normalization.method = "SCT",anchor.features = features)
}else{
  test.anchors <- FindIntegrationAnchors(object.list = PCA_test_list, normalization.method = "SCT", reference = 3,
                                         anchor.features = features, reduction = "rpca")
}
#rm(PCA_test_list)
test.combined.sct <- IntegrateData(anchorset = test.anchors, normalization.method = "SCT", verbose = F)
#test.combined.sct <- test_list[[1]]
test.combined.sct = SCTransform(test.combined.sct, return.only.var.genes = F, verbose = F)


# test.anchors <- FindIntegrationAnchors(object.list = PCA_test_list, normalization.method = "SCT", reference = 1,
#                                          anchor.features = features, reduction = "rpca", dims = 1:30, k.anchor = 20)
# #rm(PCA_test_list)
# test.combined.sct <- IntegrateData(anchorset = test.anchors, normalization.method = "SCT", verbose = F, dims = 1:30)

# Run PCA and UMAP of integrated dataset
test.combined.sct <- RunPCA(test.combined.sct, verbose = FALSE, features = VariableFeatures(test.combined.sct))
test.combined.sct <- RunUMAP(test.combined.sct, dims = 1:30, verbose = FALSE,umap.method="umap-learn",metric="correlation")
test.combined.sct <- FindNeighbors(test.combined.sct, dims = 1:30, verbose = F)
test.combined.sct <- FindClusters(test.combined.sct, resolution = seq(0.01, 1, 0.01), verbose = F, algorithm = 1)
ref <- BlueprintEncodeData()
#ref <- HumanPrimaryCellAtlasData()
test.combined.sct.diet = DietSeurat(test.combined.sct, graphs = "umap")
test.combined.sct.diet.sce = as.SingleCellExperiment(test.combined.sct.diet)
pred.hesc <- SingleR(test = test.combined.sct.diet.sce, ref = ref, assay.type.test=1, labels = ref$label.main)
test.combined.sct$sr_ct = pred.hesc$labels

Idents(test.combined.sct) = test.combined.sct$SCT_snn_res.0.09
markers <- FindAllMarkers(test.combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
top <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
test.combined.sct = label_cluster_locally(test.combined.sct, top, type = "or")

saveRDS(test.combined.sct, "all_neg_seurat_8000features.rds")



g23 = c("ATP5MC1","ATP5MC2","COA3","NDUFA5","NDUFAB1","NDUFAF3","NDUFS3","UQCRC1","MRPL12","MRPL14","
        MRPL38","MRPL47","MRPS10","MRPS18A","MRPS24","MRPS36","MRPS5","HLA-DMB","HLA-DOA","HLA-DPB1","HLA-DQA2","HLA-DQB2","HLA-DRA")
FeaturePlot(test.combined.sct, features = g23)


pred.hesc <- SingleR(test = test.combined.sct.diet.sce, ref = ref, assay.type.test=1, 
                     clusters = test.combined.sct$integrated_snn_res.0.75,
                     labels = ref$label.main)
cell_types = pred.hesc$labels
names(cell_types) = levels(unique(test.combined.sct$integrated_snn_res.0.75))
Idents(test.combined.sct) = test.combined.sct$integrated_snn_res.0.75
test.combined.sct = RenameIdents(test.combined.sct, cell_types)
test.combined.sct$sr_ct_0.75 = Idents(test.combined.sct)
DimPlot(test.combined.sct, reduction = "umap", group.by = "tissue", pt.size = 0.2, label = F)
DimPlot(test.combined.sct, reduction = "umap", group.by = "patient", pt.size = 0.2, label = F)
DimPlot(test.combined.sct, reduction = "umap", pt.size = 0.2, label = F)
DimPlot(test.combined.sct, reduction = "umap", label = F, group.by = "integrated_snn_res.0.75")
DimPlot(test.combined.sct, reduction = "umap", pt.size = 0.2, label = F, group.by = "sr_ct")

Idents(test.combined.sct) = test.combined.sct$integrated_snn_res.0.1
top = mtcars
test.combined.sct = label_cluster_locally(test.combined.sct, top, type = "gsea")




#predict cell type for each cell
pred.hesc <- SingleR(test = test.combined.sct.diet.sce, ref = ref, assay.type.test=1, labels = ref$label.main)
test.combined.sct$sr_ct = pred.hesc$labels
p_singler_highlight = DimPlot(test.combined.sct, reduction = "umap",group.by = "sr_ct", pt.size = 0.2, 
                              cells.highlight = names(test.combined.sct$sr_ct)[test.combined.sct$sr_ct=="Epithelial cells"])
p_singler = DimPlot(test.combined.sct, reduction = "umap",group.by = "sr_ct", pt.size = 0.2)
p_singler_highlight + p_singler

#split by patient
test_list_after = SplitObject(test.combined.sct, split.by = "patient")




#*********************************************************************************************
#*patient level analysis

#used for clear env
rm(list=setdiff(ls(), c("seurat_list_cd45neg", "seurat_list_cd45pos", "plot_dir", "rds_dir")))
source("single_cell_azimuth_g23_enrichment.R")
set.seed(1234)
patient_pos_seurat_objs = list()
markers_pos_list = list()
for(iter in 1:length(seurat_list_cd45pos)){
  print(iter)
  patient_cd45pos=seurat_list_cd45pos[[iter]]
  patientnumber = unique(patient_cd45pos$patient)
  patient_cd45pos <- PercentageFeatureSet(patient_cd45pos, pattern = "^MT-", col.name = "percent.mt")
  patient_cd45pos <- subset(patient_cd45pos, subset = percent.mt < 10 & nCount_RNA > 1000 & nCount_RNA < 15000)
  patient_cd45pos <- SCTransform(patient_cd45pos, return.only.var.genes = F,verbose = F)
  patient_cd45pos <- RunPCA(patient_cd45pos, features = VariableFeatures(object = patient_cd45pos))
  patient_cd45pos <- RunUMAP(patient_cd45pos, dims = 1:30, verbose = FALSE,umap.method="umap-learn",metric="correlation")
  patient_cd45pos <- FindNeighbors(patient_cd45pos, dims = 1:30, verbose = FALSE)
  res = 0.15
  # cluster_res = sapply(grep("res",colnames(patient_cd45pos@meta.data),value = TRUE),function(x)         length(unique(patient_cd45pos@meta.data[,x])))
  #for(res in seq(0.25, 1, 0.25)){
  patient_cd45pos <- FindClusters(patient_cd45pos, resolution=res, verbose = FALSE,algorithm=1)
  markers <- FindAllMarkers(patient_cd45pos, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
  patient_seurat_objs[[iter]] = patient_cd45pos
  markers_list[[iter]] = markers
}
label_cell_gsea = T
res = 0.15
for(iter in 1:length(seurat_list_cd45pos)){
  #iter = 1
  patient_cd45pos = patient_seurat_objs[[iter]]
  markers = markers_list[[iter]]
  patientnumber = unique(patient_cd45pos$patient)
  top <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  enrich_source = "local"
  if(label_cell_gsea){
    patient_cd45pos = label_cluster_locally(patient_cd45pos, top, type = "or")
  }else{
    #patient_cd45pos <- RenameIdents(patient_cd45pos, reverse_cell_type)
    patient_cd45pos = label_cluster_enrichR(patient_cd45pos, top)
    enrich_source = "enrichR"
  }
  p = DimPlot(patient_cd45pos, reduction = "umap", label = T, pt.size = 0.2)+NoLegend()
  ggsave(paste0(plot_dir, patientnumber,"_umap_res_", res,"_source_", enrich_source, "_cd45pos.png"), p, width = 12, height = 8, units = "in")
  
}

#used for clear env
rm(list=setdiff(ls(), c("seurat_list_cd45neg", "seurat_list_cd45pos", "plot_dir", "rds_dir")))
source("single_cell_azimuth_g23_enrichment.R")
set.seed(1234)
for(iter in 1:length(seurat_list_cd45neg)){
  print(iter)
  patient_cd45neg=seurat_list_cd45neg[[iter]]
  patientnumber = unique(patient_cd45neg$patient)
  patient_cd45neg <- PercentageFeatureSet(patient_cd45neg, pattern = "^MT-", col.name = "percent.mt")
  patient_cd45neg <- subset(patient_cd45neg, subset = percent.mt <= 10 & 
                              nCount_RNA > 1000 & nCount_RNA <= 25000 & 
                              nFeature_RNA<=5000 & nFeature_RNA>=500)
  patient_cd45neg <- SCTransform(patient_cd45neg, return.only.var.genes = F,verbose = F)
  patient_cd45neg <- RunPCA(patient_cd45neg, features = VariableFeatures(object = patient_cd45neg))
  patient_cd45neg <- RunUMAP(patient_cd45neg, dims = 1:30, verbose = FALSE,umap.method="umap-learn",metric="correlation")
  patient_cd45neg <- FindNeighbors(patient_cd45neg, dims = 1:30, verbose = FALSE)
  patient_cd45neg <- FindClusters(patient_cd45neg, resolution=seq(0.05,1,0.05), verbose = FALSE,algorithm=1)
  seurat_list_cd45neg[[iter]] = patient_cd45neg
}
label_cell_gsea = T
for(iter in 1:length(seurat_list_cd45neg)){
  iter = 3
  patient_cd45neg = seurat_list_cd45neg[[iter]]
  #markers = markers_neg_list[[iter]]
  patientnumber = unique(patient_cd45neg$patient)
  
  #feature plot epithelial cells
  # epithelial_markers = c("SERPINA1","CD24","NNMT","CRYAB","ANGPTL4")
  # FeaturePlot(patient_cd45neg, features = epithelial_markers)
  # markers <- FindAllMarkers(patient_cd45neg, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
  # top <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  enrich_source = "local"
  if(label_cell_gsea){
    patient_cd45neg = label_cluster_locally(patient_cd45neg, top, type = "gsea")
  }else{
    patient_cd45neg = label_cluster_enrichR(patient_cd45neg, top)
    enrich_source = "enrichR"
  }
  p = DimPlot(patient_cd45neg, reduction = "umap", label = T, pt.size = 0.2)+NoLegend()
  #Idents(patient_cd45neg) = patient_cd45neg$seurat_clusters
  ggsave(paste0(plot_dir, patientnumber,"_umap_res_", res,"_source_", enrich_source, "_cd45neg.png"), 
         p, width = 12, height = 8, units = "in")
}

# try to label cell clusters using singleR
#ref <- HumanPrimaryCellAtlasData()
ref <- BlueprintEncodeData()
patient_cd45neg.diet = DietSeurat(patient_cd45neg, graphs = "umap")
patient_cd45neg.diet.sce = as.SingleCellExperiment(patient_cd45neg.diet)
pred.hesc <- SingleR(test = patient_cd45neg.diet.sce, ref = ref, assay.type.test=1, clusters = patient_cd45neg$SCT_snn_res.0.05,
                     labels = ref$label.main)
cell_types = pred.hesc$labels
names(cell_types) = levels(unique(patient_cd45neg$SCT_snn_res.0.05))
Idents(patient_cd45neg) = patient_cd45neg$SCT_snn_res.0.05
patient_cd45neg = RenameIdents(patient_cd45neg, cell_types)
patient_cd45neg$sr_ct_0.05 = Idents(patient_cd45neg)
DimPlot(patient_cd45neg, reduction = "umap", group.by = "sr_ct_0.05", pt.size = 0.2, label = T)

pred.hesc <- SingleR(test = patient_cd45neg.diet.sce, ref = ref, assay.type.test=1, labels = ref$label.main)
patient_cd45neg$sr_ct = pred.hesc$labels
p_singler_highlight = DimPlot(patient_cd45neg, reduction = "umap",group.by = "sr_ct", pt.size = 0.2, 
                    cells.highlight = names(patient_cd45neg$sr_ct)[patient_cd45neg$sr_ct=="Epithelial cells"])
p_singler = DimPlot(patient_cd45neg, reduction = "umap",group.by = "sr_ct", pt.size = 0.2)
p_singler_highlight + p_singler

p_seurat_cluster = DimPlot(patient_cd45neg, reduction = "umap",group.by = "SCT_snn_res.0.09", pt.size = 0.2)
p_singler + p_seurat_cluster

t = plotScoreHeatmap(pred.hesc, clusters = patient_cd45neg$SCT_snn_res.0.7, order.by = "clusters", show_colnames = F)
p1 = DimPlot(patient_cd45neg, reduction = "umap", label = T, pt.size = 0.2)+NoLegend()
p2 = DimPlot(patient_cd45neg, reduction = "umap", group.by = "sr_ct_0.05", pt.size = 0.2)
p1 + p2

fold_change = FoldChange(patient_cd45pos, ident.1 = 0)
gene_list = fold_change$avg_log2FC
names(gene_list) = rownames(fold_change)
gene_list = sort(gene_list, decreasing = T)
source("single_cell_azimuth_g23_enrichment.R")
y = get_azimuth_g23_enrichment(gene_list = gene_list)
dotplot(y, split=".sign") + facet_grid(.~.sign)
