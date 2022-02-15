source("call_libraries.R")
rm(list=setdiff(ls(), "patient_raw_list"))

#define markers
hla_markers = c("HLA-DQB2","HLA-DOA","HLA-DRA","HLA-DMB","HLA-DPB1","HLA-DQA2")
mrp_markers = c("MRPL38","MRPS10","MRPL12","MRPL14","MRPL47","MRPS18A","MRPS24","MRPS36","MRPS5")
nduf_markers = c("NDUFAB1", "NDUFAF3", "NDUFA5", "NDUFS3","ATP5MC1","ATP5MC2","COA3","UQCRC1")
cancer_markers = c("CA9","NDUFA4L2","SLC17A3")
EMT_markers = c("MT2A", "SERPINE1","TM4SF1", "BIRC3", "C3", "CAV1", "TGM2", "MMP7",
                "ANXA2", "LOX", "FSTL3", "LGALS1")
PT_markers = c("NAT8", "SLC3A1", "GPX3", "CYB5A", "SLC25A5", "GSTA1")

Fibro = c("ACTA2", "TAGLN")
non_PT_epi = c("DEFB1", "UMOD")
PT = c("GPX3", "GATM")
EC = c("PLAVP", "TIMP3")

#read the raw data
all_raw_data=readRDS("./raw_stratified_rds/all_raw.rds")
#get tissue data
tissue = "Tumor"
cohort = SplitObject(all_raw_data, split.by = "tissue")[[tissue]]
#get cd45 status data
cd45_status = "CD45-"
cohort = SplitObject(cohort, split.by = "cd45")[[cd45_status]]
cohort = PercentageFeatureSet(cohort, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(cohort, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
cohort = subset(cohort, subset = percent.mt <= 30 & 
                  nCount_RNA >= 1000 & nCount_RNA <= 30000 &
                  nFeature_RNA >= 500 & nFeature_RNA <= 5000)
integrated = Integrate_seurat_obj_by_key(cohort, key = "patient",
                                               method_SCT = T, method_rpca = F)
integrated = readRDS("cd45neg_tumor_integrated_SCT_nonRPCA.rds")
integrated_azimuth = run_azimuth(integrated)
integrated$individual_anno = as.character(integrated_azimuth$predicted.annotation.l1)
temp = integrated

run_InferCnv = T
if(run_InferCnv){
  reference_pt = readRDS("reference_pt_cells.rds")
  integrated$epi_marker = 
    ifelse(grepl("Endothelial|Vascular Smooth", 
                 integrated$individual_anno), "yes", "no")
  integrated = subset(integrated, subset = epi_marker == "no")
  patient = merge(integrated, reference_pt)
  want = table(patient$individual_anno)
  want = want[want>=50]
  patient$individual_anno_new = sapply(patient$individual_anno, function(x){
    if(x %in% names(want)){
      return(x)
    }else{
      return("Miscellaneous")
    }
  })
  
  InferCnv_output_dir = paste0("./infercnv_output/",cd45_status,"_",tissue,"_","integrated", "/")
  InferCnv_output_exist_marker = paste0(InferCnv_output_dir, "infercnv.observations.txt")
  if(!file.exists(InferCnv_output_exist_marker)){
    InfercnvInputs = prepare_infercnv_data(patient, annotation = "individual_anno_new")
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=InfercnvInputs[[1]],
                                        annotations_file=InfercnvInputs[[2]],
                                        delim="\t",
                                        gene_order_file=InfercnvInputs[[3]],
                                        ref_group_names = c("Proximal Tubule")
    )
    infercnv_obj = infercnv::run(infercnv_obj,
                                 cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                 out_dir=InferCnv_output_dir, 
                                 cluster_by_groups=TRUE, denoise=TRUE, 
                                 HMM = T, num_threads = 8)
  }
}


integrated$individual_anno_score = integrated_azimuth$predicted.annotation.l1.score
integrated = run_pca_and_umap(integrated)
integrated = find_neighbours_and_clusters(integrated, res = seq(0.05,0.1,0.05))


all_types = as.character(unique(integrated$individual_anno))
palette_color = c("chocolate1", "cyan", "gold", "aquamarine", "deepskyblue", 
                  "cyan4", "darkblue", "darkolivegreen1", "darkorchid1",
                  "firebrick1", "firebrick4", "purple", "darksalmon","darkslategray", "darkmagenta")
all_colors = colorRampPalette(palette_color)(length(all_types))
names(all_colors) = all_types
DimPlot(integrated, reduction = "umap", group.by = "integrated_snn_res.0.1")
p1 = DimPlot(integrated, reduction = "umap", group.by = "individual_anno")+
  scale_color_manual(breaks = all_types, values=all_colors[all_types])
p2 = DimPlot(integrated, reduction = "umap", group.by = "patient")
# DefaultAssay(integrated) = "RNA"
# integrated = NormalizeData(integrated)
p3 = DimPlot(integrated, reduction = "umap", group.by = "integrated_snn_res.0.05")+ggtitle("res0.05_umap")
p4 = FeaturePlot(integrated, features = "CA9")

png("integrated_data_umap.png",20,15, units = "in", res = 300)
ggarrange(
  p1,p2,p3,p4,
  nrow = 2,
  ncol = 2
)
dev.off()


# selected_marker = c("HLA-DRA", "HLA-DPB1", "MRPL14", "MRPS36", "NDUFAB1", "NDUFA5", "ATP5MC1", "ATP5MC2", "COA3")



Idents(integrated) = integrated$integrated_snn_res.0.05
integrated_cancer_cell = subset(x = integrated, subset = CA9 > 0, idents = 0)
#focus on potential cancer cells and do sub clustering
temp_cancer_cell = DietSeurat(temp_cancer_cell, assays = "RNA")
integrated_cancer_cell = Integrate_seurat_obj_by_key(temp_cancer_cell, key = "patient")
#saveRDS(integrated_cancer_cell, "integrated_cancer_cells.rds")
DefaultAssay(integrated_cancer_cell) = "integrated"
integrated_cancer_cell = run_pca_and_umap(integrated_cancer_cell)
integrated_cancer_cell = find_neighbours_and_clusters(integrated_cancer_cell, res = seq(0.05,0.3, 0.05))


integrated_cancer_cell = readRDS("integrated_cancer_cells.rds")
integrated_cancer_cell = AddModuleScore(integrated_cancer_cell, 
                         features = list(hla_markers, mrp_markers, nduf_markers, 
                                         c(mrp_markers, nduf_markers), EMT_markers,
                                         PT_markers), 
                         assay = "RNA", 
                         name = c("hla_score", "mrp_score", "nduf_score", 
                                  "mrp_nduf_score", "emt_score", "pt_score"))
integrated_cancer_cell@meta.data[,c("hla_score", "mrp_score", "nduf_score",
                                    "mrp_nduf_score", "emt_score", "pt_score")] = 
integrated_cancer_cell@meta.data[,c("hla_score1", "mrp_score2", "nduf_score3", 
                                      "mrp_nduf_score4", "emt_score5", "pt_score6")]


p1 = DimPlot(integrated_cancer_cell, reduction = "umap", group.by = "integrated_snn_res.0.3")+ggtitle("seurat_cluster")
p2 = DimPlot(integrated_cancer_cell, reduction = "umap", group.by = "patient")
p3 = FeaturePlot(integrated_cancer_cell, features = c("hla_score", "mrp_score", 
                                                      "nduf_score", "mrp_nduf_score", 
                                                      "emt_score","pt_score" 
                                                      ), ncol = 3, min.cutoff = -0.5, max.cutoff = 0.5,
                                                      cols = c("grey", "red"))

png("integrated_tumor_cells_analysis.png",24,12, units = "in", res = 300)
ggarrange(
  ggarrange(p1, p2, nrow = 2), 
  p3,
  ncol = 2,
  widths = c(1,3)
)
dev.off()




# integrated_cancer_cell = AddModuleScore(integrated_cancer_cell,assay = "RNA",
#                                         features = list(c(hla_markers, mrp_markers, ndufa_markers)), 
#                                         name = "g23_score")
# FeaturePlot(integrated_cancer_cell, features = "g23_score")

DefaultAssay(integrated_cancer_cell) = "RNA"
mat = t(integrated_cancer_cell@meta.data[,c("hla_score", "mrp_score", "nduf_score", 
                                            "mrp_nduf_score", "emt_score", "pt_score")])
mat = data.frame(mat, check.names = F)
order_df = data.frame(cbind(seurat_cluster = as.character(integrated_cancer_cell$integrated_snn_res.0.3),
                            pt = integrated_cancer_cell$patient, stage = integrated_cancer_cell$stage))
order_df = arrange(order_df, seurat_cluster, pt, stage)
mat_new_order = mat[,rownames(order_df)]

seurat_cluster_colors = colorRampPalette(palette_color)(length(unique(order_df$seurat_cluster)))
names(seurat_cluster_colors) = unique(order_df$seurat_cluster)
patient_colors = colorRampPalette(palette_color)(length(unique(integrated_cancer_cell$patient)))
names(patient_colors) = unique(integrated_cancer_cell$patient)
stage_colors = colorRampPalette(c("red", "green"))(length(unique(integrated_cancer_cell$stage)))
names(stage_colors) = unique(integrated_cancer_cell$stage)
top_anno = HeatmapAnnotation(
  patient = order_df$pt,
  stage = order_df$stage,
  seurat_cluster = order_df$seurat_cluster,
  simple_anno_size = unit(1, "cm"),
  show_annotation_name = FALSE,
  col = list(patient = patient_colors,
             stage = stage_colors,
             seurat_cluster = seurat_cluster_colors),
  annotation_legend_param = list(
    seurat_cluster = list(
      grid_height = unit(6, "mm"),
      grid_width = unit(6, "mm"),
      labels_gp = gpar(fontsize = 20),
      title_gp = gpar(fontsize = 20, fontface = "bold")
      # title_gap = unit(20, "mm"),
      # title_position = "lefttop-rot",
    ),
    patient = list(
      grid_height = unit(6, "mm"),
      grid_width = unit(6, "mm"),
      labels_gp = gpar(fontsize = 20),
      title_gp = gpar(fontsize = 20, fontface = "bold")
      # title_gap = unit(20, "mm"),
      # title_position = "lefttop-rot",
    ),
    stage = list(
      grid_height = unit(6, "mm"),
      grid_width = unit(6, "mm"),
      labels_gp = gpar(fontsize = 20),
      title_gp = gpar(fontsize = 20, fontface = "bold")
      # title_gap = unit(20, "mm"),
      # title_position = "lefttop-rot",
    )
  )
)
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "black", "yellow"))

ht = Heatmap(mat_new_order, top_annotation = top_anno, cluster_rows = F, cluster_columns = F,
             show_row_names = T, show_column_names = F, col = col_fun, name = "Expression"
             , column_gap = unit(5, "mm")
             #,column_order = rownames(order_df)
             , column_split = order_df[,1],
             row_names_gp = gpar(fontsize = 25, fontface = "bold"),
             heatmap_legend_param = list(
               labels_gp = gpar(fontsize = 20),
               title_gp = gpar(fontsize = 20, fontface = "bold"),
               direction = "horizontal"
             )
)

png("integrated_tumor_cell_heatmap.png", 30,16, units = "in", res = 300)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

DimPlot(integrated_cancer_cell, reduction = "umap", group.by = "integrated_snn_res.0.3")

test = integrated_cancer_cell













cohort = temp_cancer_cell
selected_hla = hla_markers[3]
selected_mrp = mrp_markers[4]
selected_ndufa = ndufa_markers[3]
cohort$expression_cluster = apply(cohort[["RNA"]]@data, 2, function(x){
  if(x[selected_hla][[1]] > 0 & (x[selected_mrp][[1]] <= 0 & x[selected_ndufa][[1]] <= 0)){
    return("HLA_high")
  }else if(x[selected_hla][[1]] <= 0 & (x[selected_mrp][[1]] > 0 & x[selected_ndufa][[1]] > 0)){
    return("MrpNdufa_high")
  }else{
    return("Others")
  }
})

mat = cohort[["RNA"]]@data
mat = mat[selected_marker, ]
mat = data.frame(mat, check.names = F)
order_df = data.frame(cbind(pt = cohort$patient, class = cohort$expression_cluster, stage = cohort$stage))
order_df = arrange(order_df, pt, class)
mat_new_order = mat[,rownames(order_df)]
cohort_sub = subset(cohort, subset = expression_cluster != "Others")

patient_colors = colorRampPalette(palette_color)(length(unique(cohort$patient)))
names(patient_colors) = unique(cohort$patient)
expression_cluster_colors = colorRampPalette(c("blue", "purple", "pink"))(length(unique(cohort$expression_cluster)))
names(expression_cluster_colors) = unique(cohort$expression_cluster)
stage_colors = colorRampPalette(c("red", "green"))(length(unique(cohort$stage)))
names(stage_colors) = unique(cohort$stage)
top_anno = HeatmapAnnotation(
  patient = order_df$pt,
  expression_cluster = order_df$class,
  stage = order_df$stage,
  simple_anno_size = unit(1, "cm"),
  show_annotation_name = FALSE,
  col = list(patient = patient_colors, expression_cluster = expression_cluster_colors,
             stage = stage_colors),
  annotation_legend_param = list(
    patient = list(
      grid_height = unit(6, "mm"),
      grid_width = unit(6, "mm"),
      labels_gp = gpar(fontsize = 20),
      title_gp = gpar(fontsize = 20, fontface = "bold")
      # title_gap = unit(20, "mm"),
      # title_position = "lefttop-rot",
    ),
    expression_cluster = list(
      grid_height = unit(6, "mm"),
      grid_width = unit(6, "mm"),
      labels_gp = gpar(fontsize = 20),
      title_gp = gpar(fontsize = 20, fontface = "bold")
      # title_gap = unit(20, "mm"),
      # title_position = "lefttop-rot",
    ),
    stage = list(
      grid_height = unit(6, "mm"),
      grid_width = unit(6, "mm"),
      labels_gp = gpar(fontsize = 20),
      title_gp = gpar(fontsize = 20, fontface = "bold")
      # title_gap = unit(20, "mm"),
      # title_position = "lefttop-rot",
    )
  )
)
col_fun = colorRamp2(c(0, 1, 2), c("blue", "black", "yellow"))

png("test.png", 30,16, units = "in", res = 300)
Heatmap(mat_new_order, top_annotation = top_anno, cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F, col = col_fun, name = "Expression"
        #,column_gap = unit(3, "mm")
        #,column_order = rownames(order_df)
        #,column_split = order_df[,1:2]
)
dev.off()


