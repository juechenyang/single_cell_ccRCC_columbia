source("call_libraries.R")
source("Initialization.R")
#load data
integrated_cancer_cell = readRDS("integrated_cancer_cells.rds")


#######################################################
#Umaps
#######################################################

#subclustering_umap
all_types = as.character(unique(integrated_cancer_cell$integrated_snn_res.0.2))
all_types = str_sort(all_types)
palette_color = c("chocolate1", "cyan", "gold", "aquamarine", "deepskyblue",
                  "cyan4", "darkblue", "darksalmon", "darkorchid1",
                  "goldenrod4", "firebrick4", "firebrick1", "darkolivegreen1","darkslategray", "darkmagenta")
all_colors = colorRampPalette(palette_color)(length(all_types))
# all_colors = colors()[seq(18,(length(all_types)-1)*6+18, 6)]
names(all_colors) = all_types
png("subclustering_umap.png",9,9, units = "in", res = 300)
DimPlot(integrated_cancer_cell, reduction = "umap", 
        group.by = "integrated_snn_res.0.2")+
  scale_color_manual(breaks = all_types, values=all_colors[all_types])+ggtitle("Louvain cluster")
dev.off()





DefaultAssay(integrated_cancer_cell) = "RNA"
integrated_cancer_cell = NormalizeData(integrated_cancer_cell)
marker_genes = c("VHL", "PBRM1","BAP1", "CA9", "HIF1A", "EPAS1")
png("marker_genes_expression_umap.png",24,12, units = "in", res = 300)
FeaturePlot(integrated_cancer_cell, features = marker_genes, ncol = 3)
dev.off()



hla_selected_markers = c("HLA-DRA", "HLA-DPB1")
mrp_selected_markers = c("MRPS10", "MRPL12", "MRPL14", "MRPS36", "MRPS5")
nduf_selected_markers = c("NDUFAB1", "NDUFAF3", "NDUFA5", "NDUFS3", 
                          "ATP5MC1","ATP5MC2","COA3","UQCRC1")
# hla_markers = hla_selected_markers
# mrp_markers = mrp_selected_markers
# nduf_markers = nduf_selected_markers

integrated_cancer_cell = AddModuleScore(integrated_cancer_cell, 
                                        features = list(hla_markers, 
                                                        mrp_markers, 
                                                        nduf_markers, 
                                                        c(mrp_markers, nduf_markers), 
                                                        EMT_markers,
                                                        PT_markers), 
                                        assay = "RNA", 
                                        name = c("hla_score", "mrp_score", "nduf_score", 
                                                 "mrp_nduf_score", "emt_score", "pt_score"))
integrated_cancer_cell@meta.data[,c("hla_score", "mrp_score", "nduf_score",
                                    "mrp_nduf_score", "emt_score", "pt_score")] = 
integrated_cancer_cell@meta.data[,c("hla_score1", "mrp_score2", "nduf_score3", 
                                      "mrp_nduf_score4", "emt_score5", "pt_score6")]

p1 = DimPlot(integrated_cancer_cell, reduction = "umap", group.by = "integrated_snn_res.0.2")+
  scale_color_manual(breaks = all_types, values=all_colors[all_types])+ggtitle("Louvain cluster")
p2 = DimPlot(integrated_cancer_cell, reduction = "umap", group.by = "patient")
p3 = FeaturePlot(integrated_cancer_cell, features = c("hla_score", "mrp_score", 
                                                      "nduf_score", "mrp_nduf_score", 
                                                      "emt_score","pt_score"), ncol = 3, min.cutoff = -0.5, max.cutoff = 0.5)

png("module_score_umaps_all_g23_markers.png",24,12, units = "in", res = 300)
ggarrange(
  ggarrange(p1, p2, nrow = 2), 
  p3,
  ncol = 2,
  widths = c(1,3)
)
dev.off()





#######################################################
#Heatmaps
#######################################################

# module score heatmap
mat = t(integrated_cancer_cell@meta.data[,c("hla_score", "mrp_score", "nduf_score", 
                                            "mrp_nduf_score", "emt_score", "pt_score")])
mat = data.frame(mat, check.names = F)
order_df = data.frame(cbind(seurat_cluster = as.character(integrated_cancer_cell$integrated_snn_res.0.2),
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
  #show_annotation_name = FALSE,
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

png("module_score_all_g23_genes_heatmap.png", 25,12, units = "in", res = 300)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "right")
dev.off()


# gene expression heatmap
DefaultAssay(integrated_cancer_cell) = "RNA"
integrated_cancer_cell = NormalizeData(integrated_cancer_cell)
cohort = integrated_cancer_cell
mat = cohort[["RNA"]]@data
markers = c(hla_markers, mrp_markers, nduf_markers)
markers_in_data = markers %in% rownames(mat)
markers = markers[markers_in_data]
mat = mat[markers, ]
mat = data.frame(mat, check.names = F)
order_df = data.frame(cbind(pt = cohort$patient, stage = cohort$stage))
order_df = arrange(order_df, pt)
mat_new_order = mat[,rownames(order_df)]

patient_colors = colorRampPalette(palette_color)(length(unique(cohort$patient)))
names(patient_colors) = unique(cohort$patient)
stage_colors = colorRampPalette(c("red", "green"))(length(unique(cohort$stage)))
names(stage_colors) = unique(cohort$stage)
top_anno = HeatmapAnnotation(
  patient = order_df$pt,
  stage = order_df$stage,
  simple_anno_size = unit(1, "cm"),
  show_annotation_name = FALSE,
  col = list(patient = patient_colors, stage = stage_colors),
  annotation_legend_param = list(
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
col_fun = colorRamp2(c(0, 1, 2), c("blue", "black", "yellow"))

png("gene_expression_all_g23_heatmap.png", 25,16, units = "in", res = 300)
Heatmap(mat_new_order, top_annotation = top_anno, cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F, col = col_fun, name = "Expression"
        #,column_gap = unit(3, "mm")
        #,column_order = rownames(order_df)
        #,column_split = order_df[,1:2]
)
dev.off()


