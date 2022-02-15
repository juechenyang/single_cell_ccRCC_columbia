source("call_libraries.R")
rm(list=setdiff(ls(), "patient_raw_list"))

#define markers
hla_markers = c("HLA-DQB2","HLA-DOA","HLA-DRA","HLA-DMB","HLA-DPB1","HLA-DQA2")
mrp_markers = c("MRPL38","MRPS10","MRPL12","MRPL14","MRPL47","MRPS18A","MRPS24","MRPS36","MRPS5")
ndufa_markers = c("NDUFAB1", "NDUFAF3", "NDUFA5", "NDUFS3","ATP5MC1","ATP5MC2","COA3","UQCRC1")
cancer_markers = c("CA9","NDUFA4L2","SLC17A3")

#read raw data
# all_raw_data=readRDS("./raw_stratified_rds/all_raw.rds")
# patient_raw_list = SplitObject(all_raw_data, split.by = "patient")
patient_raw_list = readRDS("./raw_stratified_rds/patient_raw_list.rds")
reference_pt = readRDS("reference_pt_cells.rds")
quantile_threshold = 0.95
prediction_score_threshold = 0.8
######################################################
# Infercnv reference prepare
######################################################
# using high confident(>=0.8) PT cells from normal sample as reference
# prepare reference cells data
# neg_normal_list = list()
# for(x in 1:8){
#   single_obj = SplitObject(patient_raw_list[[x]], split.by = "cd45")[["CD45-"]]
#   single_obj = SplitObject(single_obj, split.by = "tissue")[["Normal"]]
#   single_obj = PercentageFeatureSet(single_obj, pattern = "^MT-", col.name = "percent.mt")
#   single_obj = subset(single_obj, subset = percent.mt <= 10 & nCount_RNA >= 1000 &
#                                     nCount_RNA <= quantile(single_obj$nCount_RNA, quantile_threshold) &
#                                     nFeature_RNA <= quantile(single_obj$nFeature_RNA, quantile_threshold) &
#                                     nFeature_RNA >= 500)
#   single_obj_azimuth = run_azimuth(single_obj)
#   single_obj$individual_anno = as.character(single_obj_azimuth$predicted.annotation.l1)
#   single_obj$predicted.annotation.l1.score = as.numeric(single_obj_azimuth$predicted.annotation.l1.score)
#   neg_normal_list[[x]] = single_obj
# }
# neg_normal_combine = merge(neg_normal_list[[1]], unlist(neg_normal_list[2:8]))
# neg_normal_combine_pt = subset(neg_normal_combine, subset = individual_anno=="Proximal Tubule"
#                                & predicted.annotation.l1.score >= prediction_score_threshold)

filtered_patient_seurat = list()
for(patient_number in 1:8){
  #patient_number = 4
  cd45_status = "CD45-"
  tissue = "Tumor"
  patient_seurat = patient_raw_list[[patient_number]]
  patient_seurat = SplitObject(patient_seurat, split.by = "cd45")[[cd45_status]]
  
  
  patient_seurat = SplitObject(patient_seurat, split.by = "tissue")[[tissue]]
  #add MT features and do QC
  patient_seurat = PercentageFeatureSet(patient_seurat, pattern = "^MT-", col.name = "percent.mt")
  patient_seurat = subset(patient_seurat, subset = percent.mt <= 10 & nCount_RNA >= 1000 & 
                            nCount_RNA <= quantile(patient_seurat$nCount_RNA, quantile_threshold) &
                            nFeature_RNA <= quantile(patient_seurat$nFeature_RNA, quantile_threshold) &
                            nFeature_RNA >= 500)
  
  #add azimuth results
  patient_seurat_azimuth = run_azimuth(patient_seurat)
  #process with standard pipeline
  DefaultAssay(patient_seurat) = "RNA"
  patient_seurat$individual_anno = as.character(patient_seurat_azimuth$predicted.annotation.l1)
  patient_seurat$predicted.annotation.l1.score = as.numeric(patient_seurat_azimuth$predicted.annotation.l1.score)
  patient_seurat = preprocess_standard_pipelines(patient_seurat)
  patient_seurat = run_pca_and_umap(patient_seurat, n_of_pc = 30)
  patient_seurat = find_neighbours_and_clusters(patient_seurat, n_of_pc = 30)
  patient_seurat = ScaleData(patient_seurat, features = rownames(patient_seurat))
  if(tissue == "Normal"){
    patient_seurat$cluster = as.character(patient_seurat$RNA_snn_res.0.1)
    want = table(patient_seurat$cluster)
    want = want[want>=50]
    patient_seurat$cluster_new = sapply(patient_seurat$cluster, function(x){
      if(x %in% names(want)){
        return(x)
      }else{
        return("Miscellaneous")
      }
    })
    
    selected_hla = hla_markers[3]
    selected_mrp = mrp_markers[2]
    selected_ndufa = ndufa_markers[1]
    patient_seurat$expression_cluster = apply(patient_seurat[["RNA"]]@scale.data, 2, function(x){
      if(x[selected_hla][[1]] > 0 & (x[selected_mrp][[1]] <= 0 & x[selected_ndufa][[1]] <= 0)){
        return("HLA_high")
      }else if(x[selected_hla][[1]] <= 0 & (x[selected_mrp][[1]] > 0 | x[selected_ndufa][[1]] > 0)){
        return("MrpNdufa_high")
      }else if(x[selected_hla][[1]] > 0 & (x[selected_mrp][[1]] > 0 | x[selected_ndufa][[1]] > 0)){
        return("Both_high")
      }else if(x[selected_hla][[1]] <= 0 & (x[selected_mrp][[1]] <= 0 & x[selected_ndufa][[1]] <= 0)){
        return("Both_low")
      }
    })
    
    png(paste0("Composite_marker_heatmap_Patient_normal", patient_number, ".png"),20,15, units = "in", res = 600)
    DoHeatmap(patient_seurat, features = c(hla_markers[3], mrp_markers[2], ndufa_markers[1], cancer_markers), 
              group.by = "RNA_snn_res.0.4")
    dev.off()
  }else{
    print("running InferCnv")
    ######################################################
    # Run Infercnv
    ######################################################
    # produce copy number change heatmap for tumor data(excluding unlikely cancer cell)
    #read refernce pt cells
    # patient_seurat$individual_anno = as.character(patient_seurat$predicted.annotation.l1)
    patient_seurat$epi_marker = 
      ifelse(grepl("Endothelial|Vascular Smooth", 
                   patient_seurat$individual_anno), "yes", "no")
    patient_seurat = subset(patient_seurat, subset = epi_marker == "no")
    patient = merge(patient_seurat, reference_pt)
    want = table(patient$individual_anno)
    want = want[want>=50]
    patient$individual_anno_new = sapply(patient$individual_anno, function(x){
      if(x %in% names(want)){
        return(x)
      }else{
        return("Miscellaneous")
      }
    })
    
    InferCnv_output_dir = paste0("./infercnv_output/",cd45_status,"_",tissue,"_","Patient_", patient_number, "/")
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
    
    
    print("filtering cells based on InferCnv Output")
    gene_position = read.csv("human_gene_pos.csv", header = T)
    gene_position_chr3 = gene_position[gene_position$seqname=="chr3" & gene_position$arm == "p",]
    patinet_infercnv_data = read.table(InferCnv_output_exist_marker, sep="", header = T)
    patinet_infercnv_data_chr3 = patinet_infercnv_data[rownames(patinet_infercnv_data) %in% gene_position_chr3$gene_name,]
    print(length(rownames(patinet_infercnv_data_chr3)))
    transformed_patinet_infercnv_data_chr3 = na.omit(data.frame(t(patinet_infercnv_data_chr3)))
    infercnv_signal_threshold = 0.985
    transformed_patinet_infercnv_data_chr3$keep = apply(transformed_patinet_infercnv_data_chr3, 1, function(x){
      less_than_threshold_percent = length(x[x<=infercnv_signal_threshold])/length(x)
      if (less_than_threshold_percent >= 0.2){
        return("yes")
      }else{
        return("no")
      }
    })
    transformed_filtered_patinet_infercnv_data = transformed_patinet_infercnv_data_chr3[transformed_patinet_infercnv_data_chr3$keep == "yes",]
    print(nrow(transformed_filtered_patinet_infercnv_data))
    print(nrow(transformed_filtered_patinet_infercnv_data)/nrow(transformed_patinet_infercnv_data_chr3))
    selected_cells = rownames(transformed_filtered_patinet_infercnv_data)
    selected_cells = sapply(selected_cells, function(x){
      stringr::str_replace_all(x, "\\.", "-")
    })
    
    patient_seurat$select = sapply(names(patient_seurat$orig.ident), function(x){
      if(x %in% selected_cells){
        return("RCC")
      }else{
        return("Normal_Cells")
      }
    })
    
    filtered_patient_seurat[[patient_number]] = patient_seurat
  }
}

filtered_patient_seurat = readRDS("filtered_patient_seurat.rds")
cohort = merge(filtered_patient_seurat[[1]],filtered_patient_seurat[2:8])
all_types = as.character(unique(cohort$individual_anno))
all_types = c("RCC", all_types,  "Miscellaneous")
palette_color = c("chocolate1", "cyan", "gold", "aquamarine", "deepskyblue", 
                   "darkblue", "darkolivegreen1", "darkorchid1","cyan4",
                  "firebrick1", "firebrick4", "purple", "darksalmon","darkslategray", "darkmagenta")
all_colors = colorRampPalette(palette_color)(length(all_types))
names(all_colors) = all_types

patient_seurat_list = list()
cancer_cells_patient_seurat_list = list()

for(patient_number in 1:8){
  patient_seurat = filtered_patient_seurat[[patient_number]]
  want = table(patient_seurat$individual_anno)
  want = want[want>=50]
  patient_seurat$individual_anno = sapply(patient_seurat$individual_anno, function(x){
    if(x %in% names(want)){
      return(x)
    }else{
      return("Miscellaneous")
    }
  })
  ca9 = patient_seurat[["RNA"]]@data["CA9",]
  selected_ca9 = names(ca9[ca9>0])
  selected_ch3_loss = names(patient_seurat$select[patient_seurat$select=="RCC"])
  overlap = intersect(selected_ca9, selected_ch3_loss)
  patient_seurat$individual_anno_new = sapply(names(patient_seurat$orig.ident), function(x){
    if(x %in% overlap){
      return("RCC")
    }else{
      return(as.character(patient_seurat$individual_anno[x]))
    }
  })
  cancer_cells_patient_seurat = subset(patient_seurat, subset = individual_anno_new == "RCC")
  selected_hla = hla_markers[3]
  selected_mrp = mrp_markers[4]
  selected_ndufa = ndufa_markers[3]
  cancer_cells_patient_seurat$expression_cluster = apply(cancer_cells_patient_seurat[["RNA"]]@data, 2, function(x){
    if(x[selected_hla][[1]] > 0 & (x[selected_mrp][[1]] <= 0 & x[selected_ndufa][[1]] <= 0)){
      return("HLA_high")
    }else if(x[selected_hla][[1]] <= 0 & (x[selected_mrp][[1]] > 0 & x[selected_ndufa][[1]] > 0)){
      return("MrpNdufa_high")
    }else{
      return("Others")
    }
  })
  cancer_cells_patient_seurat_list[[patient_number]] = cancer_cells_patient_seurat
  patient_seurat_list[[patient_number]] = patient_seurat
}


umap_plotlist = list()
vln_plotlist = list()
for(patient_number in 1:8){
  patient_seurat = patient_seurat_list[[patient_number]]
  cell_types = as.character(patient_seurat$individual_anno_new)
  p = DimPlot(patient_seurat, reduction = "umap", group.by = "individual_anno_new")+
    scale_color_manual(breaks = cell_types, 
                       values=all_colors[cell_types])+ggtitle(unique(patient_seurat$patient))
  umap_plotlist[[patient_number]] = p
  p1 = VlnPlot(patient_seurat, features = "CA9", group.by = "individual_anno_new")
  vln_plotlist[[patient_number]] = p1
}
png("umap_analysis.png",25,16, units = "in", res = 600)
ggarrange(
  plotlist = umap_plotlist,
  nrow = 3,
  ncol = 3
)
dev.off()
png("vln_analysis.png",25,16, units = "in", res = 600)
ggarrange(
  plotlist = vln_plotlist,
  nrow = 3,
  ncol = 3
)
dev.off()


cohort = merge(cancer_cells_patient_seurat_list[[1]],cancer_cells_patient_seurat_list[2:8])
scale_data = matrix(nrow = 19234, ncol = 1)
for(x in 1:8){
  scale_data = cbind(scale_data, cancer_cells_patient_seurat_list[[x]][["RNA"]]@scale.data)
}
# cohort = Integrate_seurat_obj_by_key(cohort, key = "patient", method_SCT = F)
# cohort = run_pca_and_umap(cohort)
# cohort = find_neighbours_and_clusters(cohort)
# DimPlot(cohort, reduction = "umap", group.by = "patient")
# DimPlot(cohort, reduction = "umap", group.by = "integrated_snn_res.0.2")
#cohort = ScaleData(cohort, features = rownames(cohort))
cohort[["RNA"]]@scale.data = scale_data
mat = cohort[["RNA"]]@scale.data
selected_marker = c(hla_markers, mrp_markers, ndufa_markers, cancer_markers)
# selected_markers_tag = sapply(selected_marker, function(x){
#   zero_count = as.numeric(table(mat[x,]==0)["TRUE"])
#   if(is.na(zero_count)){
#     return(TRUE)
#   }else{
#     zero_percent = zero_count/length(mat[x,])
#     if(zero_percent >= 0.8){
#       return(FALSE)
#     }else{
#       return(TRUE)
#     }
#   }
# })
# selected_marker = selected_marker[selected_markers_tag]
selected_marker = c("HLA-DRA", "HLA-DPB1", "MRPL14", "MRPS36", "NDUFAB1", "NDUFA5", "ATP5MC1", "ATP5MC2", "COA3")
#mat = cohort[["RNA"]]@data
mat = mat[selected_marker, ]
mat = data.frame(mat, check.names = F)
order_df = data.frame(cbind(pt = cohort$patient, class = cohort$expression_cluster, stage = cohort$stage))
order_df = arrange(order_df, class, pt)
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

png("test2.png", 30,16, units = "in", res = 300)
Heatmap(mat_new_order, top_annotation = top_anno, cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F, col = col_fun, name = "Expression"
        #,column_gap = unit(3, "mm")
        #,column_order = rownames(order_df)
        #,column_split = order_df[,1:2]
)
dev.off()









