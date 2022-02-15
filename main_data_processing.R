

rm(list=setdiff(ls(), c("neg_tumor")))

all_raw_data=readRDS("./raw_stratified_rds/all_raw.rds")
# #split by cd45
neg = SplitObject(all_raw_data, split.by = "cd45")[[2]]
#QC
neg <- PercentageFeatureSet(neg, pattern = "^MT-", col.name = "percent.mt")
neg = subset(neg, subset = percent.mt <= 10 & nCount_RNA >= 1000 & nCount_RNA <= 30000 &
                        nFeature_RNA <= 5000 & nFeature_RNA >= 500)

#get negative tumor
neg_tumor = SplitObject(neg, split.by = "tissue")[[1]]
#saveRDS(neg_tumor, "neg_tumor.rds")
# get negative normal
neg_normal = SplitObject(neg, split.by = "tissue")[[2]]
#saveRDS(neg_normal, "neg_normal.rds")


neg_tumor = readRDS("neg_tumor.rds")
number_of_component = 50
number_of_variable_features = 2500
number_of_anchors = 5
#Add azimuth results back
neg_tumor_azimuth = AddAzimuthResults(neg_tumor, "neg_tumor_azimuth_results.Rds")

####### Build color code
all_anno_individual = c()
for(x in seq(1,8,1)){
  patient = readRDS(paste0("./raw_stratified_rds/Patient", x, "_neg_tumor.rds"))
  patient = AddAzimuthResults(patient, paste0("./raw_stratified_rds/Patient", x, "_azimuth_results.Rds"))
  all_anno_individual = c(all_anno_individual, as.character(patient$predicted.annotation.l1))
}
neg_tumor_azimuth$individual_anno = all_anno_individual
all_types = as.character(unique(neg_tumor_azimuth$individual_anno))
palette_color = c("chocolate1", "cyan", "gold", "aquamarine", "deepskyblue", 
                  "cyan4", "darkblue", "darkolivegreen1", "darkorchid1",
                  "firebrick1", "firebrick4", "yellow", "darksalmon","darkslategray", "darkmagenta")
all_colors = colorRampPalette(palette_color)(length(all_types))
names(all_colors) = all_types


#############################################################
#full cell population
#############################################################

##############################without integration###############################
#process neg tumor with SCT pipeline
neg_tumor_azimuth_SCT = SCTransform(neg_tumor_azimuth, variable.features.n = 
                            number_of_variable_features, 
                            vars.to.regress = "percent.mt")
neg_tumor_azimuth_SCT = run_pca_and_umap(neg_tumor_azimuth_SCT, n_of_pc = number_of_component)
#neg tumor with standard pipeline
neg_tumor_azimuth_standard = preprocess_standard_pipelines(neg_tumor_azimuth, 
                                                  n_variable_features = number_of_variable_features)
neg_tumor_azimuth_standard = run_pca_and_umap(neg_tumor_azimuth_standard, n_of_pc = number_of_component)

# saveRDS(neg_tumor_azimuth_standard, "neg_tumor_azimuth_standard.rds")
# saveRDS(neg_tumor_azimuth_SCT, "neg_tumor_azimuth_SCT.rds")

#plot for azimuth projection
png("neg_tumor_AzimuthUmap_AzimuthL1.png",16,9, units = "in", res = 600)
DimPlot(neg_tumor_azimuth, reduction = "umap.proj", group.by = "predicted.annotation.l1", 
        pt.size = 0.2, label = F)
dev.off()

##############################with integration###############################
#integrate tumor data separated by patients with standard pipeline
neg_tumor_azimuth_integrated_by_patient_standard = Integrate_seurat_obj_by_key(neg_tumor_azimuth, key = "patient", 
                                                                n.features = number_of_variable_features, 
                                                                n.anchors = number_of_anchors,
                                                                method_SCT = F)
#run unmap on this integrated data
neg_tumor_azimuth_integrated_by_patient_standard = run_pca_and_umap(neg_tumor_azimuth_integrated_by_patient_standard,
                                                     n_of_pc = number_of_component)
#integrate tumor data separated by patients with SCT
neg_tumor_azimuth_integrated_by_patient_SCT = Integrate_seurat_obj_by_key(neg_tumor_azimuth, key = "patient", 
                                                                n.features = number_of_variable_features, n.anchors = 5,
                                                                method_SCT = T)
#run unmap on this integrated data
neg_tumor_azimuth_integrated_by_patient_SCT = run_pca_and_umap(neg_tumor_azimuth_integrated_by_patient_SCT,
                                                     n_of_pc = number_of_component)

# saveRDS(neg_tumor_azimuth_integrated_by_patient_standard, 
#         "neg_tumor_azimuth_integrated_by_patient_standard.rds")
# saveRDS(neg_tumor_azimuth_integrated_by_patient_SCT, 
#         "neg_tumor_azimuth_integrated_by_patient_SCT.rds")




####################################################
# remove cell types that are unlikely cancer cell
####################################################

#remove cells that are unlikely cancer cells
neg_tumor_azimuth$epi_marker = 
        ifelse(grepl("Endothelial|Fibroblast|Vascular Smooth|Immune", 
                     neg_tumor_azimuth$predicted.annotation.l1), "yes", "no")
neg_tumor_epi_azimuth = subset(neg_tumor_azimuth, subset = epi_marker == "no")

##############################without integration###############################
#process neg tumor with SCT pipeline
neg_tumor_epi_azimuth_SCT = SCTransform(neg_tumor_epi_azimuth, variable.features.n = 
                                            number_of_variable_features, 
                                    vars.to.regress = "percent.mt")
neg_tumor_epi_azimuth_SCT = run_pca_and_umap(neg_tumor_epi_azimuth_SCT, n_of_pc = number_of_component)
#neg tumor with standard pipeline
neg_tumor_epi_azimuth_standard = preprocess_standard_pipelines(neg_tumor_epi_azimuth, 
                                                           n_variable_features = number_of_variable_features)
neg_tumor_epi_azimuth_standard = run_pca_and_umap(neg_tumor_epi_azimuth_standard, 
                                                  n_of_pc = number_of_component)

saveRDS(neg_tumor_epi_azimuth_SCT, "neg_tumor_epi_azimuth_SCT.rds")
saveRDS(neg_tumor_epi_azimuth_standard, "neg_tumor_epi_azimuth_standard.rds")


##############################with integration##################################
#integrate tumor data separated by patients with standard pipeline
neg_tumor_epi_azimuth_integrated_by_patient_standard = Integrate_seurat_obj_by_key(neg_tumor_epi_azimuth, key = "patient", 
                                                                n.features = number_of_variable_features, 
                                                                n.anchors = number_of_anchors,
                                                                method_SCT = F)
#run unmap on this integrated data
neg_tumor_epi_azimuth_integrated_by_patient_standard = run_pca_and_umap(neg_tumor_epi_azimuth_integrated_by_patient_standard,
                                                     n_of_pc = number_of_component)
#integrate tumor data separated by patients with SCT
neg_tumor_epi_azimuth_integrated_by_patient_SCT = Integrate_seurat_obj_by_key(neg_tumor_epi_azimuth, key = "patient", 
                                                           n.features = number_of_variable_features, n.anchors = 5,
                                                           method_SCT = T)


#run unmap on this integrated data
neg_tumor_epi_azimuth_integrated_by_patient_SCT = run_pca_and_umap(neg_tumor_epi_azimuth_integrated_by_patient_SCT,
                                                n_of_pc = number_of_component)


saveRDS(neg_tumor_epi_azimuth_integrated_by_patient_standard, 
        "neg_tumor_epi_azimuth_integrated_by_patient_standard.rds")
saveRDS(neg_tumor_epi_azimuth_integrated_by_patient_SCT, 
        "neg_tumor_epi_azimuth_integrated_by_patient_SCT.rds")




DefaultAssay(integrated_data_patients_standard) = "RNA"
png("feature_plot_CA9.png",16,9, units = "in", res = 600)
FeaturePlot(integrated_data_patients_standard, features = c("CA9","NDUFA4L2"), slot = "data")
dev.off()
##############################end############################################################
##########################################################################################
##########################################################################################
##########################################################################################




#########################################################
# patient level analysis include all cells(tumor data)
#########################################################
neg_tumor_patient_list = SplitObject(neg_tumor_azimuth, split.by = "patient")
neg_tumor_patient_list_standard = list()
neg_tumor_patient_list_SCT = list()
for(i in 1:length(neg_tumor_patient_list)){
   single_patient_obj = neg_tumor_patient_list[[i]]
   patient_id = unique(single_patient_obj$patient)
   single_patient_obj = PercentageFeatureSet(single_patient_obj, pattern = "^MT-", col.name = "percent.mt")
   single_patient_obj_SCT = SCTransform(single_patient_obj, variable.features.n = 
                                                   number_of_variable_features, 
                                           vars.to.regress = "percent.mt")
   single_patient_obj_SCT = run_pca_and_umap(single_patient_obj_SCT, n_of_pc = number_of_component)
   single_patient_obj_standard = preprocess_standard_pipelines(single_patient_obj, 
                                n_variable_features = number_of_variable_features)
   single_patient_obj_standard = run_pca_and_umap(single_patient_obj_standard, 
                                                     n_of_pc = number_of_component)
   neg_tumor_patient_list_SCT[[i]] = single_patient_obj_SCT
   neg_tumor_patient_list_standard[[i]] = single_patient_obj_standard
}
saveRDS(neg_tumor_patient_list_SCT, "neg_tumor_patient_list_SCT.rds")
saveRDS(neg_tumor_patient_list_standard, "neg_tumor_patient_list_standard.rds")

#########################################################
# patient level analysis keeps likely cancer cells only
#########################################################
neg_tumor_azimuth$epi_marker = 
        ifelse(grepl("Endothelial|Fibroblast|Vascular Smooth|Immune", 
                     neg_tumor_azimuth$individual_anno), "yes", "no")
neg_tumor_epi_azimuth = subset(neg_tumor_azimuth, subset = epi_marker == "no")

neg_tumor_epi_patient_list = SplitObject(neg_tumor_epi_azimuth, split.by = "patient")
neg_tumor_epi_patient_list_standard = list()
neg_tumor_epi_patient_list_SCT = list()
for(i in 1:length(neg_tumor_epi_patient_list)){
        single_patient_obj = neg_tumor_epi_patient_list[[i]]
        patient_id = unique(single_patient_obj$patient)
        single_patient_obj = PercentageFeatureSet(single_patient_obj, pattern = "^MT-", col.name = "percent.mt")
        single_patient_obj_SCT = SCTransform(single_patient_obj, variable.features.n = 
                                                     number_of_variable_features, 
                                             vars.to.regress = "percent.mt")
        single_patient_obj_SCT = run_pca_and_umap(single_patient_obj_SCT, n_of_pc = number_of_component)
        single_patient_obj_standard = preprocess_standard_pipelines(single_patient_obj, 
                                                                    n_variable_features = number_of_variable_features)
        single_patient_obj_standard = run_pca_and_umap(single_patient_obj_standard, 
                                                       n_of_pc = number_of_component)
        neg_tumor_epi_patient_list_SCT[[i]] = single_patient_obj_SCT
        neg_tumor_epi_patient_list_standard[[i]] = single_patient_obj_standard
}
saveRDS(neg_tumor_epi_patient_list_SCT, "neg_tumor_epi_patient_list_SCT.rds")
saveRDS(neg_tumor_epi_patient_list_standard, "neg_tumor_epi_patient_list_standard.rds")


##########################################################################################
##########################################################################################


######################################################
# Infercnv analysis
######################################################
# using high confident(>=0.8) PT cells from normal sample as reference
source("call_libraries.R")
neg_normal = readRDS("neg_normal.rds")
neg_normal_list = SplitObject(neg_normal, split.by = "patient")
neg_tumor_patient_list = readRDS("neg_tumor_patient_list.rds")

# infercnv analysis data prepare
# prepare reference cells data
for(x in 1:8){
  single_obj = neg_normal_list[[x]]
  single_obj = AddAzimuthResults(single_obj, paste0("./raw_stratified_rds/Patient", x, "_normal_azimuth_results.rds"))
  single_obj$individual_anno = as.character(single_obj$predicted.annotation.l1)
  neg_normal_list[[x]] = single_obj
}
neg_normal_combine = merge(neg_normal_list[[1]], unlist(neg_normal_list[2:8]))
neg_normal_combine_pt = subset(neg_normal_combine, subset = individual_anno=="Proximal Tubule"
                        & predicted.annotation.l1.score >= 0.8)

# produce copy number change heatmap for tumor data(excluding unlikely cancer cell)
for(x in 1:8){
  single_obj_tumor = neg_tumor_patient_list[[x]]
  single_obj_tumor = AddAzimuthResults(single_obj_tumor, paste0("./raw_stratified_rds/Patient", x, "_tumor_azimuth_results.Rds"))
  single_obj_tumor$individual_anno = as.character(single_obj_tumor$predicted.annotation.l1)
  single_obj_tumor$epi_marker = 
    ifelse(grepl("Endothelial|Vascular Smooth", 
                 single_obj_tumor$individual_anno), "yes", "no")
  single_obj_tumor = subset(single_obj_tumor, subset = epi_marker == "no")
  patient = merge(single_obj_tumor, neg_normal_combine_pt)
  want = table(patient$individual_anno)
  want = want[want>=10]
  patient$individual_anno_new = sapply(patient$individual_anno, function(x){
    if(x %in% names(want)){
      return(x)
    }else{
      return("Miscellaneous")
    }
  })
  InfercnvInputs = prepare_infercnv_data(patient, annotation = "individual_anno_new")
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=InfercnvInputs[[1]],
                                      annotations_file=InfercnvInputs[[2]],
                                      delim="\t",
                                      gene_order_file=InfercnvInputs[[3]],
                                      ref_group_names = c("Proximal Tubule")
  )
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir=paste0("./infercnv_output/Patient", x, "/"), 
                               cluster_by_groups=TRUE, denoise=TRUE, 
                               HMM = T, num_threads = 8)
}

# visualization for PT cells marker verification
pt_plotlist = list()
pc_plotlist = list()
for(x in 1:8){
  single_obj_tumor = neg_tumor_patient_list_standard[[x]]
  single_obj_tumor$epi_marker = 
    ifelse(grepl("Endothelial|Vascular Smooth", 
                 single_obj_tumor$individual_anno), "yes", "no")
  single_obj_tumor = subset(single_obj_tumor, subset = epi_marker == "no")
  patient = merge(single_obj_tumor, neg_normal_combine_pt)
  want = table(patient$individual_anno)
  want = want[want>=10]
  patient$individual_anno_new = sapply(patient$individual_anno, function(x){
    if(x %in% names(want)){
      return(x)
    }else{
      return("Miscellaneous")
    }
  })
  patient = Normalize_and_scale(patient)
  proximal_tubular_marker = c("SLC13A3", "SLC16A9")
  pt_plot = VlnPlot(patient, features = proximal_tubular_marker, group.by = "individual_anno_new") + 
    plot_annotation(title = patient$patient) & theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  pt_plotlist[[x]] = pt_plot
  proximal_convoluted_marker = c("SLC17A3", "SLC22A8")
  pc_plot = VlnPlot(patient, features = proximal_convoluted_marker, group.by = "individual_anno_new") + 
    plot_annotation(title = patient$patient) & theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  pc_plotlist[[x]] = pc_plot
}
png("proximal_tubular_marker.png",20,16, units = "in", res = 600)
plt = ggarrange(
  plotlist = pt_plotlist,
  nrow = 3,
  ncol = 3
)
annotate_figure(plt, top = text_grob("Proximal Tubular", 
                                      color = "blue", face = "bold", size = 30))
dev.off()

png("proximal_convoluted_marker.png",20,16, units = "in", res = 600)
plt = ggarrange(
  plotlist = pc_plotlist,
  nrow = 3,
  ncol = 3
)
annotate_figure(plt, top = text_grob("Proximal Convoluted", 
                                     color = "red", face = "bold", size = 30))
dev.off()

#########################################################
# using infercnv results as a filter for 
# more concise cancer cell population
#########################################################
# neg_tumor = readRDS("neg_tumor.rds")
# neg_tumor_patient_list = SplitObject(neg_tumor, split.by = "patient")
# saveRDS(neg_tumor_patient_list, "neg_tumor_patient_list.rds")

#select confident cancer cells with chr3p loss
gene_position = read.csv("human_gene_pos.csv", header = T)
neg_tumor_patient_list = readRDS("neg_tumor_patient_list.rds")
hla_markers = c("HLA-DQB2","HLA-DOA","HLA-DRA","HLA-DMB","HLA-DPB1","HLA-DQA2")
mrp_markers = c("MRPL38","MRPS10","MRPL12","MRPL14","MRPL47","MRPS18A","MRPS24","MRPS36","MRPS5")
ndufa_markers = c("NDUFAB1", "NDUFAF3", "NDUFA5", "NDUFS3","ATP5MC1","ATP5MC2","COA3","UQCRC1")
cancer_markers = c("CA9","NDUFA4L2","SLC17A3")

patient_number = 2
gene_position_chr3 = gene_position[gene_position$seqname=="chr3" & gene_position$arm == "p",]
patinet_infercnv_data = read.table(paste0("./infercnv_output/Patient",
                        patient_number,"/infercnv.observations.txt"), sep="", header = T)
patinet_infercnv_data_chr3 = patinet_infercnv_data[rownames(patinet_infercnv_data) %in% gene_position_chr3$gene_name,]
transformed_patinet_infercnv_data_chr3 = na.omit(data.frame(t(patinet_infercnv_data_chr3)))
infercnv_signal_threshold = 0.95
transformed_patinet_infercnv_data_chr3$keep = apply(transformed_patinet_infercnv_data_chr3, 1, function(x){
  less_than_threshold_percent = length(x[x<=infercnv_signal_threshold])/length(x)
  if (less_than_threshold_percent >= 0.4){
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

#remove unlikely cancer cells and do umap projection
patient = neg_tumor_patient_list[[patient_number]]
patient = preprocess_standard_pipelines(patient, n_variable_features = 3000)
patient = run_pca_and_umap(patient, n_of_pc = 30)
patient = find_neighbours_and_clusters(patient, n_of_pc = 30)
patient = ScaleData(patient, features = rownames(patient))
patient = AddAzimuthResults(patient, paste0("./raw_stratified_rds/Patient", patient_number, "_tumor_azimuth_results.Rds"))
patient$individual_anno = as.character(patient$predicted.annotation.l1)
print(length(Cells(patient)))
patient$epi_marker = 
  ifelse(grepl("Endothelial|Vascular Smooth", 
               patient$individual_anno), "yes", "no")
patient = subset(patient, subset = epi_marker == "no")
print(length(Cells(patient)))

# restrict the seurat obj to selected cells
patient$select = names(patient$orig.ident) %in% selected_cells
patient_sub = subset(patient, subset = select == T)
print(length(Cells(patient_sub)))


patient_sub$cluster = as.character(patient_sub$RNA_snn_res.0.1)
want = table(patient_sub$cluster)
want = want[want>=50]
patient_sub$cluster_new = sapply(patient_sub$cluster, function(x){
  if(x %in% names(want)){
    return(x)
  }else{
    return("Miscellaneous")
  }
})

selected_hla = hla_markers[3]
selected_mrp = mrp_markers[2]
selected_ndufa = ndufa_markers[1]
patient_sub$expression_cluster = apply(patient_sub[["RNA"]]@scale.data, 2, function(x){
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

png(paste0("Composite_marker_heatmap_Patient", patient_number, ".png"),20,15, units = "in", res = 600)
DoHeatmap(patient_sub, features = c(hla_markers[3], mrp_markers[2], ndufa_markers[1], cancer_markers), 
          group.by = "expression_cluster")
dev.off()

p1 = DimPlot(patient_sub, group.by = "expression_cluster", reduction = "umap", label = F)
p2 = DimPlot(patient_sub, group.by = "cluster_new", reduction = "umap", label = F)
png(paste0("Composite_marker_umap_Patient", patient_number, ".png"),20,10, units = "in", res = 300)
ggarrange(
  p1,p2,
  ncol = 2
)
dev.off()


#########################################################
# Build composite marker for cancer cells
#########################################################


neg_tumor_patient_list_standard_clustered = list()
for(x in 1:8){
  single_obj_tumor = neg_tumor_patient_list_standard[[x]]
  single_obj_tumor$epi_marker = 
    ifelse(grepl("Endothelial|Vascular Smooth", 
                 single_obj_tumor$individual_anno), "yes", "no")
  single_obj_tumor = subset(single_obj_tumor, subset = epi_marker == "no")
  single_obj$individual_anno_new = sapply(single_obj$individual_anno,as.character)
  single_obj_tumor = find_neighbours_and_clusters(single_obj_tumor, n_of_pc = 30, res = seq(0.2,1,0.2))
  neg_tumor_patient_list_standard_clustered[[x]] = single_obj_tumor
}


g23_vln_plot_list = list()
cancer_vln_plot_list = list()
for(x in 1:8){
  single_obj_tumor = neg_tumor_patient_list_standard_clustered[[x]]
  g23_vln_plot = VlnPlot(single_obj_tumor, features = c(hla_markers[1], mrp_markers[1], ndufa_markers[2]), group.by = "RNA_snn_res.0.2") + 
    plot_annotation(title = single_obj_tumor$patient) & theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  cancer_vln_plot = VlnPlot(single_obj_tumor, features = cancer_markers, group.by = "RNA_snn_res.0.2") + 
    plot_annotation(title = single_obj_tumor$patient) & theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  g23_vln_plot_list[[x]] = g23_plot
  g23_vln_plot_list[[x]] = cancer_plot
}

feature_plot_list = list()
for(x in 1:8){
  single_obj_tumor = neg_tumor_patient_list_standard_clustered[[x]]
  feature_plot = FeaturePlot(single_obj_tumor, features = c(hla_markers[1], mrp_markers[1], ndufa_markers[3], cancer_markers), 
                             reduction = "umap", ncol = 3)+
  plot_annotation(title = single_obj_tumor$patient) & theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  feature_plot_list[[x]] = feature_plot
}


png("feature_plot_AF3.png",30,15, units = "in", res = 300)
ggarrange(
  plotlist = feature_plot_list,
  nrow = 3,
  ncol = 3
)
dev.off()

png("g23.png",22,16, units = "in", res = 300)
plt = ggarrange(
  plotlist = g23_vln_plot_list,
  nrow = 3,
  ncol = 3
)
annotate_figure(plt, top = text_grob("g23 markers", 
                                     color = "blue", face = "bold", size = 30))
dev.off()

png("cancer.png",20,16, units = "in", res = 300)
plt = ggarrange(
  plotlist = cancer_pltlist,
  nrow = 3,
  ncol = 3
)
annotate_figure(plt, top = text_grob("cancer markers", 
                                     color = "blue", face = "bold", size = 30))
dev.off()











#########################################################
# patient level analysis include all cells(normal data)
#########################################################

neg_normal_patient_list_standard = list()
neg_normal_patient_list_SCT = list()
for(i in 1:length(neg_normal_patient_list)){
  single_patient_obj = neg_normal_patient_list[[i]]
  patient_id = unique(single_patient_obj$patient)
  single_patient_obj = PercentageFeatureSet(single_patient_obj, pattern = "^MT-", col.name = "percent.mt")
  single_patient_obj_SCT = SCTransform(single_patient_obj, variable.features.n = 
                                         number_of_variable_features, 
                                       vars.to.regress = "percent.mt")
  single_patient_obj_SCT = run_pca_and_umap(single_patient_obj_SCT, n_of_pc = number_of_component)
  single_patient_obj_standard = preprocess_standard_pipelines(single_patient_obj, 
                                                              n_variable_features = number_of_variable_features)
  single_patient_obj_standard = run_pca_and_umap(single_patient_obj_standard, 
                                                 n_of_pc = number_of_component)
  neg_normal_patient_list_SCT[[i]] = single_patient_obj_SCT
  neg_normal_patient_list_standard[[i]] = single_patient_obj_standard
}
saveRDS(neg_normal_patient_list_SCT, "neg_normal_patient_list_SCT.rds")
saveRDS(neg_normal_patient_list_standard, "neg_normal_patient_list_standard.rds")













want = table(patient$individual_anno)
want = want[want>=10]
patient$individual_anno_new = sapply(patient$individual_anno, function(x){
  if(x %in% names(want)){
    return(x)
  }else{
    return("Miscellaneous")
  }
})





#using unlikely cancer cells as reference
neg_tumor_patient_list_standard = readRDS("neg_tumor_patient_list_standard.rds")
for(i in 2:8){
  single_obj = neg_tumor_patient_list_standard[[i]]
  want = table(single_obj$individual_anno)
  want = want[want>=10]
  single_obj$individual_anno_new = sapply(single_obj$individual_anno, function(x){
    if(x %in% names(want)){
      return(x)
    }else{
      return("Miscellaneous")
    }
  })
  InfercnvInputs = prepare_infercnv_data(single_obj, annotation = "individual_anno_new")
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=InfercnvInputs[[1]],
                                      annotations_file=InfercnvInputs[[2]],
                                      delim="\t",
                                      gene_order_file=InfercnvInputs[[3]],
                                      ref_group_names = c("Endothelial", "Vascular Smooth Muscle / Pericyte")
  )
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir=paste0("./infercnv_output/Patient", i, "/"), 
                               cluster_by_groups=TRUE, denoise=TRUE, 
                               HMM = T, num_threads = 8)
}


proximal_tubular_marker = c("SLC13A3", "SLC16A9")
DoHeatmap(integrated_epi_azimuth_neg, features = remaining_markers, 
          group.by = "integrated_snn_res.0.15")
DoHeatmap(integrated_epi_azimuth_neg, features = cancer_markers, 
          group.by = "integrated_snn_res.0.1")



#find all markers
Idents(integrated_epi_azimuth_neg) = integrated_epi_azimuth_neg$integrated_snn_res.0.1
cluster_markers <- FindAllMarkers(integrated_epi_azimuth_neg, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,test.use = "MAST")
top <- cluster_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)





#Stage annotations
#"pT3a"-- PatientA, Patient2, Patient3, Patient5(metastatic), Patient7
#"pT1b"-- PatientB, PatientC, Patient1, Patient4, Patient6, Patient8

#Grade annotations
#"1"-- Patient3, Patient4, 
#"2"-- PatientB, PatientC, Patient1, Patient2, Patient6, Patient8
#"3"-- PatientA, Patient7
#"4"-- Patient5
#######

# grade_anno = c("II", "II", "I", "I", "IV", "II", "III", "II", "III",
#                "II", "II")
# names(grade_anno) = c("Patient1", "Patient2", "Patient3", "Patient4", "Patient5",
#                "Patient6", "Patient7", "Patient8", "PatientA", "PatientB",
#                "PatientC")
# stage_anno = c("pT1b","pT3a","pT3a","pT1b","pT3a","pT1b","pT3a","pT1b","pT3a","pT1b","pT1b")
# names(stage_anno) = c("Patient1", "Patient2", "Patient3", "Patient4", "Patient5",
#                       "Patient6", "Patient7", "Patient8", "PatientA", "PatientB",
#                       "PatientC")
# 
# all_raw_data$grade = grade_anno[all_raw_data$patient]
# all_raw_data$stage = stage_anno[all_raw_data$patient]

print(Sys.time())
library(gplots)
break_p = seq(0.8,1.2, 0.05)
png("test.png", width = 20, height = 30, units = "in", res = 300)
heatmap.2(as.matrix(patinet_infercnv_data), # data frame a matrix
          density.info = "none", # Remove density legend lines
          trace = "none", # Remove the blue trace lines from heatmap
          Rowv = FALSE, # Do not reorder the rows
          dendrogram = "none", # Only plot column dendrogram
          col = colorRampPalette(c("blue","white","red"))(length(break_p)-1), # Make colors viridis
          breaks = break_p)
dev.off()
print(Sys.time())
