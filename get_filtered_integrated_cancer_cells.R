source("call_libraries.R")
#number_of_cores = 8
#plan("multicore", workers = number_of_cores)
tissue = "Tumor"
cd45_status = "CD45-"
integrated_cancer_cells = readRDS("integrated_cancer_cells.rds")
InferCnv_output_dir = paste0("./infercnv_output/",cd45_status,"_",tissue,"_","integrated", "/")
InferCnv_output_exist_marker = paste0(InferCnv_output_dir, "infercnv.observations.txt")
gene_position = read.csv("human_gene_pos.csv", header = T)
gene_position_chr3 = gene_position[gene_position$seqname=="chr3" & gene_position$arm == "p",]
patient_infercnv_data = read.table(InferCnv_output_exist_marker, sep="", header = T)
patient_infercnv_data_chr3 = patient_infercnv_data[rownames(patient_infercnv_data) %in% gene_position_chr3$gene_name,]
print(length(rownames(patient_infercnv_data_chr3)))
transformed_patient_infercnv_data_chr3 = na.omit(data.frame(t(patient_infercnv_data_chr3)))
infercnv_signal_threshold = 0.985
transformed_patient_infercnv_data_chr3$keep = apply(transformed_patient_infercnv_data_chr3, 1, function(x){
  less_than_threshold_percent = length(x[x<=infercnv_signal_threshold])/length(x)
  if (less_than_threshold_percent >= 0.2){
    return("yes")
  }else{
    return("no")
  }
})
transformed_filtered_patient_infercnv_data = transformed_patient_infercnv_data_chr3[transformed_patient_infercnv_data_chr3$keep == "yes",]
print(nrow(transformed_filtered_patient_infercnv_data))
print(nrow(transformed_filtered_patient_infercnv_data)/nrow(transformed_patient_infercnv_data_chr3))
selected_cells = rownames(transformed_filtered_patient_infercnv_data)
selected_cells = sapply(selected_cells, function(x){
  stringr::str_replace_all(x, "\\.", "-")
})

integrated_cancer_cells$is_tumor_cell = sapply(names(integrated_cancer_cells$orig.ident), function(x){
  if(x %in% selected_cells){
    return("Yes")
  }else{
    return("No")
  }
})

filtered_integrated_cancer_cells = subset(integrated_cancer_cells, subset = is_tumor_cell == "Yes")
# filtered_integrated_cancer_cell = readRDS("filtered_integrated_cancer_cells_0.9_0.7.rds")
# DefaultAssay(filtered_integrated_cancer_cells) = "RNA"
# filtered_integrated_cancer_cells = DietSeurat(filtered_integrated_cancer_cells, assays = "RNA")
# filtered_integrated_cancer_cells = Integrate_seurat_obj_by_key(filtered_integrated_cancer_cells, 
#                                     key = "patient", method_SCT = T, method_rpca = T)
# filtered_integrated_cancer_cells = run_pca_and_umap(filtered_integrated_cancer_cells)
# filtered_integrated_cancer_cells = find_neighbours_and_clusters(filtered_integrated_cancer_cells, res = seq(0.05,0.3, 0.05))
saveRDS(filtered_integrated_cancer_cells, "filtered_integrated_cancer_cells_0.9_0.7.rds")