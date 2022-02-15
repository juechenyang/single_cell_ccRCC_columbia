source("call_libraries.R")
number_of_cores = 8
plan("multicore", workers = number_of_cores)
tissue = "Tumor"
cd45_status = "CD45-"
integrated_cancer_cells = readRDS("integrated_cancer_cells.rds")
InferCnv_output_dir = paste0("./infercnv_output/",cd45_status,"_",tissue,"_","integrated", "/")
InferCnv_output_exist_marker = paste0(InferCnv_output_dir, "infercnv.observations.txt")
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

integrated_cancer_cells$is_tumor_cell = sapply(names(patient_seurat$orig.ident), function(x){
  if(x %in% selected_cells){
    return("Yes")
  }else{
    return("No")
  }
})

filtered_integrated_cancer_cells = subset(integrated_cancer_cells, subset = is_tumor_cell == "Yes")
saveRDS(filtered_integrated_cancer_cells, "filtered_integrated_cancer_cells.rds")