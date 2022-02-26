source("call_libraries.R")
number_of_cores = 8
plan("multicore", workers = number_of_cores)
tissue = "Tumor"
cd45_status = "CD45-"
#load data
integrated = readRDS("integrated_cancer_cells.rds")
reference_pt = readRDS("reference_pt_cells.rds")
#define grouping variable
integrated$individual_anno = as.character(integrated$integrated_snn_res.0.2)
patient = merge(integrated, reference_pt)
InferCnv_output_dir = paste0("./infercnv_output/",cd45_status,"_",tissue,"_","integrated", "/")
InferCnv_output_exist_marker = paste0(InferCnv_output_dir, "infercnv.observations.txt")
if(!file.exists(InferCnv_output_exist_marker)){
  InfercnvInputs = prepare_infercnv_data(patient, annotation = "individual_anno")
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=InfercnvInputs[[1]],
                                      annotations_file=InfercnvInputs[[2]],
                                      delim="\t",
                                      gene_order_file=InfercnvInputs[[3]],
                                      ref_group_names = c("Proximal Tubule"),
                                      max_cells_per_group = 2000
  )
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir=InferCnv_output_dir,
                               cluster_by_groups=TRUE, denoise=TRUE,
                               HMM = T, num_threads = number_of_cores)
}
