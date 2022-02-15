source("call_libraries.R")
plan("multicore", workers = 16)
#read the raw data
all_raw_data=readRDS("./all_raw.rds")
#get tissue data
tissue = "Tumor"
cohort = SplitObject(all_raw_data, split.by = "tissue")[[tissue]]
#get cd45 status data
cd45_status = "CD45-"
cohort = SplitObject(cohort, split.by = "cd45")[[cd45_status]]
cohort = PercentageFeatureSet(cohort, pattern = "^MT-", col.name = "percent.mt")
#VlnPlot(cohort, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
cohort = subset(cohort, subset = percent.mt <= 30 & 
                  nCount_RNA >= 1000 & nCount_RNA <= 30000 &
                  nFeature_RNA >= 500 & nFeature_RNA <= 5000)
integrated = Integrate_seurat_obj_by_key(cohort, key = "patient",
                                         method_SCT = T, method_rpca = F)
integrated_azimuth = run_azimuth(integrated, reference_path = "./azimuth_reference/")
integrated$individual_anno = as.character(integrated_azimuth$predicted.annotation.l1)
saveRDS(integrated, "cd45neg_tumor_integrated_SCT_nonRPCA_azimuth.rds")
