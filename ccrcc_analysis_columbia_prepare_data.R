set.seed(1234)
library(reticulate)
library(factoextra)
library(flowCore)
library(SingleR)
library(Seurat)
library(cowplot)
library(dplyr)
library(iterClust)
library(cluster)
library(umap)
library(reshape)
library(pheatmap)
library(viper)
library("org.Hs.eg.db")
library(clustree)
library(leiden)
library(MAST)
library(Hmisc)
library(ggplot2)
library(scales)
library(ggcyto)
library(infercnv)
library(ggrepel)
library(plyr)
library(celldex)
plot_dir = "./plots/"
results_stratified_rds_dir = "./results_stratified_rds/"
raw_stratified_rds_dir = "./raw_stratified_rds/"

#*********************************************************************************************
# load all data
###
###LOAD AND ANNOTATE ALL DATA
###
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/180426_DRAEK_NIVI_2_HUMAN_10X/CN004/outs/filtered_gene_bc_matrices/GRCh38")
pbmca1 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmca1$patient="PatientA"
pbmca1$tissue="Normal"
pbmca1$cd45="CD45+"
pbmca1$batch="batch1"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/180426_DRAEK_NIVI_2_HUMAN_10X/CN005/outs/filtered_gene_bc_matrices/GRCh38")
pbmca2 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmca2$patient="PatientA"
pbmca2$tissue="Tumor"
pbmca2$cd45="CD45+"
pbmca2$batch="batch1"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/180730_CHARLES_NIVI_3_HUMAN_10X/CN009/outs/filtered_gene_bc_matrices/GRCh38")
pbmcb1 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmcb1$patient="PatientB"
pbmcb1$tissue="Normal"
pbmcb1$cd45="CD45+"
pbmcb1$batch="batch1"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/180730_CHARLES_NIVI_3_HUMAN_10X/CN010/outs/filtered_gene_bc_matrices/GRCh38")
pbmcb2 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmcb2$patient="PatientB"
pbmcb2$tissue="Tumor"
pbmcb2$cd45="CD45+"
pbmcb2$batch="batch1"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/180730_CHARLES_NIVI_3_HUMAN_10X/CN011/outs/filtered_gene_bc_matrices/GRCh38")
pbmcb3 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmcb3$patient="PatientB"
pbmcb3$tissue="Tumor"
pbmcb3$cd45="CD45+"
pbmcb3$batch="batch1"

pbmc_data <- Read10X(data.dir = "singlecell_rawdata/180807_CHARLES_NIVI_3_HUMAN_10X/CN012/outs/filtered_gene_bc_matrices/GRCh38")
pbmcc1 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmcc1$patient="PatientC"
pbmcc1$tissue="Normal"
pbmcc1$cd45="CD45+"
pbmcc1$batch="batch1"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/180807_CHARLES_NIVI_3_HUMAN_10X/CN013/outs/filtered_gene_bc_matrices/GRCh38")
pbmcc2 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmcc2$patient="PatientC"
pbmcc2$tissue="Tumor"
pbmcc2$cd45="CD45+"
pbmcc2$batch="batch1"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/180807_CHARLES_NIVI_3_HUMAN_10X/CN014/outs/filtered_gene_bc_matrices/GRCh38")
pbmcc3 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmcc3$patient="PatientC"
pbmcc3$tissue="Tumor"
pbmcc3$cd45="CD45+"
pbmcc3$batch="batch1"

pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190521_CHARLES_NIVI_9_HUMAN_10X/CN1/outs/filtered_feature_bc_matrix")
pbmc1 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc1$patient="Patient1"
pbmc1$tissue="Tumor"
pbmc1$cd45="CD45+"
pbmc1$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190521_CHARLES_NIVI_9_HUMAN_10X/CN2/outs/filtered_feature_bc_matrix")
pbmc2 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc2$patient="Patient1"
pbmc2$tissue="Tumor"
pbmc2$cd45="CD45-"
pbmc2$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190521_CHARLES_NIVI_9_HUMAN_10X/CN3/outs/filtered_feature_bc_matrix")
pbmc3 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc3$patient="Patient1"
pbmc3$tissue="Normal"
pbmc3$cd45="CD45+"
pbmc3$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190521_CHARLES_NIVI_9_HUMAN_10X/CN4/outs/filtered_feature_bc_matrix")
pbmc4 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc4$patient="Patient1"
pbmc4$tissue="Normal"
pbmc4$cd45="CD45-"
pbmc4$batch="batch2"

pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190521_CHARLES_NIVI_9_HUMAN_10X/CN5/outs/filtered_feature_bc_matrix")
pbmc5 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc5$patient="Patient2"
pbmc5$tissue="Tumor"
pbmc5$cd45="CD45+"
pbmc5$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190521_CHARLES_NIVI_9_HUMAN_10X/CN6/outs/filtered_feature_bc_matrix")
pbmc6 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc6$patient="Patient2"
pbmc6$tissue="Tumor"
pbmc6$cd45="CD45-"
pbmc6$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190521_CHARLES_NIVI_9_HUMAN_10X/CN7/outs/filtered_feature_bc_matrix")
pbmc7 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc7$patient="Patient2"
pbmc7$tissue="Normal"
pbmc7$cd45="CD45+"
pbmc7$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190521_CHARLES_NIVI_9_HUMAN_10X/CN8/outs/filtered_feature_bc_matrix")
pbmc8 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc8$patient="Patient2"
pbmc8$tissue="Normal"
pbmc8$cd45="CD45-"
pbmc8$batch="batch2"

pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190524_CHARLES_NIVI_10_HUMAN_10X/CN9/outs/filtered_feature_bc_matrix")
pbmc9 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc9$patient="Patient3"
pbmc9$tissue="Tumor"
pbmc9$cd45="CD45+"
pbmc9$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190524_CHARLES_NIVI_10_HUMAN_10X/CN10/outs/filtered_feature_bc_matrix")
pbmc10 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc10$patient="Patient3"
pbmc10$tissue="Tumor"
pbmc10$cd45="CD45-"
pbmc10$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190524_CHARLES_NIVI_10_HUMAN_10X/CN11/outs/filtered_feature_bc_matrix")
pbmc11 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc11$patient="Patient3"
pbmc11$tissue="Normal"
pbmc11$cd45="CD45+"
pbmc11$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190524_CHARLES_NIVI_10_HUMAN_10X/CN12/outs/filtered_feature_bc_matrix")
pbmc12 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc12$patient="Patient3"
pbmc12$tissue="Normal"
pbmc12$cd45="CD45-"
pbmc12$batch="batch2"

pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190603_CHARLES_NIVI_9_HUMAN_10X/CN13/outs/filtered_feature_bc_matrix")
pbmc13 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc13$patient="Patient4"
pbmc13$tissue="Tumor"
pbmc13$cd45="CD45+"
pbmc13$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190603_CHARLES_NIVI_9_HUMAN_10X/CN14/outs/filtered_feature_bc_matrix")
pbmc14 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc14$patient="Patient4"
pbmc14$tissue="Tumor"
pbmc14$cd45="CD45-"
pbmc14$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190603_CHARLES_NIVI_9_HUMAN_10X/CN15/outs/filtered_feature_bc_matrix")
pbmc15 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc15$patient="Patient4"
pbmc15$tissue="Normal"
pbmc15$cd45="CD45+"
pbmc15$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190603_CHARLES_NIVI_9_HUMAN_10X/CN16/outs/filtered_feature_bc_matrix")
pbmc16 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc16$patient="Patient4"
pbmc16$tissue="Normal"
pbmc16$cd45="CD45-"
pbmc16$batch="batch2"

pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190524_CHARLES_NIVI_10_HUMAN_10X/CN21b/outs/filtered_feature_bc_matrix")
pbmc_data=pbmc_data[,sample(colnames(pbmc_data),10000)]
pbmc21 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc21$patient="Patient5"
pbmc21$tissue="Tumor"
pbmc21$cd45="CD45+"
pbmc21$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190524_CHARLES_NIVI_10_HUMAN_10X/CN22b/outs/filtered_feature_bc_matrix")
pbmc_data=pbmc_data[,sample(colnames(pbmc_data),10000)]
pbmc22 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc22$patient="Patient5"
pbmc22$tissue="Tumor"
pbmc22$cd45="CD45-"
pbmc22$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190524_CHARLES_NIVI_10_HUMAN_10X/CN23/outs/filtered_feature_bc_matrix")
pbmc_data=pbmc_data[,sample(colnames(pbmc_data),10000)]
pbmc23 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc23$patient="Patient5"
pbmc23$tissue="Normal"
pbmc23$cd45="CD45+"
pbmc23$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190524_CHARLES_NIVI_10_HUMAN_10X/CN24/outs/filtered_feature_bc_matrix")
pbmc_data=pbmc_data[,sample(colnames(pbmc_data),10000)]
pbmc24 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc24$patient="Patient5"
pbmc24$tissue="Normal"
pbmc24$cd45="CD45-"
pbmc24$batch="batch2"

pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190606_CHARLES_NIVI_8_HUMAN_10X/CN29/outs/filtered_feature_bc_matrix")
pbmc29 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc29$patient="Patient6"
pbmc29$tissue="Tumor"
pbmc29$cd45="CD45+"
pbmc29$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190606_CHARLES_NIVI_8_HUMAN_10X/CN30/outs/filtered_feature_bc_matrix")
pbmc30 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc30$patient="Patient6"
pbmc30$tissue="Tumor"
pbmc30$cd45="CD45-"
pbmc30$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190606_CHARLES_NIVI_8_HUMAN_10X/CN31/outs/filtered_feature_bc_matrix")
pbmc31 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc31$patient="Patient6"
pbmc31$tissue="Normal"
pbmc31$cd45="CD45+"
pbmc31$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190606_CHARLES_NIVI_8_HUMAN_10X/CN32/outs/filtered_feature_bc_matrix")
pbmc32 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc32$patient="Patient6"
pbmc32$tissue="Normal"
pbmc32$cd45="CD45-"
pbmc32$batch="batch2"

pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190606_CHARLES_NIVI_8_HUMAN_10X/CN33/outs/filtered_feature_bc_matrix")
pbmc33 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc33$patient="Patient7"
pbmc33$tissue="Tumor"
pbmc33$cd45="CD45+"
pbmc33$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190606_CHARLES_NIVI_8_HUMAN_10X/CN34/outs/filtered_feature_bc_matrix")
pbmc34 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc34$patient="Patient7"
pbmc34$tissue="Tumor"
pbmc34$cd45="CD45-"
pbmc34$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190606_CHARLES_NIVI_8_HUMAN_10X/CN35/outs/filtered_feature_bc_matrix")
pbmc35 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc35$patient="Patient7"
pbmc35$tissue="Normal"
pbmc35$cd45="CD45+"
pbmc35$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190606_CHARLES_NIVI_8_HUMAN_10X/CN36/outs/filtered_feature_bc_matrix")
pbmc36 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc36$patient="Patient7"
pbmc36$tissue="Normal"
pbmc36$cd45="CD45-"
pbmc36$batch="batch2"
rm(pbmc_data)

pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190825_VICTOR_VICTOR_1_HUMAN_10X/CN37/outs/filtered_feature_bc_matrix")
pbmc37 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc37$patient="Patient8"
pbmc37$tissue="Tumor"
pbmc37$cd45="CD45+"
pbmc37$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190825_VICTOR_VICTOR_1_HUMAN_10X/CN38/outs/filtered_feature_bc_matrix")
pbmc38 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc38$patient="Patient8"
pbmc38$tissue="Tumor"
pbmc38$cd45="CD45-"
pbmc38$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190825_VICTOR_VICTOR_1_HUMAN_10X/CN39/outs/filtered_feature_bc_matrix")
pbmc39 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc39$patient="Patient8"
pbmc39$tissue="Normal"
pbmc39$cd45="CD45+"
pbmc3$batch="batch2"
pbmc_data <- Read10X(data.dir = "singlecell_rawdata/190825_VICTOR_VICTOR_1_HUMAN_10X/CN40/outs/filtered_feature_bc_matrix")
pbmc40 <- CreateSeuratObject(counts = pbmc_data,min.features = 500,min.cells=50)
pbmc40$patient="Patient8"
pbmc40$tissue="Normal"
pbmc40$cd45="CD45-"
pbmc40$batch="batch2"

pbmc.big=merge(pbmc1,y=c(pbmc2,pbmc3,pbmc4,pbmc5,pbmc6,pbmc7,pbmc8,pbmc9,pbmc10),project = "RCC_SC")
pbmc.big=merge(pbmc.big,y=c(pbmc11,pbmc12,pbmc13,pbmc14,pbmc15,pbmc16,pbmc21,pbmc22,pbmc23,pbmc24),project="RCC_SC")
pbmc.big=merge(pbmc.big,y=c(pbmc29,pbmc30,pbmc31,pbmc32,pbmc33,pbmc34,pbmc35,pbmc36,pbmc37,pbmc38,pbmc39,pbmc40),project="RCC_SC")
pbmc.big<- merge(pbmc.big, y = c(pbmca1,pbmca2,pbmcb1,pbmcb2,pbmcb3,pbmcc1,pbmcc2,pbmcc3), project = "RCC_SC")
rm(pbmca1,pbmca2,pbmcb1,pbmcb2,pbmcb3,pbmcc1,pbmcc2,pbmcc3)
rm(pbmc1,pbmc2,pbmc3,pbmc4,pbmc5,pbmc6,pbmc7,pbmc8,pbmc9,pbmc10,pbmc11,pbmc12,pbmc13,pbmc14,pbmc15,pbmc16,pbmc21,pbmc22,pbmc23,pbmc24,pbmc29,pbmc30,pbmc31,pbmc32,pbmc33,pbmc34,pbmc35,pbmc36,pbmc37,pbmc38,pbmc39,pbmc40)
big_list=SplitObject(pbmc.big, split.by = "cd45")
rm(pbmc.big)
rm(pbmc_data)

#*********************************************************************************************
#*qc for CD45 pos and neg
pbmc.big=big_list[[1]]
pbmc.big <- PercentageFeatureSet(pbmc.big, pattern = "^MT-", col.name = "percent.mt")
#VlnPlot(pbmc.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0,group.by="tissue")
pbmc.big <- subset(pbmc.big, subset = percent.mt < 10 & nCount_RNA > 1500 & nCount_RNA < 15000)
VlnPlot(pbmc.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0,group.by="tissue") 

pbmc.big=big_list[[2]]
pbmc.big <- PercentageFeatureSet(pbmc.big, pattern = "^MT-", col.name = "percent.mt")
pbmc.big <- subset(pbmc.big, subset = percent.mt < 10 & nCount_RNA > 1500 & nCount_RNA < 15000)
VlnPlot(pbmc.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0,group.by="tissue")
# seurat_list_cd45pos=SplitObject(big_list[[1]],split.by="patient")
# seurat_list_cd45neg=SplitObject(big_list[[2]],split.by="patient")

seurat_list_cd45pos=SplitObject(big_list[[1]],split.by="tissue")
seurat_list_cd45neg=SplitObject(big_list[[2]],split.by="tissue")

# split seurat obj by patient and by tissue
# seurat_list_cd45neg_split_by_tissue = list()
# for(i in 1:length(seurat_list_cd45neg)){
#   inner_list = SplitObject(seurat_list_cd45neg[[i]], split.by = "tissue")
#   seurat_list_cd45neg_split_by_tissue = append(seurat_list_cd45neg_split_by_tissue, inner_list[[1]])
#   seurat_list_cd45neg_split_by_tissue = append(seurat_list_cd45neg_split_by_tissue, inner_list[[2]])
# }
# 
# seurat_list_cd45pos_split_by_tissue = list()
# for(i in 1:length(seurat_list_cd45pos)){
#   inner_list = SplitObject(seurat_list_cd45pos[[i]], split.by = "tissue")
#   seurat_list_cd45pos_split_by_tissue = append(seurat_list_cd45pos_split_by_tissue, inner_list[[1]])
#   seurat_list_cd45pos_split_by_tissue = append(seurat_list_cd45pos_split_by_tissue, inner_list[[2]])
# }

# for(i in 1:length(seurat_list_cd45neg)){
#   inner_list = SplitObject(seurat_list_cd45neg[[i]], split.by = "tissue")
#   seurat_list_cd45neg_split_by_tissue = append(seurat_list_cd45neg_split_by_tissue, inner_list[[1]])
#   seurat_list_cd45neg_split_by_tissue = append(seurat_list_cd45neg_split_by_tissue, inner_list[[2]])
# }
