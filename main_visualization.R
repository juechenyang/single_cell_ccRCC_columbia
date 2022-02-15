library(Seurat)
library(dplyr)
library(infercnv)
library(ggpubr)

######################################all cells######################################
######################################without integration######################################
#read all cell data without integration
neg_tumor_azimuth_standard = readRDS("neg_tumor_azimuth_standard.rds")
neg_tumor_azimuth_SCT = readRDS("neg_tumor_azimuth_SCT.rds")

all_types = as.character(unique(neg_tumor_azimuth_standard$individual_anno))
palette_color = c("chocolate1", "cyan", "gold", "aquamarine", "deepskyblue", 
                  "cyan4", "darkblue", "darkolivegreen1", "darkorchid1",
                  "firebrick1", "firebrick4", "purple", "darksalmon","darkslategray", "darkmagenta")
all_colors = colorRampPalette(palette_color)(length(all_types))
names(all_colors) = all_types

#plot standard pipeline umaps for azimuth annotation
data_cell_types = as.character(unique(neg_tumor_azimuth_standard$predicted.annotation.l1))
p1 = DimPlot(neg_tumor_azimuth_standard, reduction = "umap", group.by = "predicted.annotation.l1", 
             pt.size = 0.2, label = F)+ggtitle("Standard")+ 
  scale_color_manual(breaks = data_cell_types, 
                     values=all_colors[data_cell_types])
#plot SCT pipeline umaps for azimuth annotation
p2 = DimPlot(neg_tumor_azimuth_SCT, reduction = "umap", group.by = "predicted.annotation.l1", 
             pt.size = 0.2, label = F)+ggtitle("SCT")+ 
  scale_color_manual(breaks = data_cell_types, 
                     values=all_colors[data_cell_types])
png("neg_tumor_AzimuthL1Umap.png",24,9, units = "in", res = 600)
ggarrange(
  p1,p2,
  ncol = 2
)
dev.off()

#plot standard pipeline umaps for patients
p1 = DimPlot(neg_tumor_azimuth_standard, reduction = "umap", group.by = "patient", 
             pt.size = 0.2, label = F)+ggtitle("Standard")
#plot SCT pipeline umaps for patients
p2 = DimPlot(neg_tumor_azimuth_SCT, reduction = "umap", group.by = "patient", 
             pt.size = 0.2, label = F)+ggtitle("SCT")
png("neg_tumor_PatientUmap.png",20,9, units = "in", res = 600)
ggarrange(
  p1,p2,
  ncol = 2
)
dev.off()

#plot standard pipeline umaps for patients
p1 = DimPlot(neg_tumor_azimuth_standard, reduction = "umap", group.by = "batch", 
             pt.size = 0.2, label = F)+ggtitle("Standard")
#plot SCT pipeline umaps for patients
p2 = DimPlot(neg_tumor_azimuth_SCT, reduction = "umap", group.by = "batch", 
             pt.size = 0.2, label = F)+ggtitle("SCT")
png("neg_tumor_BatchUmap.png",20,9, units = "in", res = 600)
ggarrange(
  p1,p2,
  ncol = 2
)
dev.off()

######################################with integration######################################
#read all cell data with integration
neg_tumor_azimuth_integrated_by_patient_standard = readRDS("neg_tumor_azimuth_integrated_by_patient_standard.rds")
neg_tumor_azimuth_integrated_by_patient_SCT = readRDS("neg_tumor_azimuth_integrated_by_patient_SCT.rds")

data_cell_types = as.character(unique(neg_tumor_azimuth_integrated_by_patient_standard$predicted.annotation.l1))
p1 = DimPlot(neg_tumor_azimuth_integrated_by_patient_standard, reduction = "umap", group.by = "predicted.annotation.l1", 
             pt.size = 0.2, label = F)+ggtitle("Standard")+ 
  scale_color_manual(breaks = data_cell_types, 
                     values=all_colors[data_cell_types])
#plot SCT pipeline umaps for azimuth annotation
p2 = DimPlot(neg_tumor_azimuth_integrated_by_patient_SCT, reduction = "umap", group.by = "predicted.annotation.l1", 
             pt.size = 0.2, label = F)+ggtitle("SCT")+ 
  scale_color_manual(breaks = data_cell_types, 
                     values=all_colors[data_cell_types])
png("neg_tumor_IntegratedByPatient_AzimuthL1Umap.png",24,9, units = "in", res = 600)
ggarrange(
  p1,p2,
  ncol = 2
)
dev.off()

p1 = DimPlot(neg_tumor_azimuth_integrated_by_patient_standard, reduction = "umap", group.by = "patient", 
             pt.size = 0.2, label = F)+ggtitle("Standard")
#plot SCT pipeline umaps for azimuth annotation
p2 = DimPlot(neg_tumor_azimuth_integrated_by_patient_SCT, reduction = "umap", group.by = "patient", 
             pt.size = 0.2, label = F)+ggtitle("SCT")
png("neg_tumor_IntegratedByPatient_PatientUmap.png",24,9, units = "in", res = 600)
ggarrange(
  p1,p2,
  ncol = 2
)
dev.off()

######################################likely cancer cells######################################
######################################without integration######################################
neg_tumor_epi_azimuth_standard = readRDS("neg_tumor_epi_azimuth_standard.rds")
neg_tumor_epi_azimuth_SCT = readRDS("neg_tumor_epi_azimuth_SCT.rds")

#visualization
data_cell_types = as.character(unique(neg_tumor_epi_azimuth_standard$predicted.annotation.l1))
p1 = DimPlot(neg_tumor_epi_azimuth_standard, reduction = "umap", group.by = "predicted.annotation.l1", 
             pt.size = 0.2, label = F)+ggtitle("Standard")+ 
  scale_color_manual(breaks = data_cell_types, 
                     values=all_colors[data_cell_types])
#plot SCT pipeline umaps for patients
p2 = DimPlot(neg_tumor_epi_azimuth_SCT, reduction = "umap", group.by = "predicted.annotation.l1", 
             pt.size = 0.2, label = F)+ggtitle("SCT")+ 
  scale_color_manual(breaks = data_cell_types, 
                     values=all_colors[data_cell_types])
png("neg_tumor_epi_AzimuthL1Umap.png",20,9, units = "in", res = 600)
ggarrange(
  p1,p2,
  ncol = 2
)
dev.off()

#visualization
p1 = DimPlot(neg_tumor_epi_azimuth_standard, reduction = "umap", group.by = "patient", 
             pt.size = 0.2, label = F)+ggtitle("Standard")
#plot SCT pipeline umaps for patients
p2 = DimPlot(neg_tumor_epi_azimuth_SCT, reduction = "umap", group.by = "patient", 
             pt.size = 0.2, label = F)+ggtitle("SCT")
png("neg_tumor_epi_PatientUmap.png",20,9, units = "in", res = 600)
ggarrange(
  p1,p2,
  ncol = 2
)
dev.off()

######################################with integration######################################
neg_tumor_epi_azimuth_integrated_by_patient_standard = readRDS("neg_tumor_epi_azimuth_integrated_by_patient_standard.rds")
neg_tumor_epi_azimuth_integrated_by_patient_SCT = readRDS("neg_tumor_epi_azimuth_integrated_by_patient_SCT.rds")


data_cell_types = as.character(unique(neg_tumor_epi_azimuth_integrated_by_patient_standard$predicted.annotation.l1))
p1 = DimPlot(neg_tumor_epi_azimuth_integrated_by_patient_standard, reduction = "umap", group.by = "predicted.annotation.l1", 
             pt.size = 0.2, label = F)+ggtitle("Standard")+ 
  scale_color_manual(breaks = data_cell_types, 
                     values=all_colors[data_cell_types])
#plot SCT pipeline umaps for azimuth annotation
p2 = DimPlot(neg_tumor_epi_azimuth_integrated_by_patient_SCT, reduction = "umap", group.by = "predicted.annotation.l1", 
             pt.size = 0.2, label = F)+ggtitle("SCT")+ 
  scale_color_manual(breaks = data_cell_types, 
                     values=all_colors[data_cell_types])
png("neg_tumor_epi_IntegratedByPatient_AzimuthL1Umap.png",20,9, units = "in", res = 600)
ggarrange(
  p1,p2,
  ncol = 2
)
dev.off()

p1 = DimPlot(neg_tumor_epi_azimuth_integrated_by_patient_standard, reduction = "umap", group.by = "patient", 
             pt.size = 0.2, label = F)+ggtitle("Standard")
#plot SCT pipeline umaps for azimuth annotation
p2 = DimPlot(neg_tumor_epi_azimuth_integrated_by_patient_SCT, reduction = "umap", group.by = "patient", 
             pt.size = 0.2, label = F)+ggtitle("SCT")
png("neg_tumor_epi_IntegratedByPatient_PatientUmap.png",20,9, units = "in", res = 600)
ggarrange(
  p1,p2,
  ncol = 2
)
dev.off()

######################################patient level visualization######################################
######################################all cells######################################

neg_tumor_patient_list_standard = readRDS("neg_tumor_patient_list_standard.rds")
neg_tumor_patient_list_SCT = readRDS("neg_tumor_patient_list_SCT.rds")
neg_tumor_epi_patient_list_SCT = readRDS("neg_tumor_epi_patient_list_SCT.rds")
neg_tumor_epi_patient_list_standard = readRDS("neg_tumor_epi_patient_list_standard.rds")

####### plot patient level data
plotlist = list()
for(i in 1:length(neg_tumor_patient_list_standard)){
  single_obj = neg_tumor_patient_list_standard[[i]]
  single_obj_unique_annos = unique(single_obj$individual_anno)
  p = DimPlot(neg_tumor_patient_list_standard[[i]], reduction = "umap", group.by = "individual_anno", 
              pt.size = 0.2, label = F)+ggtitle(single_obj$patient) + 
    scale_color_manual(breaks = single_obj_unique_annos, 
                       values=all_colors[single_obj_unique_annos])
  plotlist[[i]] = p
}
png("neg_tumor_patient_level_standard.png",25,16, units = "in", res = 600)
ggarrange(
  plotlist = plotlist,
  nrow = 3,
  ncol = 3,
  legend = "right",
  common.legend = F
)
dev.off()



plotlist = list()
for(i in 1:length(neg_tumor_patient_list_standard)){
  single_obj = neg_tumor_patient_list_standard[[i]]
  single_obj_unique_annos = unique(single_obj$individual_anno)
  p = DimPlot(neg_tumor_patient_list_standard[[i]], reduction = "umap", group.by = "individual_anno", 
              pt.size = 0.2, label = F)+ggtitle(single_obj$patient) + 
    scale_color_manual(breaks = single_obj_unique_annos, 
                       values=all_colors[single_obj_unique_annos])
  plotlist[[i]] = p
}
png("neg_tumor_patient_level_standard_.png",25,16, units = "in", res = 600)
ggarrange(
  plotlist = plotlist,
  nrow = 3,
  ncol = 3,
  legend = "right",
  common.legend = F
)
dev.off()

plotlist = list()
for(i in 1:length(neg_tumor_patient_list_SCT)){
  single_obj = neg_tumor_patient_list_SCT[[i]]
  single_obj_unique_annos = unique(single_obj$individual_anno)
  p = DimPlot(neg_tumor_patient_list_SCT[[i]], reduction = "umap", group.by = "individual_anno", 
              pt.size = 0.2, label = F)+ggtitle(single_obj$patient)+ 
      scale_color_manual(breaks = single_obj_unique_annos, 
                       values=all_colors[single_obj_unique_annos])
  plotlist[[i]] = p
}
png("neg_tumor_patient_level_SCT.png",25,16, units = "in", res = 600)
ggarrange(
  plotlist = plotlist,
  nrow = 3,
  ncol = 3
)
dev.off()

######################################likely cancer cells######################################
plotlist = list()
for(i in 1:length(neg_tumor_epi_patient_list_standard)){
  single_obj = neg_tumor_epi_patient_list_standard[[i]]
  single_obj_unique_annos = unique(single_obj$individual_anno)
  p = DimPlot(neg_tumor_epi_patient_list_standard[[i]], reduction = "umap", group.by = "individual_anno", 
              pt.size = 0.2, label = F)+ggtitle(single_obj$patient) + 
    scale_color_manual(breaks = single_obj_unique_annos, 
                       values=all_colors[single_obj_unique_annos])
  plotlist[[i]] = p
}
png("neg_tumor_epi_patient_level_standard.png",25,16, units = "in", res = 600)
ggarrange(
  plotlist = plotlist,
  nrow = 3,
  ncol = 3
)
dev.off()

plotlist = list()
for(i in 1:length(neg_tumor_epi_patient_list_SCT)){
  single_obj = neg_tumor_epi_patient_list_SCT[[i]]
  single_obj_unique_annos = unique(single_obj$individual_anno)
  p = DimPlot(neg_tumor_epi_patient_list_SCT[[i]], reduction = "umap", group.by = "individual_anno", 
              pt.size = 0.2, label = F)+ggtitle(single_obj$patient)+ 
    scale_color_manual(breaks = single_obj_unique_annos, 
                       values=all_colors[single_obj_unique_annos])
  plotlist[[i]] = p
}
png("neg_tumor_epi_patient_level_SCT.png",25,16, units = "in", res = 600)
ggarrange(
  plotlist = plotlist,
  nrow = 3,
  ncol = 3
)
dev.off()







plotlist = list()
for(i in 1:length(neg_tumor_epi_patient_list_standard)){
  single_obj = neg_tumor_epi_patient_list_standard[[i]]
  single_obj_unique_annos = unique(single_obj$individual_anno)
  p = FeaturePlot(neg_tumor_epi_patient_list_standard[[i]], reduction = "umap", features = "NDUFA4L2", 
              pt.size = 0.2, label = F)+ggtitle(single_obj$patient)
  plotlist[[i]] = p
}
png("neg_tumor_epi_patient_level_standard_NDUFA4L2.png",25,16, units = "in", res = 600)
ggarrange(
  plotlist = plotlist,
  nrow = 3,
  ncol = 3
)
dev.off()


plotlist = list()
for(i in 1:length(neg_tumor_patient_list_standard)){
  single_obj = neg_tumor_patient_list_standard[[i]]
  single_obj_unique_annos = unique(single_obj$individual_anno)
  p = FeaturePlot(neg_tumor_patient_list_standard[[i]], reduction = "umap", features = "NDUFA4L2", 
                  pt.size = 0.2, label = F)+ggtitle(single_obj$patient) 
  plotlist[[i]] = p
}
png("neg_tumor_patient_level_standard_NDUFA4L2.png",25,16, units = "in", res = 600)
ggarrange(
  plotlist = plotlist,
  nrow = 3,
  ncol = 3,
  legend = "right",
  common.legend = F
)
dev.off()
