source("call_libraries.R")
source("Initialization.R")
integrated_cancer_cells = readRDS("integrated_cancer_cells.rds")
DefaultAssay(integrated_cancer_cells) = "RNA"
integrated_cancer_cells = NormalizeData(integrated_cancer_cells)
#add metastatic status for each cell
integrated_cancer_cells$stage = 
  ifelse(integrated_cancer_cells$patient=="Patient5", "Metastatic", integrated_cancer_cells$stage)

#add module score to selected signature
signature_list = list(
                      "Complex_I"=Complex_I
                      ,"Complex_II"=Complex_II
                      ,"Complex_III"=Complex_III
                      ,"Complex_IV"=Complex_IV
                      ,"MT" = MT
                      ,"TCA"=TCA
                      ,"glycolysis"=glycolysis
                      ,"HIF_1A" = HIF1A
                      ,"HLA"=HLA
                      ,"MRP_positive"=MRP_positive
                      ,"NRF2"=NRF2
)

# signature_list = list("Iodothyronine_deiodinases_1_3"=Iodothyronine_deiodinases_1_3,
#                       "Glutathione_peroxidases"=Glutathione_peroxidases,
#                       "Selenoproteins"=Selenoproteins,
#                       "Selenophosphate_synthetase_2"=Selenophosphate_synthetase_2
# )
signature_list_names = names(signature_list)
integrated_cancer_cells = AddModuleScore(integrated_cancer_cells, 
                                        features = signature_list, 
                                        assay = "RNA", 
                                        name = signature_list_names)
integrated_cancer_cells@meta.data[,signature_list_names] = 
integrated_cancer_cells@meta.data[,paste0(signature_list_names, 1:length(signature_list_names))]
#remove redundant features
remove_column_index = match(paste0(signature_list_names, 1:length(signature_list_names)),
                            names(integrated_cancer_cells@meta.data))
integrated_cancer_cells@meta.data = 
integrated_cancer_cells@meta.data[,-remove_column_index]

#define group
integrated_cancer_cells@meta.data[,paste0(signature_list_names, "_group")] = 
  lapply(integrated_cancer_cells@meta.data[,signature_list_names], FUN = function(vx){
    vx = sapply(vx, FUN = function(x){
      if(x < -0.25){
        return("<-0.25")
      }else if(x >= -0.25 & x < -0.1){
        return("-0.25 to -0.1")
      }else if(x >= -0.1 & x < 0.1){
        return("-0.1 to 0.1")
      }else if(x >= 0.1 & x < 0.25){
        return("0.1 to 0.25")
      }else{
        return(">=0.25")
      }
    })
  })



#draw feature plots with scaling
louvain_plots = list()
feature_plots = list()
group_order <- c(">=0.25", "0.1 to 0.25", "-0.1 to 0.1", "-0.25 to -0.1", "<-0.25")
for(x in c("pT1b","pT3a","Metastatic")){
  integrated_cancer_cells_stage = subset(integrated_cancer_cells, subset = stage == x)
  
  target_res = 0.2
  target_res_para = paste0("integrated_snn_res.", target_res)
  all_types = unique(integrated_cancer_cells_stage[[target_res_para]])[[target_res_para]]
  all_types = sort(all_types)
  all_colors = colorRampPalette(palette_color)(length(all_types))
  names(all_colors) = all_types
  cluster_maps = DimPlot(integrated_cancer_cells_stage, reduction = "umap",
                         group.by = "integrated_snn_res.0.2")+
    scale_color_manual(
      breaks = all_types, 
      values=all_colors[all_types])+ggtitle(x)+
    theme(plot.title = element_text(size=32), legend.text=element_text(size=30))+
    guides(colour = guide_legend(override.aes = list(size=12)))
  louvain_plots = append(louvain_plots, list(cluster_maps))
  
  all_types = unique(integrated_cancer_cells_stage@meta.data[,paste0(signature_list_names[1], "_group")])
  all_types = sort(all_types)
  all_colors = c("grey", "cyan", "blue", "gold", "gold3")
  names(all_colors) = all_types
  
  #sort legend
  all_types_factor = factor(all_types, levels = group_order)
  all_types = as.character(sort(all_types_factor))
  all_colors = all_colors[all_types]
  
  
  feature_group_plots = lapply(signature_list_names, FUN = function(signature){
    plot_title = signature
    if(x!="pT1b"){
      plot_title = ""
    }
    
    integrated_cancer_cells_stage_signature = subset(integrated_cancer_cells, subset = stage == x)
    p = DimPlot(integrated_cancer_cells_stage, reduction = "umap",
            group.by = paste0(signature, "_group"))+
      scale_color_manual(
        breaks = all_types, 
        values=all_colors[all_types])+
      theme(plot.title = element_text(size=32), legend.text=element_text(size=30))+
      ggtitle(plot_title)+
      guides(colour = guide_legend(override.aes = list(size=12)))
    return(p)
  })
  
  feature_plots = append(feature_plots, feature_group_plots)
}

louvain_plots_rev = lapply(louvain_plots, FUN = function(x){
  x = x + guides(colour = guide_legend(override.aes = list(size=12)))+
    theme(plot.title = element_text(size=32),legend.text=element_text(size=30))
  return(x)
})

setEPS()

# naming the eps file
postscript("module_score_u_map_by_stage.eps", height = 1080, width = 2160)

# plotting the x and y position
# vectors
png("module_score_u_map_by_stage.png",width = 60, height = 20, units = "in", res = 300)
ggarrange(
  ggarrange(plotlist=louvain_plots_rev, nrow = 3, common.legend = T, legend = "left"),
  ggarrange(plotlist=feature_plots, nrow = 3, ncol = 10, common.legend = T, legend = "right")
  ,ncol = 2
  ,widths = c(1.2, 10)
)
dev.off()

# png("human_selenoproteins_split_by_stage.png",width = 20, height = 8, units = "in", res = 400)
# ggarrange(
#   plotlist = all_plots
#   ,ncol = 5
#   ,nrow = 2
# )
# dev.off()


#function to add z-score count of super gene to object
AddScaledNormalizeCount = function(object, gene_list){
  mat = object[["RNA"]]@data
  result_list = list()
  for(x in 1:length(gene_list)){
    selected_genes = gene_list[[x]]
    marker_tag = names(gene_list[x])
    is_in_mat = selected_genes %in% rownames(mat)
    if(!all(is_in_mat)){
      not_available_genes = selected_genes[!is_in_mat]
      print(paste0(paste(not_available_genes, collapse = ","), 
                   " is not available for analysis"))
    }
    selected_genes = selected_genes[is_in_mat]
    selected_mat = mat[selected_genes, ]
    composite_marker_exp = as.numeric(apply(selected_mat, MARGIN = 2, FUN = mean))
    composite_marker_exp = scale(composite_marker_exp)
    marker_tag = paste0(marker_tag,"_Normalize_Zscore")
    result_list[[marker_tag]] = composite_marker_exp
  }
  result_list = as.data.frame(result_list)
  rownames(result_list) = colnames(x=object)
  object[[colnames(x = result_list)]] <- result_list
  return(object)
}

integrated_cancer_cells = AddScaledNormalizeCount(integrated_cancer_cells, signature_list)


hist_plts = lapply(signature_list_names, FUN = function(sig){
  p = ggplot(integrated_cancer_cells@meta.data, aes_string(x=paste0(sig, "_Normalize_Zscore")))+
    geom_histogram(color="black", fill="lightblue")+
    ggtitle(sig)
  return(p)
})

png("hist.png",width = 18, height = 12, units = "in", res = 400)
ggarrange(
  plotlist = hist_plts,
  ncol = 4,
  nrow = 3
)
dev.off()

#define group
integrated_cancer_cells@meta.data[,paste0(signature_list_names, "_Normalize_Zscore", "_group")] = 
  lapply(integrated_cancer_cells@meta.data[,paste0(signature_list_names, "_Normalize_Zscore")], FUN = function(vx){
    vx = sapply(vx, FUN = function(x){
      if(x < -1.5){
        return("<-1.5")
      }else if(x >= -1.5 & x < -0.5){
        return("-1.5 to -0.5")
      }else if(x >= -0.5 & x < 0.5){
        return("-0.5 to 0.5")
      }else if(x >= 0.5 & x < 1.5){
        return("0.5 to 1.5")
      }else{
        return(">=1.5")
      }
    })
  })

#draw feature plots with scaling
louvain_plots = list()
feature_plots = list()
group_order <- c(">=1.5", "0.5 to 1.5", "-0.5 to 0.5", "-1.5 to -0.5", "<-1.5")
for(x in c("pT1b","pT3a","Metastatic")){
  integrated_cancer_cells_stage = subset(integrated_cancer_cells, subset = stage == x)
  
  target_res = 0.2
  target_res_para = paste0("integrated_snn_res.", target_res)
  all_types = unique(integrated_cancer_cells_stage[[target_res_para]])[[target_res_para]]
  all_types = sort(all_types)
  all_colors = colorRampPalette(palette_color)(length(all_types))
  names(all_colors) = all_types
  cluster_maps = DimPlot(integrated_cancer_cells_stage, reduction = "umap",
                         group.by = "integrated_snn_res.0.2")+
    scale_color_manual(
      breaks = all_types, 
      values=all_colors[all_types])+ggtitle(x)+
    theme(plot.title = element_text(size=32))
  louvain_plots = append(louvain_plots, list(cluster_maps))
  
  all_types = unique(integrated_cancer_cells_stage@meta.data[,paste0(signature_list_names[1], "_Normalize_Zscore", "_group")])
  all_types = sort(all_types)
  all_colors = c("grey", "cyan", "blue", "gold", "gold3")
  names(all_colors) = all_types
  
  #sort legend
  all_types_factor = factor(all_types, levels = group_order)
  all_types = as.character(sort(all_types_factor))
  all_colors = all_colors[all_types]
  
  
  feature_group_plots = lapply(signature_list_names, FUN = function(signature){
    plot_title = signature
    if(x!="pT1b"){
      plot_title = ""
    }
    
    integrated_cancer_cells_stage_signature = subset(integrated_cancer_cells, subset = stage == x)
    p = DimPlot(integrated_cancer_cells_stage, reduction = "umap",
                group.by = paste0(signature, "_Normalize_Zscore", "_group"))+
      scale_color_manual(
        breaks = all_types, 
        values=all_colors[all_types])+
      theme(plot.title = element_text(size=32),legend.text=element_text(size=30))+
      guides(colour = guide_legend(override.aes = list(size=12)))+
      ggtitle(plot_title)
    return(p)
  })
  
  feature_plots = append(feature_plots, feature_group_plots)p
}

png("feature_zscore_group_map_by_stage.png",width = 70, height = 20, units = "in", res = 400)
ggarrange(
  ggarrange(plotlist=louvain_plots_rev, nrow = 3, common.legend = T, legend = "left"),
  ggarrange(plotlist=feature_plots, nrow = 3, ncol = 11, common.legend = T, legend = "right")
  ,ncol = 2
  ,widths = c(1.2, 11)
)
dev.off()

# feature_plots_rev = lapply(feature_plots, FUN = function(x){
#   x = x + theme(plot.title = element_text(size=32), legend.text=element_text(size=30))+
#     guides(colour = guide_legend(override.aes = list(size=12)))
#   return(x)
# })

louvain_plots_rev = lapply(louvain_plots, FUN = function(x){
  x = x + guides(colour = guide_legend(override.aes = list(size=12)))+
    theme(plot.title = element_text(size=32),legend.text=element_text(size=30))
  return(x)
})

png("feature_zscore_group_map_by_stage.png",width = 70, height = 20, units = "in", res = 400)
ggarrange(
  ggarrange(plotlist=louvain_plots_rev, nrow = 3, common.legend = T, legend = "left"),
  ggarrange(plotlist=feature_plots, nrow = 3, ncol = 11, common.legend = T, legend = "right")
  ,ncol = 2
  ,widths = c(1.2, 11)
)
dev.off()


# plts = FeaturePlot(integrated_cancer_cells,
#                    features = paste0(signature_list_names, "_Normalize_Zscore"),
#                    combine = F,
#                    min.cutoff = -1,
#                    max.cutoff = 1)
# plts = lapply(plts, function(x){
#   return(x + scale_colour_gradient2(
#     high = "gold"
#     ,low = "blue"
#     ,mid = "grey"
#     ,limits = c(-1, 1)))
# })



  

######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################



#plot through split view
all_plots = list()
for(x in c("metastasis", "non-metastasis")){
  integrated_cancer_cells_stage = subset(integrated_cancer_cells, subset = metastatic_status == x)
  target_res = 0.2
  target_res_para = paste0("integrated_snn_res.", target_res)
  all_types = unique(integrated_cancer_cells_stage[[target_res_para]])[[target_res_para]]
  all_types = sort(all_types)
  all_colors = colorRampPalette(palette_color)(length(all_types))
  names(all_colors) = all_types
  cluster_maps = DimPlot(integrated_cancer_cells_stage, reduction = "umap",
                         group.by = "integrated_snn_res.0.2")+
    scale_color_manual(
      breaks = all_types, 
      values=all_colors[all_types])+ggtitle(x)
  all_plots = append(all_plots, list(cluster_maps))
  
  feature_plots = FeaturePlot(integrated_cancer_cells_stage, 
                              features = signature_list_names,
                              ncol = 4
                              ,min.cutoff = -0.5
                              ,max.cutoff = 0.5
                              ,combine = F
  )
  feature_plots = lapply(feature_plots, function(x){
    return(x + scale_colour_gradient2(
      high = "gold"
      ,low = "blue"
      ,mid = "grey"
      ,limits = c(-0.5, 0.5)))
  })
  all_plots = append(all_plots, feature_plots)
}


png("split_by_metastasis.png",width = 45, height = 8, units = "in", res = 400)
ggarrange(
  plotlist = all_plots
  ,ncol = 11
  ,nrow = 2
)
dev.off()

#######################################################
#Umaps
#######################################################

#subclustering_umap
target_res = 0.2
target_res_para = paste0("integrated_snn_res.", target_res)
all_types = unique(integrated_cancer_cells[[target_res_para]])[[target_res_para]]
all_types = sort(all_types)
all_colors = colorRampPalette(palette_color)(length(all_types))
names(all_colors) = all_types
png("subclustering_umap.png",9,9, units = "in", res = 300)
DimPlot(integrated_cancer_cells, reduction = "umap", 
        group.by = "integrated_snn_res.0.2")+
  scale_color_manual(breaks = all_types, values=all_colors[all_types])+ggtitle("Louvain clusters")
dev.off()



#subclustering_umap
target_res = 0.2
target_res_para = paste0("integrated_snn_res.", target_res)
all_types = unique(integrated_cancer_cells[[target_res_para]])[[target_res_para]]
all_types = sort(all_types)
all_colors = colorRampPalette(palette_color)(length(all_types))
names(all_colors) = all_types
signature_list = list("MT" = MT
                      ,"Complex_I"=Complex_I
                      ,"Complex_II"=Complex_II
                      ,"Complex_III"=Complex_III
                      ,"Complex_IV"=Complex_IV
                      ,"Complex_V"=Complex_V
                      ,"TCA"=TCA
                      ,"glycolysis"=glycolysis
                      ,"glycolyticAndTCA"=glycolyticAndTCA
                      ,"MAS"=MAS
                      ,"HIF1A" = HIF1A
                      ,"CU"=CU
                      ,"EMT"=EMT
                      ,"HIF1A_2A"=HIF1A_2A
                      ,"HLA"=HLA
                      ,"MRP"=MRP
                      ,"MRP_positive"=MRP_positive
                      ,"MRP_negative"=MRP_negative
)
signature_list_names = names(signature_list)
integrated_cancer_cells = AddModuleScore(integrated_cancer_cells, 
                                        features = signature_list, 
                                        assay = "RNA", 
                                        name = signature_list_names)
integrated_cancer_cells@meta.data[,signature_list_names] = 
  integrated_cancer_cells@meta.data[,paste0(signature_list_names, 1:length(signature_list_names))]

p1 = DimPlot(integrated_cancer_cells, reduction = "umap", group.by = target_res_para)+
  scale_color_manual(breaks = all_types, values=all_colors[all_types])+ggtitle("Louvain cluster")
p2 = DimPlot(integrated_cancer_cells, reduction = "umap", group.by = "patient")
p3 = FeaturePlot(integrated_cancer_cells, 
                 features = signature_list_names, 
                 ncol = 4
                 ,min.cutoff = -0.5
                 ,max.cutoff = 0.5
                 ,combine = F
)

p3 = lapply(p3, function(x){
  return(x + scale_colour_gradient2(
                                    high = "gold"
                                   ,low = "blue"
                                   ,mid = "grey"
                                   ,limits = c(-0.5, 0.5)))
})
p4 = DimPlot(integrated_cancer_cells, reduction = "umap", group.by = "stage")
plot_List = append(list(p1,p2,p4), p3)
png("module_score_extension_markers.png",25,10, units = "in", res = 400)
ggarrange(
  plotlist = plot_List
  ,ncol = 7
  ,nrow = 3
)
dev.off()


#add raw count data for composite markers
AddRawCount = function(object, gene_list){
  mat = object[["RNA"]]@counts
  result_list = list()
  for(x in 1:length(gene_list)){
    selected_genes = gene_list[[x]]
    marker_tag = names(gene_list[x])
    is_in_mat = selected_genes %in% rownames(mat)
    if(!all(is_in_mat)){
      not_available_genes = selected_genes[!is_in_mat]
      print(paste0(paste(not_available_genes, collapse = ","), 
                   " is not available for analysis"))
    }
    selected_genes = selected_genes[is_in_mat]
    selected_mat = mat[selected_genes, ]
    composite_marker_exp = as.numeric(apply(selected_mat, MARGIN = 2, FUN = mean))
    marker_tag = paste0(marker_tag,"_counts")
    result_list[[marker_tag]] = composite_marker_exp
  }
  result_list = as.data.frame(result_list)
  rownames(result_list) = colnames(x=object)
  object[[colnames(x = result_list)]] <- result_list
  return(object)
}

integrated_cancer_cells = AddRawCount(integrated_cancer_cells, signature_list)

p1 = DimPlot(integrated_cancer_cells, reduction = "umap", group.by = "integrated_snn_res.0.2")+
  scale_color_manual(breaks = all_types, values=all_colors[all_types])+ggtitle("Louvain cluster")
p2 = DimPlot(integrated_cancer_cells, reduction = "umap", group.by = "patient")
p3 = FeaturePlot(integrated_cancer_cells, 
                 features = paste0(signature_list_names, "_counts") 
                 # ,ncol = 5
                 # ,min.cutoff = -0.5
                 # ,max.cutoff = 0.5
                 ,combine = F
)

p3 = lapply(p3, function(x){
  return(x + scale_colour_gradient2(
    high = "yellow"
    ,mid = "gold4"
    ,low = "gray83"
    # ,midpoint = 0.75
    # ,limits = c(0, 1.5)
  ))
})
p4 = DimPlot(integrated_cancer_cells, reduction = "umap", group.by = "stage")
plot_List = append(list(p1,p2,p4), p3)
png("raw_count_umaps_extension_markers.png",25,10, units = "in", res = 400)
ggarrange(
  plotlist = plot_List
  ,ncol = 7
  ,nrow = 3
)
dev.off()


#add normalize count to object
AddNormalizeCount = function(object, gene_list){
  mat = object[["RNA"]]@data
  result_list = list()
  for(x in 1:length(gene_list)){
    selected_genes = gene_list[[x]]
    marker_tag = names(gene_list[x])
    is_in_mat = selected_genes %in% rownames(mat)
    if(!all(is_in_mat)){
      not_available_genes = selected_genes[!is_in_mat]
      print(paste0(paste(not_available_genes, collapse = ","), 
                   " is not available for analysis"))
    }
    selected_genes = selected_genes[is_in_mat]
    selected_mat = mat[selected_genes, ]
    composite_marker_exp = as.numeric(apply(selected_mat, MARGIN = 2, FUN = mean))
    marker_tag = paste0(marker_tag,"_normalize")
    result_list[[marker_tag]] = composite_marker_exp
  }
  result_list = as.data.frame(result_list)
  rownames(result_list) = colnames(x=object)
  object[[colnames(x = result_list)]] <- result_list
  return(object)
}

integrated_cancer_cells = AddNormalizeCount(integrated_cancer_cells, signature_list)

p1 = DimPlot(integrated_cancer_cells, reduction = "umap", group.by = "integrated_snn_res.0.2")+
  scale_color_manual(breaks = all_types, values=all_colors[all_types])+ggtitle("Louvain cluster")
p2 = DimPlot(integrated_cancer_cells, reduction = "umap", group.by = "patient")
p3 = FeaturePlot(integrated_cancer_cells, 
                 features = paste0(signature_list_names, "_normalize") 
                 # ,ncol = 5
                 # ,min.cutoff = 0
                 # ,max.cutoff = 1.5
                 ,combine = F
)


p3 = lapply(p3, function(x){
  return(x + scale_colour_gradient2(
      high = "yellow"
      ,mid = "gold4"
      ,low = "gray83"
      # ,midpoint = 0.75
      # ,limits = c(0, 1.5)
    ))
})
p4 = DimPlot(integrated_cancer_cells, reduction = "umap", group.by = "stage")

plot_List = append(list(p1,p2,p4), p3)
png("normalize_no_scale_count_umaps_extension_markers.png",25,10, units = "in", res = 400)
ggarrange(
  plotlist = plot_List
  ,ncol = 7
  ,nrow = 3
)
dev.off()









integrated_cancer_cells$Complex_II_group = 
  sapply(integrated_cancer_cells$Complex_II, FUN = function(x){
    if(x < -0.25){
      return("lower than -0.25")
    }else if(x > 0.25){
      return("higher than 0.25")
    }else{
      return("between -0.25 and 0.25")
    }
  })

DimPlot(integrated_cancer_cells, group.by = "stage", 
        split.by = "Complex_II_group", reduction = "umap")+ggtitle("")
DimPlot(integrated_cancer_cells, group.by = "stage", 
        split.by = "Complex_I_group", reduction = "umap")+
        ggtitle("Complex I stratification")

FeaturePlot(integrated_cancer_cells, features = "Complex_I_group",
            split.by = "stage",
            reduction = "umap")+ggtitle("")

#split by stage umaps


png("test.png",width = 16, height = 20, units = "in", res = 400)
FeaturePlot(integrated_cancer_cells, 
            features = CoQ10
            ,ncol=2
            ,min.cutoff = -0.5
            ,max.cutoff = 0.5
            ,split.by = "stage"
)
dev.off()

condition = "pT1b"
selected_cells = integrated_cancer_cells$stage==condition
hist(integrated_cancer_cells@meta.data[selected_cells, "Complex_I"], breaks = 100,
     main = paste0(condition, "_Complex_III"))

condition = "metastasis"
selected_cells = integrated_cancer_cells$metastatic_status==condition
hist(integrated_cancer_cells@meta.data[selected_cells, "Complex_I"], breaks = 100,
     main = paste0(condition, "_Complex_III"))






#######################################################
#Heatmaps
#######################################################

# module score heatmap
mat = t(integrated_cancer_cells@meta.data[,c("hla_score", "mrp_score", "nduf_score", 
                                            "mrp_nduf_score", "emt_score", "pt_score")])
mat = data.frame(mat, check.names = F)
order_df = data.frame(cbind(seurat_cluster = as.character(integrated_cancer_cells$integrated_snn_res.0.2),
                            pt = integrated_cancer_cells$patient, stage = integrated_cancer_cells$stage))
order_df = arrange(order_df, seurat_cluster, pt, stage)
mat_new_order = mat[,rownames(order_df)]

seurat_cluster_colors = colorRampPalette(palette_color)(length(unique(order_df$seurat_cluster)))
names(seurat_cluster_colors) = unique(order_df$seurat_cluster)
patient_colors = colorRampPalette(palette_color)(length(unique(integrated_cancer_cells$patient)))
names(patient_colors) = unique(integrated_cancer_cells$patient)
stage_colors = colorRampPalette(c("red", "green"))(length(unique(integrated_cancer_cells$stage)))
names(stage_colors) = unique(integrated_cancer_cells$stage)
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
DefaultAssay(integrated_cancer_cells) = "RNA"
integrated_cancer_cells = NormalizeData(integrated_cancer_cells)
cohort = integrated_cancer_cells
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


