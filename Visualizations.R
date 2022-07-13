source("call_libraries.R")
source("Initialization.R")
integrated_cancer_cells = readRDS("integrated_cancer_cells.rds")
DefaultAssay(integrated_cancer_cells) = "RNA"
integrated_cancer_cells = NormalizeData(integrated_cancer_cells)
#add metastatic status for each cell
integrated_cancer_cells$stage = 
  ifelse(integrated_cancer_cells$patient=="Patient5", "Metastatic", integrated_cancer_cells$stage)
#define groups tag for different algs
module_score_tag = "_ModuleScore"
z_score_tag = "_ZScore"
normalized_count_tag = "_NormalizedCount"
raw_count_tag = "_RawCount"

#initialize signature list
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
# add Seurat module score based on pre-defined signature list
signature_list_names = names(signature_list)
#############################################################################################
#############################################################################################
#############################################################################################

# Module score
integrated_cancer_cells = AddModuleScore(integrated_cancer_cells, 
                                        features = signature_list, 
                                        assay = "RNA", 
                                        name = signature_list_names)
integrated_cancer_cells@meta.data[,paste0(signature_list_names, module_score_tag)] = 
integrated_cancer_cells@meta.data[,paste0(signature_list_names, 1:length(signature_list_names))]
#remove redundant features
remove_column_index = match(paste0(signature_list_names, 1:length(signature_list_names)),
                            names(integrated_cancer_cells@meta.data))
integrated_cancer_cells@meta.data = 
integrated_cancer_cells@meta.data[,-remove_column_index]

#create grouping variables for module scores of each pathway signature
integrated_cancer_cells@meta.data[,paste0(signature_list_names, module_score_tag, "_group")] = 
lapply(integrated_cancer_cells@meta.data[,paste0(signature_list_names, module_score_tag)], FUN = function(vx){
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
target_res = 0.2
target_res_para = paste0("integrated_snn_res.", target_res)
for(x in c("pT1b","pT3a","Metastatic")){
  integrated_cancer_cells_stage = subset(integrated_cancer_cells, subset = stage == x)
  all_types = unique(integrated_cancer_cells_stage[[target_res_para]])[[target_res_para]]
  all_types = sort(all_types)
  all_colors = colorRampPalette(palette_color)(length(all_types))
  names(all_colors) = all_types
  cluster_maps = DimPlot(integrated_cancer_cells_stage, reduction = "umap",
                         group.by = target_res_para)+
    scale_color_manual(
      breaks = all_types, 
      values=all_colors[all_types])+ggtitle(x)+
    theme(plot.title = element_text(size=15), legend.text=element_text(size=15))+
    guides(colour = guide_legend(override.aes = list(size=6)))
  louvain_plots = append(louvain_plots, list(cluster_maps))
  
  all_types = group_order
  all_types = sort(all_types)
  all_colors = c("grey", "cyan", "blue", "gold", "gold3")
  names(all_colors) = all_types
  
  #sort legend
  all_types_factor = factor(all_types, levels = group_order)
  all_types = as.character(sort(all_types_factor))
  all_colors = all_colors[all_types]
  
  # make feature plots
  feature_group_plots = lapply(signature_list_names, FUN = function(signature){
    plot_title = signature
    if(x!="pT1b"){
      plot_title = ""
    }
    p = DimPlot(integrated_cancer_cells_stage, reduction = "umap",
            group.by = paste0(signature, module_score_tag, "_group"))+
      scale_color_manual(
        breaks = all_types, 
        values=all_colors[all_types])+
      theme(plot.title = element_text(size=15), legend.text=element_text(size=15))+
      ggtitle(plot_title)+
      guides(colour = guide_legend(override.aes = list(size=6)))
    return(p)
  })
  feature_plots = append(feature_plots, feature_group_plots)
  
}

# plot out feature umap
png("module_score_u_map_by_stage.png",width = 100, height = 30, units = "cm", res = 300)
ggarrange(
  ggarrange(plotlist=louvain_plots, nrow = 3, common.legend = T, legend = "left"),
  ggarrange(plotlist=feature_plots, nrow = 3, ncol = 11, common.legend = T, legend = "right")
  ,ncol = 2
  ,widths = c(1.2, 11)
)
dev.off()

# plot out violin
base = 5
interval = 0.8
all_vln_plots = lapply(signature_list_names, function(x){
  p_table = compare_means(as.formula(paste0(paste0(x, module_score_tag), "~", "stage")), 
                          data = integrated_cancer_cells@meta.data, method = "t.test")
  p_table = data.frame(p_table)[,c(2,3,5)]
  p_table$y.position = seq(base,base+nrow(p_table)*interval,interval)[1:nrow(p_table)]
  p_table = as_tibble(p_table)
  reordered_stage = factor(integrated_cancer_cells$stage, levels = c("pT1b", "pT3a", "Metastatic"))
  score = integrated_cancer_cells@meta.data[,paste0(x, module_score_tag)] 
  ggplot(integrated_cancer_cells@meta.data, aes(x=reordered_stage, y=score))+
    geom_violin(aes(fill=reordered_stage), draw_quantiles = 0.5)+
    guides(fill=guide_legend(title="stage"))+
    scale_fill_manual(values = c("#00BE0E", "#489DFF", "#FF6C67"))+
    add_pvalue(p_table)+
    ylim(-1,base+2)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
      axis.title.x = element_blank(),
      plot.title = element_text(color = "blue", size = 12, face = "bold", hjust = 0.5))+
    ylab(paste0(x, module_score_tag))+
    ggtitle(x)
})

png("module_score_vln_by_stage.png",width = 25, height = 12, units = "in", res = 300)
ggarrange(plotlist = all_vln_plots,ncol = 6, nrow = 2, common.legend = T, legend = "right")
dev.off()

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

# define function to add z-score count of super gene to object
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
    marker_tag = paste0(marker_tag, z_score_tag)
    result_list[[marker_tag]] = composite_marker_exp
  }
  result_list = as.data.frame(result_list)
  rownames(result_list) = colnames(x=object)
  object[[colnames(x = result_list)]] <- result_list
  return(object)
}

integrated_cancer_cells = AddScaledNormalizeCount(integrated_cancer_cells, signature_list)

#define group
integrated_cancer_cells@meta.data[,paste0(signature_list_names, z_score_tag, "_group")] = 
  lapply(integrated_cancer_cells@meta.data[,paste0(signature_list_names, z_score_tag)], FUN = function(vx){
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
target_res = 0.2
target_res_para = paste0("integrated_snn_res.", target_res)
for(x in c("pT1b","pT3a","Metastatic")){
  integrated_cancer_cells_stage = subset(integrated_cancer_cells, subset = stage == x)
  all_types = unique(integrated_cancer_cells_stage[[target_res_para]])[[target_res_para]]
  all_types = sort(all_types)
  all_colors = colorRampPalette(palette_color)(length(all_types))
  names(all_colors) = all_types
  cluster_maps = DimPlot(integrated_cancer_cells_stage, reduction = "umap",
                         group.by = "integrated_snn_res.0.2")+
    scale_color_manual(
      breaks = all_types, 
      values=all_colors[all_types])+ggtitle(x)+
    theme(plot.title = element_text(size=32),legend.text=element_text(size=30))+
    guides(colour = guide_legend(override.aes = list(size=12)))
  louvain_plots = append(louvain_plots, list(cluster_maps))
  
  all_types = group_order
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
    p = DimPlot(integrated_cancer_cells_stage, reduction = "umap",
                group.by = paste0(signature, z_score_tag, "_group"))+
      scale_color_manual(
        breaks = all_types, 
        values=all_colors[all_types])+
      theme(plot.title = element_text(size=32),legend.text=element_text(size=30))+
      guides(colour = guide_legend(override.aes = list(size=12)))+
      ggtitle(plot_title)
    return(p)
  })
  
  feature_plots = append(feature_plots, feature_group_plots)
}

# plot out the z-score umap
png("feature_zscore_group_map_by_stage.png",width = 70, height = 20, units = "in", res = 300)
ggarrange(
  ggarrange(plotlist=louvain_plots, nrow = 3, common.legend = T, legend = "left"),
  ggarrange(plotlist=feature_plots, nrow = 3, ncol = 11, common.legend = T, legend = "right")
  ,ncol = 2
  ,widths = c(1.2, 11)
)
dev.off()

# plot out violin
base = 7
interval = 0.7
all_vln_plots = lapply(signature_list_names, function(x){
  p_table = compare_means(as.formula(paste0(paste0(x, z_score_tag), "~", "stage")), 
                          data = integrated_cancer_cells@meta.data, method = "t.test")
  p_table = data.frame(p_table)[,c(2,3,5)]
  p_table$y.position = seq(base,base+nrow(p_table)*interval,interval)[1:nrow(p_table)]
  p_table = as_tibble(p_table)
  reordered_stage = factor(integrated_cancer_cells$stage, levels = c("pT1b", "pT3a", "Metastatic"))
  score = integrated_cancer_cells@meta.data[,paste0(x, z_score_tag)] 
  ggplot(integrated_cancer_cells@meta.data, aes(x=reordered_stage, y=score))+
    geom_violin(aes(fill=reordered_stage), draw_quantiles = 0.5)+
    guides(fill=guide_legend(title="stage"))+
    scale_fill_manual(values = c("#00BE0E", "#489DFF", "#FF6C67"))+
    add_pvalue(p_table)+
    ylim(-1,base+2)+
    theme(
      axis.title.x = element_blank(),
      plot.title = element_text(color = "blue", size = 12, face = "bold", hjust = 0.5))+
    ylab(paste0(x, z_score_tag))+
    ggtitle(x)
})

png("z_score_vln_by_stage.png",width = 25, height = 12, units = "in", res = 300)
ggarrange(plotlist = all_vln_plots,ncol = 6, nrow = 2, common.legend = T, legend = "right")
dev.off()


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

#add normalize count to object
AddNormalizedCount = function(object, gene_list){
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
    marker_tag = paste0(marker_tag,normalized_count_tag)
    result_list[[marker_tag]] = composite_marker_exp
  }
  result_list = as.data.frame(result_list)
  rownames(result_list) = colnames(x=object)
  object[[colnames(x = result_list)]] <- result_list
  return(object)
}

integrated_cancer_cells = AddNormalizedCount(integrated_cancer_cells, signature_list)

#define group
integrated_cancer_cells@meta.data[,paste0(signature_list_names, normalized_count_tag, "_group")] = 
  lapply(integrated_cancer_cells@meta.data[,paste0(signature_list_names, normalized_count_tag)], FUN = function(vx){
    vx = sapply(vx, FUN = function(x){
      if(x < 0.25){
        return("<0.25")
      }else if(x >= 0.25 & x < 0.5){
        return("0.25 to 0.5")
      }else if(x >= 0.5 & x < 0.75){
        return("0.5 to 0.75")
      }else if(x >= 0.75 & x < 1){
        return("0.75 to 1")
      }else{
        return(">=1")
      }
    })
  })

# create feature plots with scaling
louvain_plots = list()
feature_plots = list()
group_order <- c(">=1", "0.75 to 1", "0.5 to 0.75", "0.25 to 0.5", "<0.25")
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
    theme(plot.title = element_text(size=32),legend.text=element_text(size=30))+
    guides(colour = guide_legend(override.aes = list(size=12)))
  louvain_plots = append(louvain_plots, list(cluster_maps))
  
  all_types = group_order
  #all_types = sort(all_types)
  all_colors = c("gold","gold3", "grey", "cyan", "blue")
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
    p = DimPlot(integrated_cancer_cells_stage, reduction = "umap",
                group.by = paste0(signature, normalized_count_tag, "_group"))+
      scale_color_manual(
        breaks = all_types, 
        values=all_colors[all_types])+
      theme(plot.title = element_text(size=32),legend.text=element_text(size=30))+
      guides(colour = guide_legend(override.aes = list(size=12)))+
      ggtitle(plot_title)
    return(p)
  })
  
  feature_plots = append(feature_plots, feature_group_plots)
}

# plot out the z-score umap
png("feature_normalized_count_group_map_by_stage.png",width = 70, height = 20, units = "in", res = 300)
ggarrange(
  ggarrange(plotlist=louvain_plots, nrow = 3, common.legend = T, legend = "left"),
  ggarrange(plotlist=feature_plots, nrow = 3, ncol = 11, common.legend = T, legend = "right")
  ,ncol = 2
  ,widths = c(1.2, 11)
)
dev.off()

# plot out violin
base = 5.5
interval = 0.6
all_vln_plots = lapply(signature_list_names, function(x){
  p_table = compare_means(as.formula(paste0(paste0(x, normalized_count_tag), "~", "stage")), 
                          data = integrated_cancer_cells@meta.data, method = "t.test")
  p_table = data.frame(p_table)[,c(2,3,5)]
  p_table$y.position = seq(base,base+nrow(p_table)*interval,interval)[1:nrow(p_table)]
  p_table = as_tibble(p_table)
  reordered_stage = factor(integrated_cancer_cells$stage, levels = c("pT1b", "pT3a", "Metastatic"))
  score = integrated_cancer_cells@meta.data[,paste0(x, normalized_count_tag)] 
  ggplot(integrated_cancer_cells@meta.data, aes(x=reordered_stage, y=score))+
    geom_violin(aes(fill=reordered_stage), draw_quantiles = 0.5)+
    guides(fill=guide_legend(title="stage"))+
    scale_fill_manual(values = c("#00BE0E", "#489DFF", "#FF6C67"))+
    add_pvalue(p_table)+
    ylim(-1,base+2)+
    theme(
      axis.title.x = element_blank(),
      plot.title = element_text(color = "blue", size = 12, face = "bold", hjust = 0.5))+
    ylab(paste0(x, normalized_count_tag))+
    ggtitle(x)
})

png("normalized_count_vln_by_stage.png",width = 25, height = 12, units = "in", res = 300)
ggarrange(plotlist = all_vln_plots,ncol = 6, nrow = 2, common.legend = T, legend = "right")
dev.off()

#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
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
    marker_tag = paste0(marker_tag,raw_count_tag)
    result_list[[marker_tag]] = composite_marker_exp
  }
  result_list = as.data.frame(result_list)
  rownames(result_list) = colnames(x=object)
  object[[colnames(x = result_list)]] <- result_list
  return(object)
}
integrated_cancer_cells = AddRawCount(integrated_cancer_cells, signature_list)

#define group
integrated_cancer_cells@meta.data[,paste0(signature_list_names, raw_count_tag, "_group")] = 
  lapply(integrated_cancer_cells@meta.data[,paste0(signature_list_names, raw_count_tag)], FUN = function(vx){
    vx = sapply(vx, FUN = function(x){
      if(x < 2.5){
        return("<2.5")
      }else if(x >= 2.5 & x < 5){
        return("2.5 to 5")
      }else if(x >= 5 & x < 7.5){
        return("5 to 7.5")
      }else if(x >= 7.5 & x < 10){
        return("7.5 to 10")
      }else{
        return(">=10")
      }
    })
  })

#draw feature plots with scaling
louvain_plots = list()
feature_plots = list()
group_order <- c(">=10", "7.5 to 10", "5 to 7.5", "2.5 to 5", "<2.5")
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
    theme(plot.title = element_text(size=32),legend.text=element_text(size=30))+
    guides(colour = guide_legend(override.aes = list(size=12)))
  louvain_plots = append(louvain_plots, list(cluster_maps))
  
  all_types = group_order
  all_types = sort(all_types)
  all_colors = c("blue", "gold", "grey", "gold3", "cyan")
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
    p = DimPlot(integrated_cancer_cells_stage, reduction = "umap",
                group.by = paste0(signature, raw_count_tag, "_group"))+
      scale_color_manual(
        breaks = all_types, 
        values=all_colors[all_types])+
      theme(plot.title = element_text(size=32),legend.text=element_text(size=30))+
      guides(colour = guide_legend(override.aes = list(size=12)))+
      ggtitle(plot_title)
    return(p)
  })
  
  feature_plots = append(feature_plots, feature_group_plots)
}

png("feature_raw_count_group_map_by_stage.png",width = 70, height = 20, units = "in", res = 300)
ggarrange(
  ggarrange(plotlist=louvain_plots, nrow = 3, common.legend = T, legend = "left"),
  ggarrange(plotlist=feature_plots, nrow = 3, ncol = 11, common.legend = T, legend = "right")
  ,ncol = 2
  ,widths = c(1.2, 11)
)
dev.off()

# plot out violin
base = 50
interval = 1
all_vln_plots = lapply(signature_list_names, function(x){
  p_table = compare_means(as.formula(paste0(paste0(x, raw_count_tag), "~", "stage")), 
                          data = integrated_cancer_cells@meta.data, method = "t.test")
  p_table = data.frame(p_table)[,c(2,3,5)]
  p_table$y.position = seq(base,base+nrow(p_table)*interval,interval)[1:nrow(p_table)]
  p_table = as_tibble(p_table)
  reordered_stage = factor(integrated_cancer_cells$stage, levels = c("pT1b", "pT3a", "Metastatic"))
  score = integrated_cancer_cells@meta.data[,paste0(x, raw_count_tag)] 
  ggplot(integrated_cancer_cells@meta.data, aes(x=reordered_stage, y=score))+
    geom_violin(aes(fill=reordered_stage), draw_quantiles = 0.5)+
    guides(fill=guide_legend(title="stage"))+
    scale_fill_manual(values = c("#00BE0E", "#489DFF", "#FF6C67"))+
    add_pvalue(p_table)+
    ylim(-1,base+2)+
    theme(
      axis.title.x = element_blank(),
      plot.title = element_text(color = "blue", size = 12, face = "bold", hjust = 0.5))+
    ylab(paste0(x, raw_count_tag))+
    ggtitle(x)
})

png("raw_count_vln_by_stage.png",width = 25, height = 15, units = "in", res = 300)
ggarrange(plotlist = all_vln_plots,ncol = 6, nrow = 2)
dev.off()

#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################


#####################################################################################################
#pairwise-correlation plot
#####################################################################################################

Pairwise_Heatmap = function(exp_df){
  col_fun = colorRampPalette(c("blue", 'white', 'red'))(20)
  p.mat <- cor.mtest(exp_df)
  cor_matrix = cor(exp_df, method="pearson")
  
  #draw the heatmap obj and return the real heatmap
  heatmap_obj = corrplot(cor_matrix, p.mat = p.mat$p, col=col_fun, type="upper", insig = "blank", tl.cex=0.8,
                         cl.cex = 1)
  return(heatmap_obj)
  
  
  #***********Used For Debug**************
  # heatmap_plot = draw(heatmap_obj)
  # return(invisible(NULL))
  # heatmap_obj = corrplot::corrplot(matrix, p.mat = res1$p, col = col_fun,
  # type="upper", pch.col = "white", insig = "label_sig",order="hclust", pch.cex=2, tl.pos="d",tl.srt=60)
  
}

mat = integrated_cancer_cells[["RNA"]]@data
plt_list = lapply(signature_list_names, FUN = function(x){
  selected_genes = signature_list[[x]]
  is_in_mat = selected_genes %in% rownames(mat)
  if(!all(is_in_mat)){
    not_available_genes = selected_genes[!is_in_mat]
    print(paste0(paste(not_available_genes, collapse = ","), 
                 " is not available for analysis"))
  }
  selected_genes = selected_genes[is_in_mat]
  is_all_zero = sapply(selected_genes, FUN = function(x){
    return(all(as.numeric(mat[x,])==0))
  })
  selected_genes = selected_genes[!is_all_zero]
  selected_mat = mat[selected_genes, ]
  selected_mat = data.frame(t(as.matrix(selected_mat)))
  col_fun = colorRampPalette(c("blue", 'white', 'red'))(20)
  p.mat <- cor_pmat(selected_mat, method="spearman")
  cor_matrix = cor(selected_mat, method="spearman")
  cor_matrix = data.frame(cor_matrix, 
                          row.names = selected_genes, check.rows = F, 
                          check.names = F)
  names(cor_matrix) = selected_genes
  return(ggcorrplot(cor_matrix, hc.order = TRUE, type = "upper",p.mat = p.mat))
})

plt_list_rev = lapply(plt_list, FUN = function(x){
  x = x + theme(plot.margin = margin(t = 0,  # Top margin
                                     r = 0,  # Right margin
                                     b = 0,  # Bottom margin
                                     l = 0))
})

png("pairwise_correlation.png", res = 300, units = "cm", height = 60, width = 100)
ggarrange(plotlist = plt_list_rev, nrow = 2, ncol = 6)
dev.off()







#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#####################################################################################################
#Heatmaps
#####################################################################################################

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


