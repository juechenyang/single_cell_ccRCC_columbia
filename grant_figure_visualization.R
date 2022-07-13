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
  ,"Glycolysis"=glycolysis
)

# add Seurat module score based on pre-defined signature list
signature_list_names = names(signature_list)
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
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),axis.title.y = element_blank(),
          plot.title = element_text(size=15), legend.text=element_text(size=15))+
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
    plot_title =str_replace_all(signature, "_", " ")
    if(x!="pT1b"){
      plot_title = ""
    }
    p = DimPlot(integrated_cancer_cells_stage, reduction = "umap",
                group.by = paste0(signature, module_score_tag, "_group"))+
      scale_color_manual(
        breaks = all_types, 
        values=all_colors[all_types])+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.title.x = element_blank(),axis.title.y = element_blank(),
            plot.title = element_text(size=15), legend.text=element_text(size=15),
            plot.margin = margin(t = 10,  # Top margin
                                 r = 10,  # Right margin
                                 b = 10,  # Bottom margin
                                 l = 10)
            )+
      ggtitle(plot_title)+
      guides(colour = guide_legend(override.aes = list(size=6)))
    return(p)
  })
  feature_plots = append(feature_plots, feature_group_plots)
  
}

# plot out feature umap
png("grant_module_score_umap.png",width = 60, height = 30, units = "cm", res = 300)
ggarrange(
  ggarrange(plotlist=louvain_plots, nrow = 3, common.legend = T, legend = "left"),
  ggarrange(plotlist=feature_plots, nrow = 3, ncol = 5, common.legend = T, legend = "right")
  ,ncol = 2
  ,widths = c(1, 5)
)
dev.off()

# plot out violin
base = 1.8
interval = 0.6
all_vln_plots = lapply(signature_list_names, function(x){
  plot_title = str_replace_all(x, "_", " ")
  p_table = compare_means(as.formula(paste0(paste0(x, module_score_tag), "~", "stage")), 
                          data = integrated_cancer_cells@meta.data, method = "t.test")
  p_table = data.frame(p_table)[,c(2,3,5)]
  if(x=="Complex_IV"){
    p_table[3,3] = round(p_table[3,3], 300)
  }
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
          axis.title.x = element_blank(),axis.title.y = element_blank(),
          legend.position="none",
          plot.title = element_blank(),
          plot.margin = margin(t = 10,  # Top margin
                               r = 20,  # Right margin
                               b = 10,  # Bottom margin
                               l = 21))+ # Left margin)+
    ylab(paste0(x, module_score_tag))
})

# png("grant_module_score_vln.png",width = 20, height = 20, units = "cm", res = 300)
# ggarrange(plotlist = all_vln_plots,ncol = 3, nrow = 2, common.legend = T, legend = "none")
# dev.off()

feature_and_vln = append(feature_plots, all_vln_plots)
png("aggregate.png",width = 40, height = 30, units = "cm", res = 300)
ggarrange(
  ggarrange(plotlist = louvain_plots, nrow = 4, ncol = 1, common.legend = T, legend = "left"),
  ggarrange(plotlist = feature_and_vln, nrow = 4, ncol = 5, common.legend = T, legend = "right"),
  ncol = 2,
  nrow = 1,
  widths = c(1,5)
)
dev.off()