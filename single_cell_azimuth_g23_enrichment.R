library(clusterProfiler)
library(tidyverse)


#*********************************************************************************************
#*load label cell funcions that generated from enrichR
# define cell type finder
get_cell_type_enrichR <- function(gene_list){
  source_python('CellTypeFinder.py')
  gene_str = paste(gene_list, collapse = "\n")
  result <- get_cell_type(gene_str)
  if(nrow(result)==0){
    return("Unknown")
  }else{
    return(result[1,"Term name"])
  }
}

label_cluster_enrichR = function(seurat_obj, top){
  top = data.frame(top)
  cell_types = c()
  top$cluster = as.character(top$cluster)
  for(i in seq(0, as.numeric(top[nrow(top), "cluster"]))){
    top_sub = top[top$cluster==as.character(i), ]
    gene_cell_type = get_cell_type_enrichR(top_sub$gene)
    x = stringr::str_split_fixed(gene_cell_type, "CL[0-9]+", 2)
    gene_cell_type = trimws(x[[1]])
    if(gene_cell_type %in% cell_types){
      cell_types = c(cell_types, paste0(gene_cell_type, "_", i))
    }else{
      cell_types = c(cell_types, gene_cell_type)
    }
  }
  names(cell_types) <- levels(seurat_obj)
  seurat_obj <- RenameIdents(seurat_obj, cell_types)
  return(seurat_obj)
}


#prepare azimuth with g23 lib
cell_type_lib = read.csv("Azimuth_kidney.csv", header = F)
#cell_type_lib = read.csv("Azimuth_Cell_Types_2021.csv", header = F)
# cell_type_lib = read.csv("CellMarker_Augmented_2021.csv", header = F)
# cell_type_lib = cell_type_lib[grepl(".+Kidney", cell_type_lib$V1), ]
rownames(cell_type_lib) = cell_type_lib[,1]
df = cell_type_lib[,-1]
df$geneID <-  do.call(Map, c(f = c, df))
df$term = rownames(df)
df = df %>%
  dplyr::select(term, geneID) %>%
  tidyr::unnest(cols = c(geneID))
df = df[df$geneID!="",]
#label cell type GSEA
get_cell_type_locally_gsea = function(gene_list){
  gsea_results <- GSEA(gene_list, TERM2GENE = df)
  gsea_results = data.frame(gsea_results)
  gsea_results = gsea_results[gsea_results$enrichmentScore>0, ]
  if(nrow(gsea_results)==0){
    return("Unknown")
  }else{
    gene_cell_type = gsea_results[1,"Description"]
    x = stringr::str_split_fixed(gene_cell_type, "CL[0-9]+", 2)
    gene_cell_type = trimws(x[[1]])
    return(gene_cell_type)
  }
}

get_cell_type_locally_or = function(gene_list){
  gsea_results <- enricher(gene_list, TERM2GENE = df)
  gsea_results = data.frame(gsea_results)
  if(nrow(gsea_results)==0){
    return("Unknown")
  }else{
    gsea_results = gsea_results[order(gsea_results[,"p.adjust"], decreasing = F), ]
    print(gsea_results[1,"Description"])
    print(gsea_results[1,"p.adjust"])
    return(gsea_results[1,"Description"])
  }
}
#GSEA enrichment
label_cluster_locally = function(seurat_obj, top, type="gsea"){
  cell_types = c()
  if(type == "gsea"){
    for(i in seq(0, max(as.numeric(Idents(seurat_obj)))-1)){
      fold_change = FoldChange(seurat_obj, ident.1 = i)
      gene_list = fold_change$avg_log2FC
      names(gene_list) = rownames(fold_change)
      gene_list = sort(gene_list, decreasing = T)
      gene_cell_type = get_cell_type_locally_gsea(gene_list)
      if(gene_cell_type %in% cell_types){
        cell_types = c(cell_types, paste0(gene_cell_type, "_", i))
      }else{
        cell_types = c(cell_types, gene_cell_type)
      }
    }
  }else{
    top = data.frame(top)
    top$cluster = as.character(top$cluster)
    for(i in seq(0, as.numeric(top[nrow(top), "cluster"]))){
      top_sub = top[top$cluster==as.character(i), ]
      gene_cell_type = get_cell_type_locally_or(top_sub$gene)
      x = stringr::str_split_fixed(gene_cell_type, "CL[0-9]+", 2)
      gene_cell_type = trimws(x[[1]])
      if(gene_cell_type %in% cell_types){
        cell_types = c(cell_types, paste0(gene_cell_type, "_", i))
      }else{
        cell_types = c(cell_types, gene_cell_type)
      }
    }
  }
  names(cell_types) <- levels(seurat_obj)
  seurat_obj <- RenameIdents(seurat_obj, cell_types)
  return(seurat_obj)
}