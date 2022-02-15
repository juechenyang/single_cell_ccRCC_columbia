library(SeuratDisk)
export_to_cellbygene = function(seurat_object){
  SaveH5Seurat(seurat_object, filename = "temp_cellby_gene.h5Seurat")
  Convert("temp_cellby_gene.h5Seurat", dest = "h5ad")
}