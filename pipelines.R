preprocess_standard_pipelines = function(obj, n_variable_features = 3000){
  # Normalize data
  obj = NormalizeData(obj)
  # FindVariableFeatures
  obj = FindVariableFeatures(obj, nfeature = n_variable_features)
  # Scale data
  obj = ScaleData(obj)
  
  return(obj)
}

preprocess_SCT = function(obj, n_variable_features = 3000){
  obj = SCTransform(obj, variable.features.n = n_variable_features, 
                    vars.to.regress = "percent.mt", verbose = F)
  return(obj)
}

run_pca_and_umap = function(obj, n_of_pc = 30){
  # Run PCA
  obj <- RunPCA(obj, npcs = n_of_pc, verbose = F)
  # Run Umap
  obj <- RunUMAP(obj, dims = 1:n_of_pc, reduction = "pca", umap.method="umap-learn",
                 metric="correlation", verbose = F)
  return(obj)
}

find_neighbours_and_clusters = function(obj, n_of_pc = 30, res = seq(0.1,1,0.1)){
  obj <- FindNeighbors(obj, dims = 1:n_of_pc, verbose = F)
  obj <- FindClusters(obj, resolution = res)
  return(obj)
}

