options(future.globals.maxSize = 20000 * 1024^2)
#function to integrate a splited seurat datasets with split key and 
#ref as the dataset who has the largest number of cells
Integrate_seurat_obj_by_key = function(obj, key, method_SCT = T, method_rpca = F,
                                       n.features=2000, n.anchors=5){
  obj_list = SplitObject(obj, split.by = key)
  if(method_SCT){
    if(method_rpca){
      obj_list <- lapply(X = obj_list, FUN = SCTransform, method = "glmGamPoi")
      features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = n.features)
      obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)
      obj_list <- lapply(X = obj_list, FUN = RunPCA, features = features)
      anchors <- FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT", 
                                        anchor.features = features, reduction = "rpca",
                                        k.anchor = n.anchors)
      result_data <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = F)
    }else{
      obj_list <- lapply(X = obj_list, FUN = SCTransform)
      features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = n.features)
      obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)
      anchors <- FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT", 
                                        anchor.features = features,
                                        k.anchor = n.anchors)
      result_data <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = F)
    }
  }else{
    
    if(method_rpca){
      obj_list <- lapply(X = obj_list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = n.features)
      })
      features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = n.features)
      obj_list <- lapply(X = obj_list, FUN = function(x) {
        x <- ScaleData(x, features = features, verbose = FALSE)
        x <- RunPCA(x, features = features, verbose = FALSE)
      })
      anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features, 
                                        reduction = "rpca", 
                                        k.anchor = n.anchors)
    }else{
      obj_list <- lapply(X = obj_list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = n.features)
      })
      features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = n.features)
      anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features, 
                                        k.anchor = n.anchors)
    }
    result_data <- IntegrateData(anchorset = anchors, verbose = F)
    DefaultAssay(result_data) <- "integrated"
    result_data <- ScaleData(result_data, verbose = FALSE)
  }
  
  return(result_data)
}
