#' annotate_clusters.R
#'
#' This function adds manual celltype annotations to a seurat object
#'
#' @param seurat_obj Input seurate object
#' @param resolution Resolution to be used for clustering
#' @param annotations Named vector of cell types for each cluster
#' @return Seurat object
#' @export


annotate_clusters <- function(seurat_obj=NULL, resolution=NULL, annotations=NULL){
  cluster_col <- paste0("ADT_snn_res.", resolution)
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop("Couldn't find cluster column: ", cluster_col)
  }
  Idents(seurat_obj) <- cluster_col
  
  # map clusters -> labels
  lab <- annotations[as.character(Idents(seurat_obj))]
  lab[is.na(lab)] <- "Unknown"              # optional default
  names(lab) <- colnames(seurat_obj)        # <<< key line
  
  seurat_obj$final_annotation <- lab
  Idents(seurat_obj) <- "final_annotation"
  return(seurat_obj)
}
