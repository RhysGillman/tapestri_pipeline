#' annotate_clusters.R
#'
#' This function adds manual celltype annotations to a seurat object
#'
#' @param seurat_obj Input seurat object
#' @param resolution Resolution to be used for clustering
#' @param annotations Named vector of cell types for each cluster
#' @param plot_directory Location to save final annotations UMAP
#' @return Seurat object
#' @export


annotate_clusters <- function(seurat_obj=NULL, 
                              resolution=NULL, 
                              annotations=NULL, 
                              plot_directory=NULL){
  
  sample_ID=seurat_obj@meta.data$sample_ID[1]
  cluster_col <- paste0("ADT_snn_res.", resolution)
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop("Couldn't find cluster column: ", cluster_col)
  }
  Idents(seurat_obj) <- cluster_col
  
  # map clusters -> labels
  lab <- annotations[as.character(Idents(seurat_obj))]
  lab[is.na(lab)] <- "Unknown"              # optional default
  names(lab) <- colnames(seurat_obj)
  
  seurat_obj$final_annotation <- lab
  Idents(seurat_obj) <- "final_annotation"
  
  if(!is.null(plot_directory)){
    DimPlot(seurat_obj,
            reduction = "umap",
            group.by  = "final_annotation",
            label     = TRUE,
            repel     = TRUE, label.box = T, pt.size=2) +
      ggtitle("Final Annotations")
    
    ggsave(file.path(plot_directory,paste0(sample_ID,"_final_annotations.png")))
  }
  
  
  
  return(seurat_obj)
  
  
  
}
