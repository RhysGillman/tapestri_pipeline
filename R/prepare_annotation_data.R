#' prepare_annotation_data.R
#'
#' This function loads annotation data from a tapestri .h5 file and generates a seurat object ready for annotation
#'
#' @param input_file Path to input .dna+protein.h5 file
#' @param sample_ID Sample ID, or NULL to derive from .h5 filname (default)
#' @param res Vector of resolutions to to be used my FindClusters
#' @param save_directory Output directory to save seurat object as .RDS, or NULL to prevent saving (default)
#' @return Seurat object
#' @export


prepare_annotation_data <- function(
  input_file=NULL,
  sample_ID=NULL,
  cluster_resolutions=c(0.2, 0.35, 0.5, 0.8, 1, 1.2),
  save_directory=NULL,
  qc_plot_directory=NULL
){
  if(is.null(sample_ID)){
    sample_ID = gsub(".dna[+]protein.h5", "", sort(basename(input_file)))
  }
  message(paste0("Preparing data for ",sample_ID," from: ", input_file))
  
  #############################
  # Reading in protein counts #
  #############################
  
  protein_count_matrix <- h5read(input_file, "/assays/protein_read_counts/layers/read_counts") 
  rownames(protein_count_matrix) <-  as.character(h5read(input_file, "/assays/protein_read_counts/ca/id"))
  colnames(protein_count_matrix) <- as.character(h5read(input_file, "/assays/protein_read_counts/ra/barcode"))
  protein_count_matrix <- Matrix(protein_count_matrix, sparse=T)
  
  obj <- CreateSeuratObject(counts = protein_count_matrix, assay = "ADT", project = "Tapestri_Protein")
  DefaultAssay(obj) <- "ADT"
  
  #################
  # Normalisation #
  #################
  
  # CLR-normalize across cells
  obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2)
  # scaled data
  obj <- ScaleData(obj, assay = "ADT", layer = "data")  # writes layer "scale.data"
  
  # Drop any zero-variance or non-finite features
  mat <- GetAssayData(obj, assay = "ADT", layer = "scale.data")
  keep <- matrixStats::rowVars(as.matrix(mat)) > 0 & Matrix::rowSums(is.finite(mat)) == ncol(mat)
  if (any(!keep)) {
    message(paste0("Warning: The following features are being removed due to having no variance or non-finite values\n",paste(names(which(keep==F)),collapse = "\n")))
    obj <- subset(obj, features = rownames(mat)[keep])
    mat <- mat[keep, , drop = FALSE]
  }
  
  ###########
  # Run PCA #
  ###########
  # calculate reasonable number of PCs
  npcs <- min(20L, nrow(mat), ncol(mat) - 1L)
  
  obj <- RunPCA(
    object  = obj,
    assay   = "ADT",
    layer   = "scale.data",
    features = rownames(mat),
    npcs    = npcs,
    approx  = FALSE,
    verbose = FALSE
  )
  
  ##############
  # Elbow Plot #
  ##############
  # Take only PCs which have a change in variation is > 0.1% from the next -- OTHERWISE, can just use 7
  elb_plot <- ElbowPlot(obj, ndims=30)$data
  pct <- elb_plot$stdev
  max_pca_dim <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  if(!is.null(qc_plot_directory)){
    elb_plot <- ggplot() + 
      geom_point(data = elb_plot %>% filter(dims <= max_pca_dim), mapping = aes(dims, stdev), size = 5, color = 'green') +
      geom_point(data = elb_plot %>% filter(dims > max_pca_dim), mapping = aes(dims, stdev), size = 5, color = 'red') +
      ggtitle(sample_ID) +
      theme_bw()
    ggsave(file.path(qc_plot_directory,sample_ID,"pca_elbow_plot.png"),plot = elb_plot)
  }
  
  ########
  # UMAP #
  ########
  
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:max_pca_dim, min.dist = 0.3)
  obj <- FindNeighbors(obj, reduction = "pca", dims = 1:max_pca_dim,k.param=20)
  
  for (res in cluster_resolutions) {
    obj <- FindClusters(obj, resolution = res,random.seed=888)
  }
  
  obj@meta.data$sample_ID <- sample_ID
  obj@meta.data <- cbind(obj@meta.data, as.data.frame(obj@reductions$pca@cell.embeddings)[rownames(obj@meta.data),c(1:npcs)])
  obj@meta.data <- cbind(obj@meta.data, as.data.frame(obj@reductions$umap@cell.embeddings)[rownames(obj@meta.data),c(1:2)])
  
  if(!is.null(save_directory)){
    saveRDS(obj,file.path(save_directory,sample_ID,"_seurat_obj.rds"))
  }
  
  return(obj)
}