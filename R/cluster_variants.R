

cluster_variants <- function(
    seurat_obj,
    genotypes_file,
    plot_directory = NULL,
    resolution=0.2
){
  
  barcodes <- colnames(seurat_obj)
  
  sample_ID <- seurat_obj@meta.data$sample_ID[1]
  
  #Read header
  
  fmt_line <- readLines(genotypes_file, n = 1L)
  format_fields <- strsplit(gsub("^#FORMAT=", "", fmt_line), ";", fixed = TRUE)[[1]]
  cols <- colnames(data.table::fread(genotypes_file, nrows = 2L, sep = "\t"))
  
  # Generate AF matrix of priority + pass variants 
  
  row_idx <- which(fread(genotypes_file, select="filter")$filter==".") + 2
  
  variants  <- fread(cmd = sprintf("awk 'NR==%s' %s", paste(row_idx, collapse="||NR=="), genotypes_file), sep = "\t") %>%
    setNames(cols) %>%
    filter(priority_flag==1) %>%
    as.data.frame()
  
  rownames(variants) <- variants$plot_ID
  
  pull_AF <- function(s) suppressWarnings(as.numeric(unlist(lapply(str_split(s, ";"), function(x) x[which(format_fields=="AF")]))))
  
  
  af_matrix <- variants[barcodes] %>%
    mutate(across(everything(), pull_AF))
  
  af_matrix <- as.matrix(af_matrix)
  af_matrix[!is.finite(af_matrix)] <- 0
  
  # filtering
  
  min_cells_with_signal <- 5L
  nonzero_counts <- rowSums(af_matrix > 0, na.rm = TRUE)
  row_var         <- matrixStats::rowVars(af_matrix)
  keep <- (nonzero_counts >= min_cells_with_signal) & (row_var > 0)
  af_matrix <- af_matrix[keep, , drop = FALSE]
  if (nrow(af_matrix) == 0) stop("No variants passed filtering for PCA")
  
  # Add assay to seurat object
  
  assay <- CreateAssayObject(counts = af_matrix)
  assay <- SetAssayData(assay, slot = "data", new.data = af_matrix)
  seurat_obj[["variants"]] <- assay
  DefaultAssay(seurat_obj) <- "variants"
  
  # scaling
  seurat_obj <- ScaleData(seurat_obj, features = rownames(af_matrix), verbose = FALSE)
  
  # calculate reasonable number of PCs
  npcs <- max(1L, min(20L, nrow(af_matrix), ncol(af_matrix) - 1L))
  
  # run PCA
  seurat_obj <- RunPCA(
    seurat_obj,
    features = rownames(af_matrix),
    npcs = npcs,
    reduction.name = "pca.variants",
    verbose = FALSE
  )
  
  # run UMAP
  
  seurat_obj <- FindNeighbors(
    seurat_obj,
    reduction = "pca.variants",
    dims = 1:npcs,
    graph.name = "variants_snn"
  )
  
  seurat_obj <- FindClusters(
    seurat_obj,
    graph.name = "variants_snn",
    resolution = resolution,
    group.singletons = TRUE
  )
  seurat_obj$variants_clusters <- seurat_obj$seurat_clusters
  
  seurat_obj <- RunUMAP(
    seurat_obj,
    reduction = "pca.variants",
    dims = 1:npcs,
    reduction.name = "umap.variants"
  )
  
  Idents(seurat_obj) <- "variants_clusters"
  p <- DimPlot(
    seurat_obj,
    reduction = "umap.variants",
    group.by  = "variants_clusters",
    label = TRUE, repel = TRUE, label.box = TRUE, pt.size = 2
  )
  print(p)
  
  if(!is.null(plot_directory)){
    ggsave(plot=p,filename = file.path(plot_directory,paste0(sample_ID,"_variant_umap.png")))
  }
  
  invisible(seurat_obj)
  
}
