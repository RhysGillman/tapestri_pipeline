#' annotate_clusters.R
#'
#' This function adds manual variant annotations to a seurat object
#'
#' @param seurat_obj Input seurate object
#' @param variant Variant to be annotated
#' @param variant_file Path to the *_all.tsv file containing cell-specific variant info
#' @return Seurat object
#' @export

annotate_variants <- function(
    seurat_obj,
    variant,
    variant_file,
    reduction = "umap",
    plot_directory = NULL,
){
  stopifnot(!is.null(seurat_obj), !is.null(variant), !is.null(variant_file))
  
  sample_ID=seurat_obj@meta.data$sample_ID[1]
  
  #Read header
  fmt_line <- readLines(variant_file, n = 1L)
  format_fields <- strsplit(gsub("^#FORMAT=", "", fmt_line), ";", fixed = TRUE)[[1]]
  cols <- colnames(data.table::fread(variant_file, nrows = 2L, sep = "\t"))
  
  # find variant row
  variant_row <- data.table::fread(
    cmd = sprintf("grep -m1 -F '%s' '%s'", variant, variant_file),
    sep = "\t", header = FALSE
  )
  if (nrow(variant_row) == 0) stop("Variant not found: ", variant)
  
  data.table::setnames(variant_row, cols)
  
  # sanity check
  if (variant_row[["variant"]][1] != variant) {
    stop("Found row does not exactly match requested variant: ", cell_row[[variant_col]][1])
  }
  
  cell_variants <- t(variant_row[,15:ncol(variant_row)]) %>%
    as.data.frame() %>%
    rownames_to_column("barcode") %>%
    deframe()
  
  # pull out data based on #FORMAT
  data <- data.table::transpose(as.data.frame(strsplit(cell_variants, ";", fixed = T)),keep.names = "barcode")
  colnames(data) <- c("barcode", format_fields)
  data <- data %>%
    mutate(across(-barcode,as.numeric))
  
  
  # Prepare safe column name suffix based on variant (for uniqueness)
  safe_var <- gsub("[^A-Za-z0-9_]+", "_", variant)
  
  seurat_cells <- colnames(seurat_obj)
  
  add_col <- function(seurat_obj, vec, base_name) {
    colname <- paste0(base_name, "_", safe_var)
    tmp <- rep(NA, length(seurat_cells))
    names(tmp) <- seurat_cells
    # subset on intersection
    keep <- intersect(seurat_cells, data$barcode)
    if (length(keep)) {
      m <- match(keep, data$barcode)
      tmp[keep] <- vec[m]
    }
    seurat_obj[[colname]] <- tmp
    return(seurat_obj)
  }
  
  for(value in format_fields){
    seurat_obj <- add_col(seurat_obj,as.vector(data[[value]]), value)
  }
  
  plot_features <- colnames(seurat_obj@meta.data)[which(str_detect(colnames(seurat_obj@meta.data), safe_var))]
  
  for(feature in plot_features){
    type <- str_extract(feature, "^([^_]+)_", group = 1)
    if (type=="NGT") {
      seurat_obj[[feature]] <- factor(
        as.vector(seurat_obj[[feature]])[[1]],
        levels = c(0,1,2,3),
        labels = c("WT","HET","HOM","Unavailable")
      )
      DimPlot(
        seurat_obj, reduction = reduction, group.by = feature,
        cols = c("WT"="black","HET"="blue","HOM"="red","Unavailable"="grey"),
        pt.size   = 2, alpha = 0.8
      ) + ggplot2::ggtitle(paste0(type,": ", variant))
      
      ggsave(file.path(plot_directory,paste0(sample_ID,"_",feature,".png")))
    }else{
    
      FeaturePlot(seurat_obj, features = feature, reduction = reduction, order = TRUE, pt.size = 2) +
        ggplot2::ggtitle(paste0(type,": ", variant))
      ggsave(file.path(plot_directory,paste0(sample_ID,"_",feature,".png")))
    }
    
  }
  
  return(seurat_obj)
}