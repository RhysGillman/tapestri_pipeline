#' annotate_mutants.R
#'
#' This function adds manual variant annotations to a seurat object
#'
#' @param seurat_obj Input seurate object
#' @param variant Variant/s to be annotated
#' @param genotypes_file Path to the *_genotypes.tsv file containing cell-specific variant info
#' @return Seurat object
#' @export

annotate_mutants <- function(
    seurat_obj,
    variant,
    name=NULL,
    combine_strategy = "intersect",
    genotypes_file,
    reduction = "umap",
    plot_directory = NULL
){
  stopifnot(!is.null(seurat_obj), !is.null(variant), !is.null(genotypes_file))
  
  if(is.null(name)){
    name <- "mutant"
  }
  
  barcodes <- colnames(seurat_obj)
  
  sample_ID=seurat_obj@meta.data$sample_ID[1]
  
  #Read header
  fmt_line <- readLines(genotypes_file, n = 1L)
  format_fields <- strsplit(gsub("^#FORMAT=", "", fmt_line), ";", fixed = TRUE)[[1]]
  cols <- colnames(data.table::fread(genotypes_file, nrows = 2L, sep = "\t"))
  
  variant_genotypes <- data.frame(barcode=barcodes)
  
  for(v in variant){
    
    message("Searching genotypes file for variant: ", v)
    
    # grep for variant row
    variant_row <- data.table::fread(
      cmd = sprintf("grep -Fw '%s' '%s'", v, genotypes_file),
      sep = "\t", header = FALSE
    ) %>% as.data.frame()
    
    if (nrow(variant_row) > 1) stop("Error: Variant search: ", v, "returned multiple rows, make sure variant ID is unique")
    
    if (nrow(variant_row) == 0) stop("Variant not found: ", v)
    
    data.table::setnames(variant_row, cols)
    
    cell_variants <- t(variant_row[,barcodes]) %>%
      as.data.frame() %>%
      rownames_to_column("barcode") %>%
      deframe()
    
    # pull out data based on #FORMAT
    data <- data.table::transpose(as.data.frame(strsplit(cell_variants, ";", fixed = T)),keep.names = "barcode")
    colnames(data) <- c("barcode", format_fields)
    data <- suppressWarnings(data %>%
      mutate(across(-barcode,as.numeric)))
    
    # Prepare safe column name suffix based on variant (for uniqueness)
    safe_var <- gsub("[^A-Za-z0-9_]+", "_", v)
    
    add_col <- function(seurat_obj, vec, base_name) {
      colname <- paste0(base_name, "_", safe_var)
      tmp <- rep(NA, length(barcodes))
      names(tmp) <- barcodes
      # subset on intersection
      keep <- intersect(barcodes, data$barcode)
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
    
    if(!is.null(plot_directory)){
      plot_features <- colnames(seurat_obj@meta.data)[which(str_detect(colnames(seurat_obj@meta.data), safe_var))]
      
      for(feature in plot_features){
        type <- str_extract(feature, "^([^_]+)_", group = 1)
        if (type=="NGT") {
          vals <- seurat_obj[[feature, drop = TRUE]]
          seurat_obj[[feature]] <- factor(
            vals,
            levels = c(0,1,2,3),
            labels = c("WT","HET","HOM","Unavailable")
          )
          p <- DimPlot(
            seurat_obj, reduction = reduction, group.by = feature,
            cols = c("WT"="black","HET"="blue","HOM"="red","Unavailable"="grey"),
            pt.size   = 2, alpha = 0.8
          ) + ggplot2::ggtitle(paste0(type,": ", safe_var))
          
          ggsave(file.path(plot_directory,paste0(sample_ID,"_",feature,".png")), plot = p)
        }else{
          
          p <- FeaturePlot(seurat_obj, features = feature, reduction = reduction, order = TRUE, pt.size = 2) +
            ggplot2::ggtitle(paste0(type,": ", safe_var))
          ggsave(file.path(plot_directory,paste0(sample_ID,"_",feature,".png")), plot = p)
        }
        
      }
    }
    
    ngt <- data %>% dplyr::select(barcode,NGT) %>% setnames(c("barcode",safe_var))
    
    variant_genotypes <- variant_genotypes %>%
      left_join(ngt, by = "barcode")
    
  }
  
  variant_genotypes <- variant_genotypes %>% column_to_rownames("barcode") %>% as.matrix
  # replace 0,1,2,3 as binary 0,1 matrix
  variant_genotypes[variant_genotypes%in%c(1,2)] <- 1
  variant_genotypes[variant_genotypes!=1] <- 0
  
  if(combine_strategy=="union"){
    mutant_cells <- names(which(rowSums(variant_genotypes) > 0))
  }else if(combine_strategy=="intersect"){
    mutant_cells <- names(which(rowSums(variant_genotypes) == ncol(variant_genotypes)))
  }
  
  seurat_obj[[name]] <- ifelse(colnames(seurat_obj) %in% mutant_cells, "mutant", "wildtype")
  
  return(seurat_obj)
}