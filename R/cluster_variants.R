

cluster_variants <- function(
    seurat_obj=NULL,
    genotypes_file,
    summary_file,
    plot_directory = NULL,
    resolution=0.2,
    n_features=100,
    sample_ID=NULL
    
){
  
  if(!is.null(seurat_obj)){
    barcodes <- colnames(seurat_obj)
    sample_ID <- seurat_obj@meta.data$sample_ID[1]
  }else{
    barcodes <- .detect_barcodes(genotypes_file)
    if(is.null(sample_ID)){
      stop("Please provide sample_ID")
    }
  }
  
  # Find features
  
  all_variants <- fread(summary_file) %>%
    mutate(rownumber=row_number()) %>%
    filter(filter==".", data_proportion_total >= 0.80) %>%
    # get NGT variance (max = 0.5, deprioritises too low OR too high frequency)
    mutate(
      p_alt   = ifelse(data_cnt_total > 0, alt_cnt_total / data_cnt_total, NA_real_),
      var_ngt = p_alt * (1 - p_alt)                             # variability score
    ) %>%
    arrange(desc(var_ngt))
  rare_variants <- all_variants %>%
    filter(
      p_alt < 0.02,
      mean_GQ_total >= 30
    ) %>%
    arrange(desc(alt_cnt_total)) %>%
    head(n_features*0.1)
  selection <- bind_rows(
    all_variants %>% head(n_features*0.9), rare_variants
  ) %>%
    pull(rownumber) %>%
    sort()
  
  row_idx <- selection + 2
  #Read header
  
  fmt_line <- readLines(genotypes_file, n = 1L)
  format_fields <- strsplit(gsub("^#FORMAT=", "", fmt_line), ";", fixed = TRUE)[[1]]
  head_file <- data.table::fread(genotypes_file, nrows = 20L, sep = "\t")
  cols <- colnames(head_file)
  coltypes <- sapply(head_file, class)
  
  
  # Generate AF matrix of pass variants 
  
  #row_idx <- which(fread(genotypes_file, select="filter")$filter==".") + 2
  pb   <- txtProgressBar(min = 0, max = length(row_idx), style = 3)
  variants  <- tryCatch({
    fread(cmd = sprintf("awk 'NR==%s' %s", paste(row_idx, collapse="||NR=="), genotypes_file), sep = "\t") %>%
    setNames(cols) %>%
    as.data.frame()}, error = function(e) {
      message("Warning: Genotype file read failed: ", conditionMessage(e), " - using fallback.")
      rbindlist(lapply(sort(unique(row_idx)), function(r) {
        i <- which(row_idx==r)
        setTxtProgressBar(pb, i)
        fread(genotypes_file, skip = r - 1L, nrows = 1L, header = FALSE, sep = "\t")
      }), use.names = FALSE) %>%
        setNames(cols) %>%
        as.data.frame()
    }
    )
  
  rownames(variants) <- variants$plot_ID
  
  pull_AF <- function(s) suppressWarnings(as.numeric(unlist(lapply(str_split(s, ";"), function(x) x[which(format_fields=="AF")]))))
  
  
  af_matrix <- variants[barcodes] %>%
    mutate(across(everything(), pull_AF))
  
  af_matrix <- as.matrix(af_matrix)
  #af_matrix[!is.finite(af_matrix)] <- 0
  
  # filtering
  
  min_cells_with_signal <- 5L
  nonzero_counts <- rowSums(af_matrix > 0, na.rm = TRUE)
  row_var         <- matrixStats::rowVars(af_matrix, na.rm = TRUE)
  keep <- (nonzero_counts >= min_cells_with_signal) & (row_var > 0)
  af_matrix <- af_matrix[keep, , drop = FALSE]
  if (nrow(af_matrix) == 0) stop("No variants passed filtering for PCA")
  
  # Add assay to seurat object
  
  assay_counts <- af_matrix
  assay_counts[!is.finite(assay_counts)] <- 0
  assay <- Seurat::CreateAssayObject(counts = assay_counts)
  assay <- Seurat::SetAssayData(assay, layer = "data", new.data = assay_counts)
  seurat_obj[["variants"]] <- assay
  Seurat::DefaultAssay(seurat_obj) <- "variants"
  
  
  
  cosine_dist_na <- function(A, min_overlap = 10L) {
    p <- nrow(A); n <- ncol(A)
    D <- matrix(0, n, n)
    for (i in seq_len(n)) {
      for (j in i:n) {
        ok <- is.finite(A[, i]) & is.finite(A[, j])
        if (sum(ok) < min_overlap) {
          d <- 1  # maximal distance if too little shared data
        } else {
          xi <- A[ok, i]; xj <- A[ok, j]
          num <- sum(xi * xj)
          den <- sqrt(sum(xi^2)) * sqrt(sum(xj^2))
          sim <- if (den > 0) num / den else 0
          d <- 1 - sim
        }
        D[i, j] <- D[j, i] <- d
      }
    }
    stats::as.dist(D)
  }
  
  D <- cosine_dist_na(af_matrix, min_overlap = 10L)
  
  # --- UMAP on distance matrix ---------------------------------------------
  set.seed(12345)
  umap_emb <- uwot::umap(
    D,
    n_neighbors = 30, min_dist = 0.3, metric = "euclidean",
    n_components = 2, verbose = FALSE
  )
  rownames(umap_emb) <- colnames(af_matrix)
  colnames(umap_emb) <- c("UMAP_1", "UMAP_2")
  
  minPts <- max(20L, round(ncol(af_matrix) * 0.01))
  db <- dbscan::dbscan(umap_emb, eps = 0.8, minPts = minPts)
  clusters <- setNames(db$cluster, rownames(umap_emb))
  cluster_labels <- ifelse(clusters == 0, "Noise", paste0("C", clusters))
  common <- intersect(colnames(seurat_obj), names(clusters))
  seurat_obj <- subset(seurat_obj, cells = common)
  umap_emb   <- umap_emb[common, , drop = FALSE]
  
  seurat_obj[["umap.variants"]] <- Seurat::CreateDimReducObject(
    embeddings = umap_emb,
    key = "umapVAF_",
    assay = "variants"
  )
  seurat_obj$variants_clusters <- clusters[colnames(seurat_obj)]  ### CHANGED: use DBSCAN labels
  Seurat::Idents(seurat_obj) <- "variants_clusters"
  
  p <- Seurat::DimPlot(
    seurat_obj, reduction = "umap.variants",
    group.by = "variants_clusters",
    label = TRUE, repel = TRUE, label.box = TRUE, pt.size = 2
  )
  print(p)
  
  if (!is.null(plot_directory)) {
    ggplot2::ggsave(plot = p, filename = file.path(plot_directory, paste0(sample_ID, "_variant_umap.png")))
  }
  
  invisible(seurat_obj)
  
}
