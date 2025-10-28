#' annotate_cells.R
#'
#' This function automatically annotates individual cells in a seurat object
#'
#' @param seurat_obj Input seurat object created using prepare_annotation_data()
#' @param markers Either "CellMarker2.0" (default) or named list of markers for cell types
#' @param plot_directory Output directory to save plots
#' @param CellMarker_path Path to CellMarker2.0 
#' @param celltypes Can be used to limit the set of possible celltypes to be considered
#' @param tissue Only required if using CellMarker2.0
#' @param scina_probability_threshold Only keep celltype assignments with probability above this thresholds (Default 0.5)
#' @return Seurat object
#' @export


annotate_cells <- function(
    seurat_obj=NULL,
    markers="CellMarker2.0",
    plot_directory=NULL,
    CellMarker_path=NULL,
    celltypes=NULL,
    tissue="normal peripheral blood",
    scina_probability_threshold=0.5
    
){
  
  if(!dir.exists(plot_directory)){
    dir.create(plot_directory, recursive = T)
  }
  
  sample_ID=seurat_obj@meta.data$sample_ID[1]
  
  ########################
  # Prepare Cell Markers #
  ########################
  
  message("Preparing cell markers...")
  
  if(is.list(markers)){
    cell_markers <- markers
    # only keep markers present in seurat object
    cell_markers <- lapply(cell_markers, function(v) v[v %in% Features(seurat_obj)])
  }else if(markers=="CellMarker2.0"){
    # Read in and filter CellMarker data
    if(tissue=="normal peripheral blood"){
      cell_markers <- read_xlsx(CellMarker_path) %>%
        filter(tissue_type=="Peripheral blood") %>%
        filter(cell_type=="Normal cell") %>%
        dplyr::select(cell_name, marker, Symbol) %>%
        group_by(cell_name) %>%
        reframe(markers=unique(append(marker,Symbol))) %>%
        na.omit()
    }
    # only keep markers present in seurat object
    cell_markers <- cell_markers %>%
      filter(markers %in% Features(seurat_obj))
    
    # create SCINA cell_markers list
    cell_markers <- split(cell_markers, cell_markers$cell_name) %>% lapply(deframe)
    
  }
  
  #############
  # Run SCINA #
  #############
  
  normalised_mat <- GetAssayData(seurat_obj, assay = "ADT", layer = "data")
  
  message("Running SCINA...")
  
  scina_result = SCINA(normalised_mat, cell_markers, max_iter = 100, convergence_n = 10, 
                       convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=FALSE, allow_unknown=TRUE, log_file='SCINA.log')
  
  # takes too long to run
  #plotheat.SCINA(normalised_mat, scina_result, cell_markers)
  
  # Pull annotations and probabilities
  
  scina_cell_annotations <- setNames(scina_result$cell_labels,colnames(normalised_mat))
  scina_probability_mat <- scina_result$probabilities
  scina_probability_mat <- scina_probability_mat[,names(scina_cell_annotations)]
  
  known <- scina_cell_annotations %in% rownames(scina_probability_mat) 
  
  df_known <- data.frame(
    cell_barcode = names(scina_cell_annotations)[known],
    cell_type    = scina_cell_annotations[known],
    probability  = scina_probability_mat[cbind(scina_cell_annotations[known], names(scina_cell_annotations)[known])],
    row.names    = NULL
  )
  
  if(any(!known)){
    # unknown assignments: keep them with NA prob
    df_unknown <- data.frame(
      cell_barcode = names(scina_cell_annotations)[!known],
      cell_type    = scina_cell_annotations[!known],
      probability  = NA,
      row.names    = NULL
    )
    df_known <- rbind(df_known, df_unknown)
  }
  
  
  final_annotations <- df_known %>%
    # filter annotations
    mutate(cell_type=ifelse(probability<scina_probability_threshold, "unknown",cell_type)) %>%
    dplyr::select(cell_barcode,cell_type) %>%
    deframe()
  
  #################################
  # Visualise Annotations on UMAP #
  #################################
  
  message(paste0("Generating ",file.path(plot_directory,paste0(sample_ID,"_SCINA_cell_annotations_UMAP.png"))))
  # add SCINA annotations to seurat object
  common <- intersect(colnames(seurat_obj), names(final_annotations))
  seurat_obj$SCINA_label <- NA
  seurat_obj$SCINA_label[common] <- final_annotations[common]
  
  DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by  = "SCINA_label",
    label     = FALSE,
    raster    = FALSE,
    pt.size   = 2
  ) + ggplot2::labs(color = "SCINA")
  
  ggsave(file.path(plot_directory,paste0(sample_ID,"_SCINA_cell_annotations_UMAP.png")))
  
  resolutions <- colnames(seurat_obj@meta.data)[str_detect(colnames(seurat_obj@meta.data),"ADT_snn_res")] %>%
    gsub(pattern="ADT_snn_res.",replacement="")
  
  ############################################
  # Investigate Cell Annotations Per Cluster #
  ############################################
  
  for(res in resolutions){
    
    res_col <- paste0("ADT_snn_res.", res)
    
    message(paste0("Generating ",file.path(plot_directory,paste0(sample_ID,"_seurat_clusters_res_",res,".png"))))
    
    # seurat clusters
    DimPlot(seurat_obj,
            reduction = "umap",
            group.by  = res_col,
            label     = TRUE,
            repel     = TRUE, label.box = T) +
      ggtitle(paste("Seurat Clusters at Resolution: ",res))
    
    ggsave(file.path(plot_directory,paste0(sample_ID,"_seurat_clusters_res_",res,".png")))
    
    message(paste0("Generating ",file.path(plot_directory,paste0(sample_ID,"_marker_expression_heatmap_res_",res,".png"))))
    
    # heatmap marker expression
    DoHeatmap(
      seurat_obj,
      features = Features(seurat_obj),
      group.by = res_col,
      raster   = TRUE,
      slot     = "data"
    )
    
    ggsave(file.path(plot_directory,paste0(sample_ID,"_marker_expression_heatmap_res_",res,".png")))
    
    
    clust_vec <- seurat_obj@meta.data[[res_col]]
    scina_vec <- seurat_obj$SCINA_label
    valid <- !is.na(clust_vec) & !is.na(scina_vec)
    
    # proportions of cell annotation per cluster
    tab  <- table(Cluster = clust_vec[valid], SCINA = scina_vec[valid])
    prop <- prop.table(tab, margin = 1)
    
    ord_rows <- order(rowSums(tab), decreasing = TRUE)
    ord_cols <- order(colSums(tab), decreasing = TRUE)
    
    message(paste0("Generating ",file.path(plot_directory,paste0(sample_ID,"_SCINA_cluster_proportions_res_",res,".png"))))
    
    pheatmap::pheatmap(
      prop[ord_rows, ord_cols, drop = FALSE],
      cluster_rows = FALSE, cluster_cols = FALSE,
      display_numbers = TRUE, number_format = "%.2f",
      main = paste0("SCINA proportions per cluster (Resolution ", res, ")"),
      filename = file.path(plot_directory,paste0(sample_ID,"_SCINA_cluster_proportions_res_",res,".png"))
    )
    
    # Build a data frame with coords + labels
    df <- Embeddings(seurat_obj, "umap")[valid, 1:2] |>
      as.data.frame() |>
      `colnames<-`(c("UMAP_1","UMAP_2"))
    df$cluster     <- clust_vec[valid]
    df$SCINA_label <- scina_vec[valid]
    
    if(length(unique(df$cluster)) <=9 ){
      pastel_cols <- brewer.pal(max(3, length(unique(df$cluster))), "Pastel1")
      
      message(paste0("Generating ",file.path(plot_directory,paste0(sample_ID,"_SCINA_annotations_per_cluster_res_",res,".png"))))
      
      ggplot(df, aes(UMAP_1, UMAP_2)) +
        # Shadow layer: opaque pastel clusters
        geom_point(
          aes(color = cluster),
          size = 10,
          alpha = 0.9,                  # fully opaque, no alpha stacking
          show.legend = TRUE
        ) +
        scale_color_manual(values = pastel_cols,
                           guide = guide_legend(
                             override.aes = list(alpha = 1, size = 10, shape = 16)
                           )) +
        # Foreground layer: SCINA colours
        ggnewscale::new_scale_color() +
        geom_point(
          aes(color = SCINA_label),
          size = 0.9,
          alpha = 0.9
        ) +
        coord_equal() +
        theme_classic() +
        labs(color = "SCINA") +
        scale_color_discrete(
          name  = "SCINA",
          guide = guide_legend(
            override.aes = list(alpha = 1, size = 2, shape = 16)
          )
        )
      
      ggsave(file.path(plot_directory,paste0(sample_ID,"_SCINA_annotations_per_cluster_res_",res,".png")))
    }
    
  }
  
  return(seurat_obj)
  
  
}