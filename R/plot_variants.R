#' annotate_clusters.R
#'
#' This function adds manual variant annotations to a seurat object
#'
#' @param seurat_object Input seurate object
#' @param variant Variant to be annotated
#' @param genotypes_file Path to the *_genotypes.tsv file containing cell-specific variant info
#' @return Seurat object
#' @export

plot_variants <- function(
    seurat_object=NULL,
    file_prefix,
    search=NULL,
    n_variants=20,
    search_col="plot_ID",
    result_dir,
    plot_directory,
    plot_prefix,
    # filters
    damaging=T,
    priority_only=T,
    p_max=NULL,
    z_min=2,
    min_varcount_total  = 20,
    max_alt_portion_total = 0.5,
    topn_col="padj",
    font_size_multiplyer=1
){
  #########
  # Files #
  #########
  
  variant_counts_file <- file.path(result_dir,paste0(file_prefix,"_variant_counts_per_celltype.tsv"))
  variant_counts_cols <- colnames(data.table::fread(variant_counts_file, nrows = 2L, sep = "\t"))
  genotype_file <- file.path(result_dir,paste0(file_prefix,"_genotypes.tsv"))
  genotype_fmt_line <- readLines(genotype_file, n = 1L)
  genotype_format_fields <- strsplit(gsub("^#FORMAT=", "", genotype_fmt_line), ";", fixed = TRUE)[[1]]
  genotype_cols <- colnames(data.table::fread(genotype_file, nrows = 2L, sep = "\t"))
  
  all_files <- list.files(
    path = result_dir
  )
  
  cell_type_files <- grep(paste0("^",file_prefix,"_(?!summary|vep_in|genotypes|variant_counts_per_celltype).*.tsv"), all_files, value = TRUE, perl = TRUE)
  
  if(priority_only){
    summary_file <- file.path(result_dir,paste0(file_prefix,"_summary_pass_only_priority.tsv"))
    cell_type_files <- cell_type_files[str_detect(cell_type_files, "pass_only_priority.tsv$")]
  }else{
    summary_file <- file.path(result_dir,paste0(file_prefix,"_summary_pass_only.tsv"))
    cell_type_files <- cell_type_files[str_detect(cell_type_files, "pass_only.tsv$")]
  }
  summary_cols <- colnames(data.table::fread(summary_file, nrows = 2L, sep = "\t"))
  
  cell_types <- unique(str_extract(cell_type_files, pattern = paste0("^",file_prefix,"_(.*)_pass.*.tsv"), group = 1))
  barcodes <- .detect_barcodes(genotype_file,genotype_format_fields)
  
  ########################
  # Variants of Interest #
  ########################
  
  # If not searching, just take top n results #
  
  if(is.null(search)){
    variant_rows <- foreach(file=file.path(result_dir,cell_type_files), .combine = "bind_rows") %do%{
      c <- str_extract(file, pattern = paste0(file_prefix,"_(.*)_pass.*.tsv"), group = 1)
      v <- fread(file) %>%
        mutate(cell_type=c)
      if (!is.null(p_max)) v <- filter(v, padj < p_max)
      if (!is.null(z_min)) v <- filter(v, Z_score > z_min)
      v
    }
  
    if(topn_col %in% c("OR","Z_score")){
      variant_rows <- variant_rows %>%
        arrange(desc(!!sym(topn_col))) %>%
        head(n_variants)
    }else if(topn_col %in% c("padj")){
      variant_rows <- variant_rows %>%
        arrange(!!sym(topn_col)) %>%
        head(n_variants)
    }
    
    variant_rows <- variant_rows %>%
      # apply filters
      filter(alt_cnt_total >= min_varcount_total) %>%
      filter(alt_proportion_total <= max_alt_portion_total)
    
    
    variants <- variant_rows$plot_ID
    names(variants) <- variant_rows$variant
    
  }else{
    
    # Let user search for variants by defining a search term and a column #
    
    if(!search_col %in% summary_cols){
      stop("search_col not found in variant data")
    }
    
    variant_rows <- foreach(s=search, .combine = "bind_rows") %do% {
      
      fread(summary_file) %>%
        filter(str_detect(.data[[search_col]], regex(s)))
    } %>%
      # apply filters
      filter(alt_cnt_total >= min_varcount_total) %>%
      filter(alt_proportion_total <= max_alt_portion_total)
    
    if (isTRUE(damaging)) variant_rows <- filter(variant_rows, str_detect(SIFT,"Deleterious") | str_detect(PolyPhen,"damaging") )
    
    
    if(nrow(variant_rows)==0){
      stop("ERROR: Search and filters returned no variants")
    }
    
    # apply cell-type result filters
    if (!is.null(p_max) || !is.null(z_min)) {
      all_results <- foreach(file=file.path(result_dir,cell_type_files), .combine = "bind_rows") %do%{
        c <- str_extract(file, pattern = paste0(file_prefix,"_(.*)_pass.*.tsv"), group = 1)
        fread(file) %>%
          mutate(cell_type=c)
      }
      
      sig_variants <- all_results
      if (!is.null(p_max)) sig_variants <- filter(sig_variants, padj < p_max)
      if (!is.null(z_min)) sig_variants <- filter(sig_variants, Z_score > z_min)
      sig_variants <- pull(sig_variants, variant)
      
      variant_rows <- variant_rows %>%
        filter(variant %in% sig_variants)
    }
    
    
    variants <- variant_rows$plot_ID
    names(variants) <- variant_rows$variant
    
    message(paste0("Search returned ", length(variants), " variants"))
    
    if(length(variants)>n_variants){
      stop("ERROR: Search returned too many variants. Use stricter settings or increase n_variants")
    }
  }
  
  ################################
  # Variant Counts Per Cell Type #
  ################################
  
  message("Retrieving variant counts per cell-type")
  
  variant_counts <- foreach(variant=variants, .combine = "bind_rows") %do% {
    suppressWarnings(variant_row <- data.table::fread(
      cmd = sprintf("grep -Fw '%s' '%s'", variant, variant_counts_file),
      sep = "\t"
    ) %>% 
      as.data.frame() %>%
      setNames(variant_counts_cols) %>%
      mutate(chromosome=as.character(chromosome),
             MAX_AF=as.double(MAX_AF)))
  } %>%
    # deal with special case duplicate plot_IDs
    mutate(plot_ID=ifelse(plot_ID %in% plot_ID[which(duplicated(plot_ID))], paste0(plot_ID,"(",variant,")"), plot_ID)) %>%
    dplyr::select(`variant`,plot_ID, ends_with("count")) %>%
    pivot_longer(-c(`variant`,plot_ID),names_to="col", values_to="count") %>%
    mutate(cell_type=gsub("_alt_count|_data_count","",col),
           feature=ifelse(str_detect(col,"alt"),"alt_count","data_count")
    )%>%
    dplyr::select(`variant`,plot_ID,cell_type,feature,count) %>%
    pivot_wider(names_from=feature, values_from=count) %>%
    mutate(variant = factor(variant, levels = names(variants))) %>%
    arrange(variant) %>%
    mutate(plot_ID=factor(plot_ID, levels=unique(plot_ID)),
           alt_frequency=alt_count/data_count)
  
  
  ggplot(variant_counts, aes(x=plot_ID, fill=cell_type, y=alt_frequency)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x="Variant ID", y="Mutant Allele Proportion") +
    guides(fill=guide_legend(title="Cell Type")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, hjust=1))
  
  ggsave(file.path(plot_directory,paste0(plot_prefix,"_barplot.png")))
  
  message(paste0("Saved plot showing variant count per cell-type to "), file.path(plot_directory,paste0(plot_prefix,"_barplot.png")))
  
  ##############
  # Upset Plot #
  ##############
  
  if(length(variants) <= 3){
    message("Skipping Upset plot due to too few variants returned")
  }else{
    is_mutated = function(s) {
      components <- unlist(str_split(s, ";"))
      NGT <- components[which(genotype_format_fields=="NGT")]
      return(suppressWarnings(as.numeric(NGT)) %in% c(1,2))
    }
    
    
    sel_idx <- match(barcodes, genotype_cols)
    
    upset_data <- foreach(v = names(variants)) %do% {
      dt <- data.table::fread(
        cmd       = sprintf("grep -Fwm1 -- %s %s", shQuote(v), shQuote(genotype_file)),
        sep       = "\t",
        header    = FALSE,
        select    = sel_idx,
        col.names = barcodes,
        showProgress = FALSE
      )
      
      if (nrow(dt) == 0L) return(character(0))
      mutated <- vapply(dt, is_mutated, logical(1))
      names(mutated)[mutated]
    }
    
    
    
    
    upset_data <- foreach(v=names(variants)) %do% {
      variant_row <- data.table::fread(
        cmd = sprintf("grep -Fw '%s' '%s'", v, genotype_file),
        sep = "\t"
      ) %>% 
        as.data.frame() %>%
        setNames(genotype_cols) %>%
        filter(variant==v)
      variant_row <- variant_row[barcodes]
      mutated <- vapply(variant_row, is_mutated, logical(1))
      names(which(mutated))
      
    }
    names(upset_data) <- variants
    
    bc <- sort(unique(unlist(upset_data)))
    upset_data <- as.data.frame(
      setNames(lapply(upset_data, function(v) bc %in% v), names(upset_data))
    )
    rownames(upset_data) <- bc
    
    sets <- colnames(upset_data)
    
    upset_data[sets] <- lapply(upset_data[sets], function(x) as.logical(as.integer(x)))
    
    p <- upset(
      n_intersections=10,
      upset_data,
      intersect = colnames(upset_data),
      base_annotations = list(
        "Intersection size" = ComplexUpset::intersection_size(
          text = list(vjust = -0.5),
          text_colors = c(on_background = "black", on_bar = "black")
        )
      ),
      set_sizes = ComplexUpset::upset_set_size()
    )
    
    ggsave(file.path(plot_directory, paste0(plot_prefix, "_upset.png")),
           p, width = 20, height = 20, units = "cm", dpi = 300, bg = "white")
    
    
    #png(filename=file.path(plot_directory,paste0(plot_prefix,"_upset.png")),res = 300, width = 20, height = 20, units = "cm", type = "cairo")
    #on.exit(grDevices::dev.off(), add = TRUE)
    #grid::grid.newpage()
    #UpSetR::upset(UpSetR::fromList(upset_data), order.by = "freq",nsets = length(variants))
    #dev.off()
    message(paste0("Saved variant upset plot to "), file.path(plot_directory,paste0(plot_prefix,"_upset.png")))
  }
    
  ##############
  # UMAP Plots #
  ##############
  
  font_scale <- seq(from=2,to=1,length.out=20)
  panel_multiplyer <- font_scale[length(variants)]
  
  umap_plot_title_size=5*font_size_multiplyer*panel_multiplyer
  umap_legend_title_size=8*font_size_multiplyer*panel_multiplyer
  umap_legend_label_size=8*font_size_multiplyer*panel_multiplyer
  umap_axis_title_size=5*font_size_multiplyer*panel_multiplyer
  umap_axis_lables_size=4*font_size_multiplyer*panel_multiplyer
  #umap_ref_title_size=8*font_size_multiplyer
  umap_ref_label_size=4*font_size_multiplyer
  umap_ref_axis_title_size=10*font_size_multiplyer
  umap_ref_axis_lables_size=8*font_size_multiplyer
  
  
  if(is.null(seurat_object)){
    paste0("Skipping UMAP plots because no seurat object supplied")
  }else{
    umap_plots <- foreach(v=names(variants)) %do% {
      
      seurat_anno_data <- data.table::fread(
        cmd = sprintf("grep -Fw '%s' '%s'", v, genotype_file),
        sep = "\t"
      ) %>% 
        as.data.frame() %>%
        setNames(genotype_cols) %>%
        filter(variant==v)
      seurat_anno_data <- seurat_anno_data[c("plot_ID",barcodes)] %>%
        column_to_rownames("plot_ID") %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("barcode") %>%
        deframe()
      
      seurat_anno_data <- data.table::transpose(as.data.frame(strsplit(seurat_anno_data, ";", fixed = T)),keep.names = "barcode")
      colnames(seurat_anno_data) <- c("barcode", genotype_format_fields)
      seurat_anno_data <- suppressWarnings(seurat_anno_data %>%
        mutate(across(-barcode,as.numeric)))
      
      
      # Prepare safe column name suffix based on variant (for uniqueness)
      safe_var <- gsub("[^A-Za-z0-9_]+", "_", v)
      
      seurat_cells <- colnames(seurat_object)
      
      add_col <- function(seurat_object, vec, base_name) {
        colname <- paste0(base_name, "_", safe_var)
        tmp <- rep(NA, length(seurat_cells))
        names(tmp) <- seurat_cells
        # subset on intersection
        keep <- intersect(seurat_cells, seurat_anno_data$barcode)
        if (length(keep)) {
          m <- match(keep, seurat_anno_data$barcode)
          tmp[keep] <- vec[m]
        }
        seurat_object[[colname]] <- tmp
        return(seurat_object)
      }
      
      for(value in genotype_format_fields){
        seurat_object <- add_col(seurat_object,as.vector(seurat_anno_data[[value]]), value)
      }
      
      plot_features <- colnames(seurat_object@meta.data)[which(str_detect(colnames(seurat_object@meta.data), safe_var))]
      
      plots <- foreach(feature=plot_features) %do% {
        type <- str_extract(feature, "^([^_]+)_", group = 1)
        if (type=="NGT") {
          
          seurat_object[[feature]] <- factor(
            as.vector(seurat_object[[feature]])[[1]],
            levels = c(3,0,1,2),
            labels = c("Unavailable","WT","HET","HOM")
          )
          p <- DimPlot(
            seurat_object, reduction = "umap", group.by = feature,
            pt.size   = 1, alpha = 0.9
          ) + ggplot2::ggtitle(paste0(v)) +
            scale_color_manual(
              name   = "NGT",
              values = c(WT="lightblue", HET="darkgreen", HOM="red", Unavailable="grey80"),
              limits = c("WT","HET","HOM","Unavailable"),
              drop   = FALSE
            ) +
            theme(plot.title = element_text(size = umap_plot_title_size),
                  legend.title = element_text(size = umap_legend_title_size),
                  legend.text = element_text(size=umap_legend_label_size),
                  axis.title = element_text(size = umap_axis_title_size),
                  axis.text = element_text(size = umap_axis_lables_size)
            ) +
            guides(fill="none")
          
          # dealing with issue generating multiple legends
          if(v!=names(variants)[1]){
            p <- p + guides(colour="none")
          }
          
        }else{
          p <- suppressMessages(FeaturePlot(seurat_object, features = feature, reduction = "umap", order = TRUE, pt.size = 1, ) +
            ggplot2::ggtitle(paste0(type,": ", variant)) +
            scale_colour_gradient("AF", limits = c(0,1), breaks = c(0,0.5,1), low = "grey85", high = "red",
                                   oob = scales::squish, na.value = "grey85") +
            theme(plot.title = element_text(size = umap_plot_title_size),
                  legend.title = element_text(size = umap_legend_title_size),
                  legend.text = element_text(size=umap_legend_label_size),
                  axis.title = element_text(size = umap_axis_title_size),
                  axis.text = element_text(size = umap_axis_lables_size)
                  )
            )
        }
        p
        
      }
      
      names(plots) <- str_extract(plot_features, "^([^_]+)_", group = 1)
      
      plots
      
    }
    
    names(umap_plots) <- variants
    
    annotated_umap <- DimPlot(seurat_object,
                              reduction = "umap",
                              group.by  = "final_annotation",
                              label     = TRUE,
                              repel     = TRUE, label.box = T, pt.size=1.5, label.size = umap_ref_label_size
    ) +
      ggtitle("") +
      theme(axis.title = element_text(size=umap_ref_axis_title_size),
            axis.text = element_text(size = umap_ref_axis_lables_size)) 

    
    NGT_plots <- .assemble_feature_grid(variant_plot_list=umap_plots, feature = "NGT", reference_umap=annotated_umap, 
                                        ncol = if(length(variants) < 3) length(variants) else 3)
    ggsave(plot = NGT_plots, file.path(plot_directory,paste0(plot_prefix,"_umap_variants_NGT.png")),  width = 40, height = 50, units = "cm")
    message(paste0("Saved genotype UMAP plot to "), file.path(plot_directory,paste0(plot_prefix,"_umap_variants_NGT.png")))
    
    AF_plots <- .assemble_feature_grid(variant_plot_list=umap_plots, feature = "AF", reference_umap=annotated_umap, 
                                       ncol = if(length(variants) < 3) length(variants) else 3)
    ggsave(plot = AF_plots, file.path(plot_directory,paste0(plot_prefix,"_umap_variants_AF.png")),  width = 40, height = 50, units = "cm")
    message(paste0("Saved allele-frequency UMAP plot to "), file.path(plot_directory,paste0(plot_prefix,"_umap_variants_AF.png")))
  }
}




.detect_barcodes <- function(gt_file,genotype_format_fields){
  head_file <- data.table::fread(gt_file, nrows = 10L, sep = "\t") %>%
    as.matrix()
  cols <- colnames(head_file)
  matches_gt_format <- grepl(head_file,pattern=paste(rep("([0-9.]+|NA|NaN)",length(genotype_format_fields)), collapse = ";"))
  dim(matches_gt_format) <- dim(head_file)
  matches_gt_format <- apply(matches_gt_format,1, rle)
  barcode_length <- lapply(matches_gt_format, function(x) tail(x$lengths,1)) %>%
    unlist()
  if(all(barcode_length==barcode_length[1])){
    barcodes <- cols[(length(cols)-barcode_length[1]+1):length(cols)]
  }else{
    stop("Error: Couldn't detect barcodes from genotype file")
  }
  return(barcodes)
}




.assemble_feature_grid <- function(variant_plot_list,
                                  feature = c("NGT","AF","DP","GQ"),
                                  reference_umap,
                                  ncol = 4,
                                  title = NULL,
                                  collect_legend = TRUE,
                                  add_titles = TRUE,
                                  variant_order = NULL) {
  
  feature <- as.character(feature)
  
  # variant names for panel titles
  vnames <- names(variant_plot_list) %||% paste0("Variant_", seq_along(variant_plot_list))
  
  # optional ordering
  if (!is.null(variant_order)) {
    keep <- vnames %in% variant_order
    variant_plot_list <- variant_plot_list[keep]
    vnames <- vnames[keep]
    ord <- match(variant_order, vnames)
    variant_plot_list <- variant_plot_list[ord[!is.na(ord)]]
    vnames <- vnames[ord[!is.na(ord)]]
  }
  
  # extract the chosen feature plot from each variant
  plots <- Map(function(el, nm) {
    p <- el[[feature]]
    if (is.null(p)) return(NULL)
    if (add_titles) p <- p + labs(title = nm)
    p
  }, variant_plot_list, vnames)
  
  # drop missing
  plots <- plots[!vapply(plots, is.null, logical(1))]
  if (!length(plots)) stop("No plots found for feature '", feature, "'.")
  
  
  # assemble patchwork grid
  wrap_plots(
    plots,
    ncol = ncol,
    guides = if (isTRUE(collect_legend)) "collect" else "keep"
  ) | (reference_umap + Seurat::NoLegend()) +
    plot_annotation(title = feature) &
    theme(plot.title = element_text(hjust = 0.5))
}


.assemble_all_features <- function(variant_plot_list,
                                  features = c("NGT","AF","DP","GQ"),
                                  ncol = 4,
                                  collect_legend = TRUE,
                                  add_titles = TRUE,
                                  variant_order = NULL) {
  res <- lapply(features, function(ftr) {
    .assemble_feature_grid(
      variant_plot_list,
      feature = ftr,
      ncol = ncol,
      title = ftr,
      collect_legend = collect_legend,
      add_titles = add_titles,
      variant_order = variant_order
    )
  })
  names(res) <- features
  res
}


.standardise_ngt_legend <- function(p) {
  p +
    scale_color_manual(
      name   = "NGT",
      values = ngt_cols,
      limits = ngt_levels,
      drop   = FALSE
    ) +
    scale_fill_manual(
      name   = "NGT",
      values = ngt_cols,
      limits = ngt_levels,
      drop   = FALSE,
      guide  = "none"
    )
}
