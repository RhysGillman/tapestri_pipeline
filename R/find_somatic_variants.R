#' find_somatic_variants.R
#'
#' This function finds variants enriched in a particular annotated cell-type
#'
#' @param h5_in Path to input .h5 file
#' @param out_file Path to output file
#' @param block Block size (number of rows) for memory-efficient file processing
#' @return nothing
#' @export

find_somatic_variants <- function(h5_in,
                                  vcf_in,
                                  cell_type="all",
                                  cell_annotations,
                                  out_dir,
                                  file_prefix,
                                  block = 1000L,
                                  # Variant-level thresholds
                                  min_varcount_total        = 1,
                                  min_varportion_total      = 0,
                                  min_datacount_total       = 1,
                                  min_dataportion_total     = 0,
                                  min_quality_mean_total    = 0,
                                  min_allelefreq_mean_total = 0,
                                  # celltype-specific thresholds
                                  min_varcount_celltype     = 1,
                                  min_varportion_celltype   = 0,
                                  min_datacount_celltype    = 1,
                                  min_varportion_in_celltype = 0,
                                  max_p_value = 1,
                                  # Somatic / statistical thresholds
                                  odds_ratio_cutoff = -1000,
                                  zscore_cutoff     = 0,
                                  # Sorting options
                                  celltype_sort_value = "p",
                                  all_sort_column      = "mean_GQ_total",
                                  overwrite=TRUE
                          ) {
  # check paths
  
  all_data_file <- file.path(out_dir,paste0(file_prefix,"_all.tsv"))
  summary_file <- file.path(out_dir,paste0(file_prefix,"_all_summary.tsv"))
  summary_pass_file <- file.path(out_dir,paste0(file_prefix,"_all_summary_pass_only.tsv"))
  
  for(path in c(all_data_file,summary_file,summary_pass_file)){
    if(file.exists(path)){
      if(overwrite){
        message(paste0("Removing ", path, "..."))
        file.remove(path)
      }else{
        stop(paste0(path," already exists!"))
      }
    }
  }
  
  
  all_cell_types <- sort(unique(cell_annotations))
  
  if (cell_type != "all") {
    analysis_cell_types <- intersect(all_cell_types, cell_type)
  }else{
    analysis_cell_types <- all_cell_types
  }
  
  if (length(analysis_cell_types)==0) stop("Selected cell type doesn't match available data")
  
  
  for (ct in analysis_cell_types) {
    
    ct_file <- file.path(out_dir,paste0(file_prefix,"_",ct,".tsv"))
    ct_pass_file <- file.path(out_dir,paste0(file_prefix,"_",ct,"_pass_only.tsv"))
    
    for(path in c(ct_file,ct_pass_file)){
      if(file.exists(path)){
        if(overwrite){
          message(paste0("Removing ", path, "..."))
          file.remove(path)
        }else{
          stop(paste0(path," already exists!"))
        }
      }
    }
  }
  
  # -------------------process variants ----------------------
  message("Processing per-variant info...")
  variant_info <- h5read(h5_in, "/assays/dna_variants/ca")[c("CHROM","POS","REF","ALT")] %>% as.data.frame()
  
  variants <- variant_info %>%
    mutate(
      # calculate lengths
      ref_len = nchar(REF),
      alt_len = nchar(ALT),
      # skip ambiguous calls
      skip = str_detect(ALT, fixed("*")) | str_detect(REF, fixed("*")) | ALT %in% c(".", "") | REF %in% c(".", ""),
      # detect INS/DEL/SNV based on ref/alt lengths
      ins = alt_len > ref_len,
      del = ref_len > alt_len,
      snv = alt_len == 1 & ref_len == 1,
      # prepare trimmed REF/ALT for ins/del
      ALT_trim = sub("(?i)^[ATGCN]", "", ALT, perl = TRUE),
      REF_trim = sub("(?i)^[ATGCN]", "", REF, perl = TRUE),
      # Apply the trims
      ALT = if_else(ins, ALT_trim, ALT),
      REF = if_else(del, REF_trim, REF),
      # classify variants
      var_type = case_when(ins ~ "INS", del ~ "DEL", snv ~ "SNV", TRUE ~ NA_character_),
      # fix pos
      POS2 = if_else(del, POS + 1L, POS),
      END  = if_else(del, POS + nchar(REF), POS),
      # generate key
      var_key  = case_when(
        ins ~ paste0(CHROM,":",POS2,"-",END,"+", ALT),
        del ~ paste0(CHROM,":",POS2,"-",END,"-", REF),
        snv ~ paste0(CHROM,":",POS2,"-",END,REF, ">", ALT),
        TRUE ~ NA_character_
      )
      
    ) %>%
    mutate(filter=skip | is.na(var_type) | (ref_len == alt_len & ref_len > 1)) %>%
    dplyr::select(chromosome=CHROM, start_position = POS2, end_position=END, reference_allele=REF, alternate_allele=ALT, variant_type=var_type, variant=var_key,filter)
  
  # store which variants are being kept
  keep_idx <- !variants$filter
  
  if(any(!keep_idx)){
    message(paste0("WARNING: ", length(which(!keep_idx)),"/",nrow(variants), " variants are being filtered out"))
  }
  
  # filter out unwanted variants
  variants <- variants %>%
    filter(!filter) %>%
    dplyr::select(-filter)
  
  
  # ---------- Process Genotypes -------------------------------
  message("Processing per-cell info...")
  # this data is all per-cell
  # working with HDF5Array for memory efficiency 
  
  ds_AF  <- "/assays/dna_variants/layers/AF"
  ds_DP  <- "/assays/dna_variants/layers/DP"
  ds_GQ  <- "/assays/dna_variants/layers/GQ"
  ds_NGT <- "/assays/dna_variants/layers/NGT"
  
  AF_raw <- HDF5Array(h5_in, ds_AF)[keep_idx, ]
  DP     <- HDF5Array(h5_in, ds_DP)[keep_idx, ]
  GQ     <- HDF5Array(h5_in, ds_GQ)[keep_idx, ]
  NGT    <- HDF5Array(h5_in, ds_NGT)[keep_idx, ]
  
  nr <- nrow(variants)
  nc <- ncol(NGT)
  
  # get cell barcodes
  barcodes <- as.vector(rhdf5::h5read(h5_in, "/assays/dna_variants/ra/barcode"))
  if (length(barcodes) != nc) stop("barcode length != number of cells")
  if (length(barcodes) != length(cell_annotations) | !identical(sort(names(cell_annotations)),sort(barcodes))) stop("variant barcodes don't match cell annotations")
  
  dimnames(AF_raw) <- dimnames(DP) <- dimnames(GQ) <- dimnames(NGT) <- list(NULL, barcodes)
  
  # Probe AF to find scale
  probe_rows <- seq_len(min(16L, nr))
  probe_cols <- seq_len(min(32L, nc))
  af_probe   <- as.matrix(AF_raw[probe_rows, probe_cols, drop = FALSE])
  af_in_percent <- any(is.finite(af_probe) & af_probe > 1)
  
  
  # write header for cell-specific file
  writeLines("#FORMAT=AF;DP;GQ;NGT", all_data_file)
  
  # work in blocks for memory efficiency
  message(paste0("Found ", nr, " variants"))
  message(paste0("Working in blocks of ", block, " variants"))
  first_block <- TRUE
  # ---------- Stream row blocks spanning ALL columns ----------
  for (i0 in seq(1L, nr, by = block)) {
    
    message(paste0("Processing variants ", i0, "-", min(i0 + block - 1L, nr), " / ", nr))
    
    rows <- i0:min(i0 + block - 1L, nr)
    
    # get variant data for the block
    block_variants <- variants[rows, , drop = FALSE]
    
    # materialize only this block into RAM
    AFb_raw <- as.matrix(AF_raw[rows, , drop = FALSE])# used for per-cell "AF="
    AFb     <- if (af_in_percent) AFb_raw / 100 else AFb_raw  # 0–1 for stats
    DPb     <- as.matrix(DP[rows, , drop = FALSE])
    GQb     <- as.matrix(GQ[rows, , drop = FALSE])
    NGTb    <- as.matrix(NGT[rows, , drop = FALSE])
    
    # Totals across all cells (for "other" computations)
    alt_cnt_total  <- rowSums(NGTb == 1L | NGTb == 2L, na.rm = TRUE) # how many cells carry alt allele
    ref_cnt_total  <- rowSums(NGTb == 0L, na.rm = TRUE) # how many cells carry ref allele
    data_cnt_total <- alt_cnt_total + ref_cnt_total # how many cells have non-missing data
    alt_proportion_total  <- ifelse(data_cnt_total > 0, alt_cnt_total / data_cnt_total, 0)
    data_proportion_total <- ifelse(data_cnt_total > 0, data_cnt_total/ nc, 0)
    
    ## Calculate means/medians only when AF>0 (mask)
    mask <- is.finite(AFb) & (AFb > 0)
    denom <- pmax(1L, rowSums(mask))
    mean_AF_total <- rowSums(AFb * mask, na.rm = TRUE) / denom
    mean_AF_total[rowSums(mask) == 0] <- NA_real_
    
    AFb_na <- AFb; AFb_na[!mask] <- NA_real_
    median_AF_total <- matrixStats::rowMedians(AFb_na, na.rm = TRUE)
    
    # GQ mean where AF>0
    GQuse <- GQb
    GQuse[!mask] <- NA_real_
    mean_GQ_total <- rowMeans(GQuse, na.rm = TRUE)
    
    # assemble lines
    block_summary <- block_variants %>%
      bind_cols(
        data.frame(
          mean_GQ_total,
          mean_AF_total,
          median_AF_total,
          alt_cnt_total,
          data_cnt_total,
          alt_proportion_total
        )
      ) %>%
    
      # add filters
      
      rowwise() %>%
      mutate(
        filter = {
          labs <- c(
            "varcount_total"      = alt_cnt_total        < min_varcount_total,
            "varportion_total"    = alt_proportion_total < min_varportion_total,
            "datacount_total"     = data_cnt_total       < min_datacount_total,
            "dataportion_total"   = data_proportion_total< min_dataportion_total,
            "quality_mean_total"  = mean_GQ_total        < min_quality_mean_total,
            "allelefreq_mean_total" = mean_AF_total      < min_allelefreq_mean_total
          )
          out <- paste(names(which(labs)), collapse = ";")
          if (out == "") "." else out
        }
      ) %>%
      ungroup()
    
    write_tsv(block_summary, summary_file, append = TRUE, col_names = first_block)
    
    write_tsv(block_summary %>%
                filter(filter==".") %>%
                dplyr::select(-filter), summary_pass_file, append = TRUE, col_names = first_block)
    
    # get per cell info
    
    block_per_cell_info <- array(
      sprintf("%s;%s;%s;%s", round(as.vector(AFb),3), as.vector(DPb), as.vector(GQb), as.vector(NGTb)),
      dim = dim(AFb), dimnames = list(NULL, barcodes)
    )
    block_per_cell_info <- as.data.frame(block_per_cell_info, check.names = FALSE)
    
    block_per_cell_info <- block_summary %>%
      bind_cols(block_per_cell_info)
    
    write_tsv(block_per_cell_info, file.path(out_dir,paste0(file_prefix,"_all.tsv")), append = TRUE, col_names = first_block)
    
    
    
    for (ct in analysis_cell_types) {
      
      ct_file <- file.path(out_dir,paste0(file_prefix,"_",ct,".tsv"))
      ct_pass_file <- file.path(out_dir,paste0(file_prefix,"_",ct,"_pass_only.tsv"))
      
      barcode_idx <- which(cell_annotations == ct)
      
      # Counts within this celltype
      NGTb_ct <- NGTb[, barcode_idx, drop = FALSE] # get NGT for celltype
      alt_cnt_ct  <- rowSums(NGTb_ct == 1L | NGTb_ct == 2L, na.rm = TRUE)
      het_cnt_ct  <- rowSums(NGTb_ct == 1L, na.rm = TRUE)
      hom_cnt_ct <- rowSums(NGTb_ct == 2L, na.rm = TRUE)
      ref_cnt_ct  <- rowSums(NGTb_ct == 0L, na.rm = TRUE)
      no_data_cnt_ct  <- rowSums(NGTb_ct == 3L, na.rm = TRUE)
      data_cnt_ct <- alt_cnt_ct + ref_cnt_ct
      alt_proportion_ct  <- ifelse(data_cnt_ct > 0, alt_cnt_ct / data_cnt_ct, 0)
      alt_proportion_in_ct <- ifelse(alt_cnt_total > 0, alt_cnt_ct / alt_cnt_total, 0)
      data_proportion_ct <- ifelse(data_cnt_ct > 0, data_cnt_ct/ length(barcode_idx), 0)
      
      # Counts for "all other" cell types
      alt_cnt_other  <- alt_cnt_total - alt_cnt_ct
      ref_cnt_other <- ref_cnt_total - ref_cnt_ct
      data_cnt_other <- data_cnt_total - data_cnt_ct
      alt_proportion_other <- ifelse(data_cnt_other > 0, alt_cnt_other / data_cnt_other, 0)
      
      
      # Odds ratio + zscore
      #OR <- rep(NA, length(alt_cnt_ct))
      #Z_score <- rep(NA, length(alt_cnt_ct))
      #ok <- (alt_cnt_ct > 0 & ref_cnt_ct > 0 & alt_cnt_other > 0 & ref_cnt_other > 0)
      #OR[ok] <- (alt_cnt_ct[ok] / alt_cnt_other[ok]) / (ref_cnt_ct[ok] / ref_cnt_other[ok])
      #SE <- sqrt(1 / alt_cnt_ct[ok] + 1 / ref_cnt_ct[ok] + 1 / alt_cnt_other[ok] + 1 / ref_cnt_other[ok])
      #Z_score[ok]  <- log(OR[ok]) / SE[ok]
      # what if 100% cells carry alt or 0 "other cells" carry alt?
      #undefined <- ((alt_cnt_ct > 0 & alt_cnt_other == 0) | (alt_cnt_ct > 0 & ref_cnt_ct ==0 ) )
      #OR[undefined] <- Inf
      #Z_score[undefined] <- Inf
      
      n <- length(alt_cnt_ct)
      
      OR        <- rep(NA_real_, n)
      logOR     <- rep(NA_real_, n)
      SE        <- rep(NA_real_, n)
      Z_score   <- rep(NA_real_, n)
      p_wald    <- rep(NA_real_, n)
      ci_low    <- rep(NA_real_, n)
      ci_high   <- rep(NA_real_, n)
      p_fisher  <- rep(NA_real_, n)
      
      # rows with all four cells > 0
      ok <- (alt_cnt_ct > 0 & ref_cnt_ct > 0 & alt_cnt_other > 0 & ref_cnt_other > 0)
      
      # Odds ratio and Wald stats
      OR[ok]    <- (alt_cnt_ct[ok] / ref_cnt_ct[ok]) / (alt_cnt_other[ok] / ref_cnt_other[ok])
      logOR[ok] <- log(OR[ok])
      SE[ok]    <- sqrt( 1/alt_cnt_ct[ok] + 1/ref_cnt_ct[ok] + 1/alt_cnt_other[ok] + 1/ref_cnt_other[ok] )
      Z_score[ok] <- logOR[ok] / SE[ok]
      p_wald[ok]  <- 2 * pnorm(-abs(Z_score[ok]))
      
      # 95% CI for the OR (Wald)
      ci_low[ok]  <- exp(logOR[ok] - 1.96 * SE[ok])
      ci_high[ok] <- exp(logOR[ok] + 1.96 * SE[ok])
      
      # ---- Handle zeros/sparse cells ----
      need_fisher <- (alt_cnt_ct + ref_cnt_ct + alt_cnt_other + ref_cnt_other) > 0 &
        (alt_cnt_ct == 0 | ref_cnt_ct == 0 | alt_cnt_other == 0 | ref_cnt_other == 0 |
           pmin(alt_cnt_ct, ref_cnt_ct, alt_cnt_other, ref_cnt_other) < 5)
      
      if (any(need_fisher)) {
        fisher_idx <- which(need_fisher)
        p_fisher[fisher_idx] <- vapply(
          fisher_idx,
          function(i) {
            m <- matrix(c(alt_cnt_ct[i], ref_cnt_ct[i],
                          alt_cnt_other[i], ref_cnt_other[i]),
                        nrow = 2, byrow = TRUE)
            fisher.test(m, alternative = "two.sided")$p.value
          },
          numeric(1)
        )
      }
      
      p_final <- ifelse(!is.na(p_fisher), p_fisher, p_wald)
      
      stats_block_ct <- block_summary %>%
        bind_cols(data.frame(
          alt_cnt_ct,
          het_cnt_ct,
          hom_cnt_ct,
          ref_cnt_ct,
          data_cnt_ct,
          alt_proportion_in_ct,
          alt_cnt_other,
          ref_cnt_other,
          data_cnt_other,
          alt_proportion_other,
          OR,
          ci_low = ci_low,
          ci_high = ci_high,
          Z_score,
          p_wald = p_wald,
          p_fisher = p_fisher,
          p = p_final
          
        )) #%>%
        # filtering
        #rowwise() %>%
        #mutate(
          # Build the new, cell-type-specific filter labels
        #  .new_filter_labels = {
        #    labs <- c(
        #      "varcount_celltype"       = alt_cnt_ct  < min_varcount_celltype,
        #      "varportion_celltype"     = alt_proportion_in_ct  < min_varportion_celltype,
        #      "datacount_celltype"      = data_cnt_ct   < min_datacount_celltype,
        #      "varportion_in_celltype"  = alt_proportion_in_ct < min_varportion_in_celltype,
        #      "odds_ratio_cutoff"       = OR < odds_ratio_cutoff,
        #      "zscore_cutoff"           = Z_score < zscore_cutoff,
        #      "p_value" = p_final > max_p_value
        #    )
        #    out <- paste(names(which(labs)), collapse = ";")
        #    if (out == "") "." else out
        #  },
          # Append to existing `filter` (which may be "." from your previous step)
        #  filter = {
        #    add <- .new_filter_labels
        #    if (add == ".") {
        #      filter
        #    } else if (filter == ".") {
        #      add
        #    } else {
        #      paste0(filter, ";", add)
        #    }
        #  }
        #) %>%
        #dplyr::select(-.new_filter_labels) %>%
        #ungroup() #%>%
        #arrange(desc(!! rlang::sym(c(celltype_sort_value))))
      

      write_tsv(stats_block_ct, ct_file, append = TRUE, col_names = first_block)
      
      #write_tsv(stats_block_ct %>%
      #            filter(filter==".") %>%
      #            dplyr::select(-filter), ct_pass_file, append = TRUE, col_names = first_block)
      
    }
    first_block <- FALSE
  }
  
  # sorting
  
  ## read in just sorting column
  ## sort
  ## get indexes
  ## rewrite file block-wise using index
  
  # padj
  ## read in p-value col
  ## adjust pvalues
  ## re-write block-wise using index
  
  message("Sorting and Adjusting P Values...")
  
  for (ct in analysis_cell_types) {
    ct_file <- file.path(out_dir,paste0(file_prefix,"_",ct,".tsv"))
    ct_pass_file <- file.path(out_dir,paste0(file_prefix,"_",ct,"_pass_only.tsv"))
    
    full_results <- fread(ct_file) %>%
      #adjust p-values
      adjust_pvalue(p.col="p",output.col = "padj",method = "BH") %>%
      # filtering
      rowwise() %>%
      mutate(
        # Build the new, cell-type-specific filter labels
        .new_filter_labels = {
          labs <- c(
            "varcount_celltype"       = alt_cnt_ct  < min_varcount_celltype,
            "varportion_celltype"     = alt_proportion_in_ct  < min_varportion_celltype,
            "datacount_celltype"      = data_cnt_ct   < min_datacount_celltype,
            "varportion_in_celltype"  = alt_proportion_in_ct < min_varportion_in_celltype,
            "odds_ratio_cutoff"       = OR < odds_ratio_cutoff,
            "zscore_cutoff"           = Z_score < zscore_cutoff,
            "p_value" = padj > max_p_value
          )
          out <- paste(names(which(labs)), collapse = ";")
          if (out == "") "." else out
        },
        # Append to existing `filter` (which may be "." from your previous step)
        filter = {
          add <- .new_filter_labels
          if (add == ".") {
            filter
          } else if (filter == ".") {
            add
          } else {
            paste0(filter, ";", add)
          }
        }
      ) %>%
      dplyr::select(-.new_filter_labels) %>%
      ungroup()
    
    
    if(celltype_sort_value=="p"){
      full_results <- full_results %>%
        arrange(padj)
    }
    
    
    write_tsv(full_results, ct_file, append = FALSE, col_names = TRUE)
    
    write_tsv(full_results %>%
                filter(filter==".") %>%
                dplyr::select(-filter), ct_pass_file, append = FALSE, col_names = TRUE)
    
    
  }
  
  
}
