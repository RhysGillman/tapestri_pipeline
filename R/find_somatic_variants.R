#' find_somatic_variants.R
#'
#' This function finds variants enriched in a particular annotated cell-type
#'
#' @param h5_in Path to input .h5 file
#' @param out_file Path to output file
#' @param block Block size (number of rows) for memory-efficient file processing
#' @return nothing
#' @export

find_somatic_variants <- function(h5_in=NULL,
                                  vcf_in=NULL,
                                  prep_vcf=F,
                                  prep_vcf_dir=NULL,
                                  cell_type="all",
                                  cell_annotations,
                                  out_dir,
                                  file_prefix,
                                  block = 1000,
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
                                  overwrite=TRUE,
                                  genome_version="hg19"
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
  
  # barcodes
  if(!is.null(vcf_in)){
    hdr <- scanVcfHeader(vcf_in)
    # genome version
    #  find_ref_genome <- meta(hdr)$reference
    #  if(is.null(find_ref_genome)){
    #    find_ref_genome <- meta(hdr)$contig 
    #  }
    
    barcodes <- samples(hdr)
    
  }else if(!is.null(h5_in)){
    barcodes <- as.vector(rhdf5::h5read(h5_in, "/assays/dna_variants/ra/barcode"))
  }
    
  if (length(barcodes) != length(cell_annotations)) stop("cell annotations are not the same length as vcf samples")
  
  if (!identical(sort(names(cell_annotations)),sort(barcodes))) {
    message("variant barcodes don't match cell annotations")
    message("attemtping to find match in substring...")
    barcode1 <- barcodes[1]
    cell1 <- names(cell_annotations)[1]
    match <- F
    for(i in 1:nchar(barcode1)){
      start_sub <- substr(barcode1,1,i)
      end_sub <- substr(barcode1,nchar(barcode1)-i,nchar(barcode1))
      if(start_sub==cell1){
        match <- T
        barcode_suffix <- substr(barcode1, i+1, nchar(barcode1))
        barcodes <- gsub(barcode_suffix,"",barcodes)
        if (identical(sort(names(cell_annotations)),sort(barcodes))){
          message(paste0("substring match found, trimming "), barcode_suffix," suffix from barcodes")
        }
      }else if (end_sub==cell1){
        match <- T
        barcode_prefix <- substr(barcode1,i, nchar(barcode1)-i-1)
        barcodes <- gsub(barcode_prefix,"",barcodes)
        if (identical(sort(names(cell_annotations)),sort(barcodes))){
          message(paste0("substring match found, trimming "), barcode_prefix," prefix from barcodes")
        }
      }
      if(match){
        break
      }
    }
    if(!match){
      stop("no substring match found")
    }
  }
    
  
  
  # -------------------process variants ----------------------
  message("Processing per-variant info...")
  
  normalise_vcf <- function(vcf_df){
    out <- vcf_df %>%
      mutate(
        # calculate lengths
        ref_len = nchar(REF),
        alt_len = nchar(ALT),
        # detect INS/DEL/SNV based on ref/alt lengths
        simple_ins = (ref_len == 1) & (alt_len > 1) & (substr(REF,1,1)==substr(ALT,1,1)),
        complex_ins = (ref_len > 1) & (alt_len > ref_len) & (substr(REF,1,1)==substr(ALT,1,1)),
        ins_len0 = (ref_len == 0) & (alt_len > ref_len),
        simple_del = (alt_len == 0) & (ref_len > alt_len),
        complex_del = (alt_len > 0) & (ref_len > alt_len) & (substr(REF,1,1)==substr(ALT,1,1)),
        simple_snv = (ref_len == 1) & (alt_len == 1) & (REF != ALT),
        complex_snv = (ref_len == alt_len) & (ref_len > 1) & (REF != ALT)) %>%
      rowwise() %>%
      mutate(
        complex_snv_pos = ifelse(complex_snv, which(str_split(ALT, pattern = "")[[1]] != str_split(REF, pattern = "")[[1]]), NA)
      ) %>%
      ungroup() %>%
      mutate(
        # skip ambiguous calls
        skip = str_detect(ALT, stringr::fixed("*")) | str_detect(REF, stringr::fixed("*")) | (!simple_ins&!complex_ins&!ins_len0&!complex_del&!simple_snv&!complex_snv),
        # handling positions
        START = case_when(
          simple_ins ~ POS,
          complex_ins ~ POS,
          ins_len0 ~ POS,
          simple_del ~ POS,
          complex_del ~ POS, # ATGAG ATG
          simple_snv ~ POS,
          complex_snv ~ POS + complex_snv_pos - 1
        ),
        END = case_when(
          simple_ins ~ START+1,
          complex_ins ~ START+1,
          ins_len0 ~ START+1,
          simple_del ~ START + ref_len - 1,
          complex_del ~ START + ref_len - alt_len,
          simple_snv ~ POS,
          complex_snv ~ POS + complex_snv_pos - 1
        ),
        var_type = case_when(
          simple_ins ~ "INS", 
          complex_ins ~ "INS", 
          ins_len0 ~ "INS", 
          simple_del ~ "DEL", 
          complex_del ~ "DEL", 
          simple_snv ~ "SNV", 
          complex_snv ~ "SNV",
          TRUE ~ NA_character_),
        new_REF = case_when(
          simple_ins ~ REF, 
          complex_ins ~ substr(REF,1,1), 
          ins_len0 ~ REF, 
          simple_del ~ REF, 
          complex_del ~ substr(REF,1,1+ref_len-alt_len), # ATGAG ATG
          simple_snv ~ REF, 
          complex_snv ~ substr(REF,complex_snv_pos,complex_snv_pos)
        ),
        new_ALT = case_when(
          simple_ins ~ ALT, 
          complex_ins ~ substr(ALT,1,1+alt_len-ref_len), 
          ins_len0 ~ ALT, 
          simple_del ~ ALT, 
          complex_del ~ substr(ALT,1,1), # ATGAG ATG
          simple_snv ~ ALT, 
          complex_snv ~ substr(ALT,complex_snv_pos,complex_snv_pos)
        ),
        var_key = case_when(
          simple_ins ~ paste0(CHROM,":",START+1,"+",substr(new_ALT,2,nchar(new_ALT))), 
          complex_ins ~ paste0(CHROM,":",START+1,"+",substr(new_ALT,2,nchar(new_ALT))), 
          ins_len0 ~ paste0(CHROM,":",START,"+",new_ALT), 
          simple_del ~ paste0(CHROM,":",START,"-",END,"-",new_REF), 
          complex_del ~ paste0(CHROM,":",START+1,"-",END,"-",substr(new_REF,2,nchar(new_REF))), 
          simple_snv ~ paste0(CHROM,":",START,new_REF,">",new_ALT), 
          complex_snv ~ paste0(CHROM,":",START,new_REF,">",new_ALT)
        )
      ) %>%
      mutate(filter=skip | is.na(var_type)) %>%
      dplyr::select(chromosome=CHROM, start_position = START, end_position=END, reference_allele=new_REF, alternate_allele=new_ALT, variant_type=var_type, variant=var_key,filter)
  }
  
  if(!is.null(vcf_in)){
    
    if(prep_vcf){
      message("Generating compressed vcf and indexing...")
      out_vcf <- file.path(prep_vcf_dir,gsub(".vcf",".vcf.gz",basename(vcf_in)))
      #module load htslib/1.16
      system(sprintf("bgzip -c %s > %s",vcf_in,out_vcf))
      system(sprintf("tabix -p vcf %s",out_vcf))
      vcf_in <- out_vcf
    }else if(!str_detect(vcf_in,".vcf.gz")){
      stop("vcf.gz required, use prep_vcf=TRUE")
    }else if(!file.exists(gsub(".vcf.gz|.vcf",".vcf.gz.tbi",vcf_in))){
      stop("vcf index required, use prep_vcf=TRUE")
    }
    
    tabix <- TabixFile(vcf_in)#, yieldSize = block*100)
    
    param <- ScanVcfParam(
      fixed = c("ALT"),
      info = NA,     
      geno = NA
    )
    
    open(tabix)
    suppressWarnings(rm(variant_info))
    repeat {
      v <- readVcf(tabix, genome = genome_version, param = param)
      if (nrow(v) == 0) break
      vcf <- expand(v)
      rr   <- rowRanges(vcf)
      
      # ---- base columns ----
      CHROM  <- as.character(seqnames(rr))
      POS    <- start(rr)
      REF    <- as.character(ref(vcf))
      ALT <- vapply(
        as.list(alt(vcf)),
        function(x) paste(as.character(x), collapse = ","),
        character(1)
      )
      
      if(exists("variant_info")){
        variant_info <- variant_info %>%
          bind_rows(data.frame(CHROM,POS,REF,ALT))
      }else{
        variant_info <- data.frame(CHROM,POS,REF,ALT)
      }
      
    }
    close(tabix)
    
    variants <- normalise_vcf(variant_info)
    
    
  }else if(!is.null(h5_in)){
    variant_info <- h5read(h5_in, "/assays/dna_variants/ca")[c("CHROM","POS","REF","ALT")] %>% as.data.frame()
    variants <- normalise_vcf(variant_info)
  }
  
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
  
  pull_expanded_ngt_from_vcf <- function(v) {
    GT <- geno(v)$GT # get unexpanded genotypes (includes multiallelic sites)
    nvar <- nrow(GT) # number of variants
    nsamp <- ncol(GT) # number of samples
    
    # count number of alt alleles per row/variant
    nALT <- elementNROWS(alt(v))
    
    # get per-sample alleles A and B
    pat <- "^([0-9]+|\\.)[\\|/]([0-9]+|\\.)$"
    # sample allele A
    A <- suppressWarnings(as.integer(sub(pat, "\\1", GT, perl = TRUE))) %>%
      matrix(nrow=nvar, ncol = nsamp)
    # sample allele B
    B <- suppressWarnings(as.integer(sub(pat, "\\2", GT, perl = TRUE))) %>%
      matrix(nrow=nvar, ncol = nsamp)
    A[GT == "./."] <- NA_integer_
    B[GT == "./."] <- NA_integer_
    
    # Build expanded row mapping (same order as expand(v))
    exp_rows <- sum(nALT) # nrows after expansion
    row_map  <- rep.int(seq_len(nvar), nALT)  # which original row for each expanded row
    alt_idx  <- unlist(lapply(nALT, seq_len), use.names = FALSE)
    
    # expand A and B rows
    Aexp <- A[row_map, , drop = FALSE]
    Bexp <- B[row_map, , drop = FALSE]
    # altid matrix map
    J <- matrix(alt_idx, nrow = exp_rows, ncol = nsamp)
    
    # Build NGT (0,1,2,3) per expanded row x sample
    NGT <- matrix(3L, nrow = exp_rows, ncol = nsamp,
                  dimnames = list(NULL, colnames(GT)))
    
    # missing stays 3
    not_na <- !(is.na(Aexp) | is.na(Bexp))
    
    # 0/0 → 0
    NGT[not_na & (Aexp == 0L) & (Bexp == 0L)] <- 0L
    # j/j → 2
    NGT[not_na & (Aexp == J)  & (Bexp == J)]  <- 2L
    # 0/j or j/0 → 1
    het <- not_na & (
      ((Aexp == 0L) & (Bexp == J)) |
        ((Aexp == J)  & (Bexp == 0L))
    )
    NGT[het] <- 1L
    
    NGT
  }
  
  
  
  calculate_variant_level_stats <- function(NGT,AF,GQ) {
    nc=ncol(NGT)
    # Totals across all cells (for "other" computations)
    alt_cnt_total  <- rowSums(NGT == 1L | NGT == 2L, na.rm = TRUE) # how many cells carry alt allele
    ref_cnt_total  <- rowSums(NGT == 0L, na.rm = TRUE) # how many cells carry ref allele
    data_cnt_total <- alt_cnt_total + ref_cnt_total # how many cells have non-missing data
    alt_proportion_total  <- ifelse(data_cnt_total > 0, alt_cnt_total / data_cnt_total, 0)
    data_proportion_total <- ifelse(data_cnt_total > 0, data_cnt_total/ nc, 0)
    
    ## Calculate means/medians only when AF>0 (mask)
    mask <- is.finite(AF) & (AF > 0)
    denom <- pmax(1L, rowSums(mask))
    mean_AF_total <- rowSums(AF * mask, na.rm = TRUE) / denom
    mean_AF_total[rowSums(mask) == 0] <- NA_real_
    
    AF_na <- AF; AF_na[!mask] <- NA_real_
    median_AF_total <- matrixStats::rowMedians(AF_na, na.rm = TRUE)
    
    # GQ mean where AF>0
    GQuse <- GQ
    GQuse[!mask] <- NA_real_
    mean_GQ_total <- rowMeans(GQuse, na.rm = TRUE)
    
    out <- data.frame(
      mean_GQ_total,
      mean_AF_total,
      median_AF_total,
      alt_cnt_total,
      ref_cnt_total,
      data_cnt_total,
      alt_proportion_total,
      data_proportion_total
    )
    return(out)
  }
  
  calculate_cell_type_stats <- function(NGT,cell_type,variant_stats){
    
    barcode_idx <- which(cell_annotations == cell_type)
    
    alt_cnt_total <- variant_stats$alt_cnt_total
    ref_cnt_total <- variant_stats$ref_cnt_total
    data_cnt_total <- variant_stats$data_cnt_total
    
    # Counts within this celltype
    NGT_ct <- NGT[, barcode_idx, drop = FALSE] # get NGT for celltype
    alt_cnt_ct  <- rowSums(NGT_ct == 1L | NGT_ct == 2L, na.rm = TRUE)
    het_cnt_ct  <- rowSums(NGT_ct == 1L, na.rm = TRUE)
    hom_cnt_ct <- rowSums(NGT_ct == 2L, na.rm = TRUE)
    ref_cnt_ct  <- rowSums(NGT_ct == 0L, na.rm = TRUE)
    no_data_cnt_ct  <- rowSums(NGT_ct == 3L, na.rm = TRUE)
    data_cnt_ct <- alt_cnt_ct + ref_cnt_ct
    alt_proportion_ct  <- ifelse(data_cnt_ct > 0, alt_cnt_ct / data_cnt_ct, 0)
    alt_proportion_in_ct <- ifelse(alt_cnt_total > 0, alt_cnt_ct / alt_cnt_total, 0)
    data_proportion_ct <- ifelse(data_cnt_ct > 0, data_cnt_ct/ length(barcode_idx), 0)
    
    # Counts for "all other" cell types
    alt_cnt_other  <- alt_cnt_total - alt_cnt_ct
    ref_cnt_other <- ref_cnt_total - ref_cnt_ct
    data_cnt_other <- data_cnt_total - data_cnt_ct
    alt_proportion_other <- ifelse(data_cnt_other > 0, alt_cnt_other / data_cnt_other, 0)
    
    n <- length(alt_cnt_ct)
    
    OR        <- rep(NA_real_, n)
    logOR     <- rep(NA_real_, n)
    SE        <- rep(NA_real_, n)
    Z_score   <- rep(NA_real_, n)
    p_wald    <- rep(NA_real_, n)
    ci_low    <- rep(NA_real_, n)
    ci_high   <- rep(NA_real_, n)
    p_fisher  <- rep(NA_real_, n)
    
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
    
    out <- data.frame(
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
      
    )
  }
  
  
  # write header for cell-specific file
  writeLines("#FORMAT=AF;DP;GQ;NGT", all_data_file)
  
  
  if(!is.null(vcf_in)){
    
    nr <- nrow(variants)
    
    
    tabix <- TabixFile(vcf_in, yieldSize = block)
    
    param <- ScanVcfParam(
      fixed = c("ALT"),
      info = NA,     
      geno = c("AD","DP","GQ","GT")
    )
    
    open(tabix)
    
    raw_row_counter <- 0
    variants_row_counter <- 0
    first_block <- TRUE
    repeat {
      
      # read in vcf block
      v <- readVcf(tabix, genome = genome_version, param = param)
      if (nrow(v) == 0) break
      # expand multi-allelic vcf
      vcf <- expand(v)
      
      # find rows to match variants of interest
      current_rows <- (raw_row_counter+1):(raw_row_counter+nrow(vcf))
      raw_row_counter <- raw_row_counter + nrow(vcf)
      keep_rows <- which(keep_idx[current_rows])
      
      # find which rows to pull from variants df
      variant_rows <- (variants_row_counter+1):(variants_row_counter+length(keep_rows))
      
      message(paste0("Processing variants ", variants_row_counter+1, "-", min(variants_row_counter + length(keep_rows), nr), " / ", nr))
      
      variants_row_counter <- variants_row_counter + length(keep_rows)
      
      # get variant data for the block
      block_variants <- variants[variant_rows, , drop = FALSE]
      
      # read in AD array and split to ref/alt matrices
      ADb <- geno(vcf)$AD[keep_rows,,1:2,drop=FALSE]
      barcodes <- colnames(ADb)
      ADb_ref <- ADb[,,1, drop=FALSE]
      dim(ADb_ref) <- dim(ADb)[1:2]
      ADb_alt <- ADb[,,2, drop=FALSE]
      dim(ADb_alt) <- dim(ADb)[1:2]
      # calculate DP from AD ref+alt
      trueDPb <- ADb_ref + ADb_alt
      # calculate AF
      AFb <- ADb_alt/trueDPb
      # read in GQ
      GQb <- geno(vcf)$GQ[keep_rows,,drop=FALSE]
      # read in GT
      GTb <- geno(vcf)$GT[keep_rows,,drop=FALSE]
      
      NGTb <- pull_expanded_ngt_from_vcf(v)[keep_rows,,drop=FALSE]
      
      # calculate variant-level stats
      var_stats <- calculate_variant_level_stats(NGT = NGTb,AF = AFb, GQ = GQb)
      
      # assemble lines
      block_summary <- block_variants %>%
        bind_cols(
          var_stats
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
        sprintf("%s;%s;%s;%s", round(as.vector(AFb),3), as.vector(trueDPb), as.vector(GQb), as.vector(NGTb)),
        dim = dim(AFb), dimnames = list(NULL, barcodes)
      )
      block_per_cell_info <- as.data.frame(block_per_cell_info, check.names = FALSE)
      
      block_per_cell_info <- block_summary %>%
        bind_cols(block_per_cell_info)
      
      write_tsv(block_per_cell_info, file.path(out_dir,paste0(file_prefix,"_all.tsv")), append = TRUE, col_names = first_block)
       
      for (ct in analysis_cell_types) {
        
        ct_file <- file.path(out_dir,paste0(file_prefix,"_",ct,".tsv"))
        ct_pass_file <- file.path(out_dir,paste0(file_prefix,"_",ct,"_pass_only.tsv"))
        
        ct_stats <- calculate_cell_type_stats(NGTb,ct,var_stats)
        
        stats_block_ct <- block_summary %>%
          bind_cols(ct_stats) 
        
        write_tsv(stats_block_ct, ct_file, append = TRUE, col_names = first_block)
        
      }
      first_block <- FALSE
    }
    close(tabix)
    
  }else if(!is.null(h5_in)){
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
    nc <- ncol(variants)
    
    dimnames(AF_raw) <- dimnames(DP) <- dimnames(GQ) <- dimnames(NGT) <- list(NULL, barcodes)
    
    # Probe AF to find scale
    probe_rows <- seq_len(min(16L, nr))
    probe_cols <- seq_len(min(32L, nc))
    af_probe   <- as.matrix(AF_raw[probe_rows, probe_cols, drop = FALSE])
    af_in_percent <- any(is.finite(af_probe) & af_probe > 1)
    
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
      
      var_stats <- calculate_variant_level_stats(NGT = NGTb,AF = AFb, GQ = GQb)
      
      # assemble lines
      block_summary <- block_variants %>%
        bind_cols(
          var_stats
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
        
        ct_stats <- calculate_cell_type_stats(NGTb,ct,var_stats)
        
        stats_block_ct <- block_summary %>%
          bind_cols(ct_stats) 
        
        write_tsv(stats_block_ct, ct_file, append = TRUE, col_names = first_block)
        
      }
      first_block <- FALSE
    }
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
