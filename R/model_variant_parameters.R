model_variant_parameters <- function(h5_in=NULL, vcf_in=NULL, plot_dir, prep_vcf=FALSE, plot_prefix="", block=1000) {
  
  
  
  if(!is.null(h5_in)&!is.null(vcf_in)){
    stop("Supply only h5_in OR vcf_in")
  }
  
  suppressWarnings(rm(stats))
  
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
    
    nr <- as.integer(system(
      sprintf("zcat -f %s | grep -cv '^#'", vcf_in), intern = T
    ))
    
    tabix <- TabixFile(vcf_in, yieldSize = block)
    
    if(isOpen(tabix)){
      close(tabix)
    }
    
    param <- ScanVcfParam(
      fixed = c("ALT"),
      info = NA,     
      geno = c("AD","DP","GQ","GT")
    )
    
    open(tabix)
    
    row_counter <- 0
    first_block <- TRUE
    repeat {
      if (row_counter >= nr) break
      # read in vcf block
      v <- readVcf(tabix, param = param)
      if (nrow(v) == 0) break
      # expand multi-allelic vcf
      vcf <- expand(v)
      
      # find rows to match variants of interest
      current_rows <- (row_counter+1):(row_counter+nrow(vcf))
      
      message(paste0("Processing variants ", row_counter+1, "-", min(max(current_rows), nr), " / ", nr))
      
      row_counter <- row_counter + nrow(vcf)
      
      # read in AD array and split to ref/alt matrices
      ADb <- geno(vcf)$AD[,,1:2,drop=FALSE]
      ADb_ref <- ADb[,,1, drop=FALSE]
      dim(ADb_ref) <- dim(ADb)[1:2]
      ADb_alt <- ADb[,,2, drop=FALSE]
      dim(ADb_alt) <- dim(ADb)[1:2]
      # calculate DP from AD ref+alt
      trueDPb <- ADb_ref + ADb_alt
      # calculate AF
      AFb <- ADb_alt/trueDPb
      # read in GQ
      GQb <- geno(vcf)$GQ[,,drop=FALSE]
      # read in GT
      #GTb <- geno(vcf)$GT[,,drop=FALSE]
      
      NGTb <- .pull_expanded_ngt_from_vcf(v)[,,drop=FALSE]
      
      #nc <- ncol(NGTb)
      # Totals across all cells (for "other" computations)
      #alt_cnt_total  <- rowSums(NGTb == 1L | NGTb == 2L, na.rm = TRUE) # how many cells carry alt allele
      #ref_cnt_total  <- rowSums(NGTb == 0L, na.rm = TRUE) # how many cells carry ref allele
      #data_cnt_total <- alt_cnt_total + ref_cnt_total # how many cells have non-missing data
      #alt_proportion_total  <- ifelse(data_cnt_total > 0, alt_cnt_total / data_cnt_total, 0)
      #data_proportion_total <- ifelse(data_cnt_total > 0, data_cnt_total/ nc, 0)
      
      ## Calculate means/medians only when AF>0 (mask)
      mask <- is.finite(AFb) & (AFb > 0)
      #mask <- NGTb %in% c(1L,2L)
      dim(mask) <- dim(NGTb)
      denom <- pmax(1L, rowSums(mask))
      mean_AF_total <- rowSums(AFb * mask, na.rm = TRUE) / denom
      mean_AF_total[rowSums(mask) == 0] <- NA_real_
      
      #AF_na <- AFb; AF_na[!mask] <- NA_real_
      #median_AF_total <- matrixStats::rowMedians(AF_na, na.rm = TRUE)
      
      # GQ mean where AF>0
      GQuse <- GQb
      GQuse[!mask] <- NA_real_
      mean_GQ_total <- rowMeans(GQuse, na.rm = TRUE)
      
      tmp <- data.frame(mean_AF=mean_AF_total, mean_GQ=mean_GQ_total)
      
      if(!exists("stats")){
        stats <- tmp
      }else{
        stats <- bind_rows(stats,tmp)
      }
      
  }
  }else if(!is.null(h5_in)){
    # working with HDF5Array for memory efficiency 
    
    ds_AF  <- "/assays/dna_variants/layers/AF"
    ds_DP  <- "/assays/dna_variants/layers/DP"
    ds_GQ  <- "/assays/dna_variants/layers/GQ"
    ds_NGT <- "/assays/dna_variants/layers/NGT"
    
    AF_raw <- HDF5Array(h5_in, ds_AF)
    DP     <- HDF5Array(h5_in, ds_DP)
    GQ     <- HDF5Array(h5_in, ds_GQ)
    NGT    <- HDF5Array(h5_in, ds_NGT)
    
    nr <- nrow(AF_raw)
    nc <- ncol(AF_raw)
    
    # Probe AF to find scale
    probe_rows <- seq_len(min(16L, nr))
    probe_cols <- seq_len(min(32L, nc))
    af_probe   <- as.matrix(AF_raw[probe_rows, probe_cols, drop = FALSE])
    af_in_percent <- any(is.finite(af_probe) & af_probe > 1)
    
    # work in blocks for memory efficiency

    first_block <- TRUE
    # ---------- Stream row blocks spanning ALL columns ----------
    for (i0 in seq(1L, nr, by = block)) {
      
      message(paste0("Processing variants ", i0, "-", min(i0 + block - 1L, nr), " / ", nr))
      
      rows <- i0:min(i0 + block - 1L, nr)
      
      # materialize only this block into RAM
      AFb_raw <- as.matrix(AF_raw[rows, , drop = FALSE])# used for per-cell "AF="
      AFb     <- if (af_in_percent) AFb_raw / 100 else AFb_raw  # 0–1 for stats
      DPb     <- as.matrix(DP[rows, , drop = FALSE])
      GQb     <- as.matrix(GQ[rows, , drop = FALSE])
      NGTb    <- as.matrix(NGT[rows, , drop = FALSE])
      
      
      
      ## Calculate means/medians only when AF>0 (mask)
      mask <- is.finite(AFb) & (AFb > 0)
      denom <- pmax(1L, rowSums(mask))
      mean_AF_total <- rowSums(AFb * mask, na.rm = TRUE) / denom
      mean_AF_total[rowSums(mask) == 0] <- NA_real_
      
      #AF_na <- AFb; AF_na[!mask] <- NA_real_
      #median_AF_total <- matrixStats::rowMedians(AF_na, na.rm = TRUE)
      
      # GQ mean where AF>0
      GQuse <- GQb
      GQuse[!mask] <- NA_real_
      mean_GQ_total <- rowMeans(GQuse, na.rm = TRUE)
      
      tmp <- data.frame(mean_AF=mean_AF_total, mean_GQ=mean_GQ_total)
      
      if(!exists("stats")){
        stats <- tmp
      }else{
        stats <- bind_rows(stats,tmp)
      }
      
  }
  
  }
  
  ggplot(stats, aes(x=mean_AF)) +
    geom_histogram(bins=30, fill="blue",colour="black") +
    scale_y_continuous(trans = "log10") +
    theme_bw() +
    labs(x="Mean Allele Frequency", y="Count")
  
  ggsave(file.path(plot_dir, paste0(plot_prefix,"_mean_AF.png")))
  
  ggplot(stats, aes(x=mean_GQ)) +
    geom_histogram(bins=30, fill="blue",colour="black") +
    scale_y_continuous(trans = "log10") +
    theme_bw() +
    labs(x="Mean Genotype Quality", y="Count")
  
  ggsave(file.path(plot_dir, paste0(plot_prefix,"_mean_GQ.png")))
  
  
}