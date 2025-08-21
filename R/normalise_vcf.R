#' normalise_vcf.R
#'
#' This function loads variant information from a .h5 file and generates a vcf-like file with per-cell variant information
#'
#' @param h5_in Path to input .h5 file
#' @param out_file Path to output file
#' @param block Block size (number of rows) for memory-efficient file processing
#' @return nothing
#' @export

normalise_vcf <- function(h5_in, 
                          out_file, 
                          block = 1000L) {
  message("Creating normalised VCF")
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
      # generate key
      var_key  = case_when(
        ins ~ paste0("+", ALT),
        del ~ paste0("-", REF),
        snv ~ paste0(REF, "->", ALT),
        TRUE ~ NA_character_
      ),
      POS2 = if_else(del, POS + 1L, POS),
      END  = if_else(del, POS + nchar(REF), POS)
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
  dimnames(AF_raw) <- dimnames(DP) <- dimnames(GQ) <- dimnames(NGT) <- list(NULL, barcodes)
  
  # Probe AF to find scale
  probe_rows <- seq_len(min(16L, nr))
  probe_cols <- seq_len(min(32L, nc))
  af_probe   <- as.matrix(AF_raw[probe_rows, probe_cols, drop = FALSE])
  af_in_percent <- any(is.finite(af_probe) & af_probe > 1)
  
  
  # write header including cell barcodes
  if (file.exists(out_file)) {
    message(paste0("WARNING: Overwriting ", out_file))
    file.remove(out_file)
  }
  writeLines("#FORMAT=AF;DP;GQ;NGT", out_file)
  #writeLines(paste(c("chromosome", "start_position", "end_position", "reference_allele", "alternate_allele", "variant_type", "variant", barcodes), collapse = "\t"), out_file)
  
  # work in blocks for memory efficiency
  message(paste0("Found ", nr, " rows"))
  message(paste0("Working in blocks of ", block, " rows"))
  first_block <- TRUE
  # ---------- Stream row blocks spanning ALL columns ----------
  for (i0 in seq(1L, nr, by = block)) {
    
    message(paste0("Processing rows ", i0, "-", min(i0 + block - 1L, nr), " / ", nr))
    
    rows <- i0:min(i0 + block - 1L, nr)
    
    # materialize only this block into RAM
    AFb_raw <- as.matrix(AF_raw[rows, , drop = FALSE])# used for per-cell "AF="
    AFb     <- if (af_in_percent) AFb_raw / 100 else AFb_raw  # 0–1 for stats
    DPb     <- as.matrix(DP[rows, , drop = FALSE])
    GQb     <- as.matrix(GQ[rows, , drop = FALSE])
    NGTb    <- as.matrix(NGT[rows, , drop = FALSE])
    
    # summary info
    no_data <- rowSums(NGTb == 3L, na.rm = TRUE)
    data_cnt <- nc - no_data
    alt_cnt  <- rowSums(NGTb == 1L | NGTb == 2L, na.rm = TRUE)
    
    ## Calculate means/medians only when AF>0 (mask)
    mask <- is.finite(AFb) & (AFb > 0)
    denom <- pmax(1L, rowSums(mask))
    mean_af <- rowSums(AFb * mask, na.rm = TRUE) / denom
    mean_af[rowSums(mask) == 0] <- NA_real_
    
    AFb_na <- AFb; AFb_na[!mask] <- NA_real_
    med_af <- matrixStats::rowMedians(AFb_na, na.rm = TRUE)
    
    # GQ mean where AF>0
    GQuse <- GQb
    GQuse[!mask] <- NA_real_
    mean_gq <- rowMeans(GQuse, na.rm = TRUE)
    
    stopifnot(identical(dim(AFb), dim(DPb)),
              identical(dim(AFb), dim(GQb)),
              identical(dim(AFb), dim(NGTb)))
    
    AFb  <- ifelse(is.finite(AFb), formatC(AFb, digits = 6, format = "fg"), ".")
    DPb  <- ifelse(is.finite(DPb), as.character(as.integer(DPb)), ".")
    GQb  <- ifelse(is.finite(GQb), as.character(as.integer(GQb)), ".")
    NGTb <- ifelse(is.finite(NGTb), as.character(as.integer(NGTb)), ".")
    
    block_per_cell_info <- array(
      sprintf("%s;%s;%s;%s", as.vector(AFb), as.vector(DPb), as.vector(GQb), as.vector(NGTb)),
      dim = dim(AFb), dimnames = list(NULL, barcodes)
    )
    block_per_cell_info <- as.data.frame(block_per_cell_info, check.names = FALSE)
    
    
    # assemble lines
    block_variants <- variants[rows, , drop = FALSE] %>%
      bind_cols(
        data.frame(
          mean_GQ=mean_gq,
          mean_AF=mean_af,
          median_AF=med_af,
          alt_count=alt_cnt,
          data_count=data_cnt
        )
      ) %>%
      bind_cols(block_per_cell_info)
    
    write_tsv(block_variants, out_file, append = TRUE, col_names = first_block)
    first_block <- FALSE
  }
  
}
