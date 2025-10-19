.detect_barcodes <- function(gt_file){
  genotype_fmt_line <- readLines(gt_file, n = 1L)
  genotype_format_fields <- strsplit(gsub("^#FORMAT=", "", genotype_fmt_line), ";", fixed = TRUE)[[1]]
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

.get_col_schema <- function(file, n_lines){
  head_file <- data.table::fread(file, nrows = n_lines)
  schema <- setNames(as.character(vapply(head_file, class, "")), names(head_file))
  
  override <- c("MAX_AF"="double","chromosome"="character")
  schema[names(override)] <- unname(override)
  
  return(schema)
}
