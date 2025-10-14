#' Run VEP and return either the output path, a data.frame, or a VCF object
#' @param input Path to VEP input file
#' @param output desired output path (.vcf recommended)
#' @param return Type of object that the function should return, one of "path", "data.frame", or "vcf"
#' @param args extra VEP CLI args (character vector)
#' @export

run_vep <- function(input,
                    output = tempfile(fileext = ".vep.vcf"),
                    return = c("path","data.frame","vcf"),
                    args = character()){
  stopifnot(file.exists(input))
  return <- match.arg(return)
  cfg <- .get_cfg()
  
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  
  mode <- if (cfg$mode == "auto") .autodetect_mode(cfg) else cfg$mode
  
  cmd <- switch(mode,
                system      = .vep_cmd_system(input, output, cfg, args),
                conda       = .vep_cmd_conda(input, output, cfg, args),
                singularity = .vep_cmd_singularity(input, output, cfg, args),
                docker      = .vep_cmd_docker(input, output, cfg, args),
                rest        = .vep_cmd_rest(input, output, cfg, args),
                stop("Unknown mode: ", mode)
  )
  
  .run_system(cmd$bin, cmd$args)
  
  switch(return,
         path = output,
         vcf  = VariantAnnotation::readVcf(output),
         "data.frame" = .parse_vep_vcf(output)
  )
}

#' Configure VEP
#' @param mode Execution backend for VEP. One of `"auto"`, `"system"`, `"conda"`, `"singularity"`, `"docker"`, or `"rest"`.
#'   `"auto"` picks a sensible mode based on what is available; `"system"` calls a local VEP binary;
#'   `"conda"` activates a conda env containing VEP; `"singularity"`/`"docker"` run VEP in a container;
#'   `"rest"` uses the Ensembl REST API (local paths like `vep_path`, `cache_dir`, `fasta` are ignored).
#'   Default: `"auto"`.
#' @param vep_path Path to the VEP executable or wrapper to call when `mode="system"` (e.g., `/usr/bin/vep`).
#'   Ignored for container/REST modes.
#' @param conda_env Name or path of the conda environment that provides VEP when `mode="conda"` (e.g., `"vep-114"` or
#'   `"/path/to/envs/vep"`). Ignored otherwise.
#' @param image Container image to use when `mode` is `"singularity"` (e.g., `/path/to/vep.sif`) or `"docker"`
#'   (e.g., `"ensemblorg/ensembl-vep:release_114.0"`). Ignored for other modes.
#' @param bind One or more bind-mounts for container modes. For Singularity, a character vector of paths or
#'   `"host:container"` mappings; for Docker, equivalent `-v host:container` semantics. Use this to expose your
#'   input/output/cache directories inside the container.
#' @param cache_dir Directory containing the VEP cache (e.g., `~/.vep`). Used by all local modes; not used by `mode="rest"`.
#'   Default comes from `Sys.getenv("VEP_CACHE")` or `~/.vep`.
#' @param fasta Absolute path to the reference FASTA used by VEP (must match `species`/`assembly` and be indexed).
#'   Optional for `mode="rest"`.
#' @param species Ensembl species identifier, e.g., `"homo_sapiens"`. Default: `"homo_sapiens"`.
#' @param assembly Genome assembly/build string, e.g., `"GRCh38"` or `"GRCh37"`. Default: `"GRCh38"`.
#' @param threads Number of CPU threads to request for VEP. Default: `max(1, parallel::detectCores() - 1L)`.
#' @param perl_path Optional character vector of directories to prepend to `PATH` before running VEP (helps locate `perl`
#'   and VEP dependencies in non-standard installs).
#' @param perl5lib Optional character vector (or single string) to prepend to `PERL5LIB` for locating Perl modules.
#' @param vep_env Optional named list or named character vector of extra environment variables to set for the VEP run
#'   (e.g., `c(HTTP_PROXY="http://...", PERL_MM_OPT="...")`).
#' @param extra_args Character vector of additional raw CLI flags to pass through to VEP (e.g., `c("--everything","--cache")`).
#' @param perl_bin Absolute path to the Perl interpreter to use (overrides whatever is found on `PATH`). Useful when VEP
#'   requires a specific Perl. Optional.
#' @param perl_inc Character vector of include directories to add to Perl's `@INC` (prepended for the VEP process).
#'   Optional.
#' @param cache_version Ensembl cache version as a string (e.g., `"114"`). If supplied, may be used to validate or
#'   construct cache paths/arguments; otherwise inferred by VEP where possible.
#' @export


vep_config <- function(
    mode = c("auto","system","conda","singularity","docker","rest"),
    vep_path = NULL,
    conda_env = NULL,
    image = NULL,
    bind = NULL,
    cache_dir = Sys.getenv("VEP_CACHE", "~/.vep"),
    fasta = NULL,
    species = "homo_sapiens",
    assembly = "GRCh38",
    threads = max(1, parallel::detectCores() - 1L),
    perl_path = NULL,      # directories to prepend to PATH (optional)
    perl5lib  = NULL,      # PERL5LIB prepend (optional)
    vep_env   = NULL,      # named list/vec of extra env (optional)
    extra_args = character(),
    perl_bin = NULL,       # <-- NEW: absolute path to perl interpreter
    perl_inc = character(),# <-- NEW: include dirs for @INC (vector)
    cache_version = NULL   # <-- NEW: e.g. "114"
){
  mode <- match.arg(mode)
  cfg <- list(
    mode=mode,
    vep_path=vep_path,
    conda_env=conda_env,
    image=image,
    bind=bind,
    cache_dir=path.expand(cache_dir),
    fasta=if (!is.null(fasta)) path.expand(fasta) else NULL,
    species=species,
    assembly=assembly,
    threads=as.integer(threads),
    perl_path=perl_path,
    perl5lib=perl5lib,
    vep_env=vep_env,
    extra_args=as.character(extra_args),
    perl_bin=perl_bin,
    perl_inc=as.character(perl_inc),
    cache_version=cache_version
  )
  options(myvep.config = cfg)
  invisible(cfg)
}


._guess_input_format <- function(path){
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("vcf","vcf.gz")) return(NULL)        # VEP auto-detects VCF
  # assume ensembl-style for tsv/ txt
  if (ext %in% c("tsv","txt")) return(c("--format","ensembl"))
  NULL
}

._infer_cache_version <- function(cfg){
  # try to infer from cache_dir/fasta path like ".../114_GRCh37/..."
  paths <- c(cfg$cache_dir, cfg$fasta)
  m <- regmatches(paths, regexpr("/(\\d{2,3})_GRCh", paths, perl=TRUE))
  v <- unique(gsub("[^0-9]", "", m))
  if (length(v) && nzchar(v[1])) return(v[1])
  NULL
}


.common_vep_args <- function(input, output, cfg){
  base <- c(
    "--offline","--cache",
    "--dir_cache", cfg$cache_dir,
    "--species", cfg$species,
    "--assembly", cfg$assembly,
    "--input_file", normalizePath(input),
    "--output_file", normalizePath(output),
    "--force_overwrite",
    "--vcf","--symbol","--canonical",
    "--fork", as.character(cfg$threads)
  )
  # cache version (explicit beats inferred)
  cv <- cfg$cache_version %||% ._infer_cache_version(cfg)
  if (!is.null(cv)) base <- c("--cache_version", cv, base)
  
  infmt <- ._guess_input_format(input)
  if (!is.null(infmt)) base <- c(infmt, base)
  if (!is.null(cfg$fasta)) base <- c(base, "--fasta", cfg$fasta)
  c(base, cfg$extra_args)
}



.build_vep_env <- function(cfg){
  if (!is.null(cfg$vep_env)) {
    # allow named vector/list; convert to KEY=VALUE strings
    kv <- unlist(cfg$vep_env, use.names = TRUE)
    return(sprintf("%s=%s", names(kv), unname(kv)))
  }
  env <- character(0)
  
  if (!is.null(cfg$perl_path) && length(cfg$perl_path)) {
    path <- paste(c(cfg$perl_path, Sys.getenv("PATH")), collapse=":")
    env <- c(env, sprintf("PATH=%s", path))
  }
  if (!is.null(cfg$perl5lib) && length(cfg$perl5lib)) {
    perl5lib <- paste(c(cfg$perl5lib, Sys.getenv("PERL5LIB")), collapse=":")
    env <- c(env, sprintf("PERL5LIB=%s", perl5lib))
  }
  # If PERL5LIB not provided but perl_inc was, synthesize PERL5LIB from perl_inc
  if (!any(startsWith(env, "PERL5LIB=")) && length(cfg$perl_inc)) {
    env <- c(env, sprintf("PERL5LIB=%s", paste(cfg$perl_inc, collapse=":")))
  }
  
  # give child a HOME to avoid $HOME expansion issues under RStudio Server
  env <- c(env, sprintf("HOME=%s", normalizePath("~")))
  env
}




# internal: fetch config (bootstraps from env vars if absent)
.get_cfg <- function(){
  cfg <- getOption("myvep.config")
  if (is.null(cfg)) {
    vep_config(
      mode = Sys.getenv("MYVEP_MODE","auto"),
      vep_path = nzchar(Sys.getenv("MYVEP_BIN")) %||% NULL,
      conda_env = nzchar(Sys.getenv("MYVEP_CONDA_ENV")) %||% NULL,
      image = nzchar(Sys.getenv("MYVEP_IMAGE")) %||% NULL,
      bind = strsplit(Sys.getenv("MYVEP_BIND",""), ",")[[1]] %||% NULL,
      cache_dir = Sys.getenv("VEP_CACHE", "~/.vep"),
      fasta = nzchar(Sys.getenv("MYVEP_FASTA")) %||% NULL,
      species = Sys.getenv("MYVEP_SPECIES","homo_sapiens"),
      assembly = Sys.getenv("MYVEP_ASSEMBLY","GRCh38"),
      threads = as.integer(Sys.getenv("MYVEP_THREADS","1")),
      extra_args = strsplit(Sys.getenv("MYVEP_EXTRA",""), " +")[[1]] %||% character()
    )
  }
  getOption("myvep.config")
}

# fallback operator
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) return(y)
  if (is.character(x)) {
    if (length(x) != 1L) return(x)
    if (is.na(x) || !nzchar(x)) return(y)
  }
  x
}

.onLoad <- function(lib, pkg){
  # initialise default config once the package loads
  .get_cfg()
}

# tiny helper that checks exit codes and prints stderr on failure
.run_system <- function(bin, args, stream=T){
  
  out_opt <- if (stream) "" else TRUE
  
  cfg <- getOption("myvep.config")
  env <- .build_vep_env(cfg)
  
  # 1) Try direct exec of vep
  out1 <- try(system2(bin, args, env = if(length(env)) env else NULL,
                      stdout = out_opt, stderr = out_opt), silent = TRUE)
  if (!inherits(out1, "try-error")) {
    status <- suppressWarnings(attr(out1, "status"))
    if (is.null(status)) status <- if (is.numeric(out1) && length(out1) == 1L) as.integer(out1) else 0L
    if (identical(status, 0L)) return(out1)
  }
  
  # 2) Fallback: run via explicit perl with include dirs (-I)
  perl <- cfg$perl_bin %||% unname(Sys.which("perl")) %||% "/usr/bin/perl"
  inc_flags <- if (length(cfg$perl_inc)) sprintf("-I%s", cfg$perl_inc) else character()
  out2 <- try(system2(perl, c(inc_flags, bin, args),
                      env = if(length(env)) env else NULL,
                      stdout = out_opt, stderr = out_opt), silent = TRUE)
  if (!inherits(out2, "try-error")) {
    status <- suppressWarnings(attr(out2, "status"))
    if (is.null(status)) status <- if (is.numeric(out2) && length(out2) == 1L) as.integer(out2) else 0L
    if (identical(status, 0L)) return(out2)
  }
  
  # 3) Fallback: run via explicit perl relying on PERL5LIB
  env3 <- env
  if (!any(startsWith(env3, "PERL5LIB=")) && length(cfg$perl_inc)) {
    env3 <- c(env3, sprintf("PERL5LIB=%s", paste(cfg$perl_inc, collapse=":")))
  }
  out3 <- try(system2(perl, c(bin, args),
                      env = if(length(env3)) env3 else NULL,
                      stdout = out_opt, stderr = out_opt), silent = TRUE)
  
  if (!inherits(out3, "try-error")) {
    status <- suppressWarnings(attr(out3, "status"))
    if (is.null(status)) status <- if (is.numeric(out3) && length(out3) == 1L) as.integer(out3) else 0L
    if (identical(status, 0L)) return(out3)
  }
  
  stop("VEP failed via direct exec and perl fallbacks.\n--- Last output ---\n",
       paste(capture.output(print(out3)), collapse = "\n"))
}

.which <- function(x){
  p <- Sys.which(x)
  if (!nzchar(p)) NULL else unname(p)
}

.autodetect_mode <- function(cfg){
  if (!is.null(cfg$vep_path) && file.exists(cfg$vep_path)) return("system")
  if (!is.null(.which("vep"))) return("system")
  if (!is.null(cfg$conda_env) && file.exists(file.path(cfg$conda_env,"bin","vep")))
    return("conda")
  if (!is.null(.which("singularity")) && !is.null(cfg$image)) return("singularity")
  if (!is.null(.which("docker")) && !is.null(cfg$image)) return("docker")
  "rest"
}


.vep_cmd_system <- function(input, output, cfg, user_args){
  vep_bin <- cfg$vep_path %||% .which("vep")
  if (is.null(vep_bin)) stop("VEP not found; set vep_config(vep_path=...).")
  if (!file.exists(vep_bin)) stop("vep_path does not exist: ", vep_bin)
  
  list(
    bin  = vep_bin,
    args = c(.common_vep_args(input, output, cfg), user_args)
  )
}


.autodetect_mode <- function(cfg){
  # prefer an explicit vep_path if it exists
  if (!is.null(cfg$vep_path) && file.exists(cfg$vep_path)) return("system")
  if (!is.null(.which("vep"))) return("system")
  if (!is.null(cfg$conda_env) && file.exists(file.path(cfg$conda_env,"bin","vep")))
    return("conda")
  if (!is.null(.which("singularity")) && !is.null(cfg$image)) return("singularity")
  if (!is.null(.which("docker")) && !is.null(cfg$image)) return("docker")
  "rest"
}


.vep_cmd_conda <- function(input, output, cfg, user_args){
  if (is.null(cfg$conda_env))
    stop("Provide conda_env for conda mode.")
  vep_bin <- cfg$vep_path %||% file.path(cfg$conda_env, "bin", "vep")
  if (!file.exists(vep_bin))
    stop("vep not found at ", vep_bin)
  list(bin = vep_bin,
       args = c(.common_vep_args(input, output, cfg), user_args))
}

.vep_cmd_singularity <- function(input, output, cfg, user_args){
  sing <- .which("singularity")
  if (is.null(sing)) stop("singularity not found on PATH")
  if (is.null(cfg$image)) stop("Set vep_config(image=...) to a .sif file")
  
  host_in  <- normalizePath(dirname(input))
  host_out <- normalizePath(dirname(output))
  binds <- unique(c(
    paste0(host_in,  ":", host_in),
    paste0(host_out, ":", host_out),
    paste0(normalizePath(cfg$cache_dir), ":/opt/vep/.vep")
  ))
  if (!is.null(cfg$fasta))
    binds <- c(binds, paste0(normalizePath(dirname(cfg$fasta)), ":", dirname(cfg$fasta)))
  if (length(cfg$bind)) binds <- unique(c(binds, cfg$bind))
  
  vep_args <- .common_vep_args(input, output, modifyList(cfg, list(cache_dir="/opt/vep/.vep")))
  list(
    bin  = sing,
    args = c("exec",
             unlist(lapply(binds, function(b) c("-B", b))),
             cfg$image, "vep",
             c(vep_args, user_args))
  )
}

.vep_cmd_docker <- function(input, output, cfg, user_args){
  d <- .which("docker")
  if (is.null(d)) stop("docker not found on PATH")
  if (is.null(cfg$image)) stop("Set image= for docker mode")
  host_in  <- normalizePath(dirname(input))
  host_out <- normalizePath(dirname(output))
  mounts <- unique(c(
    paste0(host_in,  ":", host_in),
    paste0(host_out, ":", host_out),
    paste0(normalizePath(cfg$cache_dir), ":/opt/vep/.vep")
  ))
  if (!is.null(cfg$fasta))
    mounts <- c(mounts, paste0(normalizePath(dirname(cfg$fasta)), ":", dirname(cfg$fasta)))
  if (length(cfg$bind)) mounts <- unique(c(mounts, cfg$bind))
  
  vep_args <- .common_vep_args(input, output, modifyList(cfg, list(cache_dir="/opt/vep/.vep")))
  list(
    bin  = d,
    args = c("run","--rm",
             unlist(lapply(mounts, function(m) c("-v", m))),
             cfg$image, "vep",
             c(vep_args, user_args))
  )
}

# REST fallback (stub) — small inputs only
.vep_cmd_rest <- function(input, output, cfg, user_args){
  stop("REST mode not implemented yet. Configure system/conda/singularity/docker.")
}


# Parse a VEP-annotated VCF into a long data.frame (one row per CSQ record)
.parse_vep_vcf <- function(output) {
  message("Parsing VEP-annotated VCF to data.frame...")
  stopifnot(file.exists(output))
  vcf <- VariantAnnotation::readVcf(output)
  
  # -- CSQ schema
  info_hdr <- VariantAnnotation::info(VariantAnnotation::header(vcf))
  if (!"CSQ" %in% rownames(info_hdr))
    stop("CSQ field not present in INFO header; is this a VEP-annotated VCF?")
  desc <- as.character(info_hdr["CSQ", "Description"])
  fmt  <- sub("^.*[Ff]ormat: *", "", desc)
  csq_fields <- strsplit(fmt, "\\|")[[1]]
  csq_fields <- trimws(csq_fields)
  
  # -- Per-variant block
  rr   <- rowRanges(vcf)
  n    <- length(rr)
  CHROM <- as.character(GenomeInfoDb::seqnames(rr))
  POS   <- as.integer(S4Vectors::start(rr))
  ID <- names(rr)
  REF   <- as.character(VariantAnnotation::ref(vcf))
  
 
  
  csqCL <- VariantAnnotation::info(vcf)$CSQ
  if (is.null(csqCL)) return(data.frame())
  stopifnot(length(csqCL) == n)
  
  idx <- rep.int(seq_len(n), S4Vectors::elementNROWS(csqCL))
  csq_entries <- unlist(csqCL, use.names = FALSE)
  if (!length(csq_entries)) return(data.frame())
  
  # split CSQ strings -> matrix -> data.frame
  split_csq <- function(s, nfields = length(csq_fields)) {
    parts <- strsplit(s, "\\|", fixed = FALSE)[[1]]
    length(parts) <- nfields
    stats::setNames(parts, csq_fields)
  }
  csq_mat <- do.call(rbind, lapply(csq_entries, split_csq))
  csq_df  <- as.data.frame(csq_mat, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Variant columns (expanded) — ALT comes from CSQ "Allele" (allele-specific)
  allele_col <- names(csq_df)[match("allele", tolower(names(csq_df)))]
  ALT <- if (!is.na(allele_col)) csq_df[[allele_col]] else NA_character_
  
  variant_df <- data.frame(
    CHROM = CHROM[idx],
    POS   = POS[idx],
    ID    = ID[idx],
    REF   = REF[idx],
    ALT   = ALT,
    #.vep_key = paste(CHROM[idx], POS[idx], REF[idx], ALT, sep = ":"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  
  # Final order: variant-first, then all CSQ fields
  out <- cbind(variant_df, csq_df, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}
