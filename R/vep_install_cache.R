#' Install/populate a VEP cache using your chosen backend
#' Works with singularity or system/conda. Run on a node with internet.
#' @param species e.g. "homo_sapiens"
#' @param assembly e.g. "GRCh38" or "GRCh37"
#' @param cache_version optional Ensembl cache version (e.g. 114)
#' @param plugins also fetch Plugins/
#' @param fasta also fetch reference FASTA
#' @export
vep_install_cache <- function(species = .get_cfg()$species,
                              assembly = .get_cfg()$assembly,
                              cache_version = NULL,
                              plugins = TRUE,
                              fasta = TRUE) {
  cfg  <- .get_cfg()
  mode <- if (cfg$mode == "auto") .autodetect_mode(cfg) else cfg$mode
  
  if (!dir.exists(cfg$cache_dir)) dir.create(cfg$cache_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Build AUTO letters for INSTALL.pl (-a c[f][p])
  auto <- paste0("c", if (isTRUE(fasta)) "f" else "", if (isTRUE(plugins)) "p --PLUGINS all" else "")
  
  if (mode == "singularity") {
    sing <- .which("singularity"); if (is.null(sing)) stop("singularity not on PATH")
    if (is.null(cfg$image)) stop("Set vep_config(image=...) to the .sif when using singularity")
    
    # Bind host cache to container path used by VEP
    binds <- c(paste0(normalizePath(cfg$cache_dir), ":/opt/vep/.vep"))
    
    # Use known path to INSTALL.pl in the official image
    install_pl <- "/opt/vep/src/ensembl-vep/INSTALL.pl"
    
    args <- c(
      "exec",
      unlist(lapply(binds, function(b) c("-B", b))),
      cfg$image,
      "perl", install_pl,
      "-c", "/opt/vep/.vep",
      "-a", auto,
      "-s", species,
      "-y", assembly,
      "--NO_UPDATE", "--CONVERT"#, "--QUIET"
    )
    if (!is.null(cache_version)) args <- c(args, "--CACHE_VERSION", as.character(cache_version))
    
    .run_system("singularity", args)
    
  } else if (mode %in% c("system", "conda")) {
    # Locate INSTALL.pl near the vep binary
    vep_bin <- if (mode == "system") .which("vep") else file.path(cfg$conda_env, "bin", "vep")
    if (is.null(vep_bin) || !file.exists(vep_bin))
      stop("vep binary not found for ", mode, " mode")
    
    candidates <- c(
      file.path(dirname(vep_bin), "INSTALL.pl"),
      Sys.glob(file.path(dirname(vep_bin), "..", "share", "ensembl-vep*", "INSTALL.pl")),
      Sys.glob(file.path(dirname(vep_bin), "..", "src", "ensembl-vep", "INSTALL.pl"))
    )
    install_pl <- candidates[file.exists(candidates)][1]
    if (is.na(install_pl)) stop("Could not locate INSTALL.pl; install full VEP or use singularity.")
    
    args <- c(
      basename(install_pl),
      "-c", normalizePath(cfg$cache_dir),
      "-a", auto,
      "-s", species,
      "-y", assembly,
      "--NO_UPDATE", "--CONVERT"#, "--QUIET"
    )
    if (!is.null(cache_version)) args <- c(args, "--CACHE_VERSION", as.character(cache_version))
    
    withr::with_dir(dirname(install_pl), {
      .run_system("perl", args)
    })
    
  } else {
    stop("This helper supports singularity or system/conda. Configure one of those modes.")
  }
  
  # Sanity check
  if (!vep_has_cache(cfg$cache_dir, species, assembly, cache_version)) {
    warning("Cache installation finished but expected files not detected in ", cfg$cache_dir,
            ". Check binds/permissions and re-run without --QUIET for details.")
  } else {
    message("VEP cache ready in: ", cfg$cache_dir)
  }
  invisible(TRUE)
}
