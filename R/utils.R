# Internal validation, path, and dependency helpers.

.tfregact_stan_file <- function(stan_file) {
  if (!is.character(stan_file) || length(stan_file) != 1L || !nzchar(stan_file)) {
    stop("`stan_file` must be one non-empty file name or path.", call. = FALSE)
  }
  # A bare model name denotes a package model. Prefer the development/package
  # copy before looking in the working directory, where legacy scripts may
  # contain an older Stan file with the same basename. Explicit paths retain
  # their normal override behavior.
  is_bare_name <- identical(dirname(stan_file), ".") &&
    identical(basename(stan_file), stan_file)
  if (!is_bare_name && file.exists(stan_file)) {
    return(normalizePath(stan_file, winslash = "/", mustWork = TRUE))
  }
  development_candidates <- c(
    file.path(getwd(), "inst", "stan", basename(stan_file)),
    file.path(getwd(), "TFRegAct", "inst", "stan", basename(stan_file))
  )
  development_match <- development_candidates[file.exists(development_candidates)]
  if (length(development_match)) {
    return(normalizePath(
      development_match[[1]],
      winslash = "/",
      mustWork = TRUE
    ))
  }
  packaged <- system.file("stan", basename(stan_file), package = "TFRegAct")
  if (nzchar(packaged) && file.exists(packaged)) return(normalizePath(packaged, winslash = "/", mustWork = TRUE))
  if (file.exists(stan_file)) {
    return(normalizePath(stan_file, winslash = "/", mustWork = TRUE))
  }
  stop(sprintf("Stan model `%s` was not found in the package or at the supplied path.", stan_file), call. = FALSE)
}

.tfregact_require <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("Package `%s` is required for this operation. Install it with install.packages('%s').", package, package), call. = FALSE)
  }
  invisible(TRUE)
}

.tfregact_network_path <- function(path = NULL, must_work = TRUE) {
  candidate <- path
  if (is.null(candidate) || !length(candidate) || is.na(candidate[[1]]) || !nzchar(candidate[[1]])) {
    candidate <- getOption("TFRegAct.network_path", Sys.getenv("TFREGACT_NETWORK_PATH", unset = ""))
  }
  if (is.null(candidate) || !length(candidate) || !nzchar(candidate[[1]])) {
    stop("Supply `network_edge_file`, or set options(TFRegAct.network_path = '<path>').", call. = FALSE)
  }
  candidate <- path.expand(as.character(candidate[[1]]))
  if (isTRUE(must_work) && !file.exists(candidate)) {
    stop(sprintf("TF network file was not found: %s", candidate), call. = FALSE)
  }
  normalizePath(candidate, winslash = "/", mustWork = isTRUE(must_work))
}

.tfregact_network_cache_dir <- function() {
  cache_dir <- tools::R_user_dir("TFRegAct", which = "cache")
  normalizePath(cache_dir, winslash = "/", mustWork = FALSE)
}

.tfregact_default_network_file <- function() {
  file.path(.tfregact_network_cache_dir(), "TF_Full_Map.RData")
}

.tfregact_stan_without_log_lik <- function(stan_file) {
  stan_file <- .tfregact_stan_file(stan_file)
  source_text <- paste(readLines(stan_file, warn = FALSE), collapse = "\n")
  stripped_text <- sub("\\ngenerated quantities \\{[[:space:][:print:]]*\\}[[:space:]]*$", "\n", source_text)
  if (identical(stripped_text, source_text)) return(stan_file)
  cache_dir <- file.path(.tfregact_network_cache_dir(), "stan")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  output <- file.path(cache_dir, paste0(tools::file_path_sans_ext(basename(stan_file)), "_no_log_lik.stan"))
  if (!file.exists(output) || !identical(paste(readLines(output, warn = FALSE), collapse = "\n"), stripped_text)) {
    writeLines(stripped_text, output, useBytes = TRUE)
  }
  normalizePath(output, winslash = "/", mustWork = TRUE)
}
