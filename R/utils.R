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

# Activity inference needs coefficient draws, not cell-length Stan output arrays.
# Keep the same predictor calculation as a model-local variable so neither
# log_mu nor log_lik is written for every posterior draw.
.tfregact_activity_compact_stan <- function(stan_file) {
  stripped <- .tfregact_stan_without_log_lik(stan_file)
  lines <- readLines(stripped, warn = FALSE)
  model_line <- which(trimws(lines) == "model {")
  declaration <- which(trimws(lines) == "vector[N] log_mu;")
  calculation <- grep("^  log_mu =", lines)
  if (length(model_line) != 1L || length(declaration) != 1L ||
      length(calculation) != 1L || calculation >= model_line) {
    stop("Unexpected activity Stan predictor layout.", call. = FALSE)
  }
  closing <- max(which(trimws(lines[seq_len(model_line - 1L)]) == "}"))
  predictor_lines <- lines[calculation:(closing - 1L)]
  compact <- lines[-c(declaration, calculation:(closing - 1L))]
  model_line <- which(trimws(compact) == "model {")
  compact <- append(compact, c("  vector[N] log_mu;", predictor_lines),
                    after = model_line)
  output <- sub("_no_log_lik[.]stan$", "_activity_compact.stan", stripped)
  if (identical(output, stripped)) {
    stop("Activity Stan output requires a generated-quantities block.", call. = FALSE)
  }
  if (!file.exists(output) || !identical(readLines(output, warn = FALSE), compact)) {
    writeLines(compact, output, useBytes = TRUE)
  }
  normalizePath(output, winslash = "/", mustWork = TRUE)
}

# CmdStanR 0.8 reads CSVs through external grep on Windows. MSYS startup can
# fail intermittently in parallel workers. Use the same CSV reader/parser,
# feeding it R-filtered text instead of a shell command, in this worker only.
.tfregact_native_cmdstan_csv <- function() {
  if (.Platform$OS.type != "windows" ||
      isTRUE(getOption("TFRegAct.native_csv_active"))) return(invisible(FALSE))
  ns <- asNamespace("cmdstanr")
  bridge <- new.env(parent = ns)
  bridge$.tf_native_draws <- function(file, ...) {
    lines <- readLines(file, warn = FALSE)
    data.table::fread(text = paste(lines[!startsWith(lines, "#")], collapse = "\n"), ...)
  }
  bridge$.tf_native_metadata <- function(file, ...) {
    lines <- readLines(file, warn = FALSE)
    data.table::fread(text = paste(lines[grepl("^[#a-zA-Z]", lines)], collapse = "\n"), ...)
  }
  replacements <- list(
    read_cmdstan_csv = ".tf_native_draws(output_file,",
    read_csv_metadata = ".tf_native_metadata(csv_file,"
  )
  patched <- list()
  for (name in names(replacements)) {
    fun <- get(name, envir = ns)
    code <- paste(deparse(body(fun)), collapse = "\n")
    pattern <- "data.table::fread(cmd = fread_cmd,"
    if (!grepl(pattern, code, fixed = TRUE)) return(invisible(FALSE))
    code <- gsub(pattern, replacements[[name]], code, fixed = TRUE)
    body(fun) <- parse(text = code)[[1L]]
    environment(fun) <- bridge
    patched[[name]] <- fun
  }
  for (name in names(patched)) {
    utils::assignInNamespace(name, patched[[name]], ns = "cmdstanr")
  }
  options(TFRegAct.native_csv_active = TRUE)
  invisible(TRUE)
}
