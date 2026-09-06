#!/usr/bin/env Rscript

# Complete two-stage TF-to-gene regulatory-direction workflow.
# Source this file, then call TF_regulatory_direction_computation().

.tf_regulatory_direction_dir <- local({
  source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  source_dir <- if (is.null(source_file) || !nzchar(as.character(source_file[[1]]))) getwd() else dirname(normalizePath(source_file[[1]], winslash = "/", mustWork = FALSE))
  function() source_dir
})

.tf_regulatory_direction_stop <- function(...) stop(sprintf(...), call. = FALSE)

.tf_regulatory_direction_load_dependencies <- function() invisible(NULL)

.tf_regulatory_direction_load_seurat <- function(seurat_obj, data_file) {
  if (!is.null(seurat_obj) && !is.null(data_file)) {
    .tf_regulatory_direction_stop("Supply exactly one of `seurat_obj` or `data_file`, not both.")
  }
  if (!is.null(seurat_obj)) return(seurat_obj)
  if (!is.character(data_file) || length(data_file) != 1L || !file.exists(data_file)) {
    .tf_regulatory_direction_stop("Supply `seurat_obj`, or a valid `data_file` containing an object named `pbmc`.")
  }
  loaded <- new.env(parent = emptyenv())
  load(data_file, envir = loaded)
  if (!exists("pbmc", envir = loaded, inherits = FALSE)) {
    .tf_regulatory_direction_stop("`data_file` must contain a Seurat object named `pbmc`.")
  }
  get("pbmc", envir = loaded, inherits = FALSE)
}

.tf_regulatory_direction_output_dir <- function(output, target_tf, target_gene) {
  if (isFALSE(output) || is.null(output)) return(NULL)
  if (isTRUE(output)) return(file.path(getwd(), paste0(toupper(target_tf), "_", toupper(target_gene), "_regulatory_output")))
  if (is.character(output) && length(output) == 1L && nzchar(output)) return(normalizePath(output, winslash = "/", mustWork = FALSE))
  .tf_regulatory_direction_stop("`output` must be FALSE, TRUE, or a single non-empty output path.")
}

#' Estimate the direction and strength of one TF -> target-gene regulatory edge.
#'
#' Supply exactly one input dataset: `seurat_obj` (an in-memory Seurat object)
#' or `data_file` (an .RData file containing an object named `pbmc`). The Seurat
#' object must contain the requested assay/layer and any metadata columns named
#' by `cell_column`, `condition_column`, or `batch_column`. Condition and batch
#' are optional; set the corresponding argument to `NULL` to omit it. When
#' present, both are encoded in one nuisance design matrix shared by Stage 1
#' and Stage 2.
#' The workflow searches a target-specific adjustment set, fits a sparse Stage 1
#' model by MCMC (default) or variational inference, removes unsupported
#' non-target TFs, then fits Stage 2 exclusively by MCMC. `output`
#' follows TF_activity_computation(): FALSE returns results only in memory, TRUE
#' creates `<TF>_<GENE>_regulatory_output`, and a character value is a custom path.
TF_regulatory_direction_computation <- function(
  target_tf = NULL,
  target_gene = NULL,
  input = NULL,
  seurat_obj = NULL,
  data_file = NULL,
  output = FALSE,
  cell_column = NULL,
  cell_levels = NULL,
  batch_subset = NULL,
  batch_column = NULL,
  condition_column = "sample",
  control_level = "Normal",
  disease_level = "AAA",
  assay = "RNA",
  layer = "data",
  Y_exp = NULL,
  libsize = NULL,
  network_edge_file = .tfregact_default_network_file(),
  confounder_confidence_threshold = 3,
  adjustment_search_starts = 8L,
  adjustment_search_cores = 4L,
  max_adjustment_sets = NULL,
  dagitty_beta = 2,
  gamma = 1,
  eta = 0.5,
  r_dir = 3,
  nuisance_prior_scale = 1,
  confidence_min = 1,
  confidence_max = 10,
  stage1_inference = c("mcmc", "variational"),
  stage1_chains = 3L,
  stage1_iter_warmup = 600L,
  stage1_iter_sampling = 1200L,
  stage1_adapt_delta = 0.95,
  stage1_max_treedepth = 12L,
  stage1_variational_algorithm = c("meanfield", "fullrank"),
  stage1_variational_iter = 10000L,
  stage1_variational_output_samples = 2000L,
  stage1_filter_interval = c(5, 95),
  correlation_filter = TRUE,
  correlation_threshold = 0.7,
  stage2_model = c("no_interaction", "target_interaction"),
  direction_effect = 0.2,
  beta_sd_floor = 0.5,
  stage1_sd_multiplier = 1.5,
  target_interaction_sd = 0.5,
  stage2_chains = 3L,
  stage2_iter_warmup = 500L,
  stage2_iter_sampling = 900L,
  stage2_adapt_delta = 0.95,
  stage2_max_treedepth = 12L,
  compute_loo = FALSE,
  seed = 123L,
  refresh = 50L,
  force_recompile = FALSE
) {
  stage2_model <- match.arg(stage2_model)
  stage1_inference <- match.arg(stage1_inference)
  stage1_variational_algorithm <- match.arg(stage1_variational_algorithm)
  if (!is.null(max_adjustment_sets)) {
    warning(
      "`max_adjustment_sets` is obsolete and ignored; use `adjustment_search_starts`.",
      call. = FALSE
    )
  }
  .tf_regulatory_direction_load_dependencies()
  if (!is.null(input)) {
    if (!is.null(seurat_obj) || !is.null(data_file)) .tf_regulatory_direction_stop("When `input` is supplied, do not also supply `seurat_obj` or `data_file`.")
    input_values <- tf_computation_input_values(input)
    seurat_obj <- input_values$seurat_obj
    if (!is.null(target_tf) && !identical(toupper(trimws(as.character(target_tf[[1]]))), toupper(input_values$target_tf))) .tf_regulatory_direction_stop("`target_tf` disagrees with `input@target_tf`.")
    if (!is.null(target_gene) && !is.na(input_values$target_gene) && !identical(toupper(trimws(as.character(target_gene[[1]]))), toupper(input_values$target_gene))) .tf_regulatory_direction_stop("`target_gene` disagrees with `input@target_gene`.")
    target_tf <- input_values$target_tf; target_gene <- input_values$target_gene
    condition_column <- input_values$condition_column; control_level <- input_values$control_level
    disease_level <- input_values$disease_level; assay <- input_values$assay; layer <- input_values$layer
    Y_exp <- input_values$Y_exp; libsize <- input_values$libsize
    network_edge_file <- input_values$network_edge_file
  }
  target_tf <- if (is.null(target_tf) || !length(target_tf)) NA_character_ else trimws(as.character(target_tf[[1]]))
  target_gene <- if (is.null(target_gene) || !length(target_gene)) NA_character_ else trimws(as.character(target_gene[[1]]))
  if (is.na(target_tf) || !nzchar(target_tf) || is.na(target_gene) || !nzchar(target_gene)) .tf_regulatory_direction_stop("Supply `target_tf` and `target_gene`, or a `TFComputationInput` object containing both.")

  output_dir <- .tf_regulatory_direction_output_dir(output, target_tf, target_gene)
  persist_output <- !is.null(output_dir)
  work_dir <- if (persist_output) output_dir else tempfile(paste0(toupper(target_tf), "_", toupper(target_gene), "_regulatory_"))
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  if (!persist_output) on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)

  pbmc <- .tf_regulatory_direction_load_seurat(seurat_obj, data_file)
  if (!is.null(cell_levels)) {
    if (is.null(cell_column) || !(cell_column %in% colnames(pbmc[[]]))) .tf_regulatory_direction_stop("`cell_levels` requires a valid `cell_column`.")
    pbmc <- subset(pbmc, cells = colnames(pbmc)[pbmc[[cell_column]][, 1] %in% cell_levels])
  }
  if (!is.null(batch_subset)) {
    if (is.null(batch_column) || !(batch_column %in% colnames(pbmc[[]]))) .tf_regulatory_direction_stop("`batch_subset` requires a valid `batch_column`.")
    pbmc <- subset(pbmc, cells = colnames(pbmc)[pbmc[[batch_column]][, 1] %in% batch_subset])
  }
  if (!ncol(pbmc)) .tf_regulatory_direction_stop("No cells remain after subsetting.")

  adjustment_started <- Sys.time()
  adjustment <- find_adjustment_dagitty(
    tf = target_tf, gene = target_gene, edge_file = network_edge_file, outdir = work_dir,
    beta = dagitty_beta, overall_confidence_threshold = confounder_confidence_threshold,
    search_starts = adjustment_search_starts,
    search_cores = adjustment_search_cores,
    search_seed = seed,
    write_full_outputs = persist_output,
    write_files = persist_output
  )
  adjustment_elapsed <- as.numeric(difftime(Sys.time(), adjustment_started, units = "secs"))

  analysis_object <- build_TF_analysis_object_from_seurat(
    seurat_obj = pbmc, adjustment_input = adjustment, target = c(target_tf, target_gene),
    batch = batch_column, sample = condition_column, Y_exp = Y_exp, assay = assay,
    layer = layer, libsize = libsize, confidence_col = "overall_avg_confidence",
    direction_col = "effect_on_gene_A"
  )

  stage1_started <- Sys.time()
  stage1_fit <- run_TF_directional_model(
    analysis_object = analysis_object, stan_file = NULL,
    inference = stage1_inference, gamma = gamma, eta = eta,
    r_dir = r_dir, confidence_min = confidence_min, confidence_max = confidence_max,
    chains = stage1_chains, parallel_chains = stage1_chains,
    iter_warmup = stage1_iter_warmup, iter_sampling = stage1_iter_sampling,
    seed = seed, refresh = refresh, force_recompile = force_recompile,
    control_level = control_level, disease_level = disease_level,
    nuisance_prior_scale = nuisance_prior_scale,
    adapt_delta = stage1_adapt_delta, max_treedepth = stage1_max_treedepth,
    variational_algorithm = stage1_variational_algorithm,
    variational_iter = stage1_variational_iter,
    variational_output_samples = stage1_variational_output_samples
  )
  stage1_elapsed <- as.numeric(difftime(Sys.time(), stage1_started, units = "secs"))

  filtered_analysis_object <- filter_TF_analysis_object_by_beta_ci(
    fit_result = stage1_fit, analysis_object = analysis_object,
    interval = stage1_filter_interval, force_keep_target_tf = TRUE,
    correlation_filter = correlation_filter, correlation_threshold = correlation_threshold,
    plot_correlation_heatmap = persist_output, correlation_heatmap_dir = work_dir
  )

  stage2_started <- Sys.time()
  stage2_fit <- if (identical(stage2_model, "no_interaction")) {
    run_TF_stage2_directional_model(
      analysis_object = filtered_analysis_object, stage1_fit_result = stage1_fit,
      stan_file = NULL, direction_effect = direction_effect, beta_sd_floor = beta_sd_floor,
      stage1_sd_multiplier = stage1_sd_multiplier, control_level = control_level,
      disease_level = disease_level, confidence_min = confidence_min,
      confidence_max = confidence_max, chains = stage2_chains,
      parallel_chains = stage2_chains, iter_warmup = stage2_iter_warmup,
      iter_sampling = stage2_iter_sampling, seed = seed + 1L, refresh = refresh,
      compute_loo = compute_loo, force_recompile = force_recompile,
      nuisance_prior_scale = nuisance_prior_scale,
      adapt_delta = stage2_adapt_delta, max_treedepth = stage2_max_treedepth
    )
  } else {
    run_TF_stage2_target_interaction_model(
      analysis_object = filtered_analysis_object, stage1_fit_result = stage1_fit,
      stan_file = NULL, direction_effect = direction_effect, beta_sd_floor = beta_sd_floor,
      stage1_sd_multiplier = stage1_sd_multiplier, target_interaction_sd = target_interaction_sd,
      control_level = control_level, disease_level = disease_level,
      confidence_min = confidence_min, confidence_max = confidence_max,
      chains = stage2_chains, parallel_chains = stage2_chains,
      iter_warmup = stage2_iter_warmup, iter_sampling = stage2_iter_sampling,
      seed = seed + 1L, refresh = refresh, compute_loo = compute_loo,
      force_recompile = force_recompile, adapt_delta = stage2_adapt_delta,
      nuisance_prior_scale = nuisance_prior_scale,
      max_treedepth = stage2_max_treedepth
    )
  }
  stage2_elapsed <- as.numeric(difftime(Sys.time(), stage2_started, units = "secs"))

  beta_draws <- tf_get_beta_draws_matrix(stage2_fit)
  target_index <- match(toupper(target_tf), toupper(colnames(beta_draws)))
  if (is.na(target_index)) .tf_regulatory_direction_stop("Target TF `%s` is missing from Stage 2 beta draws.", target_tf)
  target_beta <- as.numeric(beta_draws[, target_index])
  direction_summary <- data.frame(
    target_tf = target_tf, target_gene = target_gene,
    stage1_inference = stage1_inference, stage2_model = stage2_model,
    posterior_mean = mean(target_beta), posterior_median = stats::median(target_beta),
    q05 = unname(stats::quantile(target_beta, 0.05)),
    q95 = unname(stats::quantile(target_beta, 0.95)),
    probability_positive = mean(target_beta > 0),
    probability_negative = mean(target_beta < 0),
    inferred_direction = if (mean(target_beta > 0) >= 0.95) "activation" else if (mean(target_beta < 0) >= 0.95) "repression" else "uncertain",
    stringsAsFactors = FALSE
  )
  posterior_plot <- NULL
  timing <- data.frame(
    stage = c("adjustment_search", "stage1", "stage2"),
    elapsed_seconds = c(adjustment_elapsed, stage1_elapsed, stage2_elapsed),
    elapsed_minutes = c(adjustment_elapsed, stage1_elapsed, stage2_elapsed) / 60
  )

  if (persist_output) {
    utils::write.csv(direction_summary, file.path(work_dir, paste0(target_tf, "_", target_gene, "_direction_summary.csv")), row.names = FALSE)
    utils::write.csv(timing, file.path(work_dir, paste0(target_tf, "_", target_gene, "_timing.csv")), row.names = FALSE)
    saveRDS(list(adjustment = adjustment, analysis_object = analysis_object, filtered_analysis_object = filtered_analysis_object, stage1_fit = stage1_fit, stage1_inference = stage1_inference, stage2_fit = stage2_fit, direction_summary = direction_summary, timing = timing), file.path(work_dir, paste0(target_tf, "_", target_gene, "_regulatory_result.rds")))
  }

  invisible(structure(list(
    target_tf = target_tf, target_gene = target_gene,
    output_dir = if (persist_output) normalizePath(work_dir, winslash = "/", mustWork = TRUE) else NULL,
    adjustment = adjustment, analysis_object = analysis_object,
    filtered_analysis_object = filtered_analysis_object, stage1_fit = stage1_fit,
    stage2_fit = stage2_fit, stage1_inference = stage1_inference,
    direction_summary = direction_summary,
    posterior_plot = posterior_plot, timing = timing
  ), class = c("TFRegulatoryDirectionComputation", "list")))
}
