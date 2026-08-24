#!/usr/bin/env Rscript

# Complete three-stage TF activity workflow. Source this file, then call
# TF_activity_computation(). Nothing is executed on source().

.tf_activity_computation_dir <- local({
  source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  source_dir <- if (is.null(source_file) || !nzchar(as.character(source_file[[1]]))) getwd() else dirname(normalizePath(source_file[[1]], winslash = "/", mustWork = FALSE))
  function() source_dir
})

.tf_activity_computation_stop <- function(...) stop(sprintf(...), call. = FALSE)

.tf_activity_computation_load_dependencies <- function() invisible(NULL)

.tf_activity_computation_output_dir <- function(output, target_tf) {
  if (isFALSE(output) || is.null(output)) return(NULL)
  if (isTRUE(output)) return(file.path(getwd(), paste0(toupper(target_tf), "_output")))
  if (is.character(output) && length(output) == 1L && nzchar(output)) return(normalizePath(output, winslash = "/", mustWork = FALSE))
  .tf_activity_computation_stop("`output` must be FALSE, TRUE, or a single non-empty output path.")
}

.tf_activity_computation_load_seurat <- function(seurat_obj, data_file) {
  if (!is.null(seurat_obj) && !is.null(data_file)) {
    .tf_activity_computation_stop("Supply exactly one of `seurat_obj` or `data_file`, not both.")
  }
  if (!is.null(seurat_obj)) return(seurat_obj)
  if (!is.character(data_file) || length(data_file) != 1L || !file.exists(data_file)) .tf_activity_computation_stop("Supply `seurat_obj`, or a valid `data_file` containing an object named `pbmc`.")
  loaded <- new.env(parent = emptyenv())
  load(data_file, envir = loaded)
  if (!exists("pbmc", envir = loaded, inherits = FALSE)) .tf_activity_computation_stop("`data_file` must contain a Seurat object named `pbmc`.")
  get("pbmc", envir = loaded, inherits = FALSE)
}

#' Estimate cell-level activity for one target TF.
#'
#' Supply exactly one input dataset: `seurat_obj` (an in-memory Seurat object)
#' or `data_file` (an .RData file containing an object named `pbmc`). The Seurat
#' object must contain the requested assay/layer and metadata columns for
#' `cell_column`, `condition_column`, and, when used, `batch_column`.
#' `output = FALSE` returns all results in memory and removes intermediate files.
#' `output = TRUE` writes `<TF>_output` below the current working directory.
#' A character `output` is used as the output directory directly.
TF_activity_computation <- function(
  target_tf = NULL,
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
  cell_limit = NULL,
  cell_limit_seed = 123L,
  cores = 4L,
  seed = 123L,
  network_edge_file = .tfregact_default_network_file(),
  target_confidence_threshold = 4,
  target_confidence_override = NULL,
  confounder_confidence_threshold = 4,
  adjustment_search_starts = 8L,
  max_adjustment_sets = NULL,
  dagitty_beta = 2,
  beta_prior_scale = 1,
  eta = 0.5,
  r_dir = 3,
  batch_prior_scale = 1,
  target_interaction_sd = 0.5,
  predictor_sd_min = 1e-8,
  standardize_predictors = TRUE,
  prescreen_target_interval = 0.90,
  prescreen_confounder_interval = 0.90,
  prescreen_variational_iter = 10000L,
  prescreen_output_samples = 2000L,
  mcmc_target_interval = 0.90,
  mcmc_confounder_interval = 0.90,
  mcmc_chains = 3L,
  mcmc_iter_warmup = 400L,
  mcmc_iter_sampling = 600L,
  mcmc_adapt_delta = 0.95,
  mcmc_max_treedepth = 12L,
  mcmc_draws_for_em = 1500L,
  nuisance_draw_count = 100L,
  em_control = list(),
  resume = TRUE,
  force_refit = FALSE,
  force_recompile = FALSE
) {
  .tf_activity_computation_load_dependencies()
  if (!is.null(max_adjustment_sets)) {
    warning(
      "`max_adjustment_sets` is obsolete and ignored; use `adjustment_search_starts`.",
      call. = FALSE
    )
  }
  if (!is.null(input)) {
    if (!is.null(seurat_obj) || !is.null(data_file)) .tf_activity_computation_stop("When `input` is supplied, do not also supply `seurat_obj` or `data_file`.")
    input_values <- tf_computation_input_values(input)
    seurat_obj <- input_values$seurat_obj
    input_target_tf <- input_values$target_tf
    if (!is.null(target_tf) && !identical(toupper(trimws(as.character(target_tf[[1]]))), toupper(input_target_tf))) .tf_activity_computation_stop("`target_tf` disagrees with `input@target_tf`.")
    target_tf <- input_target_tf
    condition_column <- input_values$condition_column; control_level <- input_values$control_level
    disease_level <- input_values$disease_level; assay <- input_values$assay; layer <- input_values$layer
    network_edge_file <- input_values$network_edge_file
  }
  target_tf <- if (is.null(target_tf) || !length(target_tf)) NA_character_ else trimws(as.character(target_tf[[1]]))
  if (is.na(target_tf) || !nzchar(target_tf)) .tf_activity_computation_stop("Supply `target_tf`, or a `TFComputationInput` object containing `target_tf`.")
  output_dir <- .tf_activity_computation_output_dir(output, target_tf)
  persist_output <- !is.null(output_dir)
  work_dir <- if (persist_output) output_dir else tempfile(paste0(toupper(target_tf), "_activity_"))
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  if (!persist_output) on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)

  pbmc <- .tf_activity_computation_load_seurat(seurat_obj, data_file)
  if (!is.null(cell_levels)) {
    if (is.null(cell_column) || !(cell_column %in% colnames(pbmc[[]]))) .tf_activity_computation_stop("`cell_levels` requires a valid `cell_column`.")
    pbmc <- subset(pbmc, cells = colnames(pbmc)[pbmc[[cell_column]][, 1] %in% cell_levels])
  }
  if (!is.null(batch_subset)) {
    if (is.null(batch_column) || !(batch_column %in% colnames(pbmc[[]]))) .tf_activity_computation_stop("`batch_subset` requires a valid `batch_column`.")
    pbmc <- subset(pbmc, cells = colnames(pbmc)[pbmc[[batch_column]][, 1] %in% batch_subset])
  }
  if (!is.null(cell_limit) && ncol(pbmc) > cell_limit) {
    set.seed(as.integer(cell_limit_seed))
    pbmc <- subset(pbmc, cells = sample(colnames(pbmc), as.integer(cell_limit)))
  }
  if (!ncol(pbmc)) .tf_activity_computation_stop("No cells remain after subsetting.")

  graph_search_started <- Sys.time()
  direct_target_edges <- query_TF_direct_target_genes(target_tf = target_tf, edge_file = network_edge_file, confidence_threshold = target_confidence_threshold, target_confidence_override = target_confidence_override, return_all = FALSE)
  if (!nrow(direct_target_edges)) .tf_activity_computation_stop("No direct target genes passed target confidence threshold %s.", target_confidence_threshold)
  adjustment_results <- query_TF_target_adjustment_sets(target_tf = target_tf, target_genes = direct_target_edges$target_gene, edge_file = network_edge_file, outdir = work_dir, beta = dagitty_beta, confounder_confidence_threshold = confounder_confidence_threshold, search_starts = adjustment_search_starts, search_seed = seed, cores = cores, checkpoint_count = 2L, output_file = if (persist_output) file.path(work_dir, paste0(target_tf, "_adjustment_sets.rds")) else NULL, resume = resume)
  graph_search_elapsed_seconds <- as.numeric(difftime(Sys.time(), graph_search_started, units = "secs"))

  input_build_started <- Sys.time()
  screening_input <- build_TF_stage1_screening_input(seurat_obj = pbmc, direct_target_edges = direct_target_edges, adjustment_results = adjustment_results, batch_column = batch_column, condition_column = condition_column, control_level = control_level, disease_level = disease_level, assay = assay, layer = layer, confidence_col = "overall_avg_confidence", direction_col = "effect_on_gene_A", beta_prior_scale = beta_prior_scale, eta = eta, r_dir = r_dir, batch_prior_scale = batch_prior_scale, target_interaction_sd = target_interaction_sd, predictor_sd_min = predictor_sd_min, standardize_predictors = standardize_predictors, output_file = if (persist_output) file.path(work_dir, paste0(target_tf, "_stage1_screening_input.rds")) else NULL, summary_file = if (persist_output) file.path(work_dir, paste0(target_tf, "_stage1_screening_input_summary.csv")) else NULL)
  input_build_elapsed_seconds <- as.numeric(difftime(Sys.time(), input_build_started, units = "secs"))

  pipeline_started <- Sys.time()
  pipeline_result <- run_TF_three_stage_pipeline(screening_input = screening_input, output_dir = work_dir, cores = cores, seed = seed, prescreen_target_interval = prescreen_target_interval, prescreen_confounder_interval = prescreen_confounder_interval, prescreen_variational_iter = prescreen_variational_iter, prescreen_output_samples = prescreen_output_samples, mcmc_target_interval = mcmc_target_interval, mcmc_confounder_interval = mcmc_confounder_interval, mcmc_chains = mcmc_chains, mcmc_iter_warmup = mcmc_iter_warmup, mcmc_iter_sampling = mcmc_iter_sampling, mcmc_adapt_delta = mcmc_adapt_delta, mcmc_max_treedepth = mcmc_max_treedepth, mcmc_draws_for_em = mcmc_draws_for_em, nuisance_draw_count = nuisance_draw_count, em_control = em_control, resume = resume, force_refit = force_refit, force_recompile = force_recompile)
  pipeline_elapsed_seconds <- as.numeric(difftime(Sys.time(), pipeline_started, units = "secs"))

  activity_column <- paste0(target_tf, "_activity_A")
  activity_index <- match(colnames(pbmc), pipeline_result$em_fit$activity_summary$cell)
  if (anyNA(activity_index)) .tf_activity_computation_stop("Stage 3 activity cells do not match the input Seurat object.")
  activity_values <- pipeline_result$em_fit$activity_summary$activity_mean[activity_index]
  names(activity_values) <- colnames(pbmc)
  pbmc[[activity_column]] <- activity_values
  overall_timing <- rbind(data.frame(stage = c("graph_search", "stage1_input_build"), elapsed_seconds = c(graph_search_elapsed_seconds, input_build_elapsed_seconds), elapsed_minutes = c(graph_search_elapsed_seconds, input_build_elapsed_seconds) / 60), pipeline_result$timing, data.frame(stage = "three_stage_total", elapsed_seconds = pipeline_elapsed_seconds, elapsed_minutes = pipeline_elapsed_seconds / 60))

  if (persist_output) {
    utils::write.csv(direct_target_edges, file.path(work_dir, paste0(target_tf, "_direct_target_edges.csv")), row.names = FALSE)
    saveRDS(pipeline_result, file.path(work_dir, paste0(target_tf, "_three_stage_pipeline_result.rds")))
    utils::write.csv(pipeline_result$em_fit$activity_summary, file.path(work_dir, paste0(target_tf, "_stage3_activity_by_cell.csv")), row.names = FALSE)
    utils::write.csv(pipeline_result$em_fit$beta_summary, file.path(work_dir, paste0(target_tf, "_stage3_target_beta_summary.csv")), row.names = FALSE)
    utils::write.csv(pipeline_result$em_fit$convergence_summary, file.path(work_dir, paste0(target_tf, "_stage3_em_convergence.csv")), row.names = FALSE)
    utils::write.csv(pipeline_result$em_fit$activity_condition_difference_by_draw, file.path(work_dir, paste0(target_tf, "_stage3_activity_difference_by_draw.csv")), row.names = FALSE)
    utils::write.csv(pipeline_result$em_fit$activity_condition_difference_summary, file.path(work_dir, paste0(target_tf, "_stage3_activity_difference_summary.csv")), row.names = FALSE)
    saveRDS(pbmc, file.path(work_dir, paste0(target_tf, "_stage3_activity_seurat.rds")))
    utils::write.csv(overall_timing, file.path(work_dir, paste0(target_tf, "_pipeline_timing.csv")), row.names = FALSE)
  } else {
    pipeline_result$paths <- NULL
  }
  result <- list(target_tf = target_tf, output_dir = if (persist_output) normalizePath(work_dir, winslash = "/", mustWork = TRUE) else NULL, direct_target_edges = direct_target_edges, adjustment_results = adjustment_results, screening_input = screening_input, pipeline_result = pipeline_result, seurat_object = pbmc, timing = overall_timing)
  class(result) <- c("TFActivityComputation", "list")
  invisible(result)
}
