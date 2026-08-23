#!/usr/bin/env Rscript

tf_stage1_input_stop <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

.tf_stage1_input_script_dir <- local({
  source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  source_dir <- if (is.null(source_file) || !nzchar(as.character(source_file[[1]]))) {
    getwd()
  } else {
    dirname(normalizePath(source_file[[1]], winslash = "/", mustWork = FALSE))
  }
  function() source_dir
})

#' Create one condition-interaction Stan data bundle from one TF-target model.
prepare_TF_stage1_screening_stan_data <- function(
  analysis_object,
  target_confidence,
  target_direction,
  beta_prior_scale = 1,
  eta = 0.5,
  r_dir = 3,
  batch_prior_scale = 1,
  target_interaction_sd = 0.5,
  control_level = NULL,
  disease_level = NULL,
  predictor_sd_min = 1e-8,
  standardize_predictors = TRUE
) {
  expr <- as.matrix(analysis_object$expr)
  storage.mode(expr) <- "double"
  feature_names <- rownames(expr)
  target_tf <- trimws(as.character(analysis_object$target[[1]]))
  target_gene <- trimws(as.character(analysis_object$target[[2]]))
  target_match <- which(toupper(trimws(feature_names)) == toupper(target_tf))
  if (length(target_match) != 1L) {
    tf_stage1_input_stop(
      "Target TF `%s` must occur exactly once for target gene `%s`; found %d.",
      target_tf, target_gene, length(target_match)
    )
  }

  confidence <- as.numeric(analysis_object$confidence)
  confidence[target_match] <- as.numeric(target_confidence)
  if (any(!is.finite(confidence)) || any(confidence <= 0)) {
    tf_stage1_input_stop("All predictor confidence values must be finite and positive.")
  }
  direction <- tf_direction_to_int(analysis_object$direction)
  target_direction <- suppressWarnings(as.integer(target_direction[[1]]))
  if (length(direction) != nrow(expr) ||
      is.na(target_direction) || !(target_direction %in% c(-1L, 0L, 1L))) {
    tf_stage1_input_stop("Every edge direction must be -1, 0, or 1.")
  }
  direction[target_match] <- target_direction

  predictor_center <- rowMeans(expr)
  predictor_scale <- apply(expr, 1L, stats::sd)
  estimable <- is.finite(predictor_scale) & predictor_scale > predictor_sd_min
  if (!estimable[target_match]) {
    tf_stage1_input_stop(
      "Target TF `%s` has zero or near-zero variation for target gene `%s`.",
      target_tf, target_gene
    )
  }
  fitted_index <- which(estimable)
  fitted_expr <- expr[fitted_index, , drop = FALSE]
  if (isTRUE(standardize_predictors)) {
    fitted_expr <- sweep(fitted_expr, 1L, predictor_center[fitted_index], "-")
    fitted_expr <- sweep(fitted_expr, 1L, predictor_scale[fitted_index], "/")
  }

  use_batch <- !is.null(analysis_object$batch)
  batch_factor <- if (use_batch) {
    batch_values <- trimws(as.character(analysis_object$batch))
    if (any(is.na(batch_values)) || any(!nzchar(batch_values))) {
      tf_stage1_input_stop("Batch contains missing or empty labels.")
    }
    droplevels(as.factor(batch_values))
  } else {
    NULL
  }
  condition <- tf_condition_from_sample(
    analysis_object = analysis_object,
    control_level = control_level,
    disease_level = disease_level
  )
  if (is.null(condition)) {
    tf_stage1_input_stop("A two-level sample factor is required for condition interaction.")
  }
  libsize <- as.numeric(analysis_object$libsize)
  if (length(libsize) != ncol(expr) || any(!is.finite(libsize)) || any(libsize <= 0)) {
    tf_stage1_input_stop("Library sizes must be finite, positive, and cell-aligned.")
  }

  stan_data <- list(
    N = ncol(expr),
    P = length(fitted_index),
    Y = as.array(as.integer(analysis_object$Y_exp)),
    X = t(fitted_expr),
    condition = as.array(as.integer(condition)),
    log_offset = as.vector(log(libsize / mean(libsize))),
    confidence = as.vector(confidence[fitted_index]),
    direction = as.array(direction[fitted_index]),
    beta_prior_scale = as.numeric(beta_prior_scale),
    eta = as.numeric(eta),
    r_dir = as.numeric(r_dir),
    target_tf_index = as.integer(match(target_match, fitted_index)),
    target_interaction_sd = as.numeric(target_interaction_sd)
  )
  if (use_batch) {
    stan_data$K_batch <- nlevels(batch_factor)
    stan_data$batch <- as.array(as.integer(batch_factor))
    stan_data$batch_prior_scale <- as.numeric(batch_prior_scale)
  }
  list(
    status = "ready",
    target_tf = target_tf,
    target_gene = target_gene,
    target_tf_index = match(target_match, fitted_index),
    feature_names = feature_names[fitted_index],
    cell_names = colnames(expr),
    use_batch = use_batch,
    batch_levels = if (use_batch) levels(batch_factor) else character(0),
    control_level = attr(condition, "control_level"),
    disease_level = attr(condition, "disease_level"),
    predictor_center = stats::setNames(predictor_center[fitted_index], feature_names[fitted_index]),
    predictor_scale = stats::setNames(predictor_scale[fitted_index], feature_names[fitted_index]),
    target_tf_expression = as.numeric(expr[target_match, ]),
    removed_zero_variance_tfs = feature_names[!estimable],
    stan_data = stan_data
  )
}

tf_stage1_input_query_one_adjustment <- function(target_gene, target_tf, config) {
  tryCatch({
    result <- find_adjustment_dagitty(
      tf = target_tf,
      gene = target_gene,
      network = list(dagitty_bundle = .tf_stage1_input_worker_network_bundle),
      edge_file = config$edge_file,
      outdir = config$outdir,
      beta = config$beta,
      overall_confidence_threshold = config$confounder_confidence_threshold,
      max_adjustment_sets = config$max_adjustment_sets,
      include_candidate_confounders = FALSE,
      write_full_outputs = FALSE,
      write_files = FALSE
    )
    list(
      status = "ok",
      target_gene = target_gene,
      filtered_recommended_adjustment_metrics = result$filtered_recommended_adjustment_metrics,
      recommended_adjustment_set = result$recommended_adjustment_set,
      adjustment_sets_evaluated = nrow(result$adjustment_set_scores),
      selected_adjustment_set_size = nrow(result$recommended_adjustment_set),
      error = NA_character_
    )
  }, error = function(e) {
    list(
      status = "error",
      target_gene = target_gene,
      filtered_recommended_adjustment_metrics = NULL,
      recommended_adjustment_set = NULL,
      adjustment_sets_evaluated = NA_integer_,
      selected_adjustment_set_size = NA_integer_,
      error = conditionMessage(e)
    )
  })
}

#' Find the best direct-confounder adjustment set for each direct target gene.
query_TF_target_adjustment_sets <- function(
  target_tf,
  target_genes,
  edge_file = "tf_union_output/TF_Full_Map.RData",
  outdir = "adjustment_output",
  beta = 2,
  confounder_confidence_threshold = 4,
  max_adjustment_sets = 100L,
  cores = 4L,
  checkpoint_count = 2L,
  output_file = NULL,
  resume = TRUE
) {
  target_genes <- unique(trimws(as.character(target_genes)))
  target_genes <- target_genes[nzchar(target_genes)]
  if (!length(target_genes)) {
    tf_stage1_input_stop("`target_genes` must contain at least one gene.")
  }
  cores <- as.integer(cores[[1]])
  max_adjustment_sets <- as.integer(max_adjustment_sets[[1]])
  checkpoint_count <- as.integer(checkpoint_count[[1]])
  if (is.na(cores) || cores < 1L || is.na(max_adjustment_sets) || max_adjustment_sets < 1L ||
      is.na(checkpoint_count) || checkpoint_count < 0L) {
    tf_stage1_input_stop("`cores` and `max_adjustment_sets` must be positive; `checkpoint_count` must be nonnegative.")
  }
  config <- list(
    edge_file = edge_file,
    outdir = outdir,
    beta = as.numeric(beta),
    confounder_confidence_threshold = as.numeric(confounder_confidence_threshold),
    max_adjustment_sets = max_adjustment_sets
  )
  results <- stats::setNames(vector("list", length(target_genes)), target_genes)
  if (!is.null(output_file)) {
    output_file <- normalizePath(output_file, winslash = "/", mustWork = FALSE)
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    if (isTRUE(resume) && file.exists(output_file)) {
      cached <- readRDS(output_file)
      cached_names <- intersect(names(cached), target_genes)
      cached_names <- cached_names[vapply(cached[cached_names], function(x) {
        is.list(x) && identical(x$status, "ok")
      }, logical(1))]
      results[cached_names] <- cached[cached_names]
    }
  }
  pending <- names(results)[vapply(results, is.null, logical(1))]
  if (length(pending)) {
    worker_count <- min(cores, length(pending))
    if (worker_count > 1L) {
      cluster <- parallel::makeCluster(worker_count)
      on.exit(parallel::stopCluster(cluster), add = TRUE)
      workspace_dir <- getwd()
      tfregact_library <- dirname(find.package("TFRegAct"))
      parallel::clusterExport(
        cluster, c("workspace_dir", "config", "tfregact_library"), envir = environment()
      )
      parallel::clusterEvalQ(cluster, {
        setwd(workspace_dir)
        .libPaths(c(tfregact_library, .libPaths()))
        library(TFRegAct)
        .tf_stage1_input_worker_network_bundle <- TFRegAct:::find_adjustment_dagitty_extract_bundle(
          edge_file = config$edge_file
        )
        NULL
      })
      batches <- split(pending, ceiling(seq_along(pending) / worker_count))
      checkpoint_batches <- if (checkpoint_count == 0L) {
        integer(0)
      } else {
        unique(pmax(1L, pmin(length(batches), as.integer(round(seq(
          length(batches) / checkpoint_count,
          length(batches),
          length.out = checkpoint_count
        ))))))
      }
      for (batch_index in seq_along(batches)) {
        genes <- batches[[batch_index]]
        fitted <- parallel::parLapplyLB(
          cluster, genes, tf_stage1_input_query_one_adjustment,
          target_tf = target_tf, config = config
        )
        names(fitted) <- genes
        results[genes] <- fitted
        if (!is.null(output_file) && batch_index %in% checkpoint_batches) {
          saveRDS(results, output_file)
        }
      }
    } else {
      .tf_stage1_input_worker_network_bundle <- find_adjustment_dagitty_extract_bundle(
        edge_file = edge_file
      )
      on.exit(rm(.tf_stage1_input_worker_network_bundle, inherits = FALSE), add = TRUE)
      checkpoint_genes <- if (checkpoint_count == 0L) {
        integer(0)
      } else {
        unique(pmax(1L, pmin(length(pending), as.integer(round(seq(
          length(pending) / checkpoint_count,
          length(pending),
          length.out = checkpoint_count
        ))))))
      }
      for (gene_index in seq_along(pending)) {
        gene <- pending[[gene_index]]
        results[[gene]] <- tf_stage1_input_query_one_adjustment(gene, target_tf, config)
        if (!is.null(output_file) && gene_index %in% checkpoint_genes) {
          saveRDS(results, output_file)
        }
      }
    }
  }
  attr(results, "target_tf") <- target_tf
  attr(results, "confounder_confidence_threshold") <- config$confounder_confidence_threshold
  attr(results, "max_adjustment_sets") <- config$max_adjustment_sets
  attr(results, "direct_confounders_only") <- TRUE
  if (!is.null(output_file)) saveRDS(results, output_file)
  results
}

#' Build the complete nested Stage 1 input from direct TF-gene edges and Seurat.
build_TF_stage1_screening_input <- function(
  seurat_obj,
  direct_target_edges,
  adjustment_results,
  batch_column = NULL,
  condition_column,
  control_level,
  disease_level,
  assay = "RNA",
  layer = "data",
  Y_exp = NULL,
  libsize = NULL,
  confidence_col = "overall_avg_confidence",
  direction_col = "effect_on_gene_A",
  beta_prior_scale = 1,
  eta = 0.5,
  r_dir = 3,
  batch_prior_scale = 1,
  target_interaction_sd = 0.5,
  predictor_sd_min = 1e-8,
  standardize_predictors = TRUE,
  output_file = NULL,
  summary_file = NULL
) {
  required_edge_columns <- c("target_tf", "target_gene", "effect", "direction", "final_confidence")
  if (!is.data.frame(direct_target_edges) ||
      !all(required_edge_columns %in% names(direct_target_edges))) {
    tf_stage1_input_stop(
      "`direct_target_edges` must contain: %s.",
      paste(required_edge_columns, collapse = ", ")
    )
  }
  target_tf <- unique(as.character(direct_target_edges$target_tf))
  if (length(target_tf) != 1L || !nzchar(target_tf)) {
    tf_stage1_input_stop("Direct target edges must contain exactly one target TF.")
  }
  target_genes <- as.character(direct_target_edges$target_gene)
  if (anyDuplicated(target_genes) || !all(target_genes %in% names(adjustment_results))) {
    tf_stage1_input_stop("Every direct target gene must have one adjustment-set result.")
  }

  result <- stats::setNames(vector("list", length(target_genes)), target_genes)
  for (i in seq_along(target_genes)) {
    target_gene <- target_genes[[i]]
    target_edge <- direct_target_edges[i, , drop = FALSE]
    adjustment <- adjustment_results[[target_gene]]
    result[[target_gene]] <- tryCatch({
      if (!identical(adjustment$status, "ok")) stop(adjustment$error)
      analysis_object <- build_TF_analysis_object_from_seurat(
        seurat_obj = seurat_obj,
        adjustment_input = adjustment,
        target = c(target_tf, target_gene),
        batch = batch_column,
        sample = condition_column,
        Y_exp = Y_exp,
        assay = assay,
        layer = layer,
        libsize = libsize,
        confidence_col = confidence_col,
        direction_col = direction_col
      )
      prepared <- prepare_TF_stage1_screening_stan_data(
        analysis_object = analysis_object,
        target_confidence = target_edge$final_confidence[[1]],
        target_direction = target_edge$direction[[1]],
        beta_prior_scale = beta_prior_scale,
        eta = eta,
        r_dir = r_dir,
        batch_prior_scale = batch_prior_scale,
        target_interaction_sd = target_interaction_sd,
        control_level = control_level,
        disease_level = disease_level,
        predictor_sd_min = predictor_sd_min,
        standardize_predictors = standardize_predictors
      )
      metrics <- adjustment$filtered_recommended_adjustment_metrics
      metric_index <- match(toupper(prepared$feature_names), toupper(metrics$tf))
      if (anyNA(metric_index)) stop("Failed to map fitted predictors to adjustment-edge metrics.")
      edge_effect <- metrics$effect_on_gene_A[metric_index]
      edge_effect[prepared$target_tf_index] <- target_edge$effect[[1]]
      prepared$target_edge <- target_edge
      prepared$predictor_edges <- data.frame(
        tf = prepared$feature_names,
        target_gene = target_gene,
        effect = edge_effect,
        direction = as.integer(prepared$stan_data$direction),
        confidence = as.numeric(prepared$stan_data$confidence),
        distance_to_gene = metrics$distance_to_gene_A[metric_index],
        is_target_tf = seq_along(prepared$feature_names) == prepared$target_tf_index,
        stringsAsFactors = FALSE
      )
      prepared$adjustment_tf_count <- nrow(metrics)
      prepared$adjustment_sets_evaluated <- adjustment$adjustment_sets_evaluated
      prepared$selected_adjustment_set_size <- adjustment$selected_adjustment_set_size
      prepared
    }, error = function(e) {
      list(
        status = "error", target_tf = target_tf, target_gene = target_gene,
        target_tf_index = NA_integer_, feature_names = character(0),
        cell_names = character(0), batch_levels = character(0),
        use_batch = !is.null(batch_column),
        control_level = control_level, disease_level = disease_level,
        predictor_center = numeric(0), predictor_scale = numeric(0),
        target_tf_expression = numeric(0), removed_zero_variance_tfs = character(0),
        stan_data = NULL, target_edge = target_edge, predictor_edges = data.frame(),
        adjustment_tf_count = NA_integer_, adjustment_sets_evaluated = NA_integer_,
        selected_adjustment_set_size = NA_integer_, error = conditionMessage(e)
      )
    })
  }
  class(result) <- c("TFStage1ScreeningInputList", "list")
  attr(result, "target_tf") <- target_tf
  attr(result, "model_version") <- "stage1_normal_condition_interaction_stage2_ready_v2"
  attr(result, "prior_family") <- "normal"
  attr(result, "condition_interaction") <- TRUE
  attr(result, "stage2_interface_ready") <- TRUE
  attr(result, "directional_prior") <- TRUE
  attr(result, "edge_direction_included") <- TRUE
  attr(result, "direct_confounders_only") <- TRUE
  attr(result, "cell_count") <- ncol(seurat_obj)
  attr(result, "batch_column") <- batch_column
  attr(result, "use_batch") <- !is.null(batch_column)
  attr(result, "condition_column") <- condition_column
  attr(result, "confidence_threshold") <- attr(direct_target_edges, "query_summary")$confidence_threshold
  attr(result, "confounder_confidence_threshold") <- attr(adjustment_results, "confounder_confidence_threshold")
  attr(result, "max_adjustment_sets") <- attr(adjustment_results, "max_adjustment_sets")
  attr(result, "r_dir") <- as.numeric(r_dir)
  attr(result, "target_interaction_sd") <- as.numeric(target_interaction_sd)
  attr(result, "stan_file") <- if (is.null(batch_column)) {
    "TF_bayesian_screening_stage1_nb_model_nobatch.stan"
  } else {
    "TF_bayesian_screening_stage1_nb_model.stan"
  }

  summary <- do.call(rbind, lapply(result, function(entry) {
    data.frame(
      target_tf = entry$target_tf,
      target_gene = entry$target_gene,
      status = entry$status,
      cells = if (is.null(entry$stan_data)) NA_integer_ else entry$stan_data$N,
      predictors = if (is.null(entry$stan_data)) NA_integer_ else entry$stan_data$P,
      confounders = if (is.null(entry$stan_data)) NA_integer_ else entry$stan_data$P - 1L,
      target_confidence = entry$target_edge$final_confidence[[1]],
      target_direction = entry$target_edge$direction[[1]],
      adjustment_sets_evaluated = entry$adjustment_sets_evaluated,
      selected_adjustment_set_size = entry$selected_adjustment_set_size,
      error = if (is.null(entry$error)) NA_character_ else entry$error,
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  attr(result, "input_summary") <- summary
  if (!is.null(output_file)) {
    output_file <- normalizePath(output_file, winslash = "/", mustWork = FALSE)
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(result, output_file)
  }
  if (!is.null(summary_file)) {
    summary_file <- normalizePath(summary_file, winslash = "/", mustWork = FALSE)
    dir.create(dirname(summary_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(summary, summary_file, row.names = FALSE)
  }
  result
}
