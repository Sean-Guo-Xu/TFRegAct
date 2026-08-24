#!/usr/bin/env Rscript

tf_stage1_screening_stop <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

tf_stage1_screening_require_pkg <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    tf_stage1_screening_stop("Package `%s` is required.", package)
  }
}

tf_stage1_screening_script_dir <- local({
  source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  source_dir <- if (is.null(source_file) || !nzchar(as.character(source_file[[1]]))) {
    getwd()
  } else {
    dirname(normalizePath(source_file[[1]], winslash = "/", mustWork = FALSE))
  }
  function() source_dir
})

tf_stage1_screening_resolve_stan_file <- function(
  screening_input,
  stan_file,
  prior_family
) {
  input_use_batch <- attr(screening_input, "use_batch")
  use_batch <- if (is.null(input_use_batch)) {
    is.list(screening_input[[1]]) &&
      !is.null(screening_input[[1]]$stan_data$K_batch)
  } else {
    isTRUE(input_use_batch)
  }
  if (is.null(stan_file) || !nzchar(as.character(stan_file[[1]]))) {
    stan_file <- if (identical(prior_family, "laplace")) {
      if (use_batch) {
        "TF_bayesian_prescreen_laplace_nb_model.stan"
      } else {
        "TF_bayesian_prescreen_laplace_nb_model_nobatch.stan"
      }
    } else {
      input_stan_file <- attr(screening_input, "stan_file")
      if (is.null(input_stan_file) ||
          !nzchar(as.character(input_stan_file[[1]]))) {
        if (use_batch) {
          "TF_bayesian_screening_stage1_nb_model.stan"
        } else {
          "TF_bayesian_screening_stage1_nb_model_nobatch.stan"
        }
      } else {
        input_stan_file
      }
    }
  }
  stan_file <- as.character(stan_file[[1]])
  .tfregact_stan_file(stan_file)
}

tf_stage1_screening_default_output <- function(screening_input, inference) {
  target_tf <- attr(screening_input, "target_tf")
  if (is.null(target_tf) || !nzchar(as.character(target_tf[[1]]))) {
    target_tf <- "target_TF"
  }
  suffix <- if (identical(inference, "mcmc")) {
    "stage1_screening_stage2_ready_results.rds"
  } else {
    "stage1_screening_stage2_ready_results_variational.rds"
  }
  file.path(
    tf_stage1_screening_script_dir(),
    "adjustment_output",
    sprintf("%s_%s", target_tf[[1]], suffix)
  )
}

tf_stage1_screening_input_signature <- function(screening_input, ready_genes) {
  model_parts <- vapply(ready_genes, function(gene) {
    entry <- screening_input[[gene]]
    paste(
      gene,
      entry$target_tf,
      entry$stan_data$N,
      entry$stan_data$P,
      entry$stan_data$target_tf_index,
      sum(entry$stan_data$condition == 0L),
      sum(entry$stan_data$condition == 1L),
      entry$stan_data$beta_prior_scale,
      entry$stan_data$eta,
      entry$stan_data$r_dir,
      entry$stan_data$target_interaction_sd,
      paste(entry$feature_names, collapse = ","),
      sep = "|"
    )
  }, character(1))
  paste(
    attr(screening_input, "target_tf"),
    attr(screening_input, "cell_count"),
    attr(screening_input, "model_version"),
    paste(model_parts, collapse = ";"),
    sep = "::"
  )
}

tf_stage1_screening_validate_input <- function(screening_input) {
  if (!is.list(screening_input) || !length(screening_input)) {
    tf_stage1_screening_stop("`screening_input` must be a non-empty nested list.")
  }
  if (!isTRUE(attr(screening_input, "recursive_adjustment_search")) &&
      !isTRUE(attr(screening_input, "direct_confounders_only"))) {
    tf_stage1_screening_stop(
      "`screening_input` must contain a validated causal adjustment search."
    )
  }
  if (!identical(
    attr(screening_input, "model_version"),
    "stage1_normal_condition_interaction_stage2_ready_v2"
  )) {
    tf_stage1_screening_stop(
      paste0(
        "`screening_input` is not the Normal-prior condition-interaction ",
        "Stage 1 object prepared for the EM Stage 2 interface. Rebuild the Stage 1 input first."
      )
    )
  }
  ready <- names(screening_input)[vapply(
    screening_input,
    function(entry) {
      is.list(entry) && identical(entry$status, "ready") &&
        is.list(entry$stan_data) && length(entry$feature_names) > 0L
    },
    logical(1)
  )]
  if (!length(ready)) {
    tf_stage1_screening_stop("No `ready` gene models were found in `screening_input`.")
  }
  for (gene in ready) {
    entry <- screening_input[[gene]]
    if (!identical(as.integer(entry$stan_data$P), length(entry$feature_names))) {
      tf_stage1_screening_stop(
        "Feature count and Stan P differ for target gene `%s`.",
        gene
      )
    }
    if (length(entry$target_tf_index) != 1L ||
        is.na(entry$target_tf_index) ||
        entry$target_tf_index < 1L ||
        entry$target_tf_index > entry$stan_data$P) {
      tf_stage1_screening_stop("Invalid target TF index for target gene `%s`.", gene)
    }
    if (length(entry$target_tf_expression) != entry$stan_data$N ||
        any(!is.finite(entry$target_tf_expression)) ||
        any(entry$target_tf_expression < 0)) {
      tf_stage1_screening_stop(
        "Invalid nonnegative target-TF expression for target gene `%s`.",
        gene
      )
    }
    required_stan_fields <- c(
      "condition", "target_tf_index", "target_interaction_sd"
    )
    if (!all(required_stan_fields %in% names(entry$stan_data)) ||
        length(entry$stan_data$condition) != entry$stan_data$N ||
        !all(entry$stan_data$condition %in% c(0L, 1L)) ||
        !identical(
          as.integer(entry$stan_data$target_tf_index),
          as.integer(entry$target_tf_index)
        )) {
      tf_stage1_screening_stop(
        "Invalid condition-interaction Stan data for target gene `%s`.",
        gene
      )
    }
    entry_use_batch <- if (is.null(entry$use_batch)) {
      !is.null(entry$stan_data$K_batch)
    } else {
      isTRUE(entry$use_batch)
    }
    if (entry_use_batch) {
      if (!all(c("K_batch", "batch", "batch_prior_scale") %in% names(entry$stan_data)) ||
          length(entry$stan_data$batch) != entry$stan_data$N ||
          entry$stan_data$K_batch < 1L) {
        tf_stage1_screening_stop("Invalid batch data for target gene `%s`.", gene)
      }
    } else if (any(c("K_batch", "batch", "batch_prior_scale") %in% names(entry$stan_data))) {
      tf_stage1_screening_stop("No-batch input `%s` contains batch Stan data.", gene)
    }
  }
  ready
}

tf_stage1_screening_summarize_fit <- function(
  fit,
  entry,
  inference,
  stage2_draw_count,
  prior_family,
  target_interval,
  confounder_interval
) {
  P <- entry$stan_data$P
  use_batch <- if (is.null(entry$use_batch)) {
    !is.null(entry$stan_data$K_batch)
  } else {
    isTRUE(entry$use_batch)
  }
  K_batch <- if (use_batch) entry$stan_data$K_batch else 0L
  beta_variables <- sprintf("beta[%d]", seq_len(P))
  batch_variables <- if (use_batch) {
    sprintf("batch_effect[%d]", seq_len(K_batch))
  } else {
    character(0)
  }
  target_variables <- c(
    "beta_target_control",
    "beta_target_disease",
    "beta_target_delta"
  )
  scalar_variables <- c("alpha", "condition_effect", "phi")
  requested_variables <- c(
    scalar_variables,
    beta_variables,
    batch_variables,
    target_variables
  )
  joint_draws <- as.matrix(fit$draws(
    variables = c(
      "alpha", "condition_effect", "phi", "beta",
      if (use_batch) "batch_effect" else character(0),
      target_variables
    ),
    format = "draws_matrix"
  ))
  variable_match <- match(requested_variables, colnames(joint_draws))
  if (anyNA(variable_match)) {
    tf_stage1_screening_stop(
      "Missing Stage 2 interface draws for target gene `%s`: %s.",
      entry$target_gene,
      paste(requested_variables[is.na(variable_match)], collapse = ", ")
    )
  }
  joint_draws <- joint_draws[, variable_match, drop = FALSE]
  beta_draws <- joint_draws[, beta_variables, drop = FALSE]
  target_draws <- joint_draws[, target_variables, drop = FALSE]

  quantiles <- t(apply(
    beta_draws,
    2,
    stats::quantile,
    probs = c(0.025, 0.05, 0.5, 0.95, 0.975),
    names = FALSE,
    na.rm = TRUE
  ))
  if (identical(inference, "mcmc")) {
    diagnostic_summary <- fit$summary(variables = "beta")
    diagnostic_match <- match(beta_variables, diagnostic_summary$variable)
    if (anyNA(diagnostic_match)) {
      tf_stage1_screening_stop(
        "Missing beta diagnostics for target gene `%s`.",
        entry$target_gene
      )
    }
    diagnostic_summary <- diagnostic_summary[diagnostic_match, , drop = FALSE]
    beta_rhat <- as.numeric(diagnostic_summary$rhat)
    beta_ess_bulk <- as.numeric(diagnostic_summary$ess_bulk)
    beta_ess_tail <- as.numeric(diagnostic_summary$ess_tail)
  } else {
    beta_rhat <- rep(NA_real_, P)
    beta_ess_bulk <- rep(NA_real_, P)
    beta_ess_tail <- rep(NA_real_, P)
  }

  edge_match <- match(
    toupper(entry$feature_names),
    toupper(entry$predictor_edges$tf)
  )
  if (anyNA(edge_match)) {
    tf_stage1_screening_stop(
      "Failed to match predictor edges for target gene `%s`.",
      entry$target_gene
    )
  }
  edge_table <- entry$predictor_edges[edge_match, , drop = FALSE]
  is_target <- seq_len(P) == entry$target_tf_index
  interval_level <- ifelse(is_target, target_interval, confounder_interval)
  interval_column <- function(level, lower = TRUE) {
    if (isTRUE(all.equal(level, 0.90, tolerance = 1e-12))) {
      if (lower) 2L else 4L
    } else if (isTRUE(all.equal(level, 0.95, tolerance = 1e-12))) {
      if (lower) 1L else 5L
    } else {
      tf_stage1_screening_stop(
        "Only 90%% and 95%% screening intervals are currently supported."
      )
    }
  }
  interval_lower <- vapply(seq_len(P), function(j) {
    quantiles[j, interval_column(interval_level[[j]], lower = TRUE)]
  }, numeric(1))
  interval_upper <- vapply(seq_len(P), function(j) {
    quantiles[j, interval_column(interval_level[[j]], lower = FALSE)]
  }, numeric(1))

  prior_scale <- entry$stan_data$beta_prior_scale *
    (entry$stan_data$confidence / 10) ^ entry$stan_data$eta
  if (identical(prior_family, "laplace")) {
    prior_mean <- ifelse(
      entry$stan_data$direction == 0L,
      0,
      entry$stan_data$direction *
        log((1 + entry$stan_data$r_dir) / 2) * prior_scale
    )
    prior_sd <- sqrt(2) * prior_scale
  } else {
    direction_z <- stats::qnorm(
      entry$stan_data$r_dir / (1 + entry$stan_data$r_dir)
    )
    prior_mean <- ifelse(
      entry$stan_data$direction == 0L,
      0,
      entry$stan_data$direction * direction_z * prior_scale
    )
    prior_sd <- prior_scale
  }
  parameter_table <- data.frame(
    variable = beta_variables,
    beta_index = seq_len(P),
    tf = entry$feature_names,
    role = ifelse(is_target, "target_tf", "confounder"),
    effect = edge_table$effect,
    direction = as.integer(edge_table$direction),
    confidence = as.numeric(edge_table$confidence),
    prior_family = prior_family,
    prior_mean = as.numeric(prior_mean),
    prior_scale = as.numeric(prior_scale),
    prior_sd = as.numeric(prior_sd),
    beta_mean = colMeans(beta_draws),
    beta_sd = apply(beta_draws, 2, stats::sd),
    beta_median = quantiles[, 3],
    beta_q025 = quantiles[, 1],
    beta_q05 = quantiles[, 2],
    beta_q95 = quantiles[, 4],
    beta_q975 = quantiles[, 5],
    prob_beta_positive = colMeans(beta_draws > 0),
    prob_beta_negative = colMeans(beta_draws < 0),
    rhat = beta_rhat,
    ess_bulk = beta_ess_bulk,
    ess_tail = beta_ess_tail,
    screening_interval_level = interval_level,
    screening_interval_lower = interval_lower,
    screening_interval_upper = interval_upper,
    screening_interval_excludes_zero =
      interval_lower > 0 | interval_upper < 0,
    stringsAsFactors = FALSE
  )

  target_quantiles <- t(apply(
    target_draws,
    2,
    stats::quantile,
    probs = c(0.025, 0.05, 0.5, 0.95, 0.975),
    names = FALSE,
    na.rm = TRUE
  ))
  target_prefixes <- c("control", "disease", "delta")
  for (target_index in seq_along(target_prefixes)) {
    prefix <- target_prefixes[[target_index]]
    for (suffix in c("mean", "sd", "median", "q025", "q05", "q95", "q975")) {
      parameter_table[[sprintf("beta_%s_%s", prefix, suffix)]] <- NA_real_
    }
    parameter_table[[sprintf("prob_beta_%s_positive", prefix)]] <- NA_real_
    parameter_table[[sprintf("prob_beta_%s_negative", prefix)]] <- NA_real_
    values <- target_draws[, target_index]
    parameter_table[[sprintf("beta_%s_mean", prefix)]][is_target] <- mean(values)
    parameter_table[[sprintf("beta_%s_sd", prefix)]][is_target] <- stats::sd(values)
    parameter_table[[sprintf("beta_%s_median", prefix)]][is_target] <-
      target_quantiles[target_index, 3]
    parameter_table[[sprintf("beta_%s_q025", prefix)]][is_target] <-
      target_quantiles[target_index, 1]
    parameter_table[[sprintf("beta_%s_q05", prefix)]][is_target] <-
      target_quantiles[target_index, 2]
    parameter_table[[sprintf("beta_%s_q95", prefix)]][is_target] <-
      target_quantiles[target_index, 4]
    parameter_table[[sprintf("beta_%s_q975", prefix)]][is_target] <-
      target_quantiles[target_index, 5]
    parameter_table[[sprintf("prob_beta_%s_positive", prefix)]][is_target] <-
      mean(values > 0)
    parameter_table[[sprintf("prob_beta_%s_negative", prefix)]][is_target] <-
      mean(values < 0)
  }
  parameter_table$control_interval_excludes_zero <- NA
  parameter_table$disease_interval_excludes_zero <- NA
  parameter_table$beta_control_interval_lower <- NA_real_
  parameter_table$beta_control_interval_upper <- NA_real_
  parameter_table$beta_disease_interval_lower <- NA_real_
  parameter_table$beta_disease_interval_upper <- NA_real_
  target_lower_column <- interval_column(target_interval, lower = TRUE)
  target_upper_column <- interval_column(target_interval, lower = FALSE)
  parameter_table$beta_control_interval_lower[is_target] <-
    target_quantiles[1, target_lower_column]
  parameter_table$beta_control_interval_upper[is_target] <-
    target_quantiles[1, target_upper_column]
  parameter_table$beta_disease_interval_lower[is_target] <-
    target_quantiles[2, target_lower_column]
  parameter_table$beta_disease_interval_upper[is_target] <-
    target_quantiles[2, target_upper_column]
  parameter_table$control_interval_excludes_zero[is_target] <-
    target_quantiles[1, target_lower_column] > 0 |
      target_quantiles[1, target_upper_column] < 0
  parameter_table$disease_interval_excludes_zero[is_target] <-
    target_quantiles[2, target_lower_column] > 0 |
      target_quantiles[2, target_upper_column] < 0

  retained_draw_count <- min(stage2_draw_count, nrow(joint_draws))
  retained_draw_index <- unique(as.integer(round(seq(
    1,
    nrow(joint_draws),
    length.out = retained_draw_count
  ))))
  retained_draws <- joint_draws[retained_draw_index, , drop = FALSE]
  retained_beta <- retained_draws[, beta_variables, drop = FALSE]
  colnames(retained_beta) <- entry$feature_names
  retained_batch <- if (use_batch) {
    result <- retained_draws[, batch_variables, drop = FALSE]
    colnames(result) <- entry$batch_levels
    result
  } else {
    matrix(0, nrow = nrow(retained_draws), ncol = 0L)
  }
  stage2_parameter_draws <- list(
    interface_version = "tf_em_stage2_parameter_draws_v1",
    prior_family = prior_family,
    source_draw_index = retained_draw_index,
    draw_count = length(retained_draw_index),
    alpha = as.numeric(retained_draws[, "alpha"]),
    condition_effect = as.numeric(retained_draws[, "condition_effect"]),
    phi = as.numeric(retained_draws[, "phi"]),
    beta = retained_beta,
    batch_effect = retained_batch,
    beta_target_control = as.numeric(
      retained_draws[, "beta_target_control"]
    ),
    beta_target_disease = as.numeric(
      retained_draws[, "beta_target_disease"]
    ),
    beta_target_delta = as.numeric(retained_draws[, "beta_target_delta"]),
    target_tf_index = as.integer(entry$target_tf_index),
    target_interaction_sd = as.numeric(
      entry$stan_data$target_interaction_sd
    ),
    beta_prior_mean = stats::setNames(prior_mean, entry$feature_names),
    beta_prior_sd = stats::setNames(prior_sd, entry$feature_names)
  )

  finite_rhat <- parameter_table$rhat[is.finite(parameter_table$rhat)]
  finite_bulk <- parameter_table$ess_bulk[is.finite(parameter_table$ess_bulk)]
  finite_tail <- parameter_table$ess_tail[is.finite(parameter_table$ess_tail)]
  diagnostics <- list(
    inference = inference,
    posterior_draws = nrow(joint_draws),
    stage2_draws_retained = length(retained_draw_index),
    max_beta_rhat = if (length(finite_rhat)) max(finite_rhat) else NA_real_,
    min_beta_ess_bulk = if (length(finite_bulk)) min(finite_bulk) else NA_real_,
    min_beta_ess_tail = if (length(finite_tail)) min(finite_tail) else NA_real_,
    beta_rhat_over_1_05 = if (identical(inference, "mcmc")) {
      sum(parameter_table$rhat > 1.05, na.rm = TRUE)
    } else {
      NA_integer_
    },
    control_cells = sum(entry$stan_data$condition == 0L),
    disease_cells = sum(entry$stan_data$condition == 1L)
  )

  list(
    target_tf = parameter_table[is_target, , drop = FALSE],
    confounders = parameter_table[!is_target, , drop = FALSE],
    stage2_parameter_draws = stage2_parameter_draws,
    diagnostics = diagnostics
  )
}

tf_stage1_screening_fit_one <- function(job, sampling) {
  started <- Sys.time()
  entry <- job$entry
  tryCatch({
    fit <- if (identical(sampling$inference, "mcmc")) {
      .tf_stage1_worker_model$sample(
        data = entry$stan_data,
        chains = sampling$chains,
        parallel_chains = 1L,
        iter_warmup = sampling$iter_warmup,
        iter_sampling = sampling$iter_sampling,
        seed = sampling$seed + job$model_index + job$retry_seed_offset,
        refresh = sampling$refresh,
        adapt_delta = sampling$adapt_delta,
        max_treedepth = sampling$max_treedepth,
        save_warmup = FALSE
      )
    } else {
      .tf_stage1_worker_model$variational(
        data = entry$stan_data,
        seed = sampling$seed + job$model_index + job$retry_seed_offset,
        refresh = sampling$refresh,
        algorithm = sampling$variational_algorithm,
        iter = sampling$variational_iter,
        output_samples = sampling$variational_output_samples
      )
    }
    summarized <- tf_stage1_screening_summarize_fit(
      fit,
      entry,
      inference = sampling$inference,
      stage2_draw_count = sampling$stage2_draw_count,
      prior_family = sampling$prior_family,
      target_interval = sampling$target_interval,
      confounder_interval = sampling$confounder_interval
    )
    cmdstan_output_files <- fit$output_files()
    result <- list(
      status = "ok",
      inference = sampling$inference,
      target_gene = entry$target_gene,
      target_tf_name = entry$target_tf,
      target_tf = summarized$target_tf,
      confounders = summarized$confounders,
      stage2_parameter_draws = summarized$stage2_parameter_draws,
      control_level = entry$control_level,
      disease_level = entry$disease_level,
      diagnostics = summarized$diagnostics,
      sampling = sampling,
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
      error = NA_character_
    )
    rm(fit)
    gc(verbose = FALSE)

    # Remove only CmdStan files that are confirmed to live inside this worker's
    # temporary directory once compact summaries and joint draws are retained.
    worker_temp_dir <- normalizePath(
      tempdir(),
      winslash = "/",
      mustWork = TRUE
    )
    normalized_output_files <- normalizePath(
      cmdstan_output_files,
      winslash = "/",
      mustWork = FALSE
    )
    safe_temp_files <- startsWith(
      tolower(normalized_output_files),
      paste0(tolower(worker_temp_dir), "/")
    )
    if (any(safe_temp_files)) {
      unlink(normalized_output_files[safe_temp_files], force = TRUE)
    }
    result
  }, error = function(error) {
    list(
      status = "error",
      inference = sampling$inference,
      target_gene = entry$target_gene,
      target_tf_name = entry$target_tf,
      target_tf = data.frame(),
      confounders = data.frame(),
      stage2_parameter_draws = NULL,
      control_level = entry$control_level,
      disease_level = entry$disease_level,
      diagnostics = NULL,
      sampling = sampling,
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
      error = conditionMessage(error)
    )
  })
}

#' Fit all ready Stage 1 screening models using four model-level workers.
#'
#' The first level of the returned list is target gene. Each model entry keeps
#' beta summaries for the target TF and all of its direct confounders, plus a
#' compact set of aligned joint posterior draws for the EM Stage 2 interface.
#' Full CmdStan fit objects and the remaining posterior draws are discarded.
run_TF_stage1_screening_regressions <- function(
  screening_input,
  cores = 4L,
  stan_file = NULL,
  chains = 3L,
  iter_warmup = 600L,
  iter_sampling = 1200L,
  seed = 123L,
  refresh = 0L,
  adapt_delta = 0.95,
  max_treedepth = 12L,
  output_file = NULL,
  checkpoint_file = NULL,
  checkpoint_count = 2L,
  retry_seed_offset = 100000L,
  resume = TRUE,
  retry_errors = TRUE,
  force_refit = FALSE,
  force_recompile = FALSE,
  inference = c("mcmc", "variational"),
  prior_family = c("normal", "laplace"),
  target_interval = 0.90,
  confounder_interval = 0.90,
  variational_algorithm = c("meanfield", "fullrank"),
  variational_iter = 10000L,
  variational_output_samples = 2000L,
  stage2_draw_count = 50L,
  max_genes = NULL
) {
  tf_stage1_screening_require_pkg("cmdstanr")
  ready_genes <- tf_stage1_screening_validate_input(screening_input)
  if (!is.null(max_genes)) {
    max_genes <- suppressWarnings(as.integer(max_genes[[1]]))
    if (is.na(max_genes) || max_genes < 1L) {
      tf_stage1_screening_stop("`max_genes` must be a positive integer or NULL.")
    }
    ready_genes <- utils::head(ready_genes, max_genes)
  }
  skipped_genes <- setdiff(names(screening_input), ready_genes)
  input_signature <- tf_stage1_screening_input_signature(
    screening_input,
    ready_genes
  )

  cores <- suppressWarnings(as.integer(cores[[1]]))
  chains <- suppressWarnings(as.integer(chains[[1]]))
  iter_warmup <- suppressWarnings(as.integer(iter_warmup[[1]]))
  iter_sampling <- suppressWarnings(as.integer(iter_sampling[[1]]))
  seed <- suppressWarnings(as.integer(seed[[1]]))
  refresh <- suppressWarnings(as.integer(refresh[[1]]))
  max_treedepth <- suppressWarnings(as.integer(max_treedepth[[1]]))
  inference <- match.arg(as.character(inference[[1]]), c("mcmc", "variational"))
  prior_family <- match.arg(
    as.character(prior_family[[1]]),
    c("normal", "laplace")
  )
  target_interval <- as.numeric(target_interval[[1]])
  confounder_interval <- as.numeric(confounder_interval[[1]])
  variational_algorithm <- match.arg(
    as.character(variational_algorithm[[1]]),
    c("meanfield", "fullrank")
  )
  variational_iter <- suppressWarnings(as.integer(variational_iter[[1]]))
  variational_output_samples <- suppressWarnings(as.integer(
    variational_output_samples[[1]]
  ))
  stage2_draw_count <- suppressWarnings(as.integer(stage2_draw_count[[1]]))
  checkpoint_count <- suppressWarnings(as.integer(checkpoint_count[[1]]))
  retry_seed_offset <- suppressWarnings(as.integer(retry_seed_offset[[1]]))
  if (is.na(cores) || cores < 1L || is.na(chains) ||
      (identical(inference, "mcmc") && chains < 2L) ||
      is.na(iter_warmup) || iter_warmup < 0L ||
      is.na(iter_sampling) || iter_sampling < 1L ||
      is.na(seed) || is.na(refresh) || refresh < 0L ||
      !is.finite(adapt_delta) || adapt_delta <= 0 || adapt_delta >= 1 ||
      is.na(max_treedepth) || max_treedepth < 1L ||
      is.na(variational_iter) || variational_iter < 1L ||
      is.na(variational_output_samples) || variational_output_samples < 1L ||
      is.na(stage2_draw_count) || stage2_draw_count < 1L ||
      is.na(checkpoint_count) || checkpoint_count < 0L ||
      is.na(retry_seed_offset) || retry_seed_offset < 0L) {
    tf_stage1_screening_stop("Invalid parallel or sampling configuration.")
  }
  supported_intervals <- c(0.90, 0.95)
  if (!any(abs(target_interval - supported_intervals) < 1e-12) ||
      !any(abs(confounder_interval - supported_intervals) < 1e-12)) {
    tf_stage1_screening_stop(
      "`target_interval` and `confounder_interval` must be 0.90 or 0.95."
    )
  }
  worker_count <- min(cores, length(ready_genes))
  sampling <- list(
    inference = inference,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    refresh = refresh,
    adapt_delta = as.numeric(adapt_delta),
    max_treedepth = max_treedepth,
    prior_family = prior_family,
    target_interval = target_interval,
    confounder_interval = confounder_interval,
    variational_algorithm = variational_algorithm,
    variational_iter = variational_iter,
    variational_output_samples = variational_output_samples,
    stage2_draw_count = stage2_draw_count
  )

  stan_file <- tf_stage1_screening_resolve_stan_file(
    screening_input,
    stan_file,
    prior_family
  )
  if (is.null(output_file)) {
    output_file <- tf_stage1_screening_default_output(screening_input, inference)
  }
  output_file <- normalizePath(output_file, winslash = "/", mustWork = FALSE)
  if (is.null(checkpoint_file)) {
    checkpoint_file <- sub("\\.rds$", "_checkpoint.rds", output_file)
    if (identical(checkpoint_file, output_file)) {
      checkpoint_file <- paste0(output_file, "_checkpoint.rds")
    }
  }
  checkpoint_file <- normalizePath(
    checkpoint_file,
    winslash = "/",
    mustWork = FALSE
  )
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(checkpoint_file), recursive = TRUE, showWarnings = FALSE)

  is_compatible_result <- function(result) {
    result_inference <- attr(result, "inference")
    if (is.null(result_inference)) result_inference <- "mcmc"
    result_sampling <- attr(result, "sampling")
    inference_settings_match <- if (identical(inference, "variational")) {
      is.list(result_sampling) &&
        identical(
          result_sampling$variational_algorithm,
          variational_algorithm
        ) &&
        identical(result_sampling$variational_iter, variational_iter) &&
        identical(
          result_sampling$variational_output_samples,
          variational_output_samples
        )
    } else {
      is.list(result_sampling) &&
        identical(result_sampling$chains, chains) &&
        identical(result_sampling$iter_warmup, iter_warmup) &&
        identical(result_sampling$iter_sampling, iter_sampling) &&
        identical(result_sampling$seed, seed) &&
        identical(result_sampling$adapt_delta, as.numeric(adapt_delta)) &&
        identical(result_sampling$max_treedepth, max_treedepth)
    }
    stage2_draw_settings_match <- is.list(result_sampling) &&
      identical(result_sampling$stage2_draw_count, stage2_draw_count)
    prior_and_interval_settings_match <- is.list(result_sampling) &&
      identical(result_sampling$prior_family, prior_family) &&
      isTRUE(all.equal(
        result_sampling$target_interval, target_interval, tolerance = 1e-12
      )) &&
      isTRUE(all.equal(
        result_sampling$confounder_interval,
        confounder_interval,
        tolerance = 1e-12
      ))
    is.list(result) &&
      identical(attr(result, "input_signature"), input_signature) &&
      identical(attr(result, "stan_file"), stan_file) &&
      identical(
        attr(result, "model_version"),
        "stage1_normal_condition_interaction_stage2_ready_v2"
      ) &&
      identical(result_inference, inference) &&
      inference_settings_match &&
      prior_and_interval_settings_match &&
      stage2_draw_settings_match
  }
  has_error_result <- function(result) {
    any(vapply(
      result,
      function(entry) is.list(entry) && identical(entry$status, "error"),
      logical(1)
    ))
  }
  existing <- NULL
  retry_genes <- character(0)
  if (!isTRUE(force_refit) && file.exists(output_file)) {
    candidate_existing <- readRDS(output_file)
    if (is_compatible_result(candidate_existing)) {
      existing <- candidate_existing
      retry_genes <- names(existing)[vapply(
        existing,
        function(entry) is.list(entry) && identical(entry$status, "error"),
        logical(1)
      )]
    }
    if (!is.null(existing) && isTRUE(attr(existing, "complete")) &&
        (!isTRUE(retry_errors) || !has_error_result(existing))) {
      return(existing)
    }
  }

  results <- setNames(vector("list", length(ready_genes)), ready_genes)
  if (isTRUE(resume) && !is.null(existing)) {
    reusable <- intersect(names(existing), ready_genes)
    if (isTRUE(retry_errors)) {
      reusable <- reusable[vapply(
        existing[reusable],
        function(entry) is.list(entry) && identical(entry$status, "ok"),
        logical(1)
      )]
    }
    results[reusable] <- existing[reusable]
  }
  if (isTRUE(resume) && !isTRUE(force_refit) && file.exists(checkpoint_file)) {
    checkpoint <- readRDS(checkpoint_file)
    if (is_compatible_result(checkpoint)) {
      reusable <- intersect(names(checkpoint), ready_genes)
      if (isTRUE(retry_errors)) {
        reusable <- reusable[vapply(
          checkpoint[reusable],
          function(entry) is.list(entry) && identical(entry$status, "ok"),
          logical(1)
        )]
      }
      reusable <- reusable[vapply(results[reusable], is.null, logical(1))]
      results[reusable] <- checkpoint[reusable]
    }
  }

  completed <- !vapply(results, is.null, logical(1))
  pending_genes <- ready_genes[!completed]
  inference_description <- if (identical(inference, "mcmc")) {
    sprintf("%d sequential chains per model", chains)
  } else {
    sprintf(
      "%s ADVI with %d output draws per model",
      variational_algorithm,
      variational_output_samples
    )
  }
  message(sprintf(
    paste0(
      "Stage 1 regression: %d ready genes, %d completed, %d pending, ",
      "%d model workers, %s."
    ),
    length(ready_genes),
    sum(completed),
    length(pending_genes),
    worker_count,
    inference_description
  ))

  stamp_result <- function(result, complete) {
    class(result) <- c("TFStage1ScreeningRegressionList", "list")
    attr(result, "input_signature") <- input_signature
    attr(result, "model_version") <-
      "stage1_normal_condition_interaction_stage2_ready_v2"
    attr(result, "target_tf") <- attr(screening_input, "target_tf")
    attr(result, "stan_file") <- stan_file
    attr(result, "inference") <- inference
    attr(result, "prior_family") <- prior_family
    attr(result, "target_interval") <- target_interval
    attr(result, "confounder_interval") <- confounder_interval
    attr(result, "model_workers") <- worker_count
    attr(result, "chains_per_model") <- chains
    attr(result, "sampling") <- sampling
    attr(result, "input_ready_genes") <- ready_genes
    attr(result, "input_skipped_genes") <- skipped_genes
    attr(result, "input_pipeline_stage") <-
      attr(screening_input, "pipeline_stage")
    attr(result, "prescreen_prior_family") <-
      attr(screening_input, "prescreen_prior_family")
    attr(result, "prescreen_target_interval") <-
      attr(screening_input, "prescreen_target_interval")
    attr(result, "prescreen_confounder_interval") <-
      attr(screening_input, "prescreen_confounder_interval")
    attr(result, "prescreen_filter_counts") <-
      attr(screening_input, "prescreen_filter_counts")
    attr(result, "output_file") <- output_file
    attr(result, "checkpoint_file") <- checkpoint_file
    attr(result, "complete") <- isTRUE(complete)
    result
  }

  if (length(pending_genes)) {
    # Compile exactly once in the main R process before starting workers.
    compiled_model <- cmdstanr::cmdstan_model(
      stan_file = stan_file,
      force_recompile = isTRUE(force_recompile)
    )
    worker_exe_file <- compiled_model$exe_file()
    worker_stan_file <- stan_file
    rm(compiled_model)

    cluster <- parallel::makeCluster(worker_count)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterExport(
      cluster,
      varlist = c(
        "worker_stan_file",
        "worker_exe_file",
        "tf_stage1_screening_stop",
        "tf_stage1_screening_summarize_fit",
        "tf_stage1_screening_fit_one"
      ),
      envir = environment()
    )
    parallel::clusterEvalQ(cluster, {
      .tf_stage1_worker_model <- cmdstanr::cmdstan_model(
        stan_file = worker_stan_file,
        exe_file = worker_exe_file,
        compile = FALSE
      )
      NULL
    })

    batches <- split(
      pending_genes,
      ceiling(seq_along(pending_genes) / worker_count)
    )
    checkpoint_targets <- if (checkpoint_count > 0L) {
      unique(as.integer(ceiling(
        seq_len(checkpoint_count) * length(ready_genes) /
          (checkpoint_count + 1L)
      )))
    } else {
      integer(0)
    }
    last_checkpoint_completed <- sum(completed)
    model_index <- setNames(match(ready_genes, ready_genes), ready_genes)
    for (batch_index in seq_along(batches)) {
      genes <- batches[[batch_index]]
      jobs <- lapply(genes, function(gene) {
        list(
          entry = screening_input[[gene]],
          model_index = model_index[[gene]],
          retry_seed_offset = if (gene %in% retry_genes) {
            retry_seed_offset
          } else {
            0L
          }
        )
      })
      batch_results <- parallel::parLapplyLB(
        cluster,
        jobs,
        tf_stage1_screening_fit_one,
        sampling = sampling
      )
      names(batch_results) <- genes
      results[genes] <- batch_results
      completed_count <- sum(!vapply(results, is.null, logical(1)))
      crossed_target <- checkpoint_targets > last_checkpoint_completed &
        checkpoint_targets <= completed_count
      if (any(crossed_target)) {
        checkpoint <- results[!vapply(results, is.null, logical(1))]
        checkpoint <- stamp_result(checkpoint, complete = FALSE)
        saveRDS(checkpoint, checkpoint_file)
        last_checkpoint_completed <- completed_count
        message(sprintf(
          "Stage 1 regression checkpoint: %d/%d genes completed.",
          completed_count,
          length(ready_genes)
        ))
      }
    }
    parallel::stopCluster(cluster)
    on.exit(NULL, add = FALSE)
  }

  results <- results[ready_genes]
  results <- stamp_result(results, complete = TRUE)
  saveRDS(results, output_file)
  results
}
