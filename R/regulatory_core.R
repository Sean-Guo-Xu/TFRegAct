#!/usr/bin/env Rscript

tf_model_stop <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

tf_model_require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    tf_model_stop("Package '%s' is required.", pkg)
  }
}

tf_stan_model_cache <- new.env(parent = emptyenv())

tf_get_cmdstan_model <- function(stan_file, force_recompile = FALSE, include_log_lik = TRUE) {
  tf_model_require_pkg("cmdstanr")

  stan_file <- .tfregact_stan_file(stan_file)
  if (!isTRUE(include_log_lik)) stan_file <- .tfregact_stan_without_log_lik(stan_file)
  file_info <- file.info(stan_file)
  cached <- if (exists(stan_file, envir = tf_stan_model_cache, inherits = FALSE)) {
    get(stan_file, envir = tf_stan_model_cache, inherits = FALSE)
  } else {
    NULL
  }

  cache_is_current <- !is.null(cached) &&
    identical(cached$mtime, file_info$mtime) &&
    identical(cached$size, file_info$size)
  if (!isTRUE(force_recompile) && isTRUE(cache_is_current)) {
    return(cached$model)
  }

  model <- cmdstanr::cmdstan_model(
    stan_file,
    force_recompile = isTRUE(force_recompile)
  )
  assign(
    stan_file,
    list(mtime = file_info$mtime, size = file_info$size, model = model),
    envir = tf_stan_model_cache
  )
  model
}

tf_clear_stan_model_cache <- function(stan_file = NULL) {
  if (is.null(stan_file)) {
    rm(list = ls(envir = tf_stan_model_cache), envir = tf_stan_model_cache)
    return(invisible(NULL))
  }

  stan_file <- .tfregact_stan_file(stan_file)
  if (exists(stan_file, envir = tf_stan_model_cache, inherits = FALSE)) {
    rm(list = stan_file, envir = tf_stan_model_cache)
  }
  invisible(NULL)
}

tf_model_validate_object <- function(x) {
  if (inherits(x, "TFAnalysisObject")) {
    if (exists("tf_validate_analysis_object", mode = "function")) {
      tf_validate_analysis_object(x)
    }
    return(invisible(TRUE))
  }

  if (is.list(x) && exists("tf_validate_analysis_object", mode = "function")) {
    tf_validate_analysis_object(x)
    return(invisible(TRUE))
  }

  tf_model_stop("`analysis_object` must be a valid TFAnalysisObject.")
}

tf_direction_to_int <- function(direction) {
  if (is.numeric(direction) || is.integer(direction)) {
    out <- suppressWarnings(as.integer(direction))
    bad <- !(out %in% c(-1L, 0L, 1L))
    if (any(bad, na.rm = TRUE)) {
      tf_model_stop("Numeric `direction` values must be in {-1, 0, 1}.")
    }
    return(out)
  }

  direction_chr <- trimws(tolower(as.character(direction)))
  out <- rep(0L, length(direction_chr))
  out[direction_chr %in% c("activation", "activate", "positive", "up")] <- 1L
  out[direction_chr %in% c("repression", "repress", "negative", "down")] <- -1L
  out[direction_chr %in% c("unknown", "mixed", "", "na", "none")] <- 0L

  recognized <- direction_chr %in% c(
    "activation", "activate", "positive", "up",
    "repression", "repress", "negative", "down",
    "unknown", "mixed", "", "na", "none"
  )
  if (any(!recognized)) {
    bad_vals <- unique(direction_chr[!recognized])
    tf_model_stop(
      "Unrecognized direction value(s): %s",
      paste(bad_vals, collapse = ", ")
    )
  }

  out
}

tf_prepare_confidence <- function(confidence, min_value = 0, max_value = 10) {
  conf_num <- as.numeric(confidence)
  if (any(is.na(conf_num))) {
    tf_model_stop("`confidence` contains NA values after numeric coercion.")
  }
  conf_num <- pmin(pmax(conf_num, min_value), max_value)
  conf_num
}

tf_offset_from_libsize <- function(libsize) {
  libsize_num <- as.numeric(libsize)
  if (any(!is.finite(libsize_num)) || any(libsize_num <= 0)) {
    tf_model_stop("`libsize` must contain finite positive values.")
  }

  mean_lib <- mean(libsize_num)
  log(libsize_num / mean_lib)
}

tf_batch_info <- function(analysis_object, warn_single_level = TRUE) {
  batch <- analysis_object$batch
  n_cells <- ncol(analysis_object$expr)

  if (is.null(batch)) {
    return(list(use_batch = FALSE, index = NULL, levels = character(0)))
  }
  if (length(batch) != n_cells) {
    tf_model_stop("`batch` length must equal ncol(expr).")
  }

  batch_chr <- trimws(as.character(batch))
  if (any(is.na(batch_chr)) || any(batch_chr == "")) {
    tf_model_stop("`batch` must not contain NA or empty strings when provided.")
  }

  batch_factor <- droplevels(as.factor(batch_chr))
  batch_levels <- levels(batch_factor)
  if (length(batch_levels) < 2) {
    if (isTRUE(warn_single_level)) {
      warning(
        sprintf(
          "`batch` has fewer than 2 levels (%s); omitting batch columns from `W`.",
          if (length(batch_levels) == 0) "none" else paste(batch_levels, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    return(list(use_batch = FALSE, index = NULL, levels = batch_levels))
  }

  list(
    use_batch = TRUE,
    index = as.integer(batch_factor),
    levels = batch_levels
  )
}

tf_nuisance_prior_scales <- function(prior_scale, column_names) {
  Q <- length(column_names)
  if (Q == 0L) {
    return(numeric(0))
  }

  prior_scale <- as.numeric(prior_scale)
  if (any(!is.finite(prior_scale)) || any(prior_scale <= 0)) {
    tf_model_stop("`nuisance_prior_scale` must contain finite positive values.")
  }
  if (length(prior_scale) == 1L) {
    return(rep(prior_scale, Q))
  }
  if (length(prior_scale) != Q) {
    tf_model_stop(
      "`nuisance_prior_scale` must have length 1 or match the %d columns of `W`.",
      Q
    )
  }
  prior_scale
}

tf_orthonormal_batch_contrasts <- function(level_count) {
  level_count <- as.integer(level_count[[1]])
  if (is.na(level_count) || level_count < 1L) {
    tf_model_stop("`level_count` must be a positive integer.")
  }
  if (level_count == 1L) {
    return(matrix(numeric(0), nrow = 1L, ncol = 0L))
  }

  contrasts <- stats::contr.helmert(level_count)
  contrasts <- sweep(
    contrasts,
    2L,
    sqrt(colSums(contrasts ^ 2)),
    "/"
  )
  storage.mode(contrasts) <- "double"
  contrasts
}

tf_prepare_nuisance_design <- function(
  analysis_object,
  control_level = NULL,
  disease_level = NULL,
  nuisance_prior_scale = 1
) {
  n_cells <- ncol(analysis_object$expr)
  condition <- tf_condition_from_sample(
    analysis_object = analysis_object,
    control_level = control_level,
    disease_level = disease_level
  )
  if (!is.null(analysis_object$sample) && is.null(condition)) {
    tf_model_stop(
      paste(
        "A supplied condition must contain exactly two non-missing levels.",
        "Set `condition_column = NULL` to omit condition from the nuisance design."
      )
    )
  }
  has_condition <- !is.null(condition)
  batch_info <- tf_batch_info(analysis_object, warn_single_level = FALSE)

  design_parts <- list()
  if (has_condition) {
    condition_column <- matrix(as.numeric(condition) - 0.5, ncol = 1L)
    colnames(condition_column) <- "condition"
    design_parts[[length(design_parts) + 1L]] <- condition_column
  }

  if (isTRUE(batch_info$use_batch)) {
    batch_contrasts <- tf_orthonormal_batch_contrasts(
      length(batch_info$levels)
    )
    batch_matrix <- batch_contrasts[batch_info$index, , drop = FALSE]
    colnames(batch_matrix) <- paste0(
      "batch_contrast_",
      seq_len(ncol(batch_matrix))
    )
    design_parts[[length(design_parts) + 1L]] <- batch_matrix
  }

  W <- if (length(design_parts) == 0L) {
    matrix(numeric(0), nrow = n_cells, ncol = 0L)
  } else {
    do.call(cbind, design_parts)
  }
  storage.mode(W) <- "double"

  if (ncol(W) > 0L) {
    design_with_intercept <- cbind(`(Intercept)` = 1, W)
    qr_design <- qr(design_with_intercept)
    if (qr_design$rank < ncol(design_with_intercept)) {
      tf_model_stop(
        paste(
          "The nuisance design is rank deficient.",
          "Condition may be completely confounded with batch, or a supplied covariate may be redundant."
        )
      )
    }
  }

  W_names <- colnames(W)
  W_prior_scale <- tf_nuisance_prior_scales(
    prior_scale = nuisance_prior_scale,
    column_names = W_names
  )
  condition_w_index <- if (has_condition) match("condition", W_names) else 0L
  K_batch <- if (isTRUE(batch_info$use_batch)) length(batch_info$levels) else 0L
  batch_index <- if (isTRUE(batch_info$use_batch)) {
    as.integer(batch_info$index)
  } else {
    rep.int(0L, n_cells)
  }
  batch_level_design <- matrix(0, nrow = K_batch, ncol = ncol(W))
  if (K_batch > 0L) {
    batch_columns <- grep("^batch_contrast_", W_names)
    batch_level_design[, batch_columns] <-
      tf_orthonormal_batch_contrasts(K_batch)
  }
  colnames(batch_level_design) <- W_names

  list(
    Q = ncol(W),
    W = W,
    W_prior_scale = W_prior_scale,
    W_names = W_names,
    has_condition = has_condition,
    condition = if (has_condition) as.integer(condition) else rep.int(0L, n_cells),
    condition_w_index = as.integer(condition_w_index),
    control_level = attr(condition, "control_level"),
    disease_level = attr(condition, "disease_level"),
    use_batch = isTRUE(batch_info$use_batch),
    batch_levels = if (K_batch > 0L) batch_info$levels else character(0),
    K_batch = as.integer(K_batch),
    batch = batch_index,
    batch_level_design = batch_level_design
  )
}

tf_prepare_stan_data <- function(
  analysis_object,
  gamma,
  eta,
  r_dir,
  confidence_min,
  confidence_max,
  control_level = NULL,
  disease_level = NULL,
  nuisance_prior_scale = 1,
  alpha_prior_sd = 1,
  target_interaction = FALSE,
  target_interaction_sd = 0.5
) {
  tf_model_validate_object(analysis_object)

  expr <- as.matrix(analysis_object$expr)
  if (!is.numeric(expr)) {
    storage.mode(expr) <- "numeric"
  }

  y_vec <- suppressWarnings(as.integer(analysis_object$Y_exp))
  if (any(is.na(y_vec)) || any(y_vec < 0)) {
    tf_model_stop("`Y_exp` must be non-negative integers.")
  }

  nuisance <- tf_prepare_nuisance_design(
    analysis_object = analysis_object,
    control_level = control_level,
    disease_level = disease_level,
    nuisance_prior_scale = nuisance_prior_scale
  )
  direction_index <- tf_direction_to_int(analysis_object$direction)
  confidence <- tf_prepare_confidence(
    confidence = analysis_object$confidence,
    min_value = confidence_min,
    max_value = confidence_max
  )
  log_offset <- tf_offset_from_libsize(analysis_object$libsize)
  target_tf <- trimws(as.character(analysis_object$target[[1]]))
  target_tf_index <- match(toupper(target_tf), toupper(rownames(expr)))
  if (is.na(target_tf_index)) {
    tf_model_stop("Target TF `%s` was not found in `analysis_object$expr`.", target_tf)
  }
  if (!is.numeric(alpha_prior_sd) || length(alpha_prior_sd) != 1L ||
      !is.finite(alpha_prior_sd) || alpha_prior_sd <= 0) {
    tf_model_stop("`alpha_prior_sd` must be one finite positive number.")
  }
  if (!is.logical(target_interaction) || length(target_interaction) != 1L ||
      is.na(target_interaction)) {
    tf_model_stop("`target_interaction` must be TRUE or FALSE.")
  }
  if (!is.numeric(target_interaction_sd) ||
      length(target_interaction_sd) != 1L ||
      !is.finite(target_interaction_sd) || target_interaction_sd <= 0) {
    tf_model_stop("`target_interaction_sd` must be one finite positive number.")
  }
  if (isTRUE(target_interaction) && !isTRUE(nuisance$has_condition)) {
    tf_model_stop("Target interaction requires a valid two-level condition.")
  }
  use_target_interaction <- isTRUE(target_interaction) &&
    isTRUE(nuisance$has_condition)

  stan_data <- list(
    N = ncol(expr),
    P = nrow(expr),
    Y = y_vec,
    X = t(expr),
    Q = as.integer(nuisance$Q),
    W = nuisance$W,
    W_prior_scale = as.vector(nuisance$W_prior_scale),
    has_condition = as.integer(nuisance$has_condition),
    condition = as.array(nuisance$condition),
    condition_w_index = nuisance$condition_w_index,
    K_batch = nuisance$K_batch,
    batch = as.array(nuisance$batch),
    batch_level_design = nuisance$batch_level_design,
    log_offset = as.vector(log_offset),
    confidence = as.vector(confidence),
    direction = as.array(direction_index),
    gamma = as.numeric(gamma),
    eta = as.numeric(eta),
    r_dir = as.numeric(r_dir),
    alpha_prior_sd = as.numeric(alpha_prior_sd),
    target_tf_index = as.integer(target_tf_index),
    use_target_interaction = as.integer(use_target_interaction),
    target_condition_centered = if (use_target_interaction) {
      as.numeric(nuisance$condition) - 0.5
    } else {
      rep(0, ncol(expr))
    },
    target_interaction_sd = as.numeric(target_interaction_sd)
  )

  attr(stan_data, "feature_names") <- rownames(expr)
  attr(stan_data, "cell_names") <- colnames(expr)
  attr(stan_data, "target") <- analysis_object$target
  attr(stan_data, "confidence_used") <- confidence
  attr(stan_data, "W_names") <- nuisance$W_names
  attr(stan_data, "batch_levels") <- nuisance$batch_levels
  attr(stan_data, "use_batch_model") <- nuisance$use_batch
  attr(stan_data, "condition_model") <- nuisance$has_condition
  attr(stan_data, "control_level") <- nuisance$control_level
  attr(stan_data, "disease_level") <- nuisance$disease_level
  attr(stan_data, "target_tf_index") <- target_tf_index
  attr(stan_data, "target_interaction_model") <- use_target_interaction
  stan_data
}

tf_default_stan_file <- function(use_batch_model = NULL) {
  .tfregact_stan_file("TF_directional_nb_model.stan")
}

tf_select_stage1_stan_file <- function(stan_file, use_batch_model = NULL) {
  if (is.null(stan_file) || !nzchar(stan_file)) {
    return(tf_default_stan_file())
  }

  .tfregact_stan_file(stan_file)
}

run_TF_directional_model <- function(
  analysis_object,
  stan_file,
  gamma,
  eta,
  r_dir,
  confidence_min,
  confidence_max,
  chains,
  parallel_chains,
  iter_warmup,
  iter_sampling,
  seed,
  refresh,
  force_recompile,
  control_level = NULL,
  disease_level = NULL,
  nuisance_prior_scale = 1,
  alpha_prior_sd = 1,
  ...
) {
  tf_model_require_pkg("cmdstanr")

  stan_data <- tf_prepare_stan_data(
    analysis_object = analysis_object,
    gamma = gamma,
    eta = eta,
    r_dir = r_dir,
    confidence_min = confidence_min,
    confidence_max = confidence_max,
    control_level = control_level,
    disease_level = disease_level,
    nuisance_prior_scale = nuisance_prior_scale,
    alpha_prior_sd = alpha_prior_sd,
    target_interaction = FALSE
  )
  use_batch_model <- isTRUE(attr(stan_data, "use_batch_model"))
  stan_file <- tf_select_stage1_stan_file(stan_file, use_batch_model)

  model <- tf_get_cmdstan_model(
    stan_file,
    force_recompile = force_recompile,
    include_log_lik = FALSE
  )
  fit <- tf_fit_shared_laplace_stage1(
    stan_data = stan_data,
    inference = "mcmc",
    model = model,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    refresh = refresh,
    ...
  )

  list(
    fit = fit,
    stan_data = stan_data,
    stan_file = stan_file,
    target = attr(stan_data, "target"),
    feature_names = attr(stan_data, "feature_names"),
    cell_names = attr(stan_data, "cell_names"),
    use_batch_model = use_batch_model,
    condition_model = isTRUE(attr(stan_data, "condition_model")),
    control_level = attr(stan_data, "control_level"),
    disease_level = attr(stan_data, "disease_level"),
    W_names = attr(stan_data, "W_names")
  )
}

tf_get_beta_draws_matrix <- function(fit_result) {
  fit <- fit_result
  feature_names <- NULL

  if (is.list(fit_result) && !is.null(fit_result$fit)) {
    fit <- fit_result$fit
    feature_names <- fit_result$feature_names
  }

  if (is.null(fit) || !is.function(fit$draws)) {
    tf_model_stop(
      "`fit_result` must be the return value of `run_TF_directional_model()` or a cmdstanr fit object."
    )
  }

  beta_draws <- fit$draws(variables = "beta", format = "draws_matrix")
  beta_mat <- as.matrix(beta_draws)
  beta_cols <- grep("^beta\\[[0-9]+\\]$", colnames(beta_mat), value = TRUE)

  if (length(beta_cols) == 0) {
    tf_model_stop("No beta draws were found in `fit_result`.")
  }

  beta_mat <- beta_mat[, beta_cols, drop = FALSE]

  if (is.null(feature_names)) {
    beta_index <- as.integer(sub("^beta\\[([0-9]+)\\]$", "\\1", beta_cols))
    feature_names <- paste0("TF_", beta_index)
  }

  if (length(feature_names) != ncol(beta_mat)) {
    tf_model_stop(
      "Length of feature names (%d) does not match number of beta columns (%d).",
      length(feature_names),
      ncol(beta_mat)
    )
  }

  colnames(beta_mat) <- feature_names
  beta_mat
}

tf_plot_tf_correlation_heatmap <- function(
  cor_mat,
  output_dir = ".",
  filename = "filtered_tf_correlation_heatmap.pdf",
  width = 8,
  height = 7
) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    warning(
      "Package 'pheatmap' is not installed; skipping TF correlation heatmap.",
      call. = FALSE
    )
    return(NULL)
  }

  if (is.null(output_dir) || !nzchar(output_dir)) {
    output_dir <- "."
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  heatmap_file <- file.path(output_dir, filename)
  heatmap_mat <- cor_mat
  heatmap_mat[upper.tri(heatmap_mat)] <- NA_real_
  number_mat <- matrix(
    "",
    nrow = nrow(heatmap_mat),
    ncol = ncol(heatmap_mat),
    dimnames = dimnames(heatmap_mat)
  )
  number_mat[!is.na(heatmap_mat)] <- sprintf("%.2f", heatmap_mat[!is.na(heatmap_mat)])

  nature_blue_red <- grDevices::colorRampPalette(
    c("#2166AC", "#F7F7F7", "#B2182B")
  )(101)

  pheatmap::pheatmap(
    heatmap_mat,
    color = nature_blue_red,
    breaks = seq(-1, 1, length.out = length(nature_blue_red) + 1),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = number_mat,
    number_color = "black",
    fontsize_number = 8,
    na_col = "white",
    border_color = NA,
    filename = heatmap_file,
    width = width,
    height = height
  )

  normalizePath(heatmap_file, winslash = "/", mustWork = FALSE)
}

tf_prune_correlated_tfs <- function(
  expr,
  candidate_tfs,
  beta_mean,
  target_tf,
  correlation_threshold = 0.7
) {
  result <- list(
    keep_tfs = candidate_tfs,
    cor_mat = NULL,
    prune_table = data.frame(
      tf = candidate_tfs,
      correlation_pruned = FALSE,
      correlated_with = NA_character_,
      pruning_correlation = NA_real_,
      pruning_beta_mean = NA_real_,
      pruning_partner_beta_mean = NA_real_,
      stringsAsFactors = FALSE
    )
  )

  if (length(candidate_tfs) < 2) {
    return(result)
  }

  expr_sub <- as.matrix(expr[candidate_tfs, , drop = FALSE])
  cor_mat <- stats::cor(t(expr_sub), use = "pairwise.complete.obs")
  cor_mat[is.na(cor_mat)] <- 0
  diag(cor_mat) <- 1
  result$cor_mat <- cor_mat

  pair_idx <- which(
    lower.tri(cor_mat) & abs(cor_mat) >= correlation_threshold,
    arr.ind = TRUE
  )
  if (nrow(pair_idx) == 0) {
    return(result)
  }

  pair_df <- data.frame(
    tf1 = rownames(cor_mat)[pair_idx[, 1]],
    tf2 = colnames(cor_mat)[pair_idx[, 2]],
    correlation = cor_mat[pair_idx],
    stringsAsFactors = FALSE
  )
  pair_df <- pair_df[order(abs(pair_df$correlation), decreasing = TRUE), , drop = FALSE]

  keep_map <- stats::setNames(rep(TRUE, length(candidate_tfs)), candidate_tfs)
  beta_lookup <- stats::setNames(as.numeric(beta_mean[candidate_tfs]), candidate_tfs)

  for (i in seq_len(nrow(pair_df))) {
    tf1 <- pair_df$tf1[[i]]
    tf2 <- pair_df$tf2[[i]]
    if (!isTRUE(keep_map[[tf1]]) || !isTRUE(keep_map[[tf2]])) {
      next
    }

    cor_val <- pair_df$correlation[[i]]
    beta1 <- beta_lookup[[tf1]]
    beta2 <- beta_lookup[[tf2]]
    if (is.na(beta1) || is.na(beta2) || beta1 == 0 || beta2 == 0 || cor_val == 0) {
      next
    }

    direction_consistent <- sign(cor_val) == sign(beta1 * beta2)
    if (!isTRUE(direction_consistent)) {
      next
    }

    if (identical(tf1, target_tf) && !identical(tf2, target_tf)) {
      drop_tf <- tf2
      keep_tf <- tf1
    } else if (identical(tf2, target_tf) && !identical(tf1, target_tf)) {
      drop_tf <- tf1
      keep_tf <- tf2
    } else if (abs(beta1) >= abs(beta2)) {
      drop_tf <- tf2
      keep_tf <- tf1
    } else {
      drop_tf <- tf1
      keep_tf <- tf2
    }

    keep_map[[drop_tf]] <- FALSE
    drop_row <- match(drop_tf, result$prune_table$tf)
    result$prune_table$correlation_pruned[[drop_row]] <- TRUE
    result$prune_table$correlated_with[[drop_row]] <- keep_tf
    result$prune_table$pruning_correlation[[drop_row]] <- cor_val
    result$prune_table$pruning_beta_mean[[drop_row]] <- beta_lookup[[drop_tf]]
    result$prune_table$pruning_partner_beta_mean[[drop_row]] <- beta_lookup[[keep_tf]]
  }

  result$keep_tfs <- names(keep_map)[keep_map]
  result
}

filter_TF_analysis_object_by_beta_ci <- function(
  fit_result,
  analysis_object,
  interval = c(10, 90),
  force_keep_target_tf = TRUE,
  correlation_filter = TRUE,
  correlation_threshold = 0.7,
  plot_correlation_heatmap = TRUE,
  correlation_heatmap_dir = ".",
  correlation_heatmap_file = NULL
) {
  tf_model_validate_object(analysis_object)

  if (!is.numeric(interval) || length(interval) != 2 || any(is.na(interval))) {
    tf_model_stop("`interval` must be a numeric vector of length 2, for example c(10, 90).")
  }

  interval <- sort(as.numeric(interval))
  if (interval[[1]] < 0 || interval[[2]] > 100 || interval[[1]] >= interval[[2]]) {
    tf_model_stop("`interval` must contain two ordered percentile values between 0 and 100.")
  }
  if (!is.numeric(correlation_threshold) ||
      length(correlation_threshold) != 1 ||
      is.na(correlation_threshold) ||
      correlation_threshold < 0 ||
      correlation_threshold > 1) {
    tf_model_stop("`correlation_threshold` must be a single numeric value between 0 and 1.")
  }

  beta_mat <- tf_get_beta_draws_matrix(fit_result)
  object_tfs <- rownames(analysis_object$expr)
  beta_tfs <- colnames(beta_mat)
  target_label <- paste(as.character(analysis_object$target), collapse = "_")
  target_label <- gsub("[^A-Za-z0-9_.-]+", "_", target_label)
  if (is.null(correlation_heatmap_file) || !nzchar(correlation_heatmap_file)) {
    correlation_heatmap_file <- sprintf("%s_filtered_corheatmap.pdf", target_label)
  }

  missing_tfs <- setdiff(object_tfs, beta_tfs)
  if (length(missing_tfs) > 0) {
    tf_model_stop(
      "Some TFs in `analysis_object` were not found in beta draws: %s",
      paste(head(missing_tfs, 10), collapse = ", ")
    )
  }

  beta_mat <- beta_mat[, object_tfs, drop = FALSE]
  ci_mat <- apply(
    beta_mat,
    2,
    stats::quantile,
    probs = interval / 100,
    na.rm = TRUE,
    names = FALSE
  )

  ci_lower <- as.numeric(ci_mat[1, ])
  ci_upper <- as.numeric(ci_mat[2, ])
  beta_mean <- colMeans(beta_mat, na.rm = TRUE)
  beta_sd <- apply(beta_mat, 2, stats::sd, na.rm = TRUE)
  beta_median <- apply(beta_mat, 2, stats::median, na.rm = TRUE)
  ci_excludes_zero <- ci_lower > 0 | ci_upper < 0
  target_tf <- as.character(analysis_object$target[[1]])
  is_target_tf <- object_tfs == target_tf
  if (isTRUE(force_keep_target_tf) && !any(is_target_tf)) {
    warning(
      sprintf(
        "Target TF `%s` was not found in `analysis_object$expr`; it cannot be force-kept.",
        target_tf
      ),
      call. = FALSE
    )
  }
  keep_ci <- ci_excludes_zero | (isTRUE(force_keep_target_tf) & is_target_tf)
  kept_by_force <- keep_ci & !ci_excludes_zero & is_target_tf

  candidate_tfs <- object_tfs[keep_ci]
  correlation_result <- list(
    keep_tfs = candidate_tfs,
    cor_mat = NULL,
    prune_table = data.frame(
      tf = candidate_tfs,
      correlation_pruned = FALSE,
      correlated_with = NA_character_,
      pruning_correlation = NA_real_,
      pruning_beta_mean = NA_real_,
      pruning_partner_beta_mean = NA_real_,
      stringsAsFactors = FALSE
    )
  )
  correlation_heatmap_path <- NULL

  if (isTRUE(correlation_filter) && length(candidate_tfs) >= 2) {
    correlation_result <- tf_prune_correlated_tfs(
      expr = analysis_object$expr,
      candidate_tfs = candidate_tfs,
      beta_mean = beta_mean,
      target_tf = target_tf,
      correlation_threshold = correlation_threshold
    )
  } else if (length(candidate_tfs) >= 2) {
    expr_sub <- as.matrix(analysis_object$expr[candidate_tfs, , drop = FALSE])
    cor_mat <- stats::cor(t(expr_sub), use = "pairwise.complete.obs")
    cor_mat[is.na(cor_mat)] <- 0
    diag(cor_mat) <- 1
    correlation_result$cor_mat <- cor_mat
  }

  if (isTRUE(plot_correlation_heatmap) &&
      !is.null(correlation_result$cor_mat) &&
      nrow(correlation_result$cor_mat) >= 2) {
    correlation_heatmap_path <- tf_plot_tf_correlation_heatmap(
      cor_mat = correlation_result$cor_mat,
      output_dir = correlation_heatmap_dir,
      filename = correlation_heatmap_file
    )
  }

  keep <- keep_ci & object_tfs %in% correlation_result$keep_tfs
  if (isTRUE(force_keep_target_tf) && any(is_target_tf) && keep_ci[is_target_tf]) {
    keep[is_target_tf] <- TRUE
  }

  correlation_pruned <- rep(FALSE, length(object_tfs))
  correlated_with <- rep(NA_character_, length(object_tfs))
  pruning_correlation <- rep(NA_real_, length(object_tfs))
  pruning_beta_mean <- rep(NA_real_, length(object_tfs))
  pruning_partner_beta_mean <- rep(NA_real_, length(object_tfs))
  if (nrow(correlation_result$prune_table) > 0) {
    prune_match <- match(correlation_result$prune_table$tf, object_tfs)
    correlation_pruned[prune_match] <- correlation_result$prune_table$correlation_pruned
    correlated_with[prune_match] <- correlation_result$prune_table$correlated_with
    pruning_correlation[prune_match] <- correlation_result$prune_table$pruning_correlation
    pruning_beta_mean[prune_match] <- correlation_result$prune_table$pruning_beta_mean
    pruning_partner_beta_mean[prune_match] <-
      correlation_result$prune_table$pruning_partner_beta_mean
  }

  filter_table <- data.frame(
    tf = object_tfs,
    beta_mean = as.numeric(beta_mean),
    beta_sd = as.numeric(beta_sd),
    beta_median = as.numeric(beta_median),
    beta_ci_lower = ci_lower,
    beta_ci_upper = ci_upper,
    ci_excludes_zero = ci_excludes_zero,
    is_target_tf = is_target_tf,
    kept_by_force = kept_by_force,
    ci_kept = keep_ci,
    correlation_pruned = correlation_pruned,
    correlated_with = correlated_with,
    pruning_correlation = pruning_correlation,
    pruning_beta_mean = pruning_beta_mean,
    pruning_partner_beta_mean = pruning_partner_beta_mean,
    kept = keep,
    confidence = analysis_object$confidence,
    direction = analysis_object$direction,
    stringsAsFactors = FALSE
  )

  filtered_object <- create_TF_analysis_object(
    expr = analysis_object$expr[keep, , drop = FALSE],
    batch = analysis_object$batch,
    sample = analysis_object$sample,
    target = analysis_object$target,
    libsize = analysis_object$libsize,
    confidence = analysis_object$confidence[keep],
    direction = analysis_object$direction[keep],
    Y_exp = analysis_object$Y_exp
  )

  attr(filtered_object, "beta_ci_filter") <- filter_table
  attr(filtered_object, "stage1_beta_summary") <- filter_table[keep, , drop = FALSE]
  attr(filtered_object, "beta_ci_interval") <- interval
  attr(filtered_object, "tf_correlation_matrix") <- correlation_result$cor_mat
  attr(filtered_object, "tf_correlation_threshold") <- correlation_threshold
  attr(filtered_object, "tf_correlation_filter_enabled") <- isTRUE(correlation_filter)
  attr(filtered_object, "tf_correlation_heatmap_file") <- correlation_heatmap_path
  attr(filtered_object, "n_features_after_ci") <- sum(keep_ci)
  attr(filtered_object, "n_features_before") <- length(object_tfs)
  attr(filtered_object, "n_features_after") <- sum(keep)

  filtered_object
}

tf_default_stage2_stan_file <- function(condition_model = NULL, use_batch_model = NULL) {
  .tfregact_stan_file("TF_stage2_directional_nb_model.stan")
}

tf_select_stage2_stan_file <- function(
  stan_file,
  condition_model = NULL,
  use_batch_model = NULL
) {
  if (is.null(stan_file) || !nzchar(stan_file)) {
    return(tf_default_stage2_stan_file())
  }

  .tfregact_stan_file(stan_file)
}

tf_compute_loo <- function(fit, variables = "log_lik") {
  if (!requireNamespace("loo", quietly = TRUE)) {
    warning(
      "Package 'loo' is not installed; returning `loo = NULL`.",
      call. = FALSE
    )
    return(NULL)
  }
  if (is.null(fit) || !is.function(fit$draws)) {
    tf_model_stop("`fit` must be a cmdstanr fit object to compute LOO.")
  }

  log_lik <- tryCatch(
    fit$draws(variables = variables, format = "draws_matrix"),
    error = function(e) {
      tf_model_stop(
        "Failed to extract `%s` draws for LOO: %s",
        variables,
        conditionMessage(e)
      )
    }
  )
  loo::loo(as.matrix(log_lik))
}

tf_stage2_beta_summary <- function(analysis_object, stage1_fit_result = NULL) {
  object_tfs <- rownames(analysis_object$expr)

  if (!is.null(stage1_fit_result)) {
    beta_mat <- tf_get_beta_draws_matrix(stage1_fit_result)
    missing_tfs <- setdiff(object_tfs, colnames(beta_mat))
    if (length(missing_tfs) > 0) {
      tf_model_stop(
        "Some TFs in `analysis_object` were not found in `stage1_fit_result`: %s",
        paste(head(missing_tfs, 10), collapse = ", ")
      )
    }

    beta_mat <- beta_mat[, object_tfs, drop = FALSE]
    return(data.frame(
      tf = object_tfs,
      stage1_beta_mean = as.numeric(colMeans(beta_mat, na.rm = TRUE)),
      stage1_beta_sd = as.numeric(apply(beta_mat, 2, stats::sd, na.rm = TRUE)),
      stringsAsFactors = FALSE
    ))
  }

  summary_df <- attr(analysis_object, "stage1_beta_summary")
  if (is.null(summary_df)) {
    tf_model_stop(
      "`analysis_object` does not contain stage 1 beta summary. Please pass `stage1_fit_result` or create it with `filter_TF_analysis_object_by_beta_ci()` after sourcing this updated file."
    )
  }

  if (!all(c("tf", "beta_mean", "beta_sd") %in% colnames(summary_df))) {
    tf_model_stop(
      "`stage1_beta_summary` must contain columns `tf`, `beta_mean`, and `beta_sd`."
    )
  }

  summary_df <- summary_df[match(object_tfs, summary_df$tf), , drop = FALSE]
  if (any(is.na(summary_df$tf))) {
    tf_model_stop("`stage1_beta_summary` does not match all TFs in `analysis_object`.")
  }

  data.frame(
    tf = object_tfs,
    stage1_beta_mean = as.numeric(summary_df$beta_mean),
    stage1_beta_sd = as.numeric(summary_df$beta_sd),
    stringsAsFactors = FALSE
  )
}

tf_condition_from_sample <- function(
  analysis_object,
  control_level = NULL,
  disease_level = NULL
) {
  sample_vec <- analysis_object$sample
  if (is.null(sample_vec) || length(sample_vec) != ncol(analysis_object$expr)) {
    return(NULL)
  }

  sample_chr <- trimws(as.character(sample_vec))
  if (any(is.na(sample_chr)) || any(sample_chr == "")) {
    return(NULL)
  }

  sample_factor <- droplevels(as.factor(sample_chr))
  sample_levels <- levels(sample_factor)
  if (length(sample_levels) != 2) {
    return(NULL)
  }

  match_level <- function(value, argument_name) {
    value <- trimws(as.character(value[[1]]))
    exact_idx <- which(sample_levels == value)
    if (length(exact_idx) == 1) {
      return(sample_levels[[exact_idx]])
    }

    case_idx <- which(tolower(sample_levels) == tolower(value))
    if (length(case_idx) == 1) {
      return(sample_levels[[case_idx]])
    }

    tf_model_stop(
      "`%s` must match one level in `analysis_object$sample`. Available levels: %s",
      argument_name,
      paste(sample_levels, collapse = ", ")
    )
  }

  has_nonempty_level <- function(value) {
    if (is.null(value) || length(value) == 0) {
      return(FALSE)
    }
    value_chr <- trimws(as.character(value[[1]]))
    !is.na(value_chr) && nzchar(value_chr)
  }

  has_explicit_control <- has_nonempty_level(control_level)
  has_explicit_disease <- has_nonempty_level(disease_level)

  if (has_explicit_control || has_explicit_disease) {
    if (has_explicit_control) {
      control_level <- match_level(control_level, "control_level")
    }
    if (has_explicit_disease) {
      disease_level <- match_level(disease_level, "disease_level")
    }

    if (!has_explicit_control) {
      control_level <- setdiff(sample_levels, disease_level)[[1]]
    }
    if (!has_explicit_disease) {
      disease_level <- setdiff(sample_levels, control_level)[[1]]
    }

    if (identical(control_level, disease_level)) {
      tf_model_stop("`control_level` and `disease_level` must be different.")
    }
  } else {
    level_key <- tolower(sample_levels)
    control_terms <- c("control", "ctrl", "normal", "healthy", "con")
    disease_terms <- c("disease", "diseased", "case", "patient", "patients", "aaa")

    control_idx <- which(level_key %in% control_terms)
    disease_idx <- which(level_key %in% disease_terms)

    if (length(control_idx) == 1 && length(disease_idx) == 1) {
      control_level <- sample_levels[[control_idx]]
      disease_level <- sample_levels[[disease_idx]]
    } else {
      control_level <- sample_levels[[1]]
      disease_level <- sample_levels[[2]]
      warning(
        sprintf(
          "`sample` has two levels but they were not recognized as control/disease. Coding `%s` as 0 and `%s` as 1. To avoid guessing, pass `control_level` and `disease_level`.",
          control_level,
          disease_level
        ),
        call. = FALSE
      )
    }
  }

  condition <- ifelse(sample_chr == disease_level, 1L, 0L)
  attr(condition, "control_level") <- control_level
  attr(condition, "disease_level") <- disease_level
  condition
}

tf_prepare_stage2_stan_data <- function(
  analysis_object,
  stage1_fit_result,
  direction_effect,
  beta_sd_floor,
  stage1_sd_multiplier,
  control_level,
  disease_level,
  confidence_min,
  confidence_max,
  nuisance_prior_scale = 1,
  target_interaction = FALSE,
  target_interaction_sd = 0.5
) {
  tf_model_validate_object(analysis_object)

  if (!is.logical(target_interaction) || length(target_interaction) != 1L ||
      is.na(target_interaction)) {
    tf_model_stop("`target_interaction` must be TRUE or FALSE.")
  }
  if (!is.numeric(target_interaction_sd) ||
      length(target_interaction_sd) != 1L ||
      !is.finite(target_interaction_sd) ||
      target_interaction_sd <= 0) {
    tf_model_stop("`target_interaction_sd` must be one finite positive number.")
  }

  base_data <- tf_prepare_stan_data(
    analysis_object = analysis_object,
    gamma = 1,
    # Stage 2 carries confidence only for provenance; use the same default
    # scale exponent as Stage 1 so the shared data contract is coherent.
    eta = 0.5,
    r_dir = 3,
    confidence_min = confidence_min,
    confidence_max = confidence_max,
    control_level = control_level,
    disease_level = disease_level,
    nuisance_prior_scale = nuisance_prior_scale
  )
  beta_summary <- tf_stage2_beta_summary(
    analysis_object = analysis_object,
    stage1_fit_result = stage1_fit_result
  )
  use_condition_model <- isTRUE(attr(base_data, "condition_model"))
  if (isTRUE(target_interaction) && !use_condition_model) {
    tf_model_stop(
      "The target-interaction Stage 2 model requires a valid two-level condition in `analysis_object$sample`."
    )
  }

  stan_data_names <- c(
    "N",
    "P",
    "Y",
    "X",
    "Q",
    "W",
    "W_prior_scale",
    "has_condition",
    "condition",
    "condition_w_index",
    "K_batch",
    "batch",
    "batch_level_design",
    "log_offset",
    "alpha_prior_sd"
  )
  stan_data <- base_data[stan_data_names]

  stage2_prior <- tf_stage2_prior_from_stage1(
    stage1_beta_mean = beta_summary$stage1_beta_mean,
    stage1_beta_sd = beta_summary$stage1_beta_sd,
    direction = base_data$direction,
    direction_effect = direction_effect,
    beta_sd_floor = beta_sd_floor,
    stage1_sd_multiplier = stage1_sd_multiplier
  )
  stan_data$beta_prior_mean <- stage2_prior$beta_prior_mean
  stan_data$beta_prior_sd <- stage2_prior$beta_prior_sd
  stan_data$beta_init <- stage2_prior$stage1_beta_mean

  target_tf <- as.character(attr(base_data, "target")[[1]])
  target_tf_index <- match(target_tf, attr(base_data, "feature_names"))
  if (is.na(target_tf_index)) {
    tf_model_stop(
      "Target TF `%s` was not found in `analysis_object$expr`.",
      target_tf
    )
  }
  stan_data$target_tf_index <- as.integer(target_tf_index)
  stan_data$use_target_interaction <- as.integer(target_interaction)
  stan_data$target_condition_centered <- if (isTRUE(target_interaction)) {
    as.numeric(stan_data$condition) - 0.5
  } else {
    rep(0, stan_data$N)
  }
  stan_data$target_interaction_sd <- as.numeric(target_interaction_sd)

  attr(stan_data, "feature_names") <- attr(base_data, "feature_names")
  attr(stan_data, "cell_names") <- attr(base_data, "cell_names")
  attr(stan_data, "target") <- attr(base_data, "target")
  attr(stan_data, "confidence_used") <- attr(base_data, "confidence_used")
  attr(stan_data, "W_names") <- attr(base_data, "W_names")
  attr(stan_data, "batch_levels") <- attr(base_data, "batch_levels")
  attr(stan_data, "use_batch_model") <- attr(base_data, "use_batch_model")
  attr(stan_data, "condition_model") <- use_condition_model
  attr(stan_data, "control_level") <- attr(base_data, "control_level")
  attr(stan_data, "disease_level") <- attr(base_data, "disease_level")
  attr(stan_data, "target_tf") <- target_tf
  attr(stan_data, "target_tf_index") <- target_tf_index
  attr(stan_data, "target_interaction_model") <- isTRUE(target_interaction)
  attr(stan_data, "stage1_beta_summary") <- beta_summary
  attr(stan_data, "stage2_prior") <- stage2_prior

  stan_data
}

run_TF_stage2_directional_model <- function(
  analysis_object,
  stage1_fit_result,
  stan_file,
  direction_effect,
  beta_sd_floor,
  stage1_sd_multiplier,
  control_level,
  disease_level,
  confidence_min,
  confidence_max,
  chains,
  parallel_chains,
  iter_warmup,
  iter_sampling,
  seed,
  refresh,
  compute_loo,
  force_recompile,
  nuisance_prior_scale = 1,
  target_interaction = FALSE,
  target_interaction_sd = 0.5,
  ...
) {
  tf_model_require_pkg("cmdstanr")

  stan_data <- tf_prepare_stage2_stan_data(
    analysis_object = analysis_object,
    stage1_fit_result = stage1_fit_result,
    direction_effect = direction_effect,
    beta_sd_floor = beta_sd_floor,
    stage1_sd_multiplier = stage1_sd_multiplier,
    control_level = control_level,
    disease_level = disease_level,
    confidence_min = confidence_min,
    confidence_max = confidence_max,
    nuisance_prior_scale = nuisance_prior_scale,
    target_interaction = target_interaction,
    target_interaction_sd = target_interaction_sd
  )

  condition_model <- isTRUE(attr(stan_data, "condition_model"))
  use_batch_model <- isTRUE(attr(stan_data, "use_batch_model"))
  stan_file <- tf_select_stage2_stan_file(stan_file, condition_model, use_batch_model)

  model <- tf_get_cmdstan_model(
    stan_file,
    force_recompile = force_recompile,
    include_log_lik = isTRUE(compute_loo)
  )
  fit <- tf_fit_shared_stage2(
    stan_data = stan_data,
    model = model,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    refresh = refresh,
    ...
  )
  loo_result <- if (isTRUE(compute_loo)) tf_compute_loo(fit) else NULL

  list(
    fit = fit,
    loo = loo_result,
    stan_data = stan_data,
    stan_file = stan_file,
    target = attr(stan_data, "target"),
    feature_names = attr(stan_data, "feature_names"),
    cell_names = attr(stan_data, "cell_names"),
    use_batch_model = use_batch_model,
    condition_model = condition_model,
    target_interaction_model = isTRUE(attr(stan_data, "target_interaction_model")),
    target_tf = attr(stan_data, "target_tf"),
    target_tf_index = attr(stan_data, "target_tf_index"),
    control_level = attr(stan_data, "control_level"),
    disease_level = attr(stan_data, "disease_level"),
    W_names = attr(stan_data, "W_names"),
    stage1_beta_summary = attr(stan_data, "stage1_beta_summary"),
    stage2_prior = attr(stan_data, "stage2_prior")
  )
}
