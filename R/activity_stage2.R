#!/usr/bin/env Rscript

tf_em_stage2_input_stop <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

#' Convert retained Stage 2 MCMC models into the compact Stage 3 EM interface.
#'
#' Stage 2 target-TF posterior means are supplied only as EM initial values. For
#' each retained joint posterior draw, eta0 contains the intercept, condition,
#' batch, offset, and every confounder contribution, but deliberately excludes
#' the target-TF main effect and target-TF-by-condition interaction. The EM
#' M-step uses the posterior-informed Normal prior originally supplied to Stage
#' 2, while the target-TF expression is replaced by latent activity A.
build_TF_EM_stage2_input <- function(
  stage1_filtered_results,
  stage1_screening_input,
  nuisance_storage = c("all_draws", "posterior_mean"),
  nuisance_draw_ids = NULL,
  output_file = NULL
) {
  nuisance_storage <- match.arg(nuisance_storage)
  expected_model_version <- tf_shared_stage2_model_version()
  if (!is.list(stage1_filtered_results) || !length(stage1_filtered_results)) {
    tf_em_stage2_input_stop(
      "`stage1_filtered_results` must contain at least one retained gene."
    )
  }
  if (!is.list(stage1_screening_input) || !length(stage1_screening_input)) {
    tf_em_stage2_input_stop("`stage1_screening_input` must be a non-empty list.")
  }
  if (!identical(
    attr(stage1_filtered_results, "model_version"),
    expected_model_version
  ) || !identical(
    attr(stage1_screening_input, "model_version"),
    expected_model_version
  )) {
    tf_em_stage2_input_stop(
      "Stage 1 results and inputs must both use model version `%s`.",
      expected_model_version
    )
  }
  if (!identical(attr(stage1_filtered_results, "confounders_filtered"), FALSE)) {
    tf_em_stage2_input_stop(
      "The EM interface requires every Stage 1 confounder to be retained."
    )
  }
  result_stage <- attr(stage1_filtered_results, "pipeline_stage")
  if (!is.null(result_stage) && !identical(result_stage, "mcmc")) {
    tf_em_stage2_input_stop(
      "The EM interface must receive the filtered Stage 2 MCMC result."
    )
  }

  target_genes <- names(stage1_filtered_results)
  missing_input <- setdiff(target_genes, names(stage1_screening_input))
  if (length(missing_input)) {
    tf_em_stage2_input_stop(
      "Missing Stage 1 inputs for retained genes: %s.",
      paste(missing_input, collapse = ", ")
    )
  }

  first_gene <- target_genes[[1]]
  first_input <- stage1_screening_input[[first_gene]]
  if (!identical(first_input$status, "ready")) {
    tf_em_stage2_input_stop("Stage 1 input `%s` is not ready.", first_gene)
  }
  cell_names <- as.character(first_input$cell_names)
  condition <- as.integer(first_input$stan_data$condition)
  has_condition <- isTRUE(first_input$condition_model)
  target_tf_expression <- as.numeric(first_input$target_tf_expression)
  N <- length(cell_names)
  G <- length(target_genes)
  if (N < 1L || length(condition) != N ||
      length(target_tf_expression) != N ||
      any(!is.finite(target_tf_expression)) ||
      any(target_tf_expression < 0)) {
    tf_em_stage2_input_stop(
      "The target-TF expression vector must contain one finite nonnegative value per cell."
    )
  }
  if (has_condition && !identical(sort(unique(condition)), 0:1)) {
    tf_em_stage2_input_stop(
      "Condition-enabled activity inference requires both coded condition levels."
    )
  }
  if (!has_condition && any(condition != 0L)) {
    tf_em_stage2_input_stop(
      "Condition-free activity inference must use an all-zero condition vector."
    )
  }

  first_draw_bundle <-
    stage1_filtered_results[[first_gene]]$stage2_parameter_draws
  source_draw_count <- as.integer(first_draw_bundle$draw_count)
  if (length(source_draw_count) != 1L || is.na(source_draw_count) ||
      source_draw_count < 1L) {
    tf_em_stage2_input_stop("Invalid retained draw count for gene `%s`.", first_gene)
  }
  if (identical(nuisance_storage, "posterior_mean")) {
    if (!is.null(nuisance_draw_ids)) {
      tf_em_stage2_input_stop(
        "`nuisance_draw_ids` is only valid when `nuisance_storage = \"all_draws\"`."
      )
    }
    selected_draw_ids <- 0L
    S <- 1L
  } else {
    if (is.null(nuisance_draw_ids)) {
      selected_draw_ids <- seq_len(source_draw_count)
    } else {
      selected_draw_ids <- unique(as.integer(nuisance_draw_ids))
      if (!length(selected_draw_ids) || anyNA(selected_draw_ids) ||
          any(selected_draw_ids < 1L | selected_draw_ids > source_draw_count)) {
        tf_em_stage2_input_stop(
          "`nuisance_draw_ids` must be unique Stage 2 draw indices between 1 and %d.",
          source_draw_count
        )
      }
    }
    S <- length(selected_draw_ids)
  }

  Y <- matrix(NA_integer_, nrow = N, ncol = G,
              dimnames = list(cell_names, target_genes))
  eta0_draws <- array(
    NA_real_,
    dim = c(S, N, G),
    dimnames = list(
      stage1_draw = sprintf("draw_%d", seq_len(S)),
      cell = cell_names,
      target_gene = target_genes
    )
  )
  phi_draws <- matrix(
    NA_real_, nrow = S, ncol = G,
    dimnames = list(sprintf("draw_%d", seq_len(S)), target_genes)
  )
  beta_target_mean_init <- beta_target_delta_init <-
    beta_target_control_init <- beta_target_disease_init <-
      beta_prior_mean <- beta_prior_sd <- target_interaction_sd <-
        stats::setNames(rep(NA_real_, G), target_genes)
  source_draw_index <- setNames(vector("list", G), target_genes)
  confounders <- setNames(vector("list", G), target_genes)

  for (g in seq_along(target_genes)) {
    gene <- target_genes[[g]]
    model <- stage1_filtered_results[[gene]]
    input <- stage1_screening_input[[gene]]
    draws <- model$stage2_parameter_draws
    if (!identical(input$status, "ready")) {
      tf_em_stage2_input_stop("Stage 1 input `%s` is not ready.", gene)
    }
    if (!identical(draws$interface_version,
                   "tf_em_stage2_parameter_draws_v1") ||
        !identical(draws$prior_family, "normal") ||
        !identical(as.integer(draws$draw_count), source_draw_count)) {
      tf_em_stage2_input_stop(
        "Gene `%s` has an incompatible Stage 2 draw interface or draw count.",
        gene
      )
    }
    if (!identical(as.character(input$cell_names), cell_names) ||
        !identical(isTRUE(input$condition_model), has_condition) ||
        !identical(as.integer(input$stan_data$condition), condition)) {
      tf_em_stage2_input_stop(
        "Cell order or condition coding differs for retained gene `%s`.", gene
      )
    }
    if (!isTRUE(all.equal(
      as.numeric(input$target_tf_expression),
      target_tf_expression,
      tolerance = 1e-12,
      check.attributes = FALSE
    ))) {
      tf_em_stage2_input_stop(
        "Target-TF expression differs across retained gene inputs at `%s`.", gene
      )
    }

    feature_names <- as.character(input$feature_names)
    target_index <- as.integer(input$target_tf_index)
    if (!identical(colnames(draws$beta), feature_names) ||
        nrow(draws$beta) != source_draw_count ||
        ncol(draws$beta) != length(feature_names) ||
        length(target_index) != 1L || is.na(target_index) ||
        target_index < 1L || target_index > length(feature_names) ||
        !identical(as.integer(draws$target_tf_index), target_index)) {
      tf_em_stage2_input_stop(
        "Predictor mapping is inconsistent for retained gene `%s`.", gene
      )
    }
    use_batch <- if (is.null(input$use_batch)) {
      !is.null(input$stan_data$K_batch) && input$stan_data$K_batch > 0L
    } else {
      isTRUE(input$use_batch)
    }
    if (use_batch &&
        (!identical(colnames(draws$batch_effect), input$batch_levels) ||
         nrow(draws$batch_effect) != source_draw_count ||
         ncol(draws$batch_effect) != input$stan_data$K_batch)) {
      tf_em_stage2_input_stop(
        "Batch-effect mapping is inconsistent for retained gene `%s`.", gene
      )
    }
    if (!use_batch &&
        (nrow(draws$batch_effect) != source_draw_count ||
         ncol(draws$batch_effect) != 0L)) {
      tf_em_stage2_input_stop(
        "No-batch model `%s` unexpectedly contains batch-effect draws.", gene
      )
    }

    beta_baseline <- if (identical(nuisance_storage, "posterior_mean")) {
      matrix(
        colMeans(draws$beta),
        nrow = 1L,
        dimnames = list("posterior_mean", colnames(draws$beta))
      )
    } else {
      draws$beta[selected_draw_ids, , drop = FALSE]
    }
    beta_baseline[, target_index] <- 0
    eta0 <- beta_baseline %*% t(input$stan_data$X)
    alpha_nuisance <- if (identical(nuisance_storage, "posterior_mean")) {
      mean(draws$alpha)
    } else {
      as.numeric(draws$alpha[selected_draw_ids])
    }
    condition_nuisance <- if (identical(nuisance_storage, "posterior_mean")) {
      mean(draws$condition_effect)
    } else {
      as.numeric(draws$condition_effect[selected_draw_ids])
    }
    batch_nuisance <- if (use_batch) {
      if (identical(nuisance_storage, "posterior_mean")) {
        matrix(
          colMeans(draws$batch_effect),
          nrow = 1L,
          dimnames = list("posterior_mean", colnames(draws$batch_effect))
        )
      } else {
        draws$batch_effect[selected_draw_ids, , drop = FALSE]
      }
    } else {
      matrix(0, nrow = S, ncol = 0L)
    }
    phi_nuisance <- if (identical(nuisance_storage, "posterior_mean")) {
      mean(draws$phi)
    } else {
      as.numeric(draws$phi[selected_draw_ids])
    }
    eta0 <- sweep(eta0, 1L, alpha_nuisance, FUN = "+")
    if (has_condition) {
      eta0 <- eta0 + tcrossprod(
        condition_nuisance,
        as.numeric(input$stan_data$condition) - 0.5
      )
    }
    if (use_batch) {
      batch_contribution <- matrix(
        batch_nuisance[
          cbind(
            rep(seq_len(S), each = N),
            rep(as.integer(input$stan_data$batch), times = S)
          )
        ],
        nrow = S,
        ncol = N,
        byrow = TRUE
      )
      eta0 <- eta0 + batch_contribution
    }
    eta0 <- sweep(eta0, 2L, as.numeric(input$stan_data$log_offset), FUN = "+")
    if (any(!is.finite(eta0)) || length(phi_nuisance) != S ||
        any(!is.finite(phi_nuisance)) || any(phi_nuisance <= 0)) {
      tf_em_stage2_input_stop(
        "Non-finite eta0 or invalid phi draws for retained gene `%s`.", gene
      )
    }

    target_row <- model$target_tf
    target_name <- feature_names[[target_index]]
    Y[, g] <- as.integer(input$stan_data$Y)
    eta0_draws[, , g] <- eta0
    phi_draws[, g] <- phi_nuisance
    beta_target_mean_init[[g]] <- target_row$beta_mean[[1]]
    beta_target_delta_init[[g]] <- if (has_condition) {
      target_row$beta_delta_mean[[1]]
    } else {
      0
    }
    beta_target_control_init[[g]] <- if (has_condition) {
      target_row$beta_control_mean[[1]]
    } else {
      target_row$beta_mean[[1]]
    }
    beta_target_disease_init[[g]] <- if (has_condition) {
      target_row$beta_disease_mean[[1]]
    } else {
      target_row$beta_mean[[1]]
    }
    beta_prior_mean[[g]] <- draws$beta_prior_mean[[target_name]]
    beta_prior_sd[[g]] <- draws$beta_prior_sd[[target_name]]
    target_interaction_sd[[g]] <- draws$target_interaction_sd
    source_draw_index[[g]] <- if (identical(nuisance_storage, "posterior_mean")) {
      0L
    } else {
      as.integer(draws$source_draw_index[selected_draw_ids])
    }
    confounders[[g]] <- as.character(model$confounders$tf)
  }

  target_tf_center <- as.numeric(
    first_input$predictor_center[[first_input$target_tf_index]]
  )
  target_tf_scale <- as.numeric(
    first_input$predictor_scale[[first_input$target_tf_index]]
  )
  if (length(target_tf_center) != 1L || !is.finite(target_tf_center) ||
      length(target_tf_scale) != 1L || !is.finite(target_tf_scale) ||
      target_tf_scale <= 0) {
    tf_em_stage2_input_stop("Invalid target-TF centering or scaling metadata.")
  }

  result <- list(
    interface_version = "tf_em_stage2_input_v3",
    pipeline_stage = "em_input",
    model_version = expected_model_version,
    target_tf = first_input$target_tf,
    N = N,
    G = G,
    S = S,
    stage1_source_draw_count = source_draw_count,
    nuisance_source_draw_ids = as.integer(selected_draw_ids),
    nuisance_storage = nuisance_storage,
    cell_names = cell_names,
    target_genes = target_genes,
    has_condition = has_condition,
    condition = condition,
    control_level = first_input$control_level,
    disease_level = first_input$disease_level,
    target_tf_expression = target_tf_expression,
    activity_init = target_tf_expression,
    target_tf_center = target_tf_center,
    target_tf_scale = target_tf_scale,
    Y = Y,
    eta0_draws = eta0_draws,
    phi_draws = phi_draws,
    beta_target_mean_init = beta_target_mean_init,
    beta_target_delta_init = beta_target_delta_init,
    beta_target_control_init = beta_target_control_init,
    beta_target_disease_init = beta_target_disease_init,
    beta_prior_mean = beta_prior_mean,
    beta_prior_sd = beta_prior_sd,
    target_interaction_sd = target_interaction_sd,
    source_draw_index = source_draw_index,
    confounders = confounders,
    gene_weights = stats::setNames(rep(1 / G, G), target_genes)
  )
  attr(result, "prescreen_confounders_filtered") <-
    attr(stage1_filtered_results, "prescreen_confounders_filtered")
  attr(result, "prescreen_filter_counts") <-
    attr(stage1_filtered_results, "prescreen_filter_counts")
  class(result) <- c("TFEMStage2Input", "list")

  if (!is.null(output_file)) {
    output_file <- normalizePath(output_file, winslash = "/", mustWork = FALSE)
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    attr(result, "output_file") <- output_file
    saveRDS(result, output_file)
  }
  result
}
