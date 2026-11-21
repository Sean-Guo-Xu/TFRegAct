#!/usr/bin/env Rscript

tf_em_stage2_input_stop <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

# Form the Stage 3 coefficient prior from the FULL Stage 2 posterior. Compact
# nuisance draws may be thinned or selected independently and are never used
# here. Legacy v3 fits retain full posterior summaries, allowing exact reuse.
tf_stage2_posterior_normal_prior <- function(model, gene, has_condition) {
  scalar <- function(value) is.numeric(value) && length(value) == 1L && is.finite(value)
  valid_count <- function(value) scalar(value) && value >= 2 && value == trunc(value)
  fail <- function(detail) tf_em_stage2_input_stop(
    paste0("Cannot construct the Stage 3 posterior-Normal prior for `%s`: %s. ",
           "Rebuild the handoff from complete Stage 2 posterior summaries/draws; ",
           "refit Stage 2 only if its posterior is missing or degenerate. No variance floor or inflation is applied."),
    gene, detail)
  row <- model$target_tf
  required <- c("beta_mean", "beta_sd", if (has_condition) c("beta_delta_mean", "beta_delta_sd"))
  if (!is.data.frame(row) || nrow(row) != 1L || !all(required %in% names(row)) ||
      !all(vapply(row[required], scalar, logical(1)))) fail("full target-coefficient summaries are incomplete")
  summary_moments <- c(
    beta_prior_mean = row$beta_mean,
    beta_prior_variance = row$beta_sd ^ 2,
    beta_prior_sd = row$beta_sd,
    beta_delta_prior_mean = if (has_condition) row$beta_delta_mean else 0,
    beta_delta_prior_variance = if (has_condition) row$beta_delta_sd ^ 2 else 0,
    beta_delta_prior_sd = if (has_condition) row$beta_delta_sd else 0,
    beta_delta_prior_covariance = 0
  )
  if (has_condition) {
    if ("beta_mean_delta_covariance" %in% names(row) && scalar(row$beta_mean_delta_covariance)) {
      summary_moments[["beta_delta_prior_covariance"]] <- row$beta_mean_delta_covariance
    } else if (all(c("beta_control_sd", "beta_disease_sd") %in% names(row)) &&
               scalar(row$beta_control_sd) && scalar(row$beta_disease_sd)) {
      # Var(beta + delta/2) - Var(beta - delta/2) = 2 Cov(beta, delta).
      summary_moments[["beta_delta_prior_covariance"]] <- (row$beta_disease_sd ^ 2 - row$beta_control_sd ^ 2) / 2
    } else fail("full joint beta/delta covariance cannot be recovered")
  }
  moments <- summary_moments
  moment_source <- "full_stage2_posterior_summary"
  source_count <- model$diagnostics$posterior_draws
  if (!valid_count(source_count)) source_count <- NA_real_
  payload <- model$target_screening_draws
  if (!is.null(payload)) {
    count <- payload$draw_count
    if (!identical(payload$interface_version, "tf_target_screening_draws_v1") || !valid_count(count) ||
        !is.numeric(payload$overall) || length(payload$overall) != count || any(!is.finite(payload$overall))) {
      fail("full target-posterior draws are incomplete")
    }
    if (is.finite(source_count) && source_count != count) fail("full posterior draw counts disagree")
    source_count <- as.numeric(count)
    beta <- as.numeric(payload$overall)
    delta <- rep(0, count)
    delta_draws_available <- !has_condition
    if (has_condition) {
      if (all(vapply(payload[c("control", "disease")], function(x) {
        is.numeric(x) && length(x) == count && all(is.finite(x))
      }, logical(1))) && all(c("control", "disease") %in% names(payload))) {
        if (!is.null(payload$delta)) {
          delta <- as.numeric(payload$delta)
          delta_draws_available <- TRUE
        } else {
          joint <- model$stage2_parameter_draws
          complete_joint <- identical(as.integer(joint$source_draw_index), seq_len(count)) &&
            length(joint$beta_target_delta) == count && all(is.finite(joint$beta_target_delta))
          delta <- if (complete_joint) as.numeric(joint$beta_target_delta) else as.numeric(payload$disease - payload$control)
          delta_draws_available <- complete_joint
        }
      } else fail("full condition-specific target-posterior draws are incomplete")
      if (length(delta) != count || any(!is.finite(delta))) fail("full beta/delta draws are incomplete")
    }
    beta_sd <- stats::sd(beta)
    delta_sd <- if (has_condition) stats::sd(delta) else 0
    moments <- c(beta_prior_mean = mean(beta), beta_prior_variance = beta_sd ^ 2, beta_prior_sd = beta_sd,
                 beta_delta_prior_mean = if (has_condition) mean(delta) else 0,
                 beta_delta_prior_variance = delta_sd ^ 2, beta_delta_prior_sd = delta_sd,
                 beta_delta_prior_covariance = if (has_condition) stats::cov(beta, delta) else 0)
    # Older Stan CSVs round generated control/disease coefficients separately.
    # Their variance-difference covariance identity can differ at CSV precision
    # from an available full beta/delta vector. This is an audit tolerance only;
    # no posterior variance is inflated or floored.
    audit_tolerance <- if (has_condition && !("beta_mean_delta_covariance" %in% names(row))) 1e-6 else 1e-10
    if (!isTRUE(all.equal(unname(moments), unname(summary_moments), tolerance = audit_tolerance, scale = 1))) {
      fail("full posterior draws and stored summaries disagree")
    }
    if (delta_draws_available) {
      moment_source <- "full_stage2_target_draws"
    } else {
      # Complete delta means/SDs are already in the old full-posterior summary.
      # Keep them exactly if the old compact nuisance bundle omitted draws.
      moments <- summary_moments
      moment_source <- "full_stage2_posterior_summary"
    }
  }
  stored <- model$stage2_parameter_draws$target_posterior_moments
  if (!is.null(stored)) {
    stored_values <- unlist(stored[names(moments)], use.names = FALSE)
    if (!identical(stored$interface_version, "tf_stage2_target_posterior_moments_v1") ||
        !valid_count(stored$source_draw_count) || length(stored_values) != length(moments) ||
        !isTRUE(all.equal(as.numeric(stored_values), unname(moments), tolerance = 1e-8)) ||
        (is.finite(source_count) && stored$source_draw_count != source_count)) {
      fail("full posterior moment metadata disagree")
    }
    source_count <- as.numeric(stored$source_draw_count)
  }
  if (any(!is.finite(moments)) || moments[["beta_prior_sd"]] <= 0 ||
      moments[["beta_prior_variance"]] <= 0) fail("the beta posterior variance is not positive")
  if (has_condition) {
    joint <- matrix(c(moments[["beta_prior_variance"]], moments[["beta_delta_prior_covariance"]],
                      moments[["beta_delta_prior_covariance"]], moments[["beta_delta_prior_variance"]]), 2L)
    if (moments[["beta_delta_prior_sd"]] <= 0 ||
        inherits(try(chol(joint), silent = TRUE), "try-error")) fail("the full beta/delta posterior covariance is not positive definite")
  }
  c(as.list(moments), list(
    source_draw_count = source_count, moment_source = moment_source,
    source_count_status = if (is.finite(source_count)) "recorded_full_posterior_count" else "not_recorded_in_legacy_filtered_object",
    signature = paste(c("stage2_posterior_normal_approximation", gene,
                        format(source_count, digits = 17), format(moments, digits = 17)), collapse = "|")))
}

#' Convert retained Stage 2 MCMC models into the compact Stage 3 EM interface.
#'
#' Full Stage 2 target-coefficient posterior moments define a joint Normal prior
#' for Stage 3, without variance inflation. For
#' each retained joint posterior draw, eta0 contains the intercept, condition,
#' batch, offset, and every confounder contribution, but deliberately excludes
#' the target-TF main effect and target-TF-by-condition interaction. The EM
#' M-step uses this posterior Normal approximation,
#' while the target-TF expression is replaced by latent activity A.
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
      beta_prior_mean <- beta_prior_variance <- beta_prior_sd <- target_interaction_sd <-
        beta_delta_prior_mean <- beta_delta_prior_variance <- beta_delta_prior_sd <- beta_delta_prior_covariance <-
        stats::setNames(rep(NA_real_, G), target_genes)
  stage2_posterior_draw_count <- stats::setNames(rep(NA_real_, G), target_genes)
  stage2_posterior_prior_signature <- stats::setNames(rep(NA_character_, G), target_genes)
  stage2_posterior_prior_metadata <- stats::setNames(vector("list", G), target_genes)
  source_draw_index <- setNames(vector("list", G), target_genes)
  confounders <- setNames(vector("list", G), target_genes)
  knowledge_priors <- setNames(vector("list", G), target_genes)

  for (g in seq_along(target_genes)) {
    gene <- target_genes[[g]]
    model <- stage1_filtered_results[[gene]]
    input <- stage1_screening_input[[gene]]
    draws <- model$stage2_parameter_draws
    if (!identical(input$status, "ready")) {
      tf_em_stage2_input_stop("Stage 1 input `%s` is not ready.", gene)
    }
    if (!identical(draws$interface_version,
                   "tf_em_stage2_parameter_draws_v3") ||
        !identical(draws$prior_family, "normal") ||
        !identical(draws$beta_prior_family, "normal") ||
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
    expected_prior <- tf_stage2_knowledge_prior(input$stan_data)
    prior_fields <- c("beta_prior_mean", "beta_prior_variance", "beta_prior_sd")
    matching_prior <- vapply(prior_fields, function(field) {
      identical(names(draws[[field]]), feature_names) &&
        isTRUE(all.equal(as.numeric(draws[[field]]),
                         as.numeric(expected_prior[[field]]),
                         tolerance = 1e-12, check.attributes = FALSE))
    }, logical(1))
    if (!all(matching_prior)) {
      tf_em_stage2_input_stop(
        "Knowledge Normal prior metadata differ from the Stage 2 input for `%s`.", gene
      )
    }
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
    posterior_prior <- tf_stage2_posterior_normal_prior(model, gene, has_condition)
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
    beta_prior_mean[[g]] <- posterior_prior$beta_prior_mean
    beta_prior_variance[[g]] <- posterior_prior$beta_prior_variance
    beta_prior_sd[[g]] <- posterior_prior$beta_prior_sd
    beta_delta_prior_mean[[g]] <- posterior_prior$beta_delta_prior_mean
    beta_delta_prior_variance[[g]] <- posterior_prior$beta_delta_prior_variance
    beta_delta_prior_sd[[g]] <- posterior_prior$beta_delta_prior_sd
    beta_delta_prior_covariance[[g]] <- posterior_prior$beta_delta_prior_covariance
    stage2_posterior_draw_count[[g]] <- posterior_prior$source_draw_count
    stage2_posterior_prior_signature[[g]] <- posterior_prior$signature
    stage2_posterior_prior_metadata[[g]] <- posterior_prior
    knowledge_priors[[g]] <- expected_prior
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
    interface_version = "tf_em_stage2_input_v6",
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
    beta_prior_family = "normal",
    beta_prior_source = "stage2_posterior_normal_approximation",
    beta_prior_mean = beta_prior_mean,
    beta_prior_variance = beta_prior_variance,
    beta_prior_sd = beta_prior_sd,
    beta_delta_prior_mean = beta_delta_prior_mean,
    beta_delta_prior_variance = beta_delta_prior_variance,
    beta_delta_prior_sd = beta_delta_prior_sd,
    beta_delta_prior_covariance = beta_delta_prior_covariance,
    stage2_posterior_draw_count = stage2_posterior_draw_count,
    stage2_posterior_prior_signature = stage2_posterior_prior_signature,
    stage2_posterior_prior_metadata = stage2_posterior_prior_metadata,
    stage2_knowledge_priors = knowledge_priors,
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
