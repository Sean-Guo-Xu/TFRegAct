#!/usr/bin/env Rscript

tf_em_stop <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

.tf_em_script_dir <- local({
  source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  source_dir <- if (is.null(source_file) || !nzchar(as.character(source_file[[1]]))) {
    getwd()
  } else {
    dirname(normalizePath(source_file[[1]], winslash = "/", mustWork = FALSE))
  }
  function() source_dir
})

tf_em_gauss_legendre <- function(n) {
  n <- as.integer(n[[1]])
  if (is.na(n) || n < 5L) tf_em_stop("`quadrature_nodes` must be at least 5.")
  index <- seq_len(n - 1L)
  off_diagonal <- index / sqrt(4 * index ^ 2 - 1)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(index, index + 1L)] <- off_diagonal
  jacobi[cbind(index + 1L, index)] <- off_diagonal
  eig <- eigen(jacobi, symmetric = TRUE)
  order_index <- order(eig$values)
  list(
    nodes = eig$values[order_index],
    weights = 2 * eig$vectors[1L, order_index] ^ 2
  )
}

tf_em_log_phi_plus_exp_eta <- function(eta, log_phi) {
  pmax(eta, log_phi) + log1p(exp(-abs(eta - log_phi)))
}

tf_em_log_sum_exp2 <- function(a, b) {
  maximum <- pmax(a, b)
  maximum + log(exp(a - maximum) + exp(b - maximum))
}

tf_em_log_truncated_normal_density <- function(
  activity,
  mean,
  sd
) {
  stats::dnorm(activity, mean = mean, sd = sd, log = TRUE) -
    stats::pnorm(0, mean = mean, sd = sd, lower.tail = FALSE, log.p = TRUE)
}

tf_em_validate_input <- function(stage2_input) {
  if (is.character(stage2_input) && length(stage2_input) == 1L) {
    if (!file.exists(stage2_input)) {
      tf_em_stop("EM Stage2 input file not found: %s", stage2_input)
    }
    stage2_input <- readRDS(stage2_input)
  }
  if (!is.list(stage2_input) ||
      !identical(stage2_input$interface_version, "tf_em_stage2_input_v4")) {
    tf_em_stop("`stage2_input` must use interface `tf_em_stage2_input_v4`.")
  }
  required <- c(
    "N", "G", "S", "Y", "eta0_draws", "phi_draws", "has_condition", "condition",
    "target_tf_expression", "target_tf_center", "target_tf_scale",
    "beta_target_mean_init", "beta_target_delta_init", "beta_prior_family",
    "beta_lambda_positive", "beta_lambda_negative", "target_interaction_sd",
    "cell_names", "target_genes"
  )
  missing <- setdiff(required, names(stage2_input))
  if (length(missing)) {
    tf_em_stop("EM Stage2 input is missing: %s.", paste(missing, collapse = ", "))
  }
  if (!identical(stage2_input$beta_prior_family, "laplace")) {
    tf_em_stop("EM requires the Stage2 knowledge-based asymmetric Laplace beta prior.")
  }
  N <- as.integer(stage2_input$N)
  G <- as.integer(stage2_input$G)
  S <- as.integer(stage2_input$S)
  if (length(N) != 1L || length(G) != 1L || length(S) != 1L ||
      anyNA(c(N, G, S)) || any(c(N, G, S) < 1L)) {
    tf_em_stop("N, G and S must be positive integer dimensions.")
  }
  if (!identical(dim(stage2_input$Y), c(N, G)) ||
      !identical(dim(stage2_input$eta0_draws), c(S, N, G)) ||
      !identical(dim(stage2_input$phi_draws), c(S, G))) {
    tf_em_stop("Y, eta0_draws, or phi_draws has dimensions inconsistent with N/G/S.")
  }
  if (length(stage2_input$condition) != N ||
      !all(stage2_input$condition %in% c(0L, 1L)) ||
      length(stage2_input$target_tf_expression) != N ||
      any(!is.finite(stage2_input$target_tf_expression)) ||
      any(stage2_input$target_tf_expression < 0) ||
      !is.finite(stage2_input$target_tf_scale) ||
      stage2_input$target_tf_scale <= 0) {
    tf_em_stop("Invalid condition or target-TF activity-anchor metadata.")
  }
  has_condition <- isTRUE(stage2_input$has_condition)
  if (has_condition && !identical(sort(unique(stage2_input$condition)), 0:1)) {
    tf_em_stop("Condition-enabled EM requires both control and disease cells.")
  }
  if (!has_condition && any(stage2_input$condition != 0L)) {
    tf_em_stop("Condition-free EM must use an all-zero condition vector.")
  }
  if (any(stage2_input$Y < 0) || any(!is.finite(stage2_input$eta0_draws)) ||
      any(!is.finite(stage2_input$phi_draws)) || any(stage2_input$phi_draws <= 0)) {
    tf_em_stop("Y, eta0_draws, and phi_draws contain invalid values.")
  }
  gene_vectors <- c(
    "beta_target_mean_init", "beta_target_delta_init", "beta_lambda_positive",
    "beta_lambda_negative", "target_interaction_sd"
  )
  for (field in gene_vectors) {
    value <- stage2_input[[field]]
    if (length(value) != G || any(!is.finite(value))) {
      tf_em_stop("`%s` must contain one finite value per target gene.", field)
    }
  }
  if (any(stage2_input$beta_lambda_positive <= 0) ||
      any(stage2_input$beta_lambda_negative <= 0) ||
      any(stage2_input$target_interaction_sd <= 0)) {
    tf_em_stop("Beta prior rates and interaction scales must be positive.")
  }
  if (length(stage2_input$cell_names) != N ||
      length(stage2_input$target_genes) != G ||
      anyDuplicated(stage2_input$cell_names) ||
      anyDuplicated(stage2_input$target_genes)) {
    tf_em_stop("Cell and target-gene names must be complete and unique.")
  }
  stage2_input
}

tf_em_cell_log_kernel <- function(
  activity,
  y,
  eta0,
  phi,
  beta_condition,
  activity_anchor,
  sigma_activity,
  target_center,
  target_scale
) {
  activity_std <- (activity - target_center) / target_scale
  eta <- eta0 + beta_condition * activity_std
  log_denominator <- tf_em_log_phi_plus_exp_eta(eta, log(phi))
  tf_em_log_truncated_normal_density(
    activity,
    mean = activity_anchor,
    sd = sigma_activity
  ) +
    sum(y * eta - (y + phi) * log_denominator)
}

tf_em_cell_log_likelihood <- function(
  activity,
  y,
  eta0,
  phi,
  beta_condition,
  target_center,
  target_scale
) {
  activity_std <- (activity - target_center) / target_scale
  eta <- eta0 + beta_condition * activity_std
  log_denominator <- tf_em_log_phi_plus_exp_eta(eta, log(phi))
  sum(y * eta - (y + phi) * log_denominator)
}

tf_em_cell_derivatives <- function(
  activity,
  y,
  eta0,
  phi,
  beta_condition,
  activity_anchor,
  sigma_activity,
  target_center,
  target_scale
) {
  eta <- eta0 + beta_condition * (activity - target_center) / target_scale
  probability_mu <- stats::plogis(eta - log(phi))
  score_eta <- y - (y + phi) * probability_mu
  first <- -(activity - activity_anchor) / sigma_activity ^ 2 +
    sum((beta_condition / target_scale) * score_eta)
  second <- -1 / sigma_activity ^ 2 - sum(
    (beta_condition / target_scale) ^ 2 *
      (y + phi) * probability_mu * (1 - probability_mu)
  )
  c(first = first, second = second)
}

tf_em_cell_posterior <- function(
  y,
  eta0,
  phi,
  beta_condition,
  activity_anchor,
  sigma_activity,
  target_center,
  target_scale,
  active_prior_probability,
  quadrature,
  tail_log_drop = 30,
  max_upper_expansions = 12L
) {
  # Hard zero-expression anchoring has an exact spike posterior. Avoiding the
  # otherwise unnecessary mode search and quadrature is important because
  # zero-inflation is substantial in single-cell normalized expression.
  if (active_prior_probability <= 0) {
    log_spike_marginal <- tf_em_cell_log_likelihood(
      0, y, eta0, phi, beta_condition, target_center, target_scale
    )
    quadrature_count <- length(quadrature$nodes)
    return(list(
      nodes = c(0, rep(0, quadrature_count)),
      weights = c(1, rep(0, quadrature_count)),
      mean = 0,
      variance = 0,
      q05 = 0,
      median = 0,
      q95 = 0,
      mode = 0,
      slab_mode = 0,
      active_probability = 0,
      given_active_mean = 0,
      log_slab_marginal = -Inf,
      log_spike_marginal = log_spike_marginal,
      log_normalizer = log_spike_marginal,
      upper = 0,
      upper_expansions = 0L,
      tail_log_drop_achieved = Inf
    ))
  }
  derivative_zero <- tf_em_cell_derivatives(
    0, y, eta0, phi, beta_condition, activity_anchor, sigma_activity,
    target_center, target_scale
  )[["first"]]
  upper <- max(
    activity_anchor + 4 * sigma_activity,
    target_center + 4 * target_scale,
    4 * sigma_activity,
    1e-6
  )
  upper_expansions <- 0L
  if (derivative_zero <= 0) {
    mode <- 0
  } else {
    derivative_upper <- tf_em_cell_derivatives(
      upper, y, eta0, phi, beta_condition, activity_anchor, sigma_activity,
      target_center, target_scale
    )[["first"]]
    while (derivative_upper > 0 && upper_expansions < max_upper_expansions) {
      upper <- upper * 2
      upper_expansions <- upper_expansions + 1L
      derivative_upper <- tf_em_cell_derivatives(
        upper, y, eta0, phi, beta_condition, activity_anchor, sigma_activity,
        target_center, target_scale
      )[["first"]]
    }
    if (derivative_upper > 0) {
      tf_em_stop("Failed to bracket a cell-specific activity posterior mode.")
    }
    mode <- stats::uniroot(
      function(a) tf_em_cell_derivatives(
        a, y, eta0, phi, beta_condition, activity_anchor, sigma_activity,
        target_center, target_scale
      )[["first"]],
      lower = 0,
      upper = upper,
      tol = 1e-8
    )$root
  }
  mode_terms <- tf_em_cell_derivatives(
    mode, y, eta0, phi, beta_condition, activity_anchor, sigma_activity,
    target_center, target_scale
  )
  local_sd <- sqrt(1 / max(-mode_terms[["second"]], 1e-12))
  upper <- max(
    upper,
    mode + 8 * local_sd,
    activity_anchor + 8 * sigma_activity,
    8 * sigma_activity,
    1e-6
  )
  log_mode <- tf_em_cell_log_kernel(
    mode, y, eta0, phi, beta_condition, activity_anchor, sigma_activity,
    target_center, target_scale
  )
  log_upper <- tf_em_cell_log_kernel(
    upper, y, eta0, phi, beta_condition, activity_anchor, sigma_activity,
    target_center, target_scale
  )
  while (log_upper > log_mode - tail_log_drop &&
         upper_expansions < max_upper_expansions) {
    upper <- upper * 1.75
    upper_expansions <- upper_expansions + 1L
    log_upper <- tf_em_cell_log_kernel(
      upper, y, eta0, phi, beta_condition, activity_anchor, sigma_activity,
      target_center, target_scale
    )
  }

  nodes <- 0.5 * upper * (quadrature$nodes + 1)
  integration_weights <- 0.5 * upper * quadrature$weights
  log_kernel <- vapply(nodes, tf_em_cell_log_kernel, numeric(1),
    y = y, eta0 = eta0, phi = phi, beta_condition = beta_condition,
    activity_anchor = activity_anchor, sigma_activity = sigma_activity,
    target_center = target_center, target_scale = target_scale
  )
  log_weight <- log(integration_weights) + log_kernel
  max_log_weight <- max(log_weight)
  slab_posterior_weight <- exp(log_weight - max_log_weight)
  slab_posterior_weight <- slab_posterior_weight / sum(slab_posterior_weight)
  log_slab_marginal <- max_log_weight +
    log(sum(exp(log_weight - max_log_weight)))
  log_spike_marginal <- tf_em_cell_log_likelihood(
    0, y, eta0, phi, beta_condition, target_center, target_scale
  )
  log_active_joint <- log(active_prior_probability) + log_slab_marginal
  log_inactive_joint <- log1p(-active_prior_probability) + log_spike_marginal
  log_marginal <- tf_em_log_sum_exp2(log_active_joint, log_inactive_joint)
  active_probability <- exp(log_active_joint - log_marginal)
  inactive_probability <- 1 - active_probability
  posterior_weight <- c(
    inactive_probability,
    active_probability * slab_posterior_weight
  )
  posterior_nodes <- c(0, nodes)
  posterior_mean <- sum(posterior_weight * posterior_nodes)
  posterior_variance <- sum(
    posterior_weight * (posterior_nodes - posterior_mean) ^ 2
  )
  slab_mean <- sum(slab_posterior_weight * nodes)
  cumulative_weight <- cumsum(posterior_weight)
  weighted_quantile <- function(probability) {
    posterior_nodes[which(cumulative_weight >= probability)[[1]]]
  }
  list(
    nodes = posterior_nodes,
    weights = posterior_weight,
    mean = posterior_mean,
    variance = max(posterior_variance, 0),
    q05 = weighted_quantile(0.05),
    median = weighted_quantile(0.50),
    q95 = weighted_quantile(0.95),
    mode = if (inactive_probability >= 0.5) 0 else mode,
    slab_mode = mode,
    active_probability = active_probability,
    given_active_mean = slab_mean,
    log_slab_marginal = log_slab_marginal,
    log_spike_marginal = log_spike_marginal,
    log_normalizer = log_marginal,
    upper = upper,
    upper_expansions = upper_expansions,
    tail_log_drop_achieved = log_mode - log_upper
  )
}

tf_em_estep <- function(common, eta0, phi, beta_mean, beta_delta, settings) {
  N <- common$N
  K <- settings$quadrature_nodes + 1L
  nodes <- weights <- matrix(NA_real_, N, K)
  posterior_mean <- posterior_variance <- posterior_q05 <- posterior_median <-
    posterior_q95 <- posterior_mode <- active_probability <-
      given_active_mean <- log_normalizer <- numeric(N)
  upper_expansions <- integer(N)
  tail_drop <- numeric(N)
  for (i in seq_len(N)) {
    beta_condition <- beta_mean + common$condition_shift[[i]] * beta_delta
    posterior <- tf_em_cell_posterior(
      y = common$Y[i, ],
      eta0 = eta0[i, ],
      phi = phi,
      beta_condition = beta_condition,
      activity_anchor = common$activity_anchor[[i]],
      sigma_activity = common$sigma_activity,
      target_center = common$target_center,
      target_scale = common$target_scale,
      active_prior_probability = common$active_prior_probability[[i]],
      quadrature = common$quadrature,
      tail_log_drop = settings$tail_log_drop,
      max_upper_expansions = settings$max_upper_expansions
    )
    nodes[i, ] <- posterior$nodes
    weights[i, ] <- posterior$weights
    posterior_mean[[i]] <- posterior$mean
    posterior_variance[[i]] <- posterior$variance
    posterior_q05[[i]] <- posterior$q05
    posterior_median[[i]] <- posterior$median
    posterior_q95[[i]] <- posterior$q95
    posterior_mode[[i]] <- posterior$mode
    active_probability[[i]] <- posterior$active_probability
    given_active_mean[[i]] <- posterior$given_active_mean
    log_normalizer[[i]] <- posterior$log_normalizer
    upper_expansions[[i]] <- posterior$upper_expansions
    tail_drop[[i]] <- posterior$tail_log_drop_achieved
  }
  list(
    nodes = nodes,
    weights = weights,
    mean = posterior_mean,
    variance = posterior_variance,
    q05 = posterior_q05,
    median = posterior_median,
    q95 = posterior_q95,
    mode = posterior_mode,
    active_probability = active_probability,
    given_active_mean = given_active_mean,
    log_normalizer = log_normalizer,
    upper_expansions = upper_expansions,
    tail_drop = tail_drop
  )
}

tf_em_gene_expected_log_posterior <- function(
  parameters,
  gene_index,
  common,
  eta0_gene,
  phi_gene,
  estep,
  return_information = FALSE,
  beta_branch = NULL
) {
  beta_by_cell <- parameters[[1]] +
    common$condition_shift * parameters[[2]]
  activity_std <- (estep$nodes - common$target_center) / common$target_scale
  eta <- sweep(activity_std, 1L, beta_by_cell, FUN = "*")
  eta <- sweep(eta, 1L, eta0_gene, FUN = "+")
  log_denominator <- tf_em_log_phi_plus_exp_eta(eta, log(phi_gene))
  y_gene <- common$Y[, gene_index]
  log_likelihood <- sweep(eta, 1L, y_gene, FUN = "*") -
    sweep(log_denominator, 1L, y_gene + phi_gene, FUN = "*")
  expected_log_likelihood <- sum(estep$weights * log_likelihood)
  lambda_positive <- common$beta_lambda_positive[[gene_index]]
  lambda_negative <- common$beta_lambda_negative[[gene_index]]
  log_prior <- if (parameters[[1]] >= 0) {
    -lambda_positive * parameters[[1]]
  } else {
    lambda_negative * parameters[[1]]
  }
  log_prior <- log_prior - 0.5 * (
    parameters[[2]] / common$target_interaction_sd[[gene_index]]
  ) ^ 2
  # At zero the prior has two one-sided derivatives. The constrained optimizer
  # supplies its branch explicitly; there is no smooth derivative at the kink.
  positive_branch <- if (is.null(beta_branch)) parameters[[1]] >= 0 else {
    match.arg(beta_branch, c("positive", "negative")) == "positive"
  }
  prior_gradient <- if (positive_branch) -lambda_positive else lambda_negative

  probability_mu <- stats::plogis(eta - log(phi_gene))
  score_eta <- sweep(
    probability_mu,
    1L,
    y_gene + phi_gene,
    FUN = "*"
  )
  score_eta <- sweep(-score_eta, 1L, y_gene, FUN = "+")
  common_gradient <- estep$weights * score_eta * activity_std
  condition_shift <- common$condition_shift
  gradient <- c(
    sum(common_gradient) + prior_gradient,
    sum(sweep(common_gradient, 1L, condition_shift, FUN = "*")) -
      parameters[[2]] / common$target_interaction_sd[[gene_index]] ^ 2
  )
  result <- list(
    value = expected_log_likelihood + log_prior,
    gradient = gradient
  )
  if (isTRUE(return_information)) {
    curvature <- estep$weights *
      sweep(probability_mu * (1 - probability_mu), 1L,
            y_gene + phi_gene, FUN = "*") * activity_std ^ 2
    information <- matrix(c(
      sum(curvature),
      sum(sweep(curvature, 1L, condition_shift, FUN = "*")),
      sum(sweep(curvature, 1L, condition_shift, FUN = "*")),
      sum(sweep(curvature, 1L, condition_shift ^ 2, FUN = "*")) +
        1 / common$target_interaction_sd[[gene_index]] ^ 2
    ), 2L, 2L)
    result$information <- information
    result$at_prior_kink <- parameters[[1]] == 0
  }
  result
}

tf_em_mstep <- function(common, eta0, phi, estep, beta_mean, beta_delta, settings) {
  G <- common$G
  updated_mean <- updated_delta <- numeric(G)
  covariance <- array(NA_real_, dim = c(2L, 2L, G))
  convergence <- integer(G)
  objective <- numeric(G)
  at_prior_kink <- logical(G)
  branch_convergence <- matrix(NA_integer_, G, 2L,
    dimnames = list(common$target_genes, c("positive", "negative")))
  for (g in seq_len(G)) {
    branch_fits <- lapply(c("positive", "negative"), function(branch) {
      positive <- branch == "positive"
      start_beta <- if (positive) max(beta_mean[[g]], 0) else min(beta_mean[[g]], 0)
      start <- if (common$has_condition) c(start_beta, beta_delta[[g]]) else start_beta
      unpack <- function(parameters) {
        if (common$has_condition) parameters else c(parameters[[1]], 0)
      }
      objective_function <- function(parameters) {
        -tf_em_gene_expected_log_posterior(
          unpack(parameters), g, common, eta0[, g], phi[[g]], estep,
          beta_branch = branch
        )$value
      }
      gradient_function <- function(parameters) {
        gradient <- tf_em_gene_expected_log_posterior(
          unpack(parameters), g, common, eta0[, g], phi[[g]], estep,
          beta_branch = branch
        )$gradient
        -gradient[seq_along(parameters)]
      }
      stats::optim(
        par = start,
        fn = objective_function,
        gr = gradient_function,
        method = "L-BFGS-B",
        lower = c(if (positive) 0 else -Inf, if (common$has_condition) -Inf),
        upper = c(if (positive) Inf else 0, if (common$has_condition) Inf),
        control = list(
          maxit = settings$mstep_maxit,
          factr = max(1, settings$mstep_reltol / .Machine$double.eps),
          pgtol = settings$mstep_reltol
        )
      )
    })
    branch_convergence[g, ] <- vapply(branch_fits, `[[`, integer(1), "convergence")
    # The expected NB log likelihood is concave in (beta, delta), and each
    # Laplace half-axis is smooth. Comparing both constrained solutions includes
    # the kink and avoids pretending that the full prior is differentiable.
    fit <- branch_fits[[which.min(vapply(branch_fits, `[[`, numeric(1), "value"))]]
    fitted_parameters <- if (common$has_condition) fit$par else c(fit$par[[1]], 0)
    if (any(!is.finite(fit$par)) || !is.finite(fit$value)) {
      tf_em_stop("Non-finite M-step result for gene `%s`.", common$target_genes[[g]])
    }
    updated_mean[[g]] <- fitted_parameters[[1]]
    updated_delta[[g]] <- fitted_parameters[[2]]
    at_prior_kink[[g]] <- fitted_parameters[[1]] == 0
    convergence[[g]] <- fit$convergence
    objective[[g]] <- -fit$value
    information <- tf_em_gene_expected_log_posterior(
      fitted_parameters, g, common, eta0[, g], phi[[g]], estep,
      return_information = TRUE
    )$information
    if (at_prior_kink[[g]]) {
      # A Gaussian based on the one-sided curvature cannot represent the cusp.
      # Leave uncertainty unavailable here; the conditional MAP remains valid.
      covariance[, , g] <- NA_real_
    } else if (common$has_condition) {
      covariance[, , g] <- tryCatch(
        solve(information),
        error = function(e) {
          eigen_result <- eigen(information, symmetric = TRUE)
          eigen_result$vectors %*%
            diag(1 / pmax(eigen_result$values, 1e-10), 2L) %*%
            t(eigen_result$vectors)
        }
      )
    } else {
      covariance[, , g] <- matrix(0, 2L, 2L)
      covariance[1L, 1L, g] <- 1 / pmax(information[1L, 1L], 1e-10)
    }
  }
  list(
    beta_mean = updated_mean,
    beta_delta = updated_delta,
    covariance = covariance,
    convergence = convergence,
    branch_convergence = branch_convergence,
    at_prior_kink = at_prior_kink,
    uncertainty_method = ifelse(at_prior_kink, "unavailable_at_prior_kink",
      "local_complete_data_normal_not_full_posterior"),
    objective = objective
  )
}

tf_em_beta_log_prior <- function(beta_mean, beta_delta, common) {
  result <- sum(ifelse(beta_mean >= 0,
    -common$beta_lambda_positive * beta_mean,
    common$beta_lambda_negative * beta_mean))
  if (common$has_condition) {
    result <- result +
      sum(-0.5 * (beta_delta / common$target_interaction_sd) ^ 2)
  }
  result
}

tf_em_relative_change <- function(new, old, epsilon = 1e-8) {
  max(abs(new - old) / (abs(old) + epsilon))
}

tf_em_scaled_beta_change <- function(
  new_mean,
  old_mean,
  new_delta,
  old_delta,
  common
) {
  prior_sd <- sqrt(1 / common$beta_lambda_positive ^ 2 +
                     1 / common$beta_lambda_negative ^ 2)
  changes <- abs(new_mean - old_mean) / prior_sd
  if (common$has_condition) {
    changes <- c(
      changes,
      abs(new_delta - old_delta) / common$target_interaction_sd
    )
  }
  max(changes)
}

tf_em_scaled_activity_change <- function(new, old, common) {
  max(abs(new - old)) / common$target_scale
}

tf_em_sample_discrete_rows <- function(nodes, weights) {
  vapply(seq_len(nrow(nodes)), function(i) {
    sample(nodes[i, ], size = 1L, prob = weights[i, ])
  }, numeric(1))
}

tf_em_sample_bivariate_normal <- function(mean, covariance) {
  covariance <- 0.5 * (covariance + t(covariance))
  eig <- eigen(covariance, symmetric = TRUE)
  mean + eig$vectors %*% (sqrt(pmax(eig$values, 0)) * stats::rnorm(2L))
}

tf_em_sampler_ess <- function(values) {
  count <- length(values)
  if (count < 4L || !all(is.finite(values))) return(NA_real_)
  if (stats::var(values) == 0) return(NA_real_)
  correlation <- as.numeric(stats::acf(
    values, lag.max = count - 1L, plot = FALSE, demean = TRUE
  )$acf)
  pair_count <- length(correlation) %/% 2L
  paired <- correlation[seq.int(1L, 2L * pair_count, by = 2L)] +
    correlation[seq.int(2L, 2L * pair_count, by = 2L)]
  first_nonpositive <- which(paired <= 0)
  if (length(first_nonpositive)) paired <- utils::head(paired, first_nonpositive[[1]] - 1L)
  if (!length(paired)) return(as.numeric(count))
  autocorrelation_time <- max(1, -1 + 2 * sum(cummin(paired)))
  min(as.numeric(count), count / autocorrelation_time)
}

tf_em_gibbs_fit_one_draw <- function(job, common, settings) {
  started <- Sys.time()
  tryCatch({
    set.seed(settings$seed + as.integer(job$draw_id))
    samples <- tf_activity_gibbs_ars_cpp(common, job$eta0, job$phi, settings)
    activity <- samples$activity_samples
    beta <- samples$beta_samples
    delta <- samples$delta_samples
    active <- activity > 0
    active_count <- rowSums(active)
    activity_quantiles <- t(apply(activity, 1L, stats::quantile,
      probs = c(0.05, 0.5, 0.95), names = FALSE))
    coefficient_summary <- function(values) {
      quantiles <- t(apply(values, 2L, stats::quantile,
        probs = c(0.05, 0.5, 0.95), names = FALSE))
      cbind(mean = colMeans(values), sd = apply(values, 2L, stats::sd),
        q05 = quantiles[, 1L], median = quantiles[, 2L], q95 = quantiles[, 3L])
    }
    beta_covariance <- array(0, dim = c(2L, 2L, common$G))
    for (g in seq_len(common$G)) {
      beta_covariance[, , g] <- stats::cov(cbind(beta[, g], delta[, g]))
    }
    selected <- sample.int(settings$iter_sampling, 1L)
    beta_ess <- apply(beta, 2L, tf_em_sampler_ess)
    delta_ess <- if (common$has_condition) apply(delta, 2L, tf_em_sampler_ess) else {
      rep(NA_real_, common$G)
    }
    difference_samples <- if (common$has_condition) {
      colMeans(activity[common$condition == 1L, , drop = FALSE]) -
        colMeans(activity[common$condition == 0L, , drop = FALSE])
    } else numeric(0)
    list(
      status = "ok", draw_id = as.integer(job$draw_id), converged = NA,
      iterations = as.integer(settings$iter_warmup + settings$iter_sampling),
      beta_mean = colMeans(beta), beta_delta = colMeans(delta),
      beta_control = colMeans(beta - 0.5 * delta),
      beta_disease = colMeans(beta + 0.5 * delta),
      beta_covariance = beta_covariance,
      beta_at_prior_kink = rep(FALSE, common$G),
      beta_uncertainty_method = rep("gibbs_ars_empirical_posterior", common$G),
      mstep_branch_convergence = matrix(NA_integer_, common$G, 2L,
        dimnames = list(common$target_genes, c("positive", "negative"))),
      beta_estimation = "posterior_mean",
      beta_posterior_sample_mean = beta[selected, ],
      beta_posterior_sample_delta = delta[selected, ],
      coefficient_posterior_summary = list(
        beta_mean = coefficient_summary(beta), beta_delta = coefficient_summary(delta),
        beta_control = coefficient_summary(beta - 0.5 * delta),
        beta_disease = coefficient_summary(beta + 0.5 * delta)
      ),
      activity_mean = rowMeans(activity),
      activity_variance = apply(activity, 1L, stats::var),
      activity_q05 = activity_quantiles[, 1L], activity_median = activity_quantiles[, 2L],
      activity_q95 = activity_quantiles[, 3L], activity_mode = rep(NA_real_, common$N),
      activity_active_probability = rowMeans(active),
      activity_given_active_mean = ifelse(active_count > 0L,
        rowSums(activity) / pmax(active_count, 1L), NA_real_),
      activity_active_draw_count = active_count,
      activity_posterior_sample = activity[, selected],
      activity_condition_difference_samples = difference_samples,
      final_observed_objective = NA_real_,
      final_joint_log_kernel = utils::tail(samples$log_joint_trace, 1L),
      monotonicity_violations = NA_integer_, max_upper_expansions = NA_integer_,
      min_tail_log_drop = NA_real_,
      sampling_diagnostics = list(
        chains = 1L, rhat = NA_real_,
        convergence_assessment = "One chain per conditional nuisance target; convergence not established",
        ess_method = "Geyer initial positive monotone sequence for retained scalar draws",
        beta_ess = stats::setNames(beta_ess, common$target_genes),
        delta_ess = stats::setNames(delta_ess, common$target_genes),
        mean_activity_ess = tf_em_sampler_ess(samples$activity_mean_trace),
        active_count_ess = tf_em_sampler_ess(samples$active_count_trace),
        ars_proposals = samples$ars_proposals, ars_rejections = samples$ars_rejections,
        scale_acceptance_warmup = if (settings$iter_warmup > 0L) {
          samples$scale_accepted_warmup / settings$iter_warmup
        } else NA_real_,
        scale_acceptance_sampling = samples$scale_accepted_sampling / settings$iter_sampling
      ),
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
      trace = if (isTRUE(settings$save_traces)) data.frame(
        iteration = seq_along(samples$log_joint_trace),
        warmup = seq_along(samples$log_joint_trace) <= settings$iter_warmup,
        joint_log_kernel = samples$log_joint_trace
      ) else NULL,
      error = NA_character_
    )
  }, error = function(error) list(status = "error", draw_id = as.integer(job$draw_id),
    error = conditionMessage(error),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))))
}

tf_em_fit_one_draw <- function(job, common, settings) {
  if (identical(settings$algorithm, "gibbs_ars")) {
    return(tf_em_gibbs_fit_one_draw(job, common, settings))
  }
  started <- Sys.time()
  tryCatch({
    eta0 <- job$eta0
    phi <- job$phi
    beta_mean <- common$beta_mean_init
    beta_delta <- if (common$has_condition) common$beta_delta_init else rep(0, common$G)
    previous_activity <- common$activity_anchor
    previous_objective <- -Inf
    trace_rows <- vector("list", settings$max_iter)
    converged <- FALSE
    last_mstep <- NULL

    for (iteration in seq_len(settings$max_iter)) {
      estep <- tf_em_estep(
        common, eta0, phi, beta_mean, beta_delta, settings
      )
      observed_objective <- sum(estep$log_normalizer) +
        tf_em_beta_log_prior(beta_mean, beta_delta, common)
      mstep <- tf_em_mstep(
        common, eta0, phi, estep, beta_mean, beta_delta, settings
      )
      beta_change <- tf_em_scaled_beta_change(
        new_mean = mstep$beta_mean,
        old_mean = beta_mean,
        new_delta = mstep$beta_delta,
        old_delta = beta_delta,
        common = common
      )
      activity_change <- tf_em_scaled_activity_change(
        estep$mean, previous_activity, common
      )
      objective_change <- if (is.finite(previous_objective)) {
        observed_objective - previous_objective
      } else {
        NA_real_
      }
      objective_relative_change <- if (is.finite(previous_objective)) {
        abs(objective_change) / (1 + abs(previous_objective))
      } else {
        NA_real_
      }
      trace_rows[[iteration]] <- data.frame(
        iteration = iteration,
        observed_objective = observed_objective,
        objective_change = objective_change,
        objective_relative_change = objective_relative_change,
        beta_relative_change = beta_change,
        activity_relative_change = activity_change,
        mstep_nonzero_convergence = sum(mstep$convergence != 0L),
        max_upper_expansions = max(estep$upper_expansions),
        min_tail_log_drop = min(estep$tail_drop),
        stringsAsFactors = FALSE
      )
      beta_mean <- mstep$beta_mean
      beta_delta <- mstep$beta_delta
      last_mstep <- mstep
      if (iteration >= settings$min_iter &&
          beta_change < settings$beta_tolerance &&
          activity_change < settings$activity_tolerance &&
          is.finite(objective_relative_change) &&
          objective_relative_change < settings$objective_tolerance) {
        converged <- TRUE
        break
      }
      previous_activity <- estep$mean
      previous_objective <- observed_objective
    }

    final_estep <- tf_em_estep(
      common, eta0, phi, beta_mean, beta_delta, settings
    )
    final_objective <- sum(final_estep$log_normalizer) +
      tf_em_beta_log_prior(beta_mean, beta_delta, common)
    set.seed(settings$seed + as.integer(job$draw_id))
    activity_sample <- tf_em_sample_discrete_rows(
      final_estep$nodes, final_estep$weights
    )
    beta_posterior_sample <- matrix(NA_real_, common$G, 2L)
    for (g in seq_len(common$G)) {
      if (last_mstep$at_prior_kink[[g]]) next
      if (common$has_condition) {
        beta_posterior_sample[g, ] <- tf_em_sample_bivariate_normal(
          c(beta_mean[[g]], beta_delta[[g]]),
          last_mstep$covariance[, , g]
        )
      } else {
        beta_posterior_sample[g, ] <- c(
          stats::rnorm(
            1L,
            mean = beta_mean[[g]],
            sd = sqrt(pmax(last_mstep$covariance[1L, 1L, g], 0))
          ),
          0
        )
      }
    }
    trace <- do.call(rbind, trace_rows[!vapply(trace_rows, is.null, logical(1))])
    list(
      status = "ok",
      draw_id = as.integer(job$draw_id),
      converged = converged,
      iterations = nrow(trace),
      beta_mean = beta_mean,
      beta_delta = beta_delta,
      beta_control = beta_mean - 0.5 * beta_delta,
      beta_disease = beta_mean + 0.5 * beta_delta,
      beta_covariance = last_mstep$covariance,
      beta_at_prior_kink = last_mstep$at_prior_kink,
      beta_uncertainty_method = last_mstep$uncertainty_method,
      mstep_branch_convergence = last_mstep$branch_convergence,
      beta_estimation = "conditional_MAP",
      beta_posterior_sample_mean = beta_posterior_sample[, 1L],
      beta_posterior_sample_delta = beta_posterior_sample[, 2L],
      activity_mean = final_estep$mean,
      activity_variance = final_estep$variance,
      activity_median = final_estep$median,
      activity_q05 = final_estep$q05,
      activity_q95 = final_estep$q95,
      activity_mode = final_estep$mode,
      activity_active_probability = final_estep$active_probability,
      activity_given_active_mean = final_estep$given_active_mean,
      activity_posterior_sample = activity_sample,
      final_observed_objective = final_objective,
      monotonicity_violations = sum(
        trace$objective_change < -settings$monotonicity_tolerance,
        na.rm = TRUE
      ),
      max_upper_expansions = max(final_estep$upper_expansions),
      min_tail_log_drop = min(final_estep$tail_drop),
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
      trace = if (isTRUE(settings$save_traces)) trace else NULL,
      error = NA_character_
    )
  }, error = function(e) {
    list(
      status = "error",
      draw_id = as.integer(job$draw_id),
      converged = FALSE,
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
      error = conditionMessage(e)
    )
  })
}

tf_em_worker_entry <- function(job) {
  tf_em_fit_one_draw(job, .tf_em_worker_common, .tf_em_worker_settings)
}

tf_em_summarize_gibbs_results <- function(summary, draw_results, common, settings) {
  draw_results <- draw_results[order(as.integer(names(draw_results)))]
  count <- length(draw_results)
  for (prefix in c("beta_mean", "beta_delta", "beta_control", "beta_disease")) {
    components <- lapply(draw_results, function(result) {
      result$coefficient_posterior_summary[[prefix]]
    })
    component_mean <- vapply(components, function(value) value[, "mean"], numeric(common$G))
    component_variance <- vapply(components, function(value) value[, "sd"] ^ 2, numeric(common$G))
    mean <- rowMeans(matrix(component_mean, nrow = common$G))
    variance <- rowMeans(matrix(component_variance + component_mean ^ 2, nrow = common$G)) - mean ^ 2
    summary$beta_summary[[paste0(prefix, "_mean")]] <- mean
    summary$beta_summary[[paste0(prefix, "_sd")]] <- sqrt(pmax(variance, 0))
    if (count == 1L) {
      for (quantile in c("q05", "median", "q95")) {
        summary$beta_summary[[paste0(prefix, "_", quantile)]] <- components[[1]][, quantile]
      }
    }
  }
  if (!common$has_condition) {
    columns <- grep("^beta_(control|disease)_", names(summary$beta_summary), value = TRUE)
    summary$beta_summary[columns] <- NA_real_
  }
  summary$beta_summary$beta_mean_map <- NA_real_
  summary$beta_summary$beta_delta_map <- NA_real_
  summary$beta_summary$beta_at_prior_kink <- FALSE
  summary$beta_summary$beta_estimation <- "posterior_mean"
  summary$beta_summary$beta_uncertainty_method <- "gibbs_ars_empirical_posterior"
  summary$beta_mean_posterior_mean_by_draw <- summary$beta_mean_map_by_draw
  summary$beta_delta_posterior_mean_by_draw <- summary$beta_delta_map_by_draw
  summary$beta_mean_posterior_sample_by_draw <- summary$beta_mean_local_normal_sample_by_draw
  summary$beta_delta_posterior_sample_by_draw <- summary$beta_delta_local_normal_sample_by_draw
  # Preserve historical fields without labelling posterior means as MAP or
  # actual Gibbs draws as local Gaussian approximations.
  unavailable <- matrix(NA_real_, count, common$G,
    dimnames = dimnames(summary$beta_mean_map_by_draw))
  for (field in c("beta_mean_map_by_draw", "beta_delta_map_by_draw",
      "beta_mean_local_normal_sample_by_draw", "beta_delta_local_normal_sample_by_draw",
      "beta_mean_laplace_sample_by_draw", "beta_delta_laplace_sample_by_draw")) {
    summary[[field]] <- unavailable
  }
  summary$summary_method <- if (count == 1L) {
    "conditional_gibbs_ars_posterior"
  } else "equal_weight_nuisance_mixture_of_gibbs_ars_posteriors"
  summary$posterior_quantile_method <- if (count == 1L) {
    "empirical_all_retained_gibbs_draws"
  } else {
    "approximate_mixture_quantiles_from_one_joint_gibbs_draw_per_nuisance_component"
  }
  summary$sampling_diagnostics <- stats::setNames(
    lapply(draw_results, `[[`, "sampling_diagnostics"), names(draw_results))
  summary$convergence_summary$algorithm <- "gibbs_ars"
  summary$convergence_summary$iter_warmup <- settings$iter_warmup
  summary$convergence_summary$iter_sampling <- settings$iter_sampling
  summary$convergence_summary$scale_acceptance_sampling <- vapply(draw_results,
    function(result) result$sampling_diagnostics$scale_acceptance_sampling, numeric(1))
  if (common$has_condition) {
    sampled <- summary$activity_condition_difference_summary$contrast_type ==
      "one_activity_posterior_sample_per_cell"
    summary$activity_condition_difference_summary$contrast_type[sampled] <-
      "one_joint_gibbs_draw_per_nuisance_component"
    if (count == 1L) {
      values <- draw_results[[1]]$activity_condition_difference_samples
      quantiles <- stats::quantile(values, c(0.05, 0.5, 0.95), names = FALSE)
      summary$activity_condition_difference_summary[sampled,
        c("mean", "sd", "q05", "median", "q95", "probability_disease_higher")] <-
          list(mean(values), stats::sd(values), quantiles[[1]], quantiles[[2]],
            quantiles[[3]], mean(values > 0))
      summary$activity_condition_difference_summary$contrast_type[sampled] <-
        "all_retained_joint_gibbs_draws"
    }
  }
  summary
}

tf_em_summarize_draw_results <- function(draw_results, common, settings) {
  if (!length(draw_results) ||
      any(!vapply(draw_results, function(x) identical(x$status, "ok"), logical(1)))) {
    tf_em_stop("All requested nuisance-draw EM fits must be successful before summarizing.")
  }
  draw_ids <- as.integer(names(draw_results))
  order_index <- order(draw_ids)
  draw_ids <- draw_ids[order_index]
  draw_results <- draw_results[order_index]
  R <- length(draw_results)
  N <- common$N
  G <- common$G
  bind_numeric <- function(field, rows) {
    matrix(
      unlist(lapply(draw_results, `[[`, field), use.names = FALSE),
      nrow = rows,
      ncol = R
    )
  }
  activity_mean <- bind_numeric("activity_mean", N)
  activity_variance <- bind_numeric("activity_variance", N)
  activity_active_probability <- bind_numeric(
    "activity_active_probability", N
  )
  activity_given_active_mean <- bind_numeric(
    "activity_given_active_mean", N
  )
  activity_sample <- bind_numeric("activity_posterior_sample", N)
  beta_mean_map <- t(bind_numeric("beta_mean", G))
  beta_delta_map <- t(bind_numeric("beta_delta", G))
  beta_mean_sample <- t(bind_numeric("beta_posterior_sample_mean", G))
  beta_delta_sample <- t(bind_numeric("beta_posterior_sample_delta", G))
  colnames(beta_mean_map) <- colnames(beta_delta_map) <-
    colnames(beta_mean_sample) <- colnames(beta_delta_sample) <- common$target_genes
  rownames(beta_mean_map) <- rownames(beta_delta_map) <-
    rownames(beta_mean_sample) <- rownames(beta_delta_sample) <- draw_ids
  total_mean <- rowMeans(activity_mean)
  total_variance <- rowMeans(activity_variance + activity_mean ^ 2) - total_mean ^ 2
  total_active_probability <- rowMeans(activity_active_probability)
  weighted_active_mean <- activity_active_probability * activity_given_active_mean
  weighted_active_mean[activity_active_probability == 0] <- 0
  total_given_active_mean <- rowMeans(weighted_active_mean) /
    pmax(total_active_probability, .Machine$double.eps)
  total_given_active_mean[total_active_probability == 0] <- NA_real_
  activity_quantiles <- if (R == 1L) {
    cbind(
      draw_results[[1]]$activity_q05,
      draw_results[[1]]$activity_median,
      draw_results[[1]]$activity_q95
    )
  } else {
    t(apply(
      activity_sample, 1L, stats::quantile,
      probs = c(0.05, 0.5, 0.95), names = FALSE, na.rm = TRUE
    ))
  }
  activity_summary <- data.frame(
    cell = common$cell_names,
    condition = if (common$has_condition) common$condition else NA_integer_,
    condition_label = if (common$has_condition) {
      ifelse(
        common$condition == 0L, common$control_level, common$disease_level
      )
    } else {
      rep("Overall", N)
    },
    target_tf_expression = common$activity_anchor,
    activity_prior_probability = common$active_prior_probability,
    activity_active_probability = total_active_probability,
    activity_given_active_mean = total_given_active_mean,
    activity_mean = total_mean,
    activity_sd = sqrt(pmax(total_variance, 0)),
    activity_q05 = activity_quantiles[, 1L],
    activity_median = activity_quantiles[, 2L],
    activity_q95 = activity_quantiles[, 3L],
    stringsAsFactors = FALSE
  )
  summarize_beta <- function(mean_sample, delta_sample) {
    control_sample <- mean_sample - 0.5 * delta_sample
    disease_sample <- mean_sample + 0.5 * delta_sample
    summarize_matrix <- function(x, prefix) {
      quantile_matrix <- t(apply(
        x, 2L, function(value) {
          # A missing kink approximation must not silently drop that nuisance
          # component and report intervals from the remaining draws only.
          if (any(!is.finite(value))) return(rep(NA_real_, 3L))
          stats::quantile(value, probs = c(0.05, 0.5, 0.95), names = FALSE)
        }
      ))
      data.frame(
        mean = colMeans(x),
        sd = apply(x, 2L, stats::sd),
        q05 = quantile_matrix[, 1L],
        median = quantile_matrix[, 2L],
        q95 = quantile_matrix[, 3L],
        row.names = common$target_genes,
        check.names = FALSE
      ) |>
        stats::setNames(paste0(prefix, c("_mean", "_sd", "_q05", "_median", "_q95")))
    }
    cbind(
      target_gene = common$target_genes,
      summarize_matrix(mean_sample, "beta_mean"),
      summarize_matrix(delta_sample, "beta_delta"),
      summarize_matrix(control_sample, "beta_control"),
      summarize_matrix(disease_sample, "beta_disease"),
      row.names = NULL
    )
  }
  summarize_single_beta <- function(draw_result) {
    covariance <- draw_result$beta_covariance
    covariance_mean_delta <- covariance[1L, 2L, ]
    variance_mean <- covariance[1L, 1L, ]
    variance_delta <- covariance[2L, 2L, ]
    variance_control <- variance_mean + 0.25 * variance_delta -
      covariance_mean_delta
    variance_disease <- variance_mean + 0.25 * variance_delta +
      covariance_mean_delta
    summarize_normal <- function(estimate, variance, prefix) {
      standard_error <- sqrt(pmax(variance, 0))
      z05 <- stats::qnorm(0.05)
      data.frame(
        mean = estimate,
        sd = standard_error,
        q05 = estimate + z05 * standard_error,
        median = estimate,
        q95 = estimate - z05 * standard_error,
        row.names = common$target_genes,
        check.names = FALSE
      ) |>
        stats::setNames(paste0(prefix, c("_mean", "_sd", "_q05", "_median", "_q95")))
    }
    cbind(
      target_gene = common$target_genes,
      summarize_normal(draw_result$beta_mean, variance_mean, "beta_mean"),
      summarize_normal(draw_result$beta_delta, variance_delta, "beta_delta"),
      summarize_normal(draw_result$beta_control, variance_control, "beta_control"),
      summarize_normal(draw_result$beta_disease, variance_disease, "beta_disease"),
      row.names = NULL
    )
  }
  convergence_summary <- data.frame(
    draw_id = draw_ids,
    converged = vapply(draw_results, `[[`, logical(1), "converged"),
    iterations = vapply(draw_results, `[[`, integer(1), "iterations"),
    monotonicity_violations = vapply(
      draw_results, `[[`, integer(1), "monotonicity_violations"
    ),
    final_observed_objective = vapply(
      draw_results, `[[`, numeric(1), "final_observed_objective"
    ),
    elapsed_seconds = vapply(draw_results, `[[`, numeric(1), "elapsed_seconds"),
    stringsAsFactors = FALSE
  )
  activity_condition_difference_by_draw <- data.frame(
    draw_id = integer(0),
    control_activity_conditional_mean = numeric(0),
    disease_activity_conditional_mean = numeric(0),
    disease_minus_control_conditional_mean = numeric(0),
    control_activity_posterior_sample = numeric(0),
    disease_activity_posterior_sample = numeric(0),
    disease_minus_control_posterior_sample = numeric(0),
    stringsAsFactors = FALSE
  )
  summarize_difference <- function(values, contrast_type) {
    quantiles <- stats::quantile(
      values, probs = c(0.05, 0.5, 0.95), names = FALSE, na.rm = TRUE
    )
    data.frame(
      contrast = sprintf("%s minus %s", common$disease_level, common$control_level),
      contrast_type = contrast_type,
      nuisance_draws = R,
      mean = mean(values),
      sd = if (R > 1L) stats::sd(values) else 0,
      q05 = quantiles[[1]],
      median = quantiles[[2]],
      q95 = quantiles[[3]],
      probability_disease_higher = mean(values > 0),
      stringsAsFactors = FALSE
    )
  }
  activity_condition_difference_summary <- data.frame(
    contrast = character(0),
    contrast_type = character(0),
    nuisance_draws = integer(0),
    mean = numeric(0),
    sd = numeric(0),
    q05 = numeric(0),
    median = numeric(0),
    q95 = numeric(0),
    probability_disease_higher = numeric(0),
    stringsAsFactors = FALSE
  )
  if (common$has_condition) {
    control_cells <- common$condition == 0L
    disease_cells <- common$condition == 1L
    # Each column corresponds to one independently fitted nuisance posterior
    # draw. The sample contrast also propagates cell-level activity uncertainty.
    activity_condition_difference_by_draw <- data.frame(
      draw_id = draw_ids,
      control_activity_conditional_mean = colMeans(
        activity_mean[control_cells, , drop = FALSE]
      ),
      disease_activity_conditional_mean = colMeans(
        activity_mean[disease_cells, , drop = FALSE]
      ),
      disease_minus_control_conditional_mean = colMeans(
        activity_mean[disease_cells, , drop = FALSE]
      ) - colMeans(activity_mean[control_cells, , drop = FALSE]),
      control_activity_posterior_sample = colMeans(
        activity_sample[control_cells, , drop = FALSE]
      ),
      disease_activity_posterior_sample = colMeans(
        activity_sample[disease_cells, , drop = FALSE]
      ),
      disease_minus_control_posterior_sample = colMeans(
        activity_sample[disease_cells, , drop = FALSE]
      ) - colMeans(activity_sample[control_cells, , drop = FALSE]),
      stringsAsFactors = FALSE
    )
    activity_condition_difference_summary <- rbind(
      summarize_difference(
        activity_condition_difference_by_draw$disease_minus_control_conditional_mean,
        "conditional_posterior_mean"
      ),
      summarize_difference(
        activity_condition_difference_by_draw$disease_minus_control_posterior_sample,
        "one_activity_posterior_sample_per_cell"
      )
    )
  }
  beta_summary <- if (R == 1L) {
    summarize_single_beta(draw_results[[1]])
  } else {
    summarize_beta(beta_mean_sample, beta_delta_sample)
  }
  kink_by_draw <- t(bind_numeric("beta_at_prior_kink", G))
  beta_summary$beta_mean_map <- colMeans(beta_mean_map)
  beta_summary$beta_delta_map <- colMeans(beta_delta_map)
  beta_summary$beta_at_prior_kink <- colSums(kink_by_draw) > 0
  beta_summary$beta_estimation <- "conditional_MAP"
  beta_summary$beta_uncertainty_method <- ifelse(beta_summary$beta_at_prior_kink,
    "unavailable_at_prior_kink", "local_complete_data_normal_not_full_posterior")
  approximation_columns <- grep(
    "^beta_(mean|delta|control|disease)_(mean|sd|q05|median|q95)$",
    names(beta_summary), value = TRUE)
  beta_summary[beta_summary$beta_at_prior_kink, approximation_columns] <- NA_real_
  if (!common$has_condition) {
    condition_specific_columns <- grep(
      "^beta_(control|disease)_", names(beta_summary), value = TRUE
    )
    beta_summary[condition_specific_columns] <- NA_real_
  }
  result <- list(
    activity_summary = activity_summary,
    beta_summary = beta_summary,
    summary_method = if (R == 1L) {
      "conditional_activity_and_local_complete_data_beta_approximation"
    } else {
      "equal_weight_nuisance_mixture"
    },
    convergence_summary = convergence_summary,
    activity_condition_difference_by_draw = activity_condition_difference_by_draw,
    activity_condition_difference_summary = activity_condition_difference_summary,
    activity_conditional_mean = activity_mean,
    activity_conditional_variance = activity_variance,
    activity_active_probability_by_draw = activity_active_probability,
    activity_given_active_mean_by_draw = activity_given_active_mean,
    activity_mixture_sample = activity_sample,
    beta_mean_map_by_draw = beta_mean_map,
    beta_delta_map_by_draw = beta_delta_map,
    beta_at_prior_kink_by_draw = kink_by_draw,
    beta_mean_local_normal_sample_by_draw = beta_mean_sample,
    beta_delta_local_normal_sample_by_draw = beta_delta_sample,
    # Legacy field names retained as aliases; these are local normal heuristics,
    # not samples from the asymmetric-Laplace joint posterior.
    beta_mean_laplace_sample_by_draw = beta_mean_sample,
    beta_delta_laplace_sample_by_draw = beta_delta_sample
  )
  if (identical(settings$algorithm, "gibbs_ars")) {
    result <- tf_em_summarize_gibbs_results(result, draw_results, common, settings)
  }
  result
}

#' Run spike-and-slab latent-activity inference over Stage 2 nuisance draws.
#'
#' Each cell has an exact inactive state A = 0 and a positive truncated-normal
#' slab anchored to target-TF normalized expression. The negative-binomial
#' likelihood of all retained downstream genes updates the posterior active
#' probability. `active_prior_zero` is the tunable gate prior for cells without
#' detected TF expression. Cells with detected TF expression use a fixed prior
#' active probability of 0.9. The gate is soft in both groups, so zero TF
#' expression strongly downweights activity without forcing it to exactly zero.
#'
#' With the default `nuisance_draw_count = 0`, eta0 and phi use their Stage 2
#' posterior means in one fit. A positive count
#' reproducibly selects that many posterior draws without replacement
#' and fits them independently; their activity posteriors form
#' an equal-weight mixture. Explicit draw IDs remain available as an override.
#' The default Gibbs sampler uses ARS for the log-concave coefficient and
#' positive-activity conditionals, an auxiliary slab variable for the zero
#' spike, and a joint coefficient/activity scale Metropolis move. Warmup and
#' retained iterations each default to 500. Stage 2 target-beta posterior means
#' initialize the sampler and the network asymmetric-Laplace prior is retained.
#' One chain is run per nuisance component, so convergence is reported as
#' unknown, with within-chain ESS and proposal acceptance diagnostics.
#' Set `algorithm = "quadrature_em"` for the conditional-MAP EM implementation
#' and its local complete-data Gaussian uncertainty away from beta = 0.
#' Confounder effects are fixed within each nuisance
#' draw, so their Stage 2 posterior uncertainty is propagated only when more
#' than one nuisance draw is requested.
run_TF_EM_latent_activity <- function(
  stage2_input,
  nuisance_draw_count = 0L,
  nuisance_draw_ids = NULL,
  cores = 4L,
  kappa = 1,
  active_prior_zero = 0.2,
  expression_zero_tolerance = 0,
  quadrature_nodes = 21L,
  max_iter = 30L,
  min_iter = 2L,
  beta_tolerance = 1e-3,
  activity_tolerance = 1e-3,
  objective_tolerance = 1e-8,
  monotonicity_tolerance = 1e-4,
  tail_log_drop = 30,
  max_upper_expansions = 12L,
  mstep_maxit = 100L,
  mstep_reltol = 1e-8,
  seed = 123L,
  checkpoint_every = NULL,
  checkpoint_file = NULL,
  output_file = NULL,
  resume = TRUE,
  retry_errors = TRUE,
  save_traces = FALSE,
  algorithm = "gibbs_ars",
  iter_warmup = 500L,
  iter_sampling = 500L,
  scale_move_sd = 0.1
) {
  stage2_input <- tf_em_validate_input(stage2_input)
  algorithm <- match.arg(algorithm, c("gibbs_ars", "quadrature_em"))
  iter_warmup <- as.integer(iter_warmup[[1]])
  iter_sampling <- as.integer(iter_sampling[[1]])
  scale_move_sd <- as.numeric(scale_move_sd[[1]])
  if (is.na(iter_warmup) || iter_warmup < 0L ||
      is.na(iter_sampling) || iter_sampling < 2L ||
      !is.finite(scale_move_sd) || scale_move_sd <= 0 ||
      as.double(iter_warmup) + as.double(iter_sampling) > .Machine$integer.max) {
    tf_em_stop("Sampling needs nonnegative warmup, at least two retained draws and a positive scale-move SD.")
  }
  cores <- as.integer(cores[[1]])
  quadrature_nodes <- as.integer(quadrature_nodes[[1]])
  max_iter <- as.integer(max_iter[[1]])
  min_iter <- as.integer(min_iter[[1]])
  max_upper_expansions <- as.integer(max_upper_expansions[[1]])
  mstep_maxit <- as.integer(mstep_maxit[[1]])
  seed <- as.integer(seed[[1]])
  nuisance_draw_count <- as.integer(nuisance_draw_count[[1]])
  if (identical(stage2_input$nuisance_storage, "posterior_mean") &&
      nuisance_draw_count > 0L) {
    tf_em_stop(
      paste0(
        "This Stage2 input stores nuisance posterior means only. ",
        "Use `nuisance_draw_count = 0`, or rebuild the input with ",
        "`nuisance_storage = \"all_draws\"`."
      )
    )
  }
  if (is.na(nuisance_draw_count) || nuisance_draw_count < 0L ||
      nuisance_draw_count > stage2_input$S) {
    tf_em_stop(
      "`nuisance_draw_count` must be between 0 and the available draw count (%d).",
      stage2_input$S
    )
  }
  if (is.null(nuisance_draw_ids)) {
    if (nuisance_draw_count == 0L) {
      nuisance_draw_ids <- 0L
      nuisance_mode <- "posterior_mean"
    } else {
      set.seed(seed)
      nuisance_draw_ids <- sort(sample.int(
        stage2_input$S, size = nuisance_draw_count, replace = FALSE
      ))
      nuisance_mode <- "posterior_draws"
    }
  } else {
    nuisance_draw_ids <- unique(as.integer(nuisance_draw_ids))
    if (!length(nuisance_draw_ids) || anyNA(nuisance_draw_ids) ||
        any(nuisance_draw_ids < 1L | nuisance_draw_ids > stage2_input$S)) {
      tf_em_stop("`nuisance_draw_ids` contains invalid draw indices.")
    }
    nuisance_draw_count <- length(nuisance_draw_ids)
    nuisance_mode <- "explicit_posterior_draws"
  }
  numeric_positive <- c(
    kappa, beta_tolerance, activity_tolerance, objective_tolerance,
    monotonicity_tolerance, tail_log_drop, mstep_reltol
  )
  if (is.na(cores) || cores < 1L || is.na(max_iter) || max_iter < 1L ||
      is.na(min_iter) || min_iter < 1L || min_iter > max_iter ||
      is.na(max_upper_expansions) || max_upper_expansions < 1L ||
      is.na(mstep_maxit) || mstep_maxit < 1L || is.na(seed) ||
      any(!is.finite(numeric_positive)) || any(numeric_positive <= 0)) {
    tf_em_stop("Invalid EM numerical configuration.")
  }
  active_prior_positive <- 0.9
  gate_configuration <- c(active_prior_zero, expression_zero_tolerance)
  if (any(!is.finite(gate_configuration)) ||
      active_prior_zero <= 0 || active_prior_zero >= active_prior_positive ||
      expression_zero_tolerance < 0) {
    tf_em_stop(
      paste0(
        "`active_prior_zero` must lie strictly between 0 and the fixed ",
        "positive-expression prior (0.9), and `expression_zero_tolerance` ",
        "must be nonnegative."
      )
    )
  }
  quadrature <- if (identical(algorithm, "quadrature_em")) {
    tf_em_gauss_legendre(quadrature_nodes)
  } else NULL
  sigma_activity <- as.numeric(kappa) * stage2_input$target_tf_scale
  expression_is_zero <-
    stage2_input$target_tf_expression <= expression_zero_tolerance
  active_prior_probability <- ifelse(
    expression_is_zero,
    as.numeric(active_prior_zero),
    active_prior_positive
  )
  common <- list(
    N = stage2_input$N,
    G = stage2_input$G,
    Y = stage2_input$Y,
    has_condition = isTRUE(stage2_input$has_condition),
    condition = as.integer(stage2_input$condition),
    condition_shift = if (isTRUE(stage2_input$has_condition)) {
      as.numeric(stage2_input$condition) - 0.5
    } else {
      rep(0, stage2_input$N)
    },
    activity_anchor = as.numeric(stage2_input$target_tf_expression),
    active_prior_probability = active_prior_probability,
    sigma_activity = sigma_activity,
    target_center = as.numeric(stage2_input$target_tf_center),
    target_scale = as.numeric(stage2_input$target_tf_scale),
    beta_mean_init = as.numeric(stage2_input$beta_target_mean_init),
    beta_delta_init = as.numeric(stage2_input$beta_target_delta_init),
    beta_prior_family = stage2_input$beta_prior_family,
    beta_lambda_positive = as.numeric(stage2_input$beta_lambda_positive),
    beta_lambda_negative = as.numeric(stage2_input$beta_lambda_negative),
    target_interaction_sd = as.numeric(stage2_input$target_interaction_sd),
    cell_names = stage2_input$cell_names,
    target_genes = stage2_input$target_genes,
    control_level = stage2_input$control_level,
    disease_level = stage2_input$disease_level,
    quadrature = quadrature
  )
  settings <- list(
    algorithm = algorithm,
    algorithm_version = if (identical(algorithm, "gibbs_ars")) {
      "knowledge_asymmetric_laplace_gibbs_ars_v1"
    } else "knowledge_asymmetric_laplace_em_v2",
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    scale_move_sd = scale_move_sd,
    beta_prior_family = "laplace",
    beta_uncertainty = if (identical(algorithm, "gibbs_ars")) {
      "empirical_conditional_joint_posterior"
    } else "local_complete_data_normal_not_full_posterior; unavailable_at_prior_kink",
    quadrature_nodes = quadrature_nodes,
    max_iter = max_iter,
    min_iter = min_iter,
    beta_tolerance = as.numeric(beta_tolerance),
    activity_tolerance = as.numeric(activity_tolerance),
    objective_tolerance = as.numeric(objective_tolerance),
    monotonicity_tolerance = as.numeric(monotonicity_tolerance),
    tail_log_drop = as.numeric(tail_log_drop),
    max_upper_expansions = max_upper_expansions,
    mstep_maxit = mstep_maxit,
    mstep_reltol = as.numeric(mstep_reltol),
    seed = seed,
    save_traces = isTRUE(save_traces),
    kappa = as.numeric(kappa),
    sigma_activity = sigma_activity,
    active_prior_zero = as.numeric(active_prior_zero),
    active_prior_positive = active_prior_positive,
    expression_zero_tolerance = as.numeric(expression_zero_tolerance),
    hard_zero_expression = FALSE,
    activity_model = "spike_and_truncated_normal_slab"
  )
  # Dimensions and draw IDs do not identify the actual nuisance values. Hash
  # the complete handoff so a changed Stage 2 fit cannot reuse stale Stage 3.
  input_hash_file <- tempfile("tfregact-stage3-input-", fileext = ".rds")
  on.exit(unlink(input_hash_file), add = TRUE)
  saveRDS(stage2_input, input_hash_file, compress = FALSE, version = 2)
  input_hash <- unname(tools::md5sum(input_hash_file))
  signature <- paste(
    input_hash,
    stage2_input$interface_version,
    stage2_input$target_tf,
    stage2_input$N,
    stage2_input$G,
    stage2_input$S,
    stage2_input$beta_prior_family,
    paste(format(stage2_input$beta_lambda_positive, digits = 17), collapse = ","),
    paste(format(stage2_input$beta_lambda_negative, digits = 17), collapse = ","),
    paste(stage2_input$target_genes, collapse = ","),
    nuisance_mode,
    nuisance_draw_count,
    paste(nuisance_draw_ids, collapse = ","),
    paste(stage2_input$nuisance_source_draw_ids, collapse = ","),
    paste(unlist(settings, use.names = TRUE), collapse = ","),
    sep = "::"
  )
  if (is.null(output_file)) {
    output_file <- file.path(
      .tf_em_script_dir(), "adjustment_output",
      sprintf("%s_EM_latent_activity_fit.rds", stage2_input$target_tf)
    )
  }
  output_file <- normalizePath(output_file, winslash = "/", mustWork = FALSE)
  if (is.null(checkpoint_file)) {
    checkpoint_file <- sub("\\.rds$", "_checkpoint.rds", output_file)
    if (identical(checkpoint_file, output_file)) {
      checkpoint_file <- paste0(output_file, "_checkpoint.rds")
    }
  }
  checkpoint_file <- normalizePath(
    checkpoint_file, winslash = "/", mustWork = FALSE
  )
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(checkpoint_file), recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(resume) && file.exists(output_file)) {
    existing <- readRDS(output_file)
    if (identical(attr(existing, "em_signature"), signature) &&
        isTRUE(attr(existing, "complete"))) {
      return(existing)
    }
  }
  draw_results <- setNames(vector("list", length(nuisance_draw_ids)), nuisance_draw_ids)
  if (isTRUE(resume) && file.exists(checkpoint_file)) {
    checkpoint <- readRDS(checkpoint_file)
    if (identical(checkpoint$em_signature, signature)) {
      reusable <- intersect(names(checkpoint$draw_results), names(draw_results))
      if (isTRUE(retry_errors)) {
        reusable <- reusable[vapply(
          checkpoint$draw_results[reusable],
          function(x) is.list(x) && identical(x$status, "ok"), logical(1)
        )]
      }
      draw_results[reusable] <- checkpoint$draw_results[reusable]
    }
  }
  pending <- nuisance_draw_ids[vapply(draw_results, is.null, logical(1))]
  worker_count <- min(cores, max(1L, length(pending)))
  if (is.null(checkpoint_every)) checkpoint_every <- worker_count
  checkpoint_every <- as.integer(checkpoint_every[[1]])
  if (is.na(checkpoint_every) || checkpoint_every < 1L) {
    tf_em_stop("`checkpoint_every` must be a positive integer.")
  }
  save_checkpoint <- function() {
    saveRDS(list(
      interface_version = "tf_em_checkpoint_v2",
      em_signature = signature,
      draw_results = draw_results,
      settings = settings
    ), checkpoint_file)
  }
  message(sprintf(
    "Stage3 %s nuisance mode: %s; %d fit(s), %d completed, %d pending, %d workers.",
    algorithm,
    nuisance_mode,
    length(nuisance_draw_ids), sum(!vapply(draw_results, is.null, logical(1))),
    length(pending), worker_count
  ))

  posterior_mean_eta0 <- posterior_mean_phi <- NULL
  if (identical(nuisance_mode, "posterior_mean")) {
    posterior_mean_eta0 <- matrix(
      colMeans(matrix(
        stage2_input$eta0_draws,
        nrow = stage2_input$S,
        ncol = stage2_input$N * stage2_input$G
      )),
      nrow = stage2_input$N,
      ncol = stage2_input$G,
      dimnames = list(stage2_input$cell_names, stage2_input$target_genes)
    )
    posterior_mean_phi <- colMeans(stage2_input$phi_draws)
  }

  cluster <- NULL
  if (length(pending) && worker_count > 1L) {
    cluster <- parallel::makeCluster(worker_count)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    worker_common <- common
    worker_settings <- settings
    tfregact_library <- dirname(find.package("TFRegAct"))
    parallel::clusterExport(
      cluster,
      c("worker_common", "worker_settings", "tfregact_library"),
      envir = environment()
    )
    parallel::clusterEvalQ(cluster, {
      .libPaths(c(tfregact_library, .libPaths()))
      library(TFRegAct)
      .tf_em_worker_common <- worker_common
      .tf_em_worker_settings <- worker_settings
      NULL
    })
  }
  if (length(pending)) {
    batches <- split(pending, ceiling(seq_along(pending) / checkpoint_every))
    for (batch_ids in batches) {
      jobs <- lapply(batch_ids, function(draw_id) {
        eta0_job <- if (draw_id == 0L) {
          posterior_mean_eta0
        } else {
          matrix(
            stage2_input$eta0_draws[draw_id, , ],
            nrow = stage2_input$N,
            ncol = stage2_input$G,
            dimnames = list(stage2_input$cell_names, stage2_input$target_genes)
          )
        }
        phi_job <- if (draw_id == 0L) {
          posterior_mean_phi
        } else {
          as.numeric(stage2_input$phi_draws[draw_id, ])
        }
        list(
          draw_id = draw_id,
          eta0 = eta0_job,
          phi = phi_job
        )
      })
      batch_results <- if (worker_count > 1L) {
        parallel::parLapplyLB(cluster, jobs, tf_em_worker_entry)
      } else {
        lapply(jobs, tf_em_fit_one_draw, common = common, settings = settings)
      }
      names(batch_results) <- as.character(batch_ids)
      draw_results[names(batch_results)] <- batch_results
      save_checkpoint()
      message(sprintf(
        "EM checkpoint: %d/%d nuisance draws completed.",
        sum(!vapply(draw_results, is.null, logical(1))), length(draw_results)
      ))
    }
  }
  error_draws <- names(draw_results)[vapply(
    draw_results,
    function(x) !is.list(x) || !identical(x$status, "ok"),
    logical(1)
  )]
  if (length(error_draws)) {
    error_details <- vapply(error_draws, function(draw_id) {
      error_message <- draw_results[[draw_id]]$error
      if (is.null(error_message) || !nzchar(error_message)) "unknown error" else error_message
    }, character(1))
    tf_em_stop(
      paste0(
        "EM failed for nuisance draws: %s. Details: %s. ",
        "Rerun with `resume = TRUE` to retry them."
      ),
      paste(error_draws, collapse = ", "),
      paste(sprintf("%s=%s", error_draws, error_details), collapse = "; ")
    )
  }
  summarized <- tf_em_summarize_draw_results(draw_results, common, settings)
  result <- c(list(
    interface_version = if (identical(algorithm, "gibbs_ars")) {
      "tf_em_latent_activity_fit_v3"
    } else "tf_em_latent_activity_fit_v2",
    algorithm = algorithm,
    target_tf = stage2_input$target_tf,
    has_condition = isTRUE(stage2_input$has_condition),
    nuisance_mode = nuisance_mode,
    nuisance_draw_count = nuisance_draw_count,
    nuisance_draw_ids = nuisance_draw_ids,
    nuisance_source_draw_ids = stage2_input$nuisance_source_draw_ids,
    settings = settings,
    draw_results = draw_results
  ), summarized)
  class(result) <- c("TFEMLatentActivityFit", "list")
  attr(result, "em_signature") <- signature
  attr(result, "complete") <- TRUE
  attr(result, "output_file") <- output_file
  attr(result, "checkpoint_file") <- checkpoint_file
  saveRDS(result, output_file)
  result
}

# Explicit Stage2 name. Keep the shorter historical name as a compatible alias.
run_TF_stage2_EM_latent_activity <- run_TF_EM_latent_activity
