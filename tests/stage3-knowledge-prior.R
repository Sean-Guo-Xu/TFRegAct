# Run directly with Rscript; no installed package or optional dependencies needed.
source(if (file.exists("R/activity_stage3_em.R")) {
  "R/activity_stage3_em.R"
} else {
  "../R/activity_stage3_em.R"
})

common <- list(N = 4L, G = 1L, Y = matrix(c(0, 2, 4, 1), 4, 1),
  has_condition = TRUE, condition_shift = c(-0.5, -0.5, 0.5, 0.5),
  target_center = 0.3, target_scale = 1.1,
  beta_lambda_positive = 2, beta_lambda_negative = 5,
  target_interaction_sd = 0.8, target_genes = "Y")
estep <- list(nodes = matrix(c(0, 0.4, 1.2, 2, 0.7, 1, 1.8, 3), 4, 2),
  weights = matrix(rep(c(0.3, 0.7), each = 4), 4, 2))
eta0 <- c(-0.1, 0.1, 0.2, -0.2)
objective <- function(p, branch = NULL, information = FALSE) {
  tf_em_gene_expected_log_posterior(p, 1L, common, eta0, 5, estep,
    return_information = information, beta_branch = branch)
}
for (beta in c(-0.6, 0.4)) {
  p <- c(beta, 0.2)
  analytic <- objective(p)$gradient
  numerical <- vapply(1:2, function(j) {
    plus <- minus <- p
    plus[j] <- plus[j] + 1e-6
    minus[j] <- minus[j] - 1e-6
    (objective(plus)$value - objective(minus)$value) / 2e-6
  }, numeric(1))
  stopifnot(max(abs(analytic - numerical)) < 1e-6)
  # The only beta-prior dependence in the objective is the AL linear penalty.
  unpenalized <- common
  unpenalized$beta_lambda_positive <- unpenalized$beta_lambda_negative <- 0
  likelihood_delta_prior <- tf_em_gene_expected_log_posterior(
    p, 1L, unpenalized, eta0, 5, estep)$value
  expected <- if (beta >= 0) -2 * beta else 5 * beta
  stopifnot(abs(objective(p)$value - likelihood_delta_prior - expected) < 1e-10,
    identical(objective(p, information = TRUE)$information,
      tf_em_gene_expected_log_posterior(p, 1L, unpenalized, eta0, 5, estep,
        return_information = TRUE)$information))
}
right <- objective(c(0, 0.2), "positive")$gradient[[1]]
left <- objective(c(0, 0.2), "negative")$gradient[[1]]
stopifnot(abs(left - right - 7) < 1e-10,
  tf_em_beta_log_prior(-0.6, 0.2, common) == -3 - 0.5 * (0.2 / 0.8)^2)

# Known one-predictor NB optima exercise positive, negative and cusp branches.
small <- list(N = 2L, G = 1L, Y = matrix(1, 2, 1), has_condition = FALSE,
  condition_shift = c(0, 0), target_center = 0, target_scale = 1,
  beta_lambda_positive = 1, beta_lambda_negative = 1,
  target_interaction_sd = 1, target_genes = "Y")
fixed_estep <- list(nodes = matrix(1, 2, 1), weights = matrix(1, 2, 1))
settings <- list(mstep_maxit = 300L, mstep_reltol = 1e-11)
fit_small <- function(y, condition = FALSE) {
  x <- small
  x$Y[,] <- y
  x$has_condition <- condition
  x$condition_shift <- if (condition) c(-0.5, 0.5) else c(0, 0)
  tf_em_mstep(x, matrix(0, 2, 1), 5, fixed_estep, -0.2, 0, settings)
}
zero <- fit_small(1)
positive <- fit_small(10)
negative <- fit_small(0)
stopifnot(zero$beta_mean == 0, zero$at_prior_kink,
  all(is.na(zero$covariance)),
  abs(positive$beta_mean - log(5 * 19 / 11)) < 1e-5,
  abs(negative$beta_mean - log(5 / 9)) < 1e-5,
  !positive$at_prior_kink, !negative$at_prior_kink,
  fit_small(1, TRUE)$at_prior_kink)
scaled <- tf_em_scaled_beta_change(1, 0, 0, 0, small)
stopifnot(abs(scaled - 1 / sqrt(2)) < 1e-12)

# v4 handoff contains knowledge rates, never a Stage1/Stage2 posterior Normal.
input <- list(interface_version = "tf_em_stage2_input_v4", N = 2L, G = 1L,
  S = 2L, Y = matrix(c(1, 2), 2, 1), eta0_draws = array(0, c(2, 2, 1)),
  phi_draws = matrix(5, 2, 1), has_condition = FALSE, condition = c(0L, 0L),
  target_tf_expression = c(0, 1), target_tf_center = 0.5, target_tf_scale = 1,
  beta_target_mean_init = 0, beta_target_delta_init = 0,
  beta_prior_family = "laplace", beta_lambda_positive = 100,
  beta_lambda_negative = 100, target_interaction_sd = 1,
  cell_names = c("c1", "c2"), target_genes = "Y", target_tf = "T")
stopifnot(identical(tf_em_validate_input(input), input))
invalid <- input
invalid$interface_version <- "tf_em_stage2_input_v3"
stopifnot(inherits(try(tf_em_validate_input(invalid), silent = TRUE), "try-error"))
invalid <- input
invalid$beta_prior_family <- "normal"
stopifnot(inherits(try(tf_em_validate_input(invalid), silent = TRUE), "try-error"))
invalid <- input
invalid$beta_lambda_positive <- 0
stopifnot(inherits(try(tf_em_validate_input(invalid), silent = TRUE), "try-error"))

output <- tempfile(fileext = ".rds")
checkpoint <- tempfile(fileext = ".rds")
run <- function(value, draws, resume = FALSE) run_TF_EM_latent_activity(
  value, algorithm = "quadrature_em", nuisance_draw_count = draws, cores = 1L, max_iter = 2L, min_iter = 2L,
  quadrature_nodes = 5L, output_file = output, checkpoint_file = checkpoint,
  resume = resume)
single <- run(input, 0L)
multiple <- run(input, 2L)
stopifnot(identical(single$interface_version, "tf_em_latent_activity_fit_v2"),
  single$settings$beta_prior_family == "laplace",
  all(multiple$beta_at_prior_kink_by_draw),
  all(multiple$beta_mean_map_by_draw == 0),
  all(is.na(multiple$beta_mean_local_normal_sample_by_draw)),
  all(is.na(multiple$beta_summary$beta_mean_sd)),
  all(multiple$beta_summary$beta_mean_map == 0),
  all(is.finite(multiple$activity_summary$activity_mean)))
# A mixture with one unavailable cusp approximation cannot silently summarize
# the remaining nuisance component as though it represented the whole mixture.
mixed_draws <- multiple$draw_results
mixed_draws[[2]]$beta_at_prior_kink <- FALSE
mixed_draws[[2]]$beta_mean <- 0.4
mixed_draws[[2]]$beta_posterior_sample_mean <- 0.5
mixed_draws[[2]]$beta_posterior_sample_delta <- 0
summary_common <- c(input, list(activity_anchor = input$target_tf_expression,
  active_prior_probability = c(0.2, 0.9)))
mixed <- tf_em_summarize_draw_results(mixed_draws, summary_common, list())
approximation_columns <- grep(
  "^beta_(mean|delta|control|disease)_(mean|sd|q05|median|q95)$",
  names(mixed$beta_summary), value = TRUE)
stopifnot(all(is.na(mixed$beta_summary[, approximation_columns])),
  mixed$beta_summary$beta_at_prior_kink,
  mixed$beta_summary$beta_mean_map == 0.2,
  identical(as.numeric(mixed$beta_mean_map_by_draw), c(0, 0.4)))
changed <- input
changed$beta_lambda_positive <- 101
new_fit <- run(changed, 2L, resume = TRUE)
stopifnot(!identical(attr(new_fit, "em_signature"), attr(multiple, "em_signature")))
unlink(c(output, checkpoint))
message("Stage3 asymmetric-Laplace objective, branch MAP, kink and v4 EM checks passed.")
