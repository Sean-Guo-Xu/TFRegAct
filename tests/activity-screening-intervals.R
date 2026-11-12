# Load the source backend privately so this test needs no package installation.
r_dir <- c("R", "../R", "TFRegAct/R")
r_dir <- r_dir[vapply(r_dir, function(path) {
  all(file.exists(file.path(path, c("shared_stage2.R", "activity_stage1.R"))))
}, logical(1))][1L]
stopifnot(length(r_dir) == 1L, !is.na(r_dir))
ns <- new.env(parent = globalenv())
for (file in c("shared_stage2.R", "activity_stage1.R")) {
  sys.source(file.path(r_dir, file), envir = ns)
}

expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

target_values <- seq(-0.09, 0.91, length.out = 100L)
stopifnot(
  ns$tf_stage1_screening_quantile_interval(target_values, 0.9)[1L] < 0,
  ns$tf_stage1_screening_quantile_interval(target_values, 0.8)[1L] > 0
)
for (level in c(0.8, 0.81, 0.9, 0.95)) {
  expected <- stats::quantile(target_values,
    c((1 - level) / 2, (1 + level) / 2), names = FALSE, type = 7L)
  stopifnot(identical(ns$tf_stage1_screening_quantile_interval(target_values, level), expected))
}
for (level in list(0, 1, NA_real_, Inf, c(0.8, 0.9), "0.8")) {
  expect_error(ns$tf_stage1_screening_quantile_interval(target_values, level),
               "strictly between 0 and 1")
}
expect_error(ns$tf_stage1_screening_quantile_interval(c(1, NA_real_), 0.8),
             "finite vector")

# Exercise both the one-predictor case and a focal TF after a confounder.
# The EM interface deliberately keeps one draw, while rescreening must retain
# all 100 source draws and exact condition-specific target marginals.
for (P in c(1L, 2L)) {
  for (has_condition in c(FALSE, TRUE)) {
    features <- if (P == 1L) "T" else c("C", "T")
    beta <- if (P == 1L) matrix(target_values, ncol = 1L) else {
      cbind(seq(-0.2, 0.4, length.out = 100L), target_values)
    }
    colnames(beta) <- sprintf("beta[%d]", seq_len(P))
    control <- target_values - if (has_condition) 0.3 else 0
    disease <- target_values + if (has_condition) 0.3 else 0
    draws <- cbind(alpha = rep(0, 100L), condition_effect = rep(0, 100L),
      phi = rep(5, 100L), beta, beta_target_control = control,
      beta_target_disease = disease, beta_target_delta = disease - control)
    entry <- list(target_gene = "Y", target_tf = "T", target_tf_index = P,
      feature_names = features, use_batch = FALSE, batch_levels = character(),
      condition_model = has_condition,
      predictor_edges = data.frame(tf = features, effect = "activation",
        direction = 1L, confidence = 5),
      stan_data = list(N = 4L, P = P, K_batch = 0L,
        has_condition = as.integer(has_condition),
        condition = if (has_condition) c(0L, 0L, 1L, 1L) else integer(4L),
        confidence = rep(5, P), direction = rep(1L, P), gamma = 1,
        eta = 0.5, r_dir = 3, target_interaction_sd = 0.5))
    fake_fit <- list(draws = function(...) draws,
      output_files = function(...) character())
    summarized <- ns$tf_stage1_screening_summarize_fit(fake_fit, entry,
      inference = "variational", stage2_draw_count = 1L, prior_family = "laplace",
      target_interval = 0.8, confounder_interval = 0.81)
    # Match the actual CmdStanR class and its non-dropping column extraction.
    if (requireNamespace("posterior", quietly = TRUE)) {
      posterior_draws <- posterior::as_draws_matrix(draws)
      real_format_fit <- list(draws = function(...) posterior_draws,
        output_files = function(...) character())
      real_format_summary <- ns$tf_stage1_screening_summarize_fit(real_format_fit, entry,
        inference = "variational", stage2_draw_count = 1L, prior_family = "laplace",
        target_interval = 0.8, confounder_interval = 0.81)
      stopifnot(isTRUE(all.equal(summarized, real_format_summary,
                                 check.attributes = FALSE, tolerance = 1e-12)))
    }
    payload <- summarized$target_screening_draws
    stopifnot(payload$interface_version == "tf_target_screening_draws_v1",
      payload$draw_count == 100L, identical(payload$overall, target_values),
      summarized$stage2_parameter_draws$draw_count == 1L,
      nrow(summarized$stage2_parameter_draws$beta) == 1L,
      summarized$target_tf$screening_interval_level == 0.8,
      summarized$target_tf$screening_interval_lower > 0,
      isTRUE(all.equal(summarized$target_tf$beta_q025,
        unname(stats::quantile(target_values, 0.025)))),
      isTRUE(all.equal(summarized$target_tf$beta_q05,
        unname(stats::quantile(target_values, 0.05)))),
      isTRUE(all.equal(summarized$target_tf$beta_q95,
        unname(stats::quantile(target_values, 0.95)))),
      isTRUE(all.equal(summarized$target_tf$beta_q975,
        unname(stats::quantile(target_values, 0.975)))))
    if (has_condition) {
      stopifnot(identical(payload$control, control), identical(payload$disease, disease),
        isTRUE(all.equal(summarized$target_tf$beta_control_interval_lower,
          ns$tf_stage1_screening_quantile_interval(control, 0.8)[1L])),
        isTRUE(all.equal(summarized$target_tf$beta_disease_interval_upper,
          ns$tf_stage1_screening_quantile_interval(disease, 0.8)[2L])))
    } else {
      stopifnot(is.null(payload$control), is.null(payload$disease))
    }
    if (P > 1L) {
      stopifnot(summarized$confounders$screening_interval_level == 0.81)
    }

    scope <- new.env(parent = ns)
    scope$.tf_stage1_worker_model <- NULL
    scope$tf_fit_shared_laplace_stage1 <- function(...) fake_fit
    worker <- ns$tf_stage1_screening_fit_one
    environment(worker) <- scope
    sampling <- list(is_stage2 = FALSE, inference = "variational", prior_family = "laplace",
      seed = 123L, stage2_draw_count = 1L, target_interval = 0.8, confounder_interval = 0.81)
    fit <- worker(list(entry = entry, model_index = 1L, retry_seed_offset = 0L), sampling)
    stopifnot(fit$status == "ok", identical(fit$target_screening_draws, payload))
  }
}

# Stop after configuration validation, before any Stan resolution or file writes.
scope <- new.env(parent = ns)
scope$tf_stage1_screening_require_pkg <- function(...) invisible(NULL)
scope$tf_stage1_screening_validate_input <- function(x) names(x)
scope$tf_stage1_screening_input_signature <- function(...) "mock_signature"
scope$tf_stage1_screening_resolve_stan_file <- function(...) {
  scope$sampling <- parent.frame()$sampling
  stop("configuration captured")
}
runner <- ns$run_TF_stage1_screening_regressions
environment(runner) <- scope
expect_error(runner(list(Y = list()), inference = "variational",
  target_interval = 0.8, confounder_interval = 0.81), "configuration captured")
stopifnot(scope$sampling$target_interval == 0.8, scope$sampling$confounder_interval == 0.81)
expect_error(runner(list(Y = list()), target_interval = 0), "strictly between 0 and 1")
expect_error(runner(list(Y = list()), confounder_interval = 1), "strictly between 0 and 1")

cat("PASS: arbitrary screening intervals and full target draws despite the one-draw EM budget; condition and worker propagation. No fitting.\n")
