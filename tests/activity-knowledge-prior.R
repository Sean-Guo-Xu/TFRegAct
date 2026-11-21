library(TFRegAct)
ns <- asNamespace("TFRegAct")

expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"),
            grepl(pattern, conditionMessage(error)))
}

# Use nondefault Stage 1 settings and put the focal TF between adjustment TFs.
# Stage 2 must retain predictor order and use only confidence/direction for its Normal prior.
analysis <- list(
  expr = matrix(c(1, 3, 2, 5, 4, 6,
                  2, 4, 1, 3, 6, 5,
                  6, 2, 4, 1, 5, 3), 3, byrow = TRUE,
                dimnames = list(c("C1", "T", "C2"), paste0("c", 1:6))),
  target = c("T", "Y"), confidence = c(2, 8, 4),
  direction = c(-1L, 1L, 0L), Y_exp = 1:6, libsize = 101:106,
  sample = NULL, batch = NULL
)
entry <- ns$prepare_TF_stage1_screening_stan_data(
  analysis, target_confidence = 9, target_direction = 1L,
  beta_prior_scale = 2, eta = 0.7, r_dir = 4
)
entry$predictor_edges <- data.frame(
  tf = entry$feature_names,
  effect = c("repression", "activation", "unknown"),
  direction = entry$stan_data$direction,
  confidence = entry$stan_data$confidence
)
screening_input <- list(Y = entry)
attr(screening_input, "model_version") <- ns$tf_shared_stage2_model_version()
attr(screening_input, "prior_family") <- "laplace"
attr(screening_input, "direct_confounders_only") <- TRUE
attr(screening_input, "target_tf") <- "T"
attr(screening_input, "cell_count") <- 6L

prescreen <- list(Y = list(
  target_tf = data.frame(tf = "T", beta_mean = 8, beta_sd = 0.001),
  confounders = data.frame(tf = "C2", beta_mean = -9, beta_sd = 0.01)
))
attr(prescreen, "pipeline_stage") <- "prescreen"
attr(prescreen, "confounders_filtered") <- TRUE
attr(prescreen, "filter_summary") <- data.frame(
  target_gene = "Y", confounders_removed_by_interval = 1
)
mcmc_input <- ns$build_TF_mcmc_input_from_prescreen(prescreen, screening_input)

# Real Stage 1 inputs omit the Normal-only fields. Both stage signatures must
# be constructible, while a changed Stage 2 Normal prior invalidates its cache.
laplace_input <- screening_input
laplace_input$Y$stan_data$beta_prior_mean <- NULL
laplace_input$Y$stan_data$beta_prior_sd <- NULL
signature_stage1 <- ns$tf_stage1_screening_input_signature(laplace_input, "Y")
signature_stage2 <- ns$tf_stage1_screening_input_signature(mcmc_input, "Y")
changed_prior_input <- mcmc_input
changed_prior_input$Y$stan_data$beta_prior_mean[1] <-
  changed_prior_input$Y$stan_data$beta_prior_mean[1] + 0.1
stopifnot(is.character(signature_stage1), length(signature_stage1) == 1L,
          nzchar(signature_stage1), !identical(signature_stage1, signature_stage2),
          !identical(signature_stage2,
                     ns$tf_stage1_screening_input_signature(changed_prior_input, "Y")))
data <- mcmc_input$Y$stan_data
stopifnot(
  identical(mcmc_input$Y$feature_names, c("T", "C2")),
  identical(as.numeric(data$confidence), c(9, 4)),
  identical(as.integer(data$direction), c(1L, 0L)),
  data$gamma == 2, data$eta == 0.7, data$r_dir == 4,
  identical(data$beta_init, c(8, -9)),
  identical(attr(mcmc_input, "prior_family"), "normal")
)
weight <- (c(9, 4) - 1) / 9
knowledge_prior <- mcmc_input$Y$stage2_knowledge_prior
stopifnot(
  identical(knowledge_prior$prior_family, "normal"),
  isTRUE(all.equal(as.numeric(data$beta_prior_mean), c(1, 0) * (0.2 + 0.4 * weight))),
  isTRUE(all.equal(as.numeric(data$beta_prior_sd), sqrt(0.5 + 1.5 * weight))),
  isTRUE(all.equal(as.numeric(knowledge_prior$beta_prior_variance), 0.5 + 1.5 * weight)),
  is.null(knowledge_prior$beta_lambda_positive)
)

# Changing fitted posterior moments changes initialization, never the prior.
different_posterior <- prescreen
different_posterior$Y$target_tf$beta_mean <- 200
different_posterior$Y$target_tf$beta_sd <- 1000
different_input <- ns$build_TF_mcmc_input_from_prescreen(
  different_posterior, screening_input
)
stopifnot(
  identical(knowledge_prior, different_input$Y$stage2_knowledge_prior),
  !identical(data$beta_init, different_input$Y$stan_data$beta_init)
)

# Exercise posterior summaries and the complete Stage 2 -> Stage 3 handoff
# without compiling or sampling a Stan model.
draws <- cbind(
  alpha = rep(0, 100), condition_effect = rep(0, 100), phi = rep(5, 100),
  b1 = seq(0.7, 1.3, length.out = 100),
  b2 = seq(-0.3, 0.3, length.out = 100),
  beta_target_control = seq(0.7, 1.3, length.out = 100),
  beta_target_disease = seq(0.7, 1.3, length.out = 100),
  beta_target_delta = rep(0, 100)
)
colnames(draws)[4:5] <- c("beta[1]", "beta[2]")
fake_fit <- list(
  draws = function(...) draws,
  summary = function(...) data.frame(
    variable = c("beta[1]", "beta[2]"),
    rhat = 1, ess_bulk = 100, ess_tail = 100
  ),
  output_files = function(...) character(0),
  return_codes = function(...) 0L,
  num_chains = function(...) 1L
)
summarized <- ns$tf_stage1_screening_summarize_fit(
  fake_fit, mcmc_input$Y, "mcmc", NULL, "normal", 0.9, 0.9
)
bundle <- summarized$stage2_parameter_draws
stopifnot(
  bundle$interface_version == "tf_em_stage2_parameter_draws_v3",
  bundle$beta_prior_family == "normal",
  identical(names(bundle$beta_prior_mean), c("T", "C2")),
  isTRUE(all.equal(unname(bundle$beta_prior_mean), data$beta_prior_mean)),
  isTRUE(all.equal(unname(bundle$beta_prior_sd), data$beta_prior_sd)),
  is.null(bundle$beta_lambda_positive)
)
results <- list(Y = c(list(
  status = "ok", target_tf_name = "T", condition_model = FALSE
), summarized))
attr(results, "inference") <- "mcmc"
attr(results, "prior_family") <- "normal"
attr(results, "input_pipeline_stage") <- "mcmc_input"
attr(results, "model_version") <- ns$tf_shared_stage2_model_version()
filtered_file <- tempfile(fileext = ".rds")
filtered <- ns$filter_TF_mcmc_gene_results(results, output_file = filtered_file)
unlink(filtered_file)
em_input <- ns$build_TF_EM_stage2_input(
  filtered, mcmc_input, nuisance_storage = "posterior_mean"
)
stopifnot(
  em_input$interface_version == "tf_em_stage2_input_v6",
  em_input$beta_prior_family == "normal",
  em_input$beta_prior_source == "stage2_posterior_normal_approximation",
  length(em_input$beta_prior_mean) == 1L,
  isTRUE(all.equal(unname(em_input$beta_prior_mean), mean(draws[, "beta[1]"]))),
  isTRUE(all.equal(unname(em_input$beta_prior_variance), stats::sd(draws[, "beta[1]"]) ^ 2)),
  isTRUE(all.equal(unname(em_input$beta_prior_sd), stats::sd(draws[, "beta[1]"]))),
  em_input$beta_delta_prior_mean == 0, em_input$beta_delta_prior_variance == 0,
  em_input$beta_delta_prior_sd == 0, em_input$beta_delta_prior_covariance == 0,
  em_input$stage2_posterior_draw_count == 100,
  !isTRUE(all.equal(unname(em_input$beta_prior_variance), knowledge_prior$beta_prior_variance[1])),
  is.null(em_input$beta_lambda_positive)
)

# The full posterior supplies the same prior for every nuisance-storage choice.
selected <- ns$build_TF_EM_stage2_input(filtered, mcmc_input,
  nuisance_storage = "all_draws", nuisance_draw_ids = c(1L, 100L))
for (field in c("beta_prior_mean", "beta_prior_variance", "beta_prior_sd",
                "stage2_posterior_prior_signature")) stopifnot(identical(em_input[[field]], selected[[field]]))
thin_summary <- ns$tf_stage1_screening_summarize_fit(
  fake_fit, mcmc_input$Y, "mcmc", 2L, "normal", 0.9, 0.9)
thin_results <- results
thin_results$Y <- c(list(status = "ok", target_tf_name = "T", condition_model = FALSE), thin_summary)
thin_file <- tempfile(fileext = ".rds")
thin_filtered <- ns$filter_TF_mcmc_gene_results(thin_results, output_file = thin_file)
unlink(thin_file)
thin_input <- ns$build_TF_EM_stage2_input(thin_filtered, mcmc_input)
stopifnot(thin_input$S == 2L, thin_input$stage2_posterior_draw_count == 100,
          identical(thin_input$beta_prior_variance, em_input$beta_prior_variance))

# Previously filtered v3 objects need no refit: their summaries used all draws.
legacy <- filtered
legacy$Y$target_screening_draws <- NULL
legacy$Y$diagnostics <- NULL
legacy$Y$stage2_parameter_draws$target_posterior_moments <- NULL
legacy$Y$target_tf$beta_mean_delta_covariance <- NULL
legacy_input <- ns$build_TF_EM_stage2_input(legacy, mcmc_input)
stopifnot(isTRUE(all.equal(legacy_input$beta_prior_variance, em_input$beta_prior_variance)),
          is.na(legacy_input$stage2_posterior_draw_count),
          legacy_input$stage2_posterior_prior_metadata$Y$moment_source == "full_stage2_posterior_summary")
bad <- filtered
bad$Y$stage2_parameter_draws$beta_prior_mean[1] <- 1
expect_error(ns$build_TF_EM_stage2_input(
  bad, mcmc_input, nuisance_storage = "posterior_mean"
), "Knowledge Normal prior metadata differ")
old <- filtered
old$Y$stage2_parameter_draws$interface_version <- "tf_em_stage2_parameter_draws_v2"
expect_error(ns$build_TF_EM_stage2_input(old, mcmc_input), "incompatible Stage 2 draw interface")

# A conditioned posterior retains beta/delta covariance, including when compact
# nuisance draws are thinned. The legacy variance identity yields the same joint prior.
condition_input <- mcmc_input
condition_input$Y$condition_model <- TRUE
condition_input$Y$stan_data$has_condition <- 1L
condition_input$Y$stan_data$condition <- rep(0:1, each = 3)
condition_draws <- draws
beta <- draws[, "beta[1]"]
delta <- 0.3 + 0.2 * (beta - mean(beta)) + 0.05 * sin(seq_along(beta))
condition_draws[, "beta_target_delta"] <- delta
condition_draws[, "beta_target_control"] <- beta - delta / 2
condition_draws[, "beta_target_disease"] <- beta + delta / 2
condition_fit <- fake_fit
condition_fit$draws <- function(...) condition_draws
condition_summary <- ns$tf_stage1_screening_summarize_fit(
  condition_fit, condition_input$Y, "mcmc", 2L, "normal", 0.9, 0.9)
condition_model <- c(list(condition_model = TRUE), condition_summary)
joint <- ns$tf_stage2_posterior_normal_prior(condition_model, "Y", TRUE)
stopifnot(isTRUE(all.equal(joint$beta_prior_mean, mean(beta))),
          isTRUE(all.equal(joint$beta_delta_prior_mean, mean(delta))),
          isTRUE(all.equal(joint$beta_delta_prior_variance, stats::sd(delta) ^ 2)),
          isTRUE(all.equal(joint$beta_delta_prior_covariance, stats::cov(beta, delta))),
          joint$source_draw_count == 100)
condition_models <- filtered
condition_models$Y <- condition_model
condition_handoff <- ns$build_TF_EM_stage2_input(condition_models, condition_input,
                                               nuisance_storage = "posterior_mean")
stopifnot(isTRUE(all.equal(unname(condition_handoff$beta_delta_prior_covariance), stats::cov(beta, delta))),
          isTRUE(all.equal(unname(condition_handoff$beta_delta_prior_sd), stats::sd(delta))))
condition_model$target_screening_draws <- NULL
condition_model$diagnostics <- NULL
condition_model$stage2_parameter_draws$target_posterior_moments <- NULL
condition_model$target_tf$beta_mean_delta_covariance <- NULL
legacy_joint <- ns$tf_stage2_posterior_normal_prior(condition_model, "Y", TRUE)
for (field in c("beta_prior_mean", "beta_prior_variance", "beta_delta_prior_mean",
                "beta_delta_prior_variance", "beta_delta_prior_covariance")) {
  stopifnot(isTRUE(all.equal(joint[[field]], legacy_joint[[field]], tolerance = 1e-10)))
}
bad_posterior <- condition_model
bad_posterior$target_tf$beta_delta_sd <- 0
expect_error(ns$tf_stage2_posterior_normal_prior(bad_posterior, "Y", TRUE), "not positive definite")
bad_posterior <- filtered$Y
bad_posterior$target_tf$beta_sd <- 0
expect_error(ns$tf_stage2_posterior_normal_prior(bad_posterior, "Y", FALSE), "disagree|not positive")

# The stage marker selects Normal refinement versus Laplace screening.
# Mock in a private scope;
# the installed namespace is never changed.
scope <- new.env(parent = ns)
scope$.tf_stage1_worker_model <- NULL
scope$tf_fit_shared_stage2 <- function(...) {
  scope$last_route <- "stage2"
  fake_fit
}
scope$tf_fit_shared_laplace_stage1 <- function(...) {
  scope$last_route <- "stage1"
  fake_fit
}
worker <- ns$tf_stage1_screening_fit_one
environment(worker) <- scope
sampling <- list(
  is_stage2 = TRUE, inference = "mcmc", prior_family = "normal",
  seed = 123L, stage2_draw_count = NULL, chains = 1L, iter_sampling = 100L,
  target_interval = 0.9, confounder_interval = 0.9
)
job <- list(entry = mcmc_input$Y, model_index = 1L, retry_seed_offset = 0L)
fit <- worker(job, sampling)
stopifnot(fit$status == "ok", scope$last_route == "stage2")
sampling$is_stage2 <- FALSE
sampling$prior_family <- "laplace"
fit <- worker(job, sampling)
stopifnot(fit$status == "ok", scope$last_route == "stage1",
          fit$stage2_parameter_draws$beta_prior_family == "laplace",
          is.null(fit$stage2_parameter_draws$beta_prior_mean))

cat(paste0("PASS: Stage 2 network Normal prior, full-posterior v6 Normal transfer with covariance, ",
           "storage/thinning independence, reusable legacy summaries and unchanged Stage 1 Laplace dispatch.\n"))
