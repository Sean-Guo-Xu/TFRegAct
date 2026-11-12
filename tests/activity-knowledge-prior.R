library(TFRegAct)
ns <- asNamespace("TFRegAct")

expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"),
            grepl(pattern, conditionMessage(error)))
}

# Use a nondefault prior and put the focal TF between two adjustment TFs.
# Retaining T and C2 must preserve their original order and knowledge priors.
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
data <- mcmc_input$Y$stan_data
stopifnot(
  identical(mcmc_input$Y$feature_names, c("T", "C2")),
  identical(as.numeric(data$confidence), c(9, 4)),
  identical(as.integer(data$direction), c(1L, 0L)),
  data$gamma == 2, data$eta == 0.7, data$r_dir == 4,
  identical(data$beta_init, c(8, -9)),
  is.null(data$beta_prior_mean), is.null(data$beta_prior_sd)
)
scale <- 2 * (c(9, 4) / 10) ^ 0.7
knowledge_prior <- mcmc_input$Y$stage2_knowledge_prior
stopifnot(
  isTRUE(all.equal(as.numeric(knowledge_prior$beta_lambda_positive), 1 / scale)),
  isTRUE(all.equal(as.numeric(knowledge_prior$beta_lambda_negative), c(4, 1) / scale))
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
  fake_fit, mcmc_input$Y, "mcmc", NULL, "laplace", 0.9, 0.9
)
bundle <- summarized$stage2_parameter_draws
stopifnot(
  bundle$interface_version == "tf_em_stage2_parameter_draws_v2",
  bundle$beta_prior_family == "laplace",
  is.null(bundle$beta_prior_mean), is.null(bundle$beta_prior_sd)
)
results <- list(Y = c(list(
  status = "ok", target_tf_name = "T", condition_model = FALSE
), summarized))
attr(results, "inference") <- "mcmc"
attr(results, "prior_family") <- "laplace"
attr(results, "input_pipeline_stage") <- "mcmc_input"
attr(results, "model_version") <- ns$tf_shared_stage2_model_version()
filtered_file <- tempfile(fileext = ".rds")
filtered <- ns$filter_TF_mcmc_gene_results(results, output_file = filtered_file)
unlink(filtered_file)
em_input <- ns$build_TF_EM_stage2_input(
  filtered, mcmc_input, nuisance_storage = "posterior_mean"
)
stopifnot(
  em_input$interface_version == "tf_em_stage2_input_v4",
  em_input$beta_prior_family == "laplace",
  is.null(em_input$beta_prior_mean), is.null(em_input$beta_prior_sd),
  length(em_input$beta_lambda_positive) == 1L,
  em_input$beta_lambda_positive == knowledge_prior$beta_lambda_positive[1],
  em_input$beta_lambda_negative == knowledge_prior$beta_lambda_negative[1]
)
bad <- filtered
bad$Y$stage2_parameter_draws$beta_lambda_positive[1] <- 1
expect_error(ns$build_TF_EM_stage2_input(
  bad, mcmc_input, nuisance_storage = "posterior_mean"
), "Knowledge Laplace prior metadata differ")

# Both stages now have Laplace priors. Their stage marker, rather than the
# family name, must select the appropriate fitter. Mock in a private scope;
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
  is_stage2 = TRUE, inference = "mcmc", prior_family = "laplace",
  seed = 123L, stage2_draw_count = NULL, chains = 1L, iter_sampling = 100L,
  target_interval = 0.9, confounder_interval = 0.9
)
job <- list(entry = mcmc_input$Y, model_index = 1L, retry_seed_offset = 0L)
fit <- worker(job, sampling)
stopifnot(fit$status == "ok", scope$last_route == "stage2")
sampling$is_stage2 <- FALSE
fit <- worker(job, sampling)
stopifnot(fit$status == "ok", scope$last_route == "stage1")

cat(paste0("PASS: activity knowledge-prior subsets, posterior-only initialization, ",
           "v4 handoff, incompatible metadata rejection and stage dispatch.\n"))
