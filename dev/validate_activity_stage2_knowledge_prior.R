# Executable synthetic check of activity Stage 2 dispatch and Stage 3 handoff.
# Stage 1 selection is controlled; Stage 2 uses the real compact Stan model.
# The three EM iterations check execution only, not EM convergence or accuracy.
.libPaths(c("D:/Probscorer/.validation/stage2-knowledge-prior/library",
            "D:/R", .libPaths()))
library(TFRegAct)
ns <- asNamespace("TFRegAct")
outdir <- "D:/Probscorer/.validation/stage2-knowledge-prior"
fitdir <- file.path(outdir, "activity-smoke")
dir.create(fitdir, recursive = TRUE, showWarnings = FALSE)

set.seed(20261109)
N <- 64L
cells <- paste0("cell", seq_len(N))
condition <- rep(c("Normal", "AAA"), each = N / 2L)
shift <- as.numeric(condition == "AAA") - 0.5
expression <- rbind(T = runif(N, 0, 3), V = runif(N, 0, 3))
colnames(expression) <- cells
standardized <- t(scale(t(expression)))
Y <- rnbinom(N, size = 12, mu = exp(
  2 + (1.2 + 0.3 * shift) * standardized["T", ] -
    0.4 * standardized["V", ] + 0.2 * shift
))
analysis <- list(
  expr = expression, target = c("T", "Y"), confidence = c(8, 6),
  direction = c(1L, -1L), Y_exp = Y, libsize = rep(1000, N),
  sample = condition, batch = NULL
)
entry <- ns$prepare_TF_stage1_screening_stan_data(
  analysis, target_confidence = 8, target_direction = 1L,
  beta_prior_scale = 1.7, eta = 0.65, r_dir = 4,
  control_level = "Normal", disease_level = "AAA"
)
entry$predictor_edges <- data.frame(
  tf = entry$feature_names, effect = c("activation", "repression"),
  direction = entry$stan_data$direction, confidence = entry$stan_data$confidence
)
screening <- list(Y = entry)
attr(screening, "model_version") <- ns$tf_shared_stage2_model_version()
attr(screening, "prior_family") <- "laplace"
attr(screening, "direct_confounders_only") <- TRUE
attr(screening, "target_tf") <- "T"
attr(screening, "cell_count") <- N
prescreen <- list(Y = list(
  target_tf = data.frame(tf = "T", beta_mean = 1.2, beta_sd = 0.1),
  confounders = data.frame(tf = "V", beta_mean = -0.4, beta_sd = 0.1)
))
attr(prescreen, "pipeline_stage") <- "prescreen"
attr(prescreen, "confounders_filtered") <- TRUE
attr(prescreen, "target_interval") <- 0.9
attr(prescreen, "confounder_interval") <- 0.9
attr(prescreen, "filter_summary") <- data.frame(
  target_gene = "Y", confounders_removed_by_interval = 0
)
input <- ns$build_TF_mcmc_input_from_prescreen(prescreen, screening)
original_prior <- ns$tf_stage2_knowledge_prior(entry$stan_data)
stopifnot(identical(original_prior,
                    ns$tf_stage2_knowledge_prior(input$Y$stan_data)))

stage2 <- ns$run_TF_stage1_screening_regressions(
  screening_input = input, cores = 1L, inference = "mcmc",
  prior_family = "laplace", chains = 2L,
  iter_warmup = 300L, iter_sampling = 500L,
  adapt_delta = 0.95, max_treedepth = 12L, seed = 20261109L,
  target_interval = 0.9, confounder_interval = 0.9,
  stage2_draw_count = NULL, refresh = 0L,
  output_file = file.path(fitdir, "stage2-laplace-fit.rds"),
  checkpoint_file = file.path(fitdir, "stage2-laplace-checkpoint.rds"),
  resume = FALSE, force_refit = TRUE
)
stopifnot(stage2$Y$status == "ok",
          attr(stage2, "prior_family") == "laplace",
          isTRUE(stage2$Y$sampling$is_stage2),
          identical(stage2$Y$stage2_parameter_draws$knowledge_prior, original_prior))
filtered <- ns$filter_TF_mcmc_gene_results(
  stage2, output_file = file.path(fitdir, "stage2-laplace-filtered.rds")
)
stopifnot(identical(names(filtered), "Y"))
em_input <- ns$build_TF_EM_stage2_input(
  filtered, input, nuisance_storage = "posterior_mean",
  output_file = file.path(fitdir, "stage3-laplace-input.rds")
)
stopifnot(em_input$interface_version == "tf_em_stage2_input_v4",
          em_input$beta_prior_family == "laplace",
          em_input$beta_lambda_positive == original_prior$beta_lambda_positive[1],
          em_input$beta_lambda_negative == original_prior$beta_lambda_negative[1])
em <- ns$run_TF_EM_latent_activity(
  em_input, nuisance_draw_count = 0L, cores = 1L,
  quadrature_nodes = 101L, max_iter = 3L, min_iter = 2L,
  seed = 20261110L, resume = FALSE, save_traces = TRUE,
  output_file = file.path(fitdir, "stage3-laplace-em-fit.rds"),
  checkpoint_file = file.path(fitdir, "stage3-laplace-em-checkpoint.rds")
)
stopifnot(all(vapply(em$draw_results, function(x) x$status == "ok", logical(1))),
          nrow(em$activity_summary) == N,
          all(is.finite(em$activity_summary$activity_mean)),
          all(is.finite(em$activity_summary$activity_sd)),
          em$settings$beta_prior_family == "laplace")

report <- data.frame(
  check = c("activity_stage2_compact_stan", "stage3_em_execution_only"),
  status = "ok", prior_family = "laplace", cells = N,
  targets = 1L, predictors = input$Y$stan_data$P,
  gamma = original_prior$gamma, eta = original_prior$eta,
  r_dir = original_prior$r_dir,
  focal_lambda_positive = original_prior$beta_lambda_positive[1],
  focal_lambda_negative = original_prior$beta_lambda_negative[1],
  chains = c(2L, NA_integer_), retained_draws = c(1000L, NA_integer_),
  max_beta_rhat = c(stage2$Y$diagnostics$max_beta_rhat, NA_real_),
  min_beta_ess_bulk = c(stage2$Y$diagnostics$min_beta_ess_bulk, NA_real_),
  em_iterations = c(NA_integer_, em$convergence_summary$iterations),
  em_converged = c(NA, em$convergence_summary$converged),
  stringsAsFactors = FALSE
)
write.csv(report, file.path(outdir, "activity-stan-em-smoke-summary.csv"),
          row.names = FALSE)
saveRDS(list(report = report, original_prior = original_prior,
             stage2_model = attr(stage2, "stan_file"),
             stage3_interface = em_input$interface_version),
        file.path(outdir, "activity-stan-em-smoke-summary.rds"))
print(report)
cat(paste0("PASS: real activity Stage 2 compact Stan dispatch, filtering, v4 ",
           "handoff and short EM execution retain the knowledge Laplace prior.\n"))
