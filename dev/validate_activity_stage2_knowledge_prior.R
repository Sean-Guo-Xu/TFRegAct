# Executable synthetic check of activity Stage 2 dispatch and Stage 3 handoff.
# Stage 1 selection is controlled; Stage 2 uses the real compact Stan model.
# The three EM iterations check execution only, not EM convergence or accuracy.
# --reuse-stage2 reuses the already validated v3 fit; no MCMC refit is needed.
.libPaths(c(Sys.getenv("TFREGACT_VALIDATION_LIB", "D:/Probscorer/.validation/stage3-posterior-normal/library"),
            "D:/R", .libPaths()))
library(TFRegAct)
ns <- asNamespace("TFRegAct")
outdir <- Sys.getenv("TFREGACT_VALIDATION_OUTPUT", "D:/Probscorer/.validation/stage3-posterior-normal")
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

stage2 <- if ("--reuse-stage2" %in% commandArgs(TRUE)) {
  cache <- Sys.getenv("TFREGACT_STAGE2_REUSE_FILE",
    "D:/Probscorer/.validation/normal-knowledge-prior/activity-smoke/stage2-normal-fit.rds")
  stopifnot(file.exists(cache))
  value <- readRDS(cache)
  stopifnot(identical(attr(value, "model_version"), ns$tf_shared_stage2_model_version()),
            identical(attr(value, "prior_family"), "normal"))
  message("Reusing existing Stage 2 posterior: ", normalizePath(cache, winslash = "/"))
  value
} else ns$run_TF_stage1_screening_regressions(
  screening_input = input, cores = 1L, inference = "mcmc",
  prior_family = "normal", chains = 2L,
  iter_warmup = 300L, iter_sampling = 500L,
  adapt_delta = 0.95, max_treedepth = 12L, seed = 20261109L,
  target_interval = 0.9, confounder_interval = 0.9,
  stage2_draw_count = NULL, refresh = 0L,
  output_file = file.path(fitdir, "stage2-normal-fit.rds"),
  checkpoint_file = file.path(fitdir, "stage2-normal-checkpoint.rds"),
  resume = FALSE, force_refit = TRUE
)
stopifnot(stage2$Y$status == "ok",
          attr(stage2, "prior_family") == "normal",
          isTRUE(stage2$Y$sampling$is_stage2),
          identical(stage2$Y$stage2_parameter_draws$knowledge_prior, original_prior))
filtered <- ns$filter_TF_mcmc_gene_results(
  stage2, output_file = file.path(fitdir, "stage2-normal-filtered.rds")
)
stopifnot(identical(names(filtered), "Y"))
em_input <- ns$build_TF_EM_stage2_input(
  filtered, input, nuisance_storage = "posterior_mean",
  output_file = file.path(fitdir, "stage3-normal-input.rds")
)
stage2_posterior <- stage2$Y$target_tf
posterior_payload <- stage2$Y$target_screening_draws
posterior_delta <- posterior_payload$delta
if (is.null(posterior_delta) &&
    identical(as.integer(stage2$Y$stage2_parameter_draws$source_draw_index), seq_len(posterior_payload$draw_count))) {
  posterior_delta <- stage2$Y$stage2_parameter_draws$beta_target_delta
}
stage2_posterior_covariance <- if (!is.null(posterior_delta)) {
  stats::cov(posterior_payload$overall, posterior_delta)
} else if ("beta_mean_delta_covariance" %in% names(stage2_posterior)) {
  stage2_posterior$beta_mean_delta_covariance
} else (stage2_posterior$beta_disease_sd ^ 2 - stage2_posterior$beta_control_sd ^ 2) / 2
stopifnot(em_input$interface_version == "tf_em_stage2_input_v6",
          em_input$beta_prior_family == "normal",
          em_input$beta_prior_source == "stage2_posterior_normal_approximation",
          isTRUE(all.equal(unname(em_input$beta_prior_mean), stage2_posterior$beta_mean)),
          isTRUE(all.equal(unname(em_input$beta_prior_variance), stage2_posterior$beta_sd ^ 2)),
          isTRUE(all.equal(unname(em_input$beta_prior_sd), stage2_posterior$beta_sd)),
          isTRUE(all.equal(unname(em_input$beta_delta_prior_mean), stage2_posterior$beta_delta_mean)),
          isTRUE(all.equal(unname(em_input$beta_delta_prior_sd), stage2_posterior$beta_delta_sd)),
          isTRUE(all.equal(unname(em_input$beta_delta_prior_covariance), stage2_posterior_covariance)))
em <- ns$run_TF_EM_latent_activity(
  em_input, algorithm = "quadrature_em", nuisance_draw_count = 0L, cores = 1L,
  quadrature_nodes = 101L, max_iter = 3L, min_iter = 2L,
  seed = 20261110L, resume = FALSE, save_traces = TRUE,
  output_file = file.path(fitdir, "stage3-normal-em-fit.rds"),
  checkpoint_file = file.path(fitdir, "stage3-normal-em-checkpoint.rds")
)
stopifnot(all(vapply(em$draw_results, function(x) x$status == "ok", logical(1))),
          nrow(em$activity_summary) == N,
          all(is.finite(em$activity_summary$activity_mean)),
          all(is.finite(em$activity_summary$activity_sd)),
          em$settings$beta_prior_family == "normal")

# Exercise installed-namespace exports on two parallel conditional EM workers.
parallel_input <- ns$build_TF_EM_stage2_input(filtered, input,
  nuisance_storage = "all_draws", nuisance_draw_ids = 1:2)
parallel_em <- ns$run_TF_EM_latent_activity(
  parallel_input, algorithm = "quadrature_em", nuisance_draw_count = 2L, cores = 2L,
  quadrature_nodes = 21L, max_iter = 2L, min_iter = 2L, seed = 20261121L,
  resume = FALSE, output_file = file.path(fitdir, "stage3-posterior-normal-parallel-em-fit.rds"),
  checkpoint_file = file.path(fitdir, "stage3-posterior-normal-parallel-em-checkpoint.rds"))
stopifnot(length(parallel_em$draw_results) == 2L,
          all(vapply(parallel_em$draw_results, function(x) x$status == "ok", logical(1))),
          identical(parallel_input$beta_delta_prior_covariance, em_input$beta_delta_prior_covariance))
message("PASS: two installed-namespace EM workers receive the unchanged full-posterior joint prior.")

report <- data.frame(
  check = c("activity_stage2_compact_stan", "stage3_em_execution_only"),
  status = "ok", prior_family = "normal", cells = N,
  targets = 1L, predictors = input$Y$stan_data$P,
  stage1_gamma = entry$stan_data$gamma, stage1_eta = entry$stan_data$eta,
  stage1_r_dir = entry$stan_data$r_dir,
  focal_stage2_network_prior_mean = original_prior$beta_prior_mean[1],
  focal_stage2_network_prior_variance = original_prior$beta_prior_variance[1],
  focal_stage3_prior_mean = em_input$beta_prior_mean[1],
  focal_stage3_prior_variance = em_input$beta_prior_variance[1],
  focal_stage3_delta_prior_mean = em_input$beta_delta_prior_mean[1],
  focal_stage3_delta_prior_variance = em_input$beta_delta_prior_variance[1],
  focal_stage3_beta_delta_covariance = em_input$beta_delta_prior_covariance[1],
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
cat(paste0("PASS: real activity Stage 2 network Normal, filtering and v6 ",
           "handoff transfer the full posterior Normal moments/covariance without inflation.\n"))

# Check the default Stage 3 sampler on the same real Stage 2 handoff.
# This is an execution/metadata check, not a biological convergence claim.
gibbs <- ns$run_TF_EM_latent_activity(
  em_input, nuisance_draw_count = 0L, cores = 1L,
  seed = 20261120L, resume = FALSE,
  output_file = file.path(fitdir, "stage3-normal-gibbs-fit.rds"),
  checkpoint_file = file.path(fitdir, "stage3-normal-gibbs-checkpoint.rds")
)
stopifnot(gibbs$settings$algorithm == "gibbs_ars",
          gibbs$settings$iter_warmup == 1000L,
          gibbs$settings$iter_sampling == 1000L,
          gibbs$settings$beta_prior_family == "normal",
          all(vapply(gibbs$draw_results, function(x) x$status == "ok", logical(1))),
          all(is.finite(gibbs$activity_summary$activity_mean)),
          all(is.finite(gibbs$beta_summary$beta_mean_mean)))
gibbs_report <- data.frame(
  check = "stage3_default_gibbs_execution_only", status = "ok",
  prior_family = gibbs$settings$beta_prior_family, cells = em_input$N,
  targets = em_input$G, warmup = gibbs$settings$iter_warmup,
  retained_sweeps = gibbs$settings$iter_sampling,
  nuisance_storage = em_input$nuisance_storage
)
write.csv(gibbs_report, file.path(outdir, "activity-gibbs-smoke-summary.csv"), row.names = FALSE)
print(gibbs_report)
cat("PASS: default 1000/1000 Gibbs sampler executes on the real Stage 2 Normal-prior handoff.\n")
