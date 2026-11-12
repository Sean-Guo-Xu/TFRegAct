# Small executable-model check; run from the package root after installation.
# This checks prior routing and execution, not biological accuracy.
.libPaths(c("D:/Probscorer/.validation/stage2-knowledge-prior/library", "D:/R", .libPaths()))
library(TFRegAct)
ns <- asNamespace("TFRegAct")
set.seed(20261108)
N <- 96L
cells <- paste0("cell", seq_len(N))
sample_id <- rep(paste0("s", 1:8), each = 12)
condition <- rep(c("Normal", "AAA"), each = N / 2)
x <- rbind(T = runif(N, 0, 2), V = runif(N, 0, 2))
colnames(x) <- cells
shift <- as.numeric(condition == "AAA") - 0.5
sample_effect <- rep(rnorm(8, 0, 0.15), each = 12)
y <- rnbinom(N, mu = exp(1 + (0.6 + 0.2 * shift) * x[1, ] - 0.2 * x[2, ] +
  sample_effect), size = 8)
a <- ns$create_TF_analysis_object(x, target = c("T", "Y"),
  sample = condition, batch = NULL, confidence = c(4, 7),
  direction = c("activation", "unknown"), libsize = rep(1000, N), Y_exp = y)
metadata <- data.frame(donor = sample_id, row.names = cells)
stage1 <- ns$run_TF_directional_model(a,
  stan_file = normalizePath("inst/stan/TF_directional_nb_model.stan"),
  gamma = 1.6, eta = 0.7, r_dir = 4, confidence_min = 1, confidence_max = 10,
  chains = 2, parallel_chains = 2, iter_warmup = 300, iter_sampling = 500,
  seed = 11, refresh = 0, force_recompile = FALSE,
  control_level = "Normal", disease_level = "AAA", adapt_delta = 0.95,
  max_treedepth = 12)
base_args <- list(analysis_object = a, stage1_fit_result = stage1,
  stan_file = normalizePath("inst/stan/TF_stage2_directional_nb_model.stan"),
  control_level = "Normal", disease_level = "AAA", confidence_min = 1,
  confidence_max = 10, chains = 2, parallel_chains = 2,
  iter_warmup = 300, iter_sampling = 500, seed = 12, refresh = 0,
  compute_loo = FALSE, force_recompile = FALSE,
  adapt_delta = 0.95, max_treedepth = 12)
plain <- do.call(ns$run_TF_stage2_directional_model, base_args)
interaction_args <- c(base_args, list(target_interaction_sd = 0.5))
interaction <- do.call(ns$run_TF_stage2_target_interaction_model, interaction_args)
hierarchical <- ns$run_TF_stage2_hierarchical_model(a, stage1, metadata, "donor",
  control_level = "Normal", disease_level = "AAA", chains = 2, parallel_chains = 2,
  iter_warmup = 300, iter_sampling = 500, seed = 13, refresh = 0,
  adapt_delta = 0.95, max_treedepth = 12)
fits <- list(stage1 = stage1, plain = plain, interaction = interaction,
  hierarchical = hierarchical)
report <- lapply(names(fits), function(name) {
  result <- fits[[name]]
  stopifnot(all(result$fit$return_codes() == 0L))
  prior <- ns$tf_stage2_knowledge_prior(result$stan_data)
  stopifnot(identical(prior, ns$tf_stage2_knowledge_prior(stage1$stan_data)))
  summary <- result$fit$summary(c("beta", "beta_target_delta"))
  diagnostic <- result$fit$diagnostic_summary()
  data.frame(model = name, draws = nrow(result$fit$draws("beta", format = "matrix")),
    max_rhat = max(summary$rhat, na.rm = TRUE),
    divergences = sum(diagnostic$num_divergent),
    max_treedepth_hits = sum(diagnostic$num_max_treedepth))
})
report <- do.call(rbind, report)
print(report)
outdir <- "D:/Probscorer/.validation/stage2-knowledge-prior"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
write.csv(report, file.path(outdir, "stan-smoke-summary.csv"), row.names = FALSE)
saveRDS(list(report = report, prior = ns$tf_stage2_knowledge_prior(stage1$stan_data)),
  file.path(outdir, "stan-smoke-summary.rds"))
cat("PASS: Stage 1 and all regulatory Stage 2 routes execute with identical knowledge priors.\n")
