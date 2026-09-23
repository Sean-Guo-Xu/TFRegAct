# Biological-sample hierarchy for the regulatory Stage 2 branch.
# Keep X on the Stage 1 scale: within-sample centering changes its origin only.
tf_regulatory_stage2_route <- function(data_type, condition_column, stage2_model) {
  data_type <- match.arg(data_type, c("single_cell", "bulk"))
  stage2_model <- match.arg(stage2_model, c("no_interaction", "target_interaction"))
  if (is.null(condition_column)) return("no_interaction")
  if (identical(data_type, "single_cell")) return("hierarchical")
  stage2_model
}

tf_hierarchical_sample_design <- function(stan_data, metadata, sample_column) {
  if (is.null(sample_column) || length(sample_column) != 1L ||
      is.na(sample_column) || !(sample_column %in% names(metadata))) {
    tf_model_stop("Single-cell Stage 2 with condition requires `sample_column`, identifying biological samples, separately from condition and batch.")
  }
  cells <- attr(stan_data, "cell_names")
  if (anyDuplicated(rownames(metadata)) || anyNA(match(cells, rownames(metadata)))) {
    tf_model_stop("Biological-sample metadata do not align with model observations.")
  }
  md <- metadata[match(cells, rownames(metadata)), , drop = FALSE]
  samples <- trimws(as.character(md[[sample_column]]))
  if (anyNA(samples) || any(!nzchar(samples))) tf_model_stop("Biological-sample IDs must not be missing or empty.")
  if (stan_data$has_condition != 1L) tf_model_stop("Hierarchical Stage 2 requires a valid two-level condition.")
  levels <- sort(unique(samples))
  sid <- match(samples, levels)
  condition <- vapply(seq_along(levels), function(s) {
    value <- unique(stan_data$condition[sid == s])
    if (length(value) != 1L) tf_model_stop("Biological sample `%s` spans multiple conditions.", levels[s])
    as.integer(value)
  }, integer(1))
  if (any(tabulate(condition + 1L, nbins = 2L) < 2L)) {
    tf_model_stop("Hierarchical slope comparison requires at least two biological samples per condition.")
  }
  variance <- vapply(seq_along(levels), function(s) {
    x <- stan_data$X[sid == s, stan_data$target_tf_index]
    if (length(x) < 2L) 0 else stats::var(x)
  }, numeric(1))
  informed <- is.finite(variance) & variance > 1e-8
  if (any(vapply(0:1, function(k) sum(informed[condition == k]) < 2L, logical(1)))) {
    tf_model_stop("Hierarchical slope comparison requires at least two samples with within-sample focal-TF variation in each condition.")
  }
  if (any(!informed)) warning("Samples without focal-TF variation are retained through partial pooling: ",
    paste(levels[!informed], collapse = ", "), call. = FALSE)
  if (qr(cbind(1, stan_data$W))$rank != stan_data$Q + 1L) {
    tf_model_stop("Condition/batch nuisance design is not full rank.")
  }
  list(S = length(levels), sample_id = as.integer(sid), sample_condition = condition,
    sample_map = data.frame(sample = levels, condition = condition,
      cells = tabulate(sid, nbins = length(levels)), within_variance = variance,
      slope_informed_by_expression = informed, stringsAsFactors = FALSE))
}

run_TF_stage2_hierarchical_model <- function(
  analysis_object, stage1_fit_result, metadata, sample_column,
  direction_effect = 0.2, beta_sd_floor = 0.5, stage1_sd_multiplier = 1.5,
  control_level = NULL, disease_level = NULL, confidence_min = 1, confidence_max = 10,
  nuisance_prior_scale = 1, target_interaction_sd = 0.5,
  sample_intercept_sd_scale = 1, sample_slope_sd_scale = 1,
  chains = 3L, parallel_chains = chains, iter_warmup = 500L,
  iter_sampling = 900L, seed = 123L, refresh = 50L,
  adapt_delta = 0.95, max_treedepth = 12L, compute_loo = FALSE,
  force_recompile = FALSE
) {
  for (scale in list(sample_intercept_sd_scale, sample_slope_sd_scale)) {
    if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) || scale <= 0) {
      tf_model_stop("Hierarchical prior scales must be finite positive numbers.")
    }
  }
  d <- tf_prepare_stage2_stan_data(
    analysis_object, stage1_fit_result, direction_effect, beta_sd_floor,
    stage1_sd_multiplier, control_level, disease_level, confidence_min,
    confidence_max, nuisance_prior_scale, target_interaction = TRUE,
    target_interaction_sd = target_interaction_sd)
  design <- tf_hierarchical_sample_design(d, metadata, sample_column)
  d$S <- design$S
  d$sample_id <- design$sample_id
  d$sample_condition <- design$sample_condition
  d$sample_intercept_sd_scale <- sample_intercept_sd_scale
  d$sample_slope_sd_scale <- sample_slope_sd_scale
  d$save_log_lik <- as.integer(compute_loo)
  fields <- c("N", "P", "Y", "X", "log_offset", "S", "sample_id",
    "has_condition", "sample_condition", "use_target_interaction", "target_tf_index",
    "Q", "W", "W_prior_scale", "beta_prior_mean", "beta_prior_sd", "alpha_prior_sd",
    "target_interaction_sd", "sample_intercept_sd_scale", "sample_slope_sd_scale", "save_log_lik")
  stan_file <- .tfregact_stan_file("TF_stage2_hierarchical_nb_model.stan")
  model <- tf_get_cmdstan_model(stan_file, force_recompile, include_log_lik = TRUE)
  init <- function(chain_id = 1L) list(alpha = 0, beta = as.numeric(d$beta_init),
    zeta = rep(0, d$Q), beta_target_delta_free = 0, log_phi = 0,
    sample_intercept_raw = rep(0, d$S), sample_slope_raw = rep(0, d$S),
    sample_intercept_sd = 0.2, sample_slope_sd = 0.2)
  fit <- model$sample(data = d[fields], chains = chains, parallel_chains = parallel_chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling, seed = seed,
    refresh = refresh, adapt_delta = adapt_delta, max_treedepth = max_treedepth,
    init = init, save_warmup = FALSE)
  if (any(fit$return_codes() != 0L)) tf_model_stop("One or more hierarchical Stage 2 chains failed.")
  # Cell-conditional LOO is not leave-one-biological-sample-out validation.
  loo_result <- if (isTRUE(compute_loo)) tf_compute_loo(fit) else NULL
  list(fit = fit, loo = loo_result,
    loo_unit = if (isTRUE(compute_loo)) "cell_conditional_on_sample_effects" else NULL,
    stan_data = d, stan_file = stan_file, target = attr(d, "target"),
    feature_names = attr(d, "feature_names"), cell_names = attr(d, "cell_names"),
    use_batch_model = isTRUE(attr(d, "use_batch_model")), condition_model = TRUE,
    target_interaction_model = TRUE, hierarchical_model = TRUE,
    target_tf = attr(d, "target_tf"), target_tf_index = d$target_tf_index,
    control_level = attr(d, "control_level"), disease_level = attr(d, "disease_level"),
    W_names = attr(d, "W_names"), stage1_beta_summary = attr(d, "stage1_beta_summary"),
    stage2_prior = attr(d, "stage2_prior"), sample_column = sample_column,
    sample_map = design$sample_map, coefficient_scale = "within_sample_on_input_expression_scale")
}
