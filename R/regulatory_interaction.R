tf_assert_target_interaction_dependencies <- function() {
  required_functions <- c(
    "tf_model_stop",
    "run_TF_stage2_directional_model"
  )
  missing_functions <- required_functions[!vapply(
    required_functions,
    exists,
    logical(1),
    mode = "function",
    inherits = TRUE
  )]

  if (length(missing_functions) > 0) {
    stop(
      "The unified Stage 2 model dependencies are missing: ",
      paste(missing_functions, collapse = ", "),
      call. = FALSE
    )
  }
}

tf_default_stage2_target_interaction_stan_file <- function(use_batch_model = NULL) {
  .tfregact_stan_file("TF_stage2_directional_nb_model.stan")
}

tf_select_stage2_target_interaction_stan_file <- function(
  stan_file,
  use_batch_model = NULL
) {
  if (is.null(stan_file) || !nzchar(stan_file)) {
    return(tf_default_stage2_target_interaction_stan_file())
  }
  .tfregact_stan_file(stan_file)
}

run_TF_stage2_target_interaction_model <- function(
  analysis_object,
  stage1_fit_result,
  stan_file,
  direction_effect,
  beta_sd_floor,
  stage1_sd_multiplier,
  target_interaction_sd,
  control_level,
  disease_level,
  confidence_min,
  confidence_max,
  chains,
  parallel_chains,
  iter_warmup,
  iter_sampling,
  seed,
  refresh,
  compute_loo,
  force_recompile,
  nuisance_prior_scale = 1,
  ...
) {
  tf_assert_target_interaction_dependencies()

  if (!is.numeric(target_interaction_sd) ||
      length(target_interaction_sd) != 1L ||
      !is.finite(target_interaction_sd) ||
      target_interaction_sd <= 0) {
    tf_model_stop("`target_interaction_sd` must be one finite positive number.")
  }

  stan_file <- tf_select_stage2_target_interaction_stan_file(stan_file)
  result <- run_TF_stage2_directional_model(
    analysis_object = analysis_object,
    stage1_fit_result = stage1_fit_result,
    stan_file = stan_file,
    direction_effect = direction_effect,
    beta_sd_floor = beta_sd_floor,
    stage1_sd_multiplier = stage1_sd_multiplier,
    control_level = control_level,
    disease_level = disease_level,
    confidence_min = confidence_min,
    confidence_max = confidence_max,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    refresh = refresh,
    compute_loo = compute_loo,
    force_recompile = force_recompile,
    nuisance_prior_scale = nuisance_prior_scale,
    target_interaction = TRUE,
    target_interaction_sd = target_interaction_sd,
    ...
  )
  result$target_interaction_model <- TRUE
  result
}
