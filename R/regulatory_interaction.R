tf_assert_target_interaction_dependencies <- function() {
  required_functions <- c(
    "tf_model_stop",
    "tf_model_require_pkg",
    "tf_prepare_stage2_stan_data",
    "tf_get_cmdstan_model",
    "tf_compute_loo"
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
      "Source `run_TF_directional_model.R` before the target-interaction model file. Missing: ",
      paste(missing_functions, collapse = ", "),
      call. = FALSE
    )
  }
}

tf_default_stage2_target_interaction_stan_file <- function(use_batch_model = TRUE) {
  stan_file <- if (isTRUE(use_batch_model)) {
    "TF_stage2_target_interaction_directional_nb_model.stan"
  } else {
    "TF_stage2_target_interaction_directional_nb_model_nobatch.stan"
  }
  .tfregact_stan_file(stan_file)
}

tf_select_stage2_target_interaction_stan_file <- function(stan_file, use_batch_model) {
  if (is.null(stan_file) || !nzchar(stan_file)) {
    return(tf_default_stage2_target_interaction_stan_file(
      use_batch_model = use_batch_model
    ))
  }

  if (!isTRUE(use_batch_model) &&
      basename(stan_file) == "TF_stage2_target_interaction_directional_nb_model.stan") {
    return(tf_default_stage2_target_interaction_stan_file(
      use_batch_model = FALSE
    ))
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
  ...
) {
  tf_assert_target_interaction_dependencies()
  tf_model_require_pkg("cmdstanr")

  stan_data <- tf_prepare_stage2_stan_data(
    analysis_object = analysis_object,
    stage1_fit_result = stage1_fit_result,
    direction_effect = direction_effect,
    beta_sd_floor = beta_sd_floor,
    stage1_sd_multiplier = stage1_sd_multiplier,
    control_level = control_level,
    disease_level = disease_level,
    confidence_min = confidence_min,
    confidence_max = confidence_max
  )

  if (!isTRUE(attr(stan_data, "condition_model")) || is.null(stan_data$condition)) {
    tf_model_stop(
      "The target-interaction stage 2 model requires a valid two-level condition in `analysis_object$sample`."
    )
  }
  if (!is.numeric(target_interaction_sd) ||
      length(target_interaction_sd) != 1 ||
      !is.finite(target_interaction_sd) ||
      target_interaction_sd <= 0) {
    tf_model_stop("`target_interaction_sd` must be one finite positive number.")
  }
  stan_data$target_interaction_sd <- as.numeric(target_interaction_sd)

  use_batch_model <- isTRUE(attr(stan_data, "use_batch_model"))
  stan_file <- tf_select_stage2_target_interaction_stan_file(stan_file, use_batch_model)

  model <- tf_get_cmdstan_model(
    stan_file,
    force_recompile = force_recompile,
    include_log_lik = isTRUE(compute_loo)
  )
  fit <- model$sample(
    data = stan_data,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    refresh = refresh,
    ...
  )
  loo_result <- if (isTRUE(compute_loo)) tf_compute_loo(fit) else NULL

  list(
    fit = fit,
    loo = loo_result,
    stan_data = stan_data,
    stan_file = stan_file,
    target = attr(stan_data, "target"),
    target_tf = attr(stan_data, "target_tf"),
    target_tf_index = attr(stan_data, "target_tf_index"),
    feature_names = attr(stan_data, "feature_names"),
    cell_names = attr(stan_data, "cell_names"),
    use_batch_model = use_batch_model,
    condition_model = TRUE,
    target_interaction_model = TRUE,
    control_level = attr(stan_data, "control_level"),
    disease_level = attr(stan_data, "disease_level"),
    stage1_beta_summary = attr(stan_data, "stage1_beta_summary")
  )
}
