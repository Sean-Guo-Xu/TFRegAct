tf_shared_stage2_stop <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

tf_shared_stage2_model_version <- function() {
  "stage1_laplace_stage2_shared_empirical_bayes_v6"
}

tf_stage2_prior_from_stage1 <- function(
  stage1_beta_mean,
  stage1_beta_sd,
  direction,
  direction_effect = 0.2,
  beta_sd_floor = 0.5,
  stage1_sd_multiplier = 1.5
) {
  stage1_beta_mean <- as.numeric(stage1_beta_mean)
  stage1_beta_sd <- as.numeric(stage1_beta_sd)
  direction <- as.integer(direction)
  P <- length(stage1_beta_mean)
  if (P < 1L || length(stage1_beta_sd) != P || length(direction) != P ||
      any(!is.finite(stage1_beta_mean)) || any(!is.finite(stage1_beta_sd)) ||
      any(stage1_beta_sd < 0) || any(!direction %in% c(-1L, 0L, 1L))) {
    tf_shared_stage2_stop("Invalid Stage 1 beta summary or direction vector.")
  }
  hyperparameters <- c(direction_effect, beta_sd_floor, stage1_sd_multiplier)
  if (any(!is.finite(hyperparameters)) || direction_effect < 0 ||
      beta_sd_floor <= 0 || stage1_sd_multiplier <= 0) {
    tf_shared_stage2_stop(
      "`direction_effect` must be nonnegative; prior SD floor and multiplier must be positive."
    )
  }
  list(
    beta_prior_mean = stage1_beta_mean + direction_effect * direction,
    beta_prior_sd = pmax(beta_sd_floor, stage1_sd_multiplier * stage1_beta_sd),
    stage1_beta_mean = stage1_beta_mean,
    stage1_beta_sd = stage1_beta_sd,
    direction = direction,
    direction_effect = as.numeric(direction_effect),
    beta_sd_floor = as.numeric(beta_sd_floor),
    stage1_sd_multiplier = as.numeric(stage1_sd_multiplier)
  )
}

tf_shared_stage2_data_fields <- function() {
  c(
    "N", "P", "Y", "X", "Q", "W", "W_prior_scale",
    "has_condition", "condition_w_index", "condition",
    "K_batch", "batch", "batch_level_design", "log_offset",
    "beta_prior_mean", "beta_prior_sd", "alpha_prior_sd",
    "target_tf_index", "use_target_interaction",
    "target_condition_centered", "target_interaction_sd"
  )
}

tf_validate_shared_stage2_data <- function(stan_data) {
  required <- tf_shared_stage2_data_fields()
  missing <- setdiff(required, names(stan_data))
  if (length(missing)) {
    tf_shared_stage2_stop(
      "Shared Stage 2 data are missing: %s.", paste(missing, collapse = ", ")
    )
  }
  N <- as.integer(stan_data$N)
  P <- as.integer(stan_data$P)
  Q <- as.integer(stan_data$Q)
  K_batch <- as.integer(stan_data$K_batch)
  if (N < 1L || P < 1L || Q < 0L || K_batch < 0L ||
      length(stan_data$Y) != N || nrow(stan_data$X) != N ||
      ncol(stan_data$X) != P || nrow(stan_data$W) != N ||
      ncol(stan_data$W) != Q || length(stan_data$W_prior_scale) != Q ||
      length(stan_data$condition) != N || length(stan_data$batch) != N ||
      nrow(stan_data$batch_level_design) != K_batch ||
      ncol(stan_data$batch_level_design) != Q ||
      length(stan_data$log_offset) != N ||
      length(stan_data$beta_prior_mean) != P ||
      length(stan_data$beta_prior_sd) != P ||
      any(!is.finite(stan_data$beta_prior_mean)) ||
      any(!is.finite(stan_data$beta_prior_sd)) ||
      any(stan_data$beta_prior_sd <= 0) ||
      !is.finite(stan_data$alpha_prior_sd) || stan_data$alpha_prior_sd <= 0) {
    tf_shared_stage2_stop("Shared Stage 2 data have inconsistent dimensions or priors.")
  }
  has_condition <- as.integer(stan_data$has_condition)
  use_interaction <- as.integer(stan_data$use_target_interaction)
  if (!has_condition %in% 0:1 || !use_interaction %in% 0:1 ||
      any(!stan_data$condition %in% 0:1) ||
      length(stan_data$target_condition_centered) != N ||
      (has_condition == 0L && (stan_data$condition_w_index != 0L ||
        any(stan_data$condition != 0L))) ||
      (has_condition == 1L && (stan_data$condition_w_index < 1L ||
        stan_data$condition_w_index > Q)) ||
      (use_interaction == 1L && has_condition == 0L) ||
      (use_interaction == 0L && any(stan_data$target_condition_centered != 0))) {
    tf_shared_stage2_stop("Invalid optional-condition or target-interaction encoding.")
  }
  if (stan_data$target_tf_index < 1L || stan_data$target_tf_index > P ||
      !is.finite(stan_data$target_interaction_sd) ||
      stan_data$target_interaction_sd <= 0) {
    tf_shared_stage2_stop("Invalid target-TF index or interaction prior scale.")
  }
  if (!is.null(stan_data$beta_init) &&
      (length(stan_data$beta_init) != P ||
       any(!is.finite(stan_data$beta_init)))) {
    tf_shared_stage2_stop("`beta_init` must be a finite vector of length P.")
  }
  invisible(TRUE)
}

tf_fit_shared_stage2 <- function(
  stan_data,
  model,
  chains,
  parallel_chains,
  iter_warmup,
  iter_sampling,
  seed,
  refresh,
  adapt_delta,
  max_treedepth,
  init = NULL,
  ...
) {
  tf_validate_shared_stage2_data(stan_data)
  data_for_stan <- stan_data[tf_shared_stage2_data_fields()]
  beta_init <- if (is.null(stan_data$beta_init)) {
    as.numeric(data_for_stan$beta_prior_mean)
  } else {
    as.numeric(stan_data$beta_init)
  }
  if (is.null(init)) {
    init <- function(chain_id = 1L) {
      list(
        alpha = 0,
        beta = beta_init,
        zeta = rep(0, data_for_stan$Q),
        beta_target_delta_free = if (data_for_stan$use_target_interaction == 1L) 0 else numeric(0),
        log_phi = 0
      )
    }
  }
  model$sample(
    data = data_for_stan,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    refresh = refresh,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    save_warmup = FALSE,
    init = init,
    ...
  )
}
