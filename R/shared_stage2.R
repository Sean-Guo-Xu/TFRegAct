tf_shared_stage2_stop <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

tf_shared_stage2_model_version <- function() {
  "stage1_stage2_shared_knowledge_laplace_v7"
}

tf_stage2_warn_legacy_prior_arguments <- function(
  direction_effect = NULL, beta_sd_floor = NULL, stage1_sd_multiplier = NULL
) {
  supplied <- c(direction_effect = !is.null(direction_effect),
    beta_sd_floor = !is.null(beta_sd_floor),
    stage1_sd_multiplier = !is.null(stage1_sd_multiplier))
  if (any(supplied)) {
    warning(paste0("`", paste(names(supplied)[supplied], collapse = "`, `"),
      "` are deprecated and ignored; Stage 2 reuses the Stage 1 knowledge prior."),
      call. = FALSE)
  }
  invisible(NULL)
}

# Describe the original network-informed prior; posterior moments are not inputs.
tf_stage2_knowledge_prior <- function(stan_data) {
  fields <- c("P", "confidence", "direction", "gamma", "eta", "r_dir")
  missing <- setdiff(fields, names(stan_data))
  if (length(missing)) {
    tf_shared_stage2_stop("Stage 1 knowledge-prior metadata are missing: %s.",
      paste(missing, collapse = ", "))
  }
  P <- stan_data$P
  confidence <- as.numeric(stan_data$confidence)
  direction <- stan_data$direction
  hyperparameters <- unlist(stan_data[c("gamma", "eta", "r_dir")], use.names = FALSE)
  if (length(P) != 1L || !is.finite(P) || P < 1L || P != as.integer(P) ||
      length(confidence) != P || length(direction) != P ||
      any(!is.finite(confidence)) || any(confidence < 0) ||
      anyNA(direction) || any(!direction %in% c(-1L, 0L, 1L)) ||
      any(lengths(stan_data[c("gamma", "eta", "r_dir")]) != 1L) ||
      any(!is.finite(hyperparameters)) || stan_data$gamma <= 0 ||
      stan_data$eta < 0 || stan_data$r_dir < 1) {
    tf_shared_stage2_stop(
      "Invalid confidence, direction or scale in the Stage 1 knowledge prior."
    )
  }
  direction <- as.integer(direction)
  scale <- stan_data$gamma * (pmax(confidence, 1) / 10) ^ stan_data$eta
  lambda_fav <- 1 / scale
  lambda_positive <- lambda_fav * ifelse(direction == -1L, stan_data$r_dir, 1)
  lambda_negative <- lambda_fav * ifelse(direction == 1L, stan_data$r_dir, 1)
  if (any(!is.finite(scale)) || any(scale <= 0) ||
      any(!is.finite(lambda_positive)) || any(!is.finite(lambda_negative))) {
    tf_shared_stage2_stop("Knowledge-prior scales and rates must be finite and positive.")
  }
  list(
    prior_family = "laplace", confidence = confidence, direction = as.integer(direction),
    gamma = as.numeric(stan_data$gamma), eta = as.numeric(stan_data$eta),
    r_dir = as.numeric(stan_data$r_dir), beta_prior_scale = scale,
    beta_lambda_positive = lambda_positive, beta_lambda_negative = lambda_negative
  )
}

tf_shared_stage2_data_fields <- function() {
  c(
    "N", "P", "Y", "X", "Q", "W", "W_prior_scale",
    "has_condition", "condition_w_index", "condition",
    "K_batch", "batch", "batch_level_design", "log_offset",
    "confidence", "direction", "gamma", "eta", "r_dir", "alpha_prior_sd",
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
      !is.finite(stan_data$alpha_prior_sd) || stan_data$alpha_prior_sd <= 0) {
    tf_shared_stage2_stop("Shared Stage 2 data have inconsistent dimensions or priors.")
  }
  tf_stage2_knowledge_prior(stan_data)
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
    rep(0, data_for_stan$P)
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
