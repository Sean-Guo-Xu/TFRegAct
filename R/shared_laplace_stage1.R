tf_shared_laplace_stage1_stan_file <- function() {
  .tfregact_stan_file("TF_directional_nb_model.stan")
}

tf_validate_shared_laplace_stage1_data <- function(stan_data) {
  required <- c(
    "N", "P", "Y", "X", "Q", "W", "W_prior_scale",
    "has_condition", "condition", "condition_w_index", "K_batch", "batch",
    "batch_level_design", "log_offset", "confidence", "direction", "gamma",
    "eta", "r_dir", "alpha_prior_sd", "target_tf_index",
    "use_target_interaction", "target_condition_centered",
    "target_interaction_sd"
  )
  missing <- setdiff(required, names(stan_data))
  if (length(missing)) {
    stop(
      "Shared Laplace Stage 1 data are missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  N <- as.integer(stan_data$N)
  P <- as.integer(stan_data$P)
  Q <- as.integer(stan_data$Q)
  K_batch <- as.integer(stan_data$K_batch)
  has_condition <- as.integer(stan_data$has_condition)
  use_interaction <- as.integer(stan_data$use_target_interaction)
  if (N < 1L || P < 1L || Q < 0L || K_batch < 0L ||
      !(has_condition %in% 0:1) || !(use_interaction %in% 0:1) ||
      nrow(stan_data$X) != N || ncol(stan_data$X) != P ||
      nrow(stan_data$W) != N || ncol(stan_data$W) != Q ||
      length(stan_data$W_prior_scale) != Q ||
      length(stan_data$condition) != N ||
      length(stan_data$target_condition_centered) != N ||
      length(stan_data$batch) != N ||
      nrow(stan_data$batch_level_design) != K_batch ||
      ncol(stan_data$batch_level_design) != Q) {
    stop("Shared Laplace Stage 1 data have inconsistent dimensions.", call. = FALSE)
  }
  if (has_condition == 1L) {
    if (stan_data$condition_w_index < 1L ||
        stan_data$condition_w_index > Q ||
        !all(stan_data$condition %in% 0:1)) {
      stop("Invalid condition encoding in shared Laplace Stage 1 data.", call. = FALSE)
    }
  } else if (stan_data$condition_w_index != 0L ||
             any(stan_data$condition != 0L)) {
    stop("Condition-off data must use index 0 and an all-zero condition vector.", call. = FALSE)
  }
  if (use_interaction == 1L && has_condition != 1L) {
    stop("Target interaction requires a two-level condition.", call. = FALSE)
  }
  if (use_interaction == 0L && any(stan_data$target_condition_centered != 0)) {
    stop("Interaction-off data must have a zero target-condition vector.", call. = FALSE)
  }
  invisible(TRUE)
}

tf_fit_shared_laplace_stage1 <- function(
  stan_data,
  inference = c("mcmc", "variational"),
  model = NULL,
  stan_file = NULL,
  force_recompile = FALSE,
  chains = 3L,
  parallel_chains = chains,
  iter_warmup = 600L,
  iter_sampling = 1200L,
  seed = 123L,
  refresh = 0L,
  adapt_delta = 0.95,
  max_treedepth = 12L,
  variational_algorithm = c("meanfield", "fullrank"),
  variational_iter = 10000L,
  variational_output_samples = 2000L,
  ...
) {
  inference <- match.arg(inference)
  variational_algorithm <- match.arg(variational_algorithm)
  tf_validate_shared_laplace_stage1_data(stan_data)

  if (is.null(model)) {
    if (is.null(stan_file) || !nzchar(as.character(stan_file[[1]]))) {
      stan_file <- tf_shared_laplace_stage1_stan_file()
    } else {
      stan_file <- .tfregact_stan_file(as.character(stan_file[[1]]))
    }
    model <- cmdstanr::cmdstan_model(
      stan_file,
      force_recompile = isTRUE(force_recompile)
    )
  }

  if (identical(inference, "mcmc")) {
    return(model$sample(
      data = stan_data,
      chains = chains,
      parallel_chains = parallel_chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      seed = seed,
      refresh = refresh,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      save_warmup = FALSE,
      ...
    ))
  }

  model$variational(
    data = stan_data,
    seed = seed,
    refresh = refresh,
    algorithm = variational_algorithm,
    iter = variational_iter,
    output_samples = variational_output_samples,
    ...
  )
}
