#!/usr/bin/env Rscript

.tf_three_stage_script_dir <- local({
  source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  source_dir <- if (is.null(source_file) || !nzchar(as.character(source_file[[1]]))) {
    getwd()
  } else {
    dirname(normalizePath(source_file[[1]], winslash = "/", mustWork = FALSE))
  }
  function() source_dir
})

# Dependencies are loaded with the package namespace.

tf_three_stage_stop <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

tf_three_stage_read_input <- function(screening_input) {
  if (is.character(screening_input) && length(screening_input) == 1L) {
    if (!file.exists(screening_input)) {
      tf_three_stage_stop("Input file not found: %s", screening_input)
    }
    screening_input <- readRDS(screening_input)
  }
  if (!is.list(screening_input) || !length(screening_input)) {
    tf_three_stage_stop(
      "`screening_input` must be a Stage 1 input list or an RDS path."
    )
  }
  screening_input
}

# Three-stage TF activity pipeline:
#   Stage 1: Laplace-prior meanfield; filter genes and confounders.
#   Stage 2: Normal-prior MCMC; filter genes only, retain all entered confounders.
#   Stage 3: EM latent activity using the retained MCMC models.
run_TF_three_stage_pipeline <- function(
  screening_input,
  output_dir = file.path(.tf_three_stage_script_dir(), "adjustment_output"),
  cores = 4L,
  seed = 123L,
  prescreen_target_interval = 0.90,
  prescreen_confounder_interval = 0.90,
  prescreen_variational_iter = 10000L,
  prescreen_output_samples = 2000L,
  direction_effect = 0.2,
  beta_sd_floor = 0.5,
  stage1_sd_multiplier = 1.5,
  stage2_alpha_prior_sd = 1,
  mcmc_target_interval = 0.90,
  mcmc_confounder_interval = 0.90,
  mcmc_chains = 3L,
  mcmc_iter_warmup = 600L,
  mcmc_iter_sampling = 1200L,
  mcmc_adapt_delta = 0.95,
  mcmc_max_treedepth = 12L,
  mcmc_draws_for_em = 1000L,
  nuisance_draw_count = 0L,
  active_prior_zero = 0.2,
  em_control = list(),
  resume = TRUE,
  force_refit = FALSE,
  force_recompile = FALSE
) {
  screening_input <- tf_three_stage_read_input(screening_input)
  target_tf <- attr(screening_input, "target_tf")
  if (is.null(target_tf) || !nzchar(as.character(target_tf[[1]]))) {
    tf_three_stage_stop("The screening input has no `target_tf` attribute.")
  }
  target_tf <- as.character(target_tf[[1]])
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(output_dir, target_tf)

  paths <- list(
    prescreen_fit = paste0(prefix, "_stage1_laplace_meanfield_fit.rds"),
    prescreen_checkpoint = paste0(
      prefix, "_stage1_laplace_meanfield_checkpoint.rds"
    ),
    prescreen_filtered = paste0(prefix, "_stage1_prescreen_filtered.rds"),
    mcmc_input = paste0(prefix, "_stage2_mcmc_input.rds"),
    mcmc_fit = paste0(prefix, "_stage2_normal_mcmc_fit.rds"),
    mcmc_checkpoint = paste0(prefix, "_stage2_normal_mcmc_checkpoint.rds"),
    mcmc_filtered = paste0(prefix, "_stage2_mcmc_gene_filtered.rds"),
    em_input = paste0(prefix, "_stage3_em_input.rds"),
    em_fit = paste0(prefix, "_stage3_em_fit.rds"),
    em_checkpoint = paste0(prefix, "_stage3_em_checkpoint.rds"),
    em_source_draws = paste0(prefix, "_stage3_selected_stage2_draws.csv"),
    prescreen_gene_summary = paste0(
      prefix, "_stage1_prescreen_gene_summary.csv"
    ),
    prescreen_confounder_summary = paste0(
      prefix, "_stage1_prescreen_confounder_summary.csv"
    ),
    mcmc_gene_summary = paste0(prefix, "_stage2_mcmc_gene_summary.csv")
  )

  prescreen_started <- Sys.time()
  prescreen_fit <- run_TF_stage1_screening_regressions(
    screening_input = screening_input,
    cores = cores,
    inference = "variational",
    prior_family = "laplace",
    target_interval = prescreen_target_interval,
    confounder_interval = prescreen_confounder_interval,
    variational_algorithm = "meanfield",
    variational_iter = prescreen_variational_iter,
    variational_output_samples = prescreen_output_samples,
    stage2_draw_count = 1L,
    seed = seed,
    refresh = 0L,
    output_file = paths$prescreen_fit,
    checkpoint_file = paths$prescreen_checkpoint,
    resume = resume,
    retry_errors = TRUE,
    force_refit = force_refit,
    force_recompile = force_recompile
  )
  prescreen_filtered <- filter_TF_prescreen_results(
    prescreen_results = prescreen_fit,
    target_interval = prescreen_target_interval,
    confounder_interval = prescreen_confounder_interval,
    output_file = paths$prescreen_filtered
  )
  utils::write.csv(
    attr(prescreen_filtered, "filter_summary"),
    paths$prescreen_gene_summary,
    row.names = FALSE
  )
  utils::write.csv(
    attr(prescreen_filtered, "confounder_summary"),
    paths$prescreen_confounder_summary,
    row.names = FALSE
  )
  prescreen_elapsed_seconds <- as.numeric(difftime(
    Sys.time(), prescreen_started, units = "secs"
  ))
  mcmc_input <- build_TF_mcmc_input_from_prescreen(
    prescreen_filtered_results = prescreen_filtered,
    screening_input = screening_input,
    direction_effect = direction_effect,
    beta_sd_floor = beta_sd_floor,
    stage1_sd_multiplier = stage1_sd_multiplier,
    alpha_prior_sd = stage2_alpha_prior_sd,
    output_file = paths$mcmc_input
  )

  mcmc_started <- Sys.time()
  mcmc_fit <- run_TF_stage1_screening_regressions(
    screening_input = mcmc_input,
    cores = cores,
    inference = "mcmc",
    prior_family = "normal",
    target_interval = mcmc_target_interval,
    confounder_interval = mcmc_confounder_interval,
    chains = mcmc_chains,
    iter_warmup = mcmc_iter_warmup,
    iter_sampling = mcmc_iter_sampling,
    adapt_delta = mcmc_adapt_delta,
    max_treedepth = mcmc_max_treedepth,
    stage2_draw_count = mcmc_draws_for_em,
    seed = seed + 100000L,
    refresh = 0L,
    output_file = paths$mcmc_fit,
    checkpoint_file = paths$mcmc_checkpoint,
    resume = resume,
    retry_errors = TRUE,
    force_refit = force_refit,
    force_recompile = force_recompile
  )
  mcmc_filtered <- filter_TF_mcmc_gene_results(
    screening_results = mcmc_fit,
    target_interval = mcmc_target_interval,
    confounder_interval = mcmc_confounder_interval,
    output_file = paths$mcmc_filtered
  )
  utils::write.csv(
    attr(mcmc_filtered, "filter_summary"),
    paths$mcmc_gene_summary,
    row.names = FALSE
  )
  mcmc_elapsed_seconds <- as.numeric(difftime(
    Sys.time(), mcmc_started, units = "secs"
  ))
  if (!length(mcmc_filtered)) {
    tf_three_stage_stop(
      paste0(
        "No target-gene models passed Stage 2 MCMC filtering. ",
        "Inspect `%s` for fit failures or unsupported target-TF effects."
      ),
      paths$mcmc_gene_summary
    )
  }

  nuisance_draw_count <- as.integer(nuisance_draw_count[[1]])
  active_prior_zero <- as.numeric(active_prior_zero[[1]])
  if (!is.finite(active_prior_zero) || active_prior_zero <= 0 ||
      active_prior_zero >= 0.9) {
    tf_three_stage_stop(
      "`active_prior_zero` must be a finite number strictly between 0 and 0.9."
    )
  }
  nuisance_storage <- if (nuisance_draw_count == 0L) {
    "posterior_mean"
  } else {
    "all_draws"
  }
  available_em_draws <- as.integer(
    mcmc_filtered[[1]]$stage2_parameter_draws$draw_count
  )
  if (nuisance_draw_count < 0L || nuisance_draw_count > available_em_draws) {
    tf_three_stage_stop(
      "`nuisance_draw_count` must be between 0 and %d retained Stage 2 draws.",
      available_em_draws
    )
  }
  em_source_draw_ids <- if (nuisance_draw_count == 0L) {
    NULL
  } else {
    # Select upstream MCMC draws before building eta0 so memory scales with the
    # requested EM draw count rather than every retained Stage 2 draw.
    set.seed(as.integer(seed + 200000L))
    sort(sample(seq_len(available_em_draws), size = nuisance_draw_count, replace = FALSE))
  }
  if (!is.null(em_source_draw_ids)) {
    utils::write.csv(
      data.frame(
        em_draw_id = seq_along(em_source_draw_ids),
        stage2_draw_id = em_source_draw_ids,
        stringsAsFactors = FALSE
      ),
      paths$em_source_draws,
      row.names = FALSE
    )
  }
  em_input <- build_TF_EM_stage2_input(
    stage1_filtered_results = mcmc_filtered,
    stage1_screening_input = mcmc_input,
    nuisance_storage = nuisance_storage,
    nuisance_draw_ids = em_source_draw_ids,
    output_file = paths$em_input
  )

  default_em_control <- list(
    stage2_input = em_input,
    nuisance_draw_count = nuisance_draw_count,
    nuisance_draw_ids = if (nuisance_draw_count == 0L) NULL else seq_len(nuisance_draw_count),
    cores = cores,
    kappa = 1,
    active_prior_zero = active_prior_zero,
    quadrature_nodes = 21L,
    max_iter = 30L,
    min_iter = 2L,
    beta_tolerance = 1e-3,
    activity_tolerance = 1e-3,
    objective_tolerance = 1e-8,
    seed = seed + 200000L,
    checkpoint_file = paths$em_checkpoint,
    output_file = paths$em_fit,
    resume = resume,
    retry_errors = TRUE,
    save_traces = TRUE
  )
  duplicated_arguments <- intersect(names(em_control), c(
    "stage2_input", "output_file", "checkpoint_file",
    "active_prior_zero", "active_prior_positive", "hard_zero_expression"
  ))
  if (length(duplicated_arguments)) {
    tf_three_stage_stop(
      "`em_control` cannot override: %s.",
      paste(duplicated_arguments, collapse = ", ")
    )
  }
  em_started <- Sys.time()
  em_fit <- do.call(
    run_TF_EM_latent_activity,
    utils::modifyList(default_em_control, em_control)
  )
  em_elapsed_seconds <- as.numeric(difftime(
    Sys.time(), em_started, units = "secs"
  ))
  timing <- data.frame(
    stage = c("stage1_prescreen_vi", "stage2_normal_mcmc", "stage3_em"),
    elapsed_seconds = c(
      prescreen_elapsed_seconds,
      mcmc_elapsed_seconds,
      em_elapsed_seconds
    ),
    elapsed_minutes = c(
      prescreen_elapsed_seconds,
      mcmc_elapsed_seconds,
      em_elapsed_seconds
    ) / 60,
    stringsAsFactors = FALSE
  )
  timing_file <- paste0(prefix, "_stage_timing.csv")
  if (!file.exists(timing_file)) {
    utils::write.csv(timing, timing_file, row.names = FALSE)
  }

  invisible(list(
    target_tf = target_tf,
    paths = paths,
    prescreen_filter_counts = attr(prescreen_filtered, "filter_counts"),
    mcmc_filter_counts = attr(mcmc_filtered, "filter_counts"),
    em_source_draw_ids = em_source_draw_ids,
    timing = timing,
    em_fit = em_fit
  ))
}
