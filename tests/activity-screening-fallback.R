# Constructed posterior vectors only. No package installation, model or RNG call.
r_dir <- c("R", "../R", "TFRegAct/R")
r_dir <- r_dir[vapply(r_dir, function(path) {
  file.exists(file.path(path, "activity_screening_fallback.R"))
}, logical(1))][1L]
stopifnot(length(r_dir) == 1L, !is.na(r_dir))
scope <- new.env(parent = globalenv())
for (file in c("activity_stage1.R", "activity_screening_fallback.R",
               "activity_prescreen_filter.R", "activity_stage2_filter.R")) {
  source(file.path(r_dir, file), local = scope)
}

expect_error <- function(expression, pattern = NULL) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"))
  if (!is.null(pattern)) stopifnot(grepl(pattern, conditionMessage(error), ignore.case = TRUE))
}

fixed <- stats::qnorm(stats::ppoints(1001L))
interval <- function(values, level = 0.9) {
  stats::quantile(values, c((1 - level) / 2, (1 + level) / 2), names = FALSE, type = 7)
}
parameter <- function(tf, values, role, control = NULL, disease = NULL) {
  bounds <- interval(values)
  control_bounds <- if (is.null(control)) c(NA_real_, NA_real_) else interval(control)
  disease_bounds <- if (is.null(disease)) c(NA_real_, NA_real_) else interval(disease)
  data.frame(tf = tf, role = role, beta_index = 1L,
    beta_mean = mean(values), beta_sd = stats::sd(values), screening_interval_level = 0.9,
    screening_interval_lower = bounds[1L], screening_interval_upper = bounds[2L],
    screening_interval_excludes_zero = bounds[1L] > 0 || bounds[2L] < 0,
    beta_control_interval_lower = control_bounds[1L], beta_control_interval_upper = control_bounds[2L],
    beta_disease_interval_lower = disease_bounds[1L], beta_disease_interval_upper = disease_bounds[2L],
    control_interval_excludes_zero = !is.null(control) && (control_bounds[1L] > 0 || control_bounds[2L] < 0),
    disease_interval_excludes_zero = !is.null(disease) && (disease_bounds[1L] > 0 || disease_bounds[2L] < 0))
}

model <- function(overall = fixed + 1.45, control = NULL, disease = NULL) {
  n <- length(overall)
  has_condition <- !is.null(control) && !is.null(disease)
  payload <- list(interface_version = "tf_target_screening_draws_v1", draw_count = n, overall = overall)
  if (has_condition) {
    payload$control <- control
    payload$disease <- disease
  }
  list(status = "ok", target_tf_name = "T", condition_model = has_condition,
    target_tf = parameter("T", overall, "target_tf", control, disease),
    confounders = rbind(parameter("C_KEEP", fixed + 3, "confounder"),
                       parameter("C_DROP", fixed + 1.45, "confounder")),
    target_screening_draws = payload,
    stage2_parameter_draws = list(interface_version = "tf_em_stage2_parameter_draws_v2",
      draw_count = n, source_draw_index = seq_len(n), target_tf_index = 1L,
      beta = cbind(T = overall, C_KEEP = fixed + 3, C_DROP = fixed + 1.45),
      beta_target_control = if (has_condition) control else overall,
      beta_target_disease = if (has_condition) disease else overall),
    diagnostics = list(posterior_draws = n), control_level = "control", disease_level = "disease")
}

failed <- list(status = "error", target_tf_name = "T", error = "Synthetic posterior retrieval failure")
results <- function(stage, models) {
  attr(models, "inference") <- if (stage == "prescreen") "variational" else "mcmc"
  attr(models, "prior_family") <- "laplace"
  attr(models, "target_tf") <- "T"
  attr(models, "model_version") <- "constructed_test"
  if (stage == "mcmc") attr(models, "input_pipeline_stage") <- "mcmc_input"
  models
}

run_filter <- function(stage, models, fallback = NULL) {
  original <- serialize(models, NULL)
  path <- tempfile(fileext = ".rds")
  on.exit({
    unlink(path)
    stopifnot(identical(serialize(models, NULL), original))
  })
  filter <- if (stage == "prescreen") scope$filter_TF_prescreen_results else scope$filter_TF_mcmc_gene_results
  output <- filter(models, target_interval = 0.9, confounder_interval = 0.9,
                   output_file = path, fallback_target_interval = fallback)
  stopifnot(identical(serialize(models, NULL), original), file.exists(path),
            identical(readRDS(path), output))
  output
}

check_audit <- function(output, fallback, attempted, retained, initial = 0L) {
  stopifnot(attr(output, "requested_target_interval") == 0.9,
    isTRUE(all.equal(attr(output, "fallback_target_interval"), fallback)),
    identical(attr(output, "fallback_attempted"), attempted),
    identical(attr(output, "fallback_applied"), attempted),
    identical(attr(output, "fallback_succeeded"), attempted && retained > 0L),
    attr(output, "initial_retained_count") == initial,
    attr(output, "final_retained_count") == retained,
    attr(output, "target_interval") == if (attempted) fallback else 0.9,
    attr(output, "confounder_interval") == 0.9,
    attr(output, "filter_counts")[["gene_models_retained"]] == retained)
  summary <- attr(output, "filter_summary")
  stopifnot(all(c("requested_target_interval", "effective_target_interval", "fallback_target_interval",
                  "initial_target_gene_retained", "final_target_gene_retained",
                  "initial_target_supported_in") %in% names(summary)),
    sum(summary$initial_target_gene_retained) == initial,
    sum(summary$final_target_gene_retained) == retained,
    identical(summary$target_gene_retained, summary$final_target_gene_retained))
  initial_summary <- attr(output, "initial_filter_summary")
  stopifnot(sum(initial_summary$target_gene_retained) == initial,
    attr(output, "initial_filter_counts")[["gene_models_retained"]] == initial)
}

for (stage in c("prescreen", "mcmc")) {
  weak <- results(stage, list(Y = model()))
  disabled <- run_filter(stage, weak)
  stopifnot(length(disabled) == 0L)
  check_audit(disabled, NULL, FALSE, 0L)

  for (fallback in c(0.8, 0.81)) {
    output <- run_filter(stage, weak, fallback)
    stopifnot(identical(names(output), "Y"))
    check_audit(output, fallback, TRUE, 1L)
    bounds <- interval(fixed + 1.45, fallback)
    stopifnot(isTRUE(all.equal(unname(c(output$Y$target_tf$screening_interval_lower,
      output$Y$target_tf$screening_interval_upper)), bounds)),
      output$Y$target_tf$screening_interval_level == fallback)
    if (stage == "prescreen") {
      stopifnot(identical(output$Y$confounders$tf, "C_KEEP"),
        all(output$Y$confounders$screening_interval_level == 0.9),
        identical(attr(output, "confounder_summary")$interval_level, c(0.9, 0.9)))
    } else {
      stopifnot(identical(output$Y$confounders, weak$Y$confounders),
        identical(output$Y$stage2_parameter_draws, weak$Y$stage2_parameter_draws))
    }
  }

  # These fixed vectors separate an 80% interval from an 81% interval.
  boundary <- results(stage, list(Y = model(fixed + 1.3)))
  at80 <- run_filter(stage, boundary, 0.8)
  at81 <- run_filter(stage, boundary, 0.81)
  check_audit(at80, 0.8, TRUE, 1L)
  check_audit(at81, 0.81, TRUE, 0L)
  stopifnot(length(at80) == 1L, length(at81) == 0L)

  existing <- results(stage, list(STRONG = model(fixed + 2), WEAK = model()))
  baseline <- run_filter(stage, existing)
  unchanged <- run_filter(stage, existing, 0.8)
  check_audit(unchanged, 0.8, FALSE, 1L, initial = 1L)
  stopifnot(identical(names(unchanged), "STRONG"), identical(unchanged$STRONG, baseline$STRONG))

  # Failed fits remain excluded and cannot supply fake support on rescreening.
  errors_only <- run_filter(stage, results(stage, list(BAD = failed)), 0.8)
  check_audit(errors_only, 0.8, FALSE, 0L)
  partial <- run_filter(stage, results(stage, list(Y = model(), BAD = failed)), 0.8)
  check_audit(partial, 0.8, TRUE, 1L)
  audit <- attr(partial, "filter_summary")
  stopifnot(identical(names(partial), "Y"),
    audit$fit_status[audit$target_gene == "BAD"] == "error",
    !audit$target_gene_retained[audit$target_gene == "BAD"],
    audit$fit_error[audit$target_gene == "BAD"] == failed$error)

  # A successful fit still receives only one weaker screening pass if empty.
  empty <- run_filter(stage, results(stage, list(Y = model(fixed))), 0.8)
  check_audit(empty, 0.8, TRUE, 0L)
  touch_zero <- run_filter(stage, results(stage, list(Y = model(c(rep(0, 201L), rep(1, 800L))))), 0.8)
  check_audit(touch_zero, 0.8, TRUE, 0L)

  cases <- list(control_only = list(control = fixed + 1.45, disease = fixed),
                disease_only = list(control = fixed, disease = fixed - 1.45),
                both = list(control = fixed + 1.45, disease = fixed - 1.45))
  for (support in names(cases)) {
    condition <- cases[[support]]
    input <- results(stage, list(Y = model(fixed, condition$control, condition$disease)))
    output <- run_filter(stage, input, 0.8)
    check_audit(output, 0.8, TRUE, 1L)
    stopifnot(attr(output, "filter_summary")$target_supported_in == support)
  }
  overall_only <- results(stage, list(Y = model(fixed + 2, fixed, fixed)))
  output <- run_filter(stage, overall_only, 0.8)
  check_audit(output, 0.8, TRUE, 0L)

  for (invalid in list(0.9, 0.95, 0, 1, NA_real_, c(0.8, 0.81))) {
    expect_error(run_filter(stage, weak, invalid))
  }

  # Invalid successful-model intervals cannot silently join a valid fallback.
  malformed <- model()
  malformed$target_tf$screening_interval_lower <- NA_real_
  malformed$target_tf$screening_interval_excludes_zero <- FALSE
  invalid_only <- run_filter(stage, results(stage, list(INVALID = malformed)), 0.8)
  check_audit(invalid_only, 0.8, FALSE, 0L)
  expect_error(run_filter(stage, results(stage, list(Y = model(), INVALID = malformed)), 0.8),
               "invalid successful-model intervals for: INVALID")

  # New payloads must contain more than a point estimate and match diagnostics.
  one_draw <- model()
  one_draw$target_screening_draws$draw_count <- 1L
  one_draw$target_screening_draws$overall <- one_draw$target_screening_draws$overall[1L]
  one_draw$diagnostics$posterior_draws <- 1L
  expect_error(run_filter(stage, results(stage, list(Y = one_draw)), 0.8), "refit|re-fit")
  inconsistent <- model()
  inconsistent$diagnostics$posterior_draws <- length(fixed) - 1L
  expect_error(run_filter(stage, results(stage, list(Y = inconsistent)), 0.8), "refit|re-fit")
  missing_condition <- model(fixed, fixed + 1.45, fixed)
  missing_condition$target_screening_draws$control <- NULL
  expect_error(run_filter(stage, results(stage, list(Y = missing_condition)), 0.8), "refit|re-fit")
}

# MCMC screening cannot treat one completed chain as all requested chains.
for (legacy in c(FALSE, TRUE)) {
  complete <- model()
  if (legacy) complete$target_screening_draws <- NULL
  complete$sampling <- list(chains = 7L, iter_sampling = 143L)
  output <- run_filter("mcmc", results("mcmc", list(Y = complete)), 0.8)
  check_audit(output, 0.8, TRUE, 1L)
  incomplete <- complete
  incomplete$sampling <- list(chains = 3L, iter_sampling = 500L)
  expect_error(run_filter("mcmc", results("mcmc", list(Y = incomplete)), 0.8), "refit|re-fit")
  incomplete$sampling <- list(chains = 3L)
  expect_error(run_filter("mcmc", results("mcmc", list(Y = incomplete)), 0.8), "refit|re-fit")

  # Older models may carry sampling settings only on the result object.
  complete$sampling <- NULL
  inherited <- results("mcmc", list(Y = complete))
  attr(inherited, "sampling") <- list(chains = 7L, iter_sampling = 143L)
  output <- run_filter("mcmc", inherited, 0.8)
  check_audit(output, 0.8, TRUE, 1L)
  attr(inherited, "sampling") <- list(chains = 3L, iter_sampling = 500L)
  expect_error(run_filter("mcmc", inherited, 0.8), "refit|re-fit")
  inherited$Y$sampling <- list(chains = 7L, iter_sampling = 143L)
  output <- run_filter("mcmc", inherited, 0.8)
  check_audit(output, 0.8, TRUE, 1L)
}

# Restore legacy MCMC intervals only from demonstrably complete joint draws.
for (condition in c(FALSE, TRUE)) {
  legacy <- if (condition) model(fixed, fixed + 1.45, fixed) else model()
  legacy$target_screening_draws <- NULL
  restored <- run_filter("mcmc", results("mcmc", list(Y = legacy)), 0.8)
  check_audit(restored, 0.8, TRUE, 1L)
  stopifnot(restored$Y$target_tf$screening_interval_level == 0.8)
  incomplete <- legacy
  incomplete$stage2_parameter_draws$source_draw_index <- rev(seq_along(fixed))
  expect_error(run_filter("mcmc", results("mcmc", list(Y = incomplete)), 0.8), "refit|re-fit")
  if (condition) {
    incomplete <- legacy
    incomplete$stage2_parameter_draws$beta_target_control <- NULL
    expect_error(run_filter("mcmc", results("mcmc", list(Y = incomplete)), 0.8), "refit|re-fit")
  }
}

old_vi <- model()
old_vi$target_screening_draws <- NULL
old_vi$stage2_parameter_draws$draw_count <- 1L
old_vi$stage2_parameter_draws$source_draw_index <- 1L
old_vi$stage2_parameter_draws$beta <- old_vi$stage2_parameter_draws$beta[1L, , drop = FALSE]
old_vi$stage2_parameter_draws$beta_target_control <- fixed[1L]
old_vi$stage2_parameter_draws$beta_target_disease <- fixed[1L]
expect_error(run_filter("prescreen", results("prescreen", list(Y = old_vi)), 0.8), "refit|re-fit")
# Missing old draws matter only when a fallback is actually needed.
available <- results("prescreen", list(STRONG = model(fixed + 2), OLD = old_vi))
unchanged <- run_filter("prescreen", available, 0.8)
check_audit(unchanged, 0.8, FALSE, 1L, initial = 1L)

# A Stage 1 MCMC prescreen uses its full target draws for fallback, even
# though it does not retain a full Stage 2 hand-off.
prescreen_mcmc <- results("prescreen", list(Y = model()))
attr(prescreen_mcmc, "inference") <- "mcmc"
prescreen_mcmc$Y$sampling <- list(chains = 7L, iter_sampling = 143L)
prescreen_mcmc$Y$stage2_parameter_draws <- NULL
output <- run_filter("prescreen", prescreen_mcmc, 0.8)
check_audit(output, 0.8, TRUE, 1L)
attr(prescreen_mcmc, "input_pipeline_stage") <- "mcmc_input"
expect_error(run_filter("prescreen", prescreen_mcmc, 0.8), "Stage 1")

cat("PASS: explicit screening fallback, audit persistence, fixed confounder thresholds, condition OR support, failed-fit exclusion, complete-draw/chain checks and immutable source models.\n")
