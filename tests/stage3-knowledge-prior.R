# Run directly with Rscript from the source package or tests directory.
# Rcpp is needed only for the final executable Gibbs-prior check.
package_root <- if (file.exists("R/activity_stage3_em.R")) "." else ".."
source(file.path(package_root, "R", "activity_stage3_em.R"), local = TRUE)

common <- list(N = 4L, G = 1L, Y = matrix(c(0, 2, 4, 1), 4, 1),
  has_condition = TRUE, condition_shift = c(-0.5, -0.5, 0.5, 0.5),
  target_center = 0.3, target_scale = 1.1,
  beta_prior_mean = 0.3, beta_prior_variance = 0.8, beta_prior_sd = sqrt(0.8),
  beta_delta_prior_mean = -0.15, beta_delta_prior_variance = 0.64, beta_delta_prior_sd = 0.8,
  beta_delta_prior_covariance = 0.2, target_genes = "Y")
prior_covariance <- function(x) matrix(c(x$beta_prior_variance,
  x$beta_delta_prior_covariance, x$beta_delta_prior_covariance,
  x$beta_delta_prior_variance), 2L, 2L)
prior_quadratic <- function(p, x) {
  difference <- p - c(x$beta_prior_mean, x$beta_delta_prior_mean)
  -0.5 * as.numeric(crossprod(difference, solve(prior_covariance(x), difference)))
}
estep <- list(nodes = matrix(c(0, 0.4, 1.2, 2, 0.7, 1, 1.8, 3), 4, 2),
  weights = matrix(rep(c(0.3, 0.7), each = 4), 4, 2))
eta0 <- c(-0.1, 0.1, 0.2, -0.2)
objective <- function(p, information = FALSE, input = common) {
  tf_em_gene_expected_log_posterior(p, 1L, input, eta0, 5, estep,
    return_information = information)
}
# Include beta=0: the Normal objective and gradient must be smooth there.
for (beta in c(-0.6, 0, 0.4)) {
  p <- c(beta, 0.2)
  numerical <- vapply(1:2, function(j) {
    plus <- minus <- p
    plus[j] <- plus[j] + 1e-6
    minus[j] <- minus[j] - 1e-6
    (objective(plus)$value - objective(minus)$value) / 2e-6
  }, numeric(1))
  numerical_information <- -vapply(1:2, function(j) {
    plus <- minus <- p
    plus[j] <- plus[j] + 1e-5
    minus[j] <- minus[j] - 1e-5
    (objective(plus)$gradient - objective(minus)$gradient) / 2e-5
  }, numeric(2))
  stopifnot(max(abs(objective(p)$gradient - numerical)) < 1e-6,
    max(abs(objective(p, TRUE)$information - numerical_information)) < 1e-6,
    !objective(p, TRUE)$at_prior_kink)
  other <- common
  other$beta_prior_mean <- -0.2
  other$beta_prior_variance <- 1.3
  other$beta_prior_sd <- sqrt(1.3)
  other$beta_delta_prior_mean <- 0.25
  other$beta_delta_prior_variance <- 0.9
  other$beta_delta_prior_sd <- sqrt(0.9)
  other$beta_delta_prior_covariance <- -0.1
  expected <- prior_quadratic(p, common) - prior_quadratic(p, other)
  info_difference <- objective(p, TRUE)$information - objective(p, TRUE, other)$information
  stopifnot(abs(objective(p)$value - objective(p, input = other)$value - expected) < 1e-10,
    max(abs(info_difference - (solve(prior_covariance(common)) - solve(prior_covariance(other))))) < 1e-10)
}
expected_prior <- prior_quadratic(c(-0.6, 0.2), common)
stopifnot(abs(tf_em_beta_log_prior(-0.6, 0.2, common) - expected_prior) < 1e-12)

# Compare unconstrained NB+Normal MAPs against scalar score-equation roots.
small <- list(N = 2L, G = 1L, Y = matrix(1, 2, 1), has_condition = FALSE,
  condition_shift = c(0, 0), target_center = 0, target_scale = 1,
  beta_prior_mean = 0, beta_prior_variance = 0.7, beta_prior_sd = sqrt(0.7),
  beta_delta_prior_mean = 0, beta_delta_prior_variance = 0, beta_delta_prior_sd = 0,
  beta_delta_prior_covariance = 0, target_genes = "Y")
fixed_estep <- list(nodes = matrix(1, 2, 1), weights = matrix(1, 2, 1))
settings <- list(mstep_maxit = 300L, mstep_reltol = 1e-12)
for (y in c(0, 1, 10)) {
  x <- small
  x$Y[,] <- y
  fitted <- tf_em_mstep(x, matrix(0, 2, 1), 5, fixed_estep, -0.2, 0, settings)
  score <- function(beta) 2 * (y - (y + 5) * plogis(beta - log(5))) - beta / 0.7
  reference <- uniroot(score, c(-50, 50), tol = 1e-12)$root
  stopifnot(abs(fitted$beta_mean - reference) < 1e-6,
    !fitted$at_prior_kink, all(is.finite(fitted$covariance)))
}

# With zero standardized activity the likelihood contains no beta/delta
# information: the exact Normal posterior moments recover the original priors.
prior_only <- small
prior_only$has_condition <- TRUE
prior_only$condition_shift <- c(-0.5, 0.5)
prior_only$beta_prior_mean <- -0.35
prior_only$beta_prior_variance <- 0.5
prior_only$beta_prior_sd <- sqrt(0.5)
prior_only$beta_delta_prior_mean <- 0.25
prior_only$beta_delta_prior_variance <- 0.64
prior_only$beta_delta_prior_sd <- 0.8
prior_only$beta_delta_prior_covariance <- -0.18
no_signal <- list(nodes = matrix(0, 2, 1), weights = matrix(1, 2, 1))
prior_fit <- tf_em_mstep(prior_only, matrix(0, 2, 1), 5, no_signal, 1, -1, settings)
stopifnot(abs(prior_fit$beta_mean + 0.35) < 1e-6,
  abs(prior_fit$beta_delta - 0.25) < 1e-6,
  max(abs(prior_fit$covariance[, , 1] - prior_covariance(prior_only))) < 1e-10,
  abs(tf_em_scaled_beta_change(1, 0, 0, 0, small) - 1/sqrt(0.7)) < 1e-12)

# v6 supplies the Stage 2 posterior Normal approximation without inflation.
# Obsolete network-prior handoffs and variance/SD confusion are rejected.
input <- list(interface_version = "tf_em_stage2_input_v6", N = 2L, G = 1L,
  S = 2L, Y = matrix(c(1, 2), 2, 1), eta0_draws = array(0, c(2, 2, 1)),
  phi_draws = matrix(5, 2, 1), has_condition = FALSE, condition = c(0L, 0L),
  target_tf_expression = c(0, 1), target_tf_center = 0.5, target_tf_scale = 1,
  beta_target_mean_init = 0, beta_target_delta_init = 0,
  beta_prior_family = "normal", beta_prior_mean = 0.2,
  beta_prior_source = "stage2_posterior_normal_approximation",
  beta_prior_variance = 0.5, beta_prior_sd = sqrt(0.5),
  beta_delta_prior_mean = 0, beta_delta_prior_variance = 0, beta_delta_prior_sd = 0,
  beta_delta_prior_covariance = 0,
  cell_names = c("c1", "c2"), target_genes = "Y", target_tf = "T")
stopifnot(identical(tf_em_validate_input(input), input))
for (change in list(list(interface_version = "tf_em_stage2_input_v5"),
                    list(beta_prior_family = "laplace"),
                    list(beta_prior_source = "network_prior"), list(beta_delta_prior_mean = 0.1),
                    list(beta_prior_variance = 0), list(beta_prior_sd = 0.5))) {
  stopifnot(inherits(try(tf_em_validate_input(modifyList(input, change)), silent = TRUE), "try-error"))
}
tiny <- modifyList(input, list(beta_prior_variance = 1e-12, beta_prior_sd = 1e-6))
stopifnot(identical(tf_em_validate_input(tiny), tiny),
  tf_em_gene_prior_information(tiny, 1L)[1L, 1L] == 1e12)
tiny_mismatch <- tiny
tiny_mismatch$beta_prior_sd <- 1e-5
stopifnot(inherits(try(tf_em_validate_input(tiny_mismatch), silent = TRUE), "try-error"))
condition_input <- modifyList(input, list(has_condition = TRUE, condition = 0:1,
  beta_target_delta_init = 0.15, beta_delta_prior_mean = 0.15,
  beta_delta_prior_variance = 0.25, beta_delta_prior_sd = 0.5,
  beta_delta_prior_covariance = 0.1, control_level = "Control", disease_level = "Disease"))
stopifnot(identical(tf_em_validate_input(condition_input), condition_input))
for (change in list(list(beta_delta_prior_covariance = 1),
                    list(beta_delta_prior_covariance = sqrt(0.5) * 0.5),
                    list(beta_delta_prior_variance = 0), list(beta_delta_prior_sd = 0.6))) {
  stopifnot(inherits(try(tf_em_validate_input(modifyList(condition_input, change)), silent = TRUE), "try-error"))
}

output <- tempfile(fileext = ".rds")
checkpoint <- tempfile(fileext = ".rds")
run <- function(value, draws, resume = FALSE) run_TF_EM_latent_activity(
  value, algorithm = "quadrature_em", nuisance_draw_count = draws, cores = 1L,
  max_iter = 2L, min_iter = 2L, quadrature_nodes = 5L,
  output_file = output, checkpoint_file = checkpoint, resume = resume)
single <- run(input, 0L)
multiple <- run(input, 2L)
stopifnot(identical(single$interface_version, "tf_em_latent_activity_fit_v4"),
  single$settings$beta_prior_family == "normal",
  !any(multiple$beta_at_prior_kink_by_draw),
  all(is.finite(multiple$beta_mean_local_normal_sample_by_draw)),
  all(is.finite(multiple$beta_summary$beta_mean_sd)),
  all(is.finite(multiple$activity_summary$activity_mean)))
changed <- input
changed$beta_prior_mean <- -0.2
new_fit <- run(changed, 2L, resume = TRUE)
stopifnot(!identical(attr(new_fit, "em_signature"), attr(multiple, "em_signature")))
changed$beta_prior_variance <- 1
changed$beta_prior_sd <- 1
variance_fit <- run(changed, 2L, resume = TRUE)
stopifnot(!identical(attr(new_fit, "em_signature"), attr(variance_fit, "em_signature")))
condition_fit <- run(condition_input, 2L)
delta_changed <- condition_input
delta_changed$beta_delta_prior_mean <- -0.25
delta_fit <- run(delta_changed, 2L, resume = TRUE)
stopifnot(!identical(attr(condition_fit, "em_signature"), attr(delta_fit, "em_signature")))
delta_changed$beta_delta_prior_covariance <- -0.1
covariance_fit <- run(delta_changed, 2L, resume = TRUE)
stopifnot(!identical(attr(delta_fit, "em_signature"), attr(covariance_fit, "em_signature")))
unlink(c(output, checkpoint))

# Load only a matching installed implementation; otherwise compile the current
# C++ source into a temporary cache so the check cannot exercise an old DLL.
matching_installed <- requireNamespace("TFRegAct", quietly = TRUE) &&
  identical(body(get("tf_em_beta_log_prior", asNamespace("TFRegAct"))), body(tf_em_beta_log_prior))
if (matching_installed) {
  tf_activity_gibbs_ars_cpp <- get("tf_activity_gibbs_ars_cpp", asNamespace("TFRegAct"))
} else {
  if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required for the Gibbs prior check")
  Rcpp::sourceCpp(file.path(package_root, "src", "activity_gibbs_ars.cpp"),
                  env = environment(), cacheDir = tempfile("tfregact-normal-cpp-"),
                  rebuild = TRUE, showOutput = FALSE)
}
# An effectively uninformative predictor makes the scalar beta conditional
# approach its Normal prior; this checks the executable sampler's mean and SD.
prior_common <- list(N = 2L, G = 1L, Y = matrix(0, 2, 1), has_condition = FALSE,
  condition_shift = c(0, 0), activity_anchor = c(0, 1), active_prior_probability = c(0.2, 0.9),
  sigma_activity = 1, target_center = 0, target_scale = 1e10,
  beta_mean_init = 0, beta_delta_init = 0, beta_prior_mean = 0.4,
  beta_prior_variance = 1.2, beta_prior_sd = sqrt(1.2),
  beta_delta_prior_mean = 0, beta_delta_prior_variance = 0, beta_delta_prior_sd = 0,
  beta_delta_prior_covariance = 0)
set.seed(20261120)
samples <- tf_activity_gibbs_ars_cpp(prior_common, matrix(0, 2, 1), 5,
  list(iter_warmup = 100L, iter_sampling = 1500L, scale_move_sd = 1e-4))
stopifnot(abs(mean(samples$beta_samples) - 0.4) < 0.12,
  abs(var(as.numeric(samples$beta_samples)) - 1.2) < 0.2,
  all(is.finite(samples$log_joint_trace)))
joint_common <- modifyList(prior_common, list(has_condition = TRUE, condition_shift = c(-0.5, 0.5),
  beta_delta_prior_mean = -0.3, beta_delta_prior_variance = 0.7,
  beta_delta_prior_sd = sqrt(0.7), beta_delta_prior_covariance = 0.45))
set.seed(20261120)
joint_samples <- tf_activity_gibbs_ars_cpp(joint_common, matrix(0, 2, 1), 5,
  list(iter_warmup = 150L, iter_sampling = 2500L, scale_move_sd = 1e-4))
sample_covariance <- cov(cbind(as.numeric(joint_samples$beta_samples), as.numeric(joint_samples$delta_samples)))
stopifnot(abs(mean(joint_samples$beta_samples) - 0.4) < 0.12,
  abs(mean(joint_samples$delta_samples) + 0.3) < 0.1,
  max(abs(sample_covariance - prior_covariance(joint_common))) < 0.18,
  abs(sample_covariance[1, 2] - 0.45) < 0.1,
  all(is.finite(joint_samples$log_joint_trace)))

output <- tempfile(fileext = ".rds")
checkpoint <- tempfile(fileext = ".rds")
gibbs <- run_TF_EM_latent_activity(input, algorithm = "gibbs_ars", cores = 1L,
  iter_warmup = 25L, iter_sampling = 60L, resume = FALSE,
  output_file = output, checkpoint_file = checkpoint)
stopifnot(gibbs$interface_version == "tf_em_latent_activity_fit_v5",
  gibbs$settings$beta_prior_family == "normal",
  all(is.finite(gibbs$beta_summary$beta_mean_mean)),
  all(is.finite(gibbs$activity_summary$activity_mean)))
unlink(c(output, checkpoint))
message("Stage3 joint posterior-Normal gradient/curvature, MAP, v6/cache and executable Gibbs checks passed.")
