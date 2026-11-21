library(TFRegAct)
ns <- asNamespace("TFRegAct")
ns$tf_assert_target_interaction_dependencies()
expect_error <- function(expr, pattern) {
  err <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(err, "error"), grepl(pattern, conditionMessage(err)))
}
set.seed(18)
expr <- matrix(runif(90), 3, dimnames = list(c("T", "V", "W"), paste0("c", 1:30)))
a <- ns$create_TF_analysis_object(expr, batch = NULL, sample = NULL, target = c("T", "Y"),
  confidence = c(0.2, 9, 5), direction = c("activation", "repression", "unknown"),
  libsize = rep(100, 30), Y_exp = rep(3L, 30))
original <- ns$tf_prepare_stan_data(a, gamma = 1.7, eta = 0.8, r_dir = 4,
  confidence_min = 2, confidence_max = 8)
draws <- matrix(rnorm(300), 100, dimnames = list(NULL, paste0("beta[", 1:3, "]")))
fake <- list(stan_data = original, feature_names = rownames(expr),
  fit = list(draws = function(...) draws))
filtered <- ns$create_TF_analysis_object(expr[c("W", "T"), , drop = FALSE],
  batch = NULL, sample = NULL, target = c("T", "Y"),
  confidence = c(5, 0.2), direction = c("unknown", "activation"),
  libsize = rep(100, 30), Y_exp = rep(3L, 30))
make <- function(fit = fake, ...) ns$tf_prepare_stage2_stan_data(filtered, fit, ...)
d <- make()
prior <- attr(d, "stage2_prior")
stopifnot(identical(d$confidence, c(5, 2)), identical(d$direction, c(0L, 1L)),
  identical(d$gamma, 1.7), identical(d$eta, 0.8), identical(d$r_dir, 4),
  identical(prior$prior_family, "normal"),
  identical(d$beta_prior_mean, prior$beta_prior_mean),
  identical(d$beta_prior_variance, prior$beta_prior_variance),
  identical(d$beta_prior_sd, sqrt(prior$beta_prior_variance)),
  isTRUE(all.equal(d$beta_init, unname(colMeans(draws)[c(3, 1)]))))
weight <- (c(5, 2) - 1) / 9
stopifnot(isTRUE(all.equal(prior$beta_prior_mean, c(0, 0.2 + 0.4 * weight[2]))),
  isTRUE(all.equal(prior$beta_prior_variance, 0.5 + 1.5 * weight)))

# Even a radically different Stage 1 posterior only changes initialization.
changed <- fake
changed$fit <- list(draws = function(...) draws * 100 + 7)
other <- make(changed)
stopifnot(identical(attr(other, "stage2_prior"), prior),
  !identical(other$beta_init, d$beta_init))
messages <- character()
legacy <- withCallingHandlers(make(direction_effect = 50, beta_sd_floor = 99,
  stage1_sd_multiplier = 0.01), warning = function(w) {
    messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
  })
stopifnot(length(messages) == 1L, grepl("deprecated and ignored", messages),
  identical(attr(legacy, "stage2_prior"), prior))
missing_prior <- fake
missing_prior$stan_data <- NULL
expect_error(make(missing_prior), "original Stage 1 prior")
expect_error(ns$tf_stage2_knowledge_prior(modifyList(original, list(confidence = c(1, NA, 3)))),
  "Invalid")

# Array directions from Stan data produce plain, TF-aligned Normal parameters.
array_prior <- ns$tf_stage2_knowledge_prior(original)
stopifnot(is.null(dim(array_prior$beta_prior_mean)),
  array_prior$beta_prior_mean[1] > 0, array_prior$beta_prior_mean[2] < 0,
  array_prior$beta_prior_mean[3] == 0)

# The endpoints are anchored globally, independent of retained TFs and Stage 1
# shrinkage settings. These checks distinguish variance from standard deviation.
endpoint_data <- list(P = 6L, confidence = c(0, 1, 4, 7, 10, 20),
                     direction = c(1L, -1L, 0L, 1L, -1L, 1L))
endpoint <- ns$tf_stage2_knowledge_prior(endpoint_data)
stopifnot(isTRUE(all.equal(endpoint$beta_prior_mean, c(.2, -.2, 0, 7/15, -.6, .6))),
  isTRUE(all.equal(endpoint$beta_prior_variance, c(.5, .5, 1, 1.5, 2, 2))),
  isTRUE(all.equal(endpoint$beta_prior_sd ^ 2, endpoint$beta_prior_variance)),
  identical(endpoint, ns$tf_stage2_knowledge_prior(modifyList(endpoint_data,
    list(gamma = 99, eta = 9, r_dir = 50)))))
stage1_prior <- ns$tf_stage1_knowledge_prior(original)
scale <- 1.7 * (c(2, 8, 5) / 10)^0.8
stopifnot(stage1_prior$prior_family == "laplace",
  isTRUE(all.equal(stage1_prior$beta_lambda_positive, c(1, 4, 1) / scale)),
  isTRUE(all.equal(stage1_prior$beta_lambda_negative, c(4, 1, 1) / scale)))

# Exercise the actual shared fitting boundary without MCMC: the model receives
# the knowledge parameters and Stage 1 means only as starting values.
observed <- NULL
model <- list(sample = function(data, init, ...) {
  observed <<- list(data = data, init = init())
  "sampled"
})
result <- ns$tf_fit_shared_stage2(d, model, chains = 2, parallel_chains = 2,
  iter_warmup = 10, iter_sampling = 10, seed = 2, refresh = 0,
  adapt_delta = 0.9, max_treedepth = 10)
stopifnot(identical(result, "sampled"),
  identical(observed$data$confidence, d$confidence),
  identical(observed$data$direction, d$direction),
  identical(observed$init$beta, d$beta_init),
  identical(observed$data$beta_prior_mean, d$beta_prior_mean),
  identical(observed$data$beta_prior_sd, sqrt(d$beta_prior_variance)))
tampered <- d
tampered$beta_prior_sd <- tampered$beta_prior_variance
expect_error(ns$tf_validate_shared_stage2_data(tampered), "Normal prior parameters")
legacy_data <- d
legacy_data$beta_prior_mean <- NULL
legacy_data$beta_prior_sd <- NULL
expect_error(ns$tf_validate_shared_stage2_data(legacy_data), "missing")
cat("PASS: Normal knowledge prior mapping, Stage 1 preservation, posterior independence and Stan handoff.\n")
