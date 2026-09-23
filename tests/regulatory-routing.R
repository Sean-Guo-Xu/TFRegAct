library(TFRegAct)
ns <- asNamespace('TFRegAct')
expect_error <- function(expr, pattern) {
  e <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(e, 'error'), grepl(pattern, conditionMessage(e)))
}
for (type in c('single_cell', 'bulk')) for (condition in c(FALSE, TRUE)) {
  for (requested in c('no_interaction', 'target_interaction')) {
    observed <- ns$tf_regulatory_stage2_route(type, if (condition) 'group' else NULL, requested)
    expected <- if (!condition) 'no_interaction' else if (type == 'single_cell') 'hierarchical' else requested
    stopifnot(identical(observed, expected))
  }
}
set.seed(11)
counts <- matrix(rpois(3 * 40, 5), 3, dimnames = list(c('T','V','Y'), paste0('c', 1:40)))
obj <- SeuratObject::CreateSeuratObject(counts)
obj$donor <- rep(paste0('s', 1:4), each = 10)
obj$group <- rep(c('Normal','AAA'), each = 20)
input <- create_TF_computation_input(obj, 'T', 'Y', data_type = 'single_cell',
  sample_column = 'donor', condition_column = 'group')
values <- ns$tf_computation_input_values(input)
stopifnot(values$data_type == 'single_cell', values$sample_column == 'donor')
bulk <- create_TF_computation_input(obj, 'T', 'Y', data_type = 'bulk', condition_column = NULL)
stopifnot(is.null(ns$tf_computation_input_values(bulk)$condition_column))
expect_error(create_TF_computation_input(obj, 'T', data_type = 'invalid'), 'arg')
expect_error(create_TF_computation_input(obj, 'T', sample_column = 'missing'), 'not found')
x <- matrix(rnorm(80), 2, dimnames = list(c('T','V'), colnames(obj)))
a <- ns$create_TF_analysis_object(x, batch = NULL, sample = obj$group,
  target = c('T','Y'), libsize = colSums(counts), confidence = c(4,5),
  direction = c('activation','unknown'), Y_exp = counts['Y',])
draws <- matrix(rnorm(200), 100, dimnames = list(NULL,c('beta[1]','beta[2]')))
fake <- list(fit = list(draws = function(...) draws), feature_names = c('T','V'))
d <- ns$tf_prepare_stage2_stan_data(a, fake, 0.2, 0.5, 1.5, 'Normal','AAA',1,10,
  target_interaction = TRUE)
design <- ns$tf_hierarchical_sample_design(d, obj[[]], 'donor')
shuffled <- ns$tf_hierarchical_sample_design(d, obj[[]][40:1,], 'donor')
stopifnot(identical(design, shuffled), design$S == 4L,
  identical(d$X, t(x)), isTRUE(all.equal(d$beta_prior_mean[1], mean(draws[,1]) + 0.2)),
  isTRUE(all.equal(d$beta_prior_sd, pmax(0.5, 1.5 * apply(draws, 2, sd)), check.attributes = FALSE)))
expect_error(ns$tf_hierarchical_sample_design(d, obj[[]], NULL), 'requires')
bad <- obj[[]]; bad$donor[21] <- 's1'
expect_error(ns$tf_hierarchical_sample_design(d, bad, 'donor'), 'multiple conditions')
bad <- obj[[]]; bad$donor[bad$group == 'Normal'] <- 's1'
expect_error(ns$tf_hierarchical_sample_design(d, bad, 'donor'), 'two biological samples')
bad_d <- d; bad_d$X[1:10, 1] <- 0
expect_error(ns$tf_hierarchical_sample_design(bad_d, obj[[]], 'donor'), 'variation')
expect_error(ns$tf_hierarchical_sample_design(d, obj[[]][-1,], 'donor'), 'align')
# Exercise public dispatch without resampling: real object construction and
# preflight, with expensive graph/model operations replaced in a private scope.
scope <- new.env(parent = ns)
scope$find_adjustment_dagitty <- function(...) list()
scope$build_TF_analysis_object_from_seurat <- function(..., sample = NULL) {
  value <- a
  if (is.null(sample)) value$sample <- NULL
  value
}
scope$run_TF_directional_model <- function(...) fake
scope$filter_TF_analysis_object_by_beta_ci <- function(..., analysis_object) analysis_object
for (route in c('hierarchical', 'no_interaction', 'target_interaction')) {
  name <- switch(route, hierarchical = 'run_TF_stage2_hierarchical_model',
    no_interaction = 'run_TF_stage2_directional_model',
    target_interaction = 'run_TF_stage2_target_interaction_model')
  scope[[name]] <- local({ selected <- route; function(...) {
    scope$last_route <- selected
    fake
  }})
}
pipeline <- TF_regulatory_direction_computation
environment(pipeline) <- scope
for (type in c('single_cell','bulk')) for (condition in c(FALSE,TRUE)) {
  for (requested in c('no_interaction','target_interaction')) {
    input_case <- create_TF_computation_input(obj, 'T','Y', data_type = type,
      sample_column = if (type == 'single_cell') 'donor' else NULL,
      condition_column = if (condition) 'group' else NULL)
    result <- pipeline(input = input_case, stage2_model = requested)
    expected <- if (!condition) 'no_interaction' else if (type == 'single_cell') 'hierarchical' else requested
    stopifnot(scope$last_route == expected, result$stage2_model == expected,
      result$direction_summary$data_type == type)
  }
}
expect_error(pipeline(input = input, data_type = 'bulk'), 'disagrees')
cat('PASS: all routes, explicit input types, sample alignment, replicate and variation guards, Stage 2 prior transfer.\n')
