library(TFRegAct)
ns <- asNamespace('TFRegAct')

# The Stage 2 hand-off preserves every joint draw, including the last chain.
variables <- c('alpha', 'condition_effect', 'phi', 'beta[1]',
               'beta_target_control', 'beta_target_disease', 'beta_target_delta')
draws <- matrix(seq_len(1800 * length(variables)), 1800,
                dimnames = list(NULL, variables))
entry <- list(
  stan_data = list(P = 1L, K_batch = 0L, beta_prior_mean = 0,
                   beta_prior_sd = 1, target_interaction_sd = 1,
                   has_condition = 0L),
  feature_names = 'T', target_tf_index = 1L, target_gene = 'Y',
  predictor_edges = data.frame(tf = 'T', effect = 'activation',
                              direction = 1L, confidence = 5)
)
fake_fit <- list(
  draws = function(...) draws,
  summary = function(...) data.frame(variable = 'beta[1]', rhat = 1,
                                     ess_bulk = 1000, ess_tail = 1000)
)
handoff <- ns$tf_stage1_screening_summarize_fit(
  fake_fit, entry, 'mcmc', NULL, 'normal', 0.9, 0.9
)$stage2_parameter_draws
stopifnot(handoff$draw_count == 1800L,
          identical(handoff$source_draw_index, seq_len(1800)),
          identical(handoff$alpha, as.numeric(draws[, 'alpha'])))

# Exercise actual conditional EM runs on a small synthetic input.
input <- list(
  interface_version = 'tf_em_stage2_input_v3', N = 2L, G = 1L, S = 1800L,
  Y = matrix(c(1, 3), 2, 1),
  eta0_draws = array(rep(seq(-0.1, 0.1, length.out = 1800), 2), c(1800, 2, 1)),
  phi_draws = matrix(5, 1800, 1), has_condition = FALSE, condition = c(0L, 0L),
  target_tf_expression = c(0, 1), target_tf_center = 0, target_tf_scale = 1,
  beta_target_mean_init = 0.2, beta_target_delta_init = 0,
  beta_prior_mean = 0.2, beta_prior_sd = 1, target_interaction_sd = 1,
  cell_names = c('c1', 'c2'), target_genes = 'Y', target_tf = 'T'
)
run <- function(...) ns$run_TF_EM_latent_activity(
  input, cores = 1L, max_iter = 2L, quadrature_nodes = 5L,
  output_file = tempfile(fileext = '.rds'), resume = FALSE, ...
)
a <- run(seed = 9L)
b <- run(seed = 9L)
c <- run(seed = 10L)
stopifnot(a$nuisance_draw_count == 100L,
          length(unique(a$nuisance_draw_ids)) == 100L,
          max(a$nuisance_draw_ids) > 1500L,
          identical(a$nuisance_draw_ids, b$nuisance_draw_ids),
          !identical(a$nuisance_draw_ids, c$nuisance_draw_ids),
          length(a$draw_results) == 100L)
zero <- run(nuisance_draw_count = 0L)
stopifnot(zero$nuisance_mode == 'posterior_mean',
          identical(zero$nuisance_draw_ids, 0L), length(zero$draw_results) == 1L)
message('Full posterior retention, seeded sampling and zero-draw EM checks passed.')
