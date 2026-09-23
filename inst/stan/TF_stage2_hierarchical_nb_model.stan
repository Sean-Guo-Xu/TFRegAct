// Stage 2 hierarchical regulatory model for conditioned single-cell data.
// X: normalized predictors on ONE common scale across biological samples.
// Do not standardize separately within samples. Centering is done below.
// sample_id identifies biological samples, not sequencing batches.
// W: optional full-rank nuisance design (batch, condition main effect, etc.).
// No intercept or sample dummy columns in W. With condition present, include
// its main effect in W; the separate condition input defines slope groups.
// Suggested starting scales (caller must supply; Stan has no data defaults):
// sample_*_sd_scale=1; beta priors are constructed by the shared Stage 2 helper.
data {
  int<lower=1> N;
  int<lower=1> P;
  array[N] int<lower=0> Y;
  matrix[N, P] X;
  vector[N] log_offset;
  int<lower=2> S;
  array[N] int<lower=1, upper=S> sample_id;
  int<lower=0, upper=1> has_condition;
  array[S] int<lower=0, upper=1> sample_condition;
  int<lower=0, upper=1> use_target_interaction;
  int<lower=1, upper=P> target_tf_index;
  int<lower=0> Q;
  matrix[N, Q] W;
  vector<lower=0>[Q] W_prior_scale;
  vector[P] beta_prior_mean;
  vector<lower=0>[P] beta_prior_sd;
  real<lower=0> alpha_prior_sd;
  real<lower=0> target_interaction_sd;
  real<lower=0> sample_intercept_sd_scale;
  real<lower=0> sample_slope_sd_scale;
  int<lower=0, upper=1> save_log_lik;
}
transformed data {
  array[S] int sample_n = rep_array(0, S);
  matrix[S, P] X_sample_mean = rep_matrix(0, S, P);
  matrix[N, P] X_within;
  vector[S] condition_centered = rep_vector(0, S);
  int n_control = 0;
  int n_disease = 0;
  if (alpha_prior_sd <= 0 || target_interaction_sd <= 0
      || sample_intercept_sd_scale <= 0 || sample_slope_sd_scale <= 0)
    reject("All prior scale inputs must be strictly positive.");
  if (Q > 0) {
    for (q in 1:Q) {
      if (W_prior_scale[q] <= 0) reject("W_prior_scale must be positive.");
    }
  }
  if (use_target_interaction == 1 && has_condition == 0)
    reject("Target interaction requires condition.");
  for (i in 1:N) {
    sample_n[sample_id[i]] += 1;
    X_sample_mean[sample_id[i]] += X[i];
  }
  for (s in 1:S) {
    if (sample_n[s] == 0) reject("Every declared sample must contain cells.");
    X_sample_mean[s] /= sample_n[s];
    if (has_condition == 1) {
      condition_centered[s] = sample_condition[s] - 0.5;
      if (sample_condition[s] == 0) n_control += 1;
      else n_disease += 1;
    }
  }
  if (has_condition == 1 && (n_control == 0 || n_disease == 0))
    reject("Condition requires biological samples in both groups.");
  if (use_target_interaction == 1 && (n_control < 2 || n_disease < 2))
    reject("Slope comparison requires at least two samples per group.");
  for (i in 1:N) X_within[i] = X[i] - X_sample_mean[sample_id[i]];
}
parameters {
  real alpha;
  vector[P] beta; // Population within-sample slopes; target is group midpoint.
  vector[Q] zeta;
  vector[use_target_interaction] beta_target_delta_free;
  vector[S] sample_intercept_raw;
  vector[S] sample_slope_raw;
  real<lower=0> sample_intercept_sd;
  real<lower=0> sample_slope_sd;
  real log_phi;
}
transformed parameters {
  real<lower=0> phi = exp(log_phi);
  real beta_target_mean = beta[target_tf_index];
  real beta_target_delta = 0;
  real beta_target_control;
  real beta_target_disease;
  vector[S] sample_intercept;
  vector[S] beta_target_sample;
  vector[N] log_mu;
  if (use_target_interaction == 1)
    beta_target_delta = beta_target_delta_free[1];
  beta_target_control = beta_target_mean - 0.5 * beta_target_delta;
  beta_target_disease = beta_target_mean + 0.5 * beta_target_delta;
  sample_intercept = alpha + sample_intercept_sd * sample_intercept_raw;
  beta_target_sample = beta_target_mean
                      + condition_centered * beta_target_delta
                      + sample_slope_sd * sample_slope_raw;
  log_mu = X_within * beta + W * zeta + log_offset;
  for (i in 1:N) {
    int s = sample_id[i];
    log_mu[i] += sample_intercept[s]
                 + (beta_target_sample[s] - beta_target_mean)
                   * X_within[i, target_tf_index];
  }
}
model {
  alpha ~ normal(0, alpha_prior_sd);
  log_phi ~ normal(0, 1);
  sample_intercept_raw ~ std_normal();
  sample_slope_raw ~ std_normal();
  // Independent intercept/slope deviations; shared slope SD across conditions.
  sample_intercept_sd ~ normal(0, sample_intercept_sd_scale);
  sample_slope_sd ~ normal(0, sample_slope_sd_scale);
  if (Q > 0) zeta ~ normal(rep_vector(0, Q), W_prior_scale);
  if (use_target_interaction == 1)
    beta_target_delta_free[1] ~ normal(0, target_interaction_sd);
  // Stage 1 posterior-informed Normal priors on the same predictor scale.
  beta ~ normal(beta_prior_mean, beta_prior_sd);
  Y ~ neg_binomial_2_log(log_mu, phi);
}
generated quantities {
  vector[save_log_lik * N] log_lik;
  vector[save_log_lik * S] log_lik_sample;
  // Conditional on fitted sample effects. NOT a marginal leave-new-sample-out
  // likelihood: new-sample validation requires integrating/refitting effects.
  if (save_log_lik == 1) {
    log_lik_sample = rep_vector(0, S);
    for (i in 1:N) {
      log_lik[i] = neg_binomial_2_log_lpmf(Y[i] | log_mu[i], phi);
      log_lik_sample[sample_id[i]] += log_lik[i];
    }
  }
}
