data {
  int<lower=1> N;
  int<lower=1> P;
  array[N] int<lower=0> Y;
  matrix[N, P] X;

  int<lower=0> Q;
  matrix[N, Q] W;
  vector<lower=0>[Q] W_prior_scale;
  int<lower=0, upper=1> has_condition;
  int<lower=0, upper=Q> condition_w_index;
  array[N] int<lower=0, upper=1> condition;

  int<lower=0> K_batch;
  array[N] int<lower=0, upper=K_batch> batch;
  matrix[K_batch, Q] batch_level_design;

  vector[N] log_offset;
  vector[P] beta_prior_mean;
  vector<lower=0>[P] beta_prior_sd;
  real<lower=0> alpha_prior_sd;

  int<lower=1, upper=P> target_tf_index;
  int<lower=0, upper=1> use_target_interaction;
  vector[N] target_condition_centered;
  real<lower=0> target_interaction_sd;
}

parameters {
  real alpha;
  vector[P] beta;
  vector[Q] zeta;
  vector[use_target_interaction] beta_target_delta_free;
  real log_phi;
}

transformed parameters {
  vector[N] log_mu;
  vector[K_batch] batch_effect;
  real condition_effect;
  real<lower=0> phi;
  real beta_target_mean;
  real beta_target_delta;
  real beta_target_control;
  real beta_target_disease;

  batch_effect = batch_level_design * zeta;
  condition_effect = 0;
  if (has_condition == 1) {
    condition_effect = zeta[condition_w_index];
  }
  phi = exp(log_phi);

  beta_target_mean = beta[target_tf_index];
  beta_target_delta = 0;
  if (use_target_interaction == 1) {
    beta_target_delta = beta_target_delta_free[1];
  }
  beta_target_control = beta_target_mean - 0.5 * beta_target_delta;
  beta_target_disease = beta_target_mean + 0.5 * beta_target_delta;

  log_mu = alpha + X * beta + W * zeta + log_offset;
  for (i in 1:N) {
    log_mu[i] += target_condition_centered[i]
      * X[i, target_tf_index] * beta_target_delta;
  }
}

model {
  alpha ~ normal(0, alpha_prior_sd);
  if (Q > 0) {
    zeta ~ normal(rep_vector(0, Q), W_prior_scale);
  }
  log_phi ~ normal(0, 1);
  beta ~ normal(beta_prior_mean, beta_prior_sd);
  if (use_target_interaction == 1) {
    beta_target_delta_free[1] ~ normal(0, target_interaction_sd);
  }

  Y ~ neg_binomial_2_log(log_mu, phi);
}

generated quantities {
  vector[N] log_lik;
  real target_tf_regulatory_total_delta;

  for (i in 1:N) {
    log_lik[i] = neg_binomial_2_log_lpmf(Y[i] | log_mu[i], phi);
  }

  target_tf_regulatory_total_delta = 0;
  if (has_condition == 1) {
    real xbar_control_target = 0;
    real xbar_disease_target = 0;
    int N_control = 0;
    int N_disease = 0;

    for (i in 1:N) {
      if (condition[i] == 0) {
        xbar_control_target += X[i, target_tf_index];
        N_control += 1;
      } else {
        xbar_disease_target += X[i, target_tf_index];
        N_disease += 1;
      }
    }
    if (N_control > 0) {
      xbar_control_target /= N_control;
    }
    if (N_disease > 0) {
      xbar_disease_target /= N_disease;
    }

    target_tf_regulatory_total_delta =
      xbar_disease_target * beta_target_disease
      - xbar_control_target * beta_target_control;
  }
}
