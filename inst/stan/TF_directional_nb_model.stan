data {
  int<lower=1> N;
  int<lower=1> P;
  array[N] int<lower=0> Y;
  matrix[N, P] X;

  int<lower=0> Q;
  matrix[N, Q] W;
  vector<lower=0>[Q] W_prior_scale;
  int<lower=0, upper=1> has_condition;
  array[N] int<lower=0, upper=1> condition;
  int<lower=0, upper=Q> condition_w_index;

  int<lower=0> K_batch;
  array[N] int<lower=0, upper=K_batch> batch;
  matrix[K_batch, Q] batch_level_design;

  vector[N] log_offset;
  vector<lower=0>[P] confidence;
  array[P] int<lower=-1, upper=1> direction;
  real<lower=0> gamma;
  real<lower=0> eta;
  real<lower=1> r_dir;
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

  condition_effect = 0;
  if (has_condition == 1) {
    condition_effect = zeta[condition_w_index];
  }
  batch_effect = batch_level_design * zeta;
  phi = exp(log_phi);

  beta_target_mean = beta[target_tf_index];
  beta_target_delta = 0;
  if (use_target_interaction == 1) {
    beta_target_delta = beta_target_delta_free[1];
  }
  beta_target_control = beta_target_mean - 0.5 * beta_target_delta;
  beta_target_disease = beta_target_mean + 0.5 * beta_target_delta;

  log_mu = alpha + X * beta + W * zeta + log_offset;
  log_mu += target_condition_centered
    .* col(X, target_tf_index) * beta_target_delta;
}

model {
  alpha ~ normal(0, alpha_prior_sd);
  log_phi ~ normal(0, 1);
  if (Q > 0) {
    zeta ~ normal(rep_vector(0, Q), W_prior_scale);
  }
  if (use_target_interaction == 1) {
    beta_target_delta_free[1] ~ normal(0, target_interaction_sd);
  }

  // Zero-mode, confidence-weighted directional asymmetric Laplace prior.
  for (j in 1:P) {
    real b_j;
    real lambda_fav;
    real lambda_opp;
    real lambda_pos;
    real lambda_neg;

    b_j = pow(fmax(confidence[j], 1.0) / 10.0, eta) * gamma;
    lambda_fav = 1.0 / b_j;
    lambda_opp = direction[j] == 0 ? lambda_fav : r_dir * lambda_fav;

    if (direction[j] == 1) {
      lambda_pos = lambda_fav;
      lambda_neg = lambda_opp;
    } else if (direction[j] == -1) {
      lambda_pos = lambda_opp;
      lambda_neg = lambda_fav;
    } else {
      lambda_pos = lambda_fav;
      lambda_neg = lambda_fav;
    }

    target += log(lambda_pos)
              + log(lambda_neg)
              - log(lambda_pos + lambda_neg);
    if (beta[j] >= 0) {
      target += -lambda_pos * beta[j];
    } else {
      target += lambda_neg * beta[j];
    }
  }

  Y ~ neg_binomial_2_log(log_mu, phi);
}

generated quantities {
  vector[N] log_lik;

  for (i in 1:N) {
    log_lik[i] = neg_binomial_2_log_lpmf(Y[i] | log_mu[i], phi);
  }
}
