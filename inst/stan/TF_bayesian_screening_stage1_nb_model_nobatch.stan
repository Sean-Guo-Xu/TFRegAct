data {
  int<lower=1> N;
  int<lower=1> P;
  array[N] int<lower=0> Y;
  matrix[N, P] X;
  array[N] int<lower=0, upper=1> condition;
  vector[N] log_offset;
  vector<lower=0>[P] confidence;
  array[P] int<lower=-1, upper=1> direction;
  real<lower=0> beta_prior_scale;
  real<lower=0> eta;
  real<lower=1> r_dir;
  int<lower=1, upper=P> target_tf_index;
  real<lower=0> target_interaction_sd;
}

parameters {
  real alpha;
  real condition_effect;
  vector[P] beta;
  real beta_target_delta;
  real log_phi;
}

transformed parameters {
  vector<lower=0>[P] beta_prior_sd;
  vector[P] beta_prior_mean;
  real<lower=0> phi = exp(log_phi);
  real beta_target_control = beta[target_tf_index] - 0.5 * beta_target_delta;
  real beta_target_disease = beta[target_tf_index] + 0.5 * beta_target_delta;

  for (j in 1:P) {
    beta_prior_sd[j] = beta_prior_scale * pow(confidence[j] / 10.0, eta);
    beta_prior_mean[j] = direction[j] == 0
      ? 0
      : direction[j] * inv_Phi(r_dir / (1.0 + r_dir)) * beta_prior_sd[j];
  }
}

model {
  vector[N] log_mu = alpha + X * beta + log_offset;
  alpha ~ normal(0, 2);
  condition_effect ~ normal(0, 1);
  log_phi ~ normal(0, 1);
  beta ~ normal(beta_prior_mean, beta_prior_sd);
  beta_target_delta ~ normal(0, target_interaction_sd);
  for (i in 1:N) {
    log_mu[i] += condition_effect * condition[i]
      + (condition[i] - 0.5) * beta_target_delta * X[i, target_tf_index];
  }
  Y ~ neg_binomial_2_log(log_mu, phi);
}
