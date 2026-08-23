data {
  int<lower=1> N;
  int<lower=1> P;
  array[N] int<lower=0> Y;
  matrix[N, P] X;

  int<lower=1> K_batch;
  array[N] int<lower=1, upper=K_batch> batch;
  array[N] int<lower=0, upper=1> condition;
  vector[N] log_offset;

  vector<lower=0>[P] confidence;
  array[P] int<lower=-1, upper=1> direction;
  real<lower=0> beta_prior_scale;
  real<lower=0> eta;
  real<lower=1> r_dir;
  real<lower=0> batch_prior_scale;

  int<lower=1, upper=P> target_tf_index;
  real<lower=0> target_interaction_sd;
}

parameters {
  real alpha;
  real condition_effect;
  vector[P] beta;
  real beta_target_delta;
  vector[K_batch] batch_raw;
  real log_phi;
}

transformed parameters {
  vector[K_batch] batch_effect;
  vector<lower=0>[P] beta_prior_laplace_scale;
  vector[P] beta_prior_mean;
  real<lower=0> phi;
  real beta_target_control;
  real beta_target_disease;

  batch_effect = batch_prior_scale * (batch_raw - mean(batch_raw));
  phi = exp(log_phi);

  // Confidence controls shrinkage exactly as in the Normal-prior model:
  // higher confidence gives a wider prior. For a known direction, the
  // Laplace location gives prior sign odds r_dir:1.
  for (j in 1:P) {
    beta_prior_laplace_scale[j] =
      beta_prior_scale * pow(confidence[j] / 10.0, eta);
    if (direction[j] == 0) {
      beta_prior_mean[j] = 0;
    } else {
      beta_prior_mean[j] =
        direction[j]
        * log((1.0 + r_dir) / 2.0)
        * beta_prior_laplace_scale[j];
    }
  }

  beta_target_control =
    beta[target_tf_index] - 0.5 * beta_target_delta;
  beta_target_disease =
    beta[target_tf_index] + 0.5 * beta_target_delta;
}

model {
  vector[N] log_mu;

  alpha ~ normal(0, 2);
  condition_effect ~ normal(0, 1);
  batch_raw ~ normal(0, 1);
  log_phi ~ normal(0, 1);

  beta ~ double_exponential(beta_prior_mean, beta_prior_laplace_scale);
  beta_target_delta ~ normal(0, target_interaction_sd);

  log_mu =
    alpha
    + X * beta
    + batch_effect[batch]
    + log_offset;
  for (i in 1:N) {
    log_mu[i] += condition_effect * condition[i]
      + (condition[i] - 0.5)
        * beta_target_delta * X[i, target_tf_index];
  }
  Y ~ neg_binomial_2_log(log_mu, phi);
}
