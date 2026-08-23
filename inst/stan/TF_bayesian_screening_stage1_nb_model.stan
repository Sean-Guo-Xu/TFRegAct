data {
  int<lower=1> N;                         // cells
  int<lower=1> P;                         // target TF + gene-specific confounder TFs
  array[N] int<lower=0> Y;                // raw target-gene counts
  matrix[N, P] X;                         // normally standardized TF expression

  int<lower=1> K_batch;
  array[N] int<lower=1, upper=K_batch> batch;
  array[N] int<lower=0, upper=1> condition; // control = 0, disease = 1
  vector[N] log_offset;

  vector<lower=0>[P] confidence;          // original q; values are normally in [0, 10]
  array[P] int<lower=-1, upper=1> direction;
  real<lower=0> beta_prior_scale;
  real<lower=0> eta;
  real<lower=1> r_dir;                    // prior odds favoring the annotated sign
  real<lower=0> batch_prior_scale;

  int<lower=1, upper=P> target_tf_index;
  real<lower=0> target_interaction_sd;
}

parameters {
  real alpha;
  real condition_effect;
  vector[P] beta;                         // target entry is the condition-average beta
  real beta_target_delta;                 // disease minus control target-TF beta
  vector[K_batch] batch_raw;
  real log_phi;
}

transformed parameters {
  vector[K_batch] batch_effect;
  vector<lower=0>[P] beta_prior_sd;
  vector[P] beta_prior_mean;
  real<lower=0> phi;
  real beta_target_control;
  real beta_target_disease;

  // A centered batch effect also works for K_batch = 1, where it is exactly 0.
  batch_effect = batch_prior_scale * (batch_raw - mean(batch_raw));
  phi = exp(log_phi);

  // Under this Normal prior, r_dir/(1+r_dir) is the prior probability of the
  // annotated sign. Unknown directions remain zero-centered. Confidence sets
  // the amount of shrinkage: a higher q gives a wider prior.
  for (j in 1:P) {
    beta_prior_sd[j] =
      beta_prior_scale * pow(confidence[j] / 10.0, eta);
    if (direction[j] == 0) {
      beta_prior_mean[j] = 0;
    } else {
      beta_prior_mean[j] =
        direction[j]
        * inv_Phi(r_dir / (1.0 + r_dir))
        * beta_prior_sd[j];
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

  beta ~ normal(beta_prior_mean, beta_prior_sd);
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
