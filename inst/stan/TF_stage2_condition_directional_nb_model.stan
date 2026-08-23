data {
  int<lower=1> N;                         // number of cells
  int<lower=1> P;                         // number of filtered TFs
  array[N] int<lower=0> Y;                // target gene counts

  matrix[N, P] X;                         // filtered TF expression matrix
  array[N] int<lower=1> batch;            // batch index
  int<lower=1> K_batch;                   // number of batches
  array[N] int<lower=0, upper=1> condition; // 0 = control, 1 = disease

  vector[N] log_offset;                   // log library size divided by mean library size

  vector<lower=0>[P] confidence;          // confidence score, typically in [0, 10]
  array[P] int<lower=-1, upper=1> direction;
  vector[P] stage1_beta_mean;             // posterior mean from stage 1
  vector<lower=0>[P] stage1_beta_sd;      // posterior sd from stage 1

  int<lower=1, upper=P> target_tf_index;  // target TF for regulatory contribution output
  real<lower=0> direction_effect;         // absolute prior mean shift for known direction
  real<lower=0> beta_sd_floor;            // minimum prior sd
  real<lower=0> stage1_sd_multiplier;     // inflation factor for stage 1 posterior sd
}

parameters {
  real alpha;
  real condition_effect;                  // disease-control target-gene baseline shift

  vector[P] beta;                         // shared TF effects across conditions

  vector[K_batch] batch_raw;
  real<lower=0> sigma_batch;

  real<lower=0> phi;                      // NB overdispersion
}

transformed parameters {
  vector[K_batch] batch_effect;
  vector[N] log_mu;
  vector[P] beta_prior_mean;
  vector<lower=0>[P] beta_prior_sd;

  // Sum-to-zero batch effects for identifiability.
  batch_effect = sigma_batch * (batch_raw - mean(batch_raw));

  for (j in 1:P) {
    beta_prior_mean[j] =
      stage1_beta_mean[j] + direction[j] * direction_effect;
    beta_prior_sd[j] =
      fmax(beta_sd_floor, stage1_sd_multiplier * stage1_beta_sd[j]);
  }

  for (i in 1:N) {
    log_mu[i] =
      alpha
      + batch_effect[batch[i]]
      + condition_effect * condition[i]
      + X[i] * beta
      + log_offset[i];
  }
}

model {
  alpha ~ normal(0, 1);

  batch_raw ~ normal(0, 1);
  sigma_batch ~ normal(0, 1);

  condition_effect ~ normal(0, 1);
  phi ~ lognormal(0, 1);

  // Shared TF effects inherit stage 1 information and known direction.
  beta ~ normal(beta_prior_mean, beta_prior_sd);

  Y ~ neg_binomial_2_log(log_mu, phi);
}

generated quantities {
  vector[N] log_lik;
  real target_tf_regulatory_total_delta;

  {
    real xbar_control_target;
    real xbar_disease_target;
    int N_control;
    int N_disease;

    xbar_control_target = 0;
    xbar_disease_target = 0;
    N_control = 0;
    N_disease = 0;

    for (i in 1:N) {
      log_lik[i] = neg_binomial_2_log_lpmf(Y[i] | log_mu[i], phi);

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
      xbar_disease_target * beta[target_tf_index]
      - xbar_control_target * beta[target_tf_index];
  }
}
