data {
  int<lower=1> N;
  int<lower=1> P;
  array[N] int<lower=0> Y;

  matrix[N, P] X;
  array[N] int<lower=1> batch;
  int<lower=1> K_batch;
  array[N] int<lower=0, upper=1> condition;

  vector[N] log_offset;

  vector<lower=0>[P] confidence;
  array[P] int<lower=-1, upper=1> direction;
  vector[P] stage1_beta_mean;
  vector<lower=0>[P] stage1_beta_sd;

  int<lower=1, upper=P> target_tf_index;
  real<lower=0> direction_effect;
  real<lower=0> beta_sd_floor;
  real<lower=0> stage1_sd_multiplier;
  real<lower=0> target_interaction_sd;
}

parameters {
  real alpha;
  real condition_effect;
  vector[P - 1] beta_nontarget;
  real beta_target_control;
  real beta_target_disease;

  vector[K_batch] batch_raw;
  real<lower=0> sigma_batch;
  real<lower=0> phi;
}

transformed parameters {
  vector[K_batch] batch_effect;
  vector[N] log_mu;
  vector[P] beta;
  vector[P] beta_prior_mean;
  vector<lower=0>[P] beta_prior_sd;
  real beta_target_mean;
  real beta_target_delta;

  batch_effect = sigma_batch * (batch_raw - mean(batch_raw));
  beta_target_mean = 0.5 * (beta_target_control + beta_target_disease);
  beta_target_delta = beta_target_disease - beta_target_control;

  for (j in 1:P) {
    beta_prior_mean[j] =
      stage1_beta_mean[j] + direction[j] * direction_effect;
    beta_prior_sd[j] =
      fmax(beta_sd_floor, stage1_sd_multiplier * stage1_beta_sd[j]);
  }

  {
    int beta_nontarget_index;

    beta_nontarget_index = 1;
    for (j in 1:P) {
      if (j == target_tf_index) {
        beta[j] = beta_target_mean;
      } else {
        beta[j] = beta_nontarget[beta_nontarget_index];
        beta_nontarget_index += 1;
      }
    }
  }

  for (i in 1:N) {
    real beta_target_condition;

    beta_target_condition = condition[i] == 1 ?
      beta_target_disease :
      beta_target_control;
    log_mu[i] =
      alpha
      + batch_effect[batch[i]]
      + condition_effect * condition[i]
      + X[i] * beta
      + X[i, target_tf_index] * (beta_target_condition - beta[target_tf_index])
      + log_offset[i];
  }
}

model {
  alpha ~ normal(0, 1);
  batch_raw ~ normal(0, 1);
  sigma_batch ~ normal(0, 1);
  condition_effect ~ normal(0, 1);
  phi ~ lognormal(0, 1);

  {
    int beta_nontarget_index;

    beta_nontarget_index = 1;
    for (j in 1:P) {
      if (j != target_tf_index) {
        beta_nontarget[beta_nontarget_index] ~ normal(
          beta_prior_mean[j],
          beta_prior_sd[j]
        );
        beta_nontarget_index += 1;
      }
    }
  }
  beta_target_mean ~ normal(
    beta_prior_mean[target_tf_index],
    beta_prior_sd[target_tf_index]
  );
  beta_target_delta ~ normal(0, target_interaction_sd);

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
      xbar_disease_target * beta_target_disease
      - xbar_control_target * beta_target_control;
  }
}
