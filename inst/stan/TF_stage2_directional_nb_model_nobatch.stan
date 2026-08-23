data {
  int<lower=1> N;                         // number of cells
  int<lower=1> P;                         // number of filtered TFs
  array[N] int<lower=0> Y;                // target gene counts

  matrix[N, P] X;                         // filtered TF expression matrix

  vector[N] log_offset;                   // log library size divided by mean library size

  vector<lower=0>[P] confidence;          // confidence score, typically in [0, 10]
  array[P] int<lower=-1, upper=1> direction;
  vector[P] stage1_beta_mean;             // posterior mean from stage 1
  vector<lower=0>[P] stage1_beta_sd;      // posterior sd from stage 1

  real<lower=0> direction_effect;         // absolute prior mean shift for known direction
  real<lower=0> beta_sd_floor;            // minimum prior sd
  real<lower=0> stage1_sd_multiplier;     // inflation factor for stage 1 posterior sd
}

parameters {
  real alpha;

  vector[P] beta;

  real<lower=0> phi;                      // NB overdispersion
}

transformed parameters {
  vector[N] log_mu;
  vector[P] beta_prior_mean;
  vector<lower=0>[P] beta_prior_sd;

  for (j in 1:P) {
    beta_prior_mean[j] = stage1_beta_mean[j] + direction[j] * direction_effect;
    beta_prior_sd[j] = fmax(beta_sd_floor, stage1_sd_multiplier * stage1_beta_sd[j]);
  }

  for (i in 1:N) {
    log_mu[i] =
      alpha
      + X[i] * beta
      + log_offset[i];
  }
}

model {
  alpha ~ normal(0, 1);

  phi ~ lognormal(0, 1);

  // Stage 2 uses a regular Gaussian prior centered at the stage 1
  // posterior mean, with a mild direction-aware shift.
  beta ~ normal(beta_prior_mean, beta_prior_sd);

  Y ~ neg_binomial_2_log(log_mu, phi);
}

generated quantities {
  vector[N] log_lik;

  for (i in 1:N) {
    log_lik[i] = neg_binomial_2_log_lpmf(Y[i] | log_mu[i], phi);
  }
}
