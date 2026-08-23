data {
  int<lower=1> N;                         // number of cells
  int<lower=1> P;                         // number of TFs
  array[N] int<lower=0> Y;                // target gene counts

  matrix[N, P] X;                         // TF expression matrix

  vector[N] log_offset;                   // log library size divided by mean library size

  vector<lower=0>[P] confidence;          // raw confidence score
  array[P] int<lower=-1, upper=1> direction;

  real<lower=0> gamma;                    // global shrinkage scale
  real<lower=0> eta;                      // confidence exponent
  real<lower=1> r_dir;                    // opposite-direction penalty multiplier
}

parameters {
  real alpha;

  vector[P] beta;

  real<lower=0> phi;                      // NB overdispersion
}

transformed parameters {
  vector[N] log_mu;

  for (i in 1:N) {
    log_mu[i] =
      alpha
      + X[i] * beta
      + log_offset[i];
  }
}

model {
  // Hyperpriors / nuisance priors
  alpha ~ normal(0, 1);
  phi ~ lognormal(0, 1);

  // Confidence-weighted directional asymmetric Laplace prior for beta
  for (j in 1:P) {
    real b_j;
    real lambda_fav;
    real lambda_opp;
    real lambda_pos;
    real lambda_neg;

    b_j = pow(fmax(confidence[j], 1.0) / 10.0, eta) * gamma;

    lambda_fav = 1.0 / b_j;

    if (direction[j] == 0) {
      lambda_opp = lambda_fav;
    } else {
      lambda_opp = r_dir * lambda_fav;
    }

    if (direction[j] == 1) {
      // activation prior: positive beta is favored
      lambda_pos = lambda_fav;
      lambda_neg = lambda_opp;
    } else if (direction[j] == -1) {
      // repression prior: negative beta is favored
      lambda_pos = lambda_opp;
      lambda_neg = lambda_fav;
    } else {
      // unknown direction: symmetric Laplace
      lambda_pos = lambda_fav;
      lambda_neg = lambda_fav;
    }

    // Two-sided asymmetric Laplace prior centered at zero:
    // p(beta) = lambda_pos * lambda_neg / (lambda_pos + lambda_neg)
    //           * exp(-lambda_pos * beta) for beta >= 0
    //           * exp( lambda_neg * beta) for beta < 0
    target += log(lambda_pos)
              + log(lambda_neg)
              - log(lambda_pos + lambda_neg);

    if (beta[j] >= 0) {
      target += -lambda_pos * beta[j];
    } else {
      target +=  lambda_neg * beta[j];
    }
  }

  // Negative-binomial likelihood
  Y ~ neg_binomial_2_log(log_mu, phi);
}

generated quantities {
  vector[N] log_lik;

  for (i in 1:N) {
    log_lik[i] = neg_binomial_2_log_lpmf(Y[i] | log_mu[i], phi);
  }
}
