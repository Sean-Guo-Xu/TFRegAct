// Joint hierarchical model for one target TF and multiple target genes.
//
// Key design:
//   1. The target TF has a strictly positive latent activity z_i.
//   2. Non-negative Seurat-normalized target-TF expression anchors z_i through
//      a hurdle-Gamma measurement model. Zeros are explicit non-detections;
//      conditional on detection, E[x_target_i | z_i] = z_i.
//   3. Each target gene may have a different set of other TF covariates.
//      Existing gene-TF relations are stored as a sparse edge list.
//   4. Stage 1 first removes unsupported target and confounder edges. Each
//      retained Stage 2 edge receives a Normal prior derived from its Stage 1
//      posterior mean and uncertainty after conversion back to the original
//      Seurat-normalized predictor scale.
//   5. Condition changes latent target-TF activity, but does not interact
//      with any gene-TF coefficient.
//   6. Batch is included only in the target-gene outcome model, not in the
//      latent target-TF activity model.

data {
  int<lower=2> N;                           // number of cells
  int<lower=1> G;                           // number of target genes
  array[N, G] int<lower=0> Y;               // raw target-gene counts

  array[N] int<lower=0, upper=1> condition; // 0 = control, 1 = disease
  vector[N] log_offset;                     // log library-size offset

  // Non-negative target-TF expression used to anchor positive activity.
  vector<lower=0>[N] target_tf_expression;

  // Union of all non-target TF predictors. A common coefficient-prior scale
  // is easiest to interpret when predictors have comparable scales.
  int<lower=0> K_other;
  matrix[N, K_other] X_other;

  // Sparse gene-TF edge table. Every target gene must have exactly one
  // target-TF edge. Non-target edges may differ arbitrarily among genes.
  int<lower=G> E;
  array[E] int<lower=1, upper=G> edge_gene;
  array[E] int<lower=0, upper=K_other> edge_regulator;
  array[E] int<lower=0, upper=1> edge_is_target;
  array[G] int<lower=1, upper=E> target_edge;

  // Empirical-Bayes priors transferred from Stage 1. The R interface applies
  // the same conversion to target-gamma and non-target confounder edges.
  vector[E] edge_prior_mean;
  vector<lower=1e-12>[E] edge_prior_sd;

  // Batch index for the target-gene outcome model only. With one batch, use
  // K_batch = 1 and batch[i] = 1; centered batch effects are then exactly 0.
  int<lower=1> K_batch;
  array[N] int<lower=1, upper=K_batch> batch;
}

transformed data {
  int N_control = 0;
  int N_disease = 0;
  array[G] int target_edge_count = rep_array(0, G);

  for (i in 1:N) {
    if (condition[i] == 0) {
      N_control += 1;
    } else {
      N_disease += 1;
    }
  }

  if (N_control == 0 || N_disease == 0) {
    reject("Both control and disease cells are required.");
  }

  for (e in 1:E) {
    if (edge_is_target[e] == 1) {
      if (edge_regulator[e] != 0) {
        reject("A target-TF edge must have edge_regulator = 0.");
      }
      target_edge_count[edge_gene[e]] += 1;
    } else if (edge_regulator[e] == 0) {
      reject("A non-target edge must have edge_regulator >= 1.");
    }
  }

  for (g in 1:G) {
    if (target_edge_count[g] != 1) {
      reject("Each target gene must have exactly one target-TF edge.");
    }
    if (edge_is_target[target_edge[g]] != 1 ||
        edge_gene[target_edge[g]] != g) {
      reject("target_edge[g] must identify the target-TF edge for gene g.");
    }
  }
}

parameters {
  // Positive latent target-TF activity: z_i = exp(log_activity_i).
  real activity_intercept;
  real condition_log_activity_delta;
  real<lower=0> sigma_log_activity;
  vector[N] log_activity_raw;

  // Hurdle-Gamma measurement parameters. Detection must increase (or remain
  // flat) with activity, hence the non-negative slope constraint.
  real<lower=0> target_measurement_shape;
  real target_detection_intercept;
  real<lower=0> target_detection_slope;

  // Target-gene outcome model.
  vector[G] alpha;
  vector[G] condition_effect;
  vector<lower=0>[G] phi;

  matrix[G, K_batch] gene_batch_raw;
  vector<lower=0>[G] sigma_gene_batch;

  // One coefficient for each retained target or non-target TF -> gene edge.
  vector[E] edge_coef;
}

transformed parameters {
  matrix[G, K_batch] gene_batch_effect;

  vector[N] log_activity;
  vector<lower=0>[N] activity;
  real<lower=0> mean_activity_control;
  real<lower=0> mean_activity_disease;

  matrix[N, G] log_mu;

  // Sum-to-zero gene-specific batch effects.
  for (g in 1:G) {
    real batch_mean = mean(gene_batch_raw[g]);

    for (b in 1:K_batch) {
      gene_batch_effect[g, b] =
        sigma_gene_batch[g] * (gene_batch_raw[g, b] - batch_mean);
    }
  }

  // The observed target-TF expression fixes the location and scale of z_i.
  for (i in 1:N) {
    log_activity[i] =
      activity_intercept
      + condition_log_activity_delta * condition[i]
      + sigma_log_activity * log_activity_raw[i];
    activity[i] = exp(log_activity[i]);
  }

  mean_activity_control = 0;
  mean_activity_disease = 0;
  for (i in 1:N) {
    if (condition[i] == 0) {
      mean_activity_control += activity[i];
    } else {
      mean_activity_disease += activity[i];
    }
  }
  mean_activity_control /= N_control;
  mean_activity_disease /= N_disease;

  // Gene-specific intercepts, condition main effects and grouping effects.
  for (i in 1:N) {
    for (g in 1:G) {
      log_mu[i, g] =
        alpha[g]
        + log_offset[i]
        + condition_effect[g] * condition[i]
        + gene_batch_effect[g, batch[i]];
    }
  }

  // Add each sparse gene-TF edge. Centering the target activity changes only
  // the gene intercept, not gamma_g or its prior scale.
  for (e in 1:E) {
    int g = edge_gene[e];

    if (edge_is_target[e] == 1) {
      for (i in 1:N) {
        log_mu[i, g] +=
          edge_coef[e] * (activity[i] - mean_activity_control);
      }
    } else {
      int k = edge_regulator[e];
      for (i in 1:N) {
        log_mu[i, g] += edge_coef[e] * X_other[i, k];
      }
    }
  }
}

model {
  // Latent target-TF activity.
  activity_intercept ~ normal(0, 2);
  condition_log_activity_delta ~ normal(0, 1);
  sigma_log_activity ~ normal(0, 1);
  log_activity_raw ~ std_normal();

  // Hurdle-Gamma anchor. Exact zeros are non-detections; positive observations
  // follow a Gamma(shape, rate) model with conditional mean activity[i].
  target_measurement_shape ~ lognormal(0, 1);
  target_detection_intercept ~ normal(0, 2);
  target_detection_slope ~ normal(0, 1);
  for (i in 1:N) {
    real detection_logit =
      target_detection_intercept + target_detection_slope * log_activity[i];
    if (target_tf_expression[i] > 0) {
      target += bernoulli_logit_lpmf(1 | detection_logit);
      target += gamma_lpdf(
        target_tf_expression[i] |
        target_measurement_shape,
        target_measurement_shape / activity[i]
      );
    } else {
      target += bernoulli_logit_lpmf(0 | detection_logit);
    }
  }

  // Gene-level nuisance parameters.
  alpha ~ normal(0, 2);
  condition_effect ~ normal(0, 1);
  phi ~ lognormal(0, 1);

  to_vector(gene_batch_raw) ~ std_normal();
  sigma_gene_batch ~ normal(0, 1);

  // Stage 1-informed Normal priors. Direction and confidence have already
  // influenced the Stage 1 posterior and are not applied a second time here.
  edge_coef ~ normal(edge_prior_mean, edge_prior_sd);

  // Multigene negative-binomial likelihood.
  for (i in 1:N) {
    for (g in 1:G) {
      Y[i, g] ~ neg_binomial_2_log(log_mu[i, g], phi[g]);
    }
  }
}

generated quantities {
  vector[G] target_gamma;
  vector[G] target_contribution_delta;

  real expected_activity_control =
    exp(activity_intercept + 0.5 * square(sigma_log_activity));
  real expected_activity_disease =
    exp(activity_intercept + condition_log_activity_delta
        + 0.5 * square(sigma_log_activity));
  real activity_ratio_disease_vs_control =
    exp(condition_log_activity_delta);
  real overall_rms_delta = 0;

  vector[N] log_lik_target_measurement;
  matrix[N, G] log_lik_gene;
  vector[N] log_lik_cell;
  vector[K_batch] log_lik_batch = rep_vector(0, K_batch);

  for (g in 1:G) {
    int e = target_edge[g];

    target_gamma[g] = edge_coef[e];
    target_contribution_delta[g] =
      target_gamma[g]
      * (expected_activity_disease - expected_activity_control);

    overall_rms_delta += square(target_contribution_delta[g]);
  }

  overall_rms_delta = sqrt(overall_rms_delta / G);

  for (i in 1:N) {
    real detection_logit =
      target_detection_intercept + target_detection_slope * log_activity[i];
    if (target_tf_expression[i] > 0) {
      log_lik_target_measurement[i] =
        bernoulli_logit_lpmf(1 | detection_logit)
        + gamma_lpdf(
            target_tf_expression[i] |
            target_measurement_shape,
            target_measurement_shape / activity[i]
          );
    } else {
      log_lik_target_measurement[i] =
        bernoulli_logit_lpmf(0 | detection_logit);
    }
    log_lik_cell[i] = log_lik_target_measurement[i];

    for (g in 1:G) {
      log_lik_gene[i, g] = neg_binomial_2_log_lpmf(
        Y[i, g] | log_mu[i, g], phi[g]
      );
      log_lik_cell[i] += log_lik_gene[i, g];
    }

    log_lik_batch[batch[i]] += log_lik_cell[i];
  }
}
