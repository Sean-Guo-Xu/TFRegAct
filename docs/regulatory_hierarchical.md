# Regulatory Stage 2 hierarchy

## Routing and inputs

The regulatory workflow runs network adjustment search, the original Stage 1,
and predictor filtering before Stage 2. The input object's `data_type` selects
`single_cell` or `bulk`; it is not inferred from a Seurat container.
`sample_column` denotes biological sample IDs. `condition_column` denotes the
two comparison groups, and `batch_column` remains a separate nuisance term.

```text
TF–gene pair + expression + network + data_type
                      |
             Minimum adjustment set
                      |
       Original Stage 1 (no TF–condition interaction)
                      |
        Interval and optional correlation filtering
                      |
               Is condition supplied?
                 /             \
               No              Yes
               |                |
       Original Stage 2       data_type?
       without interaction    /       \
                            bulk    single_cell
                             |          |
                      Original      Hierarchical
                      Stage 2       Stage 2
                             \          /
                      Coefficient posterior,
                   direction and uncertainty
```

With condition, bulk preserves the existing `stage2_model` choice. Single-cell
uses the hierarchy regardless of that legacy choice. Without condition, both
types use the original no-interaction model. A malformed or one-level supplied
condition is an error, not equivalent to omitting condition.

## Hierarchical model

For cell i from biological sample s, all retained predictor TFs are centered
within sample: x_within[i,j] = x[i,j] - mean_s(x[,j]). They retain the same units
as Stage 1; no sample-specific standardization is performed.

```text
Y[i] ~ NB2(mu[i], phi)
log(mu[i]) = log(library_size[i] / mean_library_size)
           + alpha + sample_intercept_sd * intercept_raw[s]
           + W[i,] * zeta
           + sum_j x_within[i,j] * beta[j]
           + x_within[i,T] * ((condition[s] - 0.5) * delta
                             + sample_slope_sd * slope_raw[s])

intercept_raw[s], slope_raw[s] ~ Normal(0, 1), independently
sample_intercept_sd, sample_slope_sd ~ half-Normal(0, supplied scale)
delta ~ Normal(0, target_interaction_sd)
```

The condition main effect and optional batch contrasts remain in W. Random
intercepts and slopes use noncentered parameterization, with a shared slope SD
across conditions. Only the focal TF receives a random slope.

The coefficient prior follows the existing Stage 2 transfer:

```text
beta_prior_mean[j] = stage1_mean[j] + direction_effect * direction[j]
beta_prior_sd[j]   = max(beta_sd_floor, stage1_sd_multiplier * stage1_sd[j])
beta[j] ~ Normal(beta_prior_mean[j], beta_prior_sd[j])
```

Defaults are 0.2, 0.5 and 1.5, respectively; the default delta SD is 0.5 and
both random-effect SD scales are 1. Network confidence acts through Stage 1.
The earlier experimental `TF_stage1_hierarchical_nb_model.stan` used a Laplace
prior directly and remains an experimental reference. The production model
is `TF_stage2_hierarchical_nb_model.stan` with posterior-informed Normal priors.

Stage 1 uses uncentered expression while the hierarchy estimates within-sample
slopes. Transferred moments therefore provide empirical-Bayes regularization,
not an exact posterior transfer for an identical estimand. Reuse of the same
data across stages remains empirical Bayes. Centering does not change predictor
units; condition slopes should not be multiplied by uncentered group mean
expression and called a hierarchical total regulatory contrast.

## Validation and outputs

The hierarchy requires unique cell-to-sample alignment, nonmissing sample IDs,
one condition per sample, and at least two biological samples per group.
Each group must also have at least two samples with within-sample focal-TF
variation. Samples with no variation are retained with a warning if that guard
still passes; their slopes rely on pooling. A rank-deficient condition/batch
design is rejected. Missing sample metadata never silently selects the bulk model.

`stage2_fit$sample_map` records sample counts and within-sample variance.
`stage2_fit$stan_data` retains the aligned sample IDs, condition coding, and
prior parameters. Its `beta` draws are population within-sample coefficients;
the focal coefficient is the midpoint of the two condition slopes.
`beta_target_control = beta_target_mean - delta/2` and
`beta_target_disease = beta_target_mean + delta/2`.
`beta_target_sample` includes each sample's random deviation.

The existing direction summary uses P(beta[T] > 0) or P(beta[T] < 0) >= 0.95.
Check the condition slopes and delta separately before interpreting group
differences. The model must be assessed with R-hat, effective sample sizes,
divergences and tree-depth diagnostics. Optional cell-conditional LOO does not
measure prediction for a new biological sample.

## Example

```r
input <- create_TF_computation_input(
  seurat_obj = macrophages, target_tf = "EP300", target_gene = "MTOR",
  data_type = "single_cell", sample_column = "orig.ident",
  condition_column = "sample", control_level = "Normal", disease_level = "AAA"
)
result <- TF_regulatory_direction_computation(input = input, batch_column = "batch")
result$stage2_model
result$stage2_fit$sample_map
result$stage2_fit$fit$summary(c("beta_target_control", "beta_target_disease", "beta_target_delta"))
```

The activity/EM workflow is separate and retains its existing models. This
change routes the complete TF–gene regulatory workflow only.
