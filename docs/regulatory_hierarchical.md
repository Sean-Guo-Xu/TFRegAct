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

Stage 1 retains its confidence-weighted directional Laplace prior. Stage 2
assigns each retained predictor a network-informed Normal prior:

```text
q[j] = (min(10, max(confidence[j], 1)) - 1) / 9
prior_mean[j] = direction[j] * (0.2 + 0.4 * q[j])
prior_variance[j] = 0.5 + 1.5 * q[j]
beta[j] ~ Normal(prior_mean[j], sqrt(prior_variance[j]))
```

The Normal scale argument above is the standard deviation. Direction is +1 for
activation, -1 for repression, and 0 for unknown or mixed direction. The
variance ranges from 0.5 at confidence 1 to 2 at confidence 10; known-direction
mean magnitudes range from 0.2 to 0.6. This fixed mapping does not change when
the confidence filtering threshold changes. The default delta SD is 0.5 and
both random-effect SD scales are 1. Stage 1 posterior summaries are used for
screening and initialization, not as Stage 2 prior parameters. Network
confidence and direction determine the Normal prior directly.
The production model is `TF_stage2_hierarchical_nb_model.stan`.
`TF_stage1_hierarchical_nb_model.stan` remains an experimental reference.

Stage 1 uses expression that has not been centered within biological sample,
whereas the hierarchy estimates within-sample slopes. These are different
estimands on the same predictor units; Stage 2 constructs its network-informed Normal prior
on that scale, rather than transferring a Stage 1 posterior for a different
estimand. The same expression observations are refitted after screening; this
is not a sequential update with an independent second dataset. Centering does
not change predictor units; condition slopes should not be multiplied by
uncentered group mean expression and called a hierarchical total regulatory
contrast.

The old `direction_effect`, `beta_sd_floor` and `stage1_sd_multiplier`
arguments default to `NULL` and are deprecated: non-`NULL` values are ignored
with a warning. Saved Stage 2 fits from either the directional-Laplace or
earlier posterior-informed Normal version must be refitted before interpretation
under the current network-informed Normal model.

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

The activity workflow uses this network-informed Normal prior in its
non-hierarchical Stage 2 model. Its Stage 3 sampler or quadrature EM fit uses
a Gaussian approximation to the full Stage 2 posterior of the focal-TF
coefficient and optional interaction. That approximation retains both
posterior means and their covariance, without SD inflation or an added variance
floor. The activity workflow's routing remains separate from this
biological-sample hierarchy.
