# TFRegAct

TFRegAct is an R package for Bayesian inference of transcription-factor (TF)
regulation and cell-level TF activity from single-cell RNA-seq data.

The package provides four main steps:

1. Build a confidence-weighted human TF--gene regulatory network from multiple
   evidence sources.
2. Construct a shared `TFComputationInput` object from a Seurat object and the
   target TF (and, for edge-level inference, target gene).
3. Estimate the direction and strength of one TF-to-gene regulatory edge using
   adjustment-set selection and two-stage negative-binomial Stan models.
4. Infer relative cell-level TF activity through target-gene selection,
   screening, MCMC refinement, and latent-activity estimation.

The regulatory FullMap is stored in the TFRegAct user cache rather than the
installed package directory. Build it once with:

```r
library(TFRegAct)
build_tf_network()
```

By default this uses `tools::R_user_dir("TFRegAct", "cache")`; the network
location can be overridden with `options(TFRegAct.network_path = "<path>")`.

## Basic usage: TF-to-gene regulatory direction

Subset the Seurat object to the cell population of interest before creating the
input object. The example below estimates the direction and effect size of one
TF-to-gene edge. Supply `batch_column` only when batch information is present;
otherwise leave it as `NULL`. Condition is also optional for regulatory-edge
inference: set `condition_column = NULL` when it is unavailable or should not
be adjusted for. Internally, available condition and batch variables are
encoded in one nuisance design matrix, so Stage 1 and Stage 2 each use a single
Stan model for all four combinations (neither, either one, or both).

```r
library(TFRegAct)

edge_input <- create_TF_computation_input(
  seurat_obj = macrophages,
  target_tf = "ATF4",
  target_gene = "ATF3",
  condition_column = "sample",
  control_level = "Normal",
  disease_level = "AAA"
)

edge_fit <- TF_regulatory_direction_computation(
  input = edge_input,
  batch_column = "batch", # use NULL when no batch column exists
  output = TRUE
)

edge_fit$direction_summary
```

## Basic usage: cell-level TF activity

For activity inference, omit `target_gene`. The function identifies direct
targets, finds target-specific adjustment sets, performs screening and MCMC
refinement, then returns a Seurat object containing `<TF>_activity_A`.
By default, direct target edges and prior-network edges must have confidence
at least 4. Adjustment sets are selected from eight randomized starts, using
the smallest valid set first and path-edge confidence to break size ties.
The activity prescreen and regulatory Stage 1 share the same directional
Laplace Stan model and the same fitting function. Condition and batch terms are
optional and are encoded in one nuisance design matrix. Activity uses the
shared model with mean-field variational inference. Regulatory-edge Stage 1
uses MCMC by default, but can instead use mean-field or full-rank variational
inference for faster approximate screening. The subsequent activity Normal-prior MCMC and
regulatory Stage 2 call the same negative-binomial Stan model and MCMC fitting
function. In both workflows, the retained Stage 1 posterior is converted in R
to `beta_prior_mean` and `beta_prior_sd`; network confidence is not applied a
second time in Stage 2. This transfer is controlled by `direction_effect`,
`beta_sd_floor`, and `stage1_sd_multiplier`. Set `condition_column = NULL` to
run the entire activity pipeline without condition effects or TF-by-condition
interactions; Stage 3 then reports overall cell activity and omits the
control-versus-disease contrast tables.

Stage 3 uses a soft spike-and-slab activity gate. When the normalized target-TF
expression is zero, the prior probability that the TF is active is controlled
by `active_prior_zero` (default `0.2`); when expression is detected, this prior
probability is fixed at `0.9`. Downstream target-gene evidence then updates the
cell-specific activity probability. Thus, zero expression strongly downweights
activity but does not force activity to zero, and no cell-neighbour smoothing
or expression imputation is applied.

```r
activity_input <- create_TF_computation_input(
  seurat_obj = macrophages,
  target_tf = "IRF7",
  condition_column = "sample",
  control_level = "Normal",
  disease_level = "AAA"
)

activity_fit <- TF_activity_computation(
  input = activity_input,
  batch_column = "batch", # use NULL when no batch column exists
  active_prior_zero = 0.2, # soft activity prior when TF expression is zero
  cores = 4,
  output = TRUE
)

activity_seurat <- activity_fit$seurat_object
head(activity_fit$pipeline_result$em_fit$activity_summary)
```

Arguments omitted from the examples use the defaults below. In most analyses,
users only need to set the input data, TF/gene, condition and batch columns,
confidence thresholds, and output behaviour. Prior, sampling, and numerical
defaults should normally be changed only for a planned sensitivity analysis or
when diagnostics indicate a problem.

## Complete arguments: regulatory-direction inference

`TF_regulatory_direction_computation()` accepts the following arguments.

### Input, subsetting, and output

| Argument | Default | Description |
|---|---:|---|
| `target_tf` | `NULL` | Target TF symbol. It may instead be supplied through `input`. |
| `target_gene` | `NULL` | Target gene symbol. It may instead be supplied through `input`. |
| `input` | `NULL` | A `TFComputationInput` object. Do not also supply `seurat_obj` or `data_file`. |
| `seurat_obj` | `NULL` | In-memory Seurat object when `input` is not used. |
| `data_file` | `NULL` | `.RData` file containing an object named `pbmc`; alternative to `seurat_obj`. |
| `output` | `FALSE` | `FALSE` returns results in memory, `TRUE` creates `<TF>_<gene>_output`, and a character value specifies a directory. |
| `cell_column` | `NULL` | Metadata column used for cell-population subsetting. |
| `cell_levels` | `NULL` | Values retained from `cell_column`; `NULL` keeps all cells. |
| `batch_subset` | `NULL` | Optional batch values to retain before fitting. |
| `batch_column` | `NULL` | Batch metadata column; `NULL` fits without batch terms. |
| `condition_column` | `"sample"` | Two-level condition column; set to `NULL` for condition-free inference. |
| `control_level` | `"Normal"` | Reference condition label. |
| `disease_level` | `"AAA"` | Comparison condition label. |
| `assay` | `"RNA"` | Seurat assay used for expression extraction. |
| `layer` | `"data"` | Seurat normalized-expression layer. |
| `Y_exp` | `NULL` | Optional pre-extracted target-gene response; normally inferred from the Seurat object. |
| `libsize` | `NULL` | Optional precomputed library-size vector; normally extracted from counts. |

### Network, priors, and adjustment-set search

| Argument | Default | Description |
|---|---:|---|
| `network_edge_file` | package cache | FullMap regulatory-network file; normally left unchanged. |
| `confounder_confidence_threshold` | `3` | Minimum network confidence for candidate confounder edges. |
| `adjustment_search_starts` | `8` | Randomized recursive adjustment-set searches. |
| `adjustment_search_cores` | `4` | Cores used for adjustment-set search. |
| `max_adjustment_sets` | `NULL` | Deprecated compatibility argument; ignored in favour of `adjustment_search_starts`. |
| `dagitty_beta` | `2` | Relative weight given to confidence on paths towards the target gene when adjustment sets are scored. |
| `gamma` | `1` | Base Stage 1 beta-prior scale. This is the regulatory-interface counterpart of activity `beta_prior_scale`. |
| `eta` | `0.5` | Exponent controlling how strongly confidence changes the beta-prior scale. |
| `r_dir` | `3` | Additional shrinkage ratio on effects opposing the network direction. |
| `nuisance_prior_scale` | `1` | Normal-prior scale for available condition and batch columns in the nuisance design matrix. A scalar or one value per column is accepted. |
| `confidence_min` | `1` | Lower bound applied to network confidence. |
| `confidence_max` | `10` | Upper bound applied to network confidence. |

### Stage 1, filtering, and Stage 2

| Argument | Default | Description |
|---|---:|---|
| `stage1_inference` | `"mcmc"` | Stage 1 inference method: `"mcmc"` (default) or approximate `"variational"`. |
| `stage1_chains` | `3` | Stage 1 MCMC chains. |
| `stage1_iter_warmup` | `600` | Warmup iterations per Stage 1 chain. |
| `stage1_iter_sampling` | `1200` | Retained iterations per Stage 1 chain. |
| `stage1_adapt_delta` | `0.95` | Stage 1 HMC target acceptance probability. |
| `stage1_max_treedepth` | `12` | Stage 1 maximum HMC tree depth. |
| `stage1_variational_algorithm` | `"meanfield"` | Variational family used when `stage1_inference = "variational"`; use `"meanfield"` or `"fullrank"`. |
| `stage1_variational_iter` | `10000` | Maximum Stage 1 variational optimization iterations. |
| `stage1_variational_output_samples` | `2000` | Draws sampled from the fitted Stage 1 variational approximation. |
| `stage1_filter_interval` | `c(5, 95)` | Posterior percentile interval used to retain TF coefficients; the default is a central 90% interval. |
| `correlation_filter` | `TRUE` | Remove redundant, direction-consistent highly correlated TF predictors. |
| `correlation_threshold` | `0.7` | Absolute correlation threshold used by the redundancy filter. |
| `stage2_model` | `"no_interaction"` | Use `"no_interaction"` or `"target_interaction"`; the latter estimates a target-TF-by-condition effect. |
| `direction_effect` | `0.2` | Direction-dependent shift added to the Stage 1 beta mean when constructing the Stage 2 Normal prior. |
| `beta_sd_floor` | `0.5` | Minimum Stage 2 beta-prior standard deviation. |
| `stage1_sd_multiplier` | `1.5` | Multiplier applied to the Stage 1 beta posterior SD for the Stage 2 prior. |
| `target_interaction_sd` | `0.5` | Prior SD for the target-TF-by-condition interaction. |
| `stage2_chains` | `3` | Stage 2 MCMC chains. |
| `stage2_iter_warmup` | `500` | Warmup iterations per Stage 2 chain. |
| `stage2_iter_sampling` | `900` | Retained iterations per Stage 2 chain. |
| `stage2_adapt_delta` | `0.95` | Stage 2 HMC target acceptance probability. |
| `stage2_max_treedepth` | `12` | Stage 2 maximum HMC tree depth. |
| `compute_loo` | `FALSE` | Compute LOO from the generated pointwise log likelihood. |
| `seed` | `123` | Random seed. |
| `refresh` | `50` | CmdStan progress-printing interval. |
| `force_recompile` | `FALSE` | Recompile the Stan model even when a compiled executable is available. |

## Complete arguments: cell-level activity inference

`TF_activity_computation()` accepts the following arguments.

### Input, subsetting, and output

| Argument | Default | Description |
|---|---:|---|
| `target_tf` | `NULL` | Target TF symbol; it may instead be supplied through `input`. |
| `input` | `NULL` | A `TFComputationInput` without a target gene. Do not also supply `seurat_obj` or `data_file`. |
| `seurat_obj` | `NULL` | In-memory Seurat object when `input` is not used. |
| `data_file` | `NULL` | `.RData` file containing an object named `pbmc`; alternative to `seurat_obj`. |
| `output` | `FALSE` | `FALSE` returns results in memory, `TRUE` creates `<TF>_output`, and a character value specifies a directory. |
| `cell_column` | `NULL` | Metadata column used for cell-population subsetting. |
| `cell_levels` | `NULL` | Values retained from `cell_column`; `NULL` keeps all cells. |
| `batch_subset` | `NULL` | Optional batch values to retain before fitting. |
| `batch_column` | `NULL` | Batch metadata column; `NULL` fits without batch terms. |
| `condition_column` | `"sample"` | Two-level condition column; set to `NULL` for condition-free inference. |
| `control_level` | `"Normal"` | Reference condition label. |
| `disease_level` | `"AAA"` | Comparison condition label. |
| `assay` | `"RNA"` | Seurat assay used for expression extraction. |
| `layer` | `"data"` | Seurat normalized-expression layer. |
| `cell_limit` | `NULL` | Optional maximum cell count for pilot analyses; `NULL` uses every selected cell. |
| `cell_limit_seed` | `123` | Seed used only when cells are subsampled by `cell_limit`. |

### Network, priors, and Stage 1 prescreen

| Argument | Default | Description |
|---|---:|---|
| `cores` | `4` | Workers used for per-gene fitting and Stage 3 draw-level parallelism. |
| `seed` | `123` | Main pipeline random seed. |
| `network_edge_file` | package cache | FullMap regulatory-network file; normally left unchanged. |
| `target_confidence_threshold` | `4` | Minimum confidence for direct target-TF-to-gene edges entering activity inference. |
| `target_confidence_override` | `NULL` | Optional named values or data frame replacing selected target-edge confidence values. |
| `confounder_confidence_threshold` | `4` | Minimum network confidence for candidate confounder edges. |
| `adjustment_search_starts` | `8` | Randomized recursive adjustment-set searches per target gene. |
| `max_adjustment_sets` | `NULL` | Deprecated compatibility argument; ignored in favour of `adjustment_search_starts`. |
| `dagitty_beta` | `2` | Relative weight given to confidence on paths towards each target gene. |
| `beta_prior_scale` | `1` | Base Stage 1 beta-prior scale. |
| `eta` | `0.5` | Exponent controlling how strongly confidence changes the beta-prior scale. |
| `r_dir` | `3` | Additional shrinkage ratio on effects opposing the network direction. |
| `batch_prior_scale` | `1` | Normal-prior scale for batch contrasts. |
| `target_interaction_sd` | `0.5` | Prior SD for the target-TF-by-condition interaction. |
| `predictor_sd_min` | `1e-8` | Predictors with an SD at or below this value are treated as non-informative. |
| `standardize_predictors` | `TRUE` | Standardize TF predictors before regression. |
| `prescreen_target_interval` | `0.90` | Stage 1 interval used to retain target-gene models supported in at least one condition. |
| `prescreen_confounder_interval` | `0.90` | Stage 1 interval used to remove unsupported confounder TFs. |
| `prescreen_variational_iter` | `10000` | Maximum mean-field variational iterations per target gene. |
| `prescreen_output_samples` | `2000` | Variational posterior draws used for interval screening. |

### Stage 2 MCMC and Stage 3 activity

| Argument | Default | Description |
|---|---:|---|
| `direction_effect` | `0.2` | Direction-dependent shift added to the Stage 1 beta mean for the Stage 2 prior. |
| `beta_sd_floor` | `0.5` | Minimum Stage 2 beta-prior standard deviation. |
| `stage1_sd_multiplier` | `1.5` | Multiplier applied to Stage 1 beta posterior SDs for Stage 2 priors. |
| `stage2_alpha_prior_sd` | `1` | Stage 2 intercept-prior SD. |
| `mcmc_target_interval` | `0.90` | Stage 2 interval used to retain target-gene models. |
| `mcmc_confounder_interval` | `0.90` | Interval level recorded for confounder summaries. Stage 2 does not remove confounders. |
| `mcmc_chains` | `3` | Stage 2 MCMC chains. |
| `mcmc_iter_warmup` | `400` | Warmup iterations per Stage 2 chain. |
| `mcmc_iter_sampling` | `600` | Retained iterations per Stage 2 chain. |
| `mcmc_adapt_delta` | `0.95` | Stage 2 HMC target acceptance probability. |
| `mcmc_max_treedepth` | `12` | Stage 2 maximum HMC tree depth. |
| `mcmc_draws_for_em` | `1500` | Maximum aligned Stage 2 draws retained for possible Stage 3 uncertainty propagation. |
| `nuisance_draw_count` | `100` | Stage 2 nuisance draws fitted in Stage 3; `0` uses posterior means only. |
| `active_prior_zero` | `0.2` | Prior active probability for a cell whose normalized target-TF expression is zero. The corresponding prior zero-activity probability is `0.8`. |
| `em_control` | `list()` | Optional named list of advanced Stage 3 numerical controls described below. |
| `resume` | `TRUE` | Resume compatible pipeline checkpoints. |
| `force_refit` | `FALSE` | Refit completed stages rather than reusing compatible results. |
| `force_recompile` | `FALSE` | Recompile Stan models even when compatible executables exist. |

The active prior for cells with detected target-TF expression is fixed at
`0.9`; it is intentionally not a second user parameter. Consequently, only
one gate probability needs to be specified. The following advanced Stage 3
settings can be supplied through `em_control`:

| `em_control` entry | Default | Description |
|---|---:|---|
| `kappa` | `1` | Multiplier relating the positive-activity slab SD to the target-TF expression scale. |
| `expression_zero_tolerance` | `0` | Expression values at or below this threshold are treated as zero. |
| `quadrature_nodes` | `21` | Gauss-Legendre nodes used for each cell-specific positive-activity integral. |
| `max_iter` | `30` | Maximum EM iterations. |
| `min_iter` | `2` | Minimum EM iterations before convergence is allowed. |
| `beta_tolerance` | `1e-3` | Convergence tolerance for target beta updates. |
| `activity_tolerance` | `1e-3` | Convergence tolerance for cell-activity updates. |
| `objective_tolerance` | `1e-8` | Convergence tolerance for objective changes. |
| `monotonicity_tolerance` | `1e-4` | Permitted numerical deviation from monotone EM improvement. |
| `tail_log_drop` | `30` | Required log-density drop when truncating the positive-activity integration tail. |
| `max_upper_expansions` | `12` | Maximum expansions used to locate the integration upper bound. |
| `mstep_maxit` | `100` | Maximum optimizer iterations in each M-step. |
| `mstep_reltol` | `1e-8` | Relative optimizer tolerance in each M-step. |
| `checkpoint_every` | `NULL` | Optional checkpoint frequency across nuisance draws. |
| `save_traces` | `TRUE` in the pipeline | Retain EM iteration traces for diagnostics. |

For example, the biologically important activity-gate and slab-width settings
can be changed without specifying all other defaults:

```r
activity_fit <- TF_activity_computation(
  input = activity_input,
  active_prior_zero = 0.2,
  nuisance_draw_count = 100,
  em_control = list(kappa = 1),
  output = TRUE
)
```

Formal help pages and extended result-interpretation examples will be added as
the package interface stabilizes.
