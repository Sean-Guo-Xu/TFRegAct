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
otherwise leave it as `NULL` (the package automatically selects the no-batch
Stan model).

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
  cores = 4,
  output = TRUE
)

activity_seurat <- activity_fit$seurat_object
head(activity_fit$pipeline_result$em_fit$activity_summary)
```

Detailed installation instructions, examples, input contracts, and result
interpretation will be added as the package interface stabilizes.
