#!/usr/bin/env Rscript

# Shared, pre-computation input object for TF activity and TF -> gene analyses.

if (!methods::isClass("TFComputationInput")) {
  methods::setClass(
    "TFComputationInput",
    slots = c(
      seurat_obj = "ANY",
      target_tf = "character",
      target_gene = "character",
      condition_column = "character",
      control_level = "character",
      disease_level = "character",
      assay = "character",
      layer = "character",
      Y_exp = "ANY",
      libsize = "ANY",
      network_edge_file = "character"
    )
  )
}

tf_computation_input_stop <- function(...) stop(sprintf(...), call. = FALSE)

#' Construct shared input data and metadata for TF computation pipelines.
#'
#' `seurat_obj` is retained by reference in the S4 object. Subset the object
#' before construction; this class intentionally does not encode cell-type or
#' batch-column choices. Both downstream pipelines use the same data object,
#' target TF, and (when applicable) target gene.
create_TF_computation_input <- function(
  seurat_obj,
  target_tf,
  target_gene = NULL,
  condition_column = "sample",
  control_level = "Normal",
  disease_level = "AAA",
  assay = "RNA",
  layer = "data",
  Y_exp = NULL,
  libsize = NULL,
  network_edge_file = .tfregact_default_network_file()
) {
  if (is.null(seurat_obj)) tf_computation_input_stop("`seurat_obj` must be supplied.")
  metadata <- tryCatch(seurat_obj[[]], error = function(e) NULL)
  if (!is.data.frame(metadata)) tf_computation_input_stop("`seurat_obj` must be a Seurat-like object supporting `object[[]]` metadata extraction.")

  scalar_character <- function(value, name) {
    value <- trimws(as.character(value))
    if (length(value) != 1L || is.na(value) || !nzchar(value)) tf_computation_input_stop("`%s` must be one non-empty string.", name)
    value
  }
  target_tf <- scalar_character(target_tf, "target_tf")
  target_gene <- if (is.null(target_gene)) NA_character_ else scalar_character(target_gene, "target_gene")
  condition_column <- scalar_character(condition_column, "condition_column")
  control_level <- scalar_character(control_level, "control_level")
  disease_level <- scalar_character(disease_level, "disease_level")
  assay <- scalar_character(assay, "assay")
  layer <- scalar_character(layer, "layer")
  network_edge_file <- scalar_character(network_edge_file, "network_edge_file")
  if (!(condition_column %in% colnames(metadata))) tf_computation_input_stop("Condition metadata column `%s` was not found.", condition_column)

  methods::new(
    "TFComputationInput", seurat_obj = seurat_obj, target_tf = target_tf,
    target_gene = target_gene,
    condition_column = condition_column, control_level = control_level,
    disease_level = disease_level, assay = assay, layer = layer, Y_exp = Y_exp,
    libsize = libsize, network_edge_file = network_edge_file
  )
}

tf_computation_input_values <- function(input) {
  if (!methods::is(input, "TFComputationInput")) tf_computation_input_stop("`input` must be a `TFComputationInput` object created by `create_TF_computation_input()`.")
  list(
    seurat_obj = methods::slot(input, "seurat_obj"),
    target_tf = methods::slot(input, "target_tf"),
    target_gene = methods::slot(input, "target_gene"),
    condition_column = methods::slot(input, "condition_column"),
    control_level = methods::slot(input, "control_level"),
    disease_level = methods::slot(input, "disease_level"),
    assay = methods::slot(input, "assay"), layer = methods::slot(input, "layer"),
    Y_exp = methods::slot(input, "Y_exp"), libsize = methods::slot(input, "libsize"),
    network_edge_file = methods::slot(input, "network_edge_file")
  )
}
