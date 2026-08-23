#!/usr/bin/env Rscript

tf_stop <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

tf_as_integer_vector <- function(x, name) {
  out <- suppressWarnings(as.integer(x))
  if (length(out) != length(x) || any(is.na(out) & !is.na(x))) {
    tf_stop("`%s` must be coercible to integer without introducing NA values.", name)
  }
  out
}

tf_gene_key <- function(x) {
  toupper(trimws(as.character(x)))
}

tf_validate_analysis_object <- function(x) {
  required_fields <- c(
    "expr",
    "batch",
    "target",
    "libsize",
    "confidence",
    "direction",
    "Y_exp"
  )

  missing_fields <- setdiff(required_fields, names(x))
  if (length(missing_fields) > 0) {
    tf_stop(
      "Analysis object is missing required field(s): %s",
      paste(missing_fields, collapse = ", ")
    )
  }

  expr <- x$expr
  if (!is.matrix(expr)) {
    tf_stop("`expr` must be a matrix.")
  }
  if (!is.numeric(expr)) {
    tf_stop("`expr` must be a numeric matrix.")
  }
  if (is.null(rownames(expr)) || any(rownames(expr) == "")) {
    tf_stop("`expr` must have non-empty rownames.")
  }
  if (is.null(colnames(expr)) || any(colnames(expr) == "")) {
    tf_stop("`expr` must have non-empty colnames.")
  }

  n_rows <- nrow(expr)
  n_cols <- ncol(expr)

  if (!is.null(x$batch)) {
    if (!is.factor(x$batch)) {
      tf_stop("`batch` must be a factor when provided.")
    }
    if (length(x$batch) != n_cols) {
      tf_stop("`batch` length must equal ncol(expr).")
    }
  }

  if (!is.null(x$sample)) {
    if (!is.factor(x$sample)) {
      tf_stop("`sample` must be a factor when provided.")
    }
    if (length(x$sample) != n_cols) {
      tf_stop("`sample` length must equal ncol(expr).")
    }
  }

  if (!is.character(x$target) || length(x$target) != 2) {
    tf_stop("`target` must be a character vector of length 2.")
  }
  if (any(is.na(x$target)) || any(trimws(x$target) == "")) {
    tf_stop("`target` must not contain NA or empty strings.")
  }
  target_gene <- trimws(as.character(x$target[[2]]))
  if (tf_gene_key(target_gene) %in% tf_gene_key(rownames(expr))) {
    tf_stop(
      "Target gene `%s` must not be included in `expr` predictors. It is already modeled as `Y_exp`.",
      target_gene
    )
  }

  if (!is.numeric(x$libsize) || length(x$libsize) != n_cols) {
    tf_stop("`libsize` must be a numeric vector with length equal to ncol(expr).")
  }

  if (!is.numeric(x$confidence) || length(x$confidence) != n_rows) {
    tf_stop("`confidence` must be a numeric vector with length equal to nrow(expr).")
  }

  if (!is.character(x$direction) || length(x$direction) != n_rows) {
    tf_stop("`direction` must be a character vector with length equal to nrow(expr).")
  }

  if (!is.integer(x$Y_exp) || length(x$Y_exp) != n_cols) {
    tf_stop("`Y_exp` must be an integer vector with length equal to ncol(expr).")
  }

  invisible(TRUE)
}

# Create a standardized TF analysis object for downstream modeling.
create_TF_analysis_object <- function(
  expr,
  batch,
  sample = NULL,
  target,
  libsize,
  confidence,
  direction,
  Y_exp
) {
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"

  obj <- list(
    expr = expr,
    batch = if (is.null(batch)) NULL else as.factor(batch),
    sample = if (is.null(sample)) NULL else as.factor(sample),
    target = as.character(target),
    libsize = as.numeric(libsize),
    confidence = as.numeric(confidence),
    direction = as.character(direction),
    Y_exp = tf_as_integer_vector(Y_exp, "Y_exp")
  )

  tf_validate_analysis_object(obj)
  class(obj) <- c("TFAnalysisObject", "list")
  obj
}

print.TFAnalysisObject <- function(x, ...) {
  cat("TFAnalysisObject\n")
  cat(sprintf("  Features: %d\n", nrow(x$expr)))
  cat(sprintf("  Cells: %d\n", ncol(x$expr)))
  cat(sprintf("  Target TF: %s\n", x$target[[1]]))
  cat(sprintf("  Target gene: %s\n", x$target[[2]]))
  if (is.null(x$batch)) {
    cat("  Batches: not provided\n")
  } else {
    cat(sprintf("  Batches: %d\n", nlevels(x$batch)))
  }
  if (is.null(x$sample)) {
    cat("  Samples: not provided\n")
  } else {
    cat(sprintf("  Samples: %d\n", nlevels(x$sample)))
  }
  invisible(x)
}

tf_get_default_libsize <- function(seurat_obj, assay = NULL) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    tf_stop("Package 'SeuratObject' is required to extract libsize from a Seurat object.")
  }

  if (is.null(assay) || !nzchar(assay)) {
    assay <- SeuratObject::DefaultAssay(seurat_obj)
  }

  counts <- tryCatch(
    SeuratObject::GetAssayData(seurat_obj, assay = assay, layer = "counts"),
    error = function(e) {
      SeuratObject::GetAssayData(seurat_obj, assay = assay, slot = "counts")
    }
  )
  Matrix::colSums(counts)
}

tf_get_feature_counts <- function(seurat_obj, feature, assay = NULL) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    tf_stop("Package 'SeuratObject' is required to extract counts from a Seurat object.")
  }

  if (is.null(assay) || !nzchar(assay)) {
    assay <- SeuratObject::DefaultAssay(seurat_obj)
  }

  counts <- tryCatch(
    SeuratObject::GetAssayData(seurat_obj, assay = assay, layer = "counts"),
    error = function(e) {
      SeuratObject::GetAssayData(seurat_obj, assay = assay, slot = "counts")
    }
  )

  if (!(feature %in% rownames(counts))) {
    feature_key <- tf_gene_key(feature)
    feature_match <- which(tf_gene_key(rownames(counts)) == feature_key)
    if (length(feature_match) == 1) {
      feature <- rownames(counts)[[feature_match]]
    }
  }

  if (!(feature %in% rownames(counts))) {
    tf_stop("Feature `%s` was not found in the counts matrix.", feature)
  }

  as.vector(counts[feature, ])
}

tf_extract_metadata_column <- function(seurat_obj, column, expected_name) {
  if (is.null(column)) {
    tf_stop("`%s` must be provided.", expected_name)
  }

  if (length(column) == 1 && is.character(column) && column %in% colnames(seurat_obj[[]])) {
    return(seurat_obj[[column]][, 1])
  }

  if (length(column) == ncol(seurat_obj)) {
    return(column)
  }

  tf_stop(
    "`%s` must be either a metadata column name in the Seurat object or a vector of length ncol(seurat_obj).",
    expected_name
  )
}

tf_extract_y_exp <- function(seurat_obj, Y_exp, target_gene, assay = NULL) {
  if (is.null(Y_exp)) {
    return(tf_as_integer_vector(tf_get_feature_counts(seurat_obj, target_gene, assay = assay), "Y_exp"))
  }

  if (length(Y_exp) == 1 && is.character(Y_exp)) {
    if (Y_exp %in% colnames(seurat_obj[[]])) {
      return(tf_as_integer_vector(seurat_obj[[Y_exp]][, 1], "Y_exp"))
    }
    return(tf_as_integer_vector(tf_get_feature_counts(seurat_obj, Y_exp, assay = assay), "Y_exp"))
  }

  if (length(Y_exp) == ncol(seurat_obj)) {
    return(tf_as_integer_vector(Y_exp, "Y_exp"))
  }

  tf_stop(
    "`Y_exp` must be NULL, a metadata column name, a feature name, or an integer-like vector of length ncol(seurat_obj)."
  )
}

# Build a TF analysis object directly from a Seurat object plus a
# find_adjustment_dagitty filtered output table.
build_TF_analysis_object_from_seurat <- function(
  seurat_obj,
  adjustment_input,
  target,
  batch = NULL,
  sample = NULL,
  Y_exp = NULL,
  assay = NULL,
  layer = "data",
  libsize = NULL,
  confidence_col = "overall_avg_confidence",
  direction_col = "effect_on_gene_A"
) {
  if (!exists("extract_adjustment_tf_expression", mode = "function")) {
    tf_stop(
      "Function `extract_adjustment_tf_expression` is not available. Please source seurat_tf_correlation.R first."
    )
  }
  if (!exists("seurat_tf_read_adjustment_table", mode = "function")) {
    tf_stop(
      "Function `seurat_tf_read_adjustment_table` is not available. Please source seurat_tf_correlation.R first."
    )
  }

  adj_df <- seurat_tf_read_adjustment_table(adjustment_input)
  target_gene <- trimws(as.character(target[[2]]))
  target_gene_key <- tf_gene_key(target_gene)
  outcome_in_adj <- tf_gene_key(adj_df$tf) == target_gene_key
  if (any(outcome_in_adj, na.rm = TRUE)) {
    warning(
      sprintf(
        "Target gene `%s` was found in adjustment table `tf` and has been removed before extracting predictor expression.",
        target_gene
      ),
      call. = FALSE
    )
    adj_df <- adj_df[!outcome_in_adj, , drop = FALSE]
  }
  if (nrow(adj_df) == 0) {
    tf_stop(
      "No adjustment TFs remain after removing target gene `%s` from `adjustment_input`.",
      target_gene
    )
  }

  expr <- extract_adjustment_tf_expression(
    seurat_obj = seurat_obj,
    adjustment_input = adj_df,
    assay = assay,
    layer = layer,
    make_unique = TRUE
  )
  outcome_in_expr <- tf_gene_key(rownames(expr)) == target_gene_key
  if (any(outcome_in_expr)) {
    warning(
      sprintf(
        "Target gene `%s` was found in the adjustment TF expression matrix and has been removed from `expr` predictors.",
        target_gene
      ),
      call. = FALSE
    )
    expr <- expr[!outcome_in_expr, , drop = FALSE]
  }
  if (nrow(expr) == 0) {
    tf_stop(
      "No predictor TFs remain after removing target gene `%s` from `expr`.",
      target_gene
    )
  }

  matched_tfs <- rownames(expr)
  adj_sub <- adj_df[match(tf_gene_key(matched_tfs), tf_gene_key(adj_df$tf)), , drop = FALSE]

  if (!(confidence_col %in% colnames(adj_sub))) {
    tf_stop("Adjustment table does not contain confidence column `%s`.", confidence_col)
  }
  if (!(direction_col %in% colnames(adj_sub))) {
    tf_stop("Adjustment table does not contain direction column `%s`.", direction_col)
  }

  if (is.null(libsize)) {
    libsize <- tf_get_default_libsize(seurat_obj, assay = assay)
  }

  batch_vec <- NULL
  if (!is.null(batch)) {
    batch_vec <- tf_extract_metadata_column(seurat_obj, batch, "batch")
  }
  sample_vec <- NULL
  if (!is.null(sample)) {
    sample_vec <- tf_extract_metadata_column(seurat_obj, sample, "sample")
  }
  y_vec <- tf_extract_y_exp(seurat_obj, Y_exp, target_gene = target[[2]], assay = assay)

  obj <- create_TF_analysis_object(
    expr = expr,
    batch = batch_vec,
    sample = sample_vec,
    target = target,
    libsize = libsize,
    confidence = as.numeric(adj_sub[[confidence_col]]),
    direction = as.character(adj_sub[[direction_col]]),
    Y_exp = y_vec
  )
  obj
}
