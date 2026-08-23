#!/usr/bin/env Rscript

seurat_tf_require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required.", pkg))
  }
}

seurat_tf_read_adjustment_table <- function(adjustment_input) {
  if (is.data.frame(adjustment_input)) {
    df <- adjustment_input
  } else if (is.list(adjustment_input)) {
    if (!is.null(adjustment_input$filtered_recommended_adjustment_metrics) &&
        is.data.frame(adjustment_input$filtered_recommended_adjustment_metrics)) {
      df <- adjustment_input$filtered_recommended_adjustment_metrics
    } else if (!is.null(adjustment_input$recommended_adjustment_metrics) &&
               is.data.frame(adjustment_input$recommended_adjustment_metrics)) {
      df <- adjustment_input$recommended_adjustment_metrics
    } else if (!is.null(adjustment_input$files) &&
               !is.null(adjustment_input$files$recommended_adjustment_metrics_filtered) &&
               file.exists(adjustment_input$files$recommended_adjustment_metrics_filtered)) {
      df <- utils::read.csv(
        adjustment_input$files$recommended_adjustment_metrics_filtered,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    } else if (!is.null(adjustment_input$files) &&
               !is.null(adjustment_input$files$recommended_adjustment_metrics) &&
               file.exists(adjustment_input$files$recommended_adjustment_metrics)) {
      df <- utils::read.csv(
        adjustment_input$files$recommended_adjustment_metrics,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    } else {
      stop(
        paste(
          "`adjustment_input` list was provided, but no usable adjustment table was found.",
          "Expected one of:",
          "`filtered_recommended_adjustment_metrics`,",
          "`recommended_adjustment_metrics`,",
          "or `files$recommended_adjustment_metrics_filtered`."
        )
      )
    }
  } else if (is.character(adjustment_input) && length(adjustment_input) == 1 && file.exists(adjustment_input)) {
    df <- utils::read.csv(adjustment_input, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("`adjustment_input` must be a data.frame, a find_adjustment_dagitty result/list, or a valid CSV file path.")
  }

  if (!("tf" %in% colnames(df))) {
    stop("Adjustment table must contain a 'tf' column.")
  }

  df
}

seurat_tf_get_expression_matrix <- function(seurat_obj, assay = NULL, layer = "data") {
  seurat_tf_require_pkg("SeuratObject")

  if (is.null(assay) || !nzchar(assay)) {
    assay <- tryCatch(
      SeuratObject::DefaultAssay(seurat_obj),
      error = function(e) NULL
    )
  }
  if (is.null(assay) || !nzchar(assay)) {
    stop("Could not determine assay. Please provide `assay` explicitly.")
  }

  expr <- tryCatch(
    SeuratObject::GetAssayData(seurat_obj, assay = assay, layer = layer),
    error = function(e) {
      tryCatch(
        SeuratObject::GetAssayData(seurat_obj, assay = assay, slot = layer),
        error = function(e2) {
          stop(sprintf(
            "Failed to extract expression matrix from assay '%s' using layer/slot '%s'.",
            assay, layer
          ))
        }
      )
    }
  )

  expr
}

# Extract TF expression matrix from a Seurat object using TF names listed in a
# find_adjustment_dagitty output table. Returns a genes x cells matrix.
extract_adjustment_tf_expression <- function(
  seurat_obj,
  adjustment_input,
  assay = NULL,
  layer = "data",
  make_unique = TRUE
) {
  adj_df <- seurat_tf_read_adjustment_table(adjustment_input)
  expr_mat <- seurat_tf_get_expression_matrix(seurat_obj, assay = assay, layer = layer)

  tf_names <- trimws(as.character(adj_df$tf))
  tf_names <- tf_names[!is.na(tf_names) & tf_names != ""]
  if (make_unique) {
    tf_names <- unique(tf_names)
  }

  common_tfs <- intersect(tf_names, rownames(expr_mat))
  if (length(common_tfs) == 0) {
    stop("No overlapping TFs were found between the adjustment table and the Seurat expression matrix.")
  }

  tf_expr <- as.matrix(expr_mat[common_tfs, , drop = FALSE])
  attr(tf_expr, "requested_tfs") <- tf_names
  attr(tf_expr, "matched_tfs") <- common_tfs
  attr(tf_expr, "missing_tfs") <- setdiff(tf_names, common_tfs)
  attr(tf_expr, "assay") <- assay %||% NA_character_
  attr(tf_expr, "layer") <- layer
  tf_expr
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
