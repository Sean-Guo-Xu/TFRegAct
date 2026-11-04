#' Smooth inferred TF activity using UMI-adaptive KNN weights
#'
#' An optional post-processing step for scores that have already been inferred.
#' With `batch = NULL`, supply externally computed UMAP coordinates. When batch
#' labels are supplied, use the integrated.cca embedding. Coordinates are not
#' computed or batch-corrected by this function.
#'
#' @param score A named nonnegative numeric vector (one TF), or a numeric matrix
#'   with TFs in rows and uniquely named cells in columns.
#' @param coordinate Numeric matrix or data frame with uniquely named cells in
#'   rows. UMAP must have two columns named UMAP_1 and UMAP_2 (case and separators
#'   are ignored). For integrated.cca, use at least two columns named
#'   integratedcca_1, ...; the Seurat key integrateddr_ is also accepted.
#' @param umi Positive numeric UMI counts. Named values are aligned by cell;
#'   unnamed values must be in score-cell order.
#' @param batch NULL for UMAP, or a vector of nonmissing batch labels, one per
#'   score cell, for integrated.cca. Named labels are aligned by cell.
#' @param k Number of neighbours, excluding the query cell. Default 15.
#' @return A numeric vector or matrix matching the score's cell order and names.
#'   Original scores are not modified. Posterior uncertainty is not recomputed.
#' @details The retained fraction is U_i / (U_i + mean neighbour UMI).
#'   Remaining weight is assigned to neighbours using normalized inverse-distance
#'   weights. Exactly coincident neighbours share that weight equally.
#'   Batch labels select the coordinate type; they do not partition neighbours.
#'   Extra named coordinate rows, UMI values and batch labels are allowed.
#'   Coordinate labels are checked, but the caller is responsible for supplying
#'   the specified reduction rather than merely relabelling another embedding.
#' @export
smooth_TF_activity <- function(score, coordinate, umi, batch = NULL, k = 15L) {
  fail <- function(...) stop(sprintf(...), call. = FALSE)
  real_numeric <- function(x) is.numeric(x) && !is.complex(x)
  valid_names <- function(x) {
    !is.null(x) && !anyNA(x) && all(nzchar(trimws(x))) && !anyDuplicated(x)
  }
  is_vector <- real_numeric(score) && is.null(dim(score))
  if (is_vector) {
    cells <- names(score)
    if (!valid_names(cells)) fail("`score` must have unique, nonempty cell names.")
    scores <- matrix(score, nrow = 1L, dimnames = list(NULL, cells))
  } else {
    if (!is.matrix(score) || !real_numeric(score)) {
      fail("`score` must be a named numeric vector or a numeric TF-by-cell matrix.")
    }
    cells <- colnames(score)
    if (!valid_names(cells)) fail("`score` columns must have unique, nonempty cell names.")
    scores <- score
  }
  n <- ncol(scores)
  if (nrow(scores) < 1L || n < 2L) fail("`score` must contain at least one TF and two cells.")
  if (anyNA(scores) || any(!is.finite(scores)) || any(scores < 0)) {
    fail("`score` values must be finite and nonnegative.")
  }
  if (!real_numeric(k) || !is.null(dim(k)) || length(k) != 1L || is.na(k) || !is.finite(k) ||
      k < 1 || k != floor(k) || k >= n) {
    fail("`k` must be a positive integer smaller than the number of score cells (%d).", n)
  }
  k <- as.integer(k)

  align_vector <- function(x, label, positive = FALSE) {
    if (!is.atomic(x) || !is.null(dim(x))) fail("`%s` must be a vector.", label)
    if (!is.null(names(x))) {
      if (!valid_names(names(x))) fail("`%s` must have unique, nonempty cell names.", label)
      if (anyNA(match(cells, names(x)))) fail("`%s` is missing score cells.", label)
      x <- x[cells]
    }
    if (length(x) != n || anyNA(x)) fail("`%s` must contain one nonmissing value per score cell.", label)
    if (positive && (!real_numeric(x) || any(!is.finite(x)) || any(x <= 0))) {
      fail("`%s` values must be finite, positive numeric UMI counts.", label)
    }
    x
  }
  umi <- align_vector(umi, "umi", positive = TRUE)
  has_batch <- !is.null(batch)
  if (has_batch) {
    batch <- align_vector(batch, "batch")
    if ((!is.character(batch) && !is.factor(batch) && !real_numeric(batch)) ||
        any(!nzchar(trimws(as.character(batch)))) ||
        (is.numeric(batch) && any(!is.finite(batch)))) {
      fail("`batch` must contain nonempty, finite batch labels; use NULL for no batch.")
    }
  }

  if ((!is.matrix(coordinate) && !is.data.frame(coordinate)) ||
      (is.data.frame(coordinate) && !all(vapply(coordinate, real_numeric, logical(1))))) {
    fail("`coordinate` must be a numeric matrix or numeric data frame.")
  }
  coordinate <- as.matrix(coordinate)
  if (!real_numeric(coordinate) || nrow(coordinate) < 2L || ncol(coordinate) < 2L ||
      anyNA(coordinate) || any(!is.finite(coordinate))) {
    fail("`coordinate` must have at least two rows and columns, with all values finite and numeric.")
  }
  if (!valid_names(rownames(coordinate))) {
    fail("`coordinate` rows must have unique, nonempty cell names.")
  }
  if (!valid_names(colnames(coordinate))) {
    fail("`coordinate` columns must identify the UMAP or integrated.cca dimensions.")
  }
  dimension_names <- gsub("[^a-z0-9]", "", tolower(colnames(coordinate)))
  if (has_batch) {
    cca <- grepl("^integrated(cca|dr)[1-9][0-9]*$", dimension_names)
    if (!all(cca) || anyDuplicated(dimension_names) ||
        length(unique(sub("[0-9]+$", "", dimension_names))) != 1L) {
      fail(paste0("With `batch` supplied, `coordinate` must contain integrated.cca dimensions ",
                  "(e.g. integratedcca_1, integratedcca_2; integrateddr_ is also accepted)."))
    }
  } else if (ncol(coordinate) != 2L || !setequal(dimension_names, c("umap1", "umap2"))) {
    fail("With `batch = NULL`, `coordinate` must contain exactly the UMAP_1 and UMAP_2 dimensions.")
  }
  index <- match(cells, rownames(coordinate))
  if (anyNA(index)) fail("`coordinate` is missing score cells; row names must match score-cell names.")
  coordinate <- coordinate[index, , drop = FALSE]

  knn <- tryCatch(RANN::nn2(coordinate, coordinate, k = k + 1L, searchtype = "standard"),
                  error = function(e) fail("`coordinate` could not be used for KNN search: %s", conditionMessage(e)))
  if (anyNA(knn$nn.idx) || any(knn$nn.idx < 1L | knn$nn.idx > n) ||
      any(!is.finite(knn$nn.dists)) || any(knn$nn.dists < 0)) {
    fail("`coordinate` produced invalid KNN distances or indices; check its numerical scale.")
  }
  result <- scores
  for (i in seq_len(n)) {
    positions <- which(knn$nn.idx[i, ] != i)
    if (length(positions) < k) fail("Unable to find %d neighbours after excluding cell %s.", k, cells[i])
    positions <- positions[seq_len(k)]
    neighbours <- knn$nn.idx[i, positions]
    distance <- knn$nn.dists[i, positions]
    coincident <- distance == 0
    if (any(coincident)) {
      weights <- as.numeric(coincident) / sum(coincident)
    } else {
      # Scaling all inverse distances by the smallest distance avoids overflow.
      weights <- min(distance) / distance
      weights <- weights / sum(weights)
    }
    neighbour_umi <- mean(umi[neighbours])
    scale <- max(umi[i], neighbour_umi)
    retention <- (umi[i] / scale) / (umi[i] / scale + neighbour_umi / scale)
    result[, i] <- retention * scores[, i] +
      (1 - retention) * drop(scores[, neighbours, drop = FALSE] %*% weights)
  }
  if (any(!is.finite(result))) fail("Smoothing produced non-finite scores; check the numerical scale of `score`.")
  if (is_vector) return(stats::setNames(as.numeric(result[1L, ]), cells))
  result
}
