library(TFRegAct)
expect_error <- function(expr, pattern) {
  e <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(e, 'error'), grepl(pattern, conditionMessage(e)))
}
scores <- rbind(T1 = c(A = .9, B = .3, C = .6, D = .1),
                T2 = c(A = .2, B = .8, C = .4, D = .7))
original <- scores
coordinate <- cbind(UMAP_1 = c(A = 0, B = 1, C = 3, D = 10), UMAP_2 = 0)
umi <- c(A = 90, B = 30, C = 60, D = 100)
result <- smooth_TF_activity(scores, coordinate, umi, k = 2)
# A has neighbours B and C: distance weights 3/4 and 1/4, self retention 2/3.
expected <- (2 / 3) * scores[, 'A'] + (1 / 3) * (.75 * scores[, 'B'] + .25 * scores[, 'C'])
stopifnot(isTRUE(all.equal(unname(result[, 'A']), unname(expected), tolerance = 1e-12)),
          identical(dimnames(result), dimnames(scores)), identical(scores, original))
v <- scores['T1', ]
stopifnot(isTRUE(all.equal(smooth_TF_activity(v, coordinate, umi, k = 2), result['T1', ])))

# Independent input orders and extra named cells do not change the result.
extra <- rbind(coordinate, E = c(100, 100))
shuffled <- extra[c('E', 'D', 'C', 'A', 'B'), , drop = FALSE]
permuted_umi <- c(E = 1000, umi[c('D', 'B', 'A', 'C')])
stopifnot(isTRUE(all.equal(result, smooth_TF_activity(scores, as.data.frame(shuffled), permuted_umi, k = 2))),
          isTRUE(all.equal(result, smooth_TF_activity(scores, coordinate * 1e-12, umi, k = 2))))

# The batch route uses every integrated dimension, and allows cross-batch KNNs.
cca <- rbind(A = c(0, 0, 0), B = c(0, 0, 10), C = c(1, 0, 0))
colnames(cca) <- paste0('integratedcca_', 1:3)
batch_scores <- c(A = 0, B = 100, C = 2)
batch_umi <- c(A = 100, B = 100, C = 100)
batch <- c(C = 'b2', A = 'b1', B = 'b1')
b <- smooth_TF_activity(batch_scores, cca, batch_umi, batch = batch, k = 1)
stopifnot(b[['A']] == 1)
colnames(cca) <- paste0('integrated_dr_', 1:3)
stopifnot(isTRUE(all.equal(b, smooth_TF_activity(batch_scores, cca, batch_umi, batch, k = 1))))

# Coincident points: exclude self and share neighbour weight among ties.
ties <- matrix(0, 4, 2, dimnames = list(LETTERS[1:4], c('UMAP_1', 'UMAP_2')))
tie_scores <- c(A = 0, B = 2, C = 4, D = 6)
t <- smooth_TF_activity(tie_scores, ties, rep(100, 4), k = 3)
expected_ties <- .5 * tie_scores + .5 * (sum(tie_scores) - tie_scores) / 3
stopifnot(isTRUE(all.equal(t, expected_ties)))
two <- ties[1:2, , drop = FALSE]; two[2, 1] <- 1
stopifnot(isTRUE(all.equal(smooth_TF_activity(c(A = 0, B = 2), two, rep(1e308, 2), k = 1),
                          c(A = 1, B = 1))))

# Preserve constant signals and keep each TF within its original score range.
set.seed(42)
x <- matrix(runif(120), 3, dimnames = list(c('T1', 'T2', 'constant'), paste0('c', 1:40)))
x[3, ] <- .2
xy <- matrix(rnorm(80), 40, dimnames = list(colnames(x), c('UMAP_1', 'UMAP_2')))
y <- smooth_TF_activity(x, xy, runif(40, 10, 100), k = 15)
stopifnot(max(abs(y[3, ] - .2)) < 1e-14,
          all(y >= apply(x, 1, min) - 1e-14), all(y <= apply(x, 1, max) + 1e-14))

# Reject ambiguous or malformed coordinate inputs rather than falling back.
expect_error(smooth_TF_activity(v, NULL, umi, k = 2), 'coordinate')
bad <- coordinate; rownames(bad) <- NULL
expect_error(smooth_TF_activity(v, bad, umi, k = 2), 'row')
bad <- coordinate; rownames(bad)[2] <- 'A'
expect_error(smooth_TF_activity(v, bad, umi, k = 2), 'unique')
bad <- coordinate; rownames(bad)[2] <- 'unknown'
expect_error(smooth_TF_activity(v, bad, umi, k = 2), 'missing score cells')
bad <- coordinate; colnames(bad) <- NULL
expect_error(smooth_TF_activity(v, bad, umi, k = 2), 'columns')
bad <- coordinate; bad[1, 1] <- NA_real_
expect_error(smooth_TF_activity(v, bad, umi, k = 2), 'finite')
bad <- coordinate; bad[1, 1] <- Inf
expect_error(smooth_TF_activity(v, bad, umi, k = 2), 'finite')
expect_error(smooth_TF_activity(v, coordinate + 1i, umi, k = 2), 'numeric')
bad <- as.data.frame(coordinate); bad[, 1] <- as.character(bad[, 1])
expect_error(smooth_TF_activity(v, bad, umi, k = 2), 'numeric')
bad <- coordinate; colnames(bad) <- c('PC_1', 'PC_2')
expect_error(smooth_TF_activity(v, bad, umi, k = 2), 'UMAP')
expect_error(smooth_TF_activity(v, coordinate, umi, batch = rep('b', 4), k = 2), 'integrated.cca')
expect_error(smooth_TF_activity(batch_scores, cca, batch_umi, k = 1), 'UMAP')
bad <- cca; colnames(bad) <- c('integratedcca_1', 'integratedcca_2', 'harmony_3')
expect_error(smooth_TF_activity(batch_scores, bad, batch_umi, batch, k = 1), 'integrated.cca')

# Reject invalid scores, UMI, batch labels and k.
expect_error(smooth_TF_activity(unname(v), coordinate, umi, k = 2), 'cell names')
bad <- v; names(bad)[2] <- 'A'
expect_error(smooth_TF_activity(bad, coordinate, umi, k = 2), 'unique')
for (value in c(NA_real_, NaN, Inf, -.1)) {
  bad <- v; bad[1] <- value
  expect_error(smooth_TF_activity(bad, coordinate, umi, k = 2), 'finite and nonnegative')
}
expect_error(smooth_TF_activity(v + 1i, coordinate, umi, k = 2), 'numeric vector')
expect_error(smooth_TF_activity(scores[, 1, drop = FALSE], coordinate, umi, k = 1), 'two cells')
for (value in c(0, -1, Inf, NA_real_)) {
  bad <- umi; bad[1] <- value
  expect_error(smooth_TF_activity(v, coordinate, bad, k = 2), 'umi')
}
expect_error(smooth_TF_activity(v, coordinate, as.character(umi), k = 2), 'umi')
expect_error(smooth_TF_activity(v, coordinate, umi[-1], k = 2), 'missing score cells')
bad <- umi; names(bad)[2] <- 'A'
expect_error(smooth_TF_activity(v, coordinate, bad, k = 2), 'unique')
for (value in list(0, -1, 1.5, 4, NA_real_, Inf, '2', 1i, c(1, 2), matrix(1))) {
  expect_error(smooth_TF_activity(v, coordinate, umi, k = value), 'k')
}
for (labels in list(c('b1', NA, 'b2'), c('b1', '', 'b2'), c('b1', 'b2'), c(1, Inf, 2), TRUE)) {
  expect_error(smooth_TF_activity(batch_scores, cca, batch_umi, batch = labels, k = 1), 'batch')
}
message('UMI-adaptive smoothing, alignment, coordinate routing and invalid-input checks passed.')
