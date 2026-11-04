\name{smooth_TF_activity}
\alias{smooth_TF_activity}
\title{UMI-Adaptive KNN Smoothing of Inferred TF Activity}
\usage{
smooth_TF_activity(score, coordinate, umi, batch = NULL, k = 15L)
}
\arguments{
\item{score}{A named nonnegative numeric vector for one TF, or a numeric matrix
with TFs in rows and uniquely named cells in columns. All values must be finite.}
\item{coordinate}{Externally computed numeric coordinate matrix or data frame,
with unique cell row names. With \code{batch = NULL}, supply two UMAP columns
named \code{UMAP_1} and \code{UMAP_2}. With batch labels, supply at least two
integrated.cca columns named \code{integratedcca_1}, etc. The Seurat key
\code{integrateddr_} is also accepted. Case and separators are ignored.}
\item{umi}{Positive, finite numeric UMI counts. Named values are aligned to score
cells; unnamed values must be in score-cell order.}
\item{batch}{\code{NULL} for UMAP, or a vector of nonmissing batch labels for
integrated.cca. Named labels are aligned to score cells; unnamed labels must be
in score-cell order. A supplied label vector selects integrated.cca even if it
contains only one unique label.}
\item{k}{Number of neighbours excluding the query cell. An integer between one
and the number of score cells minus one. Default 15.}
}
\value{Smoothed numeric scores with the same vector or matrix form, cell order,
and names as \code{score}. Input scores are not modified.}
\description{An optional post-processing function called after activity inference.
It does not refit the inference model or compute coordinates.}
\details{
For cell \eqn{i}, the retained fraction is
\deqn{\alpha_i = U_i / (U_i + \overline{U}_{N_k(i)}).}
The smoothed score is
\deqn{A_i^{smooth} = \alpha_i A_i + (1-\alpha_i)\sum_{j\in N_k(i)}w_{ij}A_j.}
Neighbour weights are proportional to inverse Euclidean distance and sum to one.
Exactly coincident neighbours share all neighbour weight equally.

Extra named coordinate rows, UMI values, and batch labels are allowed and are
matched by cell name. Missing or duplicated cell identifiers, incompatible
dimension labels, non-finite inputs, and invalid \code{k} raise errors.
Coordinate labels are checked, but their computation history cannot be inferred
from a numeric matrix: callers must supply the specified reduction.

Batch labels select the coordinate type rather than restricting neighbour
search. To restrict smoothing to individual samples or other groups, subset the
inputs and call this function separately for each group. Posterior uncertainty
is not recomputed for the post-processed scores.
}
