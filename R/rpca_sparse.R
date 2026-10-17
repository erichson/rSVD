# Normalize structural, logical and pattern sparse matrices without a dense
# conversion. The resulting dgCMatrix has one stored entry per row/column pair;
# explicit zeros are permitted. Its documented p and x slots support stable
# column moments in O(nnz + n) storage.
rpca_sparse_numeric <- function(A) {
  methods::as(methods::as(methods::as(A, "generalMatrix"), "dMatrix"),
              "CsparseMatrix")
}

rpca_sparse_means <- function(A) {
  center <- Matrix::colMeans(A)
  # Matrix's sparse accumulation can lose low-order bits of a large constant.
  # Fully stored columns, and sums that overflowed in sparse accumulation,
  # can use base R's more accurate mean without filling implicit zeros.
  counts <- diff(A@p)
  stable <- which(counts == nrow(A) | !is.finite(center))
  for(j in stable) {
    center[j] <- mean(A@x[seq.int(A@p[j] + 1L, A@p[j + 1L])]) *
      (counts[j] / nrow(A))
  }
  center
}

rpca_sparse_norms <- function(A, center) {
  counts <- diff(A@p)
  if(is.logical(center)) {
    A@x <- A@x^2
    return(Matrix::colSums(A))
  }
  # Account separately for stored entries and implicit zeros. Subtracting
  # sum(x)^2/m from sum(x^2) is inaccurate for nearly constant columns.
  A@x <- (A@x - rep(center, counts))^2
  squared_norms <- Matrix::colSums(A)
  implicit <- counts < nrow(A)
  # Avoid 0 * Inf for a fully stored, very large constant column.
  squared_norms[implicit] <- squared_norms[implicit] +
    (nrow(A) - counts[implicit]) * center[implicit]^2
  squared_norms
}

rpca_sparse_products <- function(A, center, scale, replace_nan=FALSE) {
  counts <- diff(A@p)
  offset <- if(is.logical(center)) rep(0, ncol(A)) else center

  # A fully stored column can be centered explicitly without adding entries.
  # This also avoids cancellation in products with large, nearly constant
  # columns. Other column means remain a rank-one correction to each product.
  full <- counts == nrow(A)
  shift <- offset
  shift[!full] <- 0
  if(any(shift != 0)) A@x <- A@x - rep(shift, counts)
  offset[full] <- 0
  if(!is.logical(scale)) {
    A@x <- A@x / rep(scale, counts)
    offset <- offset / scale
  }
  if(replace_nan) {
    # predict.rpca historically sets NaNs to zero after scaling. An explicit
    # value equal to the implicit offset represents a zero transformed entry.
    missing <- is.nan(A@x)
    if(any(missing)) A@x[missing] <- rep(offset, counts)[missing]
  }
  A@factors <- list()
  centered <- any(offset != 0)

  forward <- function(X) {
    out <- as.matrix(A %*% X)
    if(centered) out <- sweep(out, 2, as.vector(crossprod(offset, X)), "-")
    out
  }
  adjoint <- function(X) {
    out <- as.matrix(crossprod(A, X))
    if(centered) out <- out - tcrossprod(offset, colSums(X))
    out
  }
  list(forward = forward, adjoint = adjoint,
       project = function(Q) H(adjoint(Q)),
       squared_norms = function() rpca_sparse_norms(A, offset))
}
