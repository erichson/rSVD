# Small fixtures with a prescribed spectrum keep approximation tests independent
# of singular-vector signs, the random-number generator, and the BLAS library.
nla_test_adjoint <- function(x) Conj(t(x))

nla_test_norm <- function(x) sqrt(sum(Mod(x)^2))

nla_test_fixture <- function(m, n, d, complex = FALSE) {
  basis <- function(size, rank) {
    if (complex) {
      return(exp(2i * pi * outer(seq_len(size) - 1,
                                seq_len(rank) - 1) / size) / sqrt(size))
    }
    out <- sqrt(2 / size) * cos(pi * outer(seq_len(size) - 0.5,
                                          seq_len(rank) - 1) / size)
    out[, 1L] <- out[, 1L] / sqrt(2)
    out
  }
  u <- basis(m, length(d))
  v <- basis(n, length(d))
  list(A = u %*% (d * nla_test_adjoint(v)), u = u, v = v, d = d)
}

nla_test_reconstruct <- function(fit) {
  fit$u %*% (fit$d * nla_test_adjoint(fit$v))
}

nla_test_expect_orthonormal <- function(x, info = NULL) {
  testthat::expect_lt(nla_test_norm(nla_test_adjoint(x) %*% x -
                                    diag(ncol(x))), 1e-9, label = info)
}

nla_test_expect_svd <- function(A, fit, k, info = NULL) {
  testthat::expect_s3_class(fit, "rsvd")
  testthat::expect_identical(dim(fit$u), c(nrow(A), as.integer(k)))
  testthat::expect_identical(dim(fit$v), c(ncol(A), as.integer(k)))
  testthat::expect_true(is.matrix(fit$u) && is.matrix(fit$v), info = info)
  testthat::expect_length(fit$d, k)
  testthat::expect_true(all(is.finite(fit$d)) && all(fit$d >= 0),
                        info = info)
  nla_test_expect_orthonormal(fit$u, info)
  nla_test_expect_orthonormal(fit$v, info)
}

nla_test_expect_exact <- function(A, approximation, info = NULL) {
  A <- as.matrix(A)
  testthat::expect_lt(nla_test_norm(A - approximation),
                      1e-9 * max(1, nla_test_norm(A)), label = info)
}

nla_test_expect_qb <- function(A, fit, width, info = NULL) {
  testthat::expect_s3_class(fit, "rqb")
  testthat::expect_identical(dim(fit$Q), c(nrow(A), as.integer(width)))
  testthat::expect_identical(dim(fit$B), c(as.integer(width), ncol(A)))
  testthat::expect_true(is.matrix(fit$Q) && is.matrix(fit$B), info = info)
  nla_test_expect_orthonormal(fit$Q, info)
  nla_test_expect_exact(nla_test_adjoint(fit$Q) %*% as.matrix(A),
                        fit$B, info)
}
