# Independent identities and known low-rank maps exercise numerical rank when
# the requested number of selected rows/columns exceeds the matrix rank.
testthat::test_that("pseudoinverse respects an analytic rank-one map", {
  set.seed(10)
  u <- rnorm(30)[c(22L, 24L, 12L)]
  v <- rnorm(18)
  inputs <- list(list(u = u, v = v),
                 list(u = u * (1 + 2i), v = v * (2 - 1i)))
  for (vectors in inputs) {
    # A = u v^T, so A^+ = Conj(v) Conj(u)^T / (||u||^2 ||v||^2).
    A <- outer(vectors$u, vectors$v)
    expected <- outer(Conj(vectors$v), Conj(vectors$u)) /
      (sum(Mod(vectors$u)^2) * sum(Mod(vectors$v)^2))
    for (transpose in c(FALSE, TRUE)) {
      X <- if (transpose) nla_test_adjoint(A) else A
      reference <- if (transpose) nla_test_adjoint(expected) else expected
      for (scale in c(1e-100, 1, 1e100)) {
        actual <- rsvd:::pinv(scale * X)
        target <- reference / scale
        testthat::expect_identical(dim(actual), rev(dim(X)))
        testthat::expect_lt(nla_test_norm(actual - target) /
                             nla_test_norm(target), 1e-12)
        testthat::expect_lt(nla_test_norm(scale * X %*% actual %*% (scale * X) -
                                          scale * X) /
                             nla_test_norm(scale * X), 1e-12)
        testthat::expect_lt(nla_test_norm(actual %*% (scale * X) %*% actual -
                                          actual) / nla_test_norm(actual), 1e-12)
        left <- (scale * X) %*% actual
        right <- actual %*% (scale * X)
        testthat::expect_lt(nla_test_norm(left - nla_test_adjoint(left)), 1e-12)
        testthat::expect_lt(nla_test_norm(right - nla_test_adjoint(right)), 1e-12)
      }
    }
  }
})

testthat::test_that("pseudoinverse retains small resolved singular directions", {
  # The third direction is small but still thousands of rounding-error units
  # above numerical zero. A sqrt(eps) cutoff would incorrectly discard it.
  for (complex in c(FALSE, TRUE)) {
    entries <- c(3, 2e-6, 1e-12)
    if (complex) entries <- entries * exp(1i * c(0.2, 0.7, -0.4))
    A <- matrix(if (complex) 0+0i else 0, 12, 7)
    A[cbind(1:3, 1:3)] <- entries
    expected <- matrix(if (complex) 0+0i else 0, 7, 12)
    expected[cbind(1:3, 1:3)] <- 1 / entries
    for (transpose in c(FALSE, TRUE)) {
      X <- if (transpose) nla_test_adjoint(A) else A
      reference <- if (transpose) nla_test_adjoint(expected) else expected
      for (scale in c(1e-100, 1, 1e100)) {
        actual <- rsvd:::pinv(scale * X)
        target <- reference / scale
        testthat::expect_lt(nla_test_norm(actual - target) /
                             nla_test_norm(target), 1e-12)
        testthat::expect_lt(nla_test_norm((scale * X) %*% actual %*% (scale * X) -
                                          scale * X) /
                             nla_test_norm(scale * X), 1e-12)
        # Every resolved coordinate must be retained, not just the dominant one.
        projector <- diag(c(rep(1, 3), rep(0, nrow(X) - 3L)))
        testthat::expect_lt(nla_test_norm((scale * X) %*% actual - projector),
                             1e-12)
      }
    }
  }
})

testthat::test_that("ID and CUR recover low-rank matrices when k exceeds rank", {
  # Seed 10 with rank one exposes a roundoff singular value in the selected
  # CUR rows; exact recovery is the oracle, independent of any pinv formula.
  for (seed in c(10L, 18L)) {
    for (rank in c(1L, 3L)) {
      set.seed(seed)
      left <- matrix(rnorm(30L * rank), 30, rank)
      right <- matrix(rnorm(rank * 18L), rank, 18)
      real <- left %*% right
      complex <- (left + 1i * matrix(rnorm(length(left)), 30, rank)) %*%
        (right + 1i * matrix(rnorm(length(right)), rank, 18))
      for (A in list(real, t(real), complex, nla_test_adjoint(complex))) {
        for (randomized in c(FALSE, TRUE)) {
          for (k in c(rank + 2L, min(dim(A)))) {
            set.seed(seed + 1000L)
            column <- rid(A, k = k, rand = randomized)
            testthat::expect_equal(column$C, A[, column$idx, drop = FALSE])
            testthat::expect_lt(nla_test_norm(column$Z[, column$idx, drop = FALSE] -
                                              diag(k)), 1e-12)
            testthat::expect_lt(nla_test_norm(A - column$C %*% column$Z) /
                                 nla_test_norm(A), 1e-10)
            set.seed(seed + 1000L)
            row <- rid(A, k = k, mode = "row", rand = randomized)
            testthat::expect_equal(row$R, A[row$idx, , drop = FALSE])
            testthat::expect_lt(nla_test_norm(A - row$Z %*% row$R) /
                                 nla_test_norm(A), 1e-10)
            set.seed(seed + 1000L)
            fit <- rcur(A, k = k, rand = randomized)
            testthat::expect_identical(dim(fit$U), c(k, k))
            testthat::expect_equal(fit$C, A[, fit$C.idx, drop = FALSE])
            testthat::expect_equal(fit$R, A[fit$R.idx, , drop = FALSE])
            testthat::expect_lt(nla_test_norm(A - fit$C %*% fit$U %*% fit$R) /
                                 nla_test_norm(A), 1e-10)
          }
        }
      }
    }
  }
})

testthat::test_that("over-requested sparse CUR retains scaling and dimensions", {
  set.seed(10)
  A <- matrix(rnorm(30), 30, 1) %*% matrix(rnorm(18), 1, 18)
  # Zero padding preserves rank one while making the input genuinely sparse.
  padded <- matrix(0, 42, 29)
  padded[seq_len(nrow(A)), seq_len(ncol(A))] <- A
  for (X in list(padded, t(padded))) {
    for (scale in c(1e-50, 1, 1e50)) {
      sparse <- Matrix::Matrix(scale * X, sparse = TRUE)
      original <- sparse
      for (randomized in c(FALSE, TRUE)) {
        set.seed(1010)
        fit <- rcur(sparse, k = 3, rand = randomized)
        testthat::expect_identical(dim(fit$U), c(3L, 3L))
        testthat::expect_lt(nla_test_norm(scale * X - fit$C %*% fit$U %*% fit$R) /
                             nla_test_norm(scale * X), 1e-10)
        testthat::expect_identical(sparse, original)
      }
    }
  }
})
