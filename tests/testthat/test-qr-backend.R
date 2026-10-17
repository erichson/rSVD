# The deterministic QB branch keeps the established base-QR column ordering.
# Compare with an independent full factorization, including cases where a
# truncated factorization would choose the wrong subspace without a fallback.
testthat::test_that("deterministic QB agrees with full QR on real rectangular inputs", {
  set.seed(821)
  for (shape in list(c(31L, 18L), c(18L, 31L))) {
    A <- matrix(rnorm(prod(shape)), shape[1L], shape[2L])
    dimnames(A) <- list(paste0("row", seq_len(nrow(A))),
                       paste0("col", seq_len(ncol(A))))
    for (requested in c(1L, 2.5, 5L, min(shape), min(shape) + 3L)) {
      k <- as.integer(min(requested, min(shape)))
      expected <- qr.Q(qr(A))[, seq_len(k), drop = FALSE]
      seed <- .Random.seed
      actual <- rqb(A, k = requested, rand = FALSE)
      testthat::expect_identical(.Random.seed, seed)
      testthat::expect_equal(actual$Q, expected, tolerance = 1e-12)
      testthat::expect_equal(actual$B, crossprod(expected, A), tolerance = 1e-12)
      nla_test_expect_qb(A, actual, k)
    }
  }
})

testthat::test_that("deterministic QB handles dependent leading columns and QR tolerance", {
  E <- diag(12)
  # The next useful direction occurs after the requested prefix; retaining the
  # QR of that prefix alone would return the wrong rank-two subspace.
  inputs <- list(cbind(0, E[, 1L], E[, 1L], E[, 2:8]),
                 cbind(E[, 1L], E[, 1L], E[, 2:9]),
                 matrix(0, 12, 8))
  for (epsilon in c(1e-9, 1e-5)) {
    inputs[[length(inputs) + 1L]] <-
      cbind(E[, 1L], E[, 1L] + epsilon * E[, 2L], E[, 3:10])
  }
  for (A in inputs) {
    for (k in c(1L, 2L, 4L)) {
      expected <- qr.Q(qr(A))[, seq_len(k), drop = FALSE]
      actual <- rqb(A, k = k, rand = FALSE)
      testthat::expect_equal(actual$Q, expected, tolerance = 1e-12)
      testthat::expect_equal(actual$B, crossprod(expected, A), tolerance = 1e-12)
      nla_test_expect_qb(A, actual, k)
    }
  }
})

testthat::test_that("deterministic QB preserves sparse and complex QR semantics", {
  set.seed(822)
  dense <- matrix(0, 36, 24)
  dense[1:12, 1L] <- seq_len(12)
  dense[13:24, 2L] <- seq_len(12) + 1
  dense[25:36, 3L] <- seq_len(12) + 2
  dense[, 4:24] <- dense[, 1:3] %*% matrix(rnorm(3 * 21), 3, 21)
  dimnames(dense) <- list(paste0("row", seq_len(nrow(dense))),
                         paste0("col", seq_len(ncol(dense))))
  sparse <- Matrix::Matrix(dense, sparse = TRUE)
  logical <- sparse != 0
  # Full LINPACK QR can overflow while repeatedly factoring roundoff in many
  # identical columns on some BLAS libraries. Keep this reference finite; the
  # repeated logical/pattern case has an analytic reference in the next test.
  logical[, 4:24] <- FALSE
  inputs <- list(sparse, as(sparse, "RsparseMatrix"),
                 as(sparse, "TsparseMatrix"),
                 Matrix::t(sparse),
                 Matrix::Matrix(cbind(0, dense), sparse = TRUE),
                 logical, as(logical, "nMatrix"),
                 Matrix::Diagonal(12, seq_len(12)))
  # Complex base QR uses numerical column pivoting, so columns outside a
  # requested prefix can change its leading basis vectors even at full rank.
  complex <- matrix(rnorm(28 * 16) + 1i * rnorm(28 * 16), 28, 16)
  complex[, 16L] <- 100 * complex[, 16L]
  inputs <- c(inputs, list(complex, Conj(t(complex)), matrix(0+0i, 12, 8)))
  for (A in inputs) {
    original <- A
    dense <- as.matrix(A)
    for (k in c(1L, 3L)) {
      expected <- qr.Q(qr(dense))[, seq_len(k), drop = FALSE]
      actual <- rqb(A, k = k, rand = FALSE)
      testthat::expect_equal(actual$Q, expected, tolerance = 1e-12)
      testthat::expect_equal(actual$B, nla_test_adjoint(expected) %*% dense,
                             tolerance = 1e-12)
      nla_test_expect_qb(A, actual, k)
    }
    testthat::expect_identical(A, original)
  }
})

testthat::test_that("deterministic QB handles repeated logical columns", {
  # The three disjoint indicator columns give an exact analytic basis. Many
  # repeated trailing columns expose a full-QR overflow on OpenBLAS 0.3.26;
  # compare projectors and reconstructions without using that failing oracle.
  indicators <- matrix(FALSE, 36, 3)
  indicators[cbind(seq_len(36), rep(1:3, each = 12))] <- TRUE
  dense <- cbind(indicators, matrix(TRUE, 36, 21))
  basis <- indicators / sqrt(12)
  sparse <- Matrix::Matrix(dense, sparse = TRUE)
  for (A in list(sparse, as(sparse, "nMatrix"))) {
    original <- A
    for (k in 1:3) {
      expected_projection <- tcrossprod(basis[, seq_len(k), drop = FALSE])
      actual <- rqb(A, k = k, rand = FALSE)
      testthat::expect_equal(tcrossprod(actual$Q), expected_projection,
                             tolerance = 1e-12)
      testthat::expect_equal(actual$Q %*% actual$B,
                             expected_projection %*% dense, tolerance = 1e-12)
      nla_test_expect_qb(A, actual, k)
      if (k == 3L) {
        testthat::expect_equal(actual$Q %*% actual$B, 1 * dense,
                               tolerance = 1e-12)
      }
    }
    testthat::expect_identical(A, original)
  }
})

testthat::test_that("deterministic QB rejects non-finite entries outside its prefix", {
  # Full QR has always rejected non-finite data anywhere in the matrix. A
  # prefix factorization must check later columns before using that shortcut.
  for (value in c(NA_real_, NaN, Inf, -Inf)) {
    for (column in c(1L, 6L)) {
      dense <- cbind(diag(8)[, 1:3], matrix(1, 8, 3))
      dense[8L, column] <- value
      sparse <- Matrix::Matrix(dense, sparse = TRUE)
      inputs <- list(dense, Matrix::Matrix(dense, sparse = FALSE), sparse,
                     as(sparse, "RsparseMatrix"), as(sparse, "TsparseMatrix"))
      for (A in inputs) {
        testthat::expect_error(qr(as.matrix(A)), "NA/NaN/Inf")
        testthat::expect_error(rqb(A, k = 2L, rand = FALSE), "NA/NaN/Inf")
      }
    }
  }
})
