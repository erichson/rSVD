testthat::test_that("PCA rejects undefined sample variances after row omission", {
  for (A in list(matrix(1:3, 1, 3), Matrix::Matrix(matrix(1:3, 1, 3),
                                                  sparse = TRUE))) {
    for (randomized in c(FALSE, TRUE)) {
      testthat::expect_error(rpca(A, k = 1, rand = randomized),
                             "at least two complete observations")
    }
  }
  A <- matrix(c(1, NA, 2, NA), 2, 2)
  testthat::expect_warning(
    testthat::expect_error(rpca(A, k = 1), "at least two complete observations"),
    "Missing values")
})

testthat::test_that("PCA keeps a single variable as a matrix", {
  A <- matrix(c(-3, 1, 2, 4), ncol = 1)
  for (input in list(A, Matrix::Matrix(A, sparse = TRUE))) {
    fit <- rpca(input, k = 1, scale = FALSE)
    testthat::expect_identical(dim(fit$rotation), c(1L, 1L))
    testthat::expect_identical(dim(fit$x), c(4L, 1L))
    testthat::expect_equal(fit$eigvals, stats::var(as.vector(A)))
    testthat::expect_equal(fit$var, stats::var(as.vector(A)))
    testthat::expect_equal(unname(fit$x %*% t(fit$rotation)),
                           unname(A - mean(A)), tolerance = 1e-12)
  }
})

testthat::test_that("fractional PCA and subset ranks truncate consistently", {
  A <- cbind(1:8, (1:8)^2, (-1)^(1:8), c(2, 3, 1, 8, 3, 7, 2, 9))
  for (fun in list(rpca, rid, rcur)) {
    for (randomized in c(FALSE, TRUE)) {
      set.seed(2401)
      expected <- fun(A, k = 2, rand = randomized)
      set.seed(2401)
      actual <- fun(A, k = 2.9, rand = randomized)
      testthat::expect_equal(actual, expected, tolerance = 1e-12)
    }
  }
})

testthat::test_that("invalid subset modes fail before returning incomplete factors", {
  for (mode in list("columns", "ROW", NA_character_, character(),
                    c("column", "row"))) {
    testthat::expect_error(rid(diag(3), k = 1, mode = mode), "mode")
  }
})

testthat::test_that("PCA prediction rejects a mismatched variable count", {
  fit <- rpca(matrix(seq_len(24), 6, 4), k = 1)
  for (n in c(1L, 3L, 5L)) {
    A <- matrix(seq_len(2L * n), 2, n)
    for (input in list(A, as.data.frame(A), Matrix::Matrix(A, sparse = TRUE))) {
      testthat::expect_error(predict(fit, input),
                             "same number of columns as the training data")
    }
  }
})

testthat::test_that("decomposition input checks retain finite-data semantics", {
  A <- diag(4)
  A[4, 4] <- Inf
  for (fun in list(rpca, rid, rcur, rrpca)) {
    for (input in list(A, Matrix::Matrix(A, sparse = TRUE))) {
      testthat::expect_error(fun(input), "finite|NA/NaN/Inf")
    }
    testthat::expect_error(fun(matrix(letters[1:4], 2)), "numeric")
    testthat::expect_error(fun(diag(4), rand = NA), "rand")
  }
  testthat::expect_error(rrpca(diag(3) + 1i), "real numeric")
})

testthat::test_that("invalid robust-PCA controls cannot bypass iteration", {
  for (A in list(diag(3), matrix(0, 3, 3))) {
    for (value in list(0, -1, Inf, NA_real_, numeric(), c(1, 2))) {
      testthat::expect_error(rrpca(A, lambda = value), "lambda")
      testthat::expect_error(rrpca(A, tol = value), "tol")
      testthat::expect_error(rrpca(A, maxiter = value), "maxiter")
    }
    testthat::expect_error(rrpca(A, maxiter = 1.5), "maxiter")
    testthat::expect_error(rrpca(A, trace = NA), "trace")
  }
})

testthat::test_that("a loose robust-PCA tolerance still returns a measured residual", {
  A <- diag(3)
  for (tol in c(1, 2)) {
    fit <- rrpca(A, rand = FALSE, tol = tol)
    testthat::expect_length(fit$err, 1L)
    testthat::expect_equal(fit$err, norm(A - fit$L - fit$S, "F") / norm(A, "F"))
  }
})

testthat::test_that("deterministic methods continue to ignore sampling controls", {
  A <- cbind(1:8, (1:8)^2, (-1)^(1:8))
  for (fun in list(rpca, rid, rcur, rrpca)) {
    expected <- fun(A, rand = FALSE)
    actual <- fun(A, rand = FALSE, p = NA, q = NA)
    testthat::expect_equal(actual, expected, tolerance = 1e-12)
  }
})
