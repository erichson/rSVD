# Regression tests for numerical correctness and API invariants.
# Keep random matrices small so these can run on CRAN and in CI.

testthat::test_that("deterministic QB respects the requested rank", {
  set.seed(711)
  A <- matrix(rnorm(40 * 20), 40, 20)
  out <- rqb(A, k = 5, rand = FALSE)
  testthat::expect_identical(dim(out$Q), c(40L, 5L))
  testthat::expect_identical(dim(out$B), c(5L, 20L))
  testthat::expect_equal(crossprod(out$Q), diag(5), tolerance = 1e-10)
  testthat::expect_equal(out$B, crossprod(out$Q, A), tolerance = 1e-10)
})

testthat::test_that("complex PCA scales by modulus and matches deterministic scores", {
  A <- matrix(c(1+2i, -2+3i, 4-1i, 2+0i,
                2-1i, 3+2i, -1+4i, 0-2i), 4, 2)
  fit <- rpca(A, k = 2, center = TRUE, scale = TRUE,
              rand = FALSE, retx = TRUE)
  centered <- sweep(A, 2, colMeans(A), "-")
  expected_scale <- sqrt(colSums(Mod(centered)^2) / (nrow(A)-1))
  testthat::expect_equal(fit$scale, expected_scale, tolerance = 1e-12)
  testthat::expect_true(is.numeric(fit$scale))
  testthat::expect_equal(fit$x, predict(fit, A), tolerance = 1e-10,
                         ignore_attr = TRUE)
  testthat::expect_equal(crossprod(fit$rotation), diag(2),
                         tolerance = 1e-10, ignore_attr = TRUE)
})

testthat::test_that("robust PCA handles the zero matrix without NaN or Inf", {
  A <- matrix(0, 8, 6)
  for (randomized in c(FALSE, TRUE)) {
    fit <- rrpca(A, rand = randomized, maxiter = 5)
    testthat::expect_true(all(is.finite(fit$L)))
    testthat::expect_true(all(is.finite(fit$S)))
    testthat::expect_equal(fit$L + fit$S, A, tolerance = 1e-12)
  }
})

testthat::test_that("robust PCA handles a zero-rank thresholding step", {
  # Every singular value can be below the current threshold.
  A <- matrix(0, 8, 6)
  fit <- rrpca(A, rand = FALSE, maxiter = 3)
  testthat::expect_identical(dim(fit$L), dim(A))
  testthat::expect_equal(fit$L, A, tolerance = 1e-12)
})

testthat::test_that("full-rank interpolative decomposition is an interpolation identity", {
  set.seed(712)
  A <- matrix(rnorm(12 * 7), 12, 7)
  for (randomized in c(FALSE, TRUE)) {
    fit <- rid(A, k = ncol(A), rand = randomized)
    testthat::expect_identical(dim(fit$Z), c(7L, 7L))
    testthat::expect_equal(fit$Z[, fit$idx, drop = FALSE],
                           diag(7), tolerance = 1e-9)
    testthat::expect_equal(fit$C %*% fit$Z, A, tolerance = 1e-9)
  }
})

testthat::test_that("randomized SVD returns orthonormal vectors and a competitive residual", {
  set.seed(713)
  A <- matrix(rnorm(65 * 35), 65, 35)
  k <- 8L
  fit <- rsvd(A, k = k, p = 10, q = 2)
  optimal <- svd(A)
  residual <- norm(A - fit$u %*% (fit$d * t(fit$v)), "F")
  best <- sqrt(sum(optimal$d[-seq_len(k)]^2))
  testthat::expect_equal(crossprod(fit$u), diag(k), tolerance = 1e-8)
  testthat::expect_equal(crossprod(fit$v), diag(k), tolerance = 1e-8)
  testthat::expect_gte(residual + 1e-9, best)
  testthat::expect_lt(residual, 1.3 * best)
})
