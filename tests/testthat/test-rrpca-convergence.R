# The RPCA objective and its proximal steps are invariant under transposition.
# These tests use that symmetry, closed-form shrinkage examples, and planted
# low-rank + sparse components rather than comparing singular-vector signs.

rrpca_test_corrupted_fixture <- function(m, n) {
  L <- nla_test_fixture(m, n, c(30, 15, 8))$A
  S <- matrix(0, m, n)
  i <- seq_len(50L)
  S[cbind((7L * i) %% m + 1L, (11L * i) %% n + 1L)] <- 20 * (-1)^i
  list(A = L + S, L = L, S = S)
}

testthat::test_that("deterministic robust PCA is transpose equivariant", {
  A <- outer(seq_len(12), c(0.2, -0.6, 0.3, 0.7, -0.5, 0.1)) / 12
  A[cbind(c(3, 9), c(2, 5))] <- c(3, -2)
  for (iterations in c(1L, 4L)) {
    fit <- rrpca(A, lambda = 0.3, rand = FALSE, maxiter = iterations)
    transposed <- rrpca(t(A), lambda = 0.3, rand = FALSE,
                        maxiter = iterations)
    testthat::expect_equal(fit$L, t(transposed$L), tolerance = 1e-10)
    testthat::expect_equal(fit$S, t(transposed$S), tolerance = 1e-10)
    testthat::expect_equal(fit$err, transposed$err, tolerance = 1e-10)
  }
})

testthat::test_that("randomized shrinkage expands a saturated target rank", {
  # For these matrices lambda = 3 makes S = 0 in the first proximal step.
  # All 20 singular values survive shrinkage, and the exact update is L = A.
  # Asking only for the initial 11 singular values discards nine survivors.
  # Small shapes exercise the full-SVD fallback; the last retains a small
  # randomized subproblem after increasing the requested rank.
  for (shape in list(c(80L, 60L), c(60L, 80L), c(180L, 140L))) {
    A <- matrix(0, shape[1L], shape[2L])
    diag(A)[seq_len(20L)] <- 10
    for (seed in c(107L, 211L, 307L)) {
      set.seed(seed)
      fit <- rrpca(A, lambda = 3, maxiter = 1)
      testthat::expect_lt(nla_test_norm(fit$L - A) / nla_test_norm(A), 1e-10)
      testthat::expect_equal(fit$S, matrix(0, nrow(A), ncol(A)),
                             tolerance = 1e-10)
      testthat::expect_lt(tail(fit$err, 1L), 1e-10)
    }
  }
})

testthat::test_that("robust PCA predicts a rank within both matrix dimensions", {
  A <- matrix(0, 6, 80)
  diag(A) <- 10
  progress <- capture.output(rrpca(A, lambda = 3, rand = FALSE,
                                   maxiter = 1, trace = TRUE))
  progress <- progress[grepl("target rank k =", progress, fixed = TRUE)]
  testthat::expect_length(progress, 1L)
  target <- as.integer(sub(".*target rank k = ([0-9]+).*", "\\1", progress))
  testthat::expect_gte(target, 1L)
  testthat::expect_lte(target, min(dim(A)))
})

testthat::test_that("robust PCA respects iteration limits and reports its residual", {
  fixture <- rrpca_test_corrupted_fixture(80, 60)
  for (randomized in c(FALSE, TRUE)) {
    for (iterations in c(1L, 3L)) {
      set.seed(819)
      fit <- rrpca(fixture$A, rand = randomized, maxiter = iterations,
                   tol = 1e-14)
      residual <- nla_test_norm(fixture$A - fit$L - fit$S) /
        nla_test_norm(fixture$A)
      testthat::expect_length(fit$err, iterations)
      testthat::expect_equal(tail(fit$err, 1L), residual, tolerance = 1e-12)
      testthat::expect_true(all(is.finite(fit$err) & fit$err >= 0))
      testthat::expect_s3_class(fit, "rrpca")
      testthat::expect_named(fit, c("L", "S", "err"))
    }
  }
})

testthat::test_that("robust PCA recovers low-rank and sparse components", {
  for (shape in list(c(80L, 60L), c(60L, 80L))) {
    fixture <- rrpca_test_corrupted_fixture(shape[1L], shape[2L])
    lambda <- max(shape)^-0.5
    objective <- 30 + 15 + 8 + lambda * sum(abs(fixture$S))
    for (randomized in c(FALSE, TRUE)) {
      for (seed in if (randomized) c(819L, 820L, 821L) else 819L) {
        set.seed(seed)
        fit <- rrpca(fixture$A, rand = randomized, maxiter = 100, tol = 1e-8)
        recovered_values <- svd(fit$L, nu = 0, nv = 0)$d
        testthat::expect_lt(nla_test_norm(fit$L - fixture$L) /
                             nla_test_norm(fixture$L), 1e-5)
        testthat::expect_lt(nla_test_norm(fit$S - fixture$S) /
                             nla_test_norm(fixture$S), 1e-5)
        testthat::expect_equal(sum(recovered_values > 1e-5), 3L)
        testthat::expect_equal(sum(recovered_values) + lambda * sum(abs(fit$S)),
                               objective, tolerance = 1e-6)
        testthat::expect_lte(tail(fit$err, 1L), 1e-8)
        testthat::expect_lte(length(fit$err), 100L)
        testthat::expect_true(all(is.finite(fit$err)))
      }
    }
  }
})

testthat::test_that("robust PCA scales with the input magnitude", {
  fixture <- rrpca_test_corrupted_fixture(80, 60)
  fit <- rrpca(fixture$A, rand = FALSE, maxiter = 5)
  scaled <- rrpca(2.5 * fixture$A, rand = FALSE, maxiter = 5)
  testthat::expect_equal(scaled$L, 2.5 * fit$L, tolerance = 1e-10)
  testthat::expect_equal(scaled$S, 2.5 * fit$S, tolerance = 1e-10)
  testthat::expect_equal(scaled$err, fit$err, tolerance = 1e-10)
})

testthat::test_that("feasibility alone does not terminate robust PCA", {
  # Each optimum has a simple dual certificate: for the constant matrices,
  # Y = A / sqrt(m*n); for the diagonals, Y = min(1, lambda) * A.
  # Feasible but suboptimal decompositions can have exactly zero residual.
  cases <- list(
    list(A = matrix(1, 2, 2), lambda = 2^-0.5, sparse = FALSE),
    list(A = matrix(1, 500, 4), lambda = 500^-0.5, sparse = FALSE),
    list(A = matrix(1, 100, 8), lambda = 100^-0.5, sparse = FALSE),
    list(A = diag(10), lambda = 2, sparse = FALSE),
    list(A = diag(10), lambda = 0.5, sparse = TRUE)
  )
  for (case in cases) {
    expected_L <- if (case$sparse) case$A * 0 else case$A
    expected_S <- case$A - expected_L
    optimum <- sum(svd(expected_L, nu = 0, nv = 0)$d) +
      case$lambda * sum(abs(expected_S))
    for (randomized in c(FALSE, TRUE)) {
      set.seed(912)
      fit <- rrpca(case$A, lambda = case$lambda, rand = randomized,
                   maxiter = 100, tol = 1e-9)
      objective <- sum(svd(fit$L, nu = 0, nv = 0)$d) +
        case$lambda * sum(abs(fit$S))
      testthat::expect_lt(nla_test_norm(fit$L - expected_L) /
                           nla_test_norm(case$A), 1e-8)
      testthat::expect_lt(nla_test_norm(fit$S - expected_S) /
                           nla_test_norm(case$A), 1e-8)
      testthat::expect_equal(objective, optimum, tolerance = 1e-8)
      testthat::expect_lte(tail(fit$err, 1L), 1e-9)
    }
  }
})

testthat::test_that("robust PCA convergence is not loosened by rescaling data", {
  A <- matrix(1, 2, 2)
  for (magnitude in c(1e-6, 1, 1e6)) {
    for (randomized in c(FALSE, TRUE)) {
      set.seed(913)
      fit <- rrpca(magnitude * A, rand = randomized, maxiter = 100,
                   tol = 1e-9)
      testthat::expect_equal(fit$L / magnitude, A, tolerance = 1e-8)
      testthat::expect_lt(nla_test_norm(fit$S / magnitude), 1e-8)
      testthat::expect_lte(tail(fit$err, 1L), 1e-9)
    }
  }
})

testthat::test_that("robust PCA refines a sketch when low accuracy stalls convergence", {
  fixture <- rrpca_test_corrupted_fixture(80, 60)
  set.seed(927)
  progress <- capture.output(
    fit <- rrpca(fixture$A, p = 0, q = 0, maxiter = 100, tol = 1e-8,
                  trace = TRUE)
  )
  testthat::expect_lt(nla_test_norm(fit$L - fixture$L) /
                       nla_test_norm(fixture$L), 1e-5)
  testthat::expect_lt(nla_test_norm(fit$S - fixture$S) /
                       nla_test_norm(fixture$S), 1e-5)
  testthat::expect_lte(tail(fit$err, 1L), 1e-8)
  progress <- progress[grepl("stationarity error =", progress, fixed = TRUE)]
  stationarity <- as.numeric(sub(".*stationarity error = ", "",
                                 tail(progress, 1L)))
  testthat::expect_length(stationarity, 1L)
  testthat::expect_lte(stationarity, 1e-8)
})
