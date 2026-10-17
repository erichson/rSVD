testthat::test_that("rsvd recovers exact low-rank tall and wide matrices", {
  fixture <- nla_test_fixture(32, 20, c(9, 4, 1))
  settings <- list(c(p = 0, q = 0), c(p = 4, q = 0), c(p = 4, q = 2))
  for (A in list(fixture$A, t(fixture$A))) {
    for (setting in settings) {
      set.seed(1001)
      fit <- rsvd(A, k = 3, p = setting[["p"]], q = setting[["q"]])
      info <- paste("shape", paste(dim(A), collapse = " x "),
                    "p", setting[["p"]], "q", setting[["q"]])
      nla_test_expect_svd(A, fit, 3, info)
      testthat::expect_equal(fit$d, fixture$d, tolerance = 1e-9)
      nla_test_expect_exact(A, nla_test_reconstruct(fit), info)
    }
  }
})

testthat::test_that("complex low-rank decompositions use Hermitian products", {
  fixture <- nla_test_fixture(27, 17, c(8, 3, 0.7), complex = TRUE)
  for (A in list(fixture$A, nla_test_adjoint(fixture$A))) {
    for (setting in list(c(p = 0, q = 0), c(p = 3, q = 2))) {
      set.seed(1002)
      fit <- rsvd(A, k = 3, p = setting[["p"]], q = setting[["q"]])
      nla_test_expect_svd(A, fit, 3)
      testthat::expect_equal(fit$d, fixture$d, tolerance = 1e-9)
      nla_test_expect_exact(A, nla_test_reconstruct(fit))

      set.seed(1002)
      qb <- rqb(A, k = 3, p = setting[["p"]], q = setting[["q"]])
      nla_test_expect_qb(A, qb, 3 + setting[["p"]])
      nla_test_expect_exact(A, qb$Q %*% qb$B)
    }
  }
})

testthat::test_that("decaying-spectrum residuals approach the optimal rank-k tail", {
  # Full rank alone does not preclude a useful low-rank approximation. Compare
  # to the best rank-k residual, not to zero or particular singular vectors.
  fixture <- nla_test_fixture(48, 28, 10^seq(0, -3, length.out = 28))
  k <- 6L
  optimal <- sqrt(sum(fixture$d[-seq_len(k)]^2))
  for (A in list(fixture$A, t(fixture$A))) {
    for (seed in c(12L, 72L, 908L)) {
      set.seed(seed)
      fit <- rsvd(A, k = k, p = 3, q = 2)
      residual <- nla_test_norm(A - nla_test_reconstruct(fit))
      nla_test_expect_svd(A, fit, k)
      testthat::expect_gte(residual + 1e-10, optimal)
      testthat::expect_lt(residual, 1.1 * optimal)
    }
  }
})

testthat::test_that("randomized QB retains oversampling and projection semantics", {
  fixture <- nla_test_fixture(30, 18, c(6, 2, 0.5))
  for (A in list(fixture$A, t(fixture$A))) {
    for (p in c(0L, 4L, 100L)) {
      for (q in c(0L, 2L)) {
        set.seed(1003)
        fit <- rqb(A, k = 3, p = p, q = q)
        width <- min(3 + p, dim(A))
        nla_test_expect_qb(A, fit, width)
        nla_test_expect_exact(A, fit$Q %*% fit$B)
      }
    }
  }
})

testthat::test_that("rank-one outputs retain matrix dimensions", {
  fixture <- nla_test_fixture(24, 15, 7)
  for (A in list(fixture$A, t(fixture$A))) {
    set.seed(1004)
    fit <- rsvd(A, k = 1, p = 0, q = 0)
    nla_test_expect_svd(A, fit, 1)
    nla_test_expect_exact(A, nla_test_reconstruct(fit))

    set.seed(1004)
    qb <- rqb(A, k = 1, p = 0, q = 0)
    nla_test_expect_qb(A, qb, 1)
    nla_test_expect_exact(A, qb$Q %*% qb$B)
  }
})

testthat::test_that("singular-vector requests preserve orientation and clamping", {
  fixture <- nla_test_fixture(34, 18, c(5, 2, 1))
  requests <- list(c(0L, 0L), c(1L, 0L), c(0L, 2L),
                   c(1L, 2L), c(-3L, -2L), c(99L, 99L))
  for (A in list(fixture$A, t(fixture$A))) {
    for (request in requests) {
      set.seed(1005)
      fit <- rsvd(A, k = 3, nu = request[1L], nv = request[2L], p = 4, q = 1)
      counts <- pmin(3L, pmax(0L, request))
      testthat::expect_equal(fit$d, fixture$d, tolerance = 1e-9)
      if (counts[1L] == 0L) {
        testthat::expect_null(fit$u)
      } else {
        testthat::expect_identical(dim(fit$u), c(nrow(A), counts[1L]))
        nla_test_expect_orthonormal(fit$u)
      }
      if (counts[2L] == 0L) {
        testthat::expect_null(fit$v)
      } else {
        testthat::expect_identical(dim(fit$v), c(ncol(A), counts[2L]))
        nla_test_expect_orthonormal(fit$v)
      }
    }

    # The default and an oversized k retain the historical dimension clamp.
    # Accuracy is assessed above in the intended low-rank regime.
    set.seed(1005)
    nla_test_expect_svd(A, rsvd(A), min(dim(A)))
    set.seed(1005)
    nla_test_expect_svd(A, rsvd(A, k = 100, p = 100), min(dim(A)))
  }
})

testthat::test_that("repeated singular values are compared through subspaces", {
  fixture <- nla_test_fixture(28, 19, c(7, 7, 7, 2, 1, 0.1))
  k <- 3L
  left <- fixture$u[, seq_len(k), drop = FALSE]
  right <- fixture$v[, seq_len(k), drop = FALSE]
  optimal <- sqrt(sum(fixture$d[-seq_len(k)]^2))
  for (seed in c(101L, 202L, 303L)) {
    set.seed(seed)
    fit <- rsvd(fixture$A, k = k, p = 3, q = 1)
    nla_test_expect_svd(fixture$A, fit, k)
    nla_test_expect_exact(left %*% nla_test_adjoint(left),
                          fit$u %*% nla_test_adjoint(fit$u))
    nla_test_expect_exact(right %*% nla_test_adjoint(right),
                          fit$v %*% nla_test_adjoint(fit$v))
    testthat::expect_equal(nla_test_norm(fixture$A - nla_test_reconstruct(fit)),
                           optimal, tolerance = 1e-9)
  }
})

testthat::test_that("all sampling distributions capture an exact low-rank range", {
  fixture <- nla_test_fixture(31, 21, c(5, 2, 0.4))
  for (distribution in c("normal", "unif", "rademacher")) {
    set.seed(1006)
    fit <- rsvd(fixture$A, k = 3, p = 5, q = 1, sdist = distribution)
    nla_test_expect_svd(fixture$A, fit, 3, distribution)
    nla_test_expect_exact(fixture$A, nla_test_reconstruct(fit), distribution)
  }
})

testthat::test_that("numeric data frames use the same matrix approximation API", {
  fixture <- nla_test_fixture(26, 16, c(4, 1, 0.2))
  for (A in list(fixture$A, t(fixture$A))) {
    input <- as.data.frame(A)
    set.seed(1007)
    fit <- rsvd(input, k = 3, p = 3, q = 1)
    nla_test_expect_svd(A, fit, 3)
    nla_test_expect_exact(A, nla_test_reconstruct(fit))
    set.seed(1007)
    qb <- rqb(input, k = 3, p = 3, q = 1)
    nla_test_expect_qb(A, qb, 6)
    nla_test_expect_exact(A, qb$Q %*% qb$B)
  }
})
