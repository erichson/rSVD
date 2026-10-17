testthat::test_that("robust PCA permits a zero-rank step on nonzero input", {
  # This rank-one sparse input is entirely assigned to S by soft-thresholding.
  # Every remaining singular value is below the nuclear-norm threshold.
  A <- matrix(c(1, 0, 0, 0), 2, 2)
  for (randomized in c(FALSE, TRUE)) {
    set.seed(730)
    fit <- rrpca(A, rand = randomized, maxiter = 1)
    testthat::expect_equal(fit$L, matrix(0, 2, 2), tolerance = 1e-12)
    testthat::expect_equal(fit$S, A, tolerance = 1e-12)
    testthat::expect_equal(fit$err, 0, tolerance = 1e-12)
  }
})

testthat::test_that("robust PCA omits rows with missing values before finding dimensions", {
  A <- outer(seq_len(8), seq_len(5))
  A[3, 2] <- NA_real_
  A[6, 4] <- NaN
  clean <- A[-c(3, 6), , drop = FALSE]
  for (randomized in c(FALSE, TRUE)) {
    set.seed(731)
    expected <- rrpca(clean, rand = randomized, maxiter = 5)
    set.seed(731)
    testthat::expect_warning(
      actual <- rrpca(A, rand = randomized, maxiter = 5), "Missing values")
    testthat::expect_identical(dim(actual$L), dim(clean))
    testthat::expect_equal(actual, expected, tolerance = 1e-12)
  }
  testthat::expect_warning(
    testthat::expect_error(rrpca(matrix(NA_real_, 3, 2)), "positive dimensions"),
    "Missing values")
  testthat::expect_error(rrpca(matrix(numeric(), 0, 2)), "positive dimensions")
})

testthat::test_that("ID and CUR accept numeric data frames and retain selected names", {
  set.seed(732)
  A <- matrix(rnorm(24 * 3), 24, 3) %*% matrix(rnorm(3 * 14), 3, 14)
  dimnames(A) <- list(paste0("row", seq_len(nrow(A))),
                      paste0("col", seq_len(ncol(A))))
  for (X in list(A, t(A))) {
    frame <- as.data.frame(X)
    for (randomized in c(FALSE, TRUE)) {
      for (mode in c("column", "row")) {
        fit <- rid(frame, k = 3, p = 2, q = 1,
                   mode = mode, rand = randomized)
        if (mode == "column") {
          testthat::expect_true(is.numeric(fit$C))
          testthat::expect_equal(fit$C, X[, fit$idx, drop = FALSE])
          reconstructed <- fit$C %*% fit$Z
        } else {
          testthat::expect_true(is.numeric(fit$R))
          testthat::expect_equal(fit$R, X[fit$idx, , drop = FALSE])
          reconstructed <- fit$Z %*% fit$R
        }
        testthat::expect_lt(nla_test_norm(X - reconstructed) / nla_test_norm(X), 1e-10)
      }
      fit <- rcur(frame, k = 3, p = 2, q = 1, rand = randomized)
      testthat::expect_true(is.numeric(fit$C))
      testthat::expect_true(is.numeric(fit$R))
      testthat::expect_equal(fit$C, X[, fit$C.idx, drop = FALSE])
      testthat::expect_equal(fit$R, X[fit$R.idx, , drop = FALSE])
      testthat::expect_lt(nla_test_norm(X - fit$C %*% fit$U %*% fit$R) /
                           nla_test_norm(X), 1e-10)
    }
  }
})

testthat::test_that("ID and CUR preserve sparse low-rank inputs and index-only results", {
  A <- matrix(0, 30, 20)
  A[1:8, 1:5] <- outer(seq_len(8), seq_len(5))
  A[11:18, 8:12] <- outer(seq_len(8) + 1, seq_len(5) + 2)
  A[21:28, 15:19] <- outer(seq_len(8) + 3, seq_len(5) + 1)
  sparse <- Matrix::Matrix(A, sparse = TRUE)
  original <- sparse
  for (randomized in c(FALSE, TRUE)) {
    for (mode in c("column", "row")) {
      set.seed(733)
      fit <- rid(sparse, k = 3, p = 2, q = 1,
                 mode = mode, rand = randomized)
      set.seed(733)
      index_fit <- rid(sparse, k = 3, p = 2, q = 1,
                       mode = mode, rand = randomized, idx_only = TRUE)
      testthat::expect_null(index_fit$C)
      testthat::expect_null(index_fit$R)
      testthat::expect_identical(index_fit$idx, fit$idx)
      testthat::expect_equal(index_fit$Z, fit$Z, tolerance = 1e-12)
      if (mode == "column") {
        testthat::expect_equal(fit$C, A[, fit$idx, drop = FALSE])
        reconstructed <- fit$C %*% fit$Z
      } else {
        testthat::expect_equal(fit$R, A[fit$idx, , drop = FALSE])
        reconstructed <- fit$Z %*% fit$R
      }
      testthat::expect_lt(nla_test_norm(A - reconstructed) / nla_test_norm(A), 1e-10)
    }
    set.seed(734)
    fit <- rcur(sparse, k = 3, p = 2, q = 1, rand = randomized)
    set.seed(734)
    index_fit <- rcur(sparse, k = 3, p = 2, q = 1, rand = randomized,
                      idx_only = TRUE)
    testthat::expect_null(index_fit$C)
    testthat::expect_null(index_fit$R)
    testthat::expect_identical(index_fit$C.idx, fit$C.idx)
    testthat::expect_identical(index_fit$R.idx, fit$R.idx)
    testthat::expect_equal(index_fit$U, fit$U, tolerance = 1e-12)
    testthat::expect_equal(fit$C, A[, fit$C.idx, drop = FALSE])
    testthat::expect_equal(fit$R, A[fit$R.idx, , drop = FALSE])
    testthat::expect_lt(nla_test_norm(A - fit$C %*% fit$U %*% fit$R) /
                         nla_test_norm(A), 1e-10)
  }
  testthat::expect_identical(sparse, original)
})

testthat::test_that("sparse ID and CUR omit missing rows once without losing dimensions", {
  finite <- outer(c(seq_len(5), rep(0, 7)), c(0, 1, 2, 0, 4, 0, 0, 3))
  dimnames(finite) <- list(paste0("row", seq_len(nrow(finite))),
                           paste0("col", seq_len(ncol(finite))))
  sparse <- Matrix::Matrix(finite, sparse = TRUE)
  sparse[2, 2] <- NA_real_
  sparse[10, 3] <- NaN
  original <- sparse
  cases <- list(
    list(input = sparse, clean = finite[-c(2, 10), , drop = FALSE]),
    list(input = sparse[1:2, , drop = FALSE], clean = finite[1, , drop = FALSE])
  )
  with_omission_warning <- function(expr) {
    messages <- character()
    value <- withCallingHandlers(expr, warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    testthat::expect_identical(messages, "Missing values are omitted: na.omit(A).")
    value
  }
  for (case in cases) {
    clean <- case$clean
    clean_sparse <- Matrix::Matrix(clean, sparse = TRUE)
    for (randomized in c(FALSE, TRUE)) {
      for (mode in c("column", "row")) {
        set.seed(739)
        expected <- rid(clean_sparse, k = 1, p = 2, q = 1,
                        mode = mode, rand = randomized)
        set.seed(739)
        fit <- with_omission_warning(
          rid(case$input, k = 1, p = 2, q = 1, mode = mode, rand = randomized))
        testthat::expect_equal(fit, expected, tolerance = 1e-12)
        if (mode == "column") {
          testthat::expect_identical(dim(fit$C), c(nrow(clean), 1L))
          testthat::expect_equal(fit$C, clean[, fit$idx, drop = FALSE])
          reconstructed <- fit$C %*% fit$Z
        } else {
          testthat::expect_identical(dim(fit$R), c(1L, ncol(clean)))
          testthat::expect_equal(fit$R, clean[fit$idx, , drop = FALSE])
          reconstructed <- fit$Z %*% fit$R
        }
        nla_test_expect_exact(clean, reconstructed)
      }
      set.seed(740)
      expected <- rcur(clean_sparse, k = 1, p = 2, q = 1, rand = randomized)
      set.seed(740)
      fit <- with_omission_warning(
        rcur(case$input, k = 1, p = 2, q = 1, rand = randomized))
      testthat::expect_equal(fit, expected, tolerance = 1e-12)
      testthat::expect_identical(dim(fit$C), c(nrow(clean), 1L))
      testthat::expect_identical(dim(fit$U), c(1L, 1L))
      testthat::expect_identical(dim(fit$R), c(1L, ncol(clean)))
      testthat::expect_equal(fit$C, clean[, fit$C.idx, drop = FALSE])
      testthat::expect_equal(fit$R, clean[fit$R.idx, , drop = FALSE])
      nla_test_expect_exact(clean, fit$C %*% fit$U %*% fit$R)
    }
  }
  testthat::expect_identical(sparse, original)
})

testthat::test_that("rank-one complex ID and CUR retain matrix shapes", {
  A <- outer(seq_len(12) + 1i, seq_len(7) - 2i)
  for (X in list(A, Conj(t(A)))) {
    for (randomized in c(FALSE, TRUE)) {
      for (mode in c("column", "row")) {
        set.seed(735)
        fit <- rid(X, k = 1, p = 2, q = 1, mode = mode, rand = randomized)
        if (mode == "column") {
          testthat::expect_identical(dim(fit$C), c(nrow(X), 1L))
          testthat::expect_identical(dim(fit$Z), c(1L, ncol(X)))
          reconstructed <- fit$C %*% fit$Z
        } else {
          testthat::expect_identical(dim(fit$R), c(1L, ncol(X)))
          testthat::expect_identical(dim(fit$Z), c(nrow(X), 1L))
          reconstructed <- fit$Z %*% fit$R
        }
        testthat::expect_lt(nla_test_norm(X - reconstructed) / nla_test_norm(X), 1e-10)
      }
      set.seed(736)
      fit <- rcur(X, k = 1, p = 2, q = 1, rand = randomized)
      testthat::expect_identical(dim(fit$C), c(nrow(X), 1L))
      testthat::expect_identical(dim(fit$U), c(1L, 1L))
      testthat::expect_identical(dim(fit$R), c(1L, ncol(X)))
      testthat::expect_lt(nla_test_norm(X - fit$C %*% fit$U %*% fit$R) /
                           nla_test_norm(X), 1e-10)
    }
  }
})

testthat::test_that("PCA without scores preserves eigenvalues and the rotation subspace", {
  set.seed(737)
  A <- matrix(rnorm(32 * 3), 32, 3) %*% matrix(rnorm(3 * 18), 3, 18)
  for (X in list(A, t(A))) {
    for (randomized in c(FALSE, TRUE)) {
      set.seed(738)
      with_scores <- rpca(X, k = 3, p = 2, q = 1, rand = randomized)
      set.seed(738)
      without_scores <- rpca(X, k = 3, p = 2, q = 1, rand = randomized,
                              retx = FALSE)
      testthat::expect_null(without_scores$x)
      testthat::expect_identical(names(without_scores), names(with_scores))
      testthat::expect_equal(without_scores$eigvals, with_scores$eigvals,
                             tolerance = 1e-10)
      testthat::expect_equal(tcrossprod(without_scores$rotation),
                             tcrossprod(with_scores$rotation), tolerance = 1e-10)
    }
  }
})
