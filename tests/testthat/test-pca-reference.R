# Compare invariant numerical quantities: individual singular vectors may change
# sign, and a basis within a repeated singular subspace is not unique.
pca_test_projector <- function(fit) {
  fit$rotation %*% nla_test_adjoint(fit$rotation)
}

pca_test_reconstruction <- function(fit) {
  fit$x %*% nla_test_adjoint(fit$rotation)
}

# Closed-form PCA for a zero first processed column and a second column given
# by scores. Extreme-input tests must not use a dense mean as their oracle:
# dense accumulation can also lose accuracy or overflow on some platforms.
pca_test_rank_one_reference <- function(scores, center, scale) {
  variance <- sum(scores^2) / (length(scores) - 1L)
  structure(list(
    rotation = matrix(c(0, 1), ncol = 1,
                       dimnames = list(names(center), "PC1")),
    eigvals = variance, sdev = sqrt(variance), var = variance,
    center = center, scale = scale, x = matrix(scores, ncol = 1)
  ), class = "rpca")
}

pca_test_expect_equal <- function(actual, expected, input, tolerance = 1e-9) {
  testthat::expect_s3_class(actual, "rpca")
  testthat::expect_identical(names(actual), names(expected))
  testthat::expect_equal(actual$center, expected$center, tolerance = tolerance)
  testthat::expect_equal(actual$scale, expected$scale, tolerance = tolerance)
  testthat::expect_equal(actual$eigvals, expected$eigvals, tolerance = tolerance)
  testthat::expect_equal(actual$var, expected$var, tolerance = tolerance)
  testthat::expect_equal(pca_test_projector(actual), pca_test_projector(expected),
                         tolerance = tolerance)
  testthat::expect_identical(dimnames(actual$rotation), dimnames(expected$rotation))
  testthat::expect_true(is.matrix(actual$rotation))
  if (is.null(expected$x)) {
    testthat::expect_null(actual$x)
  } else {
    testthat::expect_true(is.matrix(actual$x))
    testthat::expect_identical(dimnames(actual$x), dimnames(expected$x))
    testthat::expect_equal(pca_test_reconstruction(actual),
                           pca_test_reconstruction(expected),
                           tolerance = tolerance)
  }
  # Evaluate each fit in its own basis before comparing reconstructions.
  actual_prediction <- predict(actual, input)
  expected_prediction <- predict(expected, as.matrix(input))
  testthat::expect_true(is.matrix(actual_prediction))
  testthat::expect_identical(dimnames(actual_prediction),
                             dimnames(expected_prediction))
  testthat::expect_equal(actual_prediction %*% nla_test_adjoint(actual$rotation),
                         expected_prediction %*% nla_test_adjoint(expected$rotation),
                         tolerance = tolerance)
}

testthat::test_that("PCA matches prcomp for centered and uncentered real data", {
  # prcomp leaves data uncentered when center=FALSE. Its full singular spectrum
  # supplies the total-energy denominator even when rank. truncates the vectors.
  A <- nla_test_fixture(30, 8, c(20, 11, 6, 3, 1, 0.5, 0.2, 0.1))$A
  A <- sweep(A, 2, seq(-3, 4), "+")
  colnames(A) <- paste0("variable", seq_len(ncol(A)))
  rownames(A) <- paste0("observation", seq_len(nrow(A)))
  k <- 3L
  for (center in c(FALSE, TRUE)) {
    for (scale in c(FALSE, TRUE)) {
      reference <- stats::prcomp(A, center = center, scale. = scale, rank. = k)
      fit <- rpca(A, k = k, center = center, scale = scale, rand = FALSE)
      testthat::expect_equal(fit$center, reference$center, tolerance = 1e-12)
      testthat::expect_equal(fit$scale, reference$scale, tolerance = 1e-12)
      testthat::expect_equal(fit$sdev, reference$sdev[seq_len(k)], tolerance = 1e-12)
      testthat::expect_equal(fit$eigvals, reference$sdev[seq_len(k)]^2,
                             tolerance = 1e-12)
      testthat::expect_equal(pca_test_projector(fit), pca_test_projector(reference),
                             tolerance = 1e-10)
      testthat::expect_equal(pca_test_reconstruction(fit),
                             pca_test_reconstruction(reference), tolerance = 1e-10)
      testthat::expect_equal(fit$var, sum(reference$sdev^2), tolerance = 1e-12)
      fractions <- fit$eigvals / fit$var
      expected <- reference$sdev[seq_len(k)]^2 / sum(reference$sdev^2)
      testthat::expect_equal(fractions, expected, tolerance = 1e-12)
      testthat::expect_equal(unname(summary(fit)["Proportion of variance", ]),
                             round(expected, 3))
      testthat::expect_lt(sum(fractions), 1)
    }
  }
})

testthat::test_that("center FALSE preserves PCA of already centered data", {
  A <- nla_test_fixture(30, 8, c(20, 11, 6, 3, 1, 0.5, 0.2, 0.1))$A
  A <- scale(A, center = TRUE, scale = FALSE)
  for (scale in c(FALSE, TRUE)) {
    centered <- rpca(A, k = 3, center = TRUE, scale = scale, rand = FALSE)
    uncentered <- rpca(A, k = 3, center = FALSE, scale = scale, rand = FALSE)
    testthat::expect_identical(uncentered$center, FALSE)
    testthat::expect_equal(uncentered$scale, centered$scale, tolerance = 1e-12)
    testthat::expect_equal(uncentered$eigvals, centered$eigvals, tolerance = 1e-12)
    testthat::expect_equal(uncentered$var, centered$var, tolerance = 1e-12)
    testthat::expect_equal(pca_test_projector(uncentered), pca_test_projector(centered),
                           tolerance = 1e-10)
    testthat::expect_equal(pca_test_reconstruction(uncentered),
                           pca_test_reconstruction(centered), tolerance = 1e-10)
  }
})

testthat::test_that("randomized PCA keeps total variance independent of retained rank", {
  A <- nla_test_fixture(48, 24, 10^seq(0, -3, length.out = 24))$A
  for (center in c(FALSE, TRUE)) {
    reference <- stats::prcomp(A, center = center, scale. = FALSE)
    for (k in c(3L, 6L)) {
      set.seed(8101)
      fit <- rpca(A, k = k, p = 5, q = 2, center = center, scale = FALSE)
      testthat::expect_length(fit$sdev, k)
      testthat::expect_equal(fit$var, sum(reference$sdev^2), tolerance = 1e-12)
      testthat::expect_lte(sum(fit$eigvals) / fit$var, 1 + 1e-12)
      processed <- scale(A, center = center, scale = FALSE)
      residual <- nla_test_norm(processed - pca_test_reconstruction(fit))
      optimal <- sqrt(sum(reference$sdev[-seq_len(k)]^2) * (nrow(A) - 1))
      testthat::expect_gte(residual + 1e-12, optimal)
      testthat::expect_lt(residual, 1.05 * optimal)
    }
  }
})

testthat::test_that("randomized PCA fitted scores retain the projected SVD convention", {
  # A tall matrix and q=0 make the range-finder residual visible. Fitted scores
  # are U D, whereas prediction projects the data onto the fitted loadings.
  A <- nla_test_fixture(64, 32, 10^seq(0, -2, length.out = 32))$A
  set.seed(8102)
  fit <- rpca(A, k = 4, p = 2, q = 0, center = TRUE, scale = FALSE)
  centered <- sweep(A, 2, colMeans(A), "-")
  projected <- predict(fit, A)
  testthat::expect_equal(unname(projected), unname(centered %*% fit$rotation),
                         tolerance = 1e-12)
  testthat::expect_gt(nla_test_norm(projected - fit$x), 1e-5)
  testthat::expect_lt(nla_test_norm(crossprod(fit$x, projected - fit$x)), 1e-10)
  testthat::expect_equal(unname(crossprod(fit$x)),
                         diag(fit$eigvals * (nrow(A) - 1)), tolerance = 1e-10)
})

testthat::test_that("sparse PCA agrees with dense PCA for tall and wide low-rank data", {
  A <- matrix(0, 48, 30)
  for (block in seq_len(3L)) {
    rows <- (block - 1L) * 12L + seq_len(10L)
    columns <- (block - 1L) * 8L + seq_len(6L)
    A[rows, columns] <- outer(seq_len(10L), c(1, -2, 3, -1, 2, -3)) / block
  }
  dimnames(A) <- list(paste0("row", seq_len(nrow(A))),
                      paste0("column", seq_len(ncol(A))))
  sparse <- Matrix::Matrix(A, sparse = TRUE)
  original <- sparse
  for (transpose in c(FALSE, TRUE)) {
    input <- if (transpose) Matrix::t(sparse) else sparse
    for (center in c(FALSE, TRUE)) {
      for (scale in c(FALSE, TRUE)) {
        set.seed(8103)
        expected <- rpca(as.matrix(input), k = 3, p = 4, q = 1,
                          center = center, scale = scale)
        for (representation in c("CsparseMatrix", "RsparseMatrix", "TsparseMatrix")) {
          represented <- methods::as(input, representation)
          set.seed(8103)
          actual <- rpca(represented, k = 3, p = 4, q = 1,
                          center = center, scale = scale)
          pca_test_expect_equal(actual, expected, represented)
          set.seed(8103)
          compact <- rpca(represented, k = 3, p = 4, q = 1,
                           center = center, scale = scale, retx = FALSE)
          testthat::expect_null(compact$x)
          testthat::expect_equal(compact$eigvals, actual$eigvals, tolerance = 1e-10)
          testthat::expect_equal(pca_test_projector(compact), pca_test_projector(actual),
                                 tolerance = 1e-10)
        }
      }
    }
  }
  testthat::expect_identical(sparse, original)
})

testthat::test_that("sparse PCA respects diagonal symmetric logical and pattern matrices", {
  rows <- rep(c(seq_len(5L), 9L:13L, 17L:21L), 3L)
  columns <- rep(c(1L, 7L, 13L, 2L, 8L, 14L, 3L, 9L, 15L), each = 5L)
  symmetric <- matrix(0, 24, 24)
  for (block in seq_len(3L)) {
    index <- (block - 1L) * 7L + seq_len(5L)
    symmetric[index, index] <- outer(seq_len(5L), seq_len(5L)) / block
  }
  inputs <- list(
    diagonal = Matrix::Diagonal(x = c(8, 4, 2, rep(0, 21))),
    symmetric = Matrix::forceSymmetric(Matrix::Matrix(symmetric, sparse = TRUE),
                                       uplo = "L"),
    logical = Matrix::sparseMatrix(i = rows, j = columns, x = TRUE,
                                   dims = c(30, 18)),
    pattern = Matrix::sparseMatrix(i = rows, j = columns, dims = c(30, 18))
  )
  for (input in inputs) {
    for (center in c(FALSE, TRUE)) {
      set.seed(8104)
      expected <- rpca(as.matrix(input), k = 3, p = 4, q = 2,
                        center = center, scale = TRUE)
      set.seed(8104)
      actual <- rpca(input, k = 3, p = 4, q = 2, center = center, scale = TRUE)
      pca_test_expect_equal(actual, expected, input)
    }
  }
})

testthat::test_that("sparse PCA retains accuracy when the neglected tail is nonzero", {
  A <- matrix(0, 48, 30)
  A[cbind(seq_len(24L), seq_len(24L))] <- 10^seq(0, -3, length.out = 24)
  for (input in list(Matrix::Matrix(A, sparse = TRUE),
                     Matrix::Matrix(t(A), sparse = TRUE))) {
    for (center in c(FALSE, TRUE)) {
      processed <- scale(as.matrix(input), center = center, scale = FALSE)
      singular_values <- svd(processed, nu = 0, nv = 0)$d
      optimal <- sqrt(sum(singular_values[-seq_len(4L)]^2))
      for (q in c(0L, 2L)) {
        set.seed(8109)
        expected <- rpca(as.matrix(input), k = 4, p = 4, q = q,
                          center = center, scale = FALSE)
        set.seed(8109)
        actual <- rpca(input, k = 4, p = 4, q = q,
                        center = center, scale = FALSE)
        pca_test_expect_equal(actual, expected, input)
        residual <- nla_test_norm(processed - pca_test_reconstruction(actual))
        testthat::expect_gte(residual + 1e-12, optimal)
        if (q == 2L) testthat::expect_lt(residual, 1.05 * optimal)
      }
    }
  }
})

testthat::test_that("sparse PCA omits missing rows and retains names", {
  A <- outer(c(seq_len(12L), rep(0, 12)), c(0, 2, -1, 0, 3, 0, -2, 0))
  dimnames(A) <- list(paste0("row", seq_len(nrow(A))),
                      paste0("column", seq_len(ncol(A))))
  sparse <- Matrix::Matrix(A, sparse = TRUE)
  sparse[3, 2] <- NA_real_
  sparse[20, 4] <- NaN
  original <- sparse
  clean <- A[-c(3, 20), , drop = FALSE]
  for (randomized in c(FALSE, TRUE)) {
    set.seed(8105)
    expected <- rpca(clean, k = 1, p = 2, q = 1, rand = randomized)
    set.seed(8105)
    testthat::expect_warning(
      actual <- rpca(sparse, k = 1, p = 2, q = 1, rand = randomized),
      "Missing values")
    pca_test_expect_equal(actual, expected, Matrix::Matrix(clean, sparse = TRUE))
    testthat::expect_identical(rownames(actual$x), rownames(clean))
  }
  testthat::expect_identical(sparse, original)
})

testthat::test_that("sparse PCA prediction uses the training transformation", {
  A <- outer(c(seq_len(12L), rep(0, 12)), c(0, 2, -1, 0, 3, 0, -2, 0))
  colnames(A) <- paste0("variable", seq_len(ncol(A)))
  newdata <- A[c(2, 4, 18), , drop = FALSE]
  newdata[3, c(2, 5)] <- c(-1, 4)
  rownames(newdata) <- c("new1", "new2", "new3")
  input <- Matrix::Matrix(A, sparse = TRUE)
  for (center in c(FALSE, TRUE)) {
    for (scale in c(FALSE, TRUE)) {
      set.seed(8106)
      fit <- rpca(input, k = 1, center = center, scale = scale)
      processed <- newdata
      if (center) processed <- sweep(processed, 2, fit$center, "-")
      if (scale) processed <- sweep(processed, 2, fit$scale, "/")
      expected <- processed %*% fit$rotation
      for (represented in list(newdata, as.data.frame(newdata),
                                Matrix::Matrix(newdata, sparse = TRUE))) {
        testthat::expect_equal(predict(fit, represented), expected, tolerance = 1e-10)
      }
    }
  }
})

testthat::test_that("sparse PCA handles zero and nearly constant columns stably", {
  zero <- Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(),
                               dims = c(20, 8))
  for (center in c(FALSE, TRUE)) {
    for (scale in c(FALSE, TRUE)) {
      set.seed(8107)
      fit <- rpca(zero, k = 2, center = center, scale = scale)
      testthat::expect_equal(fit$eigvals, c(0, 0), tolerance = 1e-12)
      testthat::expect_equal(fit$var, 0, tolerance = 1e-12)
      testthat::expect_true(all(is.finite(fit$x)))
      testthat::expect_true(all(is.finite(fit$rotation)))
      if (scale) testthat::expect_equal(unname(fit$scale), rep(1, 8))
    }
  }

  # Computing variance as E[x^2] - E[x]^2 loses this small variation. Explicit
  # centering is the stable reference; the sparse products may lose a few more
  # digits when they subtract the mean contribution after multiplication.
  offset <- seq(-1, 1, length.out = 24)
  A <- cbind(zero = 0, constant = 3, near = 1e6 + 1e-2 * offset,
             trend = offset^2, sparse = c(seq_len(6L), rep(0, 18)))
  input <- Matrix::Matrix(A, sparse = TRUE)
  centered <- sweep(A, 2, colMeans(A), "-")
  expected_scale <- sqrt(colSums(centered^2) / (nrow(A) - 1))
  expected_scale[expected_scale < 1e-8] <- 1
  set.seed(8108)
  expected <- rpca(A, k = 3, center = TRUE, scale = TRUE)
  set.seed(8108)
  actual <- rpca(input, k = 3, center = TRUE, scale = TRUE)
  testthat::expect_equal(actual$scale, expected_scale, tolerance = 1e-10)
  pca_test_expect_equal(actual, expected, input, tolerance = 1e-6)
})

testthat::test_that("complex PCA total variance uses the processed modulus energy", {
  A <- nla_test_fixture(18, 8, c(12, 6, 3, 1, 0.5, 0.2, 0.1, 0.05),
                         complex = TRUE)$A
  A <- sweep(A, 2, seq_len(ncol(A)) * (1 + 2i), "+")
  for (center in c(FALSE, TRUE)) {
    for (scale in c(FALSE, TRUE)) {
      processed <- if (center) sweep(A, 2, colMeans(A), "-") else A
      if (scale) {
        scales <- sqrt(colSums(Mod(processed)^2) / (nrow(A) - 1))
        processed <- sweep(processed, 2, scales, "/")
      }
      reference <- svd(processed, nu = 0, nv = 2)
      fit <- rpca(A, k = 2, center = center, scale = scale, rand = FALSE)
      testthat::expect_equal(fit$eigvals, reference$d[1:2]^2 / (nrow(A) - 1),
                             tolerance = 1e-12)
      testthat::expect_equal(fit$var, sum(Mod(processed)^2) / (nrow(A) - 1),
                             tolerance = 1e-12)
      testthat::expect_equal(pca_test_projector(fit),
                             reference$v %*% nla_test_adjoint(reference$v),
                             tolerance = 1e-10)
    }
  }
})

testthat::test_that("sparse PCA respects implicit unit diagonals", {
  upper <- methods::as(Matrix::bandSparse(
    18, k = 1, diagonals = list(seq(0.5, 1.5, length.out = 17))),
    "triangularMatrix")
  # The diagonal consists of implicit ones, with no diagonal entries in x.
  upper@diag <- "U"
  testthat::expect_true(methods::validObject(upper))
  testthat::expect_equal(diag(as.matrix(upper)), rep(1, 18))
  for (input in list(upper, Matrix::t(upper))) {
    for (center in c(FALSE, TRUE)) {
      set.seed(8110)
      expected <- rpca(as.matrix(input), k = 3, p = 4, q = 2,
                        center = center, scale = TRUE)
      set.seed(8110)
      actual <- rpca(input, k = 3, p = 4, q = 2,
                      center = center, scale = TRUE)
      pca_test_expect_equal(actual, expected, input)
    }
  }
})

testthat::test_that("sparse PCA combines duplicate triplets before column moments", {
  A <- matrix(0, 18, 12)
  A[1:6, 1:4] <- outer(seq_len(6L), c(1, -2, 3, 1))
  A[9:14, 7:10] <- outer(seq_len(6L) / 2, c(2, 1, -1, 3))
  entries <- which(A != 0, arr.ind = TRUE)
  input <- Matrix::sparseMatrix(
    i = c(rep(entries[, 1], each = 2), 18L, 18L),
    j = c(rep(entries[, 2], each = 2), 12L, 12L),
    x = c(rep(A[entries], each = 2) * rep(c(0.75, 0.25), nrow(entries)), 1, -1),
    dims = dim(A), repr = "T")
  testthat::expect_s4_class(input, "TsparseMatrix")
  testthat::expect_gt(length(input@x), sum(A != 0))
  testthat::expect_equal(as.matrix(input), A)
  original <- input
  for (represented in list(input, Matrix::t(input))) {
    for (center in c(FALSE, TRUE)) {
      set.seed(8111)
      expected <- rpca(as.matrix(represented), k = 2, p = 3, q = 1,
                        center = center, scale = TRUE)
      set.seed(8111)
      actual <- rpca(represented, k = 2, p = 3, q = 1,
                      center = center, scale = TRUE)
      pca_test_expect_equal(actual, expected, represented)
    }
  }
  testthat::expect_identical(input, original)
})

testthat::test_that("sparse predictions retain NaN replacement and NA propagation", {
  A <- outer(seq_len(12L), c(1, -2, 3, 1, 2, -1))
  colnames(A) <- paste0("variable", seq_len(ncol(A)))
  newdata <- matrix(c(5, 6, 7, 8, 0, NaN, 0, 0, 1, 2, NaN, 4,
                      0, 0, 0, NA, 0, 0, 2, 0, 1, 0, 0, 2), 4, 6)
  dimnames(newdata) <- list(paste0("new", seq_len(nrow(newdata))), colnames(A))
  sparse <- Matrix::Matrix(newdata, sparse = TRUE)
  original <- sparse
  for (center in c(FALSE, TRUE)) {
    for (scale in c(FALSE, TRUE)) {
      set.seed(8112)
      fit <- rpca(A, k = 1, center = center, scale = scale)
      expected <- predict(fit, newdata)
      actual <- predict(fit, sparse)
      testthat::expect_equal(actual, expected, tolerance = 1e-10)
      testthat::expect_identical(is.na(actual), is.na(expected))
      testthat::expect_identical(is.nan(actual), is.nan(expected))
      if (scale) {
        testthat::expect_true(all(is.finite(actual[1:3, , drop = FALSE])))
      }
    }
  }
  testthat::expect_identical(sparse, original)
})

testthat::test_that("sparse PCA predicts a single observation without dropping dimensions", {
  A <- outer(seq_len(12L), c(1, -2, 3, 1, 2, -1))
  colnames(A) <- paste0("variable", seq_len(ncol(A)))
  for (center in c(FALSE, TRUE)) {
    for (scale in c(FALSE, TRUE)) {
      set.seed(8113)
      fit <- rpca(A, k = 1, center = center, scale = scale)
      for (values in list(c(0, 2, -1, 0, 3, 0), rep(0, ncol(A)))) {
        newdata <- matrix(values, nrow = 1,
                           dimnames = list("new", colnames(A)))
        actual <- predict(fit, Matrix::Matrix(newdata, sparse = TRUE))
        testthat::expect_identical(dim(actual), c(1L, 1L))
        testthat::expect_equal(actual, predict(fit, newdata), tolerance = 1e-10)
      }
    }
  }
})

testthat::test_that("complex PCA fits accept real sparse observations for prediction", {
  A <- nla_test_fixture(18, 8, c(12, 6, 3, 1, 0.5, 0.2, 0.1, 0.05),
                         complex = TRUE)$A
  A <- sweep(A, 2, seq_len(ncol(A)) * (1 + 2i), "+")
  colnames(A) <- paste0("variable", seq_len(ncol(A)))
  newdata <- Re(A[1:4, , drop = FALSE])
  newdata[abs(newdata) < 3] <- 0
  rownames(newdata) <- paste0("new", seq_len(nrow(newdata)))
  for (center in c(FALSE, TRUE)) {
    for (scale in c(FALSE, TRUE)) {
      fit <- rpca(A, k = 2, center = center, scale = scale, rand = FALSE)
      actual <- predict(fit, Matrix::Matrix(newdata, sparse = TRUE))
      testthat::expect_true(is.complex(actual))
      testthat::expect_equal(actual, predict(fit, newdata), tolerance = 1e-12)
    }
  }
})

testthat::test_that("omitting k retains deterministic PCA for sparse input", {
  A <- nla_test_fixture(18, 8, c(12, 6, 3, 1, 0.5, 0.2, 0.1, 0.05))$A
  for (input in list(Matrix::Matrix(A, sparse = TRUE),
                     Matrix::Matrix(t(A), sparse = TRUE))) {
    for (retx in c(FALSE, TRUE)) {
      expected <- rpca(as.matrix(input), retx = retx)
      set.seed(8114)
      before <- .Random.seed
      actual <- rpca(input, retx = retx)
      testthat::expect_identical(.Random.seed, before)
      testthat::expect_length(actual$eigvals, min(dim(input)))
      pca_test_expect_equal(actual, expected, input)
    }
  }
})

testthat::test_that("fully stored large constant columns have zero centered energy", {
  # A power of two has an exactly representable mean in both dense and sparse
  # accumulation. The absent implicit zeros must not contribute 0 * Inf.
  A <- cbind(constant = rep(2^700, 12), varying = seq_len(12L))
  input <- Matrix::Matrix(A, sparse = TRUE)
  for (scale in c(FALSE, TRUE)) {
    set.seed(8115)
    expected <- rpca(A, k = 1, center = TRUE, scale = scale)
    set.seed(8115)
    actual <- rpca(input, k = 1, center = TRUE, scale = scale)
    pca_test_expect_equal(actual, expected, input)
    if (scale) testthat::expect_identical(unname(actual$scale[1]), 1)
  }
})

testthat::test_that("fully stored constant columns use accurate column means", {
  # Sequential double-precision summation of this decimal constant can shift
  # its mean by an ulp, creating a spurious centered column of size ~1e184.
  A <- cbind(constant = rep(1e200, 12), varying = seq_len(12L))
  input <- Matrix::Matrix(A, sparse = TRUE)
  expected <- pca_test_rank_one_reference(
    scores = seq_len(12L) - 6.5,
    center = c(constant = 1e200, varying = 6.5), scale = FALSE)
  set.seed(8117)
  actual <- rpca(input, k = 1, center = TRUE, scale = FALSE)
  pca_test_expect_equal(actual, expected, input)
})

testthat::test_that("sparse column means remain finite when stored sums overflow", {
  A <- cbind(large = c(1e308, 1e308, 1e308, 0), ordinary = c(-1, 1, -1, 1))
  input <- Matrix::Matrix(A, sparse = TRUE)
  # The large column has a finite mean but its sum of squared deviations
  # overflows. Under the established scaling rule its processed column is zero.
  expected <- pca_test_rank_one_reference(
    scores = c(-1, 1, -1, 1) / sqrt(4 / 3),
    center = c(large = (3 / 4) * 1e308, ordinary = 0),
    scale = c(large = Inf, ordinary = sqrt(4 / 3)))
  set.seed(8118)
  actual <- rpca(input, k = 1, center = TRUE, scale = TRUE)
  testthat::expect_true(all(is.finite(actual$center)))
  pca_test_expect_equal(actual, expected, input)
})

testthat::test_that("PCA total energy follows the matrix after overflowing scale factors", {
  A <- cbind(large = 1e200 * c(1, -1, 2, -2, 0, 0, 0, 0),
             ordinary = c(1, -2, 3, -2, 1, -2, 3, -2), zero = rep(0, 8))
  # Preserve the established scaling convention: a sum-of-squares overflow
  # produces an infinite scale factor and a zero processed column.
  for (center in c(FALSE, TRUE)) {
    processed <- if (center) sweep(A, 2, colMeans(A), "-") else A
    scales <- sqrt(colSums(processed^2) / (nrow(A) - 1))
    scales[scales < 1e-8] <- 1
    processed <- sweep(processed, 2, scales, "/")
    expected_energy <- sum(processed^2) / (nrow(A) - 1)
    expected_eigenvalue <- svd(processed, nu = 0, nv = 0)$d[1]^2 / (nrow(A) - 1)
    for (input in list(A, Matrix::Matrix(A, sparse = TRUE))) {
      set.seed(8116)
      actual <- rpca(input, k = 1, center = center, scale = TRUE)
      testthat::expect_equal(actual$var, expected_energy, tolerance = 1e-12)
      testthat::expect_equal(actual$eigvals, expected_eigenvalue, tolerance = 1e-12)
      testthat::expect_equal(pca_test_reconstruction(actual), processed,
                             tolerance = 1e-10)
    }
  }
})
