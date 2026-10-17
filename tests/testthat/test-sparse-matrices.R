testthat::test_that("sparse C, R, and T representations preserve low-rank accuracy", {
  # Disjoint rank-one blocks give rank three and leave structural zero rows
  # and columns. Transposing also exercises the wide-input rsvd branch.
  A <- matrix(0, 36, 24)
  for (block in seq_len(3L)) {
    rows <- (block - 1L) * 8L + seq_len(8L)
    columns <- (block - 1L) * 6L + seq_len(6L)
    A[rows, columns] <- outer(seq_len(8L), c(1, -2, 3, -1, 2, -3)) / block
  }
  sparse <- Matrix::Matrix(A, sparse = TRUE)
  for (representation in c("CsparseMatrix", "RsparseMatrix", "TsparseMatrix")) {
    for (transpose in c(FALSE, TRUE)) {
      input <- methods::as(if (transpose) Matrix::t(sparse) else sparse,
                           representation)
      set.seed(2001)
      fit <- rsvd(input, k = 3, p = 4, q = 1)
      nla_test_expect_svd(input, fit, 3, representation)
      nla_test_expect_exact(input, nla_test_reconstruct(fit), representation)

      set.seed(2001)
      qb <- rqb(input, k = 3, p = 4, q = 1)
      nla_test_expect_qb(input, qb, 7, representation)
      nla_test_expect_exact(input, qb$Q %*% qb$B, representation)
    }
  }
})

testthat::test_that("symmetric and diagonal sparse inputs retain their meaning", {
  A <- matrix(0, 24, 24)
  for (block in seq_len(3L)) {
    index <- (block - 1L) * 6L + seq_len(6L)
    A[index, index] <- outer(seq_len(6L), seq_len(6L)) / block
  }
  inputs <- list(
    symmetric = Matrix::forceSymmetric(Matrix::Matrix(A, sparse = TRUE), uplo = "L"),
    diagonal = Matrix::Diagonal(x = c(6, 3, 1, rep(0, 21)))
  )
  for (input in inputs) {
    for (q in c(0L, 2L)) {
      set.seed(2002)
      fit <- rsvd(input, k = 3, p = 3, q = q)
      nla_test_expect_svd(input, fit, 3)
      nla_test_expect_exact(input, nla_test_reconstruct(fit))
      set.seed(2002)
      qb <- rqb(input, k = 3, p = 3, q = q)
      nla_test_expect_qb(input, qb, 6)
      nla_test_expect_exact(input, qb$Q %*% qb$B)
    }
  }
})

testthat::test_that("logical and pattern sparse inputs use numeric products", {
  # Each all-one rectangular block has rank one. Omitting x constructs a
  # pattern matrix, which has no numeric x slot for the implementation to use.
  rows <- integer()
  columns <- integer()
  for (block in seq_len(3L)) {
    rows <- c(rows, rep((block - 1L) * 8L + seq_len(6L), 3L))
    columns <- c(columns, rep((block - 1L) * 5L + seq_len(3L), each = 6L))
  }
  inputs <- list(
    logical = Matrix::sparseMatrix(i = rows, j = columns, x = TRUE, dims = c(30, 18)),
    pattern = Matrix::sparseMatrix(i = rows, j = columns, dims = c(30, 18))
  )
  for (input in inputs) {
    for (A in list(input, Matrix::t(input))) {
      set.seed(2003)
      fit <- rsvd(A, k = 3, p = 4, q = 2)
      nla_test_expect_svd(A, fit, 3)
      nla_test_expect_exact(A, nla_test_reconstruct(fit))
      set.seed(2003)
      qb <- rqb(A, k = 3, p = 4, q = 2)
      nla_test_expect_qb(A, qb, 7)
      nla_test_expect_exact(A, qb$Q %*% qb$B)
    }
  }
})

testthat::test_that("zero sparse inputs produce finite factors", {
  input <- Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(),
                                dims = c(14, 9))
  for (A in list(input, Matrix::t(input))) {
    for (q in c(0L, 2L)) {
      set.seed(2004)
      fit <- rsvd(A, k = 3, p = 2, q = q)
      nla_test_expect_svd(A, fit, 3)
      testthat::expect_equal(fit$d, rep(0, 3), tolerance = 1e-12)
      nla_test_expect_exact(A, nla_test_reconstruct(fit))
      set.seed(2004)
      qb <- rqb(A, k = 3, p = 2, q = q)
      nla_test_expect_qb(A, qb, 5)
      nla_test_expect_exact(A, qb$Q %*% qb$B)
    }
  }
})

testthat::test_that("one-row and one-column sparse inputs retain matrix outputs", {
  input <- Matrix::sparseMatrix(i = rep(1L, 4), j = c(1L, 3L, 5L, 7L),
                                x = c(2, -1, 3, 4), dims = c(1, 9))
  for (A in list(input, Matrix::t(input))) {
    set.seed(2005)
    fit <- rsvd(A, k = 1, p = 3, q = 2)
    nla_test_expect_svd(A, fit, 1)
    nla_test_expect_exact(A, nla_test_reconstruct(fit))
    set.seed(2005)
    qb <- rqb(A, k = 1, p = 3, q = 2)
    nla_test_expect_qb(A, qb, 1)
    nla_test_expect_exact(A, qb$Q %*% qb$B)
  }
})

testthat::test_that("deterministic QB preserves its dense QR behavior for sparse inputs", {
  # Matrix sparse QR cannot factor a wide matrix, and handles structurally
  # deficient matrices differently. This branch intentionally uses dense QR.
  input <- Matrix::sparseMatrix(i = c(2L, 4L, 7L), j = c(1L, 3L, 5L),
                                x = c(6, 3, 1), dims = c(14, 9))
  for (A in list(input, Matrix::t(input))) {
    fit <- rqb(A, k = 3, rand = FALSE)
    nla_test_expect_qb(A, fit, 3)
    nla_test_expect_exact(A, fit$Q %*% fit$B)
  }
})

testthat::test_that("sparse decaying-spectrum approximations meet the tail bound", {
  d <- 10^seq(0, -3, length.out = 24)
  A <- Matrix::Diagonal(x = d)
  k <- 4L
  optimal <- sqrt(sum(d[-seq_len(k)]^2))
  for (seed in c(17L, 81L, 407L)) {
    set.seed(seed)
    fit <- rsvd(A, k = k, p = 6, q = 2)
    nla_test_expect_svd(A, fit, k)
    residual <- nla_test_norm(as.matrix(A) - nla_test_reconstruct(fit))
    testthat::expect_gte(residual + 1e-10, optimal)
    testthat::expect_lt(residual, 1.1 * optimal)
  }
})
