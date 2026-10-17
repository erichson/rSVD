sparse_product_test_operands <- function() {
  Q <- diag(8)[, c(1L, 4L), drop = FALSE]
  dimnames(Q) <- list(paste0("row", seq_len(8)), c("first", "fourth"))
  A <- matrix(0, 8, 5, dimnames = list(rownames(Q), paste0("col", seq_len(5))))
  A[cbind(c(1L, 2L, 4L, 6L, 8L), seq_len(5))] <- c(1, -2, 3, 4, 5)
  list(Q = Q, A = A)
}

testthat::test_that("crossproducts accept row-compressed right operands", {
  fixture <- sparse_product_test_operands()
  sparse <- Matrix::Matrix(fixture$A, sparse = TRUE)
  inputs <- list(numeric = sparse, logical = sparse != 0,
                 pattern = methods::as(sparse != 0, "nMatrix"))
  dense_Q <- matrix(seq_len(16), 8, 2, dimnames = dimnames(fixture$Q))
  bases <- list(fixture$Q, fixture$Q[, 1L, drop = FALSE], dense_Q)
  for (input in inputs) {
    A <- methods::as(input, "RsparseMatrix")
    original <- A
    for (Q in bases) {
      # Coordinate bases and fully populated arrays exercise both sparse
      # and dense product dispatch with the same mathematical reference.
      expected <- t(Q) %*% as.matrix(A)
      actual <- rsvd:::crossprod_help(Q, A)
      testthat::expect_equal(as.matrix(actual), expected, tolerance = 1e-12)
      testthat::expect_identical(dim(actual), dim(expected))
    }
    testthat::expect_identical(A, original)
  }
})

testthat::test_that("crossproducts accept row-compressed left and both operands", {
  fixture <- sparse_product_test_operands()
  Q <- methods::as(Matrix::Matrix(fixture$Q, sparse = TRUE), "RsparseMatrix")
  original_Q <- Q
  sparse <- Matrix::Matrix(fixture$A, sparse = TRUE)
  dense <- fixture$A + 1
  inputs <- list(fixture$A, dense, sparse,
                 methods::as(sparse, "RsparseMatrix"),
                 methods::as(sparse, "TsparseMatrix"))
  for (A in inputs) {
    original_A <- A
    expected <- t(fixture$Q) %*% as.matrix(A)
    actual <- rsvd:::crossprod_help(Q, A)
    testthat::expect_equal(as.matrix(actual), expected, tolerance = 1e-12)
    testthat::expect_identical(A, original_A)
  }
  testthat::expect_identical(Q, original_Q)
})

testthat::test_that("crossproducts retain duplicate-triplet values and names", {
  fixture <- sparse_product_test_operands()
  # The first two triplets add; the last two cancel. A storage conversion
  # must preserve this mathematical matrix, without modifying the triplets.
  A <- Matrix::sparseMatrix(i = c(1L, 1L, 4L, 8L, 8L),
                            j = c(1L, 1L, 3L, 5L, 5L),
                            x = c(0.25, 0.75, -2, 3, -3),
                            dims = c(8L, 5L), repr = "T",
                            dimnames = dimnames(fixture$A))
  original <- A
  dense <- matrix(0, 8, 5, dimnames = dimnames(fixture$A))
  dense[1L, 1L] <- 1
  dense[4L, 3L] <- -2
  testthat::expect_equal(as.matrix(rsvd:::crossprod_help(fixture$Q, A)),
                         t(fixture$Q) %*% dense, tolerance = 1e-12)
  testthat::expect_equal(as.matrix(rsvd:::crossprod_help(A, fixture$Q)),
                         t(dense) %*% fixture$Q, tolerance = 1e-12)
  testthat::expect_identical(A, original)
})

testthat::test_that("row-compressed SVD handles rank-one structural zeros", {
  dense <- outer(c(3, rep(0, 15)), c(1, -2, 0, 3, rep(0, 6)))
  dimnames(dense) <- list(paste0("row", seq_len(nrow(dense))),
                          paste0("col", seq_len(ncol(dense))))
  for (matrix in list(dense, t(dense))) {
    A <- methods::as(Matrix::Matrix(matrix, sparse = TRUE), "RsparseMatrix")
    original <- A
    for (q in c(0L, 1L)) {
      set.seed(2301)
      fit <- rsvd(A, k = 1, p = 0, q = q)
      nla_test_expect_svd(A, fit, 1)
      testthat::expect_equal(fit$d, sqrt(126), tolerance = 1e-12)
      nla_test_expect_exact(matrix, nla_test_reconstruct(fit))
    }
    testthat::expect_identical(A, original)
  }
})

testthat::test_that("row-compressed QB handles structurally sparse bases", {
  dense <- outer(c(3, rep(0, 15)), c(1, -2, 0, 3, rep(0, 6)))
  dimnames(dense) <- list(paste0("row", seq_len(nrow(dense))),
                          paste0("col", seq_len(ncol(dense))))
  for (matrix in list(dense, t(dense))) {
    A <- methods::as(Matrix::Matrix(matrix, sparse = TRUE), "RsparseMatrix")
    original <- A
    for (randomized in c(FALSE, TRUE)) {
      set.seed(2302)
      fit <- rqb(A, k = 1, p = 0, q = 1, rand = randomized)
      nla_test_expect_qb(A, fit, 1)
      nla_test_expect_exact(matrix, fit$Q %*% fit$B)
    }
    testthat::expect_identical(A, original)
  }
})

testthat::test_that("crossproducts preserve complex conjugation and dimensions", {
  A <- matrix(c(1 + 2i, 0, -3i, 2 - 1i, -1, 4 + 2i), 3, 2,
               dimnames = list(paste0("row", 1:3), c("a", "b")))
  B <- matrix(c(2 - 1i, 3i, 1, 1 + 1i, -2, 4 - 1i), 3, 2,
               dimnames = list(rownames(A), c("c", "d")))
  original_A <- A
  original_B <- B
  for (columns in list(1L, 1:2)) {
    testthat::expect_equal(rsvd:::crossprod_help(A[, columns, drop = FALSE], B),
                           Conj(t(A[, columns, drop = FALSE])) %*% B,
                           tolerance = 1e-12)
  }
  testthat::expect_identical(A, original_A)
  testthat::expect_identical(B, original_B)
})

testthat::test_that("crossproducts preserve custom non-Matrix dispatch", {
  check_dispatch <- function() {
    # Defining methods in a local environment leaves the global generic alone.
    where <- new.env(parent = globalenv())
    where$crossprod <- Matrix::crossprod
    class_name <- "rsvd_test_product_matrix"
    generator <- methods::setClass(class_name, slots = c(data = "matrix"),
                                    where = where)
    on.exit({
      methods::removeMethod("crossprod", c(x = class_name, y = "matrix"),
                             where = where)
      methods::removeMethod("crossprod", c(x = "matrix", y = class_name),
                             where = where)
      methods::removeClass(class_name, where = where)
    })
    seen <- character()
    methods::setMethod("crossprod", c(x = class_name, y = "matrix"),
      function(x, y) {
        seen <<- c(seen, "left")
        base::crossprod(x@data, y)
      }, where = where)
    methods::setMethod("crossprod", c(x = "matrix", y = class_name),
      function(x, y) {
        seen <<- c(seen, "right")
        base::crossprod(x, y@data)
      }, where = where)
    A <- matrix(seq_len(8), 4, 2)
    wrapped <- generator(data = A)
    original <- wrapped
    testthat::expect_equal(rsvd:::crossprod_help(wrapped, A), t(A) %*% A)
    testthat::expect_equal(rsvd:::crossprod_help(A, wrapped), t(A) %*% A)
    testthat::expect_identical(seen, c("left", "right"))
    testthat::expect_identical(wrapped, original)
  }
  check_dispatch()
})
