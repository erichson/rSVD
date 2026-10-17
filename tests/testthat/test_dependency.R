testthat::test_that("sparse methods work without attaching Matrix", {
  testthat::skip_if("package:Matrix" %in% search(),
                    "Run in a fresh R session to check namespace-only use.")
  A <- Matrix::Diagonal(x = c(4, 1, rep(0, 18)))
  set.seed(700)
  fit <- rsvd(A, k = 2, p = 3, q = 1)
  testthat::expect_equal(fit$d, c(4, 1), tolerance = 1e-10)
  testthat::expect_false("package:Matrix" %in% search())
})
