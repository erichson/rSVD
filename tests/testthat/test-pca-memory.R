# Allocation checks guard the sparse-storage contract, not elapsed time.
# The numerical accuracy of these paths is checked against dense references in
# test-pca-reference.R. Warmups keep namespace loading and byte compilation out
# of the profiles. Rprofmem does not measure peak RSS or native allocations.
pca_largest_allocation <- function(run, threshold) {
  path <- tempfile("rsvd-pca-allocation-", fileext = ".out")
  on.exit(unlink(path), add = TRUE)
  run()
  utils::Rprofmem(path, threshold = threshold)
  on.exit(utils::Rprofmem(NULL), add = TRUE)
  run()
  utils::Rprofmem(NULL)
  records <- readLines(path, warn = FALSE)
  records <- records[grepl("^[0-9]+[[:space:]]*:", records)]
  sizes <- as.numeric(sub("^([0-9]+).*", "\\1", records))
  max(c(0, sizes))
}

testthat::test_that("sparse randomized PCA avoids input-sized dense allocations", {
  testthat::skip_if_not(capabilities("profmem"), "Rprofmem is unavailable")
  A <- Matrix::bdiag(lapply(seq_len(10L), function(j) {
    outer(sin(seq_len(80L) / (j + 1)), cos(seq_len(60L) / (j + 2)))
  }))
  for (X in list(A, Matrix::t(A))) {
    dense_bytes <- 8 * nrow(X) * ncol(X)
    for (scores in c(FALSE, TRUE)) {
      run <- function() {
        set.seed(746)
        rpca(X, k = 5, p = 3, q = 1, center = TRUE, scale = TRUE,
             retx = scores)
      }
      testthat::expect_lt(pca_largest_allocation(run, dense_bytes / 2),
                          dense_bytes)
    }
  }
})

testthat::test_that("sparse PCA prediction avoids an input-sized dense allocation", {
  testthat::skip_if_not(capabilities("profmem"), "Rprofmem is unavailable")
  A <- Matrix::bdiag(lapply(seq_len(10L), function(j) {
    outer(sin(seq_len(80L) / (j + 1)), cos(seq_len(60L) / (j + 2)))
  }))
  set.seed(747)
  fit <- rpca(A, k = 5, p = 3, q = 1)
  dense_bytes <- 8 * nrow(A) * ncol(A)
  run <- function() predict(fit, A)
  testthat::expect_lt(pca_largest_allocation(run, dense_bytes / 2), dense_bytes)
})
