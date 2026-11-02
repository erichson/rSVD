testthat::test_that("invalid sketches fail before producing incomplete factors", {
  A <- diag(8)
  # A negative p used to make the sketch narrower than k and returned NA in d.
  # A fractional q silently ran a different number of power iterations.
  for (fun in list(rsvd, rqb)) {
    for (args in list(list(p = -1), list(p = 0.5), list(p = NA_real_),
                      list(q = -1), list(q = 0.5), list(q = Inf),
                      list(q = c(1, 2)), list(sdist = 1L),
                      list(sdist = c("normal", "unif")))) {
      set.seed(2601)
      before <- .Random.seed
      testthat::expect_error(do.call(fun, c(list(A = A, k = 3), args)),
                             "must be")
      testthat::expect_identical(.Random.seed, before)
    }
  }
})

testthat::test_that("malformed ranks and vector counts do not draw random numbers", {
  A <- diag(8)
  for (invalid in list(numeric(), c(2, 3), NA_real_, Inf, "2", 1 + 1i)) {
    for (fun in list(rsvd, rqb)) {
      set.seed(2602)
      before <- .Random.seed
      testthat::expect_error(fun(A, k = invalid), "'k' must be")
      testthat::expect_identical(.Random.seed, before)
    }
    testthat::expect_error(rsvd(A, k = 3, nu = invalid), "'nu' must be")
    testthat::expect_error(rsvd(A, k = 3, nv = invalid), "'nv' must be")
  }
  testthat::expect_error(rqb(A, k = 3, rand = NA), "'rand' must be")
  testthat::expect_error(rqb(A, k = 3, rand = c(TRUE, FALSE)), "'rand' must be")
})

testthat::test_that("invalid matrix dimensions and types fail before sampling", {
  matrices <- list(matrix(numeric(), 0, 2), matrix(numeric(), 2, 0),
                    matrix(letters[1:4], 2, 2))
  for (A in matrices) {
    for (fun in list(rsvd, rqb)) {
      set.seed(2603)
      before <- .Random.seed
      testthat::expect_error(fun(A, k = 1), "'A' must")
      testthat::expect_identical(.Random.seed, before)
    }
  }
  testthat::expect_error(rsvd(array(1, c(2, 2, 2)), k = 1), "positive dimensions")
})

testthat::test_that("randomized products still reject non-finite data", {
  # QR checks the thin products, avoiding a redundant scan of every finite
  # input matrix. Invalid data can therefore consume random numbers.
  for (value in c(NA_real_, NaN, Inf, -Inf)) {
    A <- diag(6)
    A[6, 6] <- value
    for (input in list(A, A + 0i, Matrix::Matrix(A, sparse = TRUE))) {
      for (fun in list(rsvd, rqb)) {
        testthat::expect_error(fun(input, k = 1), "NA/NaN/Inf|infinite or missing")
      }
    }
  }
})

testthat::test_that("deterministic QB ignores randomized controls", {
  A <- matrix(seq_len(30), 6, 5)
  expected <- rqb(A, k = 2, rand = FALSE)
  set.seed(2604)
  before <- .Random.seed
  actual <- rqb(A, k = 2, p = NA_real_, q = "unused",
                sdist = "unused", rand = FALSE)
  testthat::expect_identical(actual, expected)
  testthat::expect_identical(.Random.seed, before)
})

testthat::test_that("rank clamping and legacy fractional QB ranks are preserved", {
  set.seed(2605)
  A <- matrix(rnorm(12 * 8), 12, 8)
  testthat::expect_equal(rqb(A, k = 2.5, rand = FALSE)$Q,
                         qr.Q(qr(A))[, 1:2, drop = FALSE], tolerance = 1e-12)
  set.seed(2606)
  expected <- rsvd(A, k = 8, p = 0, q = 0)
  set.seed(2606)
  actual <- rsvd(A, k = 100, p = 0, q = 0, nu = 100, nv = 100)
  testthat::expect_identical(actual, expected)
  set.seed(2606)
  values <- rsvd(A, k = 100, p = 0, q = 0, nu = -2, nv = -3)
  testthat::expect_equal(values$d, expected$d, tolerance = 1e-12)
  testthat::expect_null(values$u)
  testthat::expect_null(values$v)
})
