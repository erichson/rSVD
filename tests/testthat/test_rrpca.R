#devtools::use_package("testthat")
#library('testthat')

context("Randomized RPCA")

#Load rsvd library
library(rsvd)

#Set seed
set.seed(1234)

#Set p, q
p=5
q=1


#Accuray
atol_float64 <- 1e-8


#*************************************************************************************
# Test: Low-rank matrix with randomly corrupted entries
#*************************************************************************************
#Create real random test matrix of dimension m x n with target rank k
m = 50
n = 40
k = 5
L1 = matrix(rnorm(m*k), nrow = m, ncol = k) 
L2 = matrix(rnorm(n*k), nrow = k, ncol = n) 
L = L1 %*% L2
S = matrix(runif(m*n, -500, 500), nrow = m, ncol = n)
corruption_probability = 0.1
S = S * matrix(rbinom(m*n, size=1, prob=corruption_probability), nrow = m, ncol = n)
testMat = L + S


#*************************************************************************************
# Test 1: deterministic SVD
#*************************************************************************************
out <- rrpca(testMat, maxiter=100, tol=1e-09, trace=FALSE, rand=FALSE)
rerror <- norm(L - out$L, 'F') / norm(L, 'F')

testthat::test_that("Test 1: RPCA", {
  testthat::expect_equal(out$L, L, tolerance = 1e-02)
  testthat::expect_equal(out$S, S, tolerance = 1e-02)
  testthat::expect_lt(rerror, 0.1)
})

#*************************************************************************************
# Test 2: randomized SVD
#*************************************************************************************
out <- rrpca(testMat, p=p, q=q, maxiter=100,  tol=1e-08, trace=FALSE)
rerror <- norm(L - out$L, 'F') / norm(L, 'F')

testthat::test_that("Test 2: Randomized RPCA", {
  testthat::expect_equal(out$L, L, tolerance = 1e-02)
  testthat::expect_equal(out$S, S, tolerance = 1e-02)
  testthat::expect_lt(rerror, 0.1)
})

#*************************************************************************************
# Test: Low-rank matrix without sparse corruption
#*************************************************************************************
m = 50
n = 40
k = 5
L <- matrix(runif(m*k), m, k)
L <- L %*% t(L)
L <- L[,1:n]
S = matrix(0, nrow = m, ncol = n)
testMat = L + S

#*************************************************************************************
# Test 3: deterministic SVD
#*************************************************************************************
out <- rrpca(testMat, maxiter=100, tol=1e-08, trace=FALSE, rand=FALSE)
rerror <- norm(L - out$L, 'F') / norm(L, 'F')

testthat::test_that("Test 3: RPCA", {
  testthat::expect_equal(out$L, L, tolerance = 1e-02)
  testthat::expect_equal(out$S, S, tolerance = 1e-02)
  testthat::expect_lt(rerror, 1)
})

#*************************************************************************************
# Test 4: randomized SVD
#*************************************************************************************
out <- rrpca(testMat, p=p, q=q, maxiter=100,  tol=1e-08, trace=FALSE)
rerror <- norm(L - out$L, 'F') / norm(L, 'F')

testthat::test_that("Test 4: Randomized RPCA", {
  testthat::expect_equal(out$L, L, tolerance = 1e-02)
  testthat::expect_equal(out$S, S, tolerance = 1e-02)
  testthat::expect_lt(rerror, 1)
})
