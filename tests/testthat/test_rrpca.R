#devtools::use_package("testthat")

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
# Test : Low-rank matrix with randomly currpted entries
#*************************************************************************************
#Create real random test matrix of dimension m x n with target rank k
m = 50
n = 40
k = 5
L <- matrix(runif(m*k), m, k)
L <- L %*% t(L)
L <- L[,1:n]
S = rep(0, m*n)
S[sample(0:(m*n), 20)] = 10
S = matrix(S, m, n)
testMat = L + S


#*************************************************************************************
# Test 1: deterministic SVD
#*************************************************************************************
taurusRRPCA <- rrpca(testMat, k=1, maxiter=100,  tol=1e-08, trace=FALSE, rand=FALSE)
nrmse <- sqrt(sum((L - taurusRRPCA$L)**2) / sum(L**2)) * 100 # reconstruction error
testthat::test_that("Test 1: RPCA", {
  testthat::expect_equal(taurusRRPCA$L, L, tolerance = 1e-02)
  testthat::expect_equal(taurusRRPCA$S, S, tolerance = 1e-02)
  testthat::expect_lt(nrmse, 0.1)
})

#*************************************************************************************
# Test 2: randomized SVD
#*************************************************************************************
taurusRRPCA <- rrpca(testMat, k=1, p=p, q=q, maxiter=100,  tol=1e-08, trace=FALSE)
nrmse <- sqrt(sum((L - taurusRRPCA$L)**2) / sum(L**2)) * 100 # reconstruction error
testthat::test_that("Test 2: Randomized RPCA", {
  testthat::expect_equal(taurusRRPCA$L, L, tolerance = 1e-02)
  testthat::expect_equal(taurusRRPCA$S, S, tolerance = 1e-02)
  testthat::expect_lt(nrmse, 0.1)
})

#*************************************************************************************
# Test: Low-rank matrix with 0 patch
#*************************************************************************************
m = 50
n = 40
k = 5
L <- matrix(runif(m*k), m, k)
L <- L %*% t(L)
L <- L[,1:n]
S = matrix(0, nrow = m, ncol = n)
S[10:15,10:15] = 0
testMat = L + S

#*************************************************************************************
# Test 3: deterministtic SVD
#*************************************************************************************
taurusRRPCA <- rrpca(testMat, k=1, maxiter=100,  tol=1e-08, trace=FALSE, rand=FALSE)
nrmse <- sqrt(sum((L - taurusRRPCA$L)**2) / sum(L**2)) * 100 # reconstruction error
testthat::test_that("Test 3: RPCA", {
  testthat::expect_equal(taurusRRPCA$L, L, tolerance = 1e-02)
  testthat::expect_equal(taurusRRPCA$S, S, tolerance = 1e-02)
  testthat::expect_lt(nrmse, 1)
})

#*************************************************************************************
# Test 4: randomized SVD
#*************************************************************************************
taurusRRPCA <- rrpca(testMat, k=1, p=p, q=q, maxiter=100,  tol=1e-08, trace=FALSE)
nrmse <- sqrt(sum((L - taurusRRPCA$L)**2) / sum(L**2)) * 100 # reconstruction error
testthat::test_that("Test 4: Randomized RPCA", {
  testthat::expect_equal(taurusRRPCA$L, L, tolerance = 1e-02)
  testthat::expect_equal(taurusRRPCA$S, S, tolerance = 1e-02)
  testthat::expect_lt(nrmse, 1)
})
