#devtools::use_package("testthat")
#library('testthat')

context("Randomized SVD")

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
# Test: rSVD using real random test matrix
#*************************************************************************************
#Create real random test matrix of dimension m x n with target rank k
m = 50
n = 30
k = 10
testMat <- matrix(runif(m*k), m, k)
testMat <- testMat %*% t(testMat)
testMat <- testMat[,1:n]


#Deterministic SVD
  svd_out <- svd(testMat)

#Randomized SVD k=n
  rsvd_out <- rsvd(testMat)
  testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
  testthat::test_that("Test 1: Randomized SVD k=n", {
    testthat::expect_equal(svd_out$d, rsvd_out$d)
    testthat::expect_equal(testMat, testMat.re)
  })

#Randomized SVD k=k
  rsvd_out <- rsvd(testMat, k=k)
  testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
  testthat::test_that("Test 2: Randomized SVD k=k", {
    testthat::expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    testthat::expect_equal(testMat, testMat.re)
  })

#Randomized SVD k=k, p=0, q=0
  rsvd_out <- rsvd(testMat, k=k, p=0, q=0)
  testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
  testthat::test_that("Test 3: Randomized SVD k=k, p=0, q=0", {
    testthat::expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    testthat::expect_equal(testMat, testMat.re)
  })


#*************************************************************************************
# Test : testMat.T
#*************************************************************************************

testMat <- t(testMat)

#Deterministic SVD
  svd_out <- svd(testMat)

#Randomized SVD k=n
  rsvd_out <- rsvd(testMat)
  testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
  testthat::test_that("Test 4: Randomized SVD k=n", {
    testthat::expect_equal(svd_out$d, rsvd_out$d)
    testthat::expect_equal(testMat, testMat.re)
  })

#Randomized SVD k=k
  rsvd_out <- rsvd(testMat, k=k)
  testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
  testthat::test_that("Test 5: Randomized SVD k=k", {
    testthat::expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    testthat::expect_equal(testMat, testMat.re)
  })

#Randomized SVD k=k, p=0, q=0
  rsvd_out <- rsvd(testMat, k=k, p=0, q=0)
  testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
  testthat::test_that("Test 6: Randomized SVD k=k, p=0, q=0", {
    testthat::expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    testthat::expect_equal(testMat, testMat.re)
  })




#*************************************************************************************
# Test 3: complex random test matrix
#*************************************************************************************
  testMat <- matrix(runif(m*k), m, k) + 1i* matrix(runif(m*k), m, k)
  testMat <- testMat %*% H(testMat)
  testMat <- testMat[,1:n]


#Deterministic SVD
  svd_out <- svd(testMat)

#Randomized SVD k=n
  rsvd_out <- rsvd(testMat)
  testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% H(rsvd_out$v)
  testthat::test_that("Test 7: Randomized SVD k=n", {
    testthat::expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    testthat::expect_equal(testMat, testMat.re)
  })


#Randomized SVD k=k
  rsvd_out <- rsvd(testMat, k=k)
  testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% H(rsvd_out$v)
  testthat::test_that("Test 8: Randomized SVD k=k", {
    testthat::expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    testthat::expect_equal(testMat, testMat.re)
  })


#Randomized SVD k=k, p=0, q=0
  rsvd_out <- rsvd(testMat, k=k, p=0, q=0)
  testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% H(rsvd_out$v)
  percent_error = 100*norm(testMat - testMat.re,'2')/norm(testMat.re,'2')
  testthat::test_that("Test 9: Randomized SVD k=k, p=0, q=0", {
    testthat::expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    testthat::expect_lt(percent_error, 10)
  })


#*************************************************************************************
# Test: complex random test matrix transposed
#*************************************************************************************

  testMat <- H(testMat)

#Deterministic SVD
  svd_out <- svd(testMat)


#Randomized SVD k=n
  rsvd_out <- rsvd(testMat)
  testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% H(rsvd_out$v)
  testthat::test_that("Test 10: Randomized SVD k=n", {
    testthat::expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    testthat::expect_equal(testMat, testMat.re)
  })


#Randomized SVD k=k
  rsvd_out <- rsvd(testMat, k=k)
  testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% H(rsvd_out$v)
  testthat::test_that("Test 11: Randomized SVD k=k", {
    testthat::expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    testthat::expect_equal(testMat, testMat.re)
  })


#Randomized SVD k=k, p=0, q=0
  rsvd_out <- rsvd(testMat, k=k, p=0, q=0)
  testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% H(rsvd_out$v)
  percent_error = 100*norm(testMat - testMat.re,'2')/norm(testMat.re,'2')
  testthat::test_that("Test 12: Randomized SVD k=k, p=0, q=0", {
    testthat::expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    testthat::expect_lt(percent_error, 10)
  })


#*************************************************************************************
# Test: incorrect inputs
#*************************************************************************************

#Randomized SVD
  rsvd_out <- rsvd(testMat, k=k, nu=k+2, nv=k+2)
  testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% H(rsvd_out$v)
  testthat::test_that("Test 13: Randomized SVD", {
    testthat::expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    testthat::expect_equal(testMat, testMat.re)
    })

#Randomized SVD
  rsvd_out <- rsvd(testMat, k=k, p=100, nu=-5, nv=-3)
  testthat::test_that("Test 14: Randomized SVD", {
    testthat::expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
  })


#*************************************************************************************
# Test: random test matrices
#*************************************************************************************
  
testMat <- H(testMat)
  
#Deterministic SVD
  svd_out <- svd(testMat)

#Randomized SVD using uniform random test matrix
  rsvd_out <- rsvd(testMat, k = k, sdist = 'unif')
  testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% H(rsvd_out$v)
  testthat::test_that("Test 15: Randomized SVD using uniform random test matrix", {
    testthat::expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    testthat::expect_equal(testMat, testMat.re)
  })
  
  
#Randomized SVD using normal random test matrix
  rsvd_out <- rsvd(testMat, k = k, sdist = 'normal')
  testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% H(rsvd_out$v)
  testthat::test_that("Test 16: Randomized SVD using normal random test matrix", {
    testthat::expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    testthat::expect_equal(testMat, testMat.re)
  })
  
  
#Randomized SVD using rademacher random test matrix
  rsvd_out <- rsvd(testMat, k = k, sdist = 'rademacher')
  testMat.re = rsvd_out$u %*% diag(rsvd_out$d) %*% H(rsvd_out$v)
  testthat::test_that("Test 18: Randomized SVD using rademacher random test matrix", {
    testthat::expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    testthat::expect_equal(testMat, testMat.re)
  })

#*************************************************************************************
# Test: alternative matrix representations
#*************************************************************************************

library(Matrix)
sparseMat <- rsparsematrix(100, 200, 0.1)

set.seed(1234)
sparse_rsvd <- rsvd(sparseMat, k=k)
set.seed(1234)
ref_rsvd <- rsvd(as.matrix(sparseMat), k=k)

testthat::test_that("Test 19: Randomized SVD with sparse matrix input", {
    testthat::expect_equal(sparse_rsvd, ref_rsvd)
})
