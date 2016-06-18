#devtools::use_package("testthat", type = "Suggests")

context("Randomized Eigen")

#Load reigen library
library(rsvd)

#Set seed
set.seed(1234)

#Set p, q
  p=5
  q=1

#Accuray
atol_float64 <- 1e-8

#*************************************************************************************
# Test: reigen using real random test matrix
#*************************************************************************************
#Create real random test matrix of dimension m x n with target rank k
m = 50
n = 30
k = 10
testMat <- matrix(runif(m*k), m, k)
testMat <- testMat %*% t(testMat)
testMat <- testMat[,1:n]


#Deterministic eigen
  AtA = H(testMat)%*%testMat
  eigen.out <- eigen(AtA)


#Randomized eigen k=k
  reigen.out <- reigen(testMat, k=k)
  AtA.re = reigen.out$vectors %*% diag(reigen.out$values) %*% H(reigen.out$vectors)
  testthat::test_that("Test 2: Randomized eigen k=k", {
    testthat::expect_equal(eigen.out$values[1:k], reigen.out$values[1:k])
    testthat::expect_equal(AtA, AtA.re)
  })

#Randomized eigen k=k, p=0, q=0
  reigen.out <- reigen(testMat, k=k, p=0, q=0)
  AtA.re = reigen.out$vectors %*% diag(reigen.out$values) %*% H(reigen.out$vectors)
  testthat::test_that("Test 2: Randomized eigen k=k", {
    testthat::expect_equal(eigen.out$values[1:k], reigen.out$values[1:k])
    testthat::expect_equal(AtA, AtA.re)
  })


#*************************************************************************************
# Test : testMat.T
#*************************************************************************************
testMat = H(testMat)
AtA = H(testMat)%*%testMat

#Deterministic eigen
  eigen.out <- eigen(AtA)


#Randomized eigen k=k
reigen.out <- reigen(testMat, k=k)
AtA.re = reigen.out$vectors %*% diag(reigen.out$values) %*% H(reigen.out$vectors)
testthat::test_that("Test 2: Randomized eigen k=k", {
  testthat::expect_equal(eigen.out$values[1:k], reigen.out$values[1:k])
  testthat::expect_equal(AtA, AtA.re)
})

#Randomized eigen k=k, p=0, q=0
reigen.out <- reigen(testMat, k=k, p=0, q=0)
AtA.re = reigen.out$vectors %*% diag(reigen.out$values) %*% H(reigen.out$vectors)
testthat::test_that("Test 2: Randomized eigen k=k", {
  testthat::expect_equal(eigen.out$values[1:k], reigen.out$values[1:k])
  testthat::expect_equal(AtA, AtA.re)
})


#*************************************************************************************
# Test 3: complex random test matrix
#*************************************************************************************
testMat <- matrix(runif(m*k), m, k) + 1i* matrix(runif(m*k), m, k)
testMat <- testMat %*% H(testMat)
testMat <- testMat[,1:n]

AtA = H(testMat)%*%testMat

#Deterministic eigen
  eigen.out <- eigen(AtA)

#Randomized eigen k=k
  reigen.out <- reigen(testMat, k=k)
  AtA.re = reigen.out$vectors %*% diag(reigen.out$values) %*% H(reigen.out$vectors)
  testthat::test_that("Test 2: Randomized eigen k=k", {
    testthat::expect_equal(eigen.out$values[1:k], reigen.out$values[1:k])
    testthat::expect_equal(AtA, AtA.re)
  })


#Randomized eigen k=k, p=0, q=0
  reigen.out <- reigen(testMat, k=k, p=0, q=0)
  AtA.re = reigen.out$vectors %*% diag(reigen.out$values) %*% H(reigen.out$vectors)
  percent_error = 100*norm(AtA - AtA.re,'2')/norm(AtA,'2')
  testthat::test_that("Test 2: Randomized eigen k=k", {
    testthat::expect_equal(eigen.out$values[1:k], reigen.out$values[1:k])
    testthat::expect_equal(AtA, AtA.re)
  })


#*************************************************************************************
# Test: complex random test matrix transposed
#*************************************************************************************

  testMat <- H(testMat)
  AtA = H(testMat)%*%testMat


#Deterministic eigen
  eigen.out <- eigen(AtA)

#Randomized eigen k=k
  reigen.out <- reigen(testMat, k=k)
  AtA.re = reigen.out$vectors %*% diag(reigen.out$values) %*% H(reigen.out$vectors)
  testthat::test_that("Test 2: Randomized eigen k=k", {
    testthat::expect_equal(eigen.out$values[1:k], reigen.out$values[1:k])
    testthat::expect_equal(AtA, AtA.re)
  })


#Randomized eigen k=k, p=0, q=0
  reigen.out <- reigen(testMat, k=k, p=0, q=0)
  AtA.re = reigen.out$vectors %*% diag(reigen.out$values) %*% H(reigen.out$vectors)
  percent_error = 100*norm(AtA - AtA.re,'2')/norm(AtA,'2')
  testthat::test_that("Test 2: Randomized eigen k=k", {
    testthat::expect_equal(eigen.out$values[1:k], reigen.out$values[1:k])
    testthat::expect_equal(AtA, AtA.re)
  })
