#devtools::use_package("testthat")
#library('testthat')

context("Randomized QB decomposition")

#Load rsvd library
library(rsvd)

#Set seed
set.seed(1234)

#Accuray
atol_float64 <- 1e-8

#*************************************************************************************
# Test: real input matrix
#*************************************************************************************
#Create real random test matrix of dimension m x n with target rank k
m = 50
n = 30
k = 10
testMat <- matrix(runif(m*k), m, k)
testMat <- testMat %*% t(testMat)
testMat <- testMat[,1:n]


#Randomized QB decomposition, standard
  rqb_out <- rqb(testMat)
  testMat.re = rqb_out$Q %*% rqb_out$B
  testthat::test_that("Test 1: Randomized QB decomposition, standard", {
    testthat::expect_equal(testMat, testMat.re)
  })
  
  
#Randomized QB decomposition, low-rank
  rqb_out <- rqb(testMat, k=k)
  testMat.re = rqb_out$Q %*% rqb_out$B
  testthat::test_that("Test 2: Randomized QB decomposition, k=k", {
    testthat::expect_equal(testMat, testMat.re)
  })

#Randomized QB decomposition, low-rank
  rqb_out <- rqb(testMat, k=k, p=0, q=0)
  testMat.re = rqb_out$Q %*% rqb_out$B
  testthat::test_that("Test 3: Randomized QB decomposition, k=k, p=0, q=0", {
    testthat::expect_equal(testMat, testMat.re)
  })

#*************************************************************************************
# Test: complex input matrix
#*************************************************************************************
  testMat <- matrix(runif(m*k), m, k) + 1i* matrix(runif(m*k), m, k)
  testMat <- testMat %*% H(testMat)
  testMat <- testMat[,1:n]


#Randomized QB decomposition, standard
  rqb_out <- rqb(testMat)
  testMat.re = rqb_out$Q %*% rqb_out$B
  testthat::test_that("Test 4: Randomized QB decomposition, standard", {
    testthat::expect_equal(testMat, testMat.re)
  })
  
  
#Randomized QB decomposition, low-rank
  rqb_out <- rqb(testMat, k=k)
  testMat.re = rqb_out$Q %*% rqb_out$B
  testthat::test_that("Test 5: Randomized QB decomposition, k=k", {
    testthat::expect_equal(testMat, testMat.re)
  })
  
#Randomized QB decomposition, low-rank
  rqb_out <- rqb(testMat, k=k, p=0, q=0)
  testMat.re = rqb_out$Q %*% rqb_out$B
  testthat::test_that("Test 6: Randomized QB decomposition, k=k, p=0, q=0", {
    testthat::expect_equal(testMat, testMat.re)
  })
  
  
#*************************************************************************************
# Test: complex input matrix transposed
#*************************************************************************************
testMat <- H(testMat)
  
  
#Randomized QB decomposition, standard
  rqb_out <- rqb(testMat)
  testMat.re = rqb_out$Q %*% rqb_out$B
  testthat::test_that("Test 4: Randomized QB decomposition, standard", {
    testthat::expect_equal(testMat, testMat.re)
  })
  
  
#Randomized QB decomposition, low-rank
  rqb_out <- rqb(testMat, k=k)
  testMat.re = rqb_out$Q %*% rqb_out$B
  testthat::test_that("Test 5: Randomized QB decomposition, k=k", {
    testthat::expect_equal(testMat, testMat.re)
  })
  
#Randomized QB decomposition, low-rank
  rqb_out <- rqb(testMat, k=k, p=0, q=0)
  testMat.re = rqb_out$Q %*% rqb_out$B
  testthat::test_that("Test 6: Randomized QB decomposition, k=k, p=0, q=0", {
    testthat::expect_equal(testMat, testMat.re)
  })  