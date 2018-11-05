#devtools::use_package("testthat")
#library('testthat')

context("Randomized CUR decomposition")

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


#CUR decomposition
  cur_out <- rcur(testMat)
  testMat.re = cur_out$C %*% cur_out$U %*% cur_out$R
  testthat::test_that("Test 1: CUR decomposition c=r=NULL", {
    testthat::expect_equal(testMat, testMat.re)
  })
  
  testMat.re = testMat[,cur_out$C.idx] %*% cur_out$U %*% testMat[cur_out$R.idx,]
  testthat::test_that("Test 2: CUR decomposition (column idx)  c=r=NULL", {
    testthat::expect_equal(testMat, testMat.re)
  })  
  
  
#CUR decomposition, k
  cur_out <- rcur(testMat, k=k)
  testMat.re = cur_out$C %*% cur_out$U %*% cur_out$R
  testthat::test_that("Test 3: CUR decomposition k=k", {
    testthat::expect_equal(testMat, testMat.re)
  })

#CUR decomposition, k
  cur_out <- rcur(H(testMat), k=k)
  testMat.re = cur_out$C %*% cur_out$U %*% cur_out$R
  testthat::test_that("Test 4: CUR decomposition k=k", {
    testthat::expect_equal(H(testMat), testMat.re)
  })
    
#*************************************************************************************
# Test: complex input matrix
#*************************************************************************************
  testMat <- matrix(runif(m*k), m, k) + 1i* matrix(runif(m*k), m, k)
  testMat <- testMat %*% H(testMat)
  testMat <- testMat[,1:n]


  #CUR decomposition
  cur_out <- rcur(testMat)
  testMat.re = cur_out$C %*% cur_out$U %*% cur_out$R
  testthat::test_that("Test 5: CUR decomposition c=r=NULL", {
    testthat::expect_equal(testMat, testMat.re)
  })
  
  testMat.re = testMat[,cur_out$C.idx] %*% cur_out$U %*% testMat[cur_out$R.idx,]
  testthat::test_that("Test 6: CUR decomposition (column idx)  c=r=NULL", {
    testthat::expect_equal(testMat, testMat.re)
  })  
  
  
  #CUR decomposition, k
  cur_out <- rcur(testMat, k=k)
  testMat.re = cur_out$C %*% cur_out$U %*% cur_out$R
  testthat::test_that("Test 7: CUR decomposition k=k", {
    testthat::expect_equal(testMat, testMat.re)
  })
  
  #CUR decomposition, k
  cur_out <- rcur(H(testMat), k=k)
  testMat.re = cur_out$C %*% cur_out$U %*% cur_out$R
  testthat::test_that("Test 8: CUR decomposition k=k", {
    testthat::expect_equal(H(testMat), testMat.re)
  })
  