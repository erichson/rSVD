#devtools::use_package("testthat")
#library('testthat')

context("Randomized ID")

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


#Interpolative decomposition, column
  id_out <- rid(testMat)
  testMat.re = id_out$C %*% id_out$Z
  testthat::test_that("Test 1: Interpolative decomposition (column) k=NULL", {
    testthat::expect_equal(testMat, testMat.re)
  })
  
  testMat.re = testMat[,id_out$idx] %*% id_out$Z
  testthat::test_that("Test 2: Interpolative decomposition (column idx) k=NULL", {
    testthat::expect_equal(testMat, testMat.re)
  })  

  
#Interpolative decomposition, row
  id_out <- rid(testMat, mode='row')
  testMat.re = id_out$Z %*% id_out$R
  testthat::test_that("Test 3: Interpolative decomposition (row) k=NULL", {
    testthat::expect_equal(testMat, testMat.re)
  })
  
  testMat.re =  id_out$Z %*% testMat[id_out$idx,]
  testthat::test_that("Test 4: Interpolative decomposition (row idx) k=NULL", {
    testthat::expect_equal(testMat, testMat.re)
  })   
  
  
#Interpolative decomposition, column, k
  id_out <- rid(testMat, k=k)
  testMat.re = id_out$C %*% id_out$Z
  testthat::test_that("Test 5: Interpolative decomposition (column) k=k", {
    testthat::expect_equal(testMat, testMat.re)
  })
  
  testMat.re = testMat[,id_out$idx] %*% id_out$Z
  testthat::test_that("Test 6: Interpolative decomposition (column idx) k=k", {
    testthat::expect_equal(testMat, testMat.re)
  })  
  
  
#Interpolative decomposition, row, k
  id_out <- rid(testMat, mode='row', k=k)
  testMat.re = id_out$Z %*% id_out$R
  testthat::test_that("Test 7: Interpolative decomposition (row) k=k", {
    testthat::expect_equal(testMat, testMat.re)
  })
  
  testMat.re =  id_out$Z %*% testMat[id_out$idx,]
  testthat::test_that("Test 8: Interpolative decomposition (row idx) k=k", {
    testthat::expect_equal(testMat, testMat.re)
  })     


#*************************************************************************************
# Test: complex input matrix
#*************************************************************************************
  testMat <- matrix(runif(m*k), m, k) + 1i* matrix(runif(m*k), m, k)
  testMat <- testMat %*% H(testMat)
  testMat <- testMat[,1:n]


#Interpolative decomposition, column, k
  id_out <- rid(testMat, k=k)
  testMat.re = id_out$C %*% id_out$Z
  testthat::test_that("Test 9: Interpolative decomposition (column) k=k", {
    testthat::expect_equal(testMat, testMat.re)
  })
  
  testMat.re = testMat[,id_out$idx] %*% id_out$Z
  testthat::test_that("Test 10: Interpolative decomposition (column idx) k=k", {
    testthat::expect_equal(testMat, testMat.re)
  })  
  
  
#Interpolative decomposition, row, k
  id_out <- rid(testMat, mode='row', k=k)
  testMat.re = id_out$Z %*% id_out$R
  testthat::test_that("Test 11: Interpolative decomposition (row) k=k", {
    testthat::expect_equal(testMat, testMat.re)
  })
  
  testMat.re =  id_out$Z %*% testMat[id_out$idx,]
  testthat::test_that("Test 12: Interpolative decomposition (row idx) k=k", {
    testthat::expect_equal(testMat, testMat.re)
  })     