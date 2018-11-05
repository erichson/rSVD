#devtools::use_package("testthat")
#library('testthat')

context("Randomized PCA")

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
# Test: Randomized PCA - center = TRUE, scale. = TRUE
#*************************************************************************************
#Create real random test matrix of dimension m x n with target rank k
m = 50
n = 30
k = 10
testMat <- matrix(runif(m*k), m, k)
testMat <- testMat %*% t(testMat)
testMat <- testMat[,1:n]


#Deterministic PCA
  pca_out <- prcomp(testMat, center = TRUE, scale. = TRUE)

#Randomized PCA
  rpca_out <- rpca(testMat, k=k)
  testthat::test_that("Test 1: standard", {
    testthat::expect_equal(pca_out$sdev[1:k], rpca_out$sdev[1:k])
    testthat::expect_equal(sum(diag(1,k,k) - H(rpca_out$rotation)%*%rpca_out$rotation), 0 )
  })

#Randomized PCA
  rpca_out <- rpca(testMat, k=k-5, p=5, q=0)
  testthat::test_that("Test 2: standard 2", {
    testthat::expect_equal(pca_out$sdev[1:(k-5)], rpca_out$sdev[1:(k-5)])
    testthat::expect_equal(sum(diag(1,k-5,k-5) - H(rpca_out$rotation)%*%rpca_out$rotation), 0 )
  })

#Randomized PCA
  rpca_out <- rpca(testMat, k=k)
  cum_var = cumsum(rpca_out$sdev**2 / rpca_out$var)
  testthat::test_that("Test 3: cumsum", {
    testthat::expect_equal(cum_var[k], 1)
  })


#*************************************************************************************
# Test: Randomized PCA Complex - center = TRUE, scale. = TRUE
#*************************************************************************************
testMat <- matrix(runif(m*k), m, k) + 1i* matrix(runif(m*k), m, k)

#Deterministic PCA
  pca_out <- prcomp(testMat, center = TRUE, scale. = TRUE)

#Randomized PCA
  rpca_out <- rpca(testMat, k=k)
  testthat::test_that("Test 4: Complex", {
    testthat::expect_equal(pca_out$sdev[1:k], rpca_out$sdev[1:k])
    testthat::expect_equal(abs(sum(diag(1,k,k) - H(rpca_out$rotation)%*%rpca_out$rotation)), 0 )
  })

#Randomized PCA
  rpca_out <- rpca(testMat, k=k, p=5, q=0)
  testthat::test_that("Test 5: Complex", {
    testthat::expect_equal(pca_out$sdev[1:(k)], rpca_out$sdev[1:(k)])
    testthat::expect_equal(abs(sum(diag(1,k,k) - H(rpca_out$rotation)%*%rpca_out$rotation)), 0 )
  })

#Randomized PCA
  rpca_out <- rpca(testMat, k=k)
  cum_var = cumsum(rpca_out$sdev**2 / rpca_out$var)
  testthat::test_that("Test 6: Complex", {
    testthat::expect_equal(cum_var[k], 1)
  })


#*************************************************************************************
# Test: Randomized PCA - center = T, scale. = F
#*************************************************************************************

#Deterministic PCA
  pca_out <- prcomp(testMat, center = TRUE, scale. = FALSE)

#Randomized PCA
  rpca_out <- rpca(testMat, k=k, center = TRUE, scale = FALSE)
  testthat::test_that("Test 7: Complex", {
    testthat::expect_equal(pca_out$sdev[1:k], rpca_out$sdev[1:k])
    testthat::expect_equal(abs(sum(diag(1,k,k) - H(rpca_out$rotation)%*%rpca_out$rotation)), 0 )
  })

#Randomized PCA
  rpca_out <- rpca(testMat, k=k, p=5, q=0, center = TRUE, scale = FALSE)
  testthat::test_that("Test 8: Complex", {
    testthat::expect_equal(pca_out$sdev[1:(k)], rpca_out$sdev[1:(k)])
    testthat::expect_equal(abs(sum(diag(1,k,k) - H(rpca_out$rotation)%*%rpca_out$rotation)), 0 )
  })

#Randomized PCA
  rpca_out <- rpca(testMat, k=k, center = TRUE, scale = FALSE)
  cum_var = cumsum(rpca_out$sdev**2 / rpca_out$var)
  testthat::test_that("Test 9: Complex", {
    testthat::expect_equal(cum_var[k], 1)
  })

#*************************************************************************************
# Test: Randomized PCA - Predict complex
#*************************************************************************************

#Deterministic PCA
  pca_out <- prcomp(testMat, center = TRUE, scale. = FALSE)

#Randomized PCA
  rpca_out <- rpca(testMat, k=k, retx=TRUE, center = TRUE, scale = FALSE)

  testthat::test_that("Test 10: rPCA Predict 1", {
    testthat::expect_equal(rpca_out$x %*% H(rpca_out$rotation), pca_out$x %*% H(pca_out$rotation))
  })

#Randomized PCA
  rpca_out <- rpca(testMat, k=k, retx=TRUE, center = TRUE, scale = FALSE)
  preds <- predict(rpca_out, newdata = testMat)
  testthat::test_that("Test 11: rPCA Predict 2", {
    testthat::expect_equal(preds %*%  H(rpca_out$rotation), rpca_out$x %*% H(rpca_out$rotation) )
  })

  
#*************************************************************************************
# Test: Randomized PCA - Predict real
#*************************************************************************************
testMat <- matrix(runif(m*k), m, k)

#Deterministic PCA
pca_out <- prcomp(testMat, center = FALSE, scale. = FALSE)

#Randomized PCA
rpca_out <- rpca(testMat, k=k, retx=TRUE, center = FALSE, scale = FALSE)
testthat::test_that("Test 12: rPCA Predict 3", {
  testthat::expect_equal(rpca_out$x %*% H(rpca_out$rotation), pca_out$x %*% H(pca_out$rotation))
})

#Randomized PCA
rpca_out <- rpca(testMat, k=k, retx=FALSE, center = FALSE, scale = FALSE)
preds <- predict(rpca_out, newdata = testMat)
testthat::test_that("Test 13: rPCA Predict 5", {
  testthat::expect_equal(preds %*%  H(rpca_out$rotation), pca_out$x %*% H(pca_out$rotation) )
})

#Deterministic PCA
pca_out <- prcomp(testMat, center = TRUE, scale. = FALSE)

#Randomized PCA
rpca_out <- rpca(testMat, k=k, retx=TRUE, center = TRUE, scale = FALSE)
testthat::test_that("Test 14: rPCA Predict 6", {
  testthat::expect_equal(rpca_out$x %*% H(rpca_out$rotation), pca_out$x %*% H(pca_out$rotation))
})

#Randomized PCA
rpca_out <- rpca(testMat, k=k, retx=FALSE, center = TRUE, scale = FALSE)
preds <- predict(rpca_out, newdata = testMat)
testthat::test_that("Test 15: rPCA Predict 7", {
  testthat::expect_equal(preds %*%  H(rpca_out$rotation), pca_out$x %*% H(pca_out$rotation) )
})


#Deterministic PCA
pca_out <- prcomp(testMat, center = TRUE, scale. = TRUE)

#Randomized PCA
rpca_out <- rpca(testMat, k=k, retx=TRUE, center = TRUE, scale = TRUE)
testthat::test_that("Test 16: rPCA Predict 8", {
  testthat::expect_equal(rpca_out$x %*% H(rpca_out$rotation), pca_out$x %*% H(pca_out$rotation))
})

#Randomized PCA
rpca_out <- rpca(testMat, k=k, retx=FALSE, center = TRUE, scale = TRUE)
preds <- predict(rpca_out, newdata = testMat)
testthat::test_that("Test 17: rPCA Predict 9", {
  testthat::expect_equal(preds %*%  H(rpca_out$rotation), pca_out$x %*% H(pca_out$rotation) )
})


#*************************************************************************************
# Test: Randomized PCA - Reconstruction error - real matrix
#*************************************************************************************
testMat <- matrix(runif(m*k), m, k)

#Randomized PCA
  rpca_out <- rpca(testMat, k=k, retx=TRUE, center = FALSE, scale = FALSE)
  Re <- t(t(rpca_out$x %*% t(rpca_out$rotation)) + rpca_out$center)
  testthat::test_that("Test 18: rPCA reconstruction 1", {
    testthat::expect_equal(Re, testMat )
  })

#Deterministic PCA
  pca_out <- prcomp(testMat, center = TRUE, scale. = FALSE)
  Re <- t(t(pca_out$x %*% H(pca_out$rotation)) + pca_out$center)

  rpca_out <- rpca(testMat, k=k, retx=TRUE, center = TRUE, scale = FALSE)
  Re2 <- t(t(rpca_out$x %*% H(rpca_out$rotation)) + pca_out$center)

  testthat::test_that("Test 19: rPCA reconstruction 2", {
    testthat::expect_equal(Re, Re2  )
  })


#Deterministic PCA
  pca_out <- prcomp(testMat, center = TRUE, scale. = TRUE)
  Re <- t(t(pca_out$x %*% H(pca_out$rotation)) + pca_out$center)

  rpca_out <- rpca(testMat, k=k, retx=TRUE, center = TRUE, scale = TRUE)
  Re2 <- t(t(rpca_out$x %*% H(rpca_out$rotation)) + pca_out$center)

  testthat::test_that("Test 20: rPCA reconstruction 3", {
    testthat::expect_equal(Re, Re2  )
  })

#*************************************************************************************
# Test: Randomized PCA - Reconstruction error
#*************************************************************************************
testMat <- matrix(runif(m*k), m, k) + 1i * matrix(runif(m*k), m, k)

#Deterministic PCA
  pca_out <- prcomp(testMat, center = FALSE, scale. = FALSE)
  Re <- pca_out$x %*% H(pca_out$rotation)

#Randomized PCA
  rpca_out <- rpca(testMat, k=k, retx=TRUE, center = FALSE, scale = FALSE)
  Re2 <- rpca_out$x %*% H(rpca_out$rotation)
  testthat::test_that("Test 21: rPCA reconstruction 4", {
    testthat::expect_equal(Re, Re2 )
  })


#*************************************************************************************
# Test: Randomized PCA - Transposed input matrix
#*************************************************************************************
testMat <- H(testMat)

#Deterministic PCA
  pca_out <- prcomp(testMat, center = TRUE, scale. = TRUE)

#Randomized PCA
  rpca_out <- rpca(testMat, k=k)
  testthat::test_that("Test 22: Randomized SVD k=n", {
    testthat::expect_equal(pca_out$sdev[1:k], rpca_out$sdev[1:k])
    testthat::expect_equal(H(pca_out$rotation)%*%pca_out$rotation, H(rpca_out$rotation)%*%rpca_out$rotation)
  })

  
#*************************************************************************************
# Test: Randomized PCA - Transposed input matrix
#*************************************************************************************
testMat <- matrix(runif(m*k), m, k) + 1i* matrix(runif(m*k), m, k)
  
testMat <- H(testMat)
  
#Deterministic PCA
  pca_out <- prcomp(testMat, center = TRUE, scale. = TRUE)
  
#Randomized PCA
  rpca_out <- rpca(testMat, k=k)
  testthat::test_that("Test 23: Randomized SVD k=n", {
    testthat::expect_equal(pca_out$sdev[1:k], rpca_out$sdev[1:k])
    testthat::expect_equal(H(pca_out$rotation)%*%pca_out$rotation, H(rpca_out$rotation)%*%rpca_out$rotation)
  })  