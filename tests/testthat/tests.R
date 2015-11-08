#devtools::use_testthat()
set.seed(1234)
#Create real random test matrix of dimension m x n with target rank k
  m = 50
  n = 30
  k = 10
  testMat <- matrix(runif(m*k), m, k)
  testMat <- testMat %*% t(testMat)
  testMat <- testMat[,1:n]

#Set p, q
  p=5
  q=2

#Accuray
atol_float64 <- 1e-8

#*************************************************************************************
# Test 1
#*************************************************************************************

#Deterministic SVD
  svd_out <- svd(testMat)
  #Ak = svd_out$u %*% diag(svd_out$d) %*% t(svd_out$v)
  #print( sum((Ak-test)**2) )

#Randomized SVD k=n
  rsvd_out <- rsvd(testMat, method='standard')
  Ak = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
  test_that("Randomized SVD k=n", {
    expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    expect_equal(testMat, Ak)
  })

#Fast randomized SVD k=n
  rsvd_out <- rsvd(testMat, method='fast')
  Ak = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
  test_that("Fast Randomized SVD k=n", {
    expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    expect_equal(testMat, Ak)
  })

#Randomized SVD k=k
  rsvd_out <- rsvd(testMat, k=k, method='standard')
  Ak = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
  test_that("Randomized SVD k=k", {
    expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    expect_equal(testMat, Ak)
  })

#Fast randomized SVD k=k
  rsvd_out <- rsvd(testMat, k=k, method='fast')
  Ak = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
  test_that("Fast Randomized SVD k=k", {
    expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    expect_equal(testMat, Ak)
  })

#Randomized SVD k=k, p=0, q=0
  rsvd_out <- rsvd(testMat, k=k, p=0, q=0, method='standard')
  Ak = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
  test_that("Randomized SVD k=k, p=0, q=0", {
    expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    expect_equal(testMat, Ak)
  })

#Fast randomized SVD k=k, p=0, q=0
  rsvd_out <- rsvd(testMat, k=k, p=0, q=0, method='fast')
  Ak = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
  test_that("Fast Randomized SVD k=n, p=0, q=0", {
    expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    expect_equal(testMat, Ak)
  })



#*************************************************************************************
# Test 2: testMat.T
#*************************************************************************************

testMat <- t(testMat)
#Deterministic SVD
  svd_out <- svd(testMat)
  #Ak = svd_out$u %*% diag(svd_out$d) %*% t(svd_out$v)
  #print( sum((Ak-test)**2) )

#Randomized SVD k=n
  rsvd_out <- rsvd(testMat, method='standard')
  Ak = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
  test_that("Randomized SVD k=n", {
    expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    expect_equal(testMat, Ak)
  })

#Fast randomized SVD k=n
  rsvd_out <- rsvd(testMat, method='fast')
  Ak = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
  test_that("Fast Randomized SVD k=n", {
    expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    expect_equal(testMat, Ak)
  })

#Randomized SVD k=k
  rsvd_out <- rsvd(testMat, k=k, method='standard')
  Ak = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
  test_that("Randomized SVD k=k", {
    expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    expect_equal(testMat, Ak)
  })

#Fast randomized SVD k=k
  rsvd_out <- rsvd(testMat, k=k, method='fast')
  Ak = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
  test_that("Fast Randomized SVD k=k", {
    expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    expect_equal(testMat, Ak)
  })

#Randomized SVD k=k, p=0, q=0
  rsvd_out <- rsvd(testMat, k=k, p=0, q=0, method='standard')
  Ak = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
  test_that("Randomized SVD k=k, p=0, q=0", {
    expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    expect_equal(testMat, Ak)
  })

#Fast randomized SVD k=k, p=0, q=0
  rsvd_out <- rsvd(testMat, k=k, p=0, q=0, method='fast')
  Ak = rsvd_out$u %*% diag(rsvd_out$d) %*% t(rsvd_out$v)
  test_that("Fast Randomized SVD k=n, p=0, q=0", {
    expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    expect_equal(testMat, Ak)
  })




#*************************************************************************************
# Test 3: complex random test matrix
#*************************************************************************************
  testMat <- matrix(runif(m*k), m, k) + 1i* matrix(runif(m*k), m, k)
  testMat <- testMat %*% H(testMat)
  testMat <- testMat[,1:n]


  #Deterministic SVD
  svd_out <- svd(testMat)
  #Ak = svd_out$u %*% diag(svd_out$d) %*% t(svd_out$v)
  #print( sum((Ak-test)**2) )

  #Randomized SVD k=n
  rsvd_out <- rsvd(testMat, method='standard')
  Ak = rsvd_out$u %*% diag(rsvd_out$d) %*% H(rsvd_out$v)
  test_that("Randomized SVD k=n", {
    expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    expect_equal(testMat, Ak)
  })

  #Fast randomized SVD k=n
  rsvd_out <- rsvd(testMat, method='fast')
  Ak = rsvd_out$u %*% diag(rsvd_out$d) %*% H(rsvd_out$v)
  test_that("Fast Randomized SVD k=n", {
    expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    expect_equal(testMat, Ak)
  })

  #Randomized SVD k=k
  rsvd_out <- rsvd(testMat, k=k, method='standard')
  Ak = rsvd_out$u %*% diag(rsvd_out$d) %*% H(rsvd_out$v)
  test_that("Randomized SVD k=k", {
    expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    expect_equal(testMat, Ak)
  })

  #Fast randomized SVD k=k
  rsvd_out <- rsvd(testMat, k=k, method='fast')
  Ak = rsvd_out$u %*% diag(rsvd_out$d) %*% H(rsvd_out$v)
  test_that("Fast Randomized SVD k=k", {
    expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    expect_equal(testMat, Ak)
  })

  #Randomized SVD k=k, p=0, q=0
  rsvd_out <- rsvd(testMat, k=k, p=0, q=0, method='standard')
  Ak = rsvd_out$u %*% diag(rsvd_out$d) %*% H(rsvd_out$v)
  percent_error = 100*norm(testMat - Ak,'2')/norm(Ak,'2')
  test_that("Randomized SVD k=k, p=0, q=0", {
    expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    expect_less_than(percent_error, 10)
  })

  #Fast randomized SVD k=k, p=0, q=0
  rsvd_out <- rsvd(testMat, k=k, p=0, q=0, method='fast')
  Ak = rsvd_out$u %*% diag(rsvd_out$d) %*% H(rsvd_out$v)
  percent_error = 100*norm(testMat - Ak,'2')/norm(Ak,'2')

  test_that("Fast Randomized SVD k=n, p=0, q=0", {
    expect_equal(svd_out$d[1:k], rsvd_out$d[1:k])
    expect_less_than(percent_error, 10)
  })


