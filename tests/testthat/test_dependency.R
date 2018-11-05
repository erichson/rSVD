#devtools::use_package("testthat")
#library('testthat')

context("Dependency")

#Load rsvd library
library(rsvd)

#*************************************************************************************
# Dependency
#*************************************************************************************
testthat::test_that("Dependency", {
  testthat::expect_equal(requireNamespace("ggplot2", quietly = TRUE), TRUE)
  testthat::expect_equal(requireNamespace("plyr", quietly = TRUE), TRUE)
  testthat::expect_equal(requireNamespace("scales", quietly = TRUE), TRUE)
  testthat::expect_equal(requireNamespace("grid", quietly = TRUE), TRUE)
})
