testthat::skip_if_not_installed("ggplot2")

release_plot_fixture <- function() {
  A <- cbind(c(3, -3, 0, 0, 0, 0), c(0, 0, 2, -2, 0, 0),
             c(0, 0, 0, 0, 1, -1))
  rpca(A, k = 3, rand = FALSE, scale = FALSE)
}

testthat::test_that("factor-map labels cannot hide truncated PC indices", {
  fit <- release_plot_fixture()
  for(pcs in list(c(1.5, 2), c(1, NA), c(1, Inf), c("1", "2"),
                  c(1+0i, 2+0i), c(0, 2), c(1, 4))) {
    for(helper in list(ggbiplot, ggcorplot, ggindplot)) {
      testthat::expect_error(helper(fit, pcs = pcs), "Selected PC is not valid")
    }
    testthat::expect_error(plot(fit, pcs = pcs), "Selected PC is not valid")
  }
})

testthat::test_that("factor-map groups and labels cannot be recycled", {
  fit <- release_plot_fixture()
  for(helper in list(ggbiplot, ggindplot)) {
    testthat::expect_error(helper(fit, groups = c("a", "b")),
                           "'groups'.*one entry")
    testthat::expect_error(helper(fit, ind_labels.names = "one label"),
                           "'ind_labels.names'.*one entry")
  }
  for(helper in list(ggbiplot, ggcorplot)) {
    testthat::expect_error(helper(fit, var_labels.names = "one variable"),
                           "'var_labels.names'.*one entry")
  }
  groups <- factor(c("a", "a", "b", "b", "c", "c"))
  labels <- paste0("sample", seq_len(nrow(fit$x)))
  p <- ggindplot(fit, groups = groups, ind_labels.names = labels,
                 ellipse = FALSE)
  testthat::expect_identical(p$data$class, groups)
  testthat::expect_identical(p$data$indName, labels)
})

testthat::test_that("variable selection requires a positive integer", {
  fit <- release_plot_fixture()
  for(top.n in list(0, -1, 1.5, NA, Inf, c(1, 2))) {
    testthat::expect_error(ggcorplot(fit, top.n = top.n), "'top.n'")
  }
  p <- ggcorplot(fit, top.n = 20, var_labels = TRUE)
  points <- Filter(function(layer) inherits(layer$geom, "GeomPoint"), p$layers)
  testthat::expect_equal(nrow(points[[1]]$data), nrow(fit$rotation))
})

testthat::test_that("factor maps reject complex coordinates explicitly", {
  fit <- release_plot_fixture()
  fit$x <- fit$x + 1i * fit$x
  fit$rotation <- fit$rotation + 1i * fit$rotation
  for(helper in list(ggbiplot, ggcorplot, ggindplot)) {
    testthat::expect_error(helper(fit), "real-valued")
  }
  # Variance remains real-valued, so the scree plot is still supported.
  testthat::expect_true(inherits(ggscreeplot(fit), "ggplot"))
})
