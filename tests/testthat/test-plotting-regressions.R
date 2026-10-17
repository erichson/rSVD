testthat::skip_if_not_installed("ggplot2")

plot_test_fixture <- function(k = 3L, retx = TRUE) {
  # Three centered orthogonal directions have energies 18, 8, and 2. Their
  # known total lets summary/plot tests distinguish total from retained energy.
  A <- cbind(c(3, -3, 0, 0, 0, 0),
             c(0, 0, 2, -2, 0, 0),
             c(0, 0, 0, 0, 1, -1))
  dimnames(A) <- list(paste0("sample", 1:6), c("first", "second", "third"))
  rpca(A, k = k, rand = FALSE, center = TRUE, scale = FALSE, retx = retx)
}

plot_test_with_device <- function(expr) {
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  # Fail on deprecated plotting calls even if they have already warned in this
  # session. Ordinary warnings remain available for explicit expect_warning().
  old_options <- options(lifecycle_verbosity = "error")
  on.exit(options(old_options), add = TRUE)
  force(expr)
}

testthat::test_that("PCA plotting helpers render without deprecated aesthetics", {
  fit <- plot_test_fixture()
  plots <- list(
    function() ggbiplot(fit, ellipse = FALSE, ind_labels = FALSE),
    function() ggcorplot(fit),
    function() ggindplot(fit, ellipse = FALSE, ind_labels = FALSE),
    function() ggscreeplot(fit)
  )
  for (make_plot in plots) {
    grob <- plot_test_with_device(ggplot2::ggplotGrob(make_plot()))
    testthat::expect_s3_class(grob, "gtable")
  }
})

testthat::test_that("PCA summary and scree plots use total processed variance", {
  fit <- plot_test_fixture(k = 2L)
  summary <- summary(fit)
  testthat::expect_identical(dim(summary), c(4L, 2L))
  testthat::expect_equal(unname(summary["Proportion of variance", ]),
                         round(c(18, 8) / 28, 3))
  testthat::expect_equal(unname(summary["Cumulative proportion", ]),
                         round(cumsum(c(18, 8) / 28), 3))
  ratio <- plot_test_with_device(ggscreeplot(fit, type = "ratio"))
  cumulative <- plot_test_with_device(ggscreeplot(fit, type = "cum"))
  testthat::expect_equal(ratio$data$y, c(18, 8) / 28, tolerance = 1e-12)
  testthat::expect_equal(cumulative$data$y, cumsum(c(18, 8) / 28),
                         tolerance = 1e-12)
  testthat::expect_lt(tail(cumulative$data$y, 1), 1)
})

testthat::test_that("plot.rpca handles one component and omitted scores", {
  for (retx in c(TRUE, FALSE)) {
    fit <- plot_test_fixture(k = 1L, retx = retx)
    actual <- plot_test_with_device(withVisible(plot(fit)))
    testthat::expect_false(actual$visible)
    testthat::expect_true(inherits(actual$value, "ggplot"))
    testthat::expect_equal(actual$value$data$y, 18 / 5, tolerance = 1e-12)
    testthat::expect_identical(actual$value$labels$x, "Principal components")
  }
  fit <- plot_test_fixture(k = 2L, retx = FALSE)
  actual <- plot_test_with_device(withVisible(plot(fit)))
  testthat::expect_false(actual$visible)
  testthat::expect_true(inherits(actual$value, "ggplot"))
  testthat::expect_identical(actual$value$labels$title, "Variables factor map (PCA)")
})

testthat::test_that("plot.rpca retains its default last-plot return and score names", {
  fit <- plot_test_fixture()
  actual <- plot_test_with_device(withVisible(plot(fit)))
  testthat::expect_false(actual$visible)
  testthat::expect_true(inherits(actual$value, "ggplot"))
  testthat::expect_identical(actual$value$labels$title,
                             "Individuals factor map (PCA)")
  testthat::expect_equal(actual$value$data$a, unname(fit$x[, 1L]))
  testthat::expect_equal(actual$value$data$b, unname(fit$x[, 2L]))
  testthat::expect_identical(actual$value$data$indName, rownames(fit$x))
})

testthat::test_that("plot.rpca routes each supported named argument to its helper", {
  fit <- plot_test_fixture()
  groups <- factor(rep(c("A", "B"), each = 3))
  labels <- paste0("observation", 1:6)
  actual <- plot_test_with_device(withVisible(plot(
    fit, type = "cum", pcs = c(2L, 3L), loadings = FALSE,
    groups = groups, alpha = 0.4, ellipse = FALSE, alpha.ellipse = 0.1,
    var_labels = TRUE, var_labels.names = c("x", "y", "z"), top.n = 2L,
    ind_labels = TRUE, ind_labels.names = labels)))
  testthat::expect_false(actual$visible)
  testthat::expect_equal(actual$value$data$a, unname(fit$x[, 2L]))
  testthat::expect_equal(actual$value$data$b, unname(fit$x[, 3L]))
  testthat::expect_identical(actual$value$data$class, groups)
  testthat::expect_identical(actual$value$data$indName, labels)
  testthat::expect_match(actual$value$labels$x, "^PC2\\(")
  testthat::expect_match(actual$value$labels$y, "^PC3\\(")
  testthat::expect_equal(actual$value$layers[[1L]]$aes_params$alpha, 0.4)

  # With scores omitted, the returned loading plot exposes its routed values.
  no_scores <- plot_test_fixture(retx = FALSE)
  loading_plot <- plot_test_with_device(plot(
    no_scores, pcs = c(2L, 3L), loadings = FALSE,
    var_labels = TRUE, var_labels.names = c("x", "y", "z"), top.n = 2L))
  segments <- Filter(function(layer) inherits(layer$geom, "GeomSegment"),
                     loading_plot$layers)
  testthat::expect_length(segments, 1L)
  testthat::expect_equal(segments[[1L]]$data$a,
                         unname(no_scores$rotation[, 2L]))
  testthat::expect_equal(segments[[1L]]$data$b,
                         unname(no_scores$rotation[, 3L]))
  testthat::expect_identical(segments[[1L]]$data$varName, c("x", "y", "z"))
})

testthat::test_that("plot.rpca rejects unknown or ambiguous wrapper arguments", {
  fit <- plot_test_fixture()
  plot_test_with_device({
    testthat::expect_error(plot(fit, unknown_option = TRUE),
                           "Unknown plot argument.*unknown_option")
    testthat::expect_error(plot(fit, "cum"), "must be named")
    testthat::expect_error(plot(fit, alpha = 0.2, alpha = 0.3),
                           "supplied once")
    testthat::expect_error(plot(fit, pcs = c(1L, 4L)),
                           "Selected PC is not valid")
    testthat::expect_error(plot(plot_test_fixture(k = 1L), pcs = c(1L, 2L)),
                           "Selected PC is not valid")
    testthat::expect_error(plot(fit, pcs = c(1L, Inf)),
                           "Selected PC is not valid")
    testthat::expect_error(plot(fit, pcs = 1L), "Selected PC is not valid")
  })
})
