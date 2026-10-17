#!/usr/bin/env Rscript
# A small teaching example using the public rpca() and ggbiplot() functions.
# This example illustrates interpretation, with no claim about performance.
#
# From the repository root, with the candidate package and ggplot2 installed:
#   Rscript tools/benchmarks/illustrate-pca.R \
#     --results tools/benchmarks/example-results \
#     --figures man/figures --pdfs tools/benchmarks/figures
# Optional: --library /path/to/private/R/library --repo /path/to/rSVD

args <- commandArgs(trailingOnly = TRUE)
options <- list(results = "tools/benchmarks/example-results",
                figures = "man/figures", pdfs = "tools/benchmarks/figures",
                library = NULL, repo = ".")
if (length(args) %% 2L != 0L) stop("Supply options as --name value pairs.")
if (length(args)) {
  for (i in seq.int(1L, length(args), by = 2L)) {
    key <- sub("^--", "", args[i])
    if (!startsWith(args[i], "--") || !key %in% names(options)) {
      stop("Unknown option: ", args[i])
    }
    options[key] <- args[i + 1L]
  }
}
if (!is.null(options$library)) .libPaths(c(options$library, .libPaths()))
for (package in c("rsvd", "Matrix", "ggplot2")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Install ", package, " first.")
}
for (directory in options[c("results", "figures", "pdfs")]) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}
repo <- normalizePath(options$repo, mustWork = TRUE)

# Species is retained only for the plot, not supplied to the PCA fit.
A <- as.matrix(datasets::iris[, 1:4])
set.seed(42)
fit <- rsvd::rpca(A, k = 2, center = TRUE, scale = TRUE)
stopifnot(identical(dim(fit$x), c(150L, 2L)),
          identical(dim(fit$rotation), c(4L, 2L)),
          all(is.finite(fit$x)), all(is.finite(fit$rotation)))
variance_fraction <- fit$sdev^2 / fit$var
arrows <- sweep(fit$rotation, 2L, sqrt(fit$eigvals), "*")

# Check the quantities described in the caption against independent formulas.
reference <- stats::prcomp(A, center = TRUE, scale. = TRUE, rank. = 2)
variance_difference <- max(abs(variance_fraction -
  reference$sdev[1:2]^2 / sum(reference$sdev^2)))
correlation_difference <- max(abs(stats::cor(A, fit$x) - arrows))
stopifnot(variance_difference < 1e-12, correlation_difference < 1e-12)

# The existing public plot supplies the scores and variable arrows. Labels
# are positioned separately to distinguish the almost-parallel petal arrows.
# Suppress only ggplot2 lifecycle notices from legacy helper syntax; other
# warnings, including any rendering/data warning, remain visible.
g <- withCallingHandlers(
  rsvd::ggbiplot(fit, groups = datasets::iris$Species, ellipse = FALSE,
                 ind_labels = FALSE, var_labels = FALSE, alpha = 0.78),
  lifecycle_warning_deprecated = function(w) invokeRestart("muffleWarning")
)
labels <- data.frame(variable = gsub("\\.", " ", rownames(arrows)),
                     a = arrows[, 1], b = arrows[, 2])
labels$label_x <- ifelse(labels$a > 0, 1.28, -1.28)
labels$label_y <- labels$b
labels$hjust <- ifelse(labels$a > 0, 0, 1)
for (side in c(-1, 1)) {
  index <- which(sign(labels$a) == side)
  index <- index[order(labels$b[index])]
  if (length(index) > 1L) {
    labels$label_y[index] <- mean(labels$b[index]) +
      (seq_along(index) - (length(index) + 1) / 2) * 0.48
  }
}

g <- g +
  ggplot2::geom_segment(data = labels,
    ggplot2::aes(x = a, y = b, xend = label_x, yend = label_y),
    color = "#BC7770", linewidth = 0.35, linetype = "dotted",
    inherit.aes = FALSE) +
  ggplot2::geom_label(data = labels,
    ggplot2::aes(x = label_x, y = label_y, label = variable, hjust = hjust),
    color = "#A83C33", fill = "white", linewidth = 0, size = 3.55,
    label.padding = grid::unit(0.14, "lines"), inherit.aes = FALSE) +
  ggplot2::scale_colour_manual(values = c(setosa = "#087F8C",
    versicolor = "#BD5B29", virginica = "#5268A0")) +
  ggplot2::guides(colour = ggplot2::guide_legend(title = NULL,
    override.aes = list(alpha = 1, size = 3))) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "PCA: four measurements, two coordinates",
    subtitle = "150 iris flowers; each measurement is centered and scaled to unit variance",
    x = sprintf("PC1  (%.1f%% of total variance)", 100 * variance_fraction[1]),
    y = sprintf("PC2  (%.1f%% of total variance)", 100 * variance_fraction[2]),
    caption = sprintf(paste0(
      "Points: fitted scores. Arrows: variable–PC correlations. The two PCs retain %.1f%% of total variance.\n",
      "Species labels are used only for color. Dotted leaders connect labels to the original arrow endpoints."),
      100 * sum(variance_fraction))) +
  ggplot2::theme_minimal(base_size = 11, base_family = "DejaVu Sans") +
  ggplot2::theme(
    text = ggplot2::element_text(color = "#263C53"),
    plot.title = ggplot2::element_text(size = 16, face = "bold", margin = ggplot2::margin(b = 7)),
    plot.subtitle = ggplot2::element_text(size = 10.5, color = "#465866", margin = ggplot2::margin(b = 6)),
    plot.caption = ggplot2::element_text(size = 8.3, color = "#465866", hjust = 0,
      lineheight = 1.25, margin = ggplot2::margin(t = 12)),
    plot.title.position = "plot", plot.caption.position = "plot",
    legend.position = "top", legend.justification = "left", legend.location = "plot",
    panel.grid.major = ggplot2::element_line(color = "#E6EBEE", linewidth = 0.3),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(color = "#465866"),
    axis.title = ggplot2::element_text(size = 10.5),
    plot.margin = ggplot2::margin(15, 18, 12, 15))

# A build exercises all layers before exporting to the three file formats.
built <- ggplot2::ggplot_build(g)
stopifnot(nrow(built$data[[1]]) == 150L,
          all(is.finite(built$data[[1]]$x)), all(is.finite(built$data[[1]]$y)))
width <- 8.0
height <- 6.3
grDevices::png(file.path(options$figures, "iris-pca-biplot.png"),
               width = width, height = height, units = "in", res = 180, type = "cairo")
print(g)
grDevices::dev.off()
grDevices::svg(file.path(options$figures, "iris-pca-biplot.svg"),
               width = width, height = height, family = "DejaVu Sans")
print(g)
grDevices::dev.off()
grDevices::cairo_pdf(file.path(options$pdfs, "iris-pca-biplot.pdf"),
                     width = width, height = height, family = "DejaVu Sans")
print(g)
grDevices::dev.off()

write.csv(data.frame(component = paste0("PC", 1:2), eigenvalue = fit$sdev^2,
                     total_variance = fit$var, variance_fraction = variance_fraction),
          file.path(options$results, "iris-variance.csv"), row.names = FALSE)
write.csv(data.frame(observation = seq_len(nrow(A)), species = datasets::iris$Species,
                     PC1 = fit$x[, 1], PC2 = fit$x[, 2]),
          file.path(options$results, "iris-scores.csv"), row.names = FALSE)
write.csv(data.frame(variable = rownames(fit$rotation),
                     rotation_PC1 = fit$rotation[, 1], rotation_PC2 = fit$rotation[, 2],
                     correlation_PC1 = arrows[, 1], correlation_PC2 = arrows[, 2]),
          file.path(options$results, "iris-loadings.csv"), row.names = FALSE)
metadata <- c(
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  package = as.character(packageVersion("rsvd")), R = R.version.string,
  Matrix = as.character(packageVersion("Matrix")),
  ggplot2 = as.character(packageVersion("ggplot2")),
  BLAS = unname(extSoftVersion()["BLAS"]),
  RNG = paste(RNGkind(), collapse = "; "), seed = "42",
  call = "rpca(as.matrix(iris[, 1:4]), k = 2, center = TRUE, scale = TRUE)",
  plotting_call = "ggbiplot(fit, groups = iris$Species, ellipse = FALSE, ind_labels = FALSE, var_labels = FALSE, alpha = 0.78)",
  source_commit = paste(system2("git", c("-C", shQuote(repo), "rev-parse", "HEAD"),
                                stdout = TRUE), collapse = ""),
  md5_rpca = unname(tools::md5sum(file.path(repo, "R", "rpca.R"))),
  md5_ggbiplot = unname(tools::md5sum(file.path(repo, "R", "ggbiplot.R"))),
  retained_variance_fraction = format(sum(variance_fraction), digits = 17),
  max_variance_fraction_difference_vs_prcomp = format(variance_difference, digits = 17),
  max_arrow_difference_vs_observed_correlations = format(correlation_difference, digits = 17),
  aspect_ratio = "equal PC units on both axes",
  label_placement = "separated near-parallel labels, dotted leaders; arrow endpoints unchanged",
  PNG_pixels = paste(round(c(width, height) * 180), collapse = " x "),
  purpose = "small educational example; no performance measurement"
)
write.csv(data.frame(key = names(metadata), value = unname(metadata)),
          file.path(options$results, "iris-metadata.csv"), row.names = FALSE)
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
           file.path(options$results, "iris-session-info.txt"))
message(sprintf("PC1 %.9f%%, PC2 %.9f%%; total %.9f%%",
                100 * variance_fraction[1], 100 * variance_fraction[2],
                100 * sum(variance_fraction)))
