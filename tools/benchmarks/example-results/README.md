# PCA teaching example

This small iris example illustrates PCA output and plotting; it is not a
performance benchmark. Species is supplied only to the plot for coloring.

```r
library(rsvd)
data(iris)
set.seed(42)
pca <- rpca(as.matrix(iris[, 1:4]), k = 2, center = TRUE, scale = TRUE)
ggbiplot(pca, groups = iris$Species, ellipse = FALSE, ind_labels = FALSE)
```

Each point is a flower's fitted score in the two retained components. The
arrows use `rotation * sqrt(eigvals)`: with centered, standardized variables,
these are variable–component correlations. The observations retain their
original score scale. Axis labels report proportions of **total variance**,
including the two omitted components. A change of sign in a PC reverses both
its scores and arrows together and leaves the interpretation unchanged.

The figure extends the public `ggbiplot()` result with a
color legend, equal coordinate units, a title, and separated variable labels.
The short dotted leaders connect those labels to the original arrow endpoints;
neither scores nor arrow endpoints are rescaled. Full labels for the two almost
parallel petal arrows would otherwise overlap.

Run the complete figure script from the repository root:

```sh
Rscript tools/benchmarks/illustrate-pca.R \
  --results tools/benchmarks/example-results \
  --figures man/figures --pdfs tools/benchmarks/figures
```

Use `--library /path/to/R/library` if needed. The script requires the installed
candidate package, ggplot2, and R Cairo graphics support. It writes PNG and SVG
README assets, a PDF export, scores, loadings, variances, and environment details.
It also checks the variance fractions against `prcomp()` and arrow coordinates
against the observed variable–score correlations before rendering.

## Larger sparse-input example

The top-level README also includes an executable 10,000 by 3,000 sparse
example. Its [recorded check](sparse-example-check.txt) confirms 1,500,000
nonzero entries (5% density), 3,000 by 20 loadings, omitted training scores,
and finite 5 by 20 predictions with rsvd 2.0.0. Rerun that README code to
exercise the same fitted preprocessing and prediction path. It is an API
example; performance measurements use the separately specified benchmark
fixtures and rank settings.

The recorded source commit is the checkout HEAD at generation time. Source-file
hashes identify the numerical code independently of version or documentation
edits in the working tree.
