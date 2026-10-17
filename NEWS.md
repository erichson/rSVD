# rsvd 2.0.0

## Correctness and compatibility

* Require R 4.4 or later and Matrix 1.7-0 or later.
* Deterministic `rqb()` now honors `k` and constructs only the requested columns
  of the orthogonal factor. Randomized `rqb()` continues to return its oversampled
  basis, with width bounded by both matrix dimensions.
* Complex `rpca()` uses squared moduli to compute real scale factors.
* With `center = FALSE`, `rpca()` reports total processed squared energy,
  consistent with `prcomp()` and its uncentered SVD. Already-centered inputs
  retain the same interpretation. Prediction accepts numeric data frames when
  scaling is enabled.
* `rrpca()` handles zero input and singular-value thresholding steps with no
  retained singular values. Missing rows are omitted as documented.
  Dual initialization uses the entrywise maximum norm, and saturated
  singular-value thresholding expands the sampled rank with bounds based on
  the smaller matrix dimension.
* `rrpca()` checks both reconstruction feasibility and stationarity, and only
  increases its penalty when the stationarity test is satisfied. Stalled
  randomized updates receive bounded power-iteration refinements and then an
  exact-SVD fallback if needed. This can require more iterations and time than
  the former reconstruction-only stop. The returned `err` remains the
  reconstruction history; `trace = TRUE` also reports stationarity.
* Numeric data frames work consistently in `rsvd()`, `rid()` and `rcur()`.
  ID and CUR omit missing rows from sparse Matrix inputs without densification.
* PCA plotting handles a single retained component and omitted fitted scores,
  routes named options to the appropriate helpers, and preserves the invisible
  return of the last drawn plot.
  Use supported line-width and legend arguments to avoid deprecation warnings
  with current ggplot2, while retaining the plotted widths and public arguments.
* The internal pseudoinverse preserves matrix dimensions when one singular
  direction remains and returns the correct shape for rectangular zero matrices.
  Its numerical-rank cutoff now accounts for matrix dimensions, preventing
  CUR from amplifying roundoff singular values when `k` exceeds the true rank.
  Small but resolved singular directions remain supported.
* Public function arguments, defaults, result classes and component names are
  unchanged. The existing singular-vector count clamping and sampling
  distributions are retained.

## Efficiency

* Randomized `rqb()` keeps sparse Matrix inputs sparse during multiplication and
  converts only sketches and projected factors to dense matrices. Randomized
  ID and CUR inherit this improvement.
* Deterministic `rqb()` factors only the first `k` columns when a real,
  numerically independent prefix gives the same base-QR basis. Complex input,
  deficient prefixes and other cases retain the full dense QR fallback.
* Randomized `rpca()` with an explicit target rank centers and scales real
  sparse Matrix inputs implicitly. Prediction also avoids forming the dense
  centered input. Sparse column moments handle implicit zeros, nearly constant
  columns and overflow in mean accumulation.
* The random sketch in `rqb()` is bounded by the smaller input dimension. Very
  wide calls may therefore consume fewer random draws than earlier versions.
* `rpca(retx = FALSE)` avoids computing unused left singular vectors.
  Total-variance calculations reuse column norms or a Frobenius norm instead
  of applying sample-variance calculations to copied columns.
* Full-column interpolative decomposition avoids a discarded pseudoinverse, and
  CUR reuses the selected row and column matrices.

## Tests and documentation

* Add the testthat entry point so `R CMD check` executes the test suite.
* Check Matrix 1.7-6 on current R across Linux, Windows and macOS, plus
  Matrix 1.7-0 on R 4.4. Verify the actual Matrix version in every CI job.
* Add prescribed-spectrum residual and subspace tests, sparse representation
  coverage, rank-one and vector-count cases, and numerical bug regressions.
* Compare PCA with `prcomp()` and explicitly processed dense references,
  including sparse triangular/unit-diagonal and triplet representations.
  Allocation checks guard against restoring a full dense sparse-PCA input.
  Verify partial QR against the full base-QR factors.
* Test robust-PCA rank selection, scale behavior and recovery, including
  closed-form low-rank-only and sparse-only optima that distinguish optimality
  from a small reconstruction residual.
* Verify pseudoinverses against analytic maps and the Moore-Penrose identities,
  and test ID/CUR when the requested rank exceeds the input rank.
* Check PCA plot rendering, named options, return values, and the distinction
  between retained and total variance. Clarify the loading/correlation
  interpretation and the existing 95% multivariate-t ellipse.
  Render all four plotting helpers with lifecycle deprecations treated as errors.
* Cover sparse product storage, rank-one structural zeros, and custom matrix
  dispatch on the supported Matrix stack.
* Correct complex orthogonality and PCA reference calculations in the tests,
  preserve the robust-PCA oversampling parameter, and remove checks for
  undeclared transitive dependencies.
* Document approximation error, sparse storage boundaries, QB oversampling,
  missing-value handling, and the existing `unif` sampler's `[0, 1]` support.
  Correct spelling and mathematical notation in the documentation.
  Clarify that PCA prediction uses training column order, sample variances need
  at least two complete observations, and zero total variance gives undefined
  explained-variance proportions.
* Replace historical README performance images with reproducible comparisons
  against deterministic SVD, PCA, ID, and CUR. Include measured residuals,
  cumulative R allocations, raw timing observations, and vector figures.
* Expand the README benchmarks to matrices starting at 10,000 by 3,000,
  and restore an illustrated tiger example explaining rank, approximation
  error, and storage of the factors.
* Explain low-rank SVD, sketching costs, and interpretable ID/CUR factors in
  the README, with short display equations and executable reconstruction
  examples. Preserve the saved benchmark measurements and speedup figures.
* Clarify that ID and CUR are approximations, that ID conditioning depends on
  the selected subset, and that `idx_only = TRUE` still computes and returns
  interpolation or connector factors. Document row indices after missing-row
  omission and provide readable plain-text versions of the help equations.
* Review all help pages against the implementation, with complete argument and
  return-value descriptions, readable mathematics, cross-references, and small
  executable examples. Document the existing PCA prediction and summary
  methods, and distinguish batch robust PCA from prediction on new data.
* Correct the `digits` dataset description: its first column contains labels
  and the remaining 784 columns contain pixels. Exclude the label when
  reshaping a digit image, and document the stored scales of both datasets.
* Explain the PCA plotting helpers' coordinate scaling, label selection,
  retained-versus-total variance, and optional ellipse requirements. Remove
  obsolete contact addresses from help pages; the maintainer contact remains
  in `DESCRIPTION`.
* Expand the package and function descriptions with explanations of the
  factors, preprocessing, approximation error, and convergence diagnostics.
  Revise the README and figure captions around the mathematical ideas and
  measured results, preserving the benchmark data and executable examples.
  Qualify the variable-arrow correlation interpretation for randomized PCA.
