# rsvd 1.0.7

## Correctness and compatibility

* Deterministic `rqb()` now honors `k` and constructs only the requested columns
  of the orthogonal factor. Randomized `rqb()` continues to return its oversampled
  basis, with width bounded by both matrix dimensions.
* Complex `rpca()` uses squared moduli to compute real scale factors.
* `rrpca()` handles zero input and singular-value thresholding steps with no
  retained singular values. Missing rows are omitted as documented.
* Numeric data frames work consistently in `rsvd()`, `rid()` and `rcur()`.
  ID and CUR omit missing rows from sparse Matrix inputs without densification.
* The internal pseudoinverse preserves matrix dimensions when one singular
  direction remains and returns the correct shape for rectangular zero matrices.
* Public function arguments, defaults, result classes and component names are
  unchanged. The existing singular-vector count clamping and sampling
  distributions are retained.

## Efficiency

* Randomized `rqb()` keeps sparse Matrix inputs sparse during multiplication and
  converts only sketches and projected factors to dense matrices. Randomized
  ID and CUR inherit this improvement. Deterministic QR paths remain dense.
* The random sketch in `rqb()` is bounded by the smaller input dimension. Very
  wide calls may therefore consume fewer random draws than earlier versions.
* `rpca(retx = FALSE)` avoids computing unused left singular vectors.
* Full-column interpolative decomposition avoids a discarded pseudoinverse, and
  CUR reuses the selected row and column matrices.

## Tests and documentation

* Add the testthat entry point so `R CMD check` executes the test suite.
* Check Matrix 1.7-6 on current R across Linux, Windows and macOS, plus
  Matrix 1.6-5 on R 4.3. Verify the actual Matrix version in every CI job.
* Add prescribed-spectrum residual and subspace tests, sparse representation
  coverage, rank-one and vector-count cases, and numerical bug regressions.
* Correct complex orthogonality and PCA reference calculations in the tests,
  preserve the robust-PCA oversampling parameter, and remove checks for
  undeclared transitive dependencies.
* Document approximation error, sparse storage boundaries, QB oversampling,
  missing-value handling, and the existing `unif` sampler's `[0, 1]` support.
  Correct spelling and mathematical notation in the documentation.
