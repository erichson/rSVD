# rsvd 2.0.0

## Requirements and compatibility

* Require R 4.4 or later and Matrix 1.7-0 or later. CI checks the supported
  minimum and current R on Linux, Windows, and macOS, and verifies the loaded
  Matrix version in every job.
* Preserve public function signatures, defaults, result classes, and component
  names. Singular-vector counts retain their clamping rules, and the existing
  sampling distributions are unchanged.
* Validate ranks, flags, and sampling parameters before numerical work.
  Reject invalid randomized settings, including negative oversampling that
  could produce missing singular values. Reject empty, nonnumeric, and
  nonfinite matrix inputs with direct errors. Routines that support missing
  data still omit incomplete rows with a warning. Numeric data frames work
  consistently across the decompositions.
* Require at least two complete observations for PCA. Prediction checks that
  new data have the training number of variables and still uses their supplied
  column order. For compatibility, ordinary missing values propagate and
  `NaN` values encountered during scaling become zero.

## Numerical fixes

* Deterministic `rqb()` honors `k` and returns only the requested basis columns.
  Randomized QB retains its oversampling columns, with width bounded by both
  matrix dimensions.
* Complex PCA uses squared magnitudes to compute real scale factors. With
  `center = FALSE`, total processed energy agrees with the uncentered SVD;
  variance proportions use all processed directions, including omitted ones.
* The internal pseudoinverse preserves dimensions when one singular direction
  remains and returns the correct shape for rectangular zero matrices. Its
  rank cutoff accounts for matrix dimensions, preventing CUR from inverting
  roundoff singular values when the requested rank exceeds the true rank.
  Small, resolved singular directions remain supported.
* Robust PCA handles zero input and thresholding steps that retain no singular
  values. Dual initialization uses the entrywise maximum norm, and saturated
  singular-value thresholding increases the sampled rank within the matrix
  dimensions.
* Robust PCA checks both reconstruction feasibility and stationarity. The
  penalty increases only after the stationarity test is satisfied. Stalled
  randomized updates receive bounded power-iteration refinements and then an
  exact-SVD fallback when needed. These checks can require more time than the
  former reconstruction-only stop. `err` remains the reconstruction history;
  `trace = TRUE` also displays stationarity. A loose tolerance still completes
  an update and returns its measured residual.
* PCA plots handle one retained component and fits without stored scores.
  `plot.rpca()` routes named options to the helpers that accept them and
  invisibly returns the last plot. Factor maps require valid integer component
  indices, matching group and label lengths, and real coordinates. `top.n`
  must be a positive integer. Plotting uses supported ggplot2 line-width and
  legend arguments.

## Memory and runtime

* Randomized QB keeps sparse Matrix inputs sparse during multiplication,
  converting only sketches and projected factors to dense matrices.
  Randomized ID and CUR use the same sparse products. Missing-row omission
  in ID and CUR also preserves sparse storage.
* Deterministic QB factors only the leading `k` columns when a real,
  independent prefix gives the same basis as full base QR. Complex inputs,
  deficient prefixes, and other cases use the full dense QR fallback.
* Randomized PCA with an explicit rank centers and scales real sparse inputs
  implicitly. Prediction applies the same transformation without creating a
  full dense centered matrix. Sparse column moments account for implicit zeros,
  nearly constant columns, and overflow in mean accumulation.
* `rpca(retx = FALSE)` avoids unused left singular vectors. Total-variance
  calculations reuse column norms or a Frobenius norm instead of applying
  sample-variance calculations to copied columns.
* Full-column ID avoids a discarded pseudoinverse. CUR reuses its selected
  row and column matrices. Bounding QB sketch width by the smaller dimension
  can also reduce the number of random draws for very wide inputs.

## Documentation and tests

* Add the testthat entry point so package checks run the test suite. Cover
  prescribed spectra, invariant subspaces, complex arithmetic, vector counts,
  rank-deficient decompositions, input validation, and sparse storage formats.
  Pseudoinverse tests use analytic examples and the Moore-Penrose identities.
* Compare PCA with `prcomp()` and explicitly processed dense references. Check
  implicit sparse products, structural zeros, custom matrix dispatch, and
  allocation limits. Compare partial QR with the full base-QR factors.
* Test robust-PCA rank selection, rescaling, and recovery, including known
  low-rank-only and sparse-only optima. These distinguish an appropriate split
  from an arbitrary split with a small reconstruction residual.
* Render the PCA plotting helpers and check component selection, variance
  denominators, return values, and argument routing. Treat lifecycle
  deprecations as errors in the relevant tests.
* Edit the README, package description, and all help pages for direct language
  and concrete explanations. Add an early reconstruction example and method
  guide, explain how rank differs from oversampling and power iterations, and
  document sparse storage, preprocessing, and robust-PCA stopping diagnostics.
* Clarify PCA fitted scores versus new-data projections, total versus retained
  variance, and the loading/correlation interpretation. Explain that robust
  PCA fits a batch of observations and does not return a prediction method.
* Correct the digit dataset description and image example: the first column
  contains labels, followed by 784 pixel columns. Complete the tiger image
  attribution, and update the citation contact to match the maintainer.
* Preserve the recorded README benchmarks and their raw data, including
  reconstruction errors and cumulative R allocations. Identify their source
  revision and distinguish them from checks of later release candidates.
  Provide reproducible tiger and iris illustrations alongside the performance
  examples.
