<img src="rsvd.png" width="550" alt="rsvd: randomized matrix decompositions">

[![R-CMD-check](https://github.com/erichson/rSVD/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/erichson/rSVD/actions/workflows/R-CMD-check.yaml)
[![CRAN version](https://www.r-pkg.org/badges/version/rsvd)](https://CRAN.R-project.org/package=rsvd)

# Randomized matrix decompositions in R

Large matrices often contain much less independent information than their
dimensions suggest. A photograph has repeated patterns; many measured
variables move together. A low-rank approximation represents those patterns
with a small set of directions and weights.

**rsvd** finds such approximations for dimension reduction, compression, and
data analysis. Its central idea is to use random combinations of the data to
find the important directions, then solve a much smaller matrix problem.
The benefit is greatest when the requested rank is small relative to the
matrix dimensions and the singular values decay.

The package provides SVD, PCA, QB, interpolative (ID), CUR, and robust-PCA
decompositions. It accepts ordinary R matrices and supports sparse matrices
from **Matrix**, with implicit centering and scaling for randomized PCA.
The [package paper](https://www.jstatsoft.org/article/view/v089i11) explains
the algorithms and their applications.

## Installation

```r
install.packages("rsvd")

# Development version:
# install.packages("remotes")
# remotes::install_github("erichson/rSVD")
```

The development version is **2.0.0**; the CRAN badge shows the published version.
rsvd requires **R 4.4 or later and Matrix 1.7-0 or later**. Matrix need not be
attached with `library(Matrix)` for sparse dispatch to work.
The minimum R version follows the Matrix dependency; newer R releases are
supported and tested in CI alongside the R 4.4 series.

## What does a low-rank factorization represent?

Think of a matrix as a collection of columns: measurements, images, or
time series. Its **rank** counts the independent directions those columns
span. If every column can be assembled from ten underlying patterns, the
matrix has rank at most ten, however many columns it contains. A factorization
stores those patterns and the coefficients needed to combine them.

Measured data rarely have exactly low rank. Noise can make a matrix full rank
while most of its structure still lies in a much smaller subspace. A **low-rank
approximation** retains that structure and accepts a residual: the difference
between the original data and their reconstruction. This can make storage,
visualization, and subsequent calculations more manageable. Whether the
discarded detail is acceptable depends on the application; low rank alone
does not distinguish noise from a weak but important signal.

The **singular value decomposition (SVD)** supplies a particularly useful set
of patterns. Keeping its leading `k` terms gives the approximation

```math
A \approx U_k D_k V_k^{T}.
```

The columns of **U** and **V** are orthonormal left and right directions;
the diagonal of **D** contains their nonnegative weights, called singular
values, ordered from largest to smallest. Each paired set of directions
contributes a rank-one pattern. Larger singular values contribute more
squared energy, so retaining the leading terms captures the strongest
structure first. The exact truncated SVD is the best rank-`k` approximation
in Frobenius norm; randomized SVD seeks comparable accuracy with less work.

The equations here use **real matrices**, with superscript **T** denoting
transpose. A complex SVD uses **H**, the conjugate transpose: transpose the
matrix and conjugate its entries. The package supports this distinction;
reconstruction code must respect it too.

## See what rank means: reconstructing a tiger

A grayscale image is a matrix of pixel intensities. Its leading SVD terms
capture broad shapes and repeated patterns; additional terms recover finer
detail. The bundled tiger image lets us see the approximation directly.

```r
library(rsvd)
data(tiger)

set.seed(42)
fit <- rsvd(tiger, k = 150)

# Keep fit$u, fit$d, and fit$v for a compact representation.
# Reconstruct only when a dense image is needed; avoid an explicit diag(fit$d).
tiger_approx <- tcrossprod(sweep(fit$u, 2, fit$d, "*"), fit$v)
norm(tiger - tiger_approx, "F") / norm(tiger, "F")
```

For complex data, reconstruct with
`sweep(fit$u, 2, fit$d, "*") %*% Conj(t(fit$v))`.

![Original tiger and randomized reconstructions at ranks 10, 50, and 150, with approximation errors and factor sizes.](man/figures/tiger-rank-comparison.png)

Increasing `k` retains more detail and requires more storage. For this
1,600-by-1,200 image, rank 50 uses **7.3% of the original number of numeric
entries**, with a relative Frobenius error of **15.1%**. Rank 150 uses 21.9%
of the entries and reduces the error to 10.3%.

An `m`-by-`n` matrix stores `m * n` numbers; its rank-`k` SVD factors store
`k * (m + n + 1)`. This comparison counts numeric entries in the matrix and
factors; it is not a comparison with JPEG or PNG file compression. Reconstructing
the image creates a full dense matrix again, so retain the factors when a
compact representation is the goal.

The factors are also useful without reconstruction. To multiply the
approximation by a vector, first project the vector onto the columns of V,
scale the resulting `k` coefficients by the singular values, and combine
the columns of U. The work then grows with `k * (m + n)` instead of `m * n`.
This saving applies to each subsequent product after fitting, in addition
to any benefit from computing the factorization itself.

Image: S. Taheri (2006), edited by Fir0002,
[CC BY-SA 2.5](https://creativecommons.org/licenses/by-sa/2.5/),
via [Wikimedia Commons](https://commons.wikimedia.org/wiki/File:Sibirischer_tiger_de_edit02.jpg).
The package supplies a grayscale version; the panels add low-rank
reconstructions. [Reproduce this figure](tools/benchmarks/tiger-results/README.md).

## The core idea: find the important directions first

Suppose a tall matrix A has `m` rows and `n` columns, but we need only
`k` dominant directions. Randomized SVD first seeks a small subspace that
captures the dominant structure. Three steps turn that subspace into
an approximation.

**1. Mix the columns into a small sketch.** Draw a random matrix, denoted
by the Greek letter omega, with `n` rows and `ell = min(k + p, m, n)` columns.
The sketch includes up to `p` extra columns for oversampling; the default
uses Gaussian entries. Form the sketch:

```math
Y = A\Omega.
```

Then compute a thin QR factorization of Y and keep its orthonormal factor Q.
Each column of Y is a random weighted combination of the columns of A.
Multiplication by A gives more weight to directions with large singular
values. When a few directions dominate, their span can be captured by a
small sketch. Orthonormalization gives that span a well-conditioned basis;
oversampling helps capture directions near the truncation boundary.

**2. Express the original matrix in that basis.** Project onto the columns
of Q:

```math
B = Q^{T}A.
```

The product QB now approximates A. The matrix B has only `ell` rows, which
is where the large SVD becomes a small one. For a **20,000 × 6,000** matrix
with `k = 40` and `p = 10`, B is **50 × 6,000**: 300,000 entries instead of
120 million. The basis Q has 20,000 × 50 entries. These are the core
compressed arrays; products with the input and additional factorization
workspace still contribute to cost.

**3. Solve the small problem and map back.** Compute the ordinary SVD of
B and retain its leading `k` terms:

```math
B \approx \widetilde{U}_k D_k V_k^{T}.
```

Map the left singular vectors back to the original space:

```math
U_k = Q\widetilde{U}_k.
```

The diagonal entries of D are singular values of B and approximate
those of A; `fit$d` stores them as a vector. This produces the SVD factors
used in the tiger example. Randomness enters the subspace search; the
decomposition inside that subspace uses conventional numerical linear algebra.
The same idea applies to wide matrices through the conjugate-transposed problem.

When the spectrum decays slowly, power iterations refine the basis before
projection, sharpening the separation between dominant and trailing directions.
With `q` iterations, the idealized
sketch is

```math
Y = (AA^{T})^q A\Omega.
```

Each singular value is effectively raised to the power `2 * q + 1`,
emphasizing the larger ones. The implementation applies successive products
with A and its transpose, orthonormalizing after each product; it never
needs to build the large matrix power or the matrix of cross-products.
More iterations can improve accuracy but require more passes over the data.
The [package paper](https://www.jstatsoft.org/article/view/v089i11) develops
these ideas with worked examples; the
[Halko–Martinsson–Tropp survey](https://arxiv.org/abs/0909.4061) gives the
error analysis.

### Where the speed and memory savings come from

For a dense tall matrix, a full SVD requires work on the scale of `m * n^2`.
With a fixed number of power iterations, the main randomized products instead
cost on the scale of `m * n * ell`, followed by thin QR and SVD work on the
scale of `(m + n) * ell^2`. The advantage is substantial when the sketch
width is much smaller than the input width. Dense products also make effective
use of optimized BLAS libraries. Small inputs or ranks close to the full
dimension leave less room to recover the sketching overhead.

Randomization still requires access to the input. The default `q = 2` uses
six products with A or its transpose, including the final projection. The
saving comes from applying them to a small number of vectors together and
factorizing the compressed problem. More power iterations can improve the
subspace at the cost of additional products.

The dominant thin arrays grow with `(m + n) * ell`; a dense input still
occupies `m * n` entries. For sparse input, sketch products can use only the
stored entries, although the sketches themselves are generally dense.
Sparsity and low rank describe different properties: many zeros do not
guarantee a good low-rank approximation. The sparse benchmarks below benefit
from both compressible structure and avoiding full dense preprocessing.
Their reported allocation totals include temporary R objects, so those
numbers measure more than the sizes of the retained factors.

## PCA: from measurements to a map

Low-rank SVD is also the computational workhorse behind PCA. After
centering, and optionally scaling, the right singular vectors give principal
loadings; the left singular vectors weighted by singular values give fitted
scores. A few score columns provide a compact map of the observations.

This small example makes that map easy to interpret. The larger benchmarks
below demonstrate performance. Install `ggplot2` to use the plotting helpers.

```r
library(rsvd)
data(iris)
set.seed(42)
pca <- rpca(as.matrix(iris[, 1:4]), k = 2, center = TRUE, scale = TRUE)
ggbiplot(pca, groups = iris$Species, ellipse = FALSE,
         ind_labels = FALSE)
```

![PCA biplot of iris flowers with observation scores, variable directions, and the variance retained by the first two components.](man/figures/iris-pca-biplot.svg)

Each point is a flower's fitted score in two coordinates. For these centered,
standardized measurements, the arrows show variable–component correlations.
The nearly parallel petal arrows indicate that both measurements vary in a
similar direction. The two components retain **95.8% of total variance**;
species labels are used only for color, not for fitting the PCA.

The figure adds a legend and separates the variable labels for readability.
[Reproduce the styled figure and its numerical checks](tools/benchmarks/example-results/README.md).

## ID and CUR: keep actual columns and rows

An SVD direction can combine every variable. That is useful for compression,
but sometimes the original variables or observations are easier to interpret.
An **interpolative decomposition (ID)** represents the matrix using a subset
of its own columns. If the columns are sensors or measured features, the
selected indices identify which original measurements form the approximation.
For a column ID,

```math
A \approx CZ.
```

C contains `k` selected columns of A, and Z contains the coefficients that
approximately reconstruct all columns from that subset. At the selected
column indices, Z contains an identity matrix: those columns reproduce
themselves. `rid()` returns the subset in `C`, the coefficients in `Z`, and
the original column indices in `idx`. With `mode = "row"`, it instead returns
selected rows in `R` and reconstructs as `Z %*% R`.

**CUR** preserves both sides of the data: selected columns in C and selected
rows in R, connected by a small matrix U:

```math
A \approx CUR.
```

For an observations-by-variables matrix, C retains actual variables and R
retains actual observations. The connector U has size `k`-by-`k`; it is not
the orthonormal left factor from an SVD. `rcur()` returns `C.idx` and `R.idx`
alongside these three factors, so selections can be traced back to the data.
If rows contain missing values, both routines omit those rows with a warning;
returned row indices then refer to the remaining complete observations.

The tiger makes these definitions concrete: C and R contain selected pixel
strips from the matrix, and the factors reconstruct the entire image from them.

```r
set.seed(42)
id <- rid(tiger, k = 50)
id$idx[1:6]                     # Original column positions
tiger_id <- id$C %*% id$Z

set.seed(42)
cur <- rcur(tiger, k = 50)
tiger_cur <- cur$C %*% cur$U %*% cur$R

c(ID = norm(tiger - tiger_id, "F") / norm(tiger, "F"),
  CUR = norm(tiger - tiger_cur, "F") / norm(tiger, "F"))
```

These subset constraints usually trade some approximation accuracy for
interpretability. Selected columns need not be orthogonal, and different
valid selections can give different errors. Selection also depends on the
units and scaling of the input: ID and CUR do not automatically standardize
variables. Use their residuals to judge whether the chosen subset is adequate;
selected indices describe an approximation, not a causal feature ranking.
The [package paper](https://www.jstatsoft.org/article/view/v089i11) discusses
the pivoted-QR constructions used by these routines.

## Available decompositions

| Function | Result and typical use |
| --- | --- |
| `rsvd(A, k)` | Leading singular vectors and values; low-rank compression. |
| `rpca(A, k)` | Principal loadings, component variances, and optional scores. |
| `rqb(A, k)` | An orthonormal basis Q and a compressed matrix B. The randomized basis includes oversampling. |
| `rid(A, k)` | Selected columns and interpolation coefficients; row ID is also available. |
| `rcur(A, k)` | Selected columns, selected rows, and their connector matrix. |
| `rrpca(A)` | A low-rank component and an entrywise sparse corruption component. |

ID and CUR retain interpretable subsets of the data; their residuals need
not match a truncated SVD.

## Choosing rank and assessing accuracy

`rsvd(A, k, p = 10, q = 2)` has three main accuracy/cost controls:

- **`k`** is the requested rank. A small rank gives the largest computational
  advantage. Choose it according to the information the application needs.
- **`p`** is a nonnegative integer oversampling parameter. Extra sketch columns
  help capture the target subspace; the sketch is bounded by the matrix dimensions.
- **`q`** is a nonnegative integer number of power iterations. Increasing it
  can improve accuracy for a slowly decaying spectrum, at the cost of more
  passes over the input. The ID and CUR defaults use `q = 0`.

A full-rank matrix can still have a useful low-rank approximation. What matters
is the energy left outside the leading directions. Let A<sub>k</sub> denote
the exact truncated SVD, and let the singular values of A be ordered from
largest to smallest. Its Frobenius residual is

```math
\|A-A_k\|_F^2 = \sum_{j>k}\sigma_j^2.
```

The Frobenius norm is the square root of the sum of squared magnitudes of
matrix entries. This identity says that the best possible squared reconstruction error is
exactly the energy in the discarded singular values. A rapidly decaying
spectrum leaves little energy in that tail. When `k` is much smaller than
the matrix dimensions, a flat spectrum leaves a large
relative residual even for the optimal approximation.

Tests and benchmarks therefore compare reconstruction error with the optimal
truncation error, and check orthogonality and subspaces. Exact recovery is a
suitable expectation for an exact low-rank fixture. Singular-vector entries
can change sign, complex phase, or rotate within a repeated singular subspace.
`set.seed()` controls sampling in a fixed numerical environment; it does not
guarantee identical vectors across BLAS or package versions.

Omitting `k` in `rsvd()` requests the full available rank. This remains supported,
including downstream uses such as Procrustes alignment, but offers no small-rank
advantage. Base `svd()` is a natural reference for a full deterministic SVD.

## Measured performance

The following results compare **the same requested rank and returned factors**
with deterministic references. They were measured on Linux x86_64 with R 4.6.1,
Matrix 1.7-6, and OpenBLAS 0.3.26 using four BLAS threads. Times are medians of
five calls after warmup, with method order rotated between repetitions.
The synthetic matrices contain a low-rank signal plus noise; they are not exact
low-rank inputs. Scaling cases range from **10,000 × 3,000 to 20,000 × 6,000**;
dense PCA, ID, and CUR use **15,000 × 4,500** matrices.

**Speedup** is the deterministic median time divided by the randomized
median time, calculated from the unrounded measurements. A 10× speedup means
one tenth of the reference runtime. In the tables, paired values are always
**deterministic / randomized**. The figures retain absolute costs and add
speedup and allocation-reduction panels, so both the benefit and its scale
are visible.

**R allocations** are the cumulative bytes allocated during a separate profiled
call, expressed in GiB (`2^30` bytes). They are not peak resident memory or live
workspace, and do not include every allocation made inside native libraries. Fixture creation,
explicit garbage collection, and error calculation are outside the measurement;
preprocessing and any required sparse-to-dense conversion are inside it.

Relative error is the Frobenius norm of the reconstruction residual divided
by the Frobenius norm of the input, using the centered/scaled input in PCA.
Values below are percentages, not variance
explained. The [benchmark guide](tools/benchmarks/README.md) includes the exact
fixtures, baseline calls, raw measurements, environment, and rerun commands.

### SVD as matrix dimensions grow

The reference is `base::svd(A, nu = k, nv = k)`. Base SVD computes the full
singular spectrum even when only `k` vectors are requested; it is not a
specialized truncated-SVD solver. rsvd uses `p = 10, q = 2` here.

<!-- BEGIN BENCHMARK SVD -->

| Input | k | Time (s) | Speedup | R allocations (GiB) | Relative error |
| --- | ---: | ---: | ---: | ---: | ---: |
| Dense, 10,000 × 3,000 | 40 | 16.858 / 0.697 | **24.2×** | 1.234 / 0.133 | 9.67% / 9.68% |
| Dense, 20,000 × 6,000 | 40 | 97.157 / 2.144 | **45.3×** | 4.928 / 0.267 | 9.67% / 9.68% |

<!-- END BENCHMARK SVD -->

![Dense SVD elapsed time, speedup over base SVD, cumulative R allocations, and allocation reduction at three sizes from 10,000 by 3,000 to 20,000 by 6,000.](man/figures/benchmark-svd-scaling.svg)

At 20,000 × 6,000, the measured advantage is a **45.3× speedup** and an
**18.5× reduction in cumulative R allocations**. The relative reconstruction
errors are 9.67% for the deterministic reference and 9.68% for rsvd; most of
that residual is the unavoidable rank-40 truncation error.

The plot holds `k = 40` while the input grows from **30 million to 120 million
entries**. Its largest dense input alone occupies about 0.89 GiB in double
precision, before any factorization workspace.
Shading shows the interquartile range of the recorded times; whiskers show
their minimum and maximum. Ratio panels use a 1× reference line, with larger
values indicating a greater speedup or allocation reduction. Matrix size and
the ratio of sketch width to input width matter; the measured advantage is
not a universal speedup guarantee.

### PCA on dense and sparse input

Both methods center and scale columns and return `k = 40` loadings and score
columns. The reference is
`prcomp(as.matrix(A), rank. = k, center = TRUE, scale. = TRUE)`. The sparse
fixtures have 5% nonzero entries in twenty noisy blocks. Randomized PCA uses
`p = 10, q = 2`; for sparse inputs, it applies centering/scaling through
sparse products.

<!-- BEGIN BENCHMARK PCA -->

| Input | k | Time (s) | Speedup | R allocations (GiB) | Relative error |
| --- | ---: | ---: | ---: | ---: | ---: |
| Dense PCA, 15,000 × 4,500 | 40 | 54.085 / 5.256 | **10.3×** | 8.056 / 2.981 | 10.3% / 10.3% |
| Sparse-input PCA, 10,000 × 3,000 | 40 | 18.336 / 0.463 | **39.6×** | 3.806 / 0.213 | 0.516% / 0.516% |
| Sparse-input PCA, 20,000 × 6,000 | 40 | 110.634 / 1.976 | **56.0×** | 15.212 / 0.482 | 0.548% / 0.548% |

<!-- END BENCHMARK PCA -->

![Sparse-input PCA elapsed time, speedup over prcomp, cumulative R allocations, and allocation reduction, comparing dense preprocessing with implicit sparse centering.](man/figures/benchmark-pca-scaling.svg)

The third curve, `rpca(as.matrix(A), ...)`, isolates the effect of preserving
sparse storage from the effect of using a randomized SVD. All three curves
use the same observations, preprocessing, and requested output. The maximum
observed error relative to the optimal rank-40 error is printed on each figure.

At the largest size, implicit sparse PCA gives a **56.0× speedup** and a
**31.5× reduction in cumulative R allocations** relative to `prcomp()`.
It is also **5.1× faster than randomized PCA after densifying the same input**.
Preserving sparse storage therefore contributes a substantial part of the
benefit. The speedup for implicit sparse PCA peaks at the middle size;
runtime variability and machine behavior mean it need not rise monotonically
with dimensions.

Here “sparse PCA” means **ordinary PCA computed from sparse input**. The
principal loadings are generally dense; `rpca()` does not impose sparsity on
the loadings.

### Interpretable factors: ID and CUR

The deterministic references are `rid(rand = FALSE)` and `rcur(rand = FALSE)`,
which use LAPACK column-pivoted QR. Base R has no corresponding ID or CUR
function. Both paths return the selected data subsets and their interpolation
or connector factors; randomized calls use `k = 40, p = 10, q = 0`.

<!-- BEGIN BENCHMARK ID_CUR -->

| Input | k | Time (s) | Speedup | R allocations (GiB) | Relative error |
| --- | ---: | ---: | ---: | ---: | ---: |
| Column ID, dense, 15,000 × 4,500 | 40 | 43.727 / 0.684 | **63.9×** | 1.257 / 0.326 | 14.3% / 15.8% |
| CUR, dense, 15,000 × 4,500 | 40 | 42.141 / 1.076 | **39.2×** | 1.546 / 0.615 | 18.1% / 19.2% |
| Column ID, sparse, 15,000 × 4,500 | 40 | 46.938 / 0.287 | **163.5×** | 1.522 / 0.095 | 0.153% / 0.153% |
| CUR, sparse, 15,000 × 4,500 | 40 | 40.352 / 0.346 | **116.6×** | 1.572 / 0.145 | 0.189% / 0.189% |

<!-- END BENCHMARK ID_CUR -->

These methods select columns and rows, so different valid selections can have
different residuals. Column-pivoted QR is not a guarantee of optimal rank-`k`
approximation on every matrix. The table reports the measured errors alongside
cost, including cases where a randomized selection is less accurate.

## PCA behavior and sparse storage

Centering a sparse matrix usually turns its structural zeros into nonzero
values. The useful observation is that a centered matrix can still be
multiplied without storing that dense result. If the vector mu contains the
column means and **1** is a column of ones, then for a sketch matrix W,

```math
(A-\mathbf{1}\mu^{T})W = AW-\mathbf{1}(\mu^{T}W).
```

The first term uses the sparse input; the second is a rank-one
correction. `rpca()` uses this identity, together with implicit scaling, to
work with the centered/scaled data through products. Its sketches and
loadings are dense, but it avoids materializing the full centered input.

```r
set.seed(42)
# A 10,000 by 3,000 sparse matrix: 20 signal-plus-noise blocks, 5% density.
A <- Matrix::bdiag(replicate(20, {
  signal <- outer(seq_len(500) / 500, seq_len(150) / 150)
  signal + matrix(rnorm(500 * 150, sd = 0.001), 500, 150)
}, simplify = FALSE))
fit <- rpca(A, k = 20, center = TRUE, scale = TRUE, retx = FALSE)
scores <- predict(fit, A[1:5, , drop = FALSE])
```

An explicit target rank and `rand = TRUE` enable implicit centering/scaling for
real sparse Matrix inputs. Prediction uses the fitted transformation; supply
new columns in their **training order**. Names are not used to reorder them.
`retx = FALSE` omits training scores and avoids unused left singular vectors.

As with `prcomp()`, `center = FALSE` skips centering: this is appropriate for
already-centered data and also permits an uncentered SVD. For real columns with centering off,
`scale = TRUE` uses `sqrt(sum(x^2) / (nrow(A) - 1))`, generally not the centered
standard deviation. Unlike `prcomp()`'s default, `rpca()` defaults to scaling.
The existing rsvd convention replaces scales below `1e-8` with 1.

`fit$var` is the squared Frobenius norm of the processed observations divided
by `m - 1`, where `m` is the number of complete observations. With centering,
this is the sum of all component variances, including omitted components;
explained proportions use the same denominator regardless of `k`. Without
centering it is the total squared energy with the same normalization.
Fitted scores retain the `U D` convention, while
`predict()` projects observations onto the fitted loadings; these can differ
by the randomized approximation residual. Sample variances need at least two
complete observations, and explained proportions are undefined for zero energy.

| Routine | Sparse input and dense workspace |
| --- | --- |
| `rsvd()`, randomized `rqb()` | Keep sparse input for products; sketches, bases, and projected factors are generally dense. |
| Randomized `rid()`, `rcur()` | Use sparse sketch products and dense projected or selected factors. |
| `rpca(rand = TRUE, k = ...)` | Centers/scales real sparse input implicitly; loadings, sketches, and optional scores are dense. |
| `rpca(rand = FALSE)` or `rpca(k = NULL)` | Uses dense preprocessing and deterministic SVD. |
| Deterministic `rid()`, `rcur()` | Use dense pivoted QR, which can require a full dense copy. |
| Deterministic `rqb()` | Factors an independent real leading block when it reproduces the base-QR basis; otherwise falls back to full dense QR. |
| `rrpca()` | Uses dense input and dense working arrays. |

For sketch width `ell = min(k + p, min(dim(A)))`, thin randomized workspace
scales with `(nrow(A) + ncol(A)) * ell`, in addition to input storage and
factorization workspace. Sparse preprocessing also needs storage proportional
to the stored entries. A full reconstruction or residual matrix is dense even
when the original input is sparse. Keep factors when possible.

`rid(idx_only = TRUE)` and `rcur(idx_only = TRUE)` omit the selected data
subsets from the result, but still compute the interpolation or connector
factors. They are not index-selection-only algorithms.

## Robust PCA

`rrpca()` minimizes a nuclear norm plus `lambda` times an entrywise absolute
sum to separate a low-rank component from sparse corruption. A small
reconstruction residual alone does not establish that the split is optimal.
The solver checks reconstruction feasibility and stationarity, expands
saturated singular-value thresholding steps, and refines stalled randomized
updates before using an exact-SVD fallback when necessary.

`maxiter` remains a work limit and can stop before convergence. The returned
`err` is the reconstruction-error history; `trace = TRUE` also reports the
stationarity diagnostic. Stronger convergence checks can require more work,
and the solver still uses dense working arrays. See the
[release review](tools/validation/release-review.md) for optimality checks and
the limits of the validation.

## Tests, release notes, and reproducibility

```sh
R CMD build .
R CMD check rsvd_2.0.0.tar.gz
```

With development dependencies installed, `testthat::test_local()` runs the
suite interactively. Tests cover prescribed-spectrum approximation, invariant
subspaces, complex arithmetic, sparse representations, implicit PCA products,
rank-deficient ID/CUR, robust-PCA rank selection and convergence, and public
result shapes. CI checks the actual Matrix version before testing on Linux,
Windows, and macOS, including R 4.4 with Matrix 1.7-0 at the supported minimum.

- [NEWS.md](NEWS.md): changes in 2.0.0, including compatibility and memory behavior.
- [Release review](tools/validation/release-review.md): numerical audit,
  downstream checks, supported environments, and remaining limitations.
- [Benchmark guide and measurements](tools/benchmarks/README.md): reproduce the
  tables and figures, or substitute matrices representative of your workload.

## References and citation

- Erichson, N. B., Voronin, S., Brunton, S. L., and Kutz, J. N. (2019).
  [Randomized Matrix Decompositions Using R](https://www.jstatsoft.org/article/view/v089i11).
  *Journal of Statistical Software*, 89(11), 1–48.
- Halko, N., Martinsson, P.-G., and Tropp, J. A. (2011).
  [Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions](https://arxiv.org/abs/0909.4061).
- Voronin, S. and Martinsson, P.-G. (2015).
  [RSVDPACK: Subroutines for Computing Partial Singular Value Decompositions via Randomized Sampling](https://arxiv.org/abs/1502.05366).

```bibtex
@Article{,
  title = {Randomized Matrix Decompositions Using {R}},
  author = {N. Benjamin Erichson and Sergey Voronin and Steven L.
    Brunton and J. Nathan Kutz},
  journal = {Journal of Statistical Software},
  year = {2019},
  volume = {89},
  number = {11},
  pages = {1--48},
  doi = {10.18637/jss.v089.i11},
}
```


