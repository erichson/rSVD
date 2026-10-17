<img src="https://raw.githubusercontent.com/erichson/rSVD/master/rsvd.png" width="550">

[![R-CMD-check](https://github.com/erichson/rSVD/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/erichson/rSVD/actions/workflows/R-CMD-check.yaml)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/rsvd)](https://CRAN.R-project.org/package=rsvd)

## Fast Randomized Singular Value Decomposition using R

Randomized singular value decomposition (rsvd) is a fast probabilistic algorithm that can
be used to compute the near optimal low-rank singular value decomposition of massive data sets with high accuracy.
The key idea is to compute a compressed representation
of the data to capture the essential information. This compressed representation can then be used to obtain
the low-rank singular value decomposition. The package also provides randomized QB, interpolative (ID), and CUR decompositions.
The computational advantage is most useful when the target rank is small compared with both matrix dimensions. The following historical comparison uses target rank k=50; timings depend on matrix structure, hardware, and BLAS implementation:

![speed](https://raw.githubusercontent.com/Benli11/data/master/img/rsvd_speedups.png)

The singular value decomposition plays a central role in data analysis and scientific computing.
The SVD is also widely used for computing
(randomized) principal component analysis (PCA), a linear dimensionality reduction technique.
Randomized PCA (rpca) uses the approximated singular value decomposition
to compute the most significant principal components. This package also includes a
function to compute (randomized) robust principal component analysis (RPCA).
In addition several plot functions are provided. See for further details: [Randomized Matrix Decompositions using R](https://arxiv.org/abs/1608.02148).


## SVD example: Image compression

```R
library(rsvd)
data(tiger)

# Image compression using randomized SVD
s <- rsvd(tiger, k=150)
tiger.re = s$u %*% diag(s$d) %*% t(s$v) # reconstruct image

# Display original and reconstructed image
par(mfrow=c(1,2))
image(tiger, col = gray((0:255)/255))
image(tiger.re, col = gray((0:255)/255))
```
Here are the results:
![tiger](https://raw.githubusercontent.com/Benli11/data/master/img/reTiger.png)

and the speedup gained over the base SVD function:

```R
library(microbenchmark)

timing_svd <- microbenchmark(
  'SVD' = svd(tiger, nu=150, nv=150),
  'rSVD' = rsvd(tiger, k=150),
  times=50)

print(timing_svd, unit='s')
```
![timing](https://raw.githubusercontent.com/Benli11/data/master/img/timeing.png)


## Choosing rank and assessing accuracy

`rsvd(A, k, p = 10, q = 2)` computes an approximate rank-`k` decomposition.
Choose `k` much smaller than `min(dim(A))`. Oversampling `p` adds sketch
columns; power iterations `q` improve separation when the spectrum decays
slowly, at the cost of additional passes over `A`.

Accuracy depends on singular-value decay. A full-rank matrix with a decaying
spectrum can still have a useful low-rank approximation. For a flat spectrum,
even the optimal rank-`k` approximation leaves a substantial residual. If the
exact singular values are `d`, the best possible Frobenius residual is
`sqrt(sum(d[-seq_len(k)]^2))`; numerical tests should compare with that scale.
Exact recovery is appropriate for an exact low-rank fixture, not for every input.

For a small real matrix, a reconstruction can be formed without an explicit
diagonal matrix:

```r
set.seed(42)
fit <- rsvd(A, k = 10)
A_approx <- tcrossprod(sweep(fit$u, 2, fit$d, "*"), fit$v)
```

For complex data, use `Conj(t(fit$v))` in the reconstruction. Singular vectors
can change sign or complex phase, and repeated singular values allow rotations
within a subspace. Compare residuals or subspaces when validating results.
`set.seed()` makes sampling reproducible within a fixed numerical environment;
vector entries need not be identical across BLAS, Matrix, or package versions.

Omitting `k` requests all singular values. That remains supported for API
compatibility, but it loses the efficiency of a small target rank. Use `svd()`
when a full deterministic decomposition is the goal.

## Sparse matrices and memory

Sparse matrices from **Matrix** are accepted. `library(Matrix)` is not required;
the package imports the matrix operations it uses.

| Routine | Storage behavior |
| --- | --- |
| `rsvd()` | Retains the input representation; samples, the basis, projected problem and singular vectors are generally dense. |
| `rqb(rand = TRUE)` | Retains sparse Matrix inputs and forms dense sketches and factors. The returned basis includes oversampling. |
| `rid(rand = TRUE)`, `rcur(rand = TRUE)` | Use sparse products for the randomized sketch, then dense projected or selected factors. |
| `rqb(rand = FALSE)` | Uses only a real, independent leading column block when it gives the same base-QR basis; otherwise uses full dense QR. |
| Deterministic ID and CUR | Use dense pivoted QR; a sparse input can require a full dense copy. |
| `rpca(rand = TRUE)` with an explicit `k` | Centers and scales real sparse Matrix inputs implicitly; loadings, sketches and optional scores are dense. |
| `rpca(rand = FALSE)` or `rpca(k = NULL)` | Uses dense preprocessing and deterministic SVD. |
| `rrpca()` | Uses dense input and working arrays for its low-rank, sparse and residual updates. |

For sketch width `ell = min(k + p, min(dim(A)))`, the randomized methods need
thin dense arrays whose storage grows with `(nrow(A) + ncol(A)) * ell`, in
addition to the input and factorization workspace. Reconstructing a large
`A_approx` or computing `A - A_approx` allocates a dense matrix, even when the
input is sparse. Keep the factors when a full reconstruction is unnecessary.

`rpca(retx = FALSE)` skips unused left singular vectors and omits scores.
`rid(idx_only = TRUE)` and `rcur(idx_only = TRUE)` omit the returned selected
subsets, but still compute their interpolation or connector factors.

Matrix 1.7.x requires R 4.4 or later. rsvd retains its R 4.0 minimum and supports
older compatible Matrix releases; its Matrix imports preserve sparse dispatch
on R versions before 4.4. CI checks current R with Matrix 1.7-6 on Linux,
Windows and macOS, and R 4.3 with Matrix 1.6-5. Each job verifies the loaded
Matrix version before running the package check.

## PCA with sparse data

An explicit small target rank enables implicit centering and scaling:

```r
A <- Matrix::bdiag(outer(1:100, 1:20), outer(1:100, 20:1))
set.seed(42)
fit <- rpca(A, k = 2, center = TRUE, scale = TRUE, retx = FALSE)
scores <- predict(fit, A[1:5, , drop = FALSE])
```

The centered matrix is applied through products with thin factors. Sparse
working storage grows with the number of stored entries, and dense workspace
grows with the sketch width. Column scales use deviations from the mean
directly, including the contribution of implicit zeros. Nearly constant,
fully stored columns are centered before multiplication to reduce cancellation.

As with `prcomp()`, `center = FALSE` simply skips centering. It works for data
that are already centered and also permits an uncentered SVD. With centering
disabled, `scale = TRUE` divides each column by
`sqrt(sum(x^2) / (nrow(A) - 1))`; this is generally not its centered standard
deviation. `rpca()` retains its existing rule that scales below `1e-8` become 1.

`fit$var` includes the energy of all directions in the processed input, so
variance proportions do not depend on how many components were requested.
For uncentered input this is a total second moment. Fitted scores keep the
existing `U D` convention; `predict()` projects observations onto the fitted
loadings. These can differ by the randomized range-finder residual.

## Robust PCA convergence

`rrpca()` seeks a low-rank matrix and an entrywise sparse matrix whose sum
matches the input, minimizing the nuclear norm of the first plus `lambda`
times the entrywise absolute sum of the second. A small reconstruction
residual alone does not establish that this split is optimal.

The iteration now checks both reconstruction feasibility and stationarity.
It increases the penalty only when the stationarity test is satisfied,
expands a sampled SVD when all returned singular values survive thresholding,
and refines stalled randomized updates before falling back to exact SVD if
necessary. These checks can require more work than the earlier stopping rule.

The returned `err` keeps its existing meaning as the reconstruction-error
history. Use `trace = TRUE` to see the additional stationarity diagnostic.
`maxiter` remains a work limit and can stop an unconverged iteration. The
algorithm still uses dense working arrays, even when its sparse component
contains many zeros.

## Tests and development

The test suite runs as part of `R CMD check`. It checks low-rank reconstruction,
orthogonality, residuals relative to the optimal truncated-SVD error, complex
arithmetic, sparse representations, and API edge cases.

```sh
R CMD build .
R CMD check --no-manual rsvd_1.0.7.tar.gz
```

For an interactive development check, use `testthat::test_local()` with the
package dependencies installed. See [NEWS.md](NEWS.md) for the 1.0.7 changes
and [the numerical review](tools/validation/follow-up-review.md) for reference
tests, measured costs, compatibility decisions and recent randomized NLA
experiments. `tools/benchmark-lowrank.R` provides reproducible timing and
R-allocation checks.

## Installation

Install the rsvd package via CRAN
```R
install.packages("rsvd")
```

You can also install the development version from GitHub using [devtools](https://cran.r-project.org/package=devtools):

```r
devtools::install_github("erichson/rsvd")
```

The source package can be obtained here: [CRAN: rsvd](https://CRAN.R-project.org/package=rsvd).

## New in Version 1.0.7

See [NEWS.md](NEWS.md) for numerical fixes, sparse-memory improvements and test coverage.

## New in Version 1.0.2

* Several small issues are fixed.
* Thanks to Aaron Lun, who has fixed a bug in the rsvd function that occurred when nu=0 or nv=0.

## New in Version 1.0.0

* Support for non-default matrix types to deal with large-scale matrices that are held on file, added by Aaron Lun.
* Fixed a bug which occurred running rpca with k=1 and retx=TRUE, discovered by Will.


## References

* [Erichson NB, Voronin S, Brunton SL, Kutz JN (2019). Randomized Matrix Decompositions Using R. Journal of Statistical Software, 89(11), 1–48. doi: 10.18637/jss.v089.i11.](https://www.jstatsoft.org/article/view/v089i11)
* [Sergey Voronin, Per-Gunnar Martinsson. RSVDPACK: Subroutines for computing partial singular value decompositions via randomized sampling on single core, multi core, and GPU architectures. (2015)](https://arxiv.org/abs/1502.05366)
* [Nathan Halko, et al. Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions. (2011)](https://arxiv.org/abs/0909.4061)

## Cite as

```
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


