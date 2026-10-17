# Validation notes for rsvd 1.0.7

These notes record the local numerical validation and performance comparison of
the 1.0.7 candidate against repository commit
[`bb5bb3b`](https://github.com/erichson/rSVD/commit/bb5bb3b), whose package version
is 1.0.6. Validation was performed on 2026-10-16, America/Los_Angeles.

## Compatibility and correctness

The installed-package test suite completed with the same results in both
environments:

| R | Matrix | Test blocks | Passing assertions | Failures | Warnings | Skips |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 4.6.1 | 1.7-6 | 128 | 1,267 | 0 | 0 | 0 |
| 4.3.3 | 1.6-5 | 128 | 1,267 | 0 | 0 | 0 |

The tests cover exact low-rank reconstruction, decaying-spectrum approximation
error relative to the optimal rank-k tail, Hermitian orthogonality for complex
data, and projector comparisons for repeated singular values. They also check
tall and wide inputs, rank-one and zero cases, oversampling, power iterations,
singular-vector suppression, data frames, and several sparse representations.
Timing assertions are kept out of the unit tests.

The API comparison found unchanged formals for all 25 existing functions and
unchanged names for all 10 exports. This includes argument order and defaults.
Numerical bug fixes can change previously incorrect results; the validation
does not require identical singular-vector signs or bases across BLAS and
Matrix versions.

### Package check status

The ordinary staged local `R CMD check --as-cran --no-manual` run completed
with **0 errors, 0 warnings, and 2 NOTEs**. Tests, examples, and documentation
checks passed. The notes concern a maintainer email change that predates this
patch and a `00LOCK-rsvd` directory left by installation in this environment.

A diagnostic rerun using fresh output and `--install-args=--no-lock` also
completed with **0 errors, 0 warnings, and 2 NOTEs**. Tests, examples, and
documentation passed, and the check for non-standard check-directory items
passed. Its notes concern the inherited maintainer email change and the
expected non-staged installation caused by `--no-lock`. This diagnostic option
is not a proposed change to the normal installation workflow. No further local
reruns are needed to distinguish this environment issue. Hosted CI results for
ordinary installation are tracked in the draft pull request.

### Matrix version and dispatch

The [CRAN Matrix page](https://cran.r-project.org/web/packages/Matrix/index.html)
listed **Matrix 1.7-6**, requiring **R >= 4.4**, at validation time. The number
4.4 is the R requirement. [Matrix NEWS](https://cran.r-project.org/web/packages/Matrix/news.html)
records that this requirement began with Matrix 1.7-0.

Base R's `crossprod()` and `tcrossprod()` became S3/S4 generics in R 4.4.0; see
the [R cross-product documentation](https://stat.ethz.ch/R-manual/R-devel/library/base/html/crossprod.html).
The package retains its imports of Matrix's transpose and cross-product
generics, which also provide dispatch under older R. The R 4.3.3 / Matrix 1.6-5
run exercises that compatibility path. The package's declared R >= 4.0.0
minimum was retained, but R 4.0 itself was not run locally.

## Benchmark method

The reproducible driver is [`../benchmark-lowrank.R`](../benchmark-lowrank.R).
The sanitized measurements are in
[`benchmark-comparison.csv`](benchmark-comparison.csv). The CSV retains numeric
ratios even where a small difference in elapsed time is not evidence of an
improvement. Machine-specific absolute paths have been removed.

Baseline and candidate ran in separate, sequential R processes using R 4.6.1,
Matrix 1.7-6, Linux x86_64, and the system reference BLAS
(`libblas.so.3.12.0`). `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, and
`MKL_NUM_THREADS` were each set to 1. Each case used one warmup, three elapsed-time
repetitions, and one separate `Rprofmem` run. Garbage collection and RNG setup
were outside elapsed-time measurements. Fixture construction and numerical
checks were outside both timing and allocation profiles. Every randomized call
used seed 731.

The sparse fixtures contain 20 disjoint rank-one blocks, with 400,000 nonzeros
in 8,000,000 entries, or 5% density. Their singular values are 20 through 1.
The dense fixtures use orthonormal cosine bases and a prescribed spectrum.
All benchmark inputs have exact low rank except the deliberately full-column-rank
deterministic ID case.

| Case | Dimensions | Target rank | Other settings |
| --- | ---: | ---: | --- |
| Sparse QB, tall | 4,000 x 2,000 | 20 | `p=10`, `q=2`, randomized |
| Sparse QB, wide | 2,000 x 4,000 | 20 | `p=10`, `q=2`, randomized |
| Dense SVD control | 2,200 x 1,000 | 20 | `p=10`, `q=2` |
| PCA without scores | 6,000 x 600 | 40 | `p=10`, `q=2`, randomized, `center=FALSE`, `scale=FALSE`, `retx=FALSE` |
| Full-column ID | 600 x 250 | 250 | `rand=FALSE`, column mode |

### Elapsed time and cumulative allocations

Times are the medians of three repetitions. Allocation sizes use MiB
(1,048,576 bytes).

| Case | Baseline time | Candidate time | Baseline cumulative R allocation | Candidate cumulative R allocation |
| --- | ---: | ---: | ---: | ---: |
| Sparse QB, tall | 1.403 s | 0.062 s | 95.382 MiB | 38.468 MiB |
| Sparse QB, wide | 1.286 s | 0.064 s | 92.635 MiB | 35.722 MiB |
| Dense SVD control | 0.328 s | 0.326 s | 20.019 MiB | 20.019 MiB |
| PCA without scores | 1.197 s | 1.195 s | 179.978 MiB | 178.132 MiB |
| Full-column ID | 0.085 s | 0.025 s | 14.739 MiB | 5.929 MiB |

For these sparse fixtures, removing the full-input dense conversion gives
observed time ratios of 22.6 and 20.1. The largest recorded R allocation falls
from 64,000,048 to 960,048 bytes, a 98.5% reduction. The baseline profile locates
the large allocation in `as.matrix.Matrix`. Cumulative recorded allocations
fall 59.7% for the tall case and 61.4% for the wide case.

The dense SVD control shows the same allocation volume and essentially the same
elapsed time. No dense SVD speed improvement is claimed. PCA without scores
saves about 1.846 MiB of cumulative allocation, or 1.0%, while its largest
allocation remains 27.466 MiB. Its elapsed times are also effectively unchanged;
no PCA speed improvement is claimed. Removing the discarded full-column ID
pseudoinverse gives an observed time ratio of 3.4 and a 59.8% reduction in
cumulative recorded allocation.

These ratios describe the measured fixtures and runtime. Sparsity pattern,
density, target rank, spectrum, and BLAS implementation can change the results.
Three repetitions do not establish a general timing guarantee.

### Numerical results and memory interpretation

Baseline and candidate had the same recorded relative Frobenius residuals:

| Case | Baseline and candidate residual |
| --- | ---: |
| Sparse QB, tall | 3.87097e-15 |
| Sparse QB, wide | 2.00637e-15 |
| Dense SVD control | 3.43995e-15 |
| PCA projection | 3.31970e-15 |
| Full-column ID | 0 |

Orthogonality checks passed. The PCA check forms the projection using the
returned rotation outside the timed call and verifies that `x` remains `NULL`.
The ID check also verifies the interpolation identity on the selected columns.

`Rprofmem` reports R-managed allocation events. Cumulative allocation is the
sum of their sizes; largest allocation is the largest individual recorded
event. Neither quantity is peak RSS or retained memory, and native-library
allocations may be absent. Peak resident memory was not measured. The input
fixtures already existed when profiling began, so their initial construction
is excluded from the allocation totals.

### Reproduction

Install the baseline checkout and candidate checkout into different libraries
using the same R executable and dependencies. Run the versions sequentially:

```sh
mkdir -p before-library after-library
R CMD INSTALL --library=before-library /path/to/baseline-checkout
R CMD INSTALL --library=after-library /path/to/candidate-checkout

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
Rscript tools/benchmark-lowrank.R \
  --library before-library --output before.csv --reps 3

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
Rscript tools/benchmark-lowrank.R \
  --library after-library --output after.csv --reps 3
```

The driver writes each timing sample, runtime versions, algorithm settings,
allocation statistics, and numerical checks to CSV. Raw allocation profiles go
into directories derived from the output filenames. `--cases` selects a subset
and `--profile-dir` changes the profile destination. R must have memory
profiling enabled for `Rprofmem`.

## Sparse boundaries and deferred work

Randomized SVD and QB preserve sparse inputs during multiplication by thin
dense factors. Their sketches, orthogonal bases, projected matrices, and
returned singular vectors still require dense storage. Randomized ID and CUR
also need dense reduced or selected factors. PCA's preprocessing and robust
PCA's matrix updates retain dense storage requirements.

Deterministic QB continues to use dense QR while constructing only the
requested columns of Q. Matrix sparse QR has different requirements and
permutation behavior, including a restriction to matrices with at least as
many rows as columns and special handling of structural rank deficiency; see
the [Matrix QR documentation](https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/qr-methods.html).
Using sparse input does not remove the dense storage requirement for this
deterministic branch.

An existing PCA summary inconsistency remains deferred: with `center=FALSE`,
eigenvalues come from the uncentered matrix, while `var` is calculated from
column sample variances. Consequently, reported proportions can exceed one
for data with a substantial mean. Changing that denominator changes established
summary semantics and needs a separate compatibility decision.

Release validation still needs an actual R 4.0 run, PDF manual generation,
and checks against downstream reverse dependencies. The local results above
do not cover Windows, macOS, or alternative BLAS implementations. Hosted CI
status is tracked separately in the draft pull request.
