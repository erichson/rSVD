# PCA, QR and robust-PCA follow-up review

This review accompanies rsvd 1.0.7 and compares the follow-up changes with
commit [`ad8f05a`](https://github.com/erichson/rSVD/commit/ad8f05a08491f4a75f12e8d7f14a524ce025d80d),
the first candidate in draft PR #17. Local validation was completed on
2026-10-17 UTC. The [initial report](README.md) used a
different baseline, the original 1.0.6 repository. Measurements in the two
reports should not be combined as if they shared a baseline.

## Compatibility decisions

The 25 existing function signatures, their argument order and defaults, and
the 10 exported names are unchanged. Public result component names and
classes are retained. The new sparse-PCA helpers and shared SVD computation
function are private. `methods`, supplied with R, is added to Imports for
explicit Matrix coercions; no additional external dependency is introduced.

The randomized SVD core was separated from its public method so PCA can
supply multiplication operators. Across 180 comparisons with the previous
implementation, both complete return objects and the final RNG states were
identical. These comparisons cover tall and wide real/complex matrices,
real sparse matrices, all three samplers, zero/two power iterations, and
omitted or truncated singular vectors. This is a same-runtime regression
check, not a promise of bitwise equality across BLAS implementations.
The [driver](rsvd-core-compatibility.R),
[per-case results](rsvd-core-compatibility-outputs.csv),
[formals comparison](rsvd-core-compatibility-formals.csv) and
[summary](rsvd-core-compatibility-summary.txt) are included.

Correcting a numerical bug necessarily changes the affected result. In
particular, uncentered PCA summaries and previously premature robust-PCA
stops are deliberately corrected. The public default sampling distribution,
oversampling, initial power count and robust-PCA iteration limit are retained.

## PCA semantics and sparse implementation

### `center = FALSE` and the variance denominator

R's [`prcomp()` documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/prcomp.html),
[`scale()` documentation](https://stat.ethz.ch/R-manual/R-devel/library/base/html/scale.html)
and [`prcomp()` source](https://svn.r-project.org/R/trunk/src/library/stats/R/prcomp.R)
establish that disabling centering skips subtraction of the mean. This works
for precentered input, but also permits an uncentered SVD. If scaling remains
enabled, the divisor is `sqrt(sum(x^2)/(m-1))`, not the centered sample
standard deviation. The package's pre-existing rule that scales below `1e-8`
become 1 is retained, including for constant columns.

The prior `rpca()` singular values used the processed matrix, but its total
`var` applied `stats::var()` to each processed column, implicitly centering
again. On a small uncentered example, the old denominator was 15 while the
processed squared energy divided by `m-1` was 104.6667. The first two reported
fractions summed to 6.96. The corrected denominator is
`sum(Mod(processed_A)^2)/(m-1)`, including the energy outside the requested
components. Already-centered input retains the same interpretation.

The fitted score convention remains `U D`. Prediction computes processed
observations times the loadings. On approximate fits these can differ by the
range-finder residual; the tests explicitly retain this distinction. Numeric
data frames now work in prediction when scaling is enabled.

### Implicit centering and scaling

For real sparse Matrix input, an explicit `k` and `rand = TRUE` select the
sparse path. Write the processed matrix as `M = C - 1 g^T`, where `C` is a
sparse, scaled working copy and `g` holds the scaled column offsets. Products
are evaluated as:

```
M X   = C X - 1 (g^T X)
M^T Y = C^T Y - g (1^T Y)
```

Only the thin products become dense. The shared range finder handles tall
and wide orientation through these operators. Prediction uses the same
representation. Deterministic PCA and `k = NULL` retain their dense SVD path.
This is PCA of sparse input, not a method imposing sparsity on the loadings.

Structural, logical and pattern representations are converted to a general
numeric compressed-column matrix without a full dense copy. This makes unit
diagonals explicit and combines duplicate triplet entries before working with
the documented `p` and `x` slots of
[`dgCMatrix`](https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/dgCMatrix-class.html).
The implementation does not assume that every stored entry is nonzero.

Column sums of squares use stored deviations plus the contribution from
implicit zeros. They avoid subtracting two large raw moments. Fully stored
columns can be centered explicitly without adding entries; this also reduces
cancellation in multiplication for large offsets with small variance. Their
means use base R's more accurate accumulation. Mean sums that overflow are
recomputed from the stored entries with an appropriate count ratio. Guards
avoid `0 * Inf` for fully stored large constants and recompute processed
energy when cached squared norms lead to `Inf / Inf` after scaling.

The implementation preserves the existing row-omission warning, names,
one-column dimensions, small-scale floor, and prediction's post-scaling NaN
policy. It does not mutate the caller's sparse input.

### Reference tests written before the production changes

The first 846 reference assertions produced 14 failing denominator checks
and one data-frame prediction error on both supported local runtimes. The
expanded PCA reference file now contains 1,157 assertions. It compares with
`prcomp()` and explicitly centered/scaled dense references, using low-rank
and decaying-spectrum inputs, and includes complex data, every center/scale
combination, compressed-column/row/triplet inputs, triangular unit diagonals,
duplicate triplets, logical/pattern matrices, constant and nearly constant
columns, extreme means, missing rows, `retx = FALSE`, and prediction.

Five allocation assertions also failed against the previous implementation:
each profiled sparse fit or prediction allocated at least one full input-sized
dense matrix. They pass with the implicit implementation. These are warmed
allocation guards, not timing tests, and skip on R builds without memory
profiling. They do not claim to measure native-library storage or peak RSS.

## Deterministic QR

The safe optimization is to avoid factoring columns that cannot affect the
requested basis. For real finite input, integer `k < min(m,n)`, and a
numerically independent leading `k`-column block, base R's default QR uses
the same first `k` Householder vectors whether it receives that block or the
whole matrix. Only the block is factored, and `qr.qy()` constructs the
requested columns of `Q`; projection still uses the complete input.

This follows the limited pivoting rule documented for
[`qr()`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/qr.html) and
implemented in [`dqrdc2`](https://svn.r-project.org/R/trunk/src/appl/dqrdc2.f).
Complex QR pivots across all columns, and a deficient real prefix may need
later columns, so both retain full dense QR. Full-rank requests and legacy
fractional-rank behavior also use the fallback. Separate minimum/maximum
checks preserve full QR's rejection of non-finite values in later columns
without a full dense logical mask.

The QR regression file has 529 assertions against full base QR and analytic
projectors. Coverage
includes tolerance-boundary prefixes, zero/dependent columns, tall/wide and
complex data, rank-one/full/oversized/fractional requests, sparse structural
formats, non-finite trailing columns, orthogonality, names, input immutability
and unchanged RNG state.

A general replacement with Matrix sparse QR is not included. Its documented
row/column permutations, `m >= n` requirement and structural-rank augmentation
are material compatibility differences; see the
[Matrix QR documentation](https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/qr-methods.html).
Executed counterexamples changed the selected column ordering relative to
dense pivoted QR and lost the desired orthonormal-column interpretation after
augmented rows were removed. Deterministic ID/CUR therefore retain their
existing pivoted backend. The deterministic QB basis is the base-QR basis,
not necessarily the best rank-k subspace.

## Robust PCA

### Norm, rank selection and stopping

The initialization now uses `norm(A, "M") / lambda`, the entrywise maximum
required by the dual of the entrywise l1 norm. The former `"I"` norm is a
maximum row sum, making the initial scaling depend on orientation. It was
a conservative bound, so this difference alone does not prove that the
solver could not converge; the premature stopping defect below is decisive.
Rank prediction and growth use `min(m,n)`. If every returned singular value
survives shrinkage, the randomized rank increases before that proximal step
is accepted; a sufficiently large request uses deterministic SVD.

These changes follow the formulation and rank guidance in
[Lin, Chen and Ma (2010)](https://arxiv.org/abs/1009.5055).
Checking only `||A-L-S||_F/||A||_F` can accept a feasible but nonoptimal split.
Closed-form constant and diagonal examples expose this defect even though
the reconstruction residual is essentially zero.

The existing sparse-then-low-rank update order is retained. With exact
proximal steps, the two subgradient conditions differ by
`mu * (L_new - L_old)`. The implementation requires both its normalized norm
and the existing reconstruction residual to be at most `tol`. It normalizes
the former by `max(1, ||Z||_F)` so input rescaling does not loosen convergence,
and increases the penalty only when this stationarity test is satisfied.
The role of primal and dual residuals is described in
[Boyd et al. (2011), section 3.3](https://stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf);
the normalization here is an implementation choice, not a claim to reproduce
all of that paper's absolute/relative tolerance conventions.

Randomized proximal accuracy can itself stall this test. After three
randomized steps without a 10% improvement, the implementation increases
the power count by two, at most twice, and then uses exact SVD if another
stall occurs. This bounded refinement is an engineering safeguard, not a
general convergence or optimality certificate. Likewise, a nonsaturated
approximate SVD is not a rigorous bound on an unobserved singular value.

The public `q` remains the initial power count. `maxiter` can still stop an
unconverged iteration. The returned `err` retains only the reconstruction
history; `trace = TRUE` also prints stationarity. Dense working arrays remain
necessary. Temporary threshold and shifted-input arrays are released earlier
to reduce their overlap with the new previous-iterate comparison.

### Numerical evidence and cost

The 178 new assertions and 12 existing robust-PCA assertions pass on both
local runtimes. They check transpose symmetry, saturated threshold ranks,
iteration limits, planted component recovery, objective values, output rank,
input scaling, and low-rank-only/sparse-only optima with known dual
certificates. A deliberately weak `p = 0, q = 0` starting approximation checks
that refinement can recover instead of stagnating.

The stronger stopping rule changes the amount of work. The benchmark compares
default calls and records component errors, objective, both available stopping
diagnostics, iterations and actual SVD backends. It is not a speed comparison
between equally accurate solutions. The controlled results are recorded below.

## Controlled performance measurements

All measurements below used R 4.6.1, Matrix 1.7-6, Linux x86_64 and the system
reference BLAS (`libblas.so.3.12.0`). Numerical jobs ran sequentially. Times
are medians of three repetitions; allocation sizes use MiB (1,048,576 bytes).
These are observations on selected fixtures, not general timing guarantees.

### PCA and SVD control

The sparse fixture has 20 disjoint rank-one blocks, 400,000 stored nonzeros
and dimensions 4,000 by 2,000 (5% density). Its transpose supplies the wide
case. Sparse PCA uses centering and scaling, `k = 20`, `p = 10`, `q = 2`;
the prediction model is fitted outside the measurement. The dense PCA case
is 6,000 by 600 with rank 40, no centering or scaling and no fitted scores.
The SVD control is 2,200 by 1,000 with rank 20.

| Case | Before time | After time | Before cumulative R allocation | After cumulative R allocation |
| --- | ---: | ---: | ---: | ---: |
| Sparse PCA, tall | 2.026 s | 0.088 s | 650.23 MiB | 59.13 MiB |
| Sparse PCA, wide | 2.205 s | 0.089 s | 710.92 MiB | 57.58 MiB |
| Sparse PCA, no scores | 2.265 s | 0.086 s | 647.77 MiB | 56.66 MiB |
| Sparse PCA prediction | 0.786 s | 0.008 s | 916.14 MiB | 7.16 MiB |
| Dense PCA, no scores | 1.205 s | 1.172 s | 178.13 MiB | 81.93 MiB |
| Dense SVD control | 0.343 s | 0.354 s | 20.02 MiB | 20.02 MiB |

The largest recorded R allocation for the sparse PCA cases fell from
64,000,048 bytes to 3,200,048 bytes: 61.04 to 3.05 MiB. Cumulative
allocations fell by about 91% for fitting and 99% for prediction on these
fixtures. Dense PCA saved about 54% of cumulative allocations by avoiding
column-wise variance copies, but its small timing difference does not
establish a speed improvement. The dense SVD control has unchanged allocation
volume and similar elapsed time; no core dense-SVD speedup is claimed.

All relative Frobenius residuals were below `4.4e-14`, and the numerical
orthogonality checks passed. The sparse PCA residuals changed only at ordinary
roundoff scale. Full settings, individual timing samples, allocation sizes
and numerical results are in [pca-benchmark-comparison.csv](pca-benchmark-comparison.csv).

### Guarded deterministic QR

All cases have rank 20 and independent leading columns. Sparse fixtures use
20 disjoint rank-one blocks and arrange the independent columns first so
the guarded path applies. The reference performs the full base QR; both
versions construct only the requested Q columns, making this a comparison
against the previous draft rather than the original full-Q implementation.

| Case | Full QR time | Guarded time | Full QR cumulative R allocation | Guarded cumulative R allocation |
| --- | ---: | ---: | ---: | ---: |
| Dense, 2,500 by 1,000 | 2.997 s | 0.052 s | 78.03 MiB | 3.60 MiB |
| Sparse, 2,000 by 1,000 | 2.181 s | 0.004 s | 77.73 MiB | 3.08 MiB |
| Sparse, 1,000 by 2,000 | 3.666 s | 0.003 s | 77.33 MiB | 2.00 MiB |

Both Q and B matched the full-QR reference exactly in every measured case.
Reconstruction residuals were at most `3.1e-15`. Millisecond-scale sparse
times are especially sensitive to timer resolution and system load; the
large savings arise because the selected fixtures permit skipping most of
the factorization. Cases requiring the fallback do not receive this gain.
See [qr-prefix-benchmark.csv](qr-prefix-benchmark.csv) for unrounded results.

### Robust-PCA costs

| Input and planted rank | Before time | After time | Before iterations | After iterations |
| --- | ---: | ---: | ---: | ---: |
| 80 by 60, rank 3 | 0.014 s | 0.064 s | 9 | 27 |
| 600 by 400, rank 5 | 0.384 s | 1.990 s | 8 | 37 |
| 1,200 by 120, rank 3 | 0.253 s | 1.201 s | 10 | 43 |

All three revised default calls satisfied both diagnostics at `tol = 1e-5`
within `maxiter = 50`. None used an exact-SVD update. The medium case refined
from `q = 2` to `q = 4`; the other cases retained `q = 2`. Runtime increased
by roughly 4-5 times, with correspondingly more cumulative allocation. For
the medium case, cumulative R allocations rose from 215.4 to 1,035.0 MiB,
while the largest allocation stayed at 1.83 MiB. Earlier temporary releases
do not offset the allocation volume of many more iterations. No robust-PCA
speedup or total-allocation reduction is claimed.

The stronger stopping rule is not a promise that each component error on
every already-successful fixture becomes smaller. At the default tolerance,
relative L errors changed from `7.74e-6` to `3.25e-5` on the small fixture,
`4.51e-5` to `1.17e-5` on the medium fixture, and `1.98e-6` to `1.18e-5`
on the tall fixture. Feasibility and stationarity have different meanings
from relative component recovery error. Both components and recovered ranks
are recorded in [rrpca-benchmark-summary.csv](rrpca-benchmark-summary.csv).

The decisive correctness evidence comes from analytic counterexamples. For
an all-ones 500 by 4 matrix with default `lambda`, the previous split had
relative L error 0.1111 and objective 49.6904, despite a reconstruction
residual around `1.8e-14`; the known optimum is 44.7214. The revised split
reaches that optimum with relative L error around `5.7e-15`. On a 10 by 10
identity with `lambda = 2`, the previous objective was 14.6667 instead of 10,
again with negligible reconstruction residual; the revised objective is 10.
The sparse-only identity case with `lambda = 0.5` remains correct. The
rank-20 single-update example separately verifies that surviving singular
directions are not discarded; its earlier iterate was infeasible and should
not be judged by objective alone. See the
[correctness results](rrpca-correctness-summary.csv) and
[cross-version test summary](rrpca-test-summary.csv).

### Final package validation

The complete installed-package test suite passes in both local runtimes
with both the system reference BLAS and the hosted Linux runner's OpenBLAS:

| R | Matrix | BLAS | Test blocks | Passing assertions | Failures | Errors | Warnings | Skips |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 4.6.1 | 1.7-6 | Reference 3.12.0 | 165 | 3,136 | 0 | 0 | 0 | 0 |
| 4.3.3 | 1.6-5 | Reference 3.12.0 | 165 | 3,136 | 0 | 0 | 0 | 0 |
| 4.6.1 | 1.7-6 | OpenBLAS 0.3.26 | 165 | 3,136 | 0 | 0 | 0 | 0 |
| 4.3.3 | 1.6-5 | OpenBLAS 0.3.26 | 165 | 3,136 | 0 | 0 | 0 | 0 |

The [test summary](full-test-summary.csv) records these totals. The suite
retains the earlier SVD/QB/ID/CUR, complex, sparse and API regressions as well
as the new PCA, QR and robust-PCA cases. Exact low-rank fixtures are checked
for reconstruction. Full-rank inputs with decaying spectra are checked
against the best rank-k tail error, not an unattainable zero residual.

The first hosted follow-up run exposed numerical defects in three test
references. On OpenBLAS 0.3.26, full base LINPACK QR of three independent
indicator columns followed by 21 identical columns repeatedly factors
roundoff until subnormal scaling produces non-finite entries. This was
reproduced without Matrix or rsvd. The full-QR comparison now uses a finite
reference, while a new 52-assertion test retains the original repeated-column
case and checks its known projector and reconstruction at ranks one, two
and three. On macOS, dense column means lost accuracy on a constant `1e200`
column and overflowed on partially stored `1e308` data. Those two sparse-PCA
tests now use the analytically known processed matrices, eigenvalues and
loadings. These corrections preserve the edge cases and tolerances and do
not change production code. The actual loaded BLAS library was verified for
each local configuration; changing `LD_LIBRARY_PATH` alone is insufficient
when R's startup library path overrides it.

The source package built successfully. The final local
`R CMD check --as-cran --no-manual` completed with **0 errors, 0 warnings and
2 NOTEs**. Installation, examples, documentation, code checks and all 3,136
installed-package assertions passed. The notes concern the maintainer email
change already present in the original repository and an installation-lock
directory retained by this execution environment. The ordinary staged
installation procedure is unchanged. The first restricted-network attempt
stopped at CRAN's incoming-index lookup; the reported complete run used
network access for that standard check.

Hosted CI checks current R/Matrix 1.7-6 on Linux, Windows and macOS, and R 4.3
with Matrix 1.6-5 on Linux. Each job asserts the actual Matrix version.
Final job links and statuses are tracked in
[draft PR #17](https://github.com/erichson/rSVD/pull/17).

## Recent randomized NLA: what to retain and what to defer

The existing stabilized Gaussian range finder remains a good default. The
[2023 RandNLA software monograph](https://arxiv.org/abs/2302.11474) explains
intermediate orthogonalization and pass-efficient range finding. The current
implementation uses `2*q + 2` products with the input or its adjoint: six at
the default `q = 2`. Avoiding dense preprocessing is therefore valuable on
sparse input without changing the algorithm. Existing `p` and `q` already
expose the principal accuracy/work tradeoff.

[Camaño, Epperly, Meyer and Tropp (2025)](https://arxiv.org/abs/2508.21189)
give a particularly relevant structured sparse-sign construction,
SparseStack. The developer prototype implements their Definition 1.7: four
independent equal-width hashing blocks, one signed entry per row in each
block, normalized by `1/sqrt(4)`. It supports widths divisible by four and
does not conflate this construction with a single CountSketch or SparseIID.
The paper's constant-four setting is supported empirically/conjecturally;
it is not a general high-probability theorem at every sketch width.

The isolated prototype completed 1,280 comparisons across 20 seeds,
coordinate/smooth bases, exact/geometric/slow/gapped spectra, widths 20/32,
and zero/two power iterations. Exact-rank residuals were below `1e-12`, and
the decaying-spectrum results were broadly comparable to Gaussian sketches.
However, a separate rank-two/width-four test failed to reconstruct in
60 of 500 SparseStack trials at `q = 0`, versus zero Gaussian failures
using the same error threshold. Neither method failed in the tested `q = 2`
trials, but numerical completion of a deficient basis is not a general cure
for lost sketch rank. Matrix 1.7-6 also rejected the tested complex sparse
construction and a complex-dense times real-sparse sketch product.

The scripts and compact summaries are retained so these observations are
reproducible. No new public sampler is added: small widths, arbitrary widths
and complex behavior need a coherent design before users can rely on it.
The kernel timing experiment measures sketch construction and the initial
product only, and cannot establish whole-decomposition speedups.
At width 80, observed combined sketch/product times were 0.100 to 0.010
seconds for dense tall input, 0.093 to 0.011 for dense wide input, 0.006 to
0.002 for sparse tall input, and 0.007 to 0.003 for sparse wide input
(Gaussian to SparseStack). Dense-tall cumulative R allocations slightly
increased, from 2.319 to 2.402 MiB. The mixed memory result and limited kernel
scope reinforce the need to evaluate the complete algorithm before exposing
an option. The [kernel measurements](sparsestack-kernel-benchmark.csv),
[accuracy summary](sparsestack-accuracy-summary.csv),
[small-width failures](sparsestack-small-width-summary.csv) and
[complex-support checks](sparsestack-complex-support.csv) accompany the
[experimental driver](sparsestack-experiment.R).

[Epperly and Webber (2026)](https://arxiv.org/abs/2605.19096) give sharp results
for randomized projections and sketched least squares. Their minimax result
concerns the untruncated estimator `projection(range(A*Omega))*A` under its
stated loss and assumptions; it does not establish universal optimality for
this package's truncated, power-iterated output. It supports caution about
discrete sketches on coherent, small-width problems. Orthogonalizing a
Gaussian sketch does not change its exact range and would add work here.

[Chenakkod and Dereziński (2026)](https://arxiv.org/abs/2605.09755) accelerate
power iteration through an intermediate sketch and a second sampling stage.
That requires additional sketch-size/error decisions and is a larger design
change than this maintenance revision. It is not used as a drop-in backend.

A future pass-budget option or residual-based accuracy diagnostic could fit
the existing design, but would need a separate return-value/API decision and
tests on slow spectra. In particular, subtracting nearly equal squared norms
to estimate a tiny residual is numerically unsafe. No new public option is
needed for the measured improvements in this revision.

## Reproduction and validation scope

Run from the repository root with separate installed baseline and candidate
libraries. Use the same R/Matrix/BLAS environment for comparisons, and keep
other numerical processes idle. The benchmark drivers use a warmup, three
elapsed-time runs and a separate `Rprofmem` run. RNG setup, garbage collection,
fixture construction and result validation are outside the timed operation.
Recorded R allocations exclude fixture construction, may exclude native
allocations, and are neither peak RSS nor retained memory.

```sh
Rscript tools/benchmark-lowrank.R --library before-library \
  --output pca-before.csv --reps 3 \
  --cases rpca_sparse_tall,rpca_sparse_wide,rpca_sparse_no_scores,rpca_sparse_predict,rpca_no_scores,rsvd_dense_tall
Rscript tools/benchmark-lowrank.R --library after-library \
  --output pca-after.csv --reps 3 \
  --cases rpca_sparse_tall,rpca_sparse_wide,rpca_sparse_no_scores,rpca_sparse_predict,rpca_no_scores,rsvd_dense_tall
Rscript tools/validation/qr-experiment.R after-library qr.csv
Rscript tools/validation/rrpca-experiment.R after-library robust.csv
Rscript tools/validation/sparsestack-experiment.R after-library sampler-results accuracy
Rscript tools/validation/sparsestack-experiment.R after-library sampler-timing timing
Rscript tools/validation/rsvd-core-compatibility.R baseline-checkout . compatibility
```

Set `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1` and `MKL_NUM_THREADS=1` for
these comparisons. The robust driver reads the reference function from
`ad8f05a`; an optional final argument can instead name a saved reference R
file. The sampler is a developer experiment and is not exported by rsvd.

The package retains R >= 4.0.0. Executed environments cover R 4.3.3 with
Matrix 1.6-5 and current R with Matrix 1.7-6; R 4.0 itself and a full
downstream reverse-dependency campaign remain outside this validation.
Checks use `--no-manual`, so they do not establish PDF-manual generation.
No release or merge is performed by this review.
