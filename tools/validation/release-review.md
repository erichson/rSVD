# Release review for rsvd 1.0.7

This review combines an independent mathematical audit, targeted numerical
experiments, regression tests, downstream-interface checks, and reproducible
performance measurements. Passing tests is evidence about the cases tested;
it is not a proof for every floating-point input or a guarantee about future
R, Matrix, or downstream releases.

## Assessment

The supported stack is **R 4.4 or later with Matrix 1.7-0 or later**.
The package requirements and the minimum-version CI job enforce this policy.

The reviewed range finder, implicit PCA operator, ID/CUR construction, and
robust-PCA updates are consistent with their intended algorithms. The second
release review found a material numerical-rank problem in CUR and valid-input
failures in the PCA plotting wrapper. Both received focused corrections and
tests established against the previous implementation. No new public
decomposition, sampling distribution, or solver option was added.

Public decomposition signatures, defaults, classes, and component names are
preserved. Some numerical results intentionally change where the old code was
incorrect, or where stronger robust-PCA convergence checks require more work.
An explicit comparison with the original 1.0.6 source confirms identical
argument lists and defaults for all 21 exported functions and registered S3
methods.
The supported full-width SVD contract is retained: real downstream packages
use it for Procrustes alignment even though the package's speed advantage is
primarily in low-rank computations.

Final execution results are recorded in the linked downstream and test
summaries below.

## Findings and conservative decisions

### SVD, QB, ID, and CUR

The SVD range finder stabilizes both sides of the power iteration and uses
the correct conjugate adjoints. Wide-matrix orientation, omitted singular
vectors, and the implicit-operator interface are consistent. Randomized QB
retains its oversampled basis; an SVD fit returns the requested truncated
triplets. Tests use residuals, singular values, orthogonality, and subspaces
instead of requiring arbitrary singular-vector signs to agree.

The deterministic QB optimization is restricted to a real, finite,
independent leading block whose Householder basis agrees with full base QR.
Complex, deficient-prefix, and boundary cases use the full fallback. A general
sparse-QR replacement was not introduced: pivoting, rank policy, and basis
ordering are observable behavior, and a different backend is not automatically
interchangeable.

The new CUR reproduction requested three columns/rows from an exactly
rank-one 30-by-18 matrix. The former single-epsilon pseudoinverse cutoff
inverted a roundoff singular direction. Its relative reconstruction residual
was **0.1873**, with a connector norm about **2.14e14**. A dimension-scaled
machine-precision cutoff reduced the residual to **1.85e-16**, with connector
norm **0.1098**, in the reference environment.

This cutoff follows the numerical-rank scale documented by
[Matrix::rankMatrix](https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/rankMatrix.html).
It is not the much larger `sqrt(eps)` cutoff. Analytic pseudoinverses,
Moore–Penrose identities, scaling tests, and directions down to `1e-12` verify
both removal of roundoff directions and retention of resolvable small ones.
The four new rank tests contain 720 assertions; the previous implementation
failed 19 on each of two R/Matrix environments.

Details: [core audit](release-core-audit.md).

### PCA and sparse products

The centered/scaled sparse operator and its adjoint apply the same matrix.
The review checked their bilinear identity, explicit dense references, wide
orientation, implicit zeros, fully stored columns with large means, and
complex dense PCA. In 64 additional small sparse configurations, the largest
forward/adjoint identity error was **1.78e-14**, and the largest processed-energy
error was **2.85e-14**. Complex rank-three reconstruction had relative error
at most **1.07e-15** in these probes.

`center = FALSE` skips centering, consistently with R's `prcomp()`/`scale()`
convention. Scaling without centering uses a root mean square with divisor
`m - 1`; it is not generally the centered standard deviation. `rpca()` keeps
its existing default `scale = TRUE` and small-scale floor. The total variance
includes all processed directions, including those not retained in the fit.

Several existing contracts are now explicit rather than silently changed:
prediction uses training column order; fitted scores use `U D`; prediction
projects onto the loadings and can differ by the randomized residual; sample
variance requires at least two observations; zero total energy gives undefined
explained proportions. PCA on sparse input does not impose sparse loadings.

Details and reproducible probes: [PCA audit](release-pca-audit.md).

### Robust PCA

The update order gives different feasibility and stationarity conditions.
The solver now checks both, controls penalty growth, expands saturated
singular-value thresholding, and refines stalled inexact steps before an
exact-SVD fallback. The existing `err` field remains a reconstruction-history
vector, and `maxiter` remains a cap on work.

The independent review also evaluated feasible primal upper bounds and dual
lower bounds, rather than trusting only the solver's own stopping quantities.
Twenty small deterministic fixtures had maximum normalized gap **2.68e-9**
at requested tolerance `1e-9`; six larger calls through the randomized solver
had maximum relative gap **6.16e-7** at tolerance `1e-7`. These runs deliberately
allowed more iterations than the default. They include deterministic SVD
updates when the solver's rank or refinement policy calls for them.

An approximate singular-value thresholding step does not provide a universal
certificate that no omitted singular value exceeds the threshold. Nor does
the default 50-iteration cap guarantee convergence. The documentation retains
these limits and explains that stronger checks can increase runtime.

Details, equations, and independent bounds: [PCA/robust-PCA audit](release-pca-audit.md).

### Plotting

The former generic wrapper always tried to draw two-component maps and sent
every extra argument to all helpers, despite their different signatures.
Consequently, a one-component fit, a fit without scores, or a valid
helper-specific named option could fail after partial output.

The wrapper now draws the available plots, routes exact named options to the
appropriate helpers, and keeps the usual plot order and invisible last-plot
return. Tests construct and render plots, check names and selected components,
and verify that variance proportions use total energy rather than only the
retained components. Documentation qualifies the correlation interpretation
of loadings and describes the existing 95% multivariate-t ellipse accurately.

## Validation evidence

The complete package suite passes **4,013 assertions in 181 test blocks** on
R 4.6.1 / Matrix 1.7-6 with both reference BLAS 3.12.0 and single-threaded
OpenBLAS 0.3.26. Each run has zero failures, errors, warnings, and skips.
The added sparse-product cases cover compressed-row, compressed-column,
coordinate, logical, and pattern storage, structural zeros, complex adjoints,
and custom matrix dispatch without changing the product helper.

- [Release test summary](release-test-summary.csv): actual test blocks,
  assertions, and outcomes by R/Matrix/BLAS environment.
- [Core audit](release-core-audit.md): tests established before the CUR and
  plotting corrections, plus algorithm and interface review.
- [PCA audit](release-pca-audit.md): independent sparse-operator, complex-PCA,
  and robust-PCA optimality probes with reproducible scripts and CSV results.
- [Downstream review](release-downstream-review.md): official source inventory,
  package checks, actual adapter probes, baseline comparisons, and coverage
  limits. An empty or disconnected downstream suite is not counted as proof
  that its rsvd integration was tested.
- [GitHub checks](https://github.com/erichson/rSVD/actions/workflows/R-CMD-check.yaml):
  Linux, Windows, and macOS package checks with an explicit loaded-Matrix
  version check, including R 4.4 / Matrix 1.7-0 at the supported minimum.

The regenerated 25-page PDF manual also builds successfully. Its pages were
rendered and checked for clipping. The conventional Rd layout retains some
long signature lines in margins and external-package PDF destinations that
cannot be resolved within this package's standalone manual.

The earlier [follow-up review](follow-up-review.md) retains the implementation
experiments, same-environment SVD output/RNG comparisons, and decisions about
more recent randomized-NLA options. Those historical benchmarks compare an
earlier draft; they are separate from the new README comparisons.

## Performance evidence

The README uses 19 fixtures and 220 timed calls against deterministic
references, with one OpenBLAS thread, warmups, rotated method positions, and
separate allocation profiling. Results include achieved residuals. Dense
fixtures have a geometric signal plus noise, sparse fixtures have noisy
blocks at 5% density, and the bundled tiger image supplies a nonsynthetic case.

Representative measured outcomes:

| Case | Deterministic time | Randomized time | Relative residual, deterministic / randomized |
| --- | ---: | ---: | ---: |
| Dense tall SVD, 8,000 by 2,000, rank 40 | 7.924 s | 0.364 s | 9.67% / 9.67% |
| Dense square SVD, 2,000 by 2,000, rank 40 | 3.123 s | 0.094 s | 9.67% / 9.68% |
| Sparse-input PCA, 6,400 by 1,600, rank 40 | 4.818 s | 0.144 s | 0.495% / 0.495% |
| Sparse-input PCA, 4,000 by 2,000, rank 40 | 5.765 s | 0.115 s | 0.487% / 0.487% |
| Dense column ID, 4,000 by 1,000, rank 40 | 0.502 s | 0.040 s | 14.9% / 16.4% |
| Dense CUR, 4,000 by 1,000, rank 40 | 0.556 s | 0.055 s | 19.2% / 19.9% |

The small 500-by-125 SVD case took 0.005 s with either method and allocated
more with rsvd. It remains in the scaling plot and raw measurements; the
tables emphasize larger matrices. ID/CUR retain their measured
accuracy tradeoff. The sparse-PCA figure includes a dense-input randomized
control to distinguish algorithmic speed from sparse-storage savings.

Base `svd()`/`prcomp()` compute the full singular spectrum even when only
`k` vectors are requested; these are not specialized partial-SVD baselines.
ID/CUR use the package's deterministic pivoted-QR paths because base R has
no corresponding routines. Memory figures report **cumulative R allocation**,
not peak RSS or complete native-library memory. Hardware details unavailable
in the environment are identified rather than guessed.

Methods, raw observations, and editable/exportable figures:
[benchmark guide](../benchmarks/README.md).

## Remaining limits

The validation covers specific package versions, platforms, matrices, and
downstream branches. Full application workflows with optional external
services, executables, or large training datasets are not all exercised.
Numerical-rank decisions near floating-point precision are inherently sensitive
to scale, dimensions, and perturbations. Approximation quality depends on the
spectrum and chosen rank, and performance depends on the BLAS and workload.

Some permissive behavior for malformed arguments and minor legacy plot
highlighting policies is intentionally left unchanged; the core audit names
those cases. This maintenance release does not redesign validation policy or
make exact-vector reproducibility promises. Release confidence should be
based on the concrete evidence above and maintained through future CI and
downstream checks.
