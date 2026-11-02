# November review of rsvd 2.0.0

This pass starts from `759ef01b4e196486ec3b308c11eac97218e987f9` and reviews
all package R code, public help topics, the README, metadata, tests, and the
source-package workflow. It adds focused fixes to the existing 2.0.0
candidate. Earlier validation and large-matrix benchmarks remain historical
evidence for the source revisions named in their reports.

## Correctness fixes

| Case | Previous behavior | Revised behavior |
| --- | --- | --- |
| Negative oversampling, for example `rsvd(diag(8), k = 3, p = -1)` | Could construct fewer sampled directions than requested and return an `NA` singular value. | Rejects invalid sketch controls before sampling. |
| Fractional or malformed power-iteration counts | Could silently run a different number of iterations or fail inside a low-level operation. | Requires a finite non-negative integer. |
| Invalid ID selection mode | Could return an object with neither the requested columns nor rows. | Requires `"column"` or `"row"`. |
| Fractional PCA, ID, and CUR ranks | Could produce inconsistent factor dimensions. | Truncates a positive fractional rank before constructing factors. |
| Fewer than two complete PCA observations | Divided by zero in sample-variance calculations. | Reports that at least two complete observations are required. |
| Prediction with the wrong variable count | Could recycle preprocessing vectors before failing in matrix multiplication. | Checks the training variable count before preprocessing. |
| Robust PCA with `tol >= 1` | Could return untouched zero components and no measured residual. | Performs an update before applying the tolerance. |
| Invalid robust-PCA controls or complex input | Failed late or entered an unsupported real soft-thresholding calculation. | Checks controls and rejects complex input explicitly. |
| Complex QB with non-finite data | LAPACK QR could return non-finite factors without an error. | Checks the thin complex factors and reports the failure. |
| Fractional plot component indices | Could plot a truncated index under the wrong component label. | Requires finite integer component indices. |
| Short plot groups or labels | Could silently recycle entries to different observations or variables. | Requires one entry per relevant row. |
| Non-positive or fractional `top.n` | Could select an unintended variable subset. | Requires a positive integer. |

Public function signatures, defaults, result names, and valid-input RNG
behavior are retained. Missing training rows still follow each routine's
documented omission policy. Deterministic routines ignore unused randomized
controls. Prediction retains its existing missing-value behavior, including
replacing post-scaling `NaN` entries with zero; the help page now describes it.

## Performance decisions

The review checked memory use in sparse matrix products and PCA preprocessing,
the thin QR path, and repeated work in robust PCA. Sparse-input PCA continues
to center and scale through its product operators rather than allocating the
full processed matrix. No replacement sparse QR backend or new factorization
algorithm was introduced.

An early validation draft scanned every dense matrix before randomized SVD
and QB. A bounded comparison exposed a substantial cost for small requested
ranks, so that extra scan was removed from those paths. Scalar, shape, and type
checks remain; the factorization checks non-finite products. Deterministic QB
still checks the entire matrix when its leading-column shortcut needs that
protection. PCA and randomized ID also rely on the factorization's checks;
CUR avoids duplicating the input check performed by its column ID. PCA and
robust PCA use the private computation engine after their own checks,
avoiding repeated wrapper validation inside iterative work.

The README timing tables retain their original measured source and date.
They are not presented as new measurements of this candidate. The release
bundle contains the focused before/after comparison from this pass.

## Documentation

The README now opens with a runnable example and a function-selection table.
It explains what reconstruction error measures, when to change `k`, and when
oversampling or power iterations can help. The mathematical descriptions,
sparse-memory qualifications, and limitations of randomized approximations
remain in place.

The DESCRIPTION, NEWS, dataset attribution, and plotting help were edited
using direct language and concrete explanations. The citation contact now
matches the maintainer contact in DESCRIPTION. The optional stale Date field
was removed. Function help is regenerated from the roxygen comments.

## Verification and scope

The added tests are `test-release-core-review.R`,
`test-release-numerical-review.R`, and `test-release-plots-review.R`. They run
with the existing numerical, sparse, complex, memory, and plotting tests in
the ordinary package check. Regression cases were also run against the
starting source to distinguish actual failures from new expected behavior.

Independent probes compare reconstruction and variance identities with
analytic or dense references. They cover small tall and wide shapes,
one-variable matrices, rank-deficient data, real and complex matrices,
sparse storage, deterministic and randomized paths, and robust-PCA input
scales from `1e-300` to `1e300`. Paired comparisons also check outputs and RNG
state for supported SVD/QB calls.

The release bundle records the final source archive checksum, exact R and
dependency versions, test totals, complete `R CMD check --as-cran` output,
and hosted CI results. These identify the tested artifact; the version number
alone does not distinguish successive candidate archives.

This pass reviews the recorded downstream call sites and rechecks public
interfaces. Earlier complete reverse-dependency campaigns are described in
`release-downstream-review.md`; those historical results should not be read
as a new full downstream campaign on every candidate. The final report
identifies which checks were actually rerun.

Passing checks supports the cases tested. It does not prove correctness for
every floating-point matrix. Approximation error still depends on the
spectrum and chosen rank. Robust PCA's iteration cap remains a work limit,
and its randomized thresholding safeguard is not a general optimality
certificate.
