# Independent release audit: PCA and robust PCA

Reviewed production commit: `eb301bf8e90043ef664b5c089c4ce9929dcedc0e`.

This is a second, independent review of `R/rpca.R`, `R/rpca_sparse.R`,
`R/rrpca.R`, the operator interface in `R/rsvd.R`, and their associated tests
and documentation. The review starts from the matrix operations and proximal
optimality conditions. Passing the existing test suite is supporting evidence,
rather than the definition of correct behavior.

**Assessment:** no new release-blocking defect was found in these paths for
their documented input domains. The implicit PCA operations are consistent
adjoints; real and complex PCA use the intended processed matrix; robust PCA
uses the correct shrinkage and stationarity quantities for its update order.
The limited-iteration and inexact-SVD qualifications in the robust-PCA
documentation should remain visible. No production runtime code was changed
as a result of this review.

## PCA algebra and numerical behavior

For the real sparse path, the working representation is

\[
M=C-\mathbf{1}g^\mathsf{T}.
\]

The forward operation `C %*% X - 1 %*% (t(g) %*% X)` and adjoint operation
`t(C) %*% Y - g %*% colSums(Y)` are the two adjoint products of this same
matrix. The projected matrix is the transpose of the adjoint product. For
wide inputs, swapping forward and adjoint and transposing the projected
operation gives the corresponding operations on the transposed matrix.
The code follows these identities.

Sparse column normalization first makes structural entries explicit and
converts to a numeric general compressed-column representation. This is
necessary for the subsequent slot-based moment calculations: a stored
row/column pair occurs once, including after duplicate triplets are combined,
and an implicit unit diagonal participates in the moments.

The centered squared column norm is computed as a sum over stored deviations
plus the contribution from absent entries. It does not subtract two large raw
moments. Fully stored columns can be centered before multiplication without
increasing storage; this avoids cancellation of a large offset in a matrix
product. Partially stored columns retain their mean correction implicitly.
The `squared_norms` closure applies the same representation after scaling, so
the fallback total energy describes the matrix actually decomposed.

`center = FALSE` leaves the mean in the decomposition. In that case the
reported total is a second moment using the package's divisor `m - 1`.
Scaling uses squared moduli for complex input and squared entries for real
input. The established small-scale floor is retained. The variance summary
includes discarded singular directions, so changing `k` does not change its
denominator.

Dense complex products use conjugate adjoints through the existing helper
functions. The complex mean, real-valued column scales, eigenvalues, and
predictions are mutually consistent. Sparse complex Matrix storage is not
advertised as supported.

Fitted scores retain `U %*% D`; prediction applies the trained transformation
and loadings. For a tall randomized decomposition, the difference between
these quantities is the residual of the sampled range. Requiring them to be
identical would test a different estimator. The existing tests correctly
check this distinction and compare reconstruction or projectors rather than
arbitrary vector signs.

Missing training rows, prediction's existing post-scaling NaN replacement,
ordinary NA propagation, single-observation prediction dimensions, and
sparse input preservation are covered in the current tests. An additional
probe verified NA/NaN propagation when the loading matrix itself contains
zeros.

## Additional small numerical probes

The following probes ran with reference BLAS 3.12.0 on both:

- R 4.6.1 and Matrix 1.7-6;
- R 4.3.3 and Matrix 1.6-5.

There were 64 sparse operator cases: four small shapes (including a single
column and a wide matrix), four densities including zero and fully stored,
and all four center/scale combinations. The comparison matrix was explicitly
processed dense data. Complex PCA used independent random complex bases in
tall and wide rank-three matrices, under each center/scale combination.

| Probe | Maximum observed error, both environments |
|---|---:|
| Sparse forward product, absolute Frobenius error | `7.35e-15` |
| Sparse adjoint product, absolute Frobenius error | `6.36e-15` |
| Sparse forward/adjoint inner-product identity, absolute error | `1.78e-14` |
| Sparse processed squared energy, absolute error | `2.85e-14` |
| Complex PCA eigenvalues, absolute error | `3.56e-14` |
| Complex rank-three PCA reconstruction, relative Frobenius error | `1.07e-15` |
| Complex prediction versus explicit processed product | `0` |

These values describe small independent probes, not worst-case bounds or
performance measurements. The existing tests supply broader structural,
extreme-value, approximation-tail, and platform coverage.

The reproducible driver is [release-pca-audit.R](release-pca-audit.R).
Compact results are in [release-pca-audit-r46.csv](release-pca-audit-r46.csv)
and [release-pca-audit-r43.csv](release-pca-audit-r43.csv).

## Robust PCA: proximal conditions and independent bounds

The solver updates the sparse component before the low-rank component.
Writing `Z_old` for the multiplier at the start of a step, the exact sparse
proximal update gives

\[
Z_{old}+\mu(A-L_{old}-S_{new})
\in\lambda\partial\|S_{new}\|_1.
\]

With exact singular-value thresholding, the low-rank update and subsequent
multiplier update give

\[
Z_{new}=Z_{old}+\mu(A-L_{new}-S_{new})
\in\partial\|L_{new}\|_*.
\]

The difference between these two subgradient candidates is
`mu * (L_new - L_old)`. Thus feasibility and the implemented stationarity
quantity address different optimality conditions. Dividing the latter by
`max(1, ||Z_new||_F)` preserves invariance under nonzero rescaling of the
input, subject to floating-point effects. Controlled penalty increases avoid
reducing changes in the primal variables solely by raising the penalty.

The entrywise maximum in the initialization is the appropriate dual norm of
the entrywise one-norm. Bounding the rank by `min(m, n)` is necessary for
rectangular inputs. Expanding a saturated sampled thresholding update
addresses insufficient requested rank, while the bounded power refinement
and eventual deterministic fallback address observed stagnation. A
nonsaturated approximate SVD can still underestimate an omitted singular
value; neither the rank check nor the stopping tests prove an exact proximal
update. The existing documentation states this limitation appropriately.

For an independent numerical check, a private copy of the production solver
returned its final multiplier as well as its existing internal diagnostics.
Only the final return expression was changed in that copy. Neither the
production solver nor the iteration formulas were modified.

For each result, let

\[
c=\max\left(1,\|Z\|_2,\frac{\|Z\|_{\max}}{\lambda}\right),
\qquad Y=Z/c.
\]

Then `Y` satisfies the dual constraints to numerical precision. The dual
lower bound is `sum(A * Y)`. The exactly feasible decomposition
`(L, A - L)` gives the primal upper bound

\[
\|L\|_*+\lambda\|A-L\|_1.
\]

These are independently evaluated spectral and entrywise norms, not a
restatement of the solver's convergence diagnostic. As floating-point
calculations, they are numerical certificates rather than interval-arithmetic
proofs.

Twenty arbitrary real `6 x 4` fixtures used deterministic SVD, four lambda
values, `tol = 1e-9`, and a maximum of 3,000 iterations. In both R/Matrix
environments, the largest normalized primal/dual gap was **`2.68e-9`** and
the largest returned reconstruction residual was **`9.88e-10`**. The gap was
normalized by `max(1, abs(primal_upper_bound))`.

Six `60 x 40` fixtures used `rand = TRUE`, `tol = 1e-7`, and a maximum of
1,000 iterations on R 4.6.1 / Matrix 1.7-6. Their largest gap divided by the
positive primal upper bound was **`6.16e-7`**. All completed with both
diagnostics at most the requested tolerance. These calls exercised the public
randomized solver, including its rank-based deterministic SVD path when
applicable; they should not be described as using randomized thresholding at
every iteration. The CSV's `stagnation_exact_fallback` flag refers only to
the additional fallback selected after repeated accuracy refinement.

These probes deliberately allow more work than the default iteration cap;
they establish solution behavior when the requested tolerances are reached.
They do not establish that 50 iterations are enough for every input.

The randomized driver and results are in
[release-pca-randomized-audit.R](release-pca-randomized-audit.R) and
[release-pca-randomized-audit.csv](release-pca-randomized-audit.csv).

## Longstanding edge contracts made explicit

Two user-facing details merit documentation without a late change to
established runtime behavior:

1. **Prediction uses column positions.** `predict.rpca()` expects the same
   variables in the same order as training input; it does not match names.
   A small named-column permutation changed predictions by `53.78025`, while
   `predict.prcomp()` reordered those names and returned the original scores.
   R's [current `prcomp` source](https://svn.r-project.org/R/trunk/src/library/stats/R/prcomp.R)
   makes that name matching explicit. The `rpca` documentation now states its
   positional contract.
2. **Sample variance needs two complete observations.** The package uses
   `m - 1`; a one-row uncentered unscaled fit can consequently report infinite
   eigenvalues. For a processed zero matrix, the zero eigenvalues are valid,
   but the zero-over-zero explained-variance proportions are undefined. The
   documentation now distinguishes these cases. No new validation or change
   to returned values was introduced.

## Reproduction

Use a library containing the reviewed package and its dependencies:

```sh
Rscript tools/validation/release-pca-audit.R LIBRARY audit-results
Rscript tools/validation/release-pca-randomized-audit.R LIBRARY randomized-results.csv
```

The first command writes `audit-results.csv` and prints session details.
These scripts perform small accuracy probes, not controlled timing or memory
benchmarks.

The mathematical references checked in this review were Lin, Chen, and Ma's
[augmented-Lagrangian robust-PCA formulation](https://arxiv.org/pdf/1009.5055)
and the primal/dual stopping discussion in Boyd et al.'s
[ADMM review](https://stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf).
