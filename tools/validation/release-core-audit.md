# Release audit: SVD, QB, ID, and CUR

This independent review began from commit
`eb301bf8e90043ef664b5c089c4ce9929dcedc0e` and examined the randomized range
finder, deterministic QB optimization, ID/CUR construction, conjugate-product
helpers, numerical-rank policy, and related tests. The review found one material
numerical issue that existing tests did not exercise: CUR with a requested rank
larger than the input's numerical rank. A focused private-helper fix and
independent regression tests address it. The remaining observations below do
not require a new algorithm or a public API change.

## Finding fixed: pseudoinverse rank detection

The pseudoinverse retained every singular value greater than
`d[1] * .Machine$double.eps`. Roundoff singular values of a rank-deficient
rectangular matrix can exceed that threshold. Inverting them gives enormous
connector coefficients and an inaccurate CUR reconstruction, even when the
selected columns and the ID interpolation are already accurate.

A deterministic reproduction on both R 4.6.1 / Matrix 1.7-6 and
R 4.3.3 / Matrix 1.6-5, with reference BLAS:

```r
set.seed(10)
A <- matrix(rnorm(30), 30, 1) %*% matrix(rnorm(18), 1, 18)
fit <- rcur(A, k = 3, rand = FALSE)
norm(A - fit$C %*% fit$U %*% fit$R, "F") / norm(A, "F")
```

Although `A` has rank one, the previous result had a relative residual of
`0.1872657331`. The selected row matrix had singular values approximately
`12.68684`, `4.15e-15`, and `2.16e-16`. The second value was only
`1.47 * .Machine$double.eps * d[1]`, yet was inverted. The connector's Frobenius
norm was approximately `2.14e14`. Directly checking the first Moore–Penrose
identity for the selected rows also showed a relative residual of about `0.075`.

The private helper now retains singular values satisfying

```r
d > d[1] * (max(dim(A)) * .Machine$double.eps)
```

This uses the dimension-scaled machine-precision criterion documented by
[Matrix::rankMatrix](https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/rankMatrix.html).
The dimension is multiplied by epsilon before multiplying by the leading
singular value. It avoids the substantially larger `sqrt(eps)` cutoff, which
would discard small but resolvable directions.

With this change, the reproduction's relative CUR residual is `1.85e-16` and
the connector norm is `0.1098`. The deficient pseudoinverse also satisfies the
independent analytic reference and Moore–Penrose identities. The original
single-epsilon threshold predates this maintenance revision; this is a legacy
numerical-rank defect exposed by the release review.

The cutoff changes results for singular directions near machine precision.
That is intentional numerical-rank handling. Public function signatures,
classes, fields, target-rank dimensions, sampling, and RNG usage are unchanged.

## Tests established before the fix

[`test-pseudoinverse-rank.R`](../../tests/testthat/test-pseudoinverse-rank.R)
adds four test blocks and 720 assertions. Its references are analytic maps,
projector identities, and exact reconstruction of independently formed
low-rank products; the tests do not duplicate the pseudoinverse implementation.

Coverage includes:

- The known pseudoinverse of a rank-one outer product and all four
  Moore–Penrose identities, in real and complex arithmetic, both orientations,
  and at scales `1e-100`, `1`, and `1e100`.
- A diagonal embedding with nonzero singular directions down to `1e-12`,
  including complex phases and the same scales. The full coordinate projector
  verifies that the small direction remains present.
- Rank-one and rank-three inputs with `k = rank + 2` and
  `k = min(dim(A))`, for tall/wide real/complex matrices, column and row ID,
  and deterministic/randomized CUR.
- Sparse padded rank-one inputs, scaling, result dimensions, selected-row and
  selected-column identities, and preservation of the caller's input.

| Runtime | BLAS | Before: failed / assertions | After: passed / assertions |
|---|---|---:|---:|
| R 4.6.1, Matrix 1.7-6 | Reference BLAS 3.12.0 | 19 / 720 | 720 / 720 |
| R 4.3.3, Matrix 1.6-5 | Reference BLAS 3.12.0 | 19 / 720 | 720 / 720 |
| R 4.6.1, Matrix 1.7-6 | OpenBLAS 0.3.26 | Not rerun | 720 / 720 |
| R 4.3.3, Matrix 1.6-5 | OpenBLAS 0.3.26 | Not rerun | 720 / 720 |

All after-fix runs had zero failures, errors, warnings, or skipped assertions.
The loaded OpenBLAS library path was checked explicitly. These are focused
validation results; the complete package and downstream release checks are
reported separately.

## Mathematical and behavioral review

### Randomized SVD and QB

The range finder computes `A %*% Omega`, stabilizes both sides of each power
iteration, forms the final orthonormal basis, and projects the input. The
conjugate-product helpers are consistent with complex adjoints. The SVD's wide
input branch reverses the orientation and requested vector counts consistently.
The PCA operator hook preserves these operations without requiring a centered
input copy.

SVD returns the requested truncated singular triplets. Randomized QB deliberately
retains the oversampled width, so its factor width generally equals `k + p`
after dimension clamping. These are different contracts and are correctly
distinguished in the updated documentation and tests. Tests based on spectral
tails, residuals, orthogonality, and repeated-singular-value subspaces are
appropriate; individual singular-vector entries are not reliable references.

### Deterministic QB

Factoring only an independent leading real block produces the same requested
Householder basis as full base QR. The rank-deficient-prefix, complex-input,
and full-width fallbacks preserve the necessary broader column selection.
Checking finiteness across the whole input prevents the optimization from
silently accepting non-finite entries outside the prefix.

This basis follows base QR column ordering. It is not in general the best
rank-`k` approximation, so it should not be used as an optimal low-rank reference
in performance claims.

### ID and CUR

ID's column-pivoted QR, permutation inversion, identity submatrix, and row-mode
adjoint are mathematically consistent. CUR selects rows from the column
skeleton and computes the connector from the column interpolation matrix and
the selected rows' pseudoinverse. The new rank-deficiency tests specifically
cover the distinction between the number of selected rows/columns and the
input's numerical rank.

The deterministic ID/CUR paths use dense LAPACK column-pivoted QR. They provide
appropriate same-output baselines for randomized ID/CUR. ID/CUR equations in
user-facing text should use approximation unless exact-rank recovery is being
discussed. Unqualified claims that ordinary pivoted QR always supplies a strongly
rank-revealing or optimally conditioned interpolation are stronger than this
implementation guarantees.

## Nonblocking legacy observations

Argument validation remains permissive and inconsistent outside documented
inputs. For example, `rsvd(diag(8), k = 3, p = -1)` can return an `NA` singular
value because the sketch has fewer columns than the requested output rank;
ID/CUR with that negative oversampling fail later during indexing. An invalid
`rid(mode = "bogus")` can return neither the column nor row factor. These are
not valid documented configurations. A coordinated future validation policy is
preferable to changing rounding, clamping, and error behavior piecemeal during
this maintenance release.

For benchmarks, specify valid integer ranks and nonnegative oversampling and
power counts, report achieved residuals, and include any required sparse-to-dense
conversion in the deterministic baseline. `base::svd(A, nu = k, nv = k)` saves
unrequested singular vectors but still computes the complete singular spectrum;
describe it accordingly rather than as a specialized partial-SVD algorithm.
