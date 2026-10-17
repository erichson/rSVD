# Downstream compatibility review for rsvd 1.0.7

This review checks the interfaces that current downstream packages actually use.
It supplements rsvd's numerical regression suite; a successful downstream
package check is useful evidence, but it does not establish that the downstream
package tested its rsvd integration.

## Source inventory

The inventory was retrieved on **2026-10-17 at 03:07:38 UTC** from the official
[CRAN package metadata](https://cran.r-project.org/src/contrib/PACKAGES.gz).
The downloaded metadata has MD5 `962eb9fbfb99c15e04877e2981f16344`.
There are **10 direct Imports users and 6 direct Suggests users**, with no direct
Depends, LinkingTo, or Enhances entries in this snapshot. All 16 official source
archives were downloaded and verified against their metadata MD5 checksums.
The exact versions, dates, checksums, and source URLs are retained in
[release-reverse-dependency-inventory.csv](release-reverse-dependency-inventory.csv).
R's parser identified 20 production call sites across these sources; the
file names, line numbers, and expressions are retained in
[release-downstream-call-sites.csv](release-downstream-call-sites.csv).
Flashier accepts the resulting SVD object without directly calling `rsvd()`.

CRAN's published rsvd version in this snapshot is **1.0.5**. That is the
downstream regression baseline; it is distinct from the repository's former
1.0.6 version. The README's performance comparisons use deterministic methods.

| Package | Version | Relationship | Interface found in current source |
|---|---|---|---|
| FPCdpca | 0.4.0 | Imports | Positional `rsvd(X, m, p, q)`; consumes `v` |
| harf | 0.1.0 | Imports | Centered, scaled `rpca()`; consumes fitted scores |
| jackstraw | 1.3.21 | Imports | Repeated `rsvd()` calls; consumes right singular vectors |
| MUGS | 0.1.0 | Imports | Two SVD initializers and consumption of `u`, `d`, `v` |
| nlpembeds | 1.0.0 | Imports | SVD-based embeddings; consumes `u`, `d`, `v` |
| sparsepca | 0.1.2 | Imports | Randomized `rqb(X, k=l, p=0, q=2)`; consumes `Q`, `B` |
| SPECK | 1.0.1 | Imports | SVD rank selection and reduced-rank reconstruction |
| TCA | 1.2.1 | Imports | Centered, scaled `rpca()`; consumes `rotation` and `x` |
| text2map | 0.4.0 | Imports | Full-width Procrustes SVD and truncated embedding SVD |
| wordvector | 0.6.5 | Imports | Sparse SVD with `nu=0`, `nv=k` |
| eigencore | 1.0.3 | Suggests | Public `rsvd()` as a validation/benchmark adapter |
| flashier | 1.0.7 | Suggests | Accepts an SVD result as a low-rank data representation |
| LSX | 1.5.2 | Suggests | Sparse SVD with `nu=0`, `nv=k`, optionally cached |
| multivarious | 0.3.2 | Suggests | Public SVD wrapper forwarding `k`, `q`, `p` |
| Seurat | 5.6.0 | Suggests | SVD scores for the FIt-SNE PCA initializer |
| stm | 1.3.8 | Suggests | Centered, unscaled `rpca()` scores for spectral initialization |

The source scan found no calls to `rid()`, `rcur()`, `id()`, `cur()`, `rrpca()`,
or private rsvd helpers in these packages. The sole `rqb()` consumer leaves
`rand=TRUE`, so the deterministic QB width correction has no detected direct
exposure. The `rpca()` consumers use scores/loadings with centering enabled;
none of the inspected call sites consumes its corrected `$var` field.

Two details matter for maintenance:

- `text2map` omits `k` in its Procrustes call and constructs `V U^T` from the
  full-width result. `multivarious` also permits `k=min(dim(X))`. These are real
  compatibility requirements even though rsvd's performance advantage concerns
  low-rank work.
- Seurat's term "RPCA" normally refers to reciprocal PCA. Its detected use of
  this package is the separate `rsvd()` call in the FIt-SNE initializer.

## Bioconductor adapter

The review additionally includes the official
[BiocSingular 1.28.1 release](https://bioconductor.org/packages/3.23/bioc/html/BiocSingular.html).
Its source archive MD5 is `e81f52e88e8fceee9a3362c63ecad3f5`, matching the
Bioconductor 3.23 metadata. This is a focused adapter review, not a complete
Bioconductor reverse-dependency campaign.

`runRandomSVD()` calls the public `rsvd()` interface on `DelayedMatrix` and
`ScaledMatrix` objects. Its `safe_rsvd()` wrapper restores zero-column matrices
when `nu=0` or `nv=0`. `runPCA(..., BSPARAM=RandomParam())` uses this path and
constructs PCA scores from `U D`. Preserving custom matrix multiplication and
transpose methods is therefore important; coercing every input to a base matrix
would undermine this adapter's deferred computation.

## Validation method

The test environment is R 4.6.1, Matrix 1.7-6, Ubuntu 24.04 x86_64, and the
reference BLAS 3.12.0. This review concerns the modern supported environment;
the release minimum is R 4.4 and Matrix 1.7-0. Candidate and CRAN rsvd installations are isolated from
each other and from the shared runtime. All 16 inventoried CRAN packages load
at the exact versions above. Their runtime dependency closure was installed
from CRAN-derived Posit Ubuntu binaries, with source builds where needed.
Private Ubuntu GLPK/SuiteSparse libraries satisfy Seurat's igraph dependency;
the private OpenBLAS library satisfies harf's float dependency. R itself remains
on the reference BLAS.

The tested numerical core is commit `b4c64a2`, with the working plotting fix
included in the candidate installation. The exact installed source hashes are
recorded in [release-downstream-source-manifest.json](release-downstream-source-manifest.json).
Subsequent edits to the R/Matrix minimum and plot presentation do not change
the numerical paths covered here.

Source-package checks use `R CMD check --no-manual --no-vignettes` with
`_R_CHECK_FORCE_SUGGESTS_=false`. Thus the examples and connected tests are
executed, but PDF manuals, vignettes, and unavailable optional integrations are
outside the check scope. Standard `donttest`/`dontrun` exclusions and the
downstream tests' own skip conditions remain active. Failures are compared with the published rsvd 1.0.5
under the same R/Matrix/dependency environment before assigning a regression
to the candidate.

The focused probes exercise the actual installed downstream functions or
unchanged expressions extracted from their bodies. The latter are explicitly
labeled: they cover the rsvd-facing branches of harf, MUGS, Seurat, and stm
without claiming to run forest training, the complete MUGS optimizer, an
external FIt-SNE executable, or t-SNE. Other probes execute whole public or
internal adapter functions. Analytic low-rank reconstruction, orthogonality,
dimensions, finite values, names, normalized embeddings, and projection
identities provide checks beyond successful return alone.

### Limitations in existing downstream tests

- FPCdpca's only connected unit test checks `2 * 2 == 4`. Its package examples
  do exercise the distributed SVD, and this review adds an explicit adapter
  probe with finite results and the expected result contract.
- sparsepca and flashier ship testthat files without a `tests/testthat.R`
  runner. Their ordinary package checks do not execute those files.
- sparsepca's shipped tests still use the obsolete `library(spca)` name. An
  unchanged run and a compatibility run with only that library name corrected
  are reported separately. Existing numerical assertions are not weakened.
- Several sparsepca reconstruction assertions use a signed sum of residuals,
  which can conceal cancellation. The added adapter probe independently checks
  a Frobenius reconstruction error and PCA standard deviations on an analytic
  low-rank fixture.

## Execution results

**No downstream regression was found in this campaign.** All 16 CRAN packages
completed source checks with **zero errors and zero warnings** after the local
dependency/build environment was completed. Fifteen checks report `Status: OK`.
TCA reports two notes that are identical with the released rsvd 1.0.5 baseline:
unused declared imports (`Matrix`, `quadprog`) and malformed braces in its own
`tca.Rd` documentation.

### Source-package checks

| Package | Final check | Connected assertions passed | Declared test skips |
|---|---|---:|---:|
| FPCdpca | OK | 1 | 0 |
| harf | OK | 6 | 0 |
| jackstraw | OK | 861 | 0 |
| MUGS | OK | 68 | 0 |
| nlpembeds | OK | 14 | 0 |
| sparsepca | OK | No connected runner | — |
| SPECK | OK | No unit tests | — |
| TCA | 2 baseline-matched notes | No unit tests | — |
| text2map | OK | 1,008 | 24 |
| wordvector | OK | 195 | 16 |
| eigencore | OK | 3,811 | 46 |
| flashier | OK | No connected runner | — |
| LSX | OK | 292 | 11 |
| multivarious | OK | 381 | 0 |
| Seurat | OK | 623 | 13 |
| stm | OK | 9 | 0 |

These connected suites contain **7,269 passing assertions and 110 declared
skips**. This is a count of the downstream packages' checks, not a claim that
all 7,269 assertions directly exercise rsvd. In particular, FPCdpca's assertion
is the placeholder noted above. Seurat's skipped tests concern absent optional
`ape`, `BPCells`, `hdf5r`, `leidenbase`, and `sf` integrations, plus four tests
that explicitly skip on CRAN. The other packages also retain their own skip
conditions.

The final selection is in
[release-downstream-check-summary.csv](release-downstream-check-summary.csv).
[release-downstream-check-results.csv](release-downstream-check-results.csv)
retains every completed attempt, including the resolved setup failures and
TCA's baseline check. Initial build failures for eigencore, MUGS, Seurat, and
stm came from missing development linker names in this minimal container.
Initial example/test failures for flashier, multivarious, and stm came from
optional dependencies used without a skip guard. Supplying the private linker
aliases and the required small packages resolved those failures. No downstream
production code was changed.

### Explicit execution of otherwise missed tests

| Package / test scope | Candidate 1.0.7 | CRAN baseline 1.0.5 |
|---|---:|---:|
| sparsepca, shipped files unchanged | 0 assertions; all 3 files skipped | 0 assertions; all 3 files skipped |
| sparsepca, only `library(spca)` renamed in a separate copy | 25 passed, 0 failed | 25 passed, 0 failed |
| flashier, all 7 shipped testthat files | 44 passed, 0 failed | Not rerun; its focused input adapter was compared below |
| BiocSingular, `test-random-svd.R` | 61 passed, 0 failed | 61 passed, 0 failed |

The sparsepca compatibility copy changes exactly three library-name statements;
all existing assertions remain intact. Its 19 test blocks and flashier's 22
test blocks provide coverage that a plain `R CMD check` would have missed.
BiocSingular's eight randomized-SVD test blocks cover vector counts, zero-column
outputs, thin/wide matrices, centering/scaling, deferred computation, parallel
execution, names, and invalid inputs. Detailed totals are in
[release-downstream-testthat-results.csv](release-downstream-testthat-results.csv).

### Focused numerical and interface probes

**All 30 probes pass with both candidate 1.0.7 and released 1.0.5.** All 30
paired result comparisons agree under `all.equal(..., tolerance=1e-8)`, with
attributes checked. The probes include 19 CRAN cases covering all 16 packages,
and 11 additional BiocSingular cases. Results and explicitly labeled execution
scope are in [release-downstream-probe-results.csv](release-downstream-probe-results.csv);
paired comparisons are in
[release-downstream-baseline-comparison.csv](release-downstream-baseline-comparison.csv).

The fixtures are small matrices with known low rank, including a rank-four
matrix with singular values 20, 10, 5, and 2, a positive semidefinite matrix,
and sparse block count data. Selected workflow probes add small noise. Tests
use reconstruction errors, orthogonality, projection identities, PCA standard
deviations, dimensions, names, finite values, and normalized embeddings as
appropriate. Comparisons use reconstructions or Gram matrices where singular
vector signs should not affect the result.

The CRAN cases cover positional arguments, omitted singular vectors, tall and
wide sparse inputs, the randomized QB path in sparsepca, full-width Procrustes
alignment, truncated retrofitting, PCA scores/loadings, resampling, and low-rank
input conversion. The BiocSingular cases cover dense, sparse, delayed dense,
and delayed sparse representations with both deferred settings, plus omitted
vectors, the wide-matrix cross-product adapter, and `RandomParam` PCA.

This gives substantial evidence for the currently used public interfaces and
their expected numerical behavior. It does not test every transitive downstream
package or optional workflow. No HDF5-backed data set, external FIt-SNE process,
full harf training run, complete MUGS optimizer, or full Bioconductor
reverse-dependency campaign is claimed. The isolated environment's 313 available
package versions are recorded in
[release-downstream-environment.csv](release-downstream-environment.csv).

## Reproduction

The scripts intentionally keep downloaded sources, installations, and check
logs outside the package source tree:

```sh
Rscript tools/validation/release-downstream-inventory.R "$WORK"
R CMD INSTALL --library="$CANDIDATE_LIB" rsvd_1.0.7.tar.gz
R CMD INSTALL --library="$BASELINE_LIB" rsvd_1.0.5.tar.gz
Rscript tools/validation/release-downstream-install.R "$WORK" "$CANDIDATE_LIB"
Rscript tools/validation/release-downstream-bioconductor.R "$WORK" "$CANDIDATE_LIB"
Rscript tools/validation/release-downstream-check.R "$WORK" "$CANDIDATE_LIB" candidate ALL
Rscript tools/validation/release-downstream-check.R "$WORK" "$BASELINE_LIB" baseline-review TCA
Rscript tools/validation/release-downstream-probes.R "$WORK" "$CANDIDATE_LIB" candidate
Rscript tools/validation/release-downstream-probes.R "$WORK" "$BASELINE_LIB" baseline
Rscript tools/validation/release-downstream-testthat.R "$WORK" "$CANDIDATE_LIB" candidate-unchanged sparsepca
Rscript tools/validation/release-downstream-testthat.R "$WORK" "$CANDIDATE_LIB" candidate-name-fix sparsepca --fix-spca-name
Rscript tools/validation/release-downstream-testthat.R "$WORK" "$CANDIDATE_LIB" candidate flashier
Rscript tools/validation/release-downstream-testthat.R "$WORK" "$CANDIDATE_LIB" candidate BiocSingular --filter=random-svd
Rscript tools/validation/release-downstream-testthat.R "$WORK" "$BASELINE_LIB" baseline-unchanged sparsepca
Rscript tools/validation/release-downstream-testthat.R "$WORK" "$BASELINE_LIB" baseline-name-fix sparsepca --fix-spca-name
Rscript tools/validation/release-downstream-testthat.R "$WORK" "$BASELINE_LIB" baseline BiocSingular --filter=random-svd
Rscript tools/validation/release-downstream-summarize.R "$WORK" tools/validation
```

The Bioconductor installer uses the official 3.23 metadata, verifies each
source archive checksum, and installs dependencies serially. Exact versions
are in [release-bioconductor-dependency-versions.csv](release-bioconductor-dependency-versions.csv).
The probe script includes its BiocSingular cases only when that package is
available; count the resulting rows to confirm that this part ran.

Use a normal R development installation with BLAS/LAPACK/Fortran development
libraries. In this minimal container, private linker aliases to the already
installed reference libraries were needed for source builds. Private GLPK,
SuiteSparse, and OpenBLAS runtime libraries were also added for the downstream
dependencies, without modifying system or shared R libraries. The install
script includes the small optional packages found necessary for unguarded
examples/tests (`RcppML`, `bench`, `reformulas`, `GPArotation`, `tm`, `SnowballC`,
and `wordcloud`), while other optional dependency closures remain excluded.

Do not edit a driver while it is running. The first two long-running driver
logs ended with an orchestration parse error after all their package results
had already been written, because the driver file was edited to add the local
linker workaround. The corrected driver parses successfully; later invocations
used frozen copies. All individual `00check.log` files and package exit statuses
were preserved and consolidated. This orchestration issue did not trigger any
package reruns or count as a downstream failure.
