# README benchmarks

This directory contains the driver, raw measurements, and renderer for the
performance tables and figures in the top-level README. The comparisons use
deterministic references at the same requested rank. They do not compare with
an older rsvd version.

## Recorded experiment

The complete run on 2026-10-17 measured 19 fixtures, 44 fixture/method pairs,
and five timed repetitions per pair: 220 timed calls. Every method also had
one warmup and one separate allocation profile. No observations were removed.
Other local numerical jobs were paused during measurement.

The run used `--suite all`, so every case was measured before the final table
selection. The README tables now feature the already measured 8,000 by 2,000
SVD case and 6,400 by 1,600 sparse PCA case. Where present, the `table` flag in
the original CSV files records the initial selection and is retained unchanged;
all numeric measurements are also unchanged. A new `--suite tables` run uses
the current selection in `readme_cases()`.

| Setting | Recorded value |
| --- | --- |
| Numerical source | `b4c64a291aae9aa70a05206573e2d69c96704969`, rsvd 1.0.7 |
| Runtime | R 4.6.1, Matrix 1.7-6 |
| BLAS | OpenBLAS 0.3.26, one thread |
| Platform | Linux x86_64, virtualized environment |
| Hardware information | CPU model and physical RAM were not exposed by the environment |
| Synthetic fixture seed | `20261016` |
| Timed sampling seeds | `8101` through `8105` |
| Allocation sampling seed | `8101` |

The subsequent requirement, plotting-wrapper, and documentation changes do not alter any
decomposition used in these measurements. Exact package/BLAS paths and session
details are preserved in [metadata.csv](results/metadata.csv) and
[session-info.txt](results/session-info.txt). Absolute paths in those files
describe the measurement environment; use your own library and output paths
when rerunning.

## Baselines and output

| Family | Deterministic reference | Randomized call | Matched output |
| --- | --- | --- | --- |
| SVD | `svd(as.matrix(A), nu = k, nv = k)` | `rsvd(A, k, p = 10, q = 2)` | `k` left and right singular vectors and the leading singular values |
| PCA | `prcomp(as.matrix(A), rank. = k, center = TRUE, scale. = TRUE, retx = TRUE)` | `rpca(A, k, center = TRUE, scale = TRUE, retx = TRUE, p = 10, q = 2)` | `k` loadings and score columns, with the same preprocessing |
| Column ID | `rid(A, k, rand = FALSE)` | `rid(A, k, p = 10, q = 0)` | Selected columns and the interpolation factor |
| CUR | `rcur(A, k, rand = FALSE)` | `rcur(A, k, p = 10, q = 0)` | Selected columns, rows, and connector |

Base SVD and `prcomp()` compute the full singular spectrum even when fewer
vectors are requested. They are familiar deterministic references, not
specialized truncated-SVD competitors. There are no base-R ID/CUR functions;
those references are the package's dense LAPACK column-pivoted-QR paths.
Their selected columns or rows need not agree with randomized selections.

For sparse PCA, a third method calls `rpca(as.matrix(A), ...)` with otherwise
identical arguments. This separates the benefit of implicit preprocessing
from the benefit of the randomized decomposition. The dense conversion is
inside the measured call for both dense-input comparators. All methods
compute scores; the sparse method is not given a `retx = FALSE` advantage.

## Fixtures

The code in [benchmark-readme.R](benchmark-readme.R) is the authoritative
specification. Shapes, density, object size, and target rank are recorded in
[fixtures.csv](results/fixtures.csv).

**Dense synthetic data:** independent random orthonormal left and right
bases carry 80 geometrically spaced singular values from 1 to 0.01. Dense
Gaussian noise with expected Frobenius norm 0.01 makes the matrices full rank
with probability one. PCA adds nonzero column offsets before preprocessing.
The SVD scaling grid has 500, 1,000, 2,000, 4,000, and 8,000 rows, with one
quarter as many columns and `k = 40`. Separate cases cover wide and square
inputs. Dense ID, CUR, and PCA use 4,000 by 1,000 matrices.

**Sparse synthetic data:** twenty diagonal blocks each contain a rank-two
signal: a constant and a zero-mean slope on the left, random orthonormal
directions on the right, and two positive singular values. Gaussian noise
within each block produces a full-rank block while retaining exactly 5%
density. This represents data with clear low-dimensional structure and
nonzero column means; it is not intended as a universal sparse-data model.
The sparse PCA table uses 6,400 by 1,600, 4,000 by 2,000, and 2,000 by 4,000
matrices; sparse ID and CUR use 4,000 by 2,000. The PCA scaling grid has 800,
1,600, 3,200, and 6,400 rows, one quarter as many columns, and `k = 40`.
Noise and centering mean that rank-40 reconstruction is approximate.

**Image:** the package's 1,600 by 1,200 `tiger` grayscale image uses `k = 150`.
It provides a nonsynthetic spectrum and a less pronounced memory advantage.

These examples deliberately target compressible matrices. The unit suite
also checks difficult spectra and numerical edge cases. Performance curves
here do not claim to cover every spectral shape, sparsity pattern, rank,
hardware configuration, or application.

## Timing, allocations, and accuracy

The driver creates an input before timing, warms each method, and rotates
method positions across repetitions. It seeds each call separately and runs
explicit garbage collection before starting the timer. Fixture generation,
seed setup, explicit garbage collection, and numerical validation are excluded;
garbage collection triggered inside a call remains part of its elapsed time.
The README shows medians. The CSV also contains all five observations,
interquartile ranges, and extrema. Millisecond-scale measurements should not
be interpreted more precisely than their timer and run-to-run variation allow.

Allocation profiling is a separate, untimed run using `Rprofmem()`.
The driver sums numeric allocation records and retains their largest value
and count. It does not count unquantified `new page` markers. This is
**cumulative R allocation volume**, not peak RSS, a maximum-live-heap estimate,
or a complete accounting of allocations inside BLAS/LAPACK. Reused memory
can be counted repeatedly. The input fixture already exists; conversion,
preprocessing, temporary arrays, and returned factors created by the call
are included. See R's [Rprofmem documentation](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/Rprofmem.html).
Raw profiles are regenerated on every rerun; their compact numeric summaries
are committed here.

Accuracy is measured from an actual reconstruction in blocks of 128 rows,
outside timing and profiling. This avoids allocating a full reconstruction
just to validate a large input. SVD and PCA also check dimensions,
orthogonality, and the error relative to the exact optimal rank-`k` Frobenius
tail. The oracle comes from the full deterministic singular spectrum and is
checked for completeness and a finite positive tail. PCA errors refer to the
processed observations, not to the uncentered data.

ID checks the selected identity submatrix as well as reconstruction; CUR
checks factor dimensions and reconstruction. Neither comparison assumes that
pivoted QR is an optimal rank-`k` approximation. Dense ID/CUR results show
the measured increase in residual from randomized selection at `q = 0`.
Robust PCA is not included in these speed tables: its iterative convergence
and optimality checks are documented in the separate release audit.

## Reproduce the run

Install the candidate and dependencies into a separate R library. Set thread
limits **before** starting R and verify the BLAS reported by `sessionInfo()`.
Use a fresh output directory; the driver refuses to overwrite prior results.

```sh
mkdir -p benchmark-library
R CMD INSTALL --library=benchmark-library .

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  Rscript tools/benchmarks/benchmark-readme.R \
  --library benchmark-library --output benchmark-results \
  --suite all --reps 5 --source-ref YOUR_COMMIT \
  --machine 'Your CPU, RAM, operating system, and virtualization details'
```

The driver requires only base R, Matrix, rsvd, and an R build with memory
profiling enabled. `--suite smoke --reps 3` runs smaller functional checks;
it is not a substitute for the published performance run. `--suite tables`,
`--suite scaling`, and `--cases name,name` permit focused reruns. For timing,
run numerical work serially and avoid competing workloads.

Render a complete result set with Python 3, NumPy, and Matplotlib:

```sh
python tools/benchmarks/render-readme.py \
  --results benchmark-results --figures man/figures \
  --pdfs tools/benchmarks/figures --previews benchmark-previews \
  --update-readme README.md
```

The renderer reads measurements only; it does not run benchmarks or adjust
values. It replaces three marked table sections in the README and produces
two SVG/PDF figures. Curves use both shape and color, logarithmic axes, and
interquartile timing bands. SVG text is editable; PDF fonts are embedded.
PNG previews are used for visual inspection. Renderer versions are recorded
in [figure-environment.json](results/figure-environment.json).

## Files

- [summary.csv](results/summary.csv): 44 method summaries, including timings,
  allocations, approximation errors, and orthogonality.
- [timings.csv](results/timings.csv): all 220 timed observations and seeds.
- [tables.md](results/tables.md): generated Markdown tables.
- [SVD PDF](figures/benchmark-svd-scaling.pdf) and
  [PCA PDF](figures/benchmark-pca-scaling.pdf): vector figures for export.
- SVG figures live in `man/figures` and are displayed in the main README.

The source data retain unrounded numbers. Display rounding can make two
slightly different approximation errors look equal; use the CSV ratios when
assessing the small gap from optimal truncation.
