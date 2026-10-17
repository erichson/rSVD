Paired values are **deterministic / randomized**. **Speedup** is the ratio of the deterministic and randomized median elapsed times, computed from unrounded measurements. R allocations are cumulative allocation volume in GiB, not peak RSS.

## SVD

| Input | k | Time (s) | Speedup | R allocations (GiB) | Relative error |
| --- | ---: | ---: | ---: | ---: | ---: |
| Dense, 10,000 × 3,000 | 40 | 16.858 / 0.697 | **24.2×** | 1.234 / 0.133 | 9.67% / 9.68% |
| Dense, 20,000 × 6,000 | 40 | 97.157 / 2.144 | **45.3×** | 4.928 / 0.267 | 9.67% / 9.68% |

## PCA

| Input | k | Time (s) | Speedup | R allocations (GiB) | Relative error |
| --- | ---: | ---: | ---: | ---: | ---: |
| Dense PCA, 15,000 × 4,500 | 40 | 54.085 / 5.256 | **10.3×** | 8.056 / 2.981 | 10.3% / 10.3% |
| Sparse-input PCA, 10,000 × 3,000 | 40 | 18.336 / 0.463 | **39.6×** | 3.806 / 0.213 | 0.516% / 0.516% |
| Sparse-input PCA, 20,000 × 6,000 | 40 | 110.634 / 1.976 | **56.0×** | 15.212 / 0.482 | 0.548% / 0.548% |

## ID_CUR

| Input | k | Time (s) | Speedup | R allocations (GiB) | Relative error |
| --- | ---: | ---: | ---: | ---: | ---: |
| Column ID, dense, 15,000 × 4,500 | 40 | 43.727 / 0.684 | **63.9×** | 1.257 / 0.326 | 14.3% / 15.8% |
| CUR, dense, 15,000 × 4,500 | 40 | 42.141 / 1.076 | **39.2×** | 1.546 / 0.615 | 18.1% / 19.2% |
| Column ID, sparse, 15,000 × 4,500 | 40 | 46.938 / 0.287 | **163.5×** | 1.522 / 0.095 | 0.153% / 0.153% |
| CUR, sparse, 15,000 × 4,500 | 40 | 40.352 / 0.346 | **116.6×** | 1.572 / 0.145 | 0.189% / 0.189% |
