## SVD

| Input | k | Time (s), deterministic / randomized | R allocations (GiB), deterministic / randomized | Relative error, deterministic / randomized |
| --- | ---: | ---: | ---: | ---: |
| Dense, 10,000 × 3,000 | 40 | 16.858 / 0.697 | 1.234 / 0.133 | 9.67% / 9.68% |
| Dense, 20,000 × 6,000 | 40 | 97.157 / 2.144 | 4.928 / 0.267 | 9.67% / 9.68% |

## PCA

| Input | k | Time (s), deterministic / randomized | R allocations (GiB), deterministic / randomized | Relative error, deterministic / randomized |
| --- | ---: | ---: | ---: | ---: |
| Dense PCA, 15,000 × 4,500 | 40 | 54.085 / 5.256 | 8.056 / 2.981 | 10.3% / 10.3% |
| Sparse-input PCA, 10,000 × 3,000 | 40 | 18.336 / 0.463 | 3.806 / 0.213 | 0.516% / 0.516% |
| Sparse-input PCA, 20,000 × 6,000 | 40 | 110.634 / 1.976 | 15.212 / 0.482 | 0.548% / 0.548% |

## ID_CUR

| Input | k | Time (s), deterministic / randomized | R allocations (GiB), deterministic / randomized | Relative error, deterministic / randomized |
| --- | ---: | ---: | ---: | ---: |
| Column ID, dense, 15,000 × 4,500 | 40 | 43.727 / 0.684 | 1.257 / 0.326 | 14.3% / 15.8% |
| CUR, dense, 15,000 × 4,500 | 40 | 42.141 / 1.076 | 1.546 / 0.615 | 18.1% / 19.2% |
| Column ID, sparse, 15,000 × 4,500 | 40 | 46.938 / 0.287 | 1.522 / 0.095 | 0.153% / 0.153% |
| CUR, sparse, 15,000 × 4,500 | 40 | 40.352 / 0.346 | 1.572 / 0.145 | 0.189% / 0.189% |
