## SVD

| Input | k | Time (s), deterministic / randomized | R allocations (MiB), deterministic / randomized | Relative error, deterministic / randomized |
| --- | ---: | ---: | ---: | ---: |
| Dense, large, 8,000 × 2,000 | 40 | 7.924 / 0.364 | 644.75 / 105.28 | 9.67% / 9.67% |
| Dense, tall, 4,000 × 1,000 | 40 | 0.981 / 0.098 | 162.16 / 52.71 | 9.67% / 9.68% |
| Dense, wide, 1,000 × 4,000 | 40 | 1.129 / 0.111 | 163.08 / 83.23 | 9.67% / 9.68% |
| Dense, square, 2,000 × 2,000 | 40 | 3.123 / 0.094 | 246.17 / 41.65 | 9.67% / 9.68% |
| Tiger image, 1,600 × 1,200 | 150 | 0.755 / 0.186 | 107.24 / 95.18 | 10.2% / 10.3% |

## PCA

| Input | k | Time (s), deterministic / randomized | R allocations (MiB), deterministic / randomized | Relative error, deterministic / randomized |
| --- | ---: | ---: | ---: | ---: |
| Dense PCA, 4,000 × 1,000 | 40 | 1.147 / 0.202 | 483.00 / 224.58 | 10.2% / 10.2% |
| Sparse PCA, large, 6,400 × 1,600 | 40 | 4.818 / 0.144 | 1,312.34 / 127.00 | 0.495% / 0.495% |
| Sparse PCA, tall, 4,000 × 2,000 | 40 | 5.765 / 0.115 | 1,102.10 / 94.31 | 0.487% / 0.487% |
| Sparse PCA, wide, 2,000 × 4,000 | 40 | 5.825 / 0.110 | 1,103.53 / 91.23 | 0.511% / 0.511% |

## ID_CUR

| Input | k | Time (s), deterministic / randomized | R allocations (MiB), deterministic / randomized | Relative error, deterministic / randomized |
| --- | ---: | ---: | ---: | ---: |
| Column ID, dense, 4,000 × 1,000 | 40 | 0.502 / 0.040 | 73.82 / 34.89 | 14.9% / 16.4% |
| CUR, dense, 4,000 × 1,000 | 40 | 0.556 / 0.055 | 98.81 / 59.88 | 19.2% / 19.9% |
| Column ID, sparse, 4,000 × 2,000 | 40 | 2.663 / 0.033 | 227.95 / 28.01 | 0.152% / 0.152% |
| CUR, sparse, 4,000 × 2,000 | 40 | 2.739 / 0.042 | 241.71 / 41.77 | 0.189% / 0.189% |
