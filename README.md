[![Build Status](https://travis-ci.org/Benli11/rSVD.svg?branch=master)](https://travis-ci.org/Benli11/rSVD)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/rsvd)](http://cran.r-project.org/package=rsvd)

Randomized Singular Value Decomposition using R
***********************************************
Version: 0.5

Depends: R (>= 3.2.2)

License: GPL (>= 2)

Date: 2016-06-15

Author: N. Benjamin Erichson { nbe () st-andrews.ac.uk }

Randomized singular value decomposition (rsvd) is a very fast
probabilistic algorithm to compute the near optimal low-rank singular value
decomposition of massive data sets with high accuracy. SVD plays a central role
in data analysis and scientific computing. SVD is also widely used for computing
(randomized) principal component analysis (PCA), a linear dimensionality reduction technique.
Randomized PCA (rpca) is using the approximated singular value decomposition
to compute the most significant principal components. The package includes also an
function for computing (randomized) robust principal component analysis (RPCA).
In addition several plot functions are provided.


Get started
*************
Install the rsvd packages via CRAN
```R
install.packages("rsvd")
```

Run a motivating example: Image compression
```R
library(rsvd)
data(tiger)

# Image compression using randomized SVD
s <- rsvd(tiger, k=150)
tiger.re = s$u %*% diag(s$d) %*% t(s$v) # reconstruct image

# Display orginal and reconstrucuted image
par(mfrow=c(1,2))
image(tiger, col = gray((0:255)/255))
image(tiger.re, col = gray((0:255)/255))
```
Here are the results:
![tiger](https://raw.githubusercontent.com/Benli11/data/master/img/reTiger.png)

and the speedup gained over the base SVD function:

```R
library(microbenchmark)

timing_svd <- microbenchmark(
  'SVD' = svd(tiger, nu=150, nv=150),
  'rSVD' = rsvd(tiger, k=150),
  times=50)

print(timing_svd, unit='s')
```
![timing](https://raw.githubusercontent.com/Benli11/data/master/img/timeing.png)