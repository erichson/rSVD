[![Build Status](https://travis-ci.org/Benli11/rSVD.svg?branch=master)](https://travis-ci.org/Benli11/rSVD)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/rsvd)](http://cran.r-project.org/package=rsvd)

Fast Randomized Singular Value Decomposition using R
****************************************************
Randomized singular value decomposition (rsvd) is a very fast
probabilistic algorithm that can be used to compute the near optimal low-rank singular value
decomposition of massive data sets with high accuracy. SVD plays a central role
in data analysis and scientific computing. SVD is also widely used for computing
(randomized) principal component analysis (PCA), a linear dimensionality reduction technique.
Randomized PCA (rpca) uses the approximated singular value decomposition
to compute the most significant principal components. This package also includes a
function to compute (randomized) robust principal component analysis (RPCA).
In addition several plot functions are provided.


Get started
*************
Install the rsvd package via CRAN
```R
install.packages("rsvd")
```

You can also install the development version from GitHub using [devtools](https://cran.r-project.org/package=devtools):

```r
devtools::install_github("benli11/rsvd")
```

SVD example: Image compression
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


References
*************
* [N. Benjamin Erichson, Sergey Voronin, Steven L. Brunton, J. Nathan Kutz. “Randomized Matrix Decompositions using R.” (2016)](http://arxiv.org/abs/1608.02148)
* [Sergey Voronin, Per-Gunnar Martinsson. “RSVDPACK: Subroutines for computing partial singular value decompositions via randomized sampling on single core, multi core, and GPU architectures.” (2015)](https://arxiv.org/abs/1502.05366)
* [Nathan Halko, Per-Gunnar Martinsson, Joel A. Tropp. “Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions.” (2011)](https://arxiv.org/abs/0909.4061)