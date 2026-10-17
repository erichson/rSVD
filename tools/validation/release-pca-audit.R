# Independent PCA operator and robust-PCA release audit.
# Usage: Rscript tools/validation/release-pca-audit.R LIBRARY OUTPUT_PREFIX
args <- commandArgs(TRUE)
library(rsvd, lib.loc = args[1L])
library(Matrix)
set.seed(48271)
fn <- function(name) getFromNamespace(name, 'rsvd')
frob <- function(A) sqrt(sum(Mod(A)^2))
results <- list()
record <- function(label, value) {
  results[[length(results) + 1L]] <<- data.frame(label = label, value = value)
}
max_forward <- max_adjoint <- max_inner <- max_energy <- 0
for (shape in list(c(2,1), c(7,11), c(11,7), c(20,12))) {
  for (density in c(0,0.1,0.7,1)) {
    A <- matrix(rnorm(prod(shape)),shape[1],shape[2])
    A[matrix(runif(prod(shape)),shape[1],shape[2]) > density] <- 0
    for (center in c(FALSE,TRUE)) for (scale in c(FALSE,TRUE)) {
      sparse <- fn('rpca_sparse_numeric')(Matrix(A,sparse=TRUE))
      ctr <- if (center) fn('rpca_sparse_means')(sparse) else FALSE
      norms <- fn('rpca_sparse_norms')(sparse,ctr)
      scl <- if (scale) sqrt(norms/(nrow(A)-1)) else FALSE
      if (scale) scl[scl<1e-8] <- 1
      operators <- fn('rpca_sparse_products')(sparse,ctr,scl)
      dense <- if(center) sweep(A,2,ctr,'-') else A
      if(scale) dense <- sweep(dense,2,scl,'/')
      X <- matrix(rnorm(ncol(A)*3),ncol(A),3)
      Y <- matrix(rnorm(nrow(A)*3),nrow(A),3)
      max_forward <- max(max_forward,frob(operators$forward(X)-dense%*%X))
      max_adjoint <- max(max_adjoint,frob(operators$adjoint(Y)-t(dense)%*%Y))
      max_inner <- max(max_inner,abs(sum(operators$forward(X)*Y)-sum(X*operators$adjoint(Y))))
      max_energy <- max(max_energy,abs(sum(operators$squared_norms())-sum(dense^2)))
    }
  }
}
record('maximum sparse forward absolute error',max_forward)
record('maximum sparse adjoint absolute error',max_adjoint)
record('maximum sparse inner-product identity absolute error',max_inner)
record('maximum sparse processed squared-energy absolute error',max_energy)

# Dense complex PCA, arbitrary complex bases and wide/tall shapes.
max_complex_eigen <- max_complex_reconstruction <- max_complex_prediction <- 0
for(shape in list(c(24,12),c(12,24))) {
  U <- qr.Q(qr(matrix(rnorm(shape[1]*3)+1i*rnorm(shape[1]*3),shape[1],3)))
  V <- qr.Q(qr(matrix(rnorm(shape[2]*3)+1i*rnorm(shape[2]*3),shape[2],3)))
  A <- U %*% (c(20,8,2)*Conj(t(V)))
  for(center in c(FALSE,TRUE)) for(scale in c(FALSE,TRUE)) {
    d <- if(center) sweep(A,2,colMeans(A),'-') else A
    if(scale) d <- sweep(d,2,sqrt(colSums(Mod(d)^2)/(nrow(A)-1)),'/')
    fit <- rpca(A,k=3,center=center,scale=scale,p=3,q=2)
    ref <- svd(d)
    max_complex_eigen <- max(max_complex_eigen,max(abs(fit$eigvals-ref$d[1:3]^2/(nrow(A)-1))))
    max_complex_reconstruction <- max(max_complex_reconstruction,frob(d-fit$x%*%Conj(t(fit$rotation)))/frob(d))
    max_complex_prediction <- max(max_complex_prediction,frob(predict(fit,A)-d%*%fit$rotation))
  }
}
record('maximum complex PCA eigenvalue absolute error',max_complex_eigen)
record('maximum complex PCA rank-three reconstruction relative error',max_complex_reconstruction)
record('maximum complex PCA prediction absolute error',max_complex_prediction)

# Observe primal/dual state in a private copy; production methods remain unchanged.
# The proximal subgradient conditions permit an independent numerical dual bound.
solver <- fn('rrpca.default')
b <- body(solver)
b[[length(b)]] <- quote(list(L=L,S=S,err=rrpcaObj$err,Z=Z,dualErr=dualErr,iterations=niter-1L))
body(solver) <- b
max_gap <- max_dual_violation <- max_primal <- 0
for(j in 1:20) {
  A <- matrix(rnorm(24),6,4)
  lambda <- c(0.2,0.5,0.9,1.3)[(j-1)%%4+1]
  z <- solver(A,lambda=lambda,rand=FALSE,tol=1e-9,maxiter=3000)
  spectral <- svd(z$Z,nu=0,nv=0)$d[1]
  violation <- max(spectral,max(abs(z$Z))/lambda)
  dual_bound <- sum(A*z$Z)/max(1,violation)
  # L and A-L form an exactly feasible decomposition for a primal upper bound.
  primal_bound <- sum(svd(z$L,nu=0,nv=0)$d)+lambda*sum(abs(A-z$L))
  gap <- (primal_bound-dual_bound)/max(1,abs(primal_bound))
  max_gap <- max(max_gap,gap)
  max_dual_violation <- max(max_dual_violation,max(0,violation-1))
  max_primal <- max(max_primal,tail(z$err,1))
}
record('maximum deterministic robust-PCA normalized primal-dual gap, 20 fixtures',max_gap)
record('maximum deterministic robust-PCA dual infeasibility, 20 fixtures',max_dual_violation)
record('maximum deterministic robust-PCA relative reconstruction residual, 20 fixtures',max_primal)

# Sparse prediction: preserve dense NA/NaN behavior when loadings contain zeros.
fit <- structure(list(rotation=diag(3),center=c(1,2,3),scale=c(1,2,3)),class='rpca')
newdata <- matrix(c(1,NA,0,NaN,0,1,0,0,0),3,3)
a <- predict(fit,newdata)
b <- predict(fit,Matrix(newdata,sparse=TRUE))
record('sparse prediction with zero loadings equals dense including NA masks',isTRUE(all.equal(a,b)))

write.csv(do.call(rbind,results),paste0(args[2L],'.csv'),row.names=FALSE)
print(do.call(rbind,results),row.names=FALSE)
cat('\nSession:\n')
print(sessionInfo())
