#' @title  Randomized robust principal component analysis (rrpca).
#
#' @description Robust principal components analysis separates a matrix into a low-rank plus sparse component.
#
#' @details
#' Robust principal component analysis (RPCA) is a method for the robust seperation of a
#' a rectangular \eqn{(m,n)} matrix \eqn{A} into a low-rank component \eqn{L} and a
#' sparse comonent \eqn{S}: 
#'
#' \deqn{A = L + S}
#'
#' To decompose the matrix, we use the inexact augmented Lagrange multiplier
#' method (IALM). The algorithm can be used in combination with either the randomized or deterministic SVD. 
#'
#'
#' @param A       array_like; \cr
#'                a real \eqn{(m, n)} input matrix (or data frame) to be decomposed. \cr
#'                na.omit is applied, if the data contain \eqn{NA}s.
#'
#' @param lambda  scalar, optional; \cr
#'                tuning parameter (default \eqn{lambda = max(m,n)^-0.5}).
#'
#' @param maxiter integer, optional; \cr
#'                maximum number of iterations (default \eqn{maxiter = 50}).
#'
#' @param tol     scalar, optional; \cr
#'                precision parameter (default \eqn{tol = 1.0e-5}).
#'
#' @param p       integer, optional; \cr
#'                oversampling parameter for \eqn{rsvd} (default \eqn{p=10}), see \code{\link{rsvd}}.
#'
#' @param q       integer, optional; \cr
#'                number of additional power iterations for \eqn{rsvd} (default \eqn{q=2}), see \code{\link{rsvd}}.
#'
#' @param trace   bool, optional; \cr
#'                print progress.
#'                
#' @param rand    bool, optional; \cr
#'                if (\eqn{TRUE}), the \eqn{rsvd} routine is used, otherwise \eqn{svd} is used.
#'
#'
#' @return \code{rrpca} returns a list containing the following components:
#' \describe{
#'    \item{L}{  array_like; \cr
#'              low-rank component; \eqn{(m, n)} dimensional array.
#'    }
#'    \item{S}{  array_like \cr
#'               sparse component; \eqn{(m, n)} dimensional array.
#'    }
#'}
#'
#' @author N. Benjamin Erichson, \email{erichson@berkeley.edu}
#'
#'
#' @references
#' \itemize{
#'  \item [1] N. B. Erichson, S. Voronin, S. L. Brunton and J. N. Kutz. 2019.
#'        Randomized Matrix Decompositions Using {R}. 
#'        Journal of Statistical Software, 89(11), 1-48.
#'         \doi{10.18637/jss.v089.i11}.
#' 
#'   \item  [2] Lin, Zhouchen, Minming Chen, and Yi Ma.
#'          "The augmented lagrange multiplier method for exact
#'          recovery of corrupted low-rank matrices." (2010).
#'          (available at arXiv \url{https://arxiv.org/abs/1009.5055}).
#' }
#'
#' @examples
#' library('rsvd')
#'
#' # Create toy video
#' # background frame
#' xy <- seq(-50, 50, length.out=100)
#' mgrid <- list( x=outer(xy*0,xy,FUN="+"), y=outer(xy,xy*0,FUN="+") )
#' bg <- 0.1*exp(sin(-mgrid$x**2-mgrid$y**2))
#' toyVideo <- matrix(rep(c(bg), 100), 100*100, 100)
#'
#' # add moving object
#' for(i in 1:90) {
#'   mobject <- matrix(0, 100, 100)
#'   mobject[i:(10+i), 45:55] <- 0.2
#'   toyVideo[,i] =  toyVideo[,i] + c( mobject )
#' }
#'
#' # Foreground/Background separation
#' out <- rrpca(toyVideo, trace=TRUE)
#'
#' # Display results of the seperation for the 10th frame
#' par(mfrow=c(1,4))
#' image(matrix(bg, ncol=100, nrow=100)) #true background
#' image(matrix(toyVideo[,10], ncol=100, nrow=100)) # frame
#' image(matrix(out$L[,10], ncol=100, nrow=100)) # seperated background
#' image(matrix(out$S[,10], ncol=100, nrow=100)) #seperated foreground

#' @export
rrpca <- function(A, lambda=NULL, maxiter=50, tol=1.0e-5, p=10, q=2, trace=FALSE, rand=TRUE) UseMethod("rrpca")

#' @export
rrpca.default <- function(A, lambda=NULL, maxiter=50, tol=1.0e-5, p=10, q=2, trace=FALSE, rand=TRUE) {
  #*************************************************************************
  #***        Author: N. Benjamin Erichson <nbe@st-andrews.ac.uk>        ***
  #***                              <2016>                               ***
  #***                       License: BSD 3 clause                       ***
  #*************************************************************************
  A <- as.matrix(A)
  m <- nrow(A)
  n <- ncol(A)

  rrpcaObj = list(L = NULL,
                  S = NULL,
                  err = NULL)


  # Set target rank
  k <- 1
  if(k > min(m,n)) rrpcaObj$k <- min(m,n)

  # Deal with missing values
  is.na(A) <- 0
  
  # Set lambda, gamma, rho
  if(is.null(lambda)) lambda <- max(m,n)**-0.5
  gamma <- 1.25
  rho <- 1.5

  if(rand == TRUE) {
    svdalg = 'rsvd'
  }else { 
    svdalg = 'svd' 
  }  
  
  # Compute matrix norms
  spectralNorm <- switch(svdalg,
                    svd = norm(A, "2"),
                    rsvd = rsvd(A, k=1, p=10, q=1, nu=0, nv=0)$d,
                    stop("Selected SVD algorithm is not supported!")
  )

  infNorm <- norm( A , "I") / lambda
  dualNorm <- max( spectralNorm , infNorm)
  froNorm <- norm( A , "F")

  # Initalize Lagrange multiplier
  Z <- A / dualNorm

  # Initialize tuning parameter
  mu <- gamma / spectralNorm
  mubar <- mu * 1e7
  mu <- min( mu * rho , mubar )
  muinv <- 1 / mu

  # Init low-rank and sparse matrix
  L = matrix(0, nrow = m, ncol = n)
  S = matrix(0, nrow = m, ncol = n)  
  
  niter <- 1
  err <- 1
  while(err > tol && niter <= maxiter) {

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Update S using soft-threshold
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      epsi = lambda / mu
      temp_S = A - L + Z / mu
      
      S = matrix(0, nrow = m, ncol = n)

      idxL <- which(temp_S < -epsi)
      idxH <- which(temp_S > epsi)
      S[idxL] <- temp_S[idxL] + epsi
      S[idxH] <- temp_S[idxH] - epsi

      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Singular Value Decomposition
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      R <- A - S + Z / mu

      if(svdalg == 'svd') svd_out <- svd(R)
      if(svdalg == 'rsvd') {
        if(k > min(m,n)/5 ) auto_svd = 'svd' else auto_svd = 'rsvd' 
      
        svd_out <- switch(auto_svd,
                          svd = svd(R),
                          rsvd = rsvd(R, k=k+10, p=p, q=q)) 
      }  
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Predict optimal rank and update
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      svp = sum(svd_out$d > 1/mu)
      
      if(svp <= k){
        k = min(svp + 1, n)
      } else {
        k = min(svp + round(0.05 * n), n)
      }

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Truncate SVD and update L
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # rrpcaObj$L =  svd_out$u[,1:rrpcaObj$k] %*% diag(svd_out$d[1:rrpcaObj$k] - muinv, nrow=rrpcaObj$k, ncol=rrpcaObj$k)  %*% t(svd_out$v[,1:rrpcaObj$k])
      L =  t(t(svd_out$u[,1:svp, drop=FALSE]) * (svd_out$d[1:svp] - 1/mu)) %*% t(svd_out$v[,1:svp, drop=FALSE])

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Compute error
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Astar = A - L - S
      Z = Z + Astar * mu

      err = norm( Astar , 'F') / froNorm
      rrpcaObj$err <- c(rrpcaObj$err, err)

      if(trace==TRUE){
        cat('\n', paste0('Iteration: ', niter ), paste0(' predicted rank = ', svp ), paste0(' target rank k = ', k ),  paste0(' Fro. error = ', rrpcaObj$err[niter] ))
        }

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Update mu
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      mu = min(mu * rho, mubar);
      muinv = 1 / mu

      niter =  niter + 1

  }# End while loop
  rrpcaObj$L <- L
  rrpcaObj$S <- S
  
  class(rrpcaObj) <- "rrpca"
  return( rrpcaObj )

}



