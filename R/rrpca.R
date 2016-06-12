#' @title  Randomized robust principal component analysis (rRPCA).
#
#' @description Robust principal components analysis using randomized singular value decomposition.
#
#' @details
#' Robust principal component analysis is via PCP using the IALM (Lin et. al) algorithm.
#'
#'
#' @param A       array_like \cr
#'                a numeric input matrix (or data frame), with dimensions \eqn{(m, n)}. \cr
#'                If the data contain \eqn{NA}s na.omit is applied.
#'
#' @param k       int, optional \cr
#'                determines the number of principle components to compute. It is required that \eqn{k} is smaller or equal to
#'                \eqn{n}, but it is recommended that \eqn{k << min(m,n)}.
#'
#' @param lamb    real, optional \cr
#'                tuning paramter (default \eqn{gamma=m^-0.5}).
#'
#' @param gamma   real, optional \cr
#'                tuning paramter (default \eqn{gamma=1.25}).
#'
#' @param rho     real, optional \cr
#'                tuning paramter (default \eqn{rho=1.5}).
#'
#' @param maxiter int, optional \cr
#'                determines the maximal numbers of iterations.
#'
#' @param tol     real, optional \cr
#'                tolarance paramter for the desired convergence of the algorithm.
#'
#' @param svdalg  str c('auto', 'rsvd', 'svd'), optional \cr
#'                Determines which algorithm should be used for computing the singular value decomposition.
#'                By default 'auto' is used, which decides whether to use \code{\link{rsvd}} or \code{\link{svd}},
#'                depending on the number of principle components. If \eqn{k < min(n,m)/1.5} randomized svd is used.
#'
#' @param p       int, optional \cr
#'                oversampling parameter for \eqn{rsvd}  (default \eqn{p=5}), see \code{\link{rsvd}}.
#'
#' @param q       int, optional \cr
#'                number of power iterations  for \eqn{rsvd} (default \eqn{q=2}), see \code{\link{rsvd}}.
#'
#' @param trace   bool, optional \cr
#'                print progress.
#'
#' @param ...     arguments passed to or from other methods, see \code{\link{rsvd}}.
#'
#' @param ................. .
#'
#' @return \code{rpca} returns a list with class \eqn{rpca} containing the following components:
#'    \item{L}{  array_like \cr
#'              Low-rank component, array of shape \eqn{(m, n)}.
#'    }
#'    \item{S}{  array_like \cr
#'               Sparse component, array of shape \eqn{(m, n)}.
#'    }
#'
#'    \item{k}{  int \cr
#'               target-rank used for the final iteration.
#'    }
#'
#'    \item{err}{  vector \cr
#'               Frobenious error archieved by each iteration.
#'    }
#'
#'    \item{.................}{.}
#'
#'
#'
#' @note  ...
#'
#' @author N. Benjamin Erichson, \email{nbe@st-andrews.ac.uk}
#'
#'
#' @examples
#'
#' library(rsvd)
#' data(highway)
#'
#' # Foreground/Background separation
#' out <- rrpca(highway, k=1, p=0, q=0, maxiter=20, svdalg='rsvd')
#'
#' # Display results for the 200th video frame
#' par(mfrow=c(1,3))
#' image(matrix(highway[,200], ncol=144, nrow=176), col = gray((0:255)/255))
#' image(matrix(out$L[,200], ncol=144, nrow=176), col = gray((0:255)/255))
#' image(matrix(out$S[,200], ncol=144, nrow=176), col = gray((0:255)/255))
#'

#' @export
rrpca <- function(A, k=NULL, lamb=NULL, gamma=1.25, rho=1.5, maxiter=10, tol=1.0e-3, svdalg='auto', p=10, q=1, trace=FALSE, ...) UseMethod("rrpca")

#' @export
rrpca.default <- function(A, k=NULL, lamb=NULL, gamma=1.25, rho=1.5, maxiter=10, tol=1.0e-3, svdalg='auto', p=10, q=1, trace=FALSE, ...) {
  #*************************************************************************
  #***        Author: N. Benjamin Erichson <nbe@st-andrews.ac.uk>        ***
  #***                              <2016>                               ***
  #***                       License: BSD 3 clause                       ***
  #*************************************************************************
  m <- nrow(A)
  n <- ncol(A)

  rrpcaObj = list(L = matrix(0, nrow = m, ncol = n),
                  S = matrix(0, nrow = m, ncol = n),
                  k = k,
                  lamb = lamb,
                  gamma = gamma,
                  rho = rho,
                  err = NULL)


  #Set target rank
  if(is.null(rrpcaObj$k)) rrpcaObj$k <- 2
  if(rrpcaObj$k>n) rrpcaObj$k <- n
  if(rrpcaObj$k<1) stop("Target rank is not valid!")

  A <- as.matrix(stats::na.omit(A))

  # Set lambda, gamma, rho
  if(is.null(rrpcaObj$lamb)) rrpcaObj$lamb <- m^-0.5
  if(is.null(rrpcaObj$gamma)) rrpcaObj$gamma <- 1.25
  if(is.null(rrpcaObj$rho)) rrpcaObj$rho <- 1.5

  # Compute matrix norms
  spectralNorm <- switch(svdalg,
                    svd = norm(A, "2"),
                    rsvd = rsvd(A, k=1, p=5, q=0, nu=0, nv=0)$d,
                    auto = rsvd(A, k=1, p=5, q=0, nu=0, nv=0)$d,
                    stop("Selected SVD algorithm is not supported!")
  )

  infNorm <- norm( A , "I") / rrpcaObj$lamb
  dualNorm <- max( spectralNorm , infNorm)
  froNorm <- norm( A , "F")

  # Normalize A
  Y <- A / dualNorm


  # Computing further tuning parameter
  mu <- rrpcaObj$gamma / spectralNorm
  mubar <- mu * 1e7
  mu <- min( mu*rrpcaObj$rho , mubar )
  muinv <- mu**-1


  rrpcaObj$niter <- 0
  err <- 1
  while(err > tol && rrpcaObj$niter <= maxiter) {

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Update S using soft-threshold
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      epsi = rrpcaObj$lamb*muinv
      tempS = A - rrpcaObj$L + muinv*Y
      rrpcaObj$S = matrix(0, nrow = m, ncol = n)
      #rrpcaObj$S = ifelse( tempS > epsi, tempS - epsi, ifelse( tempS < (- epsi), tempS + epsi, 0) )

      idxL <-which(tempS < -epsi)
      idxH <-which(tempS > epsi)
      rrpcaObj$S[idxL] <- tempS[idxL]+epsi
      rrpcaObj$S[idxH] <- tempS[idxH]-epsi


      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Singular Value Decomposition
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(svdalg=='auto'){
        if(rrpcaObj$k < (n/1.5)) {svdalg='rsvd'} else svdalg='svd'
      }
      svd_out <- switch(svdalg,
                        svd = svd(A - rrpcaObj$S + muinv*Y, nu = rrpcaObj$k, nv = rrpcaObj$k),
                        rsvd = rsvd(A - rrpcaObj$S + muinv*Y, k=rrpcaObj$k, p=p, q=q, ...),
                        stop("Selected SVD algorithm is not supported!")
      )


      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Update L
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      rrpcaObj$L =  svd_out$u %*% diag(svd_out$d - muinv, nrow=rrpcaObj$k, ncol=rrpcaObj$k)  %*% t(svd_out$v)


      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Predict optimal rank and update
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      kopt = sum(svd_out$d > muinv)
      if(kopt <= rrpcaObj$k){
        rrpcaObj$k = min(kopt+1, n)
      } else {
        rrpcaObj$k = min(kopt + round(0.05*n), n)
      }


      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Compute error
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Z = A - rrpcaObj$L - rrpcaObj$S
      Y = Y + mu * Z

      err = norm( Z , 'F') / froNorm
      rrpcaObj$err <- c(rrpcaObj$err, err)

      if(trace==TRUE){
        print(paste0('Iteration: ', rrpcaObj$niter ))
        print(paste0('******************' ))
        print(paste0('Fro. error = ', rrpcaObj$err ))
        print(paste0('k = ', rrpcaObj$k ))
      }

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Update mu
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      mu = min(mu*rrpcaObj$rho, mubar);
      muinv = 1 / mu

      rrpcaObj$niter = rrpcaObj$niter + 1

  }# End while loop

  class(rrpcaObj) <- "rrpca"
  return( rrpcaObj )

}



