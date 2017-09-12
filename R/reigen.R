#' @title  Randomized Spectral Decomposition of a matrix (reigen).
#
#' @description Computes the approximate low-rank eigendecomposition of a symmetric matrix.
#
#' @details
#' The eigenvalue decomposition plays a central role in data analysis and scientific computing.
#' Randomized eigen (reigen) is a fast algorithm to compute the the approximate
#' low-rank eigenvalue decomposition of \eqn{A'A} given the rectangular
#' \eqn{(m,n)} matrix \eqn{A} using a probablistic algorithm.
#' Given a target rank \eqn{k << n}, \code{reigen} factors the input matrix \eqn{A} as
#' \eqn{A'A = V * diag(d) * V'}. The eigenvectors are the columns of the
#' real or complex unitary matrix \eqn{V}. The eigenvalues \eqn{d} are
#' non-negative and real numbers.
#'
#' The parameter \eqn{p} is a oversampling parameter to improve the approximation.
#' A value between 2 and 10 is recommended and \eqn{p=10} is set as default.
#'
#' The parameter \eqn{q} specifies the number of normalized power iterations
#' (subspace iterations) to reduce the approximation error. This is recommended
#' if the the singular values decay slowly. In practice 1 or 2 iterations
#' archive good results, however, computing power iterations increases the
#' computational time. The number of power iterations is set to \eqn{q=1} by default.
#'
#' If \eqn{k > (min(n,m)/1.5)}, a deterministic partial or truncated \code{\link{eigen}}
#' algorithm might be faster.
#'
#'
#' @param A       array_like \cr
#'                a real/complex input matrix (or data frame), with dimensions \eqn{(m, n)}.
#'
#' @param k       int, optional \cr
#'                determines the target rank of the low-rank decomposition and should satisfy \eqn{k << min(m,n)}.
#'
#' @param p       int, optional \cr
#'                oversampling parameter for (default \eqn{p=10}).
#'
#' @param q       int, optional \cr
#'                number of power iterations (default \eqn{q=1}).
#'
#' @param sdist  str c('normal', 'unif', 'spixel'), optional \cr
#'               Specifies the sampling distribution. \cr
#'               'unif' : (default) Uniform `[-1,1]`. \cr
#'               'normal' : Normal `~N(0,1)`. \cr
#'               'col' : Random column sampling. \cr
#'
#' @param ............. .
#'
#'
#'@return \code{reigen} returns a list containing the following two components:
#'\item{values}{  array_like \cr
#'           Eigenvalues; 1-d array of length \eqn{(k)}.
#'}
#'
#'\item{vectors}{  array_like \cr
#'           Eigenvectors; array with dimensions \eqn{(k, k)}. \cr
#'}
#'\item{.............}{.}
#'
#' @note The eigenvectors are not unique and only defined up to sign
#' (a constant of modulus one in the complex case).
#'
#'
#' @references
#' \itemize{
#'   \item  [1] N. Halko, P. Martinsson, and J. Tropp.
#'          "Finding structure with randomness: probabilistic
#'          algorithms for constructing approximate matrix
#'          decompositions" (2009).
#'          (available at arXiv \url{http://arxiv.org/abs/0909.4061}).
#'   \item  [2] S. Voronin and P.Martinsson.
#'          "RSVDPACK: Subroutines for computing partial singular value
#'          decompositions via randomized sampling on single core, multi core,
#'          and GPU architectures" (2015).
#'          (available at `arXiv \url{http://arxiv.org/abs/1502.05366}).
#' }
#'
#' @author N. Benjamin Erichson, \email{nbe@st-andrews.ac.uk}
#' @seealso \code{\link{rsvd}}, \code{\link{rpca}}, \code{\link{eigen}}
#' @examples
#'library(rsvd)
#'set.seed(123)
#'
#'#Create real random test matrix with dimension (m, n) and rank k
#'m = 10
#'n = 5
#'k = 3
#'A <- matrix(runif(m*k), m, k)
#'A <- A %*% t(A)
#'A <- A[,1:n]
#'
#'AtA = t(A) %*% A
#'
#'# Randomized low-rank eigenvalue decomposition k=3
#'reigen.out <- reigen(A, k=3)
#'AtA.re = reigen.out$vectors %*% diag(reigen.out$values) %*% t(reigen.out$vectors)
#'100 * norm( AtA - AtA.re, 'F') / norm( AtA,'F') #Percentage reconstruction error << 1e-8
#'print(reigen.out$values) # print eigenvalues
#'
#'# Compare with the deterministic eigenvalue decomposition
#'eigen.out <- eigen(AtA)
#'AtA.re2 = eigen.out$vectors %*% diag(eigen.out$values) %*% t(eigen.out$vectors)
#'100 * norm( AtA - AtA.re2, 'F') / norm( AtA,'F') #Percentage reconstruction error << 1e-8
#'print(eigen.out$values) # print eigenvalues
#'all.equal(eigen.out$values[1:k], reigen.out$values)


#' @export
reigen <- function(A, k=NULL, p=10, q=1, sdist="unif") UseMethod("reigen")

#' @export
reigen.default <- function(A, k=NULL, p=10, q=1, sdist="unif") {
  #*************************************************************************
  #***        Author: N. Benjamin Erichson <nbe@st-andrews.ac.uk>        ***
  #***                              <2015>                               ***
  #***                       License: BSD 3 clause                       ***
  #*************************************************************************

  A <- as.matrix(A)

  #Dim of input matrix
  m <- nrow(A)
  n <- ncol(A)

  #Flipp matrix, if wide
  if(m<n){
    A <- H(A)
    m <- nrow(A)
    n <- ncol(A)
    flipped <- TRUE
  } else flipped <- FALSE

  #Set target rank
  if(is.null(k)) k=n
  if(k>n) k <- n
  if(is.character(k)) stop("Target rank is not valid!")
  if(k<1) stop("Target rank is not valid!")

  #Set oversampling parameter
  l <- round(k)+round(p)
  if(l>n) l <- n
  if(l<1) stop("Target rank is not valid!")

  #Check if array is real or complex
  if(is.complex(A)) {
    isreal <- FALSE
  } else {
    isreal <- TRUE
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Generate a random sampling matrix O
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  O <- switch(sdist,
              normal = matrix(stats::rnorm(l*n), n, l),
              unif = matrix(stats::runif(l*n), n, l),
              col = sample.int(n, size = l),
              stop("Selected sampling distribution is not supported!"))

  if(isreal==FALSE) {
    O <- O + switch(sdist,
                    normal = 1i * matrix(stats::rnorm(l*n), n, l),
                    unif = 1i * matrix(stats::runif(l*n), n, l),
                    col = NULL,
                    stop("Selected sampling distribution is not supported!"))
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Build sample matrix Y : Y = A * O
  #Note: Y should approximate the range of A
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(sdist=='col'){
    Y = A[,O]
  }else{
    Y <- A %*% O
  }
  remove(O)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Orthogonalize Y using economic QR decomposition: Y=QR
  #If q > 0 perfrom q subspace iterations
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if( q > 0 ) {
    for( i in 1:q) {
      if( ((2*i-2) %% 1) == 0 ) {
        Y <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
      }

      Z = crossprod_help( A , Y )

      if( ((2*i-1) %% 1) == 0 ) {
        Z <- qr.Q( qr(Z, complete = FALSE) , complete = FALSE )
      }

      Y <- A %*% Z
      remove(Z)
    }#End for
  }#End if

  Q <- qr.Q( qr(Y) , complete = FALSE )
  remove(Y)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Project the data matrix a into a lower dimensional subspace
  #B = Q.T * A
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  B <- crossprod_help( Q , A )


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Compute BB*
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BBt = B %*% H(B)
  BBt=0.5*(BBt+H(BBt)) # ensure symmetry


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Compute eigenvalue decomposition
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  out <- eigen(BBt, symmetric=TRUE, only.values = FALSE, EISPACK = FALSE)

  #Create rsvd class
  eigenObj <- list(values = NULL,
                  vectors = NULL)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Recover eigenvectors / eigenvalues
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(flipped==TRUE) {
      eigenObj$vectors <- (Q %*% out$vectors[,1:k])
  }else {
      eigenObj$vectors <- crossprod_help(B, t(t(out$vectors[,1:k]) * out$values[1:k]**-0.5 ) )
  }

  eigenObj$values <- out$values[1:k]

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Compute eigenvalue decomposition
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  return(eigenObj)

} # End rsvd
