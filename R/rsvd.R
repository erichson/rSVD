#' @title  Randomized Singular Value Decomposition (rsvd).
#
#' @description Compute the near-optimal low-rank singular value decomposition (SVD) of a rectangular matrix.
#
#' @details
#' The singular value decomposition (SVD) plays a central role in data analysis and scientific computing.
#' Randomized SVD (rSVD) is a fast algorithm to compute the approximate
#' low-rank SVD of a rectangular \eqn{(m,n)} matrix \eqn{A}
#' using a probablistic algorithm.
#' Given a target rank \eqn{k << min(m,n)}, \code{rsvd} factors the input matrix \eqn{A} as
#' \eqn{A = U * diag(d) * V'}. The right singluar vectors are the columns of the
#' real or complex unitary matrix \eqn{U} . The left singular vectors are the columns
#' of the real or complex unitary matrix \eqn{V}. The singular values \eqn{d} are
#' non-negative and real numbers.
#'
#' \eqn{p} is an oversampling parameter to improve the approximation.
#' A value between 5 and 10 is recommended, and \eqn{p=10} is set by default.
#'
#' The parameter \eqn{q} specifies the number of power (subspace) iterations
#' (subspace iterations) to reduce the approximation error. This is recommended
#' if the the singular values decay slowly. In practice 1 or 2 iterations
#' achieve good results, however, computing power iterations increases the
#' computational time. The number of power iterations is set to \eqn{q=2} by default.
#'
#' If \eqn{k > (min(n,m)/2)}, a deterministic partial or truncated \code{\link{svd}}
#' algorithm might be faster.
#'
#'
#' @param A       Array_like. \cr
#'                A real/complex input matrix (or data frame), with dimensions \eqn{(m, n)}.
#'
#' @param k       Int, optional. \cr
#'                Determines the target rank of the low-rank decomposition. It should satisfy \eqn{k << min(m,n)}.
#'
#' @param nu       Int, optional. \cr
#'                 The number of left singular vectors to be computed. This must be between \eqn{0}
#'                 and \eqn{k}.
#'
#' @param nv       Int, optional. \cr
#'                 The number of right singular vectors to be computed. This must be between \eqn{0}
#'                 and \eqn{k}.
#'
#' @param p       Int, optional. \cr
#'                Oversampling parameter for (default \eqn{p=10}).
#'
#' @param q       Int, optional. \cr
#'                Number of power iterations (default \eqn{q=2}).
#'
#' @param sdist   String \eqn{c( 'unif', 'normal', 'rademacher')}, optional. \cr
#'                Specifies the sampling distribution. \cr
#'                \eqn{'unif'} :  Uniform `[-1,1]`. \cr
#'                \eqn{'normal}' (default) : Normal `~N(0,1)`. \cr
#'                \eqn{'rademacher'} : Rademacher random variates. \cr
#'
#' @param ............. .
#'
#'
#'@return \code{rsvd} returns a list containing the following three components:
#'\item{d}{  Array_like. \cr
#'           Singular values; 1-d array of length \eqn{(k)}.
#'}
#'
#'\item{u}{  Array_like. \cr
#'           Left singular values; array with dimensions \eqn{(m, k)} or \eqn{(m, nu)}.
#'}
#'
#'\item{v}{  Array_like. \cr
#'           Right singular values; array with dimensions \eqn{(n, k)} or \eqn{(n, nv)}. \cr
#'}
#'\item{.............}{.}
#'
#' @note The singular vectors are not unique and only defined up to sign
#' (a constant of modulus one in the complex case). If a left singular vector
#' has its sign changed, changing the sign of the corresponding right vector
#' gives an equivalent decomposition.
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
#' @author N. Benjamin Erichson, \email{erichson@uw.edu}
#' @seealso \code{\link{svd}}, \code{\link{rpca}}
#' @examples
#'
#'# Create a n by n Hilbert matrix of order n,
#'# with entries H[i,j] = 1 / (i + j + 1).
#'hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
#'H <- hilbert(n=50)
#'
#'# Low-rank (k=10) matrix approximation using rsvd
#'k=10
#'s <- rsvd(H, k=k)
#'Hre <- s$u %*% diag(s$d) %*% t(s$v) # matrix approximation
#'print(100 * norm( H - Hre, 'F') / norm( H,'F')) # percentage error

#'# Compare to truncated base svd
#'s <- svd(H)
#'Hre <- s$u[,1:k] %*% diag(s$d[1:k]) %*% t(s$v[,1:k]) # matrix approximation
#'print(100 * norm( H - Hre, 'F') / norm( H,'F')) # percentage error
#'

#' @export
rsvd <- function(A, k=NULL, nu=NULL, nv=NULL, p=10, q=2, sdist="normal") UseMethod("rsvd")

#' @export
rsvd.default <- function(A, k=NULL, nu=NULL, nv=NULL, p=10, q=2, sdist="normal") {
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
    if(m < n){
      A <- H(A)
      m <- nrow(A)
      n <- ncol(A)
      flipped <- TRUE
    } else flipped <- FALSE

    #Set target rank
    if(is.null(k)) k = n
    if(k > n) k <- n
    if(is.character(k)) stop("Target rank is not valid!")
    if(k < 1) stop("Target rank is not valid!")

    #Set oversampling parameter
    l <- round(k) + round(p)
    if(l > n) l <- n
    if(l < 1) stop("Target rank is not valid!")

    #Check if array is real or complex
    if(is.complex(A)) {
      isreal <- FALSE
    } else {
      isreal <- TRUE
    }

    #Set number of singular vectors
    if(is.null(nu)) nu <- k
    if(is.null(nv)) nv <- k
    if(nu < 0) nu <- 0
    if(nv < 0) nv <- 0
    if(nu > k) nu <- k
    if(nv > k) nv <- k
    if(flipped==TRUE) {
      temp <- nu
      nu <- nv
      nv <- temp
    }


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Generate a random sampling matrix O
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    O <- switch(sdist,
                normal = matrix(stats::rnorm(l*n), n, l),
                unif = matrix(stats::runif(l*n), n, l),
                rademacher = matrix(sample(c(-1,1), (l*n), replace = TRUE, prob = c(0.5,0.5)), n, l),
                stop("Selected sampling distribution is not supported!"))

    if(isreal==FALSE) {
      O <- O + switch(sdist,
                      normal = 1i * matrix(stats::rnorm(l*n), n, l),
                      unif = 1i * matrix(stats::runif(l*n), n, l),
                      rademacher = 1i * matrix(sample(c(-1,1), (l*n), replace = TRUE, prob = c(0.5,0.5)), n, l),
                      stop("Selected sampling distribution is not supported!"))
    }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Build sample matrix Y : Y = A * O
    #Note: Y should approximate the range of A
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Y <- A %*% O
    remove(O)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Orthogonalize Y using economic QR decomposition: Y=QR
    #If q > 0 perfrom q subspace iterations
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if( q > 0 ) {
        for( i in 1:q) {
            Y <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
            Z <- crossprod_help(A , Y )
            Z <- qr.Q( qr(Z, complete = FALSE) , complete = FALSE )
            Y <- A %*% Z
        }#End for
        remove(Z)
    }#End if

    Q <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
    remove(Y)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Project the data matrix a into a lower dimensional subspace
    #B := Q.T * A
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    B <- crossprod_help(Q , A )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Singular Value Decomposition
    #Note: B =: U * S * Vt
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rsvdObj <- svd(B, nu=nu, nv=nv) # Compute SVD
    rsvdObj$d <- rsvdObj$d[1:k] # Truncate singular values
    
    if(nu != 0) rsvdObj$u <- Q %*% rsvdObj$u # Recover left singular vectors

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Flipp SVD back
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(flipped == TRUE) {
        u_temp <- rsvdObj$u
        rsvdObj$u <- rsvdObj$v
        rsvdObj$v <- u_temp
    }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Return
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(nu == 0){ rsvdObj$u <- NULL}
    if(nv == 0){ rsvdObj$v <- NULL}
    class(rsvdObj) <- "rsvd"
    return(rsvdObj) 
    
} # End rsvd

