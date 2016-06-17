#' @title  Randomized Singular Value Decomposition (rSVD).
#
#' @description Compute the approximate low-rank singular value decomposition (SVD) of a rectangular matrix.
#
#' @details
#' The singular value decomposition (SVD) plays a central role in data analysis and scientific computing.
#' Randomized SVD (rSVD) is a fast algorithm to compute the the approximate
#' low-rank SVD of a rectangular \eqn{(m,n)} matrix \eqn{A}
#' using a probablistic algorithm.
#' Given a target rank \eqn{k << n}, \code{rsvd} factors the input matrix \eqn{A} as
#' \eqn{A = U * diag(d) * V'}. The right singluar vectors are the columns of the
#' real or complex unitary matrix \eqn{U} . The left singular vectors are the columns
#' of the real or complex unitary matrix \eqn{V}. The singular values \eqn{d} are
#' non-negative and real numbers.
#'
#' The parameter \eqn{p} is a oversampling parameter to improve the approximation.
#' A value between 2 and 10 is recommended and \eqn{p=5} is set as default.
#'
#' The parameter \eqn{q} specifies the number of normalized power iterations
#' (subspace iterations) to reduce the approximation error. This is recommended
#' if the the singular values decay slowly. In practice 1 or 2 iterations
#' archive good results, however, computing power iterations increases the
#' computational time. The number of power iterations is set to \eqn{q=2} by default.
#'
#' If \eqn{k > (min(n,m)/1.5)}, a deterministic partial or truncated \code{\link{svd}}
#' algorithm might be faster.
#'
#'
#' @param A       array_like \cr
#'                a real/complex input matrix (or data frame), with dimensions \eqn{(m, n)}.
#'
#' @param k       int, optional \cr
#'                determines the target rank of the low-rank decomposition and should satisfy \eqn{k << min(m,n)}.
#'
#' @param nu       int, optional \cr
#'                 the number of left singular vectors to be computed. This must be between \eqn{0}
#'                 and \eqn{k}.
#'
#' @param nv       int, optional \cr
#'                 the number of right singular vectors to be computed. This must be between \eqn{0}
#'                 and \eqn{k}.
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
#' @param vt      bool (\eqn{TRUE}, \eqn{FALSE}), optional \cr
#'                \eqn{TRUE} : returns the transposed right singular vectors \eqn{vt}. \cr
#'                \eqn{FALSE} : (default) returns the right singular vectors as \eqn{v}, this is the format
#'                as \code{\link{svd}} returns \eqn{v}.
#'
#' @param ............. .
#'
#'
#'@return \code{rsvd} returns a list containing the following three components:
#'\item{d}{  array_like \cr
#'           Singular values; 1-d array of length \eqn{(k)}.
#'}
#'
#'\item{u}{  array_like \cr
#'           Right singular values; array with dimensions \eqn{(m, k)}.
#'}
#'
#'\item{v (or vt)}{  array_like \cr
#'           Left singular values; array with dimensions \eqn{(n, k)}. \cr
#'           Or if \eqn{vt=TRUE}, array with dimensions \eqn{(k, n)}.
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
#' @author N. Benjamin Erichson, \email{nbe@st-andrews.ac.uk}
#' @seealso \code{\link{svd}}, \code{\link{rpca}}
#' @examples
#'library(rsvd)
#'data(tiger)
#'
#'# Randomized SVD, low-rank approximation k=100 for image compression
#'s <- rsvd(tiger, k=100)
#'tiger.re = s$u %*% diag(s$d) %*% t(s$v) # reconstruct image
#'print(100 * norm( tiger - tiger.re, 'F') / norm( tiger,'F')) # percentage error
#'
#'# Display orginal and reconstrucuted image
#'par(mfrow=c(1,2))
#'image(tiger, col = gray((0:255)/255))
#'image(tiger.re, col = gray((0:255)/255))
#'


#' @export
rsvd <- function(A, k=NULL, nu=NULL, nv=NULL, p=10, q=1, sdist="unif", vt=FALSE) UseMethod("rsvd")

#' @export
rsvd.default <- function(A, k=NULL, nu=NULL, nv=NULL, p=10, q=1, sdist="unif", vt=FALSE) {
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

    #Set number of singular vectors
    if(is.null(nu)) nu <- k
    if(is.null(nv)) nv <- k
    if(nu<0) nu <- 0
    if(nv<0) nv <- 0
    if(nu>k) nu <- k
    if(nv>k) nv <- k
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
    #Singular Value Decomposition
    #Note: B = U" * S * Vt
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Compute SVD
    rsvdObj <- La.svd(B, nu=nu, nv=nv)

    #Recover right singular vectors
    if(nu!=0) {
        rsvdObj$u <- Q %*% rsvdObj$u
    }else{ rsvdObj$u <- matrix(0)}

    if(nv==0) {
        rsvdObj$v <- matrix(0)
    }

    #Return
    if(flipped==TRUE) {
        u_temp <- rsvdObj$u
        rsvdObj$u <- H(rsvdObj$v)
        rsvdObj$d <- rsvdObj$d[1:k]
        if(vt==FALSE) {
          rsvdObj$v <- u_temp
          rsvdObj$vt <- NULL
        } else { rsvdObj$vt <- H(u_temp) }
        if(nu==0){ rsvdObj$v <- NULL}
        if(nv==0){ rsvdObj$u <- NULL}
        return(rsvdObj)

    } else {
        rsvdObj$d <- rsvdObj$d[1:k]
        if(vt==FALSE) {
          rsvdObj$v <- H(rsvdObj$v)
          rsvdObj$vt <- NULL
        }
        if(nu==0){ rsvdObj$u <- NULL}
        if(nv==0){ rsvdObj$v <- NULL}
        return(rsvdObj)
      }
} # End rsvd

