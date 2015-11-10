#' @title  Randomized Singular Value Decomposition .
#
#' @description Compute the approximate low-rank singular value decomposition (SVD) of a rectangular matrix.
#
#' @details
#' The singular value decomposition plays a central role in data analysis and scientific computing.
#' Randomized SVD (rSVD) is a very fast algorithm to compute the the approximate
#' low-rank singular value decompositon (SVD) of a rectangular \eqn{(m,n)} matrix \eqn{A}
#' using a probablistic algorithm.
#' Given a target rank \eqn{k << n}, \code{rsvd} factors the input matrix \eqn{A} as
#' \eqn{A = U * diag(d) * V'}. The right singluar vectors are the columns of the
#' real or complex unitary matrix \eqn{U} . The left singular vectors are the columns
#' of the real or complex unitary matrix \eqn{V}. The singular values \eqn{d} are
#' non-negative and real numbers.
#'
#' The paramter \eqn{p} is a oversampling parameter to improve the approximation.
#' A value between 2 and 10 is recommended and \eqn{p=5} is set as default.
#'
#' The paramter \eqn{q} specifies the number of normlized power iterations
#' (subspace iterations) to reduce the approximation error. This is recommended
#' if the the singular values decay slowly and in practice 1 or 2 iterations
#' achive good results. However, computing power iterations is increasing the
#' computational time. The number of power iterations is set to \eqn{q=2} by default.
#'
#' If \eqn{k > (min(n,m)/1.5)}, a deterministic partial or trancated \code{\link{svd}}
#' algorithm might be faster.
#'
#'
#' @param A       array_like \cr
#'                a real/complex input matrix (or data frame), with dimensions \eqn{(m, n)}.
#'
#' @param k       int, optional \cr
#'                determines the target rank of the low-rank decomposition and should satisfy \eqn{k << min(m,n)}.
#'
#'
#' @param p       int, optional \cr
#'                oversampling parameter (default \eqn{p=5}).
#'
#' @param q       int, optional \cr
#'                number of power iterations (default \eqn{q=5}).
#'
#' @param method  str c('standard', 'fast'), optional \cr
#'                'standard' : (default): Standard algorithm as described in [1, 2]. \cr
#'                'fast' : Version II algorithm as described in [2].
#'
#' @param sdist  str c('normal', 'unif'), optional \cr
#'               Specifies the distribution to draw the random samples from.
#'               'unif' : (default) Uniform `[-1,1]`. \cr
#'               'norm' : Normal `~N(0,1)`. \cr
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
#'\item{u}{  array_like \cr
#'           Right singular values; array with dimensions \eqn{(m, k)}.
#'}
#'\item{d}{  array_like \cr
#'           Singular values; 1-d array of length \eqn{(k)}.
#'}
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
#' @seealso \code{\link{svd}}
#' @examples
#' Will be added.

rsvd <- function(A, k=NULL, nu=NULL, nv=NULL, p=5, q=2, method='standard', sdist="unif", vt=FALSE) UseMethod("rsvd")

rsvd.default <- function(A, k=NULL, nu=NULL, nv=NULL, p=5, q=2, method='standard', sdist="unif", vt=FALSE) {
    #*************************************************************************
    #***        Author: N. Benjamin Erichson <nbe@st-andrews.ac.uk>        ***
    #***                              <2015>                               ***
    #***                       License: BSD 3 clause                       ***
    #*************************************************************************

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
      a <- nu
      nu <- nv
      nv <- a
    }


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Generate a random sampling matrix O
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    O <- switch(sdist,
                normal = matrix(rnorm(l*n), n, l),
                unif = matrix(runif(l*n), n, l),
                stop("Selected sampling distribution is not supported!"))

    if(isreal==FALSE) {
      O <- O + switch(sdist,
                normal = 1i * matrix(rnorm(l*n), n, l),
                unif = 1i * matrix(runif(l*n), n, l),
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
    #Note: check_finite=False may give a performance gain
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if( q > 0 ) {
        for( i in 1:q) {
          if( ((2*i-2) %% q) == 0 ) {
            Y <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
          }

          Z = crossprod_help( A , Y )

          if( ((2*i-1) %% q) == 0 ) {
            Z <- qr.Q( qr(Z, complete = FALSE) , complete = FALSE )
          }

          Y <- A %*% Z
          remove(Z)
        }#End for
    }#End if

    Q <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
    remove(Y)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Project the data matrix a into a lower dimensional subspace
    #B = Q.T * A
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    B <- crossprod_help( Q , A )

    if(method=='standard') {
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Singular Value Decomposition
      #Note: B = U" * S * Vt
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Compute SVD
      rsvdObj <- La.svd(B, nu=nu, nv=nv)

      #Recover right singular vectors
      if(nu==0) {
        rsvdObj$u <- matrix(0)
      } else{
        rsvdObj$u <- Q %*% rsvdObj$u
      }

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
        } else {
          rsvdObj$vt <- H(u_temp)
        }
        return(rsvdObj)


      } else {
        rsvdObj$d <- rsvdObj$d[1:k]
        if(vt==FALSE) {
          rsvdObj$v <- H(rsvdObj$v)
          rsvdObj$vt <- NULL
        }
        return(rsvdObj)
      }

    }else if(method=='fast'){
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Orthogonalize B.T using reduced QR decomposition: B.T = Q" * R"
      #Note: reduced QR returns Q and R, and destroys B_gpu
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Compute QR
      qr_out <- qr( H(B) , complete = FALSE)
      Qstar <- qr.Q(qr_out , complete = FALSE)
      Rstar <- qr.R(qr_out , complete = FALSE)

      #Compute singular value decomposition
      rsvdObj <- La.svd(Rstar, nu=nv, nv=nu)

      #Recover right and left singular vectors
      if(nu==0) {
        U <- matrix(0)
      } else{
        U <- tcrossprod_help( Q, rsvdObj$v )
      }

      if(nv==0) {
        V <- matrix(0)
      } else{
        V <- Qstar %*% rsvdObj$u
      }



      #Return
      if(flipped==TRUE) {
        rsvdObj$u <- V
        rsvdObj$d <- rsvdObj$d[1:k]
        if(vt==FALSE) {
          rsvdObj$v <- U
          rsvdObj$vt <- NULL
        } else {
          rsvdObj$vt <- H(U)
        }
        return(rsvdObj)

        } else {
        rsvdObj$u <- U
        rsvdObj$d <- rsvdObj$d[1:k]
        if(vt==FALSE) {
          rsvdObj$v <- V
          rsvdObj$vt <- NULL
        } else {
          rsvdObj$vt <- H(V)
        }
        return(rsvdObj)
      }

    }

}#End rsvd
