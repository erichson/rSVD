#' @title  Randomized Singular Value Decomposition (rsvd).
#
#' @description The randomized SVD computes the near-optimal low-rank approximation of a rectangular matrix
#' using a fast probablistic algorithm.
#
#' @details
#' The singular value decomposition (SVD) plays an important role in data analysis, and scientific computing.
#' Given a rectangular \eqn{(m,n)} matrix \eqn{A}, and a target rank \eqn{k << min(m,n)}, 
#' the SVD factors the input matrix \eqn{A} as
#' 
#' \deqn{ A  =  U_{k} diag(d_{k}) V_{k}^\top }{ A = U diag(d) t(V)}
#' 
#' The \eqn{k} left singular vectors are the columns of the
#' real or complex unitary matrix \eqn{U}. The \eqn{k} right singular vectors are the columns
#' of the real or complex unitary matrix \eqn{V}. The \eqn{k} dominant singular values are the 
#' entries of \eqn{d}, and non-negative and real numbers.
#'
#' \eqn{p} is an oversampling parameter to improve the approximation.
#' A value of at least 10 is recommended, and \eqn{p=10} is set by default.
#'
#' The parameter \eqn{q} specifies the number of power (subspace) iterations
#' to reduce the approximation error. The power scheme is recommended,
#' if the singular values decay slowly. In practice, 2 or 3 iterations
#' achieve good results, however, computing power iterations increases the
#' computational costs. The power scheme is set to \eqn{q=2} by default.
#'
#' If \eqn{k > (min(n,m)/4)}, a deterministic partial or truncated \code{\link{svd}}
#' algorithm might be faster.
#'
#'
#' @param A       array_like; \cr
#'                a real/complex \eqn{(m, n)} input matrix (or data frame) to be decomposed.
#'
#' @param k       integer; \cr
#'                specifies the target rank of the low-rank decomposition. \eqn{k} should satisfy \eqn{k << min(m,n)}.
#'
#' @param nu      integer, optional; \cr
#'                number of left singular vectors to be returned. \eqn{nu} must be between \eqn{0} and \eqn{k}.
#'
#' @param nv      integer, optional; \cr
#'                number of right singular vectors to be returned. \eqn{nv} must be between \eqn{0} and \eqn{k}.
#'
#' @param p       integer, optional; \cr
#'                oversampling parameter (by default \eqn{p=10}).
#'
#' @param q       integer, optional; \cr
#'                number of additional power iterations (by default \eqn{q=2}).
#'
#' @param sdist   string \eqn{c( 'unif', 'normal', 'rademacher')}, optional; \cr
#'                specifies the sampling distribution of the random test matrix: \cr
#'                		\eqn{'unif'} :  Uniform `[-1,1]`. \cr
#'                		\eqn{'normal}' (default) : Normal `~N(0,1)`. \cr
#'                		\eqn{'rademacher'} : Rademacher random variates. \cr
#'
#'@return \code{rsvd} returns a list containing the following three components:
#'\describe{
#'\item{d}{  array_like; \cr
#'           singular values; vector of length \eqn{(k)}.
#'}
#'
#'\item{u}{  array_like; \cr
#'           left singular vectors; \eqn{(m, k)} or \eqn{(m, nu)} dimensional array.
#'}
#'
#'\item{v}{  array_like; \cr
#'           right singular vectors; \eqn{(n, k)} or \eqn{(n, nv)} dimensional array. \cr
#'}
#'}
#'
#' @note The singular vectors are not unique and only defined up to sign
#' (a constant of modulus one in the complex case). If a left singular vector
#' has its sign changed, changing the sign of the corresponding right vector
#' gives an equivalent decomposition.
#'
#'
#' @references
#' \itemize{
#'  \item [1] N. B. Erichson, S. Voronin, S. L. Brunton and J. N. Kutz. 2019.
#'        Randomized Matrix Decompositions Using {R}. 
#'        Journal of Statistical Software, 89(11), 1-48.
#'         \doi{10.18637/jss.v089.i11}.
#' 
#'  \item  [2] N. Halko, P. Martinsson, and J. Tropp.
#'          "Finding structure with randomness: probabilistic
#'          algorithms for constructing approximate matrix
#'          decompositions" (2009).
#'          (available at arXiv \url{https://arxiv.org/abs/0909.4061}).
#' }
#'
#' @author N. Benjamin Erichson, \email{erichson@berkeley.edu}
#' 
#' @seealso \code{\link{svd}}, \code{\link{rpca}}
#' 
#' @examples
#'library('rsvd')
#'
#'# Create a n x n Hilbert matrix of order n,
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
    # Flip SVD back
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(nu == 0){ rsvdObj$u <- NULL}
    if(nv == 0){ rsvdObj$v <- NULL}
    if(flipped == TRUE) {
        u_temp <- rsvdObj$u
        rsvdObj$u <- rsvdObj$v
        rsvdObj$v <- u_temp
    }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Return
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    class(rsvdObj) <- "rsvd"
    return(rsvdObj) 
    
} # End rsvd

