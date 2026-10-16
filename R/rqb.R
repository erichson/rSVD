#' @title  Randomized QB Decomposition (rqb).
#
#' @description Compute the near-optimal QB decomposition of a rectangular matrix.
#
#' @details
#' The randomized QB decomposition factors a rectangular \eqn{(m,n)} matrix \eqn{A} as
#' \eqn{A = Q * B}. \eqn{Q} is an \eqn{(m,k)} matrix with orthogonal columns, and \eqn{B} a \eqn{(k,n)} matrix.
#' The target rank is assumed to be \eqn{k << min(m,n)}.   
#'
#' \eqn{p} is an oversampling parameter to improve the approximation.
#' A value between 5 and 10 is recommended, and \eqn{p=10} is set by default.
#'
#' The parameter \eqn{q} specifies the number of power (subspace) iterations
#' to reduce the approximation error. This is recommended
#' if the the singular values decay slowly. In practice 1 or 2 iterations
#' achieve good results, however, computing power iterations increases the
#' computational time. The number of power iterations is set to \eqn{q=2} by default.
#'
#'
#' @param A       array_like; \cr
#'                real/complex \eqn{(m, n)} input matrix (or data frame).
#'
#' @param k       integer, optional; \cr
#'                target rank of the low-rank decomposition. It should satisfy \eqn{k << min(m,n)}.
#'
#' @param p       integer, optional; \cr
#'                oversampling parameter (default \eqn{p=10}).
#'
#' @param q       integer, optional; \cr
#'                number of power iterations (default \eqn{q=2}).
#'
#' @param sdist   string \eqn{c( 'unif', 'normal', 'rademacher')}, optional; \cr
#'                specifies the sampling distribution: \cr
#'                \eqn{'unif'} :  Uniform `[-1,1]`. \cr
#'                \eqn{'normal}' (default) : Normal `~N(0,1)`. \cr
#'                \eqn{'rademacher'} : Rademacher random variates. \cr
#'
#' @param rand  bool, optional; \cr
#'              If (\eqn{TRUE}), a probabilistic strategy is used, otherwise a deterministic algorithm is used.
#'
#'
#'@return \code{rqb} returns a list containing the following components:
#' \describe{
#'\item{Q}{  array_like; \cr
#'           matrix with orthogonal columns; \eqn{(m, k)} dimensional array.
#'}
#'
#'\item{B}{  array_like; \cr
#'           smaller matrix; \eqn{(k, n)} dimensional array.
#'}
#'}
#'
#'
#' @references
#' \itemize{
#'  \item [1] N. B. Erichson, S. Voronin, S. L. Brunton and J. N. Kutz. 2019.
#'        Randomized Matrix Decompositions Using {R}. 
#'        Journal of Statistical Software, 89(11), 1-48.
#'         \doi{10.18637/jss.v089.i11}.
#' 
#'   \item  [2] N. Halko, P. Martinsson, and J. Tropp.
#'          "Finding structure with randomness: probabilistic
#'          algorithms for constructing approximate matrix
#'          decompositions" (2009).
#'          (available at arXiv \url{https://arxiv.org/abs/0909.4061}).
#'}
#'
#' @author N. Benjamin Erichson, \email{erichson@berkeley.edu}
#' 
#' @seealso \code{\link{svd}}
#' 
#'

#' @export
rqb <- function(A, k=NULL, p=10, q=2, sdist="normal", rand = TRUE) UseMethod("rqb")

#' @export
rqb.default <- function(A, k=NULL, p=10, q=2, sdist="normal", rand = TRUE) {

    A <- as.matrix(A)

    # Create rqb object
    rqbObj <- list()
    
    #Dim of input matrix
    m <- nrow(A)
    n <- ncol(A)

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

    if(rand == TRUE) {
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

        rqbObj$Q <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
    }else{
        rqbObj$Q <- qr.Q( qr(A, complete = FALSE) , complete = FALSE )
    }
      
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Project the data matrix a into a lower dimensional subspace
    #B := Q.T * A
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rqbObj$B <- crossprod_help(rqbObj$Q , A )


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Return
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    class(rqbObj) <- "rqb"
    return(rqbObj) 
    
} # End rqb

