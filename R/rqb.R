#' @title  Randomized QB Decomposition (rqb).
#
#' @description Compute the near-optimal QB decomposition of a rectangular matrix.
#
#' @details
#' The randomized QB decomposition factors a rectangular \eqn{(m,n)} matrix \eqn{A} as
#' \eqn{A \approx Q * B}. In the randomized case, \eqn{Q} has \eqn{l} orthonormal
#' columns and \eqn{B} has dimensions \eqn{(l,n)}, where
#' \eqn{l = \min(k+p,m,n)} for integer \eqn{k} and \eqn{p}. The oversampled
#' basis is retained. In the deterministic case, \eqn{l = \min(k,m,n)}.
#' The target rank is assumed to be \eqn{k << min(m,n)}.   
#'
#' \eqn{p} is an oversampling parameter to improve the approximation.
#' A value between 5 and 10 is recommended, and \eqn{p=10} is set by default.
#'
#' The parameter \eqn{q} specifies the number of power (subspace) iterations
#' to reduce the approximation error. This is recommended
#' if the singular values decay slowly. In practice 1 or 2 iterations
#' achieve good results, however, computing power iterations increases the
#' computational time. The number of power iterations is set to \eqn{q=2} by default.
#'
#' In the randomized case, sparse \pkg{Matrix} inputs are kept sparse during
#' multiplication. Only the samples, basis and projected matrix are dense.
#' The deterministic case returns the first \eqn{k} basis vectors of base R's
#' QR factorization. For real input with a numerically independent leading
#' block, only those \eqn{k} columns need a dense factorization. Other cases,
#' including complex input and deficient leading columns, use the full dense
#' QR fallback. The deterministic basis follows base QR's column ordering and
#' is not necessarily an optimal rank-\eqn{k} subspace. Both cases return
#' ordinary R matrices for \code{Q} and \code{B}.
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
#'                \eqn{'unif'} :  Uniform `[0,1]`. \cr
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
#'           matrix with orthonormal columns; \eqn{(m, l)} dimensional array.
#'}
#'
#'\item{B}{  array_like; \cr
#'           smaller matrix; \eqn{(l, n)} dimensional array.
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

    # Keep Matrix inputs for products with thin dense factors.
    if(!inherits(A, "Matrix")) A <- as.matrix(A)

    # Create rqb object
    rqbObj <- list()
    
    #Dim of input matrix
    m <- nrow(A)
    n <- ncol(A)

    #Set target rank
    if(is.null(k)) k = min(m, n)
    if(k > min(m, n)) k <- min(m, n)
    if(is.character(k)) stop("Target rank is not valid!")
    if(k < 1) stop("Target rank is not valid!")

    #Set oversampling parameter
    l <- round(k) + round(p)
    if(l > min(m, n)) l <- min(m, n)
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
        Y <- as.matrix(A %*% O)
        remove(O)
    
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Orthogonalize Y using economic QR decomposition: Y=QR
        #If q > 0 perform q subspace iterations
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if( q > 0 ) {
            for( i in 1:q) {
                Y <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
                Z <- as.matrix(crossprod_help(A , Y ))
                Z <- qr.Q( qr(Z, complete = FALSE) , complete = FALSE )
                Y <- as.matrix(A %*% Z)
            }#End for
            remove(Z)
        }#End if

        rqbObj$Q <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
    }else{
        # With real input, base QR keeps the column order unless a column is
        # numerically dependent. If the first k columns have full rank, their
        # Householder vectors are also the first k vectors of the full QR.
        # This avoids factoring unused columns while retaining the same basis.
        qrA <- NULL
        # Preserve full QR's rejection of non-finite values in later columns.
        # Summary methods check sparse inputs without a dense logical mask.
        if(isreal && k < min(m, n) && k == as.integer(k) &&
           is.finite(min(A)) && is.finite(max(A))) {
            qrA <- qr(as.matrix(A[, seq_len(k), drop = FALSE]))
            if(qrA$rank < k) qrA <- NULL
        }
        # Complex QR pivots across every column. Dependent real prefixes may
        # also need later columns, so both cases retain the full factorization.
        if(is.null(qrA)) {
            A <- as.matrix(A)
            qrA <- qr(A)
        }
        # Apply Q to only k identity columns instead of forming the full Q.
        identity <- diag(if(isreal) 1 else 1+0i, nrow = m, ncol = k)
        rqbObj$Q <- qr.qy(qrA, identity)
        remove(qrA, identity)
    }
      
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Project the data matrix a into a lower dimensional subspace
    #B := Q.T * A
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rqbObj$B <- as.matrix(crossprod_help(rqbObj$Q , A ))


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Return
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    class(rqbObj) <- "rqb"
    return(rqbObj) 
    
} # End rqb

