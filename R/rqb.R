#' @title Randomized QB decomposition
#'
#' @description Compute an orthonormal basis and the projection of a matrix onto that
#' basis, using either randomized range finding or a deterministic QR
#' factorization.
#'
#' @details
#' The QB decomposition approximates a matrix by
#' \deqn{A \approx Q B, \qquad B = Q^{H} A.}{A is approximately Q \%*\% B, where B = Conj(t(Q)) \%*\% A.}
#' The superscript \eqn{H} denotes the conjugate transpose, which is the ordinary
#' transpose for real input. Thus \code{Q \%*\% B} is the orthogonal projection
#' of \code{A} onto the column space of \code{Q}. The input is not centered
#' or scaled.
#'
#' When \code{rand = TRUE}, the basis width is
#' \eqn{\ell = \min(k+p,m,n)}{ell = min(k + p, m, n)} for integer \code{k}
#' and \code{p}. The oversampled columns are retained. Consequently, the result
#' generally has \code{k + p} basis vectors, rather than the \code{k} singular
#' triplets returned by \code{\link{rsvd}}. An SVD of the smaller matrix \code{B}
#' can be used to truncate the approximation further.
#'
#' Randomized range finding samples \code{A}, applies \code{q} power iterations
#' with intermediate orthogonalization, and projects onto the resulting basis.
#' Oversampling and power iterations can improve approximation quality,
#' especially when singular values decay slowly, but increase computational
#' work. The parameters are fixed: there is no adaptive error tolerance or
#' automatic rank selection. Low-rank approximation is useful when the retained
#' subspace captures the information needed by the application; full algebraic
#' rank by itself does not rule this out.
#'
#' When \code{rand = FALSE}, the basis width is
#' \eqn{\ell = \min(k,m,n)}{ell = min(k, m, n)}. The result uses the first
#' \eqn{\ell}{ell} columns of \code{qr.Q(qr(as.matrix(A)))}, with base R's
#' default QR conventions. For real input with a numerically independent leading
#' block, only that block needs dense factorization. Complex input, dependent
#' leading columns and other cases use the full dense QR fallback. This preserves
#' base QR's ordering, including its pivoting behavior; it does not guarantee an
#' optimal rank-\code{k} subspace.
#'
#' For randomized decompositions, real-valued sparse \pkg{Matrix} inputs are
#' multiplied by thin dense factors without converting the entire input to a
#' dense matrix. Logical and pattern matrices are treated as zero/one data.
#' Complex input is supported for ordinary R matrices. Samples, orthonormal
#' bases and projected matrices are dense. The deterministic fallback can
#' densify the entire input, so it is not a general sparse QR implementation.
#' Constructing \code{Q \%*\% B} also forms a dense matrix.
#'
#' Setting \code{k = NULL} uses the full available basis width and gives up the
#' small-rank storage advantage. The function does not reduce that width based
#' on an estimated numerical rank.
#'
#' @param A A real or complex numeric matrix or a numeric data frame. Sparse and
#' dense matrices from the \pkg{Matrix} package are supported; see Details.
#' Other matrix-like inputs are converted with \code{as.matrix}. The input
#' must have positive dimensions and finite entries.
#'
#' @param k A positive integer giving the target rank before randomized
#' oversampling. Requests larger than \code{min(dim(A))} are capped at that
#' value. The default, \code{NULL}, requests the full available rank. Use a
#' rank much smaller than the matrix dimensions for the main computational
#' savings.
#'
#' @param p A non-negative integer oversampling parameter for the randomized
#' algorithm. The default is 10. Unlike \code{rsvd}, \code{rqb} retains the
#' oversampled basis; see Details. This does not change the basis width when
#' \code{rand = FALSE}.
#'
#' @param q A non-negative integer giving the number of randomized power
#' (subspace) iterations. The default is 2; zero skips these iterations.
#' Unused when \code{rand = FALSE}.
#'
#' @param sdist A character string selecting the sampling distribution.
#' \code{"normal"} (the default) uses
#' standard normal entries, \code{"unif"} uses uniform entries on \eqn{[0,1]},
#' and \code{"rademacher"} uses equally likely values of -1 and 1. For complex
#' input, the real and imaginary parts are sampled independently. Unused when
#' \code{rand = FALSE}.
#'
#' @param rand A logical value. If \code{TRUE} (the default), use randomized range
#' finding; otherwise use the deterministic QR basis described below.
#'
#' @return
#' A list of class \code{"rqb"} containing two ordinary R matrices:
#' \describe{
#'   \item{Q}{An \eqn{m \times \ell}{m by ell} matrix with orthonormal columns.}
#'   \item{B}{An \eqn{\ell \times n}{ell by n} matrix containing the projection
#'   \eqn{Q^{H} A}{Conj(t(Q)) \%*\% A}.}
#' }
#' The input dimensions are \eqn{m \times n}{m by n}. The basis width \eqn{\ell}{ell}
#' is defined in Details and is not automatically reduced to the numerical rank.
#'
#' @note
#' The basis is not unique: signs, complex phases or rotations within the
#' same subspace can change \code{Q} and \code{B} while preserving their product.
#' Assess the reconstruction or subspace instead of individual basis entries.
#' Use \code{set.seed} for reproducible random sampling within a fixed numerical
#' environment. The deterministic branch does not draw random numbers.
#'
#' @references
#' \itemize{
#'   \item N. B. Erichson, S. Voronin, S. L. Brunton and J. N. Kutz. 2019.
#'   Randomized Matrix Decompositions Using R.
#'   Journal of Statistical Software, 89(11), 1-48.
#'   \doi{10.18637/jss.v089.i11}.
#'
#'   \item N. Halko, P. Martinsson and J. Tropp. 2009.
#'   Finding structure with randomness: probabilistic algorithms for
#'   constructing approximate matrix decompositions.
#'   \url{https://arxiv.org/abs/0909.4061}.
#' }
#'
#' @author N. Benjamin Erichson
#'
#' @seealso \code{\link{rsvd}}, \code{\link{svd}}, \code{\link{qr}}
#'
#' @examples
#' # An exact rank-three signal.
#' set.seed(42)
#' A <- matrix(rnorm(60 * 3), 60, 3) %*% matrix(rnorm(3 * 40), 3, 40)
#' qb <- rqb(A, k = 3, p = 2, q = 1)
#' dim(qb$Q)  # 60 by 5: randomized QB retains the two oversampling columns.
#' norm(A - qb$Q %*% qb$B, "F") / norm(A, "F")
#'
#' # Further truncation is performed in the projected matrix.
#' small <- svd(qb$B, nu = 3, nv = 3)
#' U <- qb$Q %*% small$u
#' A_approx <- sweep(U, 2, small$d[1:3], "*") %*% Conj(t(small$v))
#' norm(A - A_approx, "F") / norm(A, "F")
#'
#' # The deterministic branch returns k columns, without oversampling.
#' qb_deterministic <- rqb(A, k = 3, rand = FALSE)
#' dim(qb_deterministic$Q)
#' norm(A - qb_deterministic$Q %*% qb_deterministic$B, "F") / norm(A, "F")
#'
#' # Sparse randomized QB still returns ordinary dense factors.
#' S <- Matrix::Diagonal(x = c(9, 4, 1, rep(0, 17)))
#' qb_sparse <- rqb(S, k = 3, p = 2, q = 1)
#' dim(qb_sparse$B)

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

