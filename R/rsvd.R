#' @title Randomized singular value decomposition
#'
#' @description
#' Approximate the leading singular values and singular vectors of a matrix.
#' Random matrix products identify a subspace in which the singular value
#' decomposition (SVD) is computed. This can reduce the factorization work when
#' the requested rank is small relative to the matrix dimensions.
#'
#' @details
#' A low-rank approximation represents the columns of a matrix using a small
#' set of shared directions. When these directions capture the structure needed
#' for an analysis, storing and working with the factors can save both memory
#' and computation. For an \eqn{m \times n}{m by n} input matrix \eqn{A}, the
#' rank-\eqn{k} reconstruction has the form
#' \deqn{A \approx U_k D_k V_k^{H}.}{A is approximately U \%*\% diag(d) \%*\% Conj(t(V)).}
#' The columns of \eqn{U_k} and \eqn{V_k} are orthonormal and the diagonal of
#' \eqn{D_k} contains the returned singular values. Each pair of left and right
#' vectors contributes one matrix component, weighted by its singular value.
#' The superscript \eqn{H} denotes the conjugate transpose, which is the ordinary
#' transpose for real data. Reconstruction requires both sets of vectors and
#' uses \code{Conj(t(v))} for complex input.
#'
#' To find these directions, the algorithm samples a subspace of width
#' \eqn{\ell = \min(k+p,m,n)}{ell = min(k + p, m, n)}. It refines the sample
#' with \code{q} power iterations and intermediate orthogonalization, then
#' computes an SVD of the matrix projected onto that subspace. The large input
#' is accessed through matrix products. When the subspace is small, the final
#' SVD acts on a smaller dense matrix.
#'
#' Oversampling gives the algorithm extra directions from which to identify
#' the dominant rank-\code{k} subspace; only \code{k} singular values are
#' returned. Increasing \code{p} or \code{q} can improve accuracy, at the cost
#' of more computation and, for \code{p}, more intermediate storage. The
#' function uses the supplied parameters throughout. It does not adapt them
#' to an error tolerance or estimate a numerical rank.
#'
#' The accuracy achievable at a given rank depends on the singular values of
#' the input. A full-rank matrix can have a useful low-rank approximation when
#' these values decay rapidly. With a flat spectrum and a small target rank,
#' even an optimal approximation leaves a large residual. If
#' \eqn{\sigma_j}{sigma[j]} are the singular values of the original matrix in
#' decreasing order, its optimal rank-\eqn{k} Frobenius error is
#' \deqn{\left(\sum_{j>k}\sigma_j^2\right)^{1/2}.}{The square root of the sum of sigma[j]^2 over j > k.}
#' This expression measures the contribution of the discarded components.
#' It is zero when \code{k} reaches the full available rank. Use this
#' truncation error as a reference for the randomized reconstruction residual,
#' with the same preprocessing and target rank for both approximations. An
#' arbitrary full-rank matrix cannot in general be recovered exactly at a
#' smaller rank.
#'
#' The decomposition acts directly on the supplied matrix, without centering
#' or scaling. See \code{\link{rpca}} for principal component analysis with
#' those preprocessing options.
#'
#' Real-valued sparse \pkg{Matrix} inputs are multiplied by thin dense factors
#' without converting the entire input to a dense matrix. Logical and pattern
#' matrices are treated as zero/one data. Complex input is supported for ordinary
#' R matrices. Dense storage is still needed for the random samples,
#' intermediate bases, projected matrices and returned singular vectors.
#' Explicitly reconstructing the approximation also creates a dense matrix.
#'
#' Omitting \code{k} requests all available singular values, even if \code{nu}
#' or \code{nv} requests only a few vectors. The resulting subspace no longer
#' provides the storage savings of a small-rank approximation. When the target
#' rank is a substantial fraction of the matrix dimensions, a deterministic
#' method may be faster. Use \code{\link{svd}} as a reference for a full SVD.
#'
#' @param A A real or complex numeric matrix or a numeric data frame. Sparse and
#' dense matrices from the \pkg{Matrix} package are supported; see Details.
#' The input must have positive dimensions and finite entries.
#'
#' @param k A positive integer giving the target rank. The default,
#' \code{NULL}, uses the smaller matrix dimension; larger requests are capped
#' at that value. Use a rank much smaller than both dimensions for the main
#' computational savings.
#'
#' @param nu An integer giving the number of left singular vectors to return.
#' The default, \code{NULL}, returns \code{k} vectors. Integer requests are
#' clamped to the range from zero to the effective target rank; zero omits
#' \code{u}.
#'
#' @param nv An integer giving the number of right singular vectors to return.
#' The default, \code{NULL}, returns \code{k} vectors. Integer requests are
#' clamped to the range from zero to the effective target rank; zero omits
#' \code{v}.
#'
#' @param p A non-negative integer oversampling parameter. The default is 10.
#' The sketch uses up to \code{k + p} columns, bounded by the matrix dimensions.
#'
#' @param q A non-negative integer giving the number of power (subspace)
#' iterations. The default is 2; zero skips the power iterations.
#'
#' @param sdist A character string selecting the sampling distribution.
#' \code{"normal"} (the default) uses standard normal entries,
#' \code{"unif"} uses uniform entries on \eqn{[0,1]},
#' and \code{"rademacher"} uses equally likely values of -1 and 1. For complex
#' input, the real and imaginary parts are sampled independently.
#'
#' @return
#' A list of class \code{"rsvd"} with the following components. Here
#' \code{k} denotes the target rank after clamping to the matrix dimensions.
#' \describe{
#'   \item{d}{A numeric vector of length \code{k} containing the approximate
#'   singular values in decreasing order. The values are non-negative.}
#'   \item{u}{An ordinary R matrix with \code{nrow(A)} rows and \code{nu}
#'   orthonormal columns, after applying the default and clamping rules above.
#'   Omitted when no left singular vectors are requested.}
#'   \item{v}{An ordinary R matrix with \code{ncol(A)} rows and \code{nv}
#'   orthonormal columns, after applying the default and clamping rules above.
#'   Omitted when no right singular vectors are requested.}
#' }
#'
#' @note
#' Singular vectors are defined only up to sign, or a constant of modulus one
#' for complex data; multiplying a corresponding left and right vector by the
#' same factor leaves the reconstruction unchanged. Repeated singular values
#' also allow rotations within their singular subspace, and closely spaced
#' values can make individual directions sensitive. Compare reconstruction
#' residuals or subspaces, since different valid singular vectors can represent
#' the same approximation.
#'
#' Use \code{set.seed} for reproducibility within a fixed numerical environment;
#' identical vectors are not guaranteed across BLAS or package versions.
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
#' @seealso \code{\link{svd}}, \code{\link{rqb}}, \code{\link{rpca}}
#'
#' @examples
#' # A low-rank signal with small dense noise.
#' set.seed(42)
#' A <- matrix(rnorm(60 * 5), 60, 5) %*% matrix(rnorm(5 * 40), 5, 40)
#' A <- A + 0.01 * matrix(rnorm(60 * 40), 60, 40)
#' k <- 5L
#' fit <- rsvd(A, k = k, p = 5, q = 2)
#' A_approx <- sweep(fit$u, 2, fit$d, "*") %*% Conj(t(fit$v))
#'
#' # Compare residuals at the same target rank on this small matrix.
#' d <- svd(A, nu = 0, nv = 0)$d
#' best <- sqrt(sum(d[-seq_len(k)]^2))
#' c(randomized = norm(A - A_approx, "F"), optimal = best)
#'
#' # Complex data require a conjugate transpose for reconstruction.
#' Z <- A + 1i * A[, rev(seq_len(ncol(A)))]
#' fit_complex <- rsvd(Z, k = k, p = 5, q = 2)
#' Z_approx <- sweep(fit_complex$u, 2, fit_complex$d, "*") %*%
#'   Conj(t(fit_complex$v))
#' sqrt(sum(Mod(Z - Z_approx)^2) / sum(Mod(Z)^2))
#'
#' # Sparse input stays sparse; the returned factors are dense.
#' S <- Matrix::Diagonal(x = c(9, 4, 1, rep(0, 17)))
#' fit_sparse <- rsvd(S, k = 3, p = 2, q = 1)
#' fit_sparse$d

#' @export
rsvd <- function(A, k=NULL, nu=NULL, nv=NULL, p=10, q=2, sdist="normal") UseMethod("rsvd")

#' @export
rsvd.default <- function(A, k=NULL, nu=NULL, nv=NULL, p=10, q=2, sdist="normal") {
    rsvd_compute(A, k=k, nu=nu, nv=nv, p=p, q=q, sdist=sdist)
}

# Shared range finder. Internal products let PCA apply centering and scaling
# implicitly; ordinary rsvd calls retain the original matrix operations.
rsvd_compute <- function(A, k=NULL, nu=NULL, nv=NULL, p=10, q=2,
                         sdist="normal", operators=NULL) {
    #*************************************************************************
    #***        Author: N. Benjamin Erichson <nbe@st-andrews.ac.uk>        ***
    #***                              <2015>                               ***
    #***                       License: BSD 3 clause                       ***
    #*************************************************************************

    # Preserve custom matrix representations; numeric data frames need coercion.
    if(is.data.frame(A)) A <- as.matrix(A)

    #Dim of input matrix
    m <- nrow(A)
    n <- ncol(A)

    #Flip matrix, if wide
    if(m < n){
      if(is.null(operators)) {
        A <- H(A)
      } else {
        original_operators <- operators
        operators <- list(
          forward = original_operators$adjoint,
          adjoint = original_operators$forward,
          project = function(Q) H(original_operators$forward(Q)))
      }
      temp <- m
      m <- n
      n <- temp
      flipped <- TRUE
    } else flipped <- FALSE

    if(is.null(operators)) {
      operators <- list(forward = function(X) A %*% X,
                        adjoint = function(X) crossprod_help(A, X),
                        project = function(Q) crossprod_help(Q, A))
    }

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
    Y <- operators$forward(O)
    remove(O)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Orthogonalize Y using economic QR decomposition: Y=QR
    #If q > 0 perform q subspace iterations
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if( q > 0 ) {
        for( i in 1:q) {
            Y <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
            Z <- operators$adjoint(Y)
            Z <- qr.Q( qr(Z, complete = FALSE) , complete = FALSE )
            Y <- operators$forward(Z)
        }#End for
        remove(Z)
    }#End if

    Q <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
    remove(Y)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Project the data matrix a into a lower dimensional subspace
    #B := Q.T * A
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    B <- operators$project(Q)

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

