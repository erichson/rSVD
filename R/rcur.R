#' @title Randomized CUR matrix decomposition
#'
#' @description
#' Approximate a matrix using selected columns and rows of the original data,
#' joined by a small matrix of coefficients. The selected data retain their
#' correspondence to measured variables and observations. Randomized column
#' selection is useful when the input has low-rank structure and the requested
#' rank is small.
#'
#' @param A A real or complex numeric matrix, a numeric data frame, or a matrix
#' from the \pkg{Matrix} package. Rows containing missing values are omitted
#' with a warning; the remaining entries must be finite.
#' @param k A positive integer giving the number of columns and rows to select.
#' The default, \code{NULL}, uses \code{min(dim(A))} after missing rows are
#' omitted. Larger requests are capped at this value.
#' @param p A non-negative integer oversampling parameter for randomized
#' column selection. The default is 10. Ignored when \code{rand = FALSE}.
#' @param q A non-negative integer number of power iterations for randomized
#' column selection. The default is zero. Ignored when \code{rand = FALSE}.
#' @param idx_only Logical. If \code{TRUE}, leave \code{C} and \code{R} as
#' \code{NULL} in the returned object. Their indices and the connector
#' \code{U} are still computed and returned. The default is \code{FALSE}.
#' @param rand Logical. If \code{TRUE} (the default), use randomized column
#' selection; otherwise use a deterministic column-pivoted QR decomposition.
#' Row selection from the chosen column matrix is deterministic in both cases.
#'
#' @details
#' A column interpolative decomposition (ID) represents a matrix using selected
#' columns and coefficients for combining them. CUR also selects rows, so
#' both sides of the approximation contain subsets of the original data.
#' For an \eqn{m \times n}{m by n} matrix, these factors are combined as
#' \deqn{A \approx C U R.}{A is approximately C \%*\% U \%*\% R.}
#' The \code{m} by \code{k} matrix \code{C} contains actual columns of the
#' input, and the \code{k} by \code{n} matrix \code{R} contains actual rows.
#' These subsets can be examined as measured variables and observations in
#' their original units. The \code{k} by \code{k} matrix \code{U} contains
#' the coefficients that connect the two subsets. Its role differs from that
#' of the orthonormal left singular-vector matrix in a singular value
#' decomposition (SVD).
#'
#' The algorithm first computes a column interpolative decomposition with
#' \code{\link{rid}}. It then chooses rows using column-pivoted QR of the
#' conjugate transpose of \code{C}. To form the connector, it uses these rows
#' to approximate the interpolation coefficients from the column ID,
#' \deqn{U = Z R^{+},}{U = Z times the numerical Moore-Penrose inverse of R,}
#' where \code{Z} is the interpolation matrix from the column ID.
#' The numerical pseudoinverse excludes negligible singular values using a
#' threshold that depends on matrix dimensions and floating-point precision.
#' This allows the calculation to proceed when the selected subsets are rank
#' deficient, as can happen when \code{k} exceeds the input's numerical rank.
#'
#' Oversampling and power iterations can improve randomized column selection
#' at the cost of more computation. The supplied settings remain fixed;
#' the algorithm does not adjust them to satisfy an error tolerance. Row
#' selection from the chosen column matrix is deterministic in both modes.
#'
#' The approximation error depends on the input's spectrum and on which
#' subsets are selected, including their conditioning. Column and row
#' selection does not guarantee the smallest possible rank-\code{k} error,
#' so inspect the reconstruction residual when choosing a rank. The QR scores
#' describe the selection procedure and have no interpretation as statistical
#' importance scores.
#'
#' Centering or rescaling changes the matrix used for selection and can
#' therefore change which variables and observations are chosen. The function
#' uses the supplied data without automatic centering or scaling.
#'
#' Sparse \pkg{Matrix} input is supported. Randomized column selection uses
#' sparse matrix products, but the selected column and row matrices and the
#' connector are dense. With \code{rand = FALSE}, column selection also
#' requires a dense input matrix. The option \code{idx_only = TRUE} reduces
#' the returned object's size by omitting the selected data. Those data are
#' still used internally to compute \code{U}, so this option does not remove
#' their temporary storage. Explicit reconstruction creates a dense matrix.
#'
#' Set the random seed before a randomized call for a repeatable result in
#' the same computational environment.
#'
#' @return
#' A list of class \code{"rcur"} with the following components. Dimensions
#' refer to the input after missing-row omission and the effective \code{k}
#' after clamping.
#' \describe{
#'   \item{C}{The dense \code{m} by \code{k} matrix containing columns
#'   \code{C.idx} of \code{A}. It is \code{NULL} when \code{idx_only = TRUE}.}
#'   \item{R}{The dense \code{k} by \code{n} matrix containing rows
#'   \code{R.idx} of \code{A}. It is \code{NULL} when \code{idx_only = TRUE}.}
#'   \item{U}{The dense \code{k} by \code{k} connector matrix.}
#'   \item{C.idx}{An integer vector of \code{k} selected column indices,
#'   in selection order.}
#'   \item{R.idx}{An integer vector of \code{k} selected row indices,
#'   in selection order. These index the matrix after rows with missing
#'   values have been omitted.}
#'   \item{C.scores}{Column QR scores in input column order, as returned by
#'   \code{\link{rid}}. Entries without a QR diagonal are \code{NA}.}
#'   \item{R.scores}{The \code{k} selected-row QR scores, in the same order
#'   as \code{R.idx}. Scores for unselected input rows are not returned.}
#'   \item{rand}{The randomization setting.}
#' }
#'
#' @references
#' N. B. Erichson, S. Voronin, S. L. Brunton, and J. N. Kutz (2019).
#' Randomized Matrix Decompositions Using R.
#' \emph{Journal of Statistical Software}, 89(11), 1--48.
#' \doi{10.18637/jss.v089.i11}.
#'
#' N. Halko, P. G. Martinsson, and J. A. Tropp (2009).
#' Finding structure with randomness: probabilistic algorithms for constructing
#' approximate matrix decompositions. \url{https://arxiv.org/abs/0909.4061}.
#'
#' @author N. Benjamin Erichson
#' @seealso \code{\link{rid}}, \code{\link{rsvd}}
#'
#' @examples
#' set.seed(42)
#' A <- matrix(rnorm(60 * 5), 60, 5) %*% matrix(rnorm(5 * 35), 5, 35)
#' A <- A + 0.01 * matrix(rnorm(length(A)), nrow(A), ncol(A))
#'
#' fit <- rcur(A, k = 5)
#' A_cur <- fit$C %*% fit$U %*% fit$R
#' norm(A - A_cur, "F") / norm(A, "F")
#' fit$C.idx
#' fit$R.idx
#'
#' # idx_only still returns the connector, so reconstruction remains possible.
#' small <- rcur(A, k = 5, idx_only = TRUE, rand = FALSE)
#' A_small <- A[, small$C.idx, drop = FALSE] %*% small$U %*%
#'   A[small$R.idx, , drop = FALSE]
#' norm(A - A_small, "F") / norm(A, "F")

#' @export
rcur <- function(A, k = NULL, p = 10, q = 0, idx_only = FALSE, rand = TRUE) UseMethod("rcur")

#' @export
rcur.default <- function(A, k = NULL, p = 10, q = 0, idx_only = FALSE, rand = TRUE) {

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Checks
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.data.frame(A)) A <- as.matrix(A)
  if (any(is.na(A))) {
    warning("Missing values are omitted: na.omit(A).")
    if(inherits(A, "Matrix")) {
      A <- A[Matrix::rowSums(is.na(A)) == 0, , drop=FALSE]
    } else {
      A <- stats::na.omit(A)
    }
  }     

  m <- nrow(A)
  n <- ncol(A)  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Init id object
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rcurObj = list( C = NULL,
                R = NULL,
                U = NULL,
                C.idx = NULL,
                R.idx = NULL,
                C.scores = NULL,
                R.scores = NULL,
                rand = rand)  
  

  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Set target rank
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.null(k)) k <- min(n,m)
  if(k > min(m,n)) k <- min(m,n)
  if(k < 1) stop("Target rank is not valid!")
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Compute column ID
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  out_cid <- rid(A, k = k, p = p, q = q, mode = 'column', idx_only = TRUE, rand = rand)

  # Only the selected columns need to be dense for the subsequent QR.
  C <- as.matrix(A[, out_cid$idx, drop=FALSE])
  
  # Select column subset
  if(idx_only == FALSE) {
    rcurObj$C <- C
  } 
  
  rcurObj$C.idx <- out_cid$idx
  rcurObj$C.scores <- out_cid$scores

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Compute row ID of C
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #out_rid <- rid(A[, rcurObj$C.idx], k = k, p = p, q = q, mode = 'row', idx_only = TRUE, rand = FALSE)
  out <- qr(H(C), LAPACK=TRUE)
  
  # Get index set
  rcurObj$R.idx <- out$pivot[1:k] # Get row set

  # Get R =: S
  S <- qr.R( out ) 
  #Q <- qr.Q( out ) 
  
  # Select row subset
  R <- as.matrix(A[rcurObj$R.idx, , drop=FALSE])
  if(idx_only == FALSE) {
    rcurObj$R <- R
  } 
  
  rcurObj$R.scores <- abs(diag(S))[1:k]

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Compute U
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #rcurObj$U = Z %*% MASS::ginv(matrix(A[rcurObj$R.idx, ], nrow = k, ncol = n) )    
  rcurObj$U = out_cid$Z %*% pinv(R)

  #colinv = matrix(0, 1 , n)
  #colinv[out_cid$pivot] = 1:n
  
  #V = rbind(diag(k), H(out_cid$Z))[colinv,]
  
  # Rt*Ut = V via RRt*Ut = R*V
  #RV = rcurObj$R %*% V
  #rcurObj$U = H(pinv(rcurObj$R %*% H(rcurObj$R)) %*% RV)

  # Rt*Ut = V
  #rcurObj$U = H(pinv(H(rcurObj$R)) %*% V)

  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Return
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  class(rcurObj) <- "rcur"
  return( rcurObj )
}  



