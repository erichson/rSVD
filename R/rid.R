#' @title Randomized interpolative decomposition (ID)
#'
#' @description
#' Approximate a matrix using selected columns or rows of the original data and
#' an interpolation matrix. A randomized sketch can reduce the work needed to
#' select the subset when the target rank is small.
#'
#' @usage
#' rid(A, k = NULL, mode = "column", p = 10, q = 0,
#'     idx_only = FALSE, rand = TRUE)
#'
#' @param A A real or complex numeric matrix, a numeric data frame, or a matrix
#' from the \pkg{Matrix} package. Rows containing missing values are omitted
#' with a warning; the remaining entries must be finite.
#' @param k A positive integer giving the number of columns or rows to select.
#' The default, \code{NULL}, uses \code{min(dim(A))} after missing rows are
#' omitted. Larger requests are capped at this value. Use a substantially
#' smaller rank to obtain a compressed approximation.
#' @param mode Either \code{"column"} (the default) or \code{"row"}, specifying
#' which type of subset to select.
#' @param p A non-negative integer oversampling parameter, used only when
#' \code{rand = TRUE}. The default is 10. The sketch width is bounded by the
#' matrix dimensions.
#' @param q A non-negative integer number of power iterations, used only when
#' \code{rand = TRUE}. The default is zero.
#' @param idx_only Logical. If \code{TRUE}, leave both selected-data components
#' \code{C} and \code{R} as \code{NULL}. The interpolation matrix \code{Z},
#' indices and scores are still computed and returned. The default is
#' \code{FALSE}.
#' @param rand Logical. If \code{TRUE} (the default), select the subset using a
#' randomized QB sketch. If \code{FALSE}, use column-pivoted QR on the full
#' input, or its conjugate transpose for a row ID.
#'
#' @details
#' For an \eqn{m \times n}{m by n} matrix, a column ID has the form
#' \deqn{A \approx C Z.}{A is approximately C \%*\% Z.}
#' Here \code{C = A[, idx, drop = FALSE]} contains \code{k} actual columns,
#' and \code{Z} expresses all columns in terms of that subset. Its selected
#' columns satisfy \code{Z[, idx, drop = FALSE] = diag(k)}. For a row ID,
#' \deqn{A \approx Z R,}{A is approximately Z \%*\% R,}
#' \code{R = A[idx, , drop = FALSE]} contains the selected rows, and
#' \code{Z[idx, , drop = FALSE] = diag(k)}. Dimensions and indices in these
#' expressions refer to the input after missing-row omission.
#'
#' Unlike the orthonormal directions returned by \code{\link{rsvd}}, the
#' selected columns or rows retain their original units and can correspond to
#' measured variables or observations. Selection is sensitive to preprocessing:
#' centering or rescaling the data can change which subset is chosen. This
#' function does not center or scale the input automatically.
#'
#' The randomized algorithm applies column-pivoted QR to the projected factor
#' from \code{\link{rqb}}. Oversampling and power iterations can improve the
#' subspace used for selection, at additional time and memory cost. They are
#' fixed controls, not an adaptive error tolerance. The interpolation
#' coefficients are computed using a numerical pseudoinverse of the selected
#' QR block, so requesting more columns than the numerical rank is allowed.
#' Neither column-pivoted QR nor this pseudoinverse guarantees an optimal
#' rank-\code{k} approximation or uniformly small coefficients. Inspect the
#' reconstruction error when approximation quality matters; QR scores are
#' selection diagnostics, not statistical importance scores or probabilities.
#'
#' With sparse \pkg{Matrix} input, the randomized algorithm forms its dense
#' sketch using sparse matrix products. Returned factors are ordinary dense R
#' matrices. With \code{rand = FALSE}, the input is converted to a dense matrix
#' before QR. Forming \code{C \%*\% Z} or \code{Z \%*\% R} also produces a
#' dense reconstruction. Use \code{set.seed()} before a randomized call for a
#' repeatable result in the same computational environment; pivot choices can
#' differ across numerical libraries when columns have similar scores.
#'
#' @return
#' A list of class \code{"rid"} with the following components. Here \code{m}
#' is the number of complete input rows and \code{k} is the effective rank
#' after clamping.
#' \describe{
#'   \item{C}{The selected \code{m} by \code{k} column matrix for a column ID,
#'   or \code{NULL} for a row ID or when \code{idx_only = TRUE}.}
#'   \item{R}{The selected \code{k} by \code{n} row matrix for a row ID,
#'   or \code{NULL} for a column ID or when \code{idx_only = TRUE}.}
#'   \item{Z}{The dense interpolation matrix: \code{k} by \code{n} for a
#'   column ID, or \code{m} by \code{k} for a row ID.}
#'   \item{idx}{An integer vector of length \code{k} containing the selected
#'   column or row indices, in selection order. Row indices refer to the
#'   matrix after rows with missing values have been omitted.}
#'   \item{pivot}{The full QR column-pivot permutation, or the corresponding
#'   row permutation for a row ID. Its first \code{k} entries are \code{idx}.}
#'   \item{scores}{Absolute diagonal entries of the QR factor, mapped back to
#'   input column or row order. Entries without a corresponding QR diagonal
#'   are \code{NA}, which can occur for rectangular inputs or reduced sketches.}
#'   \item{scores.idx}{The QR scores of the selected subset, in the same
#'   order as \code{idx}.}
#'   \item{mode, rand}{The selection mode and randomization setting.}
#' }
#'
#' @references
#' N. Halko, P. G. Martinsson, and J. A. Tropp (2009).
#' Finding structure with randomness: probabilistic algorithms for constructing
#' approximate matrix decompositions. \url{https://arxiv.org/abs/0909.4061}.
#'
#' N. B. Erichson, S. Voronin, S. L. Brunton, and J. N. Kutz (2019).
#' Randomized Matrix Decompositions Using R.
#' \emph{Journal of Statistical Software}, 89(11), 1--48.
#' \doi{10.18637/jss.v089.i11}.
#'
#' @author N. Benjamin Erichson
#' @seealso \code{\link{rcur}}, \code{\link{rqb}}, \code{\link{rsvd}}
#'
#' @examples
#' set.seed(42)
#' A <- matrix(rnorm(60 * 5), 60, 5) %*% matrix(rnorm(5 * 35), 5, 35)
#' A <- A + 0.01 * matrix(rnorm(length(A)), nrow(A), ncol(A))
#'
#' # Column ID: all columns are approximated using five original columns.
#' fit <- rid(A, k = 5)
#' A_id <- fit$C %*% fit$Z
#' norm(A - A_id, "F") / norm(A, "F")
#' fit$idx
#' fit$Z[, fit$idx, drop = FALSE]
#'
#' # Row ID has the opposite factor order.
#' row_fit <- rid(A, k = 5, mode = "row")
#' norm(A - row_fit$Z %*% row_fit$R, "F") / norm(A, "F")
#'
#' # Retain indices and coefficients without returning the selected data.
#' small <- rid(A, k = 5, idx_only = TRUE, rand = FALSE)
#' is.null(small$C)
#' dim(small$Z)

#' @export
rid <- function(A, k = NULL, mode = 'column', p = 10, q = 0, idx_only = FALSE, rand = TRUE) UseMethod("rid")

#' @export
rid.default <- function(A, k = NULL, mode = 'column', p = 10, q = 0, idx_only = FALSE, rand = TRUE) {

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

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Init id object
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  idObj = list( C = NULL,
                R = NULL,
                Z = NULL,
                idx = NULL,
                pivot = NULL,
                scores = NULL,
                scores.idx = NULL,
                mode = mode,
                rand = rand)  
    
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Transpose input matrix if mode == 'row'
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  if(mode == 'row') A = H(A)
  
  m <- nrow(A)
  n <- ncol(A)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Set target rank
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.null(k)) k <- n
  if(k > min(m,n)) k <- min(m,n)
  if(k < 1) stop("Target rank is not valid!")


  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Compute interpolative decomposition
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Pivoted QR decomposition
  if(rand == TRUE) { 
    out_rqb <- rqb(A, k = k, p = p, q = q)
    out <- qr(as.matrix(out_rqb$B), LAPACK=TRUE)
  } else { 
    out <- qr(as.matrix(A), LAPACK=TRUE)
  }
    
    
  
  # Get index set
  idObj$idx <- out$pivot[1:k] # Get row set
  idObj$pivot <- out$pivot # Get row set
  
  ordered.pivot <- order(idObj$pivot)
  
  # Get R =: S
  S <- qr.R( out ) 
  
  # Compute scores
  idObj$scores <- abs(diag(S))
  idObj$scores.idx <- idObj$scores[1:k]
  idObj$scores <- idObj$scores[ordered.pivot]
  
  
  # Compute Z
  if(k == n) {
    idObj$Z <- diag(k)
  }else{
    V <- pinv(S[1:k, 1:k, drop=FALSE]) %*%
      S[1:k, (k+1):n, drop=FALSE]
    idObj$Z <- cbind(diag(k), V)
  }
  idObj$Z = matrix(idObj$Z[, ordered.pivot, drop=FALSE], nrow = k, ncol = n)

  if(mode == 'row') {
    idObj$Z = H(idObj$Z)
    
  }  
  
  # Create column / row subset
  if(idx_only == FALSE) {
    if(mode == 'column') {
      idObj$C = as.matrix(A[, idObj$idx, drop=FALSE])
    }  
    if(mode == 'row')  {
      idObj$R = H(as.matrix(A[, idObj$idx, drop=FALSE]))
      
    }
  }

  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Return
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  class(idObj) <- "rid"
  return( idObj )
}  
