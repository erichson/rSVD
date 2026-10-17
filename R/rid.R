#' @title  Randomized interpolative decomposition (ID)
#
#' @description Randomized interpolative decomposition.
#
#' @details
#' Algorithm for computing the ID of a rectangular \eqn{(m, n)} matrix \eqn{A}, with target rank 
#' \eqn{k << min(m,n)}. The input matrix is approximated as
#' 
#' \deqn{A \approx C Z}{A is approximately C \%*\% Z}
#' 
#' using the column pivoted QR decomposition. The factor matrix \eqn{C} is formed as a subset of 
#' columns of \eqn{A}, also called the partial column skeleton. 
#' If \code{mode='row'}, then the input matrix is approximated as
#' 
#' \deqn{A \approx Z R}{A is approximately Z \%*\% R}
#' 
#' using the row pivoted QR decomposition. The factor matrix \eqn{R} is now formed as
#' a subset of rows of \eqn{A}, also called the partial row skeleton. 
#' The factor matrix \eqn{Z} contains a \eqn{(k, k)} identity matrix as a submatrix.
#' Its conditioning and the approximation error depend on the input and the
#' selected subset; column-pivoted QR does not guarantee the optimal rank-\eqn{k}
#' approximation.
#' 
#' If \eqn{rand='TRUE'} a probabilistic strategy is used to compute the decomposition, otherwise a
#' deterministic algorithm is used.
#'
#' Sparse matrices from the \pkg{Matrix} package are supported. The randomized
#' algorithm uses sparse matrix products to form a dense sketch; the returned
#' factors are dense matrices. The deterministic algorithm converts the input
#' to a dense matrix before column-pivoted QR factorization.
#'
#'
#' @param A   array_like; \cr
#'            numeric \eqn{(m, n)} input matrix (or data frame). \cr
#'            Rows containing missing values are omitted with a warning.
#'
#' @param k   integer, optional; \cr
#'            number of rows/columns to be selected. 
#'            It is required that \eqn{k} is smaller or equal to \eqn{min(m,n)}.
#'
#' @param mode  string c('column', 'row'), optional; \cr
#'              columns or rows ID.
#'  
#' @param p    integer, optional; \cr
#'             oversampling parameter (default \eqn{p=10}).
#'
#' @param q    integer, optional. \cr
#'             number of additional power iterations (default \eqn{q=0}).
#'
#' @param idx_only  bool, optional; \cr
#'              if (\eqn{TRUE}), omit the selected matrix \code{C} or \code{R}
#'              from the result. The indices and interpolation matrix \code{Z}
#'              are still computed and returned.
#'              
#' @param rand  bool, optional; \cr
#'              if (\eqn{TRUE}), a probabilistic strategy is used, otherwise a deterministic algorithm is used.                
#'
#'
#'
#' @return \code{rid} returns a list containing the following components:
#' \describe{
#'    \item{C}{ array_like; \cr
#'              column subset \eqn{C = A[,idx]}, if \code{mode='column'}; array with dimensions \eqn{(m, k)}.
#'    }
#'
#'    \item{R}{ array_like; \cr
#'              row subset \eqn{R = A[idx, ]}, if \code{mode='row'}; array with dimensions \eqn{(k, n)}.
#'    }
#'    
#'    \item{Z}{ array_like; \cr
#'             interpolation coefficients; depending on the selected mode, this is an
#'             array with dimensions \eqn{(k,n)} or \eqn{(m,k)}.
#'    }
#'    
#'    \item{idx}{ array_like; \cr
#'               index set of the \eqn{k} selected columns or rows used to form \eqn{C} or \eqn{R}. 
#'               Row indices refer to the input after rows with missing values are omitted.
#'    }   
#'   
#'    \item{pivot}{ array_like; \cr
#'                information on the pivoting strategy used during the decomposition. 
#'    }  
#'      
#'    \item{scores}{ array_like; \cr
#'                  absolute diagonal entries of the pivoted QR factor, mapped to
#'                  the original columns or rows of \eqn{A}. Entries without a
#'                  corresponding QR diagonal are \code{NA}; this can occur for
#'                  a rectangular input or a reduced randomized sketch.
#'    } 
#'                   
#'    \item{scores.idx}{ array_like; \cr
#'                  scores of the \eqn{k} selected columns or rows in \eqn{C} or \eqn{R}.                    
#'    }
#'}
#'
#' @references
#' \itemize{
#'   \item  [1] N. Halko, P. Martinsson, and J. Tropp.
#'          "Finding structure with randomness: probabilistic
#'          algorithms for constructing approximate matrix
#'          decompositions" (2009).
#'          (available at arXiv \url{https://arxiv.org/abs/0909.4061}).
#'  \item [2] N. B. Erichson, S. Voronin, S. L. Brunton and J. N. Kutz. 2019.
#'        Randomized Matrix Decompositions Using R.
#'        Journal of Statistical Software, 89(11), 1-48.
#'         \doi{10.18637/jss.v089.i11}.
#'}
#'
#'
#' @author N. Benjamin Erichson, \email{erichson@uw.edu}
#'
#' @seealso \code{\link{rcur}},
#'
#'
#' @examples
#' \dontrun{
#' # Load test image
#' data("tiger")
#'
#' # Compute (column) randomized interpolative decomposition
#' # Note that the image needs to be transposed for correct plotting
#' out <- rid(t(tiger), k = 150)
#'
#' # Show selected columns 
#' tiger.partial <- matrix(0, 1200, 1600)
#' tiger.partial[,out$idx] <- t(tiger)[,out$idx]
#' image(t(tiger.partial), col = gray((0:255)/255), useRaster = TRUE)
#'
#' # Reconstruct image
#' tiger.re <- t(out$C %*% out$Z)
#'
#' # Compute relative error
#' print(norm(tiger-tiger.re, 'F') / norm(tiger, 'F')) 
#'
#' # Plot approximated image
#' image(tiger.re, col = gray((0:255)/255))
#' }

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
