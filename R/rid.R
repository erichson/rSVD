#' @title  Randomized interpolative decomposition (ID).
#
#' @description Randomized interpolative matrix decomposition.
#
#' @details
#' Algorithm for computing the ID of a rectangular \eqn{(m, n)} matrix \eqn{A}, with target rank 
#' \eqn{k << min(m,n)}. The input matrix is factored as \eqn{A = C * Z}, 
#' using the column pivoted QR decomposition. The factor matrix \eqn{C} is formed as a subset of 
#' columns of \eqn{A}, also called the partial column skeleton. 
#' If \code{mode='row'}, then the input matrix is factored as \eqn{A = Z * R}, using the 
#' row pivoted QR decomposition. The factor matrix \eqn{R} is now formed as
#' a subset of rows of \eqn{A}, also called the partial row skeleton. 
#' The factor matrix \eqn{Z} contains a \eqn{(k, k)} identity matrix as a submatrix, 
#' and is well-conditioned. 
#' 
#' If \eqn{rand='TRUE'} a probabilistic strategy is used to compute the decomposition, otherwise a
#' deterministic algorithm is used. 
#'
#'
#' @param A   Array_like. \cr
#'            A numeric input matrix (or data frame), with dimensions \eqn{(m, n)}. \cr
#'            If the data contain \eqn{NA}s na.omit is applied.
#'
#' @param k   Int, optional. \cr
#'            Determines the number of rows/columns to be selected. 
#'            It is required that \eqn{k} is smaller or equal to \eqn{min(m,n)}.
#'
#' @param mode  String c('column', 'row'). \cr
#'              Determines whether a subset of columns or rows is selected.
#'  
#' @param p       Int, optional. \cr
#'                Oversampling parameter (default \eqn{p=10}).
#'
#' @param q       Int, optional. \cr
#'                Number of power iterations (default \eqn{q=0}).
#'
#' @param idx_only  Bool (\eqn{TRUE}, \eqn{FALSE}), optional. \cr
#'              If (\eqn{TRUE}), the index set \code{idx} is returned, but not the matrix \code{C} or \code{R}. 
#'              This is more memory efficient, when dealing with large-scale data. 
#'              
#' @param rand  Bool (\eqn{TRUE}, \eqn{FALSE}). \cr
#'              If (\eqn{TRUE}), a probabilistic strategy is used, otherwise a deterministic algorithm is used.                
#'
#' @param ................. .
#'
#' @return \code{id} returns a list with class \eqn{id} containing the following components:
#'    \item{C}{ Array_like. \cr
#'              Column subset \eqn{C = A[,idx]}, if \code{mode='column'}; array with dimensions \eqn{(m, k)}.
#'    }
#'
#'    \item{R}{ Array_like. \cr
#'               Row subset \eqn{R = A[idx, ]}, if \code{mode='row'}; array with dimensions \eqn{(k, n)}.
#'    }
#'    
#'    \item{Z}{ Array_like \cr
#'             Well conditioned matrix; Dependng on the slected mode, this is an 
#'             array with dimensions \eqn{(k,n)} or \eqn{(m,k)}.
#'    }
#'    
#'    \item{idx}{ Array_like \cr
#'                The index set of the \eqn{k} selcted columns or rows used to form \eqn{C} or \eqn{R}. 
#'    }   
#'   
#'    \item{pivot}{ Array_like \cr
#'                Information on the pivoting strategy used during the decomposition. 
#'    }  
#'      
#'    \item{scores}{ Array_like .\cr
#'                   The scores (importancies) of the columns or rows of the input matrix \eqn{A}.
#'    } 
#'                   
#'    \item{scores.idx}{ Array_like .\cr
#'                   The scores (importancies) of the \eqn{k} selected columns or rows in \eqn{C} or \eqn{R}.                    
#'    }
#'    \item{.................}{.}
#'
#'
#'
#' @author N. Benjamin Erichson, \email{erichson@uw.edu}
#'
#' @seealso \code{\link{rcur}},
#'
#'

#' @export
rid <- function(A, k = NULL, mode = 'column', p = 10, q = 0, idx_only = FALSE, rand = 'TRUE') UseMethod("rid")

#' @export
rid.default <- function(A, k = NULL, mode = 'column', p = 10, q = 0, idx_only = FALSE, rand = 'TRUE') {

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Checks
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (any(is.na(A))) {
    warning("Missing values are omitted: na.omit(A).")
    A <- stats::na.omit(A)
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
  # Compute interpolative decompositon
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Pivoted QR decomposition
  if(rand == 'TRUE') { 
    out_rqb <- rqb(A, k = k, p = p, q = q)
    out <- qr(out_rqb$B, LAPACK=TRUE) 
  } else { 
    out <- qr(A, LAPACK=TRUE) 
  }
    
    
  
  # Get index set
  idObj$idx <- out$pivot[1:k] # Get row set
  idObj$pivot <- out$pivot # Get row set
  
  ordered.pivtos <- order(idObj$pivot)
  
  # Get R
  R <- qr.R( out ) 
  
  # Compute scores
  idObj$scores <- abs(diag(R))
  idObj$scores.idx <- idObj$scores[1:k]
  idObj$scores <- idObj$scores[ordered.pivtos]
  
  
  # Compute Z
  if(k == n) {
    V = pinv(R[1:k, 1:k])
  }else{
    V = pinv(R[1:k, 1:k]) %*% R[1:k, (k+1):n] 
  }
  idObj$Z = cbind(diag(k), V) 
  idObj$Z = matrix(idObj$Z[, ordered.pivtos], nrow = k, ncol = n)  

  if(mode == 'row') {
    idObj$Z = H(idObj$Z)
    
  }  
  
  # Create column / row subset
  if(idx_only == FALSE) {
    if(mode == 'column') {
      idObj$C = matrix(A[, idObj$idx], nrow = m, ncol = k)
      colnames(idObj$C) <- colnames(A)[idObj$idx]
      rownames(idObj$C) <- rownames(A)
    }  
    if(mode == 'row')  {
      idObj$R = H(matrix(A[, idObj$idx], nrow = m, ncol = k))  
      rownames(idObj$R) <- colnames(A)[idObj$idx]
      colnames(idObj$R) <- rownames(A)
      
    }
  }

  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Return
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  class(idObj) <- "rid"
  return( idObj )
  
  
}  