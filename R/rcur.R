#devtools::use_package("MASS")
#library("MASS")

#' @title  Randomized CUR matrix decomposition.
#
#' @description Randomized CUR matrix decomposition.
#
#' @details
#' Algorithm for computing the CUR of a rectangular \eqn{(m, n)} matrix \eqn{A}, with target rank 
#' \eqn{k << min(m,n)}. The input matrix is factored as \eqn{A = C * U * R}, 
#' using the \code{\link{rid}} decomposition. The factor matrix \eqn{C} is formed as a subset of 
#' columns of \eqn{A}, also called the partial column skeleton. The factor matrix \eqn{R} is formed as
#' a subset of rows of \eqn{A}, also called the partial row skeleton. The factor matrix \eqn{U} 
#' is well-conditioned. 
#' 
#' If \eqn{rand='TRUE'} a probabilistic strategy is used to compute the decomposition, otherwise a
#' deterministic algorithm is used. 
#'
#'
#' @param A   Array_like. \cr
#'            A numeric input matrix (or data frame), with dimensions \eqn{(m, n)}. \cr
#'            If the data contain \eqn{NA}s na.omit is applied.
#' @param k   Int, optional. \cr
#'            Sets the target rank of the low-rank approximation, i.e., the number of columns/rows
#'             to be selected. It is required that \eqn{k} is smaller or equal to \eqn{min(m,n)}.
#'        
#' @param p       Int, optional. \cr
#'                Oversampling parameter (default \eqn{p=10}).
#'
#' @param q       Int, optional. \cr
#'                Number of power iterations (default \eqn{q=0}).
#'                                            
#' @param rand  Bool (\eqn{TRUE}, \eqn{FALSE}). \cr
#'              If (\eqn{TRUE}), a probabilistic strategy is used, otherwise a deterministic algorithm is used.               
#'
#' @param idx_only  Bool (\eqn{TRUE}, \eqn{FALSE}), optional. \cr
#'              If (\eqn{TRUE}), the index set \code{C.idx} and \code{R.idx} is returned, but not 
#'              the matrices \code{C} and \code{R}. 
#'              This is more memory efficient, when dealing with large-scale data. 
#'              
#' @param rand  Bool (\eqn{TRUE}, \eqn{FALSE}). \cr
#'              If (\eqn{TRUE}), a probabilistic strategy is used, otherwise a deterministic algorithm is used.                
#'
#' @param ................. .
#'
#' @return \code{id} returns a list with class \eqn{id} containing the following components:
#'    \item{C}{ Array_like. \cr
#'              Column subset \eqn{C = A[,C.idx]}; array with dimensions \eqn{(m, k)}.
#'    }
#'
#'    \item{R}{ Array_like. \cr
#'               Row subset \eqn{R = A[R.idx, ]}; array with dimensions \eqn{(k, n)}.
#'    }
#'    
#'    \item{U}{ Array_like \cr
#'             Well conditioned matrix; array with dimensions \eqn{(k,k)}.
#'    }
#'    
#'    \item{C.idx}{ Array_like \cr
#'                The index set of the \eqn{k} selcted columns used to form \eqn{C}. 
#'    }   
#' 
#'    \item{R.idx}{ Array_like \cr
#'                The index set of the \eqn{k} selcted rows used to form \eqn{R}. 
#'    }   
#'       
#'      
#'    \item{C.scores}{ Array_like .\cr
#'                   The scores (importancies) of the columns of the input matrix \eqn{A}.
#'    } 
#' 
#'    \item{R.scores}{ Array_like .\cr
#'                   The scores (importancies) of the rows of the input matrix \eqn{A}.
#'    }                   
#'                                                       
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
rcur <- function(A, k = NULL, p = 10, q = 0, idx_only = FALSE, rand = 'TRUE') UseMethod("rcur")

#' @export
rcur.default <- function(A, k = NULL, p = 10, q = 0, idx_only = FALSE, rand = 'TRUE') {

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Checks
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (any(is.na(A))) {
    warning("Missing values are omitted: na.omit(A).")
    A <- stats::na.omit(A)
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
  
  # Select column subset
  if(idx_only == FALSE) {
    rcurObj$C <- matrix(A[, out_cid$idx], nrow = m, ncol = k)    
    colnames(rcurObj$C) <- colnames(A)[out_cid$idx]
    rownames(rcurObj$C) <- rownames(A)
    
  } 
  
  rcurObj$C.idx <- out_cid$idx
  rcurObj$C.scores <- out_cid$scores
  Z <- out_cid$Z
  remove(out_cid)
    
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Compute row ID of C
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  out_rid <- rid(A[, rcurObj$C.idx], k = k, p = p, q = q, mode = 'row', idx_only = TRUE, rand = FALSE)

  # Select column subset
  if(idx_only == FALSE) {
    rcurObj$R <- matrix(A[out_rid$idx, ], nrow = k, ncol = n) 
    rownames(rcurObj$R) <- rownames(A)[out_rid$idx]
    colnames(rcurObj$R) <- colnames(A)
    
  } 
  
  rcurObj$R.idx <- out_rid$idx
  rcurObj$R.scores <- out_rid$scores
  remove(out_rid)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Compute U
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rcurObj$U = Z %*% MASS::ginv(matrix(A[rcurObj$R.idx, ], nrow = k, ncol = n) )    
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Return
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  class(rcurObj) <- "rcur"
  return( rcurObj )
  
  
}  