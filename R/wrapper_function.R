#Helper function for conjugate transpose
#' @importFrom Matrix t
H <- function( X ) {
  if(is.complex(X)) {
    return( Conj(t(X)) )
  } else {
    return( t(X) )
  }
}

#Helper function for conjugate crossprod
#' @importFrom Matrix crossprod
crossprod_help <- function( A , B ) {
  if(is.complex(A)) {
    return( crossprod( Conj(A) , B) )
  } else {
    return( crossprod( A , B ) )
  }
}

#Helper function for conjugate tcrossprod
#' @importFrom Matrix tcrossprod
tcrossprod_help <- function( A , B ) {
  if(is.complex(B)) {
    return( tcrossprod( A , Conj(B) ) )
  } else {
    return( tcrossprod( A , B ) )
  }
}

#Helper function for Moore Penrose pseudoinverse
pinv <- function(A){
  s <- svd(A)
  nz <- s$d > s$d[1] * .Machine$double.eps
  if(any(nz)){
    return(s$v[, nz] %*% (H(s$u[, nz]) / s$d[nz]))
  } else {
    return(A)
  }
}
