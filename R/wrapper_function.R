#Helper function for conjugate transpose
H <- function( X ) {
  if(is.complex(X)) {
    return( Conj(t(X)) )
  } else {
    return( t(X) )
  }
}

#Helper function for conjugate crossprod
crossprod_help <- function( A , B ) {
  if(is.complex(A)) {
    return( crossprod( Conj(A) , B) )
  } else {
    return( crossprod( A , B ) )
  }
}

#Helper function for conjugate tcrossprod
tcrossprod_help <- function( A , B ) {
  if(is.complex(B)) {
    return( tcrossprod( A , Conj(B) ) )
  } else {
    return( tcrossprod( A , B ) )
  }
}
