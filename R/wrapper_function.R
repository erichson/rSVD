# Check controls before allocating a sketch or advancing the random stream.
# Ranks retain their historical rounding/truncation in the callers; iteration
# counts and sketch oversampling must be non-negative integers.
check_scalar_number <- function(x, name, lower = -Inf, integer = FALSE) {
  if(length(x) != 1L || !is.numeric(x) || is.complex(x) ||
     !is.finite(x) || x < lower || (integer && x != floor(x))) {
    requirement <- if(integer) "a finite integer" else "a finite number"
    if(is.finite(lower)) requirement <- paste0(requirement, " >= ", lower)
    stop("'", name, "' must be ", requirement, ".", call. = FALSE)
  }
  invisible(NULL)
}

check_flag <- function(x, name) {
  if(!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("'", name, "' must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(NULL)
}

check_matrix_input <- function(A, finite = TRUE) {
  dimensions <- dim(A)
  if(length(dimensions) != 2L || any(dimensions < 1L)) {
    stop("'A' must be a matrix with positive dimensions.", call. = FALSE)
  }
  if(is.matrix(A) && !(is.numeric(A) || is.complex(A) || is.logical(A))) {
    stop("'A' must contain numeric, complex, or logical values.", call. = FALSE)
  }
  # Summary methods inspect sparse and structured matrices without creating a
  # full logical mask or treating unstored entries as arbitrary slot contents.
  # Leave other matrix representations to their own product methods: forcing
  # an eager scan here can materialize a delayed matrix or change its dispatch.
  if(finite && (is.matrix(A) || inherits(A, "Matrix"))) {
    all_finite <- if(is.complex(A)) all(is.finite(A)) else
      is.finite(min(A)) && is.finite(max(A))
    if(!all_finite) {
      stop("'A' must contain finite values (no NA/NaN/Inf).", call. = FALSE)
    }
  }
  invisible(NULL)
}

check_sampling <- function(p, q, sdist) {
  check_scalar_number(p, "p", lower = 0, integer = TRUE)
  check_scalar_number(q, "q", lower = 0, integer = TRUE)
  if(!is.character(sdist) || length(sdist) != 1L || is.na(sdist) ||
     !sdist %in% c("normal", "unif", "rademacher")) {
    stop("'sdist' must be 'normal', 'unif', or 'rademacher'.", call. = FALSE)
  }
  invisible(NULL)
}

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
  # Roundoff in a rank-deficient rectangular matrix can exceed eps * d[1].
  # Use a dimension-aware relative cutoff, as in Matrix::rankMatrix, so those
  # directions do not produce a very large, inaccurate inverse. Keep the
  # threshold near machine precision to retain small resolved singular values.
  nz <- s$d > s$d[1] * (max(dim(A)) * .Machine$double.eps)
  if(any(nz)){
    return(s$v[, nz, drop = FALSE] %*%
             (H(s$u[, nz, drop = FALSE]) / s$d[nz]))
  } else {
    return(matrix(if(is.complex(A)) 0+0i else 0,
                  nrow = ncol(A), ncol = nrow(A)))
  }
}
