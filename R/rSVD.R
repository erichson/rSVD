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

rsvd <- function(A, k=NULL, p=5, q=2, method='standard', sdist="unif", vt=FALSE) {
    # Randomized Singular Value Decomposition.
    #
    # Randomized algorithm for computing the approximate low-rank singular value
    # decomposition of a rectangular (m, n) matrix `a` with target rank `k << n`.
    # The input matrix a is factored as `a = U * diag(s) * Vt`. The right singluar
    # vectors are the columns of the real or complex unitary matrix `U`. The left
    # singular vectors are the columns of the real or complex unitary matrix `V`.
    # The singular values `s` are non-negative and real numbers.
    #
    # The paramter `p` is a oversampling parameter to improve the approximation.
    # A value between 2 and 10 is recommended.
    #
    # The paramter `q` specifies the number of normlized power iterations
    # (subspace iterations) to reduce the approximation error. This is recommended
    # if the the singular values decay slowly and in practice 1 or 2 iterations
    # achive good results. However, computing power iterations is increasing the
    # computational time.
    #
    # If k > (n/1.5), partial SVD or trancated SVD might be faster.
    #
    #
    # Arguments
    # ----------
    # A : array_like
    #   Real/complex input matrix  `a` with dimensions `(m, n)`.
    #
    # k : int
    #   `k` is the target rank of the low-rank decomposition, k << min(m,n).
    #
    # p : int, optional
    #   `p` sets the oversampling parameter (default `p=5`).
    #
    # q : int, optional
    #   `q` sets the number of power iterations (default `q=2`).
    #
    # method : str `{'standard', 'fast'}`, optional
    #   'standard' (default): Standard algorithm as described in [1, 2].
    #
    #   'fast' : Version II algorithm as described in [2].
    #
    # sdist : str `{'normal', 'unif'}`, optional
    #   'unif' (default): Uniform `[-1,1]`.
    #
    #   'norm' : Normal `~N(0,1)`.
    #
    # vt : bool `{TRUE, FALSE}`, optional
    #   'TRUE' : returns the transpose `vt` of the right singular vectors
    #
    #   'FALSE' (default): returns the right singular vectors `v`. This is the format
    #    as svd {base} returns `v`.
    #
    # Value
    # -------
    # The returned value is a list with components
    #
    # u:  array_like
    #   Right singular values, array of shape `(m, k)`.
    #
    # d : array_like
    #   Singular values, 1-d array of length `k`.
    #
    # v : array_like
    #   Left singular values, array of shape `(n, k)`. Or
    #   if vt=TRUE, array of shape `(k, n)` is returned.
    #
    # Notes
    # -----
    #
    #
    #
    # References
    # ----------
    # N. Halko, P. Martinsson, and J. Tropp.
    # "Finding structure with randomness: probabilistic
    # algorithms for constructing approximate matrix
    # decompositions" (2009).
    # (available at `arXiv <http://arxiv.org/abs/0909.4061>`_).
    #
    # S. Voronin and P.Martinsson.
    # "RSVDPACK: Subroutines for computing partial singular value
    # decompositions via randomized sampling on single core, multi core,
    # and GPU architectures" (2015).
    # (available at `arXiv <http://arxiv.org/abs/1502.05366>`_).
    #
    #
    # Examples
    # --------
    #
    #
    #*************************************************************************
    #***        Author: N. Benjamin Erichson <nbe@st-andrews.ac.uk>        ***
    #***                              <2015>                               ***
    #***                       License: BSD 3 clause                       ***
    #*************************************************************************

    m <- nrow(A)
    n <- ncol(A)
    if(is.null(k)) k=n
    l <- k+p
    if(l>n){
      l <- n
      k <- n
    }

    #Check if array is real or complex
    if(is.complex(A)) {
      isreal <- FALSE
    } else {
      isreal <- TRUE
    }

    #Flipp matrix, if wide
    if(m<n){
      A <- H(A)
      m <- nrow(A)
      n <- ncol(A)
      if(l>n){
        l <- n
        k <- n
      }
      flipped <- TRUE
    } else flipped <- FALSE

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Generate a random sampling matrix O
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(sdist=='normal') {
      O <- matrix(rnorm(l*n), n, l)
      if(isreal==FALSE) O <- O + 1i * matrix(rnorm(l*n), n, l)
    } else if(sdist=='unif') {
      O <- matrix(runif(l*n), n, l)
      if(isreal==FALSE) O <- O + 1i * matrix(runif(l*n), n, l)
    }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Build sample matrix Y : Y = A * O
    #Note: Y should approximate the range of A
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Y <- A %*% O
    remove(O)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Orthogonalize Y using economic QR decomposition: Y=QR
    #If q > 0 perfrom q subspace iterations
    #Note: check_finite=False may give a performance gain
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if( q > 0 ) {
        for( i in 1:q) {
          if( ((2*i-2) %% q) == 0 ) {
            Y <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
          }

          Z = crossprod_help( A , Y )

          if( ((2*i-1) %% q) == 0 ) {
            Z <- qr.Q( qr(Z, complete = FALSE) , complete = FALSE )
          }

          Y <- A %*% Z
          remove(Z)
        }#End for
    }#End if

    Q <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
    remove(Y)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Project the data matrix a into a lower dimensional subspace
    #B = Q.T * A
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    B <- crossprod_help( Q , A )

    if(method=='standard') {
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Singular Value Decomposition
      #Note: B = U" * S * Vt
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Compute SVD
      svd_out <- La.svd(B, nu=k, nv=k)

      #Recover right singular vectors
      svd_out$u <- Q %*% svd_out$u

      #Return
      if(flipped==TRUE) {
        u_temp <- svd_out$u
        svd_out$u <- H(svd_out$v)
        svd_out$d <- svd_out$d[1:k]
        if(vt==FALSE) {
          svd_out$v <- u_temp
        } else {
          svd_out$v <- H(u_temp)
        }
        return(svd_out)

      } else {
        svd_out$d <- svd_out$d[1:k]
        if(vt==FALSE) {
          svd_out$v <- H(svd_out$v)
        }
        return(svd_out)
      }

    }else if(method=='fast'){
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Orthogonalize B.T using reduced QR decomposition: B.T = Q" * R"
      #Note: reduced QR returns Q and R, and destroys B_gpu
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Compute QR
      qr_out <- qr( H(B) , complete = FALSE)
      Qstar <- qr.Q(qr_out , complete = FALSE)
      Rstar <- qr.R(qr_out , complete = FALSE)

      #Compute right singular vectors
      svd_out <- La.svd(Rstar, nu=k, nv=k)

      U <- tcrossprod_help( Q, svd_out$v )
      V <- Qstar %*% svd_out$u

      #Return
      if(flipped==TRUE) {
        svd_out$u <- V
        svd_out$d <- svd_out$d[1:k]
        if(vt==FALSE) {
          svd_out$v <- U
        } else {
          svd_out$v <- H(U)
        }
        return(svd_out)

        } else {
        svd_out$u <- U
        svd_out$d <- svd_out$d[1:k]
        if(vt==FALSE) {
          svd_out$v <- V
        } else {
          svd_out$v <- H(V)
        }
        return(svd_out)
      }

    }

}#End rsvd
