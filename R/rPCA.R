rpca <- function(A, k=NULL, center=TRUE, scale=TRUE, whiten=FALSE, p=5, q=2, ...) UseMethod("rpca")

rpca.default <- function(A, k=NULL, center=TRUE, scale=TRUE, whiten=FALSE, p=5, q=2, ...) {
    # Principal component analysis (PCA) using randomized SVD (rSVD)
    #
    # Linear dimensionality reduction using approximated Singular Value
    # Decomposition of the data and keeping only the most significant
    # singular vectors to project the data to a lower dimensional space.
    #
    # Arguments
    # ----------
    # k : int, optional
    #   `k` is the number of principle components to compute. The number of
    #    components should be k << min(m,n). If k > (n/1.5), partial SVD
    #    or trancated SVD might be faster. However, when k is not given
    #    k is set to min(m,n)
    #
    # p : int, optional
    #   `p` sets the oversampling parameter (default k=5).
    #
    # q : int, optional
    #   `q` sets the number of power iterations (default=2).
    #
    #
    # whiten : bool `{TRUE, FALSE}`, optional
    #   When True (False by default) the `components_` vectors are divided
    #   by the singular values to ensure uncorrelated outputs with unit
    #   component-wise variances.
    #   Whitening will remove some information from the transformed signal
    #   (the relative variance scales of the components) but can sometime
    #   improve the predictive accuracy of the downstream estimators by
    #   making their data respect some hard-wired assumptions.
    #
    # method : str `{'standard', 'fast'}`, optional
    #   'standard' (default) : Standard algorithm as described in [1, 2].
    #
    #   'fast' : Version II algorithm as described in [2].
    #
    # sdist : str `{'normal', 'unif'}`, optional
    #   'unif'  (default) : Uniform `[-1,1]`.
    #
    #   'norm' : Normal `~N(0,1)`.
    #
    # Value
    # -------
    # U:  array_like
    #   Right singular values, array of shape `(m, k)`.
    #
    # s : array_like
    #   Singular values, 1-d array of length `k`.
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
    rpcaObj = list(rotation = NULL,
                   eigvals = NULL,
                   sdev = NULL,
                   center=center,
                   scale=scale,
                   x=NULL)

    m <- nrow(A)
    n <- ncol(A)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Center/Scale data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    A <- scale( A , center = center, scale = scale )
    rpcaObj$center <- attributes(A)$`scaled:center`
    rpcaObj$scale <- attributes(A)$`scaled:scale`
    attributes(A)$`scaled:center`<- NULL
    attributes(A)$`scaled:scale`<- NULL

    rpcaObj$x <- A
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Compute randomized svd
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    svd_out = rsvd(A, k=k, p=p, q=q, ...)

    if(whiten==TRUE){
      svd_out$v <- svd_out$v / svd_out$d * sqrt(m)
    }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Explained variance and explained variance ratio
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rpcaObj$eigvals <- sqrt( svd_out$d )
    rpcaObj$sdev <- ( svd_out$d ) / sqrt( m-1 )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Add row and col names
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rownames(svd_out$v) <- colnames(A)
    colnames(svd_out$v) <- paste(rep('PC', length(svd_out$d)), 1:length(svd_out$d), sep = "")

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Return
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rpcaObj$rotation <- svd_out$v
    class(rpcaObj) <- "rpca"
    return( rpcaObj )

}


print.rpca <- function(rpcaObj) {
  cat("Standard deviations:\n")
  print(round(rpcaObj$sdev,3))
  cat("\nEigenvalues:\n")
  print(round(rpcaObj$eigvals,3))
  cat("\nRotation:\n")
  print(round(rpcaObj$rotation,3))
}


summary.rpca <- function( rpcaObj )
{

  variance = rpcaObj$sdev**2
  explained_var_ratio = variance / sum( apply( rpcaObj$x , 2, var ) )
  cum_explained_var_ratio = cumsum( explained_var_ratio )

  summaryObj <- t(data.frame( var = variance,
                              sdev = rpcaObj$sdev,
                              prob = explained_var_ratio,
                              cum= cum_explained_var_ratio,
                              eigv = rpcaObj$eigvals))

  rownames( summaryObj ) <- c('Explained variance',
                              'Standard deviations',
                              'Proportion of Variance',
                              'Cumulative Proportion',
                              'Eigenvalues')

  colnames( summaryObj ) <- paste(rep('PC', length(rpcaObj$sdev)), 1:length(rpcaObj$sdev), sep = "")

  return( summaryObj )
}


print.summary.rpca <- function( summaryObj )
{
  cat( "Importance of components:\n" )
  print(round(summaryObj,3))
  cat("\n")
}
