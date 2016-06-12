#' @title  Randomized principal component analysis (rPCA).
#
#' @description Principal components analysis using randomized singular value decomposition.
#
#' @details
#' Principal component analysis is a linear dimensionality reduction technique,
#' aiming to keep only the most significant principal components to allow
#' a better interpretation of the data and to project the data to a lower dimensional space.
#'
#' Traditionally, the computation is done by a (deterministic) singular value decomposition.
#' Randomized PCA is computed using a fast randomized algorithm (\code{\link{rsvd}})
#' to compute the approximate low-rank SVD decomposition.
#' The computational gain is high if the desired number of principal components
#' is small, i.e. \eqn{k << n}.
#'
#' \code{\link{rsvd}} expects a numeric (real/complex) input matrix with dimensions \eqn{(m, n)}.
#' Given a target rank \eqn{k}, \code{rsvd} factors the input matrix \eqn{A} as
#' \eqn{A = W * diag(s) * W'}. The columns of the real or complex unitary matrix \eqn{W}
#' contain the eigenvectors (i.e. principal components). The vector \eqn{s} contains the corresponding
#' eigenvalues. Following \code{\link{prcomp}} we denote this matrix \eqn{W} as
#' rotation matrix (commonly also called loadings).
#'
#' The print and summary method can be used to present the results in a nice format.
#' A scree plot can be produced with the plot function or as recommended with
#' \code{\link{ggscreeplot}}. A biplot can be produced with \code{\link{ggbiplot}},
#' and a correlation plot with \code{\link{ggcorplot}}.
#'
#' The predict function can be used to compute the scores of new observations. The data
#' will automatically be centred (and scaled if requested). This is not fully supported for
#' complex input matrices.
#'
#'
#' @param A       array_like \cr
#'                a numeric input matrix (or data frame), with dimensions \eqn{(m, n)}. \cr
#'                If the data contain \eqn{NA}s na.omit is applied.
#'
#' @param k       int, optional \cr
#'                determines the number of principle components to compute. It is required that \eqn{k} is smaller or equal to
#'                \eqn{n}, but it is recommended that \eqn{k << min(m,n)}.
#'
#' @param center  bool (\eqn{TRUE}, \eqn{FALSE}), optional \cr
#'                a logical value (\eqn{TRUE} by default) indicating whether the variables should be
#'                shifted to be zero centered. Alternatively, a vector of length equal the number of
#'                columns of \eqn{A} can be supplied. The value is passed to scale.
#'
#' @param scale   bool (\eqn{TRUE}, \eqn{FALSE}), optional \cr
#'                a logical value (\eqn{TRUE} by default) indicating whether the variables should
#'                be scaled to have unit variance. Alternatively, a vector of length equal the number of
#'                columns of \eqn{A} can be supplied. The value is passed to scale.
#'
#' @param retx    bool (\eqn{TRUE}, \eqn{FALSE}), optional \cr
#'                a logical value (\eqn{FALSE} by default) indicating whether the rotated variables / scores
#'                should be returned.
#'
#' @param loading   bool (\eqn{TRUE}, \eqn{FALSE}), optional \cr
#'                  When \eqn{TRUE} (by default \eqn{FALSE}) the eigenvectors
#'                  are unit scaled by the square root of the eigenvalues \eqn{W = W * diag(sqrt(eigvals))}.
#'
#' @param svdalg  str c('auto', 'rsvd', 'svd'), optional \cr
#'                Determines which algorithm should be used for computing the singular value decomposition.
#'                By default 'auto' is used, which decides whether to use \code{\link{rsvd}} or \code{\link{svd}},
#'                depending on the number of principle components. If \eqn{k < min(n,m)/1.5} randomized svd is used.
#'
#' @param p       int, optional \cr
#'                oversampling parameter for \eqn{rsvd}  (default \eqn{p=5}), see \code{\link{rsvd}}.
#'
#' @param q       int, optional \cr
#'                number of power iterations  for \eqn{rsvd} (default \eqn{q=2}), see \code{\link{rsvd}}.
#'
#' @param ...     arguments passed to or from other methods, see \code{\link{rsvd}}.
#'
#' @param ................. .
#'
#' @return \code{rpca} returns a list with class \eqn{rpca} containing the following components:
#'    \item{rotation}{  array_like \cr
#'                      matrix containing the rotation (eigenvectors),
#'                      or the variable loadings if \eqn{loadings=TRUE}; array with dimensions \eqn{(n, k)}.
#'    }
#'    \item{loading}{  array_like \cr
#'                      matrix containing the loadings (scaled eigenvectors),
#'                      if \eqn{loadings=TRUE}; array with dimensions \eqn{(n, k)}.
#'    }
#'
#'    \item{eigvals}{  array_like \cr
#'                     the eigenvalues; 1-d array of length \eqn{k}.
#'    }
#'    \item{sdev}{     array_like \cr
#'                     the standard deviations of the principal components.
#'    }
#'    \item{x}{        array_like \cr
#'                     if \eqn{retx} is true a matrix containing the scores / rotated data
#'                     (centred and scaled if requested) is returned.
#'    }
#'    \item{center, scale}{  array_like \cr
#'                     the centering and scaling used, or \eqn{FALSE}.
#'    }
#'    \item{.................}{.}
#'
#' @note  The principal components are not unique and only defined up to sign
#' (a constant of modulus one in the complex case) and so may differ between different
#'  PCA implementations.
#'
#' Similar to \code{\link{prcomp}} the variances are computed with the usual divisor N - 1.
#'
#' Note also that \eqn{scale = TRUE} cannot be used if there are zero or constant (for \eqn{center = TRUE} ) variables.
#'
#'
#' @author N. Benjamin Erichson, \email{nbe@st-andrews.ac.uk}
#'
#' @seealso \code{\link{ggscreeplot}}, \code{\link{ggbiplot}},
#'          \code{\link{ggcorplot}}, \code{\link{plot.rpca}},
#'          \code{\link{predict}},   \code{\link{rsvd}}
#'
#' @examples
#'
#'library(rsvd)
#'#
#'# Load Edgar Anderson's Iris Data
#'#
#'data(iris)
#'
#'#
#'# log transform
#'#
#'log.iris <- log( iris[ , 1:4] )
#'iris.species <- iris[ , 5]
#'
#'#
#'# Perform rPCA and compute all PCs, similar to \code{\link{prcomp}}
#'#
#'iris.rpca <- rpca(log.iris, retx=TRUE)
#'summary(iris.rpca) # Summary
#'print(iris.rpca) # Prints the rotations
#'
#'# You can compare the results with prcomp
#'# iris.pca <- prcomp(log.iris, center = TRUE, scale. = TRUE)
#'# summary(iris.pca) # Summary
#'# print(iris.pca) # Prints the rotations
#'
#'#
#'# Plot functions
#'#
#'ggscreeplot(iris.rpca) # Screeplot
#'ggscreeplot(iris.rpca, 'cum') # Screeplot
#'ggscreeplot(iris.rpca, type='eigenvals') # Screeplot of the eigenvalues
#'
#'ggcorplot(iris.rpca, pcs=c(1,2)) # The correlation of the original variable with the PCs
#'
#'ggbiplot(iris.rpca, groups = iris.species, circle = FALSE) #Biplot
#'
#'#
#'# Perform rPCA and compute only the first two PCs
#'#
#'iris.rpca <- rpca(log.iris, k=2,  svdalg = 'rsvd')
#'summary(iris.rpca) # Summary
#'print(iris.rpca) # Prints the rotations
#'
#'#
#'# Compute the scores of new observations
#'#
#'preds <- predict(iris.rpca, newdata=data.frame(log.iris))
#'


#' @export
rpca <- function(A, k=NULL, center=TRUE, scale=TRUE, loading=FALSE, retx=FALSE,  svdalg='auto', p=5, q=1, ...) UseMethod("rpca")

#' @export
rpca.default <- function(A, k=NULL, center=TRUE, scale=TRUE, loading=FALSE, retx=FALSE,  svdalg='auto', p=5, q=1, ...) {
    #*************************************************************************
    #***        Author: N. Benjamin Erichson <nbe@st-andrews.ac.uk>        ***
    #***                              <2015>                               ***
    #***                       License: BSD 3 clause                       ***
    #*************************************************************************
    rpcaObj = list(rotation = NULL,
                   loading = NULL,
                   eigvals = NULL,
                   sdev = NULL,
                   var = NULL,
                   center=center,
                   scale=scale,
                   x=NULL)

    m <- nrow(A)
    n <- ncol(A)

    #Set target rank
    if(is.null(k)) k=n
    if(k>n) k <- n
    if(k<1) stop("Target rank is not valid!")

    A <- stats::na.omit(A)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Center/Scale data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    A <- scale( A , center = center, scale = scale )
    rpcaObj$center <- attributes(A)$`scaled:center`
    rpcaObj$scale <- attributes(A)$`scaled:scale`
    attributes(A)$`scaled:center`<- NULL
    attributes(A)$`scaled:scale`<- NULL

    if(center == FALSE ) rpcaObj$center <- FALSE
    if(scale == FALSE ) rpcaObj$scale <- FALSE

    if(is.complex(A)) {
      rpcaObj$center <- center
      rpcaObj$scale <- scale
    }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Compute randomized svd
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(svdalg=='auto'){
      if(k < (n/1.5)) {svdalg='rsvd'} else svdalg='svd'
    }

    svd_out <- switch(svdalg,
                      svd = svd(A, nu = 0, nv = k),
                      rsvd = rsvd(A, k=k, p=p, q=q, nu = 0, ...),
                      stop("Selected SVD algorithm is not supported!")
    )


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Explained variance and explained variance ratio
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rpcaObj$singvals <- svd_out$d
    rpcaObj$eigvals <- svd_out$d**2 / (m-1)
    rpcaObj$sdev <- svd_out$d / sqrt( m-1 )
    rpcaObj$var <- sum( apply( Re(A) , 2, stats::var ) )
    if(is.complex(A)) rpcaObj$var <- Re(rpcaObj$var + sum( apply( Im(A) , 2, stats::var ) ))

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # loadings
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(loading==TRUE){
      rpcaObj$loading <- svd_out$v %*% diag(rpcaObj$eigvals**0.5)
    }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Add row and col names
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rownames(svd_out$v) <- colnames(A)
    colnames(svd_out$v) <- paste(rep('PC', length(svd_out$d)), 1:length(svd_out$d), sep = "")

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Compute rotated data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(retx==TRUE) rpcaObj$x <- A %*% svd_out$v

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Return
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rpcaObj$rotation <- svd_out$v
    class(rpcaObj) <- "rpca"
    return( rpcaObj )

}#End rPCA


#' @export
print.rpca <- function(x , ...) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Print rpca
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cat("Standard deviations:\n")
  print(round(x$sdev,3))
  cat("\nEigenvalues:\n")
  print(round(x$eigvals,3))
  cat("\nRotation:\n")
  print(round(x$rotation,3))
}

#' @export
summary.rpca <- function( object , ... )
{
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Summary rpca
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  variance = object$sdev**2
  explained_variance_ratio = variance / object$var
  cum_explained_variance_ratio = cumsum( explained_variance_ratio )

  x <- t(data.frame( var = round(variance, 3),
                              sdev = round(object$sdev, 3),
                              prob = round(explained_variance_ratio, 3),
                              cum = round(cum_explained_variance_ratio, 3),
                              eigv = round(object$eigvals, 3)))

  rownames( x ) <- c( 'Explained variance',
                      'Standard deviations',
                      'Proportion of variance',
                      'Cumulative proportion',
                      'Eigenvalues')

  colnames( x ) <- paste(rep('PC', length(object$sdev)), 1:length(object$sdev), sep = "")

  x <- as.matrix(x)

  return( x )
}

#' @export
print.summary.rpca <- function( x , ... )
{
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Print summary rpca
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cat( "Importance of components:\n" )
  print( x )
  cat("\n")
}

#' @export
predict.rpca <- function( object, newdata, ...)
{
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Predict
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  newdata <- scale(newdata, center = object$center, scale = object$scale)
  x <- as.matrix(newdata) %*% as.matrix(object$rotation)

  return( x )
}
