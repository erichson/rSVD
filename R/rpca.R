#' @title  Randomized principal component analysis (rpca).
#
#' @description Principal components analysis using the randomized singular value decomposition.
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
#' is small, i.e. \eqn{k << min(m,n)}.
#'
#' \code{\link{rsvd}} expects a numeric (real/complex) input matrix with dimensions \eqn{(m, n)}.
#' Given a target rank \eqn{k}, \code{rsvd} factors the input matrix \eqn{A} as
#' \eqn{A = W * diag(s) * W'}. The columns of the real or complex unitary matrix \eqn{W}
#' contain the eigenvectors (i.e. principal components). The vector \eqn{s} contains the corresponding
#' eigenvalues. Following \code{\link{prcomp}} we denote this matrix \eqn{W} as
#' rotation matrix (commonly also called loadings).
#'
#' The print and summary method can be used to present the results in a nice format.
#' A scree plot can be produced with
#' \code{\link{ggscreeplot}}. The individuals factor map can be produced with \code{\link{ggindplot}},
#' and a correlation plot with \code{\link{ggcorplot}}.
#'
#' The predict function can be used to compute the scores of new observations. The data
#' will automatically be centred (and scaled if requested). This is not fully supported for
#' complex input matrices.
#'
#'
#' @param A       Array_like. \cr
#'                A numeric input matrix (or data frame), with dimensions \eqn{(m, n)}. \cr
#'                If the data contain \eqn{NA}s na.omit is applied.
#'
#' @param k       Int, optional. \cr
#'                Determines the number of principle components to compute. It is required that \eqn{k} is smaller or equal to
#'                \eqn{n}, but it is recommended that \eqn{k << min(m,n)}.
#'
#' @param center  Bool (\eqn{TRUE}, \eqn{FALSE}), optional. \cr
#'                A logical value (\eqn{TRUE} by default) indicating whether the variables should be
#'                shifted to be zero centered. Alternatively, a vector of length equal the number of
#'                columns of \eqn{A} can be supplied. The value is passed to scale.
#'
#' @param scale   Bool (\eqn{TRUE}, \eqn{FALSE}), optional. \cr
#'                A logical value (\eqn{TRUE} by default) indicating whether the variables should
#'                be scaled to have unit variance. Alternatively, a vector of length equal the number of
#'                columns of \eqn{A} can be supplied. The value is passed to scale.
#'
#' @param retx    Bool (\eqn{TRUE}, \eqn{FALSE}), optional \cr
#'                A logical value (\eqn{TRUE} by default) indicating whether the rotated variables / scores
#'                should be returned.
#'
#' @param p       Int, optional. \cr
#'                Oversampling parameter for \eqn{rsvd}  (default \eqn{p=10}), see \code{\link{rsvd}}.
#'
#' @param q       Int, optional. \cr
#'                Number of power iterations  for \eqn{rsvd} (default \eqn{q=1}), see \code{\link{rsvd}}.
#'
#' @param rand  Bool (\eqn{TRUE}, \eqn{FALSE}). \cr
#'              If (\eqn{TRUE}), a probabilistic strategy is used, otherwise a deterministic algorithm is used.
#'
#' @param ................. .
#'
#' @return \code{rpca} returns a list with class \eqn{rpca} containing the following components:
#'    \item{rotation}{  Array_like. \cr
#'                      Matrix containing the rotation (eigenvectors); array with dimensions \eqn{(n, k)}.
#'    }
#'
#'    \item{eigvals}{  Array_like. \cr
#'                     The eigenvalues; 1-d array of length \eqn{k}.
#'    }
#'    \item{sdev}{     Array_like \cr
#'                     The standard deviations of the principal components.
#'    }
#'    \item{x}{        Array_like \cr
#'                     If \eqn{retx} is true a matrix containing the scores / rotated data
#'                     (centred and scaled if requested) is returned.
#'    }
#'    \item{center, scale}{  Array_like .\cr
#'                     The centering and scaling used, or \eqn{FALSE}.
#'    }
#'    \item{.................}{.}
#'
#' @note  The principal components are not unique and only defined up to sign
#' (a constant of modulus one in the complex case) and so may differ between different
#'  PCA implementations.
#'
#' Similar to \code{\link{prcomp}} the variances are computed with the usual divisor N - 1.
#'
#' Note, that \eqn{scale = TRUE} cannot be used, if some columns have only constant entries.
#'
#'
#' @author N. Benjamin Erichson, \email{erichson@uw.edu}
#'
#' @seealso \code{\link{ggscreeplot}}, \code{\link{ggindplot}},
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
#'# Perform rPCA and compute only the first two PCs
#'#
#'iris.rpca <- rpca(log.iris, k=2)
#'summary(iris.rpca) # Summary
#'print(iris.rpca) # Prints the rotations
#'
#'#
#'# Use rPCA to compute all PCs, similar to \code{\link{prcomp}}
#'#
#'iris.rpca <- rpca(log.iris)
#'summary(iris.rpca) # Summary
#'print(iris.rpca) # Prints the rotations
#'plot(iris.rpca) # Produce screeplot, variable and individuls factor maps.
#'
#'
#'# You can compare the results with prcomp
#'# iris.pca <- prcomp(log.iris, center = TRUE, scale. = TRUE)
#'# summary(iris.pca) # Summary
#'# print(iris.pca) # Prints the rotations

#' @export
rpca <- function(A, k=NULL, center=TRUE, scale=TRUE, retx=TRUE, p=10, q=2, rand = TRUE) UseMethod("rpca")

#' @export
rpca.default <- function(A, k=NULL, center=TRUE, scale=TRUE, retx=TRUE, p=10, q=2, rand = TRUE) {

    A <- as.matrix(A)
  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Checks
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (any(is.na(A))) {
      warning("Missing values are omitted: na.omit(A).")
      A <- stats::na.omit(A)
    }   

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Init rpca object
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rpcaObj = list(rotation = NULL,
                   eigvals = NULL,
                   sdev = NULL,
                   var = NULL,
                   center = center,
                   scale = scale,
                   x=NULL)

    m <- nrow(A)
    n <- ncol(A)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Set target rank
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(is.null(k)) rand <- FALSE
    if(is.null(k)) k <- min(n,m)
    if(k > min(n,m)) k <- min(n,m)
    if(k<1) stop("Target rank is not valid!")

   

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Center/Scale data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(center == TRUE) {
      rpcaObj$center <- colMeans(A)
      A <- sweep(A, MARGIN = 2, STATS = rpcaObj$center, FUN = "-", check.margin = TRUE)
      #A <- H(H(A) - rpcaObj$center)
    } else { rpcaObj$center <- FALSE }
    
    if(scale == TRUE) {
      rpcaObj$scale <- sqrt(colSums(A**2) / (m-1))
      if(is.complex(rpcaObj$scale)) { rpcaObj$scale[Re(rpcaObj$scale) < 1e-8 ] <- 1+0i  
      } else {rpcaObj$scale[rpcaObj$scale < 1e-8] <- 1}
      A <- sweep(A, MARGIN = 2, STATS = rpcaObj$scale, FUN = "/", check.margin = TRUE)
      #A <- H(H(A) / rpcaObj$scale)
      #A[is.nan(A)] <- 0
    } else { rpcaObj$scale <- FALSE }
    
    

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Compute randomized svd / eigen decomposition
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(rand == TRUE) {
      svdalg = 'rsvd'
    }else { 
      svdalg = 'svd' 
    }

    out <- switch(svdalg,
                      svd = svd(A, nu = k, nv = k),
                      rsvd = rsvd(A, k = k, p = p, q = q),
                      stop("Selected SVD algorithm is not supported!")
    )

    rpcaObj$eigvals <- switch(svdalg,
                              svd = out$d[1:k]**2 / (m-1),
                              rsvd = out$d**2 / (m-1)
    )

    rpcaObj$rotation <- switch(svdalg,
                              svd = out$v,
                              rsvd = out$v
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Explained variance and explained variance ratio
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rpcaObj$sdev <-  sqrt( rpcaObj$eigvals )
    rpcaObj$var <- sum( apply( Re(A) , 2, stats::var ) )
    if(is.complex(A)) rpcaObj$var <- Re(rpcaObj$var + sum( apply( Im(A) , 2, stats::var ) ))

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Add row and col names
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rownames(rpcaObj$rotation) <- colnames(A)
    colnames(rpcaObj$rotation) <- paste(rep('PC', k), 1:k, sep = "")

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Compute rotated data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(retx==TRUE) {
      #rpcaObj$x <- A %*% rpcaObj$rotation # slow
      #rpcaObj$x <- H(H(out$u[,1:k]) * out$d[1:k])
      rpcaObj$x <- sweep(out$u[,1:k], MARGIN = 2, STATS = out$d[1:k], FUN = "*", check.margin = TRUE)
      rownames(rpcaObj$x) <- rownames(A)
    }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Return
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    class(rpcaObj) <- "rpca"
    return( rpcaObj )

}#End rPCA


#' @export
print.rpca <- function(x , ...) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Print rpca
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cat("Standard deviations:\n")
  print(round(x$sdev, 3))
  cat("\nEigenvalues:\n")
  print(round(x$eigvals, 3))
  cat("\nRotation:\n")
  print(round(x$rotation, 3))
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
                              cum = round(cum_explained_variance_ratio, 3)))

  rownames( x ) <- c( 'Explained variance',
                      'Standard deviations',
                      'Proportion of variance',
                      'Cumulative proportion')

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
  if(!is.logical(object$center)) {
    #newdata <- H(H(newdata) - object$center)
    newdata <- sweep(newdata, MARGIN = 2, STATS = object$center, FUN = "-", check.margin = TRUE)
  }
  
  if(!is.logical(object$scale)) {
    #newdata <- H(H(newdata) / object$scale)
    newdata <- sweep(newdata, MARGIN = 2, STATS = object$scale, FUN = "/", check.margin = TRUE)
    newdata[is.nan(newdata)] <- 0
  }
  
  x <- as.matrix(newdata) %*% as.matrix(object$rotation)
  rownames(x) <- rownames(newdata)
  
  return( x )
}
