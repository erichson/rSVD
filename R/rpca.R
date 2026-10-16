#' @title  Randomized principal component analysis (rpca).
#
#' @description Fast computation of the principal components analysis using the randomized singular value decomposition.
#
#' @details
#' Principal component analysis is an important linear dimension reduction technique.
#'
#' Randomized PCA is computed via the randomized SVD algorithm (\code{\link{rsvd}}).
#' The computational gain is substantial, if the desired number of principal components
#' is relatively small, i.e. \eqn{k << min(m,n)}.
#'
#' The print and summary method can be used to present the results in a nice format.
#' A scree plot can be produced with \code{\link{ggscreeplot}}. 
#' The individuals factor map can be produced with \code{\link{ggindplot}},
#' and a correlation plot with \code{\link{ggcorplot}}.
#'
#' The predict function can be used to compute the scores of new observations. The data
#' will automatically be centered (and scaled if requested). This is not fully supported for
#' complex input matrices.
#'
#'
#' @param A      array_like; \cr
#'               a numeric \eqn{(m, n)} input matrix (or data frame) to be analyzed. \cr
#'               If the data contain \eqn{NA}s na.omit is applied.
#'
#' @param k      integer; \cr
#'               number of dominant principle components to be computed. It is required that \eqn{k} is smaller or equal to
#'               \eqn{min(m,n)}, but it is recommended that \eqn{k << min(m,n)}.
#'
#' @param center  bool, optional; \cr
#'                logical value which indicates whether the variables should be
#'                shifted to be zero centered (\eqn{TRUE} by default).
#'
#' @param scale   bool, optional; \cr
#'                logical value which indicates whether the variables should
#'                be scaled to have unit variance (\eqn{TRUE} by default).
#'
#' @param retx    bool, optional; \cr
#'                logical value indicating whether the rotated variables / scores
#'                should be returned (\eqn{TRUE} by default).
#'
#' @param p       integer, optional; \cr
#'                oversampling parameter for \eqn{rsvd} (default \eqn{p=10}), see \code{\link{rsvd}}.
#'
#' @param q       integer, optional; \cr
#'                number of additional power iterations for \eqn{rsvd} (default \eqn{q=1}), see \code{\link{rsvd}}.
#'
#' @param rand   bool, optional; \cr
#'               if (\eqn{TRUE}), the \eqn{rsvd} routine is used, otherwise \eqn{svd} is used.
#'
#'
#' @return \code{rpca} returns a list with class \eqn{rpca} containing the following components:
#' \describe{
#'    \item{rotation}{  array_like; \cr
#'                      the rotation (eigenvectors); \eqn{(n, k)} dimensional array.
#'    }
#'
#'    \item{eigvals}{  array_like; \cr
#'                     eigenvalues; \eqn{k} dimensional vector.
#'    }
#'    \item{sdev}{     array_like; \cr
#'                     standard deviations of the principal components; \eqn{k} dimensional vector.
#'    }
#'    \item{x}{        array_like; \cr
#'                     the scores / rotated data; \eqn{(m, k)} dimensional array.
#'    }
#'    \item{center, scale}{  array_like; \cr
#'                     the centering and scaling used.
#'    }
#'}
#'
#'
#' @note  The principal components are not unique and only defined up to sign
#' (a constant of modulus one in the complex case) and so may differ between different
#'  PCA implementations.
#'
#' Similar to \code{\link{prcomp}} the variances are computed with the usual divisor N - 1.
#'
#'
#' @references
#' \itemize{
#'  \item [1] N. B. Erichson, S. Voronin, S. L. Brunton and J. N. Kutz. 2019.
#'        Randomized Matrix Decompositions Using {R}. 
#'        Journal of Statistical Software, 89(11), 1-48.
#'         \doi{10.18637/jss.v089.i11}.
#' 
#'   \item  [2] N. Halko, P. Martinsson, and J. Tropp.
#'          "Finding structure with randomness: probabilistic
#'          algorithms for constructing approximate matrix
#'          decompositions" (2009).
#'          (available at arXiv \url{https://arxiv.org/abs/0909.4061}).
#'}
#'
#'
#' @author N. Benjamin Erichson, \email{erichson@berkeley.edu}
#'
#' @seealso \code{\link{ggscreeplot}}, \code{\link{ggindplot}},
#'          \code{\link{ggcorplot}}, \code{\link{plot.rpca}},
#'          \code{\link{predict}},   \code{\link{rsvd}}
#'
#' @examples
#'
#'library('rsvd')
#'#
#'# Load Edgar Anderson's Iris Data
#'#
#'data('iris')
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
    rpcaObj$x <- sweep(out$u[, 1:k, drop=FALSE], MARGIN = 2, STATS = out$d[1:k], FUN = "*", check.margin = TRUE)
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

