#' @title  Randomized principal component analysis (PCA).
#
#' @description Performs a approximated principal components analysis using randomized SVD \code{rsvd}.
#
#' @details
#' Principle component analysis is a linear dimensionality reduction technique,
#' aiming to keep only the most significant principle components to allow a
#' a better interpretation of the data and to project the data to a lower dimensional space.
#'
#' Traditionally, the computation is done by a (deterministic) singular value decomposition,
#' randomized PCA is using a very fast randomized SVD algorithm \code{\link{rsvd}}
#' to compute the approximate low-rank SVD decomposition, instead.
#' The computational gain is in particular high, if the desired number of principle components
#' is small, i.e. \eqn{k << n}.
#'
#' \code{\link{rsvd}} expects a numeric (real/complex) input matrix with dimensions \eqn{(m, n)}.
#' Given a target rank \eqn{k}, \code{rsvd} factors the input matrix \eqn{A} as
#' \eqn{A = W * diag(s) * W'}. The the columns of the real or complex unitary matrix \eqn{W}
#' contain the eigenvectors (i.e. principle components). The vector \eqn{s} contains the corresponding
#' eigenvalues. Following \code{\link{prcomp}} we denote this matrix \eqn{W} as
#' rotation matrix (but also called loadings).
#'
#' The print and summary method can be used to present the results in a nice format.
#' A scree plot can be produced with the plot function or as rommended with
#' \code{\link{ggscreeplot}}. A biplot can be proeduced with \code{\link{ggbiplot}},
#' and a rotation plot with \code{\link{ggrotplot}}.
#'
#' The predict function can be used to compute the scores of new observations. The data
#' will automatically be centred (and scaled if requested).
#'
#'
#' @param A       array_like \cr
#'                a real/complex input matrix (or data frame), with dimensions \eqn{(m, n)}. \cr
#'                If the data contain \eqn{NA}s na.omit is applied.
#'
#' @param k       int, optional \cr
#'                determines the number of principle components to compute. It is required that \eqn{k} is smaller or equal to
#'                \eqn{min(n)}, but it is recommended that \eqn{k << min(m,n)}.
#'
#' @param center  bool (\eqn{TRUE}, \eqn{FALSE}), optional \cr
#'                a logical value (\eqn{TRUE} by default) indicating whether the variables should be
#'                shifted to be zero centered. Alternately, a vector of length equal the number of
#'                columns of \eqn{A} can be supplied. The value is passed to scale.
#'
#' @param scale   bool (\eqn{TRUE}, \eqn{FALSE}), optional \cr
#'                a logical value (\eqn{TRUE} by default) indicating whether the variables should
#'                be scaled to have unit variance.Alternately, a vector of length equal the number of
#'                columns of \eqn{A} can be supplied. The value is passed to scale.
#'
#' @param retx    bool (\eqn{TRUE}, \eqn{FALSE}), optional \cr
#'                a logical value (\eqn{FALSE} by default) indicating whether the rotated variables / scores
#'                should be returned.
#'
#' @param whiten  bool (\eqn{TRUE}, \eqn{FALSE}), optional \cr
#'                When True (\eqn{FALSE} by default) the eigenvectors
#'                are divided by the the square root of the singular values \eqn{W = W * diag(1/sqrt(s))}.
#'                Whitening can sometime improve the predictive accuracy.
#'
#' @param svdalg  str c('auto', 'rsvd', 'svd'), optional \cr
#'                Determines which algorithm should be used for computing the singular value decomposition.
#'                By default 'auto' is used, which decides weather to use \code{\link{rsvd}} or \code{\link{svd}},
#'                depending on the number of principle components. If \eqn{k < min(n,m)/1.5} randomized svd is used.
#'
#' @param ...     arguments passed to or from other methods, see \code{\link{rsvd}}.
#'
#' @param ................. .
#'
#' @return \code{rpca} returns a list with class \eqn{rpca} containing the following components:
#'    \item{rotation}{  array_like \cr
#'                      the matrix containing the roations (eigenvectors),
#'                      i.e., the variable loadings; array with dimensions \eqn{(n, k)}.
#'    }
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
#' @note  The principle components are not unique and only defined up to sign
#' (a constant of modulus one in the complex case) and so may differ between different
#'  PCA implementations.
#'
#' Similar to \code{\link{prcomp}} the variances are computed with the usual divisor N - 1.
#' using approximated Singular Value
#'
#' Note also that \eqn{scale = TRUE} cannot be used if there are zero or constant (for \eqn{center = TRUE} ) variables.
#'
#'
#' @author N. Benjamin Erichson, \email{nbe@st-andrews.ac.uk}
#'
#' @seealso \code{\link{ggscreeplot}}, \code{\link{ggbiplot}},
#'          \code{\link{ggrotplot}}, \code{\link{plot}},
#'          \code{\link{predict}},   \code{\link{rsvd}}
#'
#' @examples
#'
#'library()
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
#'iris.rpca <- rpca(log.iris, retx=T,  svdalg = 'rsvd')
#'summary(iris.rpca) # Summary
#'print(iris.rpca) # Prints the loadings/ rotations
#'
#'# You can compare the results with prcomp
#'# iris.pca <- prcomp(log.iris, center = T, scale. = T)
#'# summary(iris.pca) # Summary
#'# print(iris.pca) # Prints the loadings/ rotations
#'
#'#
#'# Plot functions
#'#
#'ggscreeplot(iris.rpca) # Screeplot
#'ggscreeplot(iris.rpca, 'cum_explained_var_ratio') # Screeplot
#'ggscreeplot(ir.rpca, type='raw') # Screeplot of the eigenvalues
#'
#'ggrotplot(iris.rpca) # The correlation of the original variable with the PCs
#'
#'ggbiplot(iris.rpca, groups = ir.species, circle = F) #Biplot
#'
#'#
#'# Perform rPCA and compute only the first two PCs
#'#
#'iris.rpca <- rpca(log.iris, k=2,  svdalg = 'rsvd')
#'summary(iris.rpca) # Summary
#'print(iris.rpca) # Prints the loadings/ rotations
#'
#'#
#'# Compute the scores of new observations
#'#
#'preds <- predict(iris.rpca, newdata=data.frame(log.iris))
#'


rpca <- function(A, k=NULL, center=TRUE, scale=TRUE, whiten=FALSE, retx=FALSE, svdalg='auto', p=5, q=2, ...) UseMethod("rpca")

rpca.default <- function(A, k=NULL, center=TRUE, scale=TRUE, whiten=FALSE, retx=FALSE,  svdalg='auto', p=5, q=2, ...) {
    #*************************************************************************
    #***        Author: N. Benjamin Erichson <nbe@st-andrews.ac.uk>        ***
    #***                              <2015>                               ***
    #***                       License: BSD 3 clause                       ***
    #*************************************************************************
    rpcaObj = list(rotation = NULL,
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

    A <- na.omit(A)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Center/Scale data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    A <- scale( A , center = center, scale = scale )
    rpcaObj$center <- attributes(A)$`scaled:center`
    rpcaObj$scale <- attributes(A)$`scaled:scale`
    attributes(A)$`scaled:center`<- NULL
    attributes(A)$`scaled:scale`<- NULL


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

    if(whiten==TRUE){
      svd_out$v <- svd_out$v / sqrt(svd_out$d)
    }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Explained variance and explained variance ratio
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rpcaObj$eigvals <- sqrt( svd_out$d )
    rpcaObj$sdev <- ( svd_out$d ) / sqrt( m-1 )
    rpcaObj$var <- sum( apply( Re(A) , 2, var ) )
    if(is.complex(A)) rpcaObj$var <- Re(rpcaObj$var + sum( apply( Im(A) , 2, var ) ))

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



print.rpca <- function(rpcaObj) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Print rpca
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cat("Standard deviations:\n")
  print(round(rpcaObj$sdev,3))
  cat("\nEigenvalues:\n")
  print(round(rpcaObj$eigvals,3))
  cat("\nRotation:\n")
  print(round(rpcaObj$rotation,3))
}


summary.rpca <- function( rpcaObj, ... )
{
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Summary rpca
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  variance = rpcaObj$sdev**2
  explained_var_ratio = variance / rpcaObj$var
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

  summaryObj <- as.matrix(summaryObj)

  return( summaryObj )
}


print.summary.rpca <- function( summaryObj, ... )
{
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Print summary rpca
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cat( "Importance of components:\n" )
  print( summaryObj )
  cat("\n")
}


predict.rpca <- function( rpcaObj, newdata )
{
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Predict
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  newdata <- scale(newdata, center = rpcaObj$center, scale = rpcaObj$scale)
  x <- as.matrix(newdata) %*% as.matrix(rpcaObj$rotation)

}
