#' @title Randomized principal component analysis
#'
#' @description Compute principal component loadings, scores, and variances from an
#' optionally centered and scaled matrix, using a randomized or deterministic
#' singular value decomposition (SVD).
#'
#' @details The randomized algorithm uses \code{\link{rsvd}} to approximate the
#' dominant singular subspace. Its computational advantage is greatest for
#' large matrices and a small requested rank. The retained scores and
#' loadings describe a low-rank approximation to the processed data. The
#' \code{p} and \code{q} arguments trade additional work for approximation
#' accuracy; use \code{set.seed} before fitting for reproducible random draws.
#'
#' Both centering and scaling are enabled by default. As in
#' \code{\link[stats]{prcomp}}, \code{center = FALSE} simply skips centering;
#' it supports both already-centered and uncentered data. With
#' \code{scale = TRUE}, each column is divided by the square root of its sum
#' of squared magnitudes divided by \eqn{m-1}, after any requested centering.
#' This is the sample standard deviation for centered data. Scale factors
#' below \code{1e-8} are replaced by one for compatibility, so constant and
#' nearly constant columns are retained. Complex columns use squared
#' magnitudes, giving real scale factors. Unlike the default in
#' \code{prcomp}, scaling is enabled by default here.
#'
#' Total variance is the squared Frobenius norm of the processed matrix
#' divided by \eqn{m-1}, including directions omitted from the approximation.
#' For uncentered data it is squared energy divided by \eqn{m-1} and includes
#' the contribution of the column means. At least two complete observations
#' are needed for meaningful sample variances. Variance
#' proportions are undefined when the processed matrix has zero total
#' variance. The proportion attributed to a retained component is its
#' \code{eigvals} entry divided by \code{var}, rather than divided by the sum
#' of the retained eigenvalues.
#'
#' \code{predict(object, newdata)} applies the stored training centering and
#' scaling, then multiplies the processed new observations by
#' \code{object$rotation}. It does not estimate new means or scale factors,
#' reorder columns by name, or omit incomplete new rows automatically.
#' Prediction also works when the fit was created with \code{retx = FALSE}.
#' Fitted scores retain the projected-SVD convention \eqn{U_k D_k}{U_k D_k}.
#' For a randomized approximation, projecting the original training data
#' with \code{predict} can differ from \code{object$x} by the range-finder
#' residual. To reconstruct observations in the input variable coordinates,
#' multiply fitted scores by the conjugate transpose of the loadings
#' and then undo any scaling and centering. Scores alone are expressed in
#' component coordinates.
#'
#' With \code{rand = TRUE} and a supplied \code{k}, real sparse inputs are
#' centered and scaled through implicit matrix products. The entire centered
#' matrix is never formed; storage includes the sparse entries and dense
#' working factors of width up to \code{k + p}, capped by the matrix
#' dimensions. Loadings, fitted scores, and predicted scores are dense.
#' With \code{rand = FALSE} or \code{k = NULL}, a sparse input is converted
#' to a dense matrix for the deterministic SVD. This is ordinary PCA of
#' sparse input, not a method imposing sparsity on the loadings.
#'
#' The \code{print} and \code{summary} methods describe a fit.
#' \code{\link{plot.rpca}} draws the available scree and factor plots;
#' the individual plotting helpers can also be called directly.
#'
#' @param A A numeric matrix or data frame with observations in rows and variables
#' in columns. Dense real or complex data and real sparse matrices from
#' \pkg{Matrix} are supported. Rows containing missing values are omitted
#' with a warning. In the dimensions below, \eqn{m} is the number of retained
#' observations and \eqn{n} is the number of variables.
#'
#' @param k A positive integer giving the number of principal components to return.
#' Values greater than \code{min(m, n)} are reduced to that limit. The default,
#' \code{NULL}, returns \code{min(m, n)} components using the deterministic
#' SVD. The randomized method is most useful when \code{k} is much smaller
#' than both matrix dimensions.
#'
#' @param center Logical; subtract the training column means before computing the SVD
#' (default \code{TRUE}). \code{FALSE} leaves the column means unchanged.
#'
#' @param scale Logical; divide each processed column by its scale factor (default
#' \code{TRUE}). This gives unit sample variance for centered columns whose
#' scale factors are at least \code{1e-8}; see Details.
#'
#' @param retx Logical; return the fitted component scores in \code{x} (default
#' \code{TRUE}). With \code{FALSE}, \code{x} is \code{NULL}; loadings and
#' new-data prediction remain available.
#'
#' @param p A nonnegative integer oversampling parameter for \code{\link{rsvd}}
#' (default \code{10}). Ignored for the deterministic SVD.
#'
#' @param q A nonnegative integer number of power iterations for
#' \code{\link{rsvd}} (default \code{2}). Ignored for the deterministic SVD.
#'
#' @param rand Logical; use the randomized SVD when \code{TRUE} and \code{k} is
#' supplied. With \code{FALSE}, or \code{k = NULL}, use the deterministic
#' SVD. Requesting fewer components does not make the deterministic backend
#' a partial SVD algorithm.
#'
#' @param object,x An object returned by \code{rpca}.
#'
#' @param newdata A numeric matrix, data frame, or real sparse \pkg{Matrix} containing
#' new observations. Supply the same variables in the same column order as
#' the training data, with no missing values; see Details.
#'
#' @param ... Additional arguments; currently ignored by the \code{print},
#' \code{summary}, and \code{predict} methods.
#'
#' @return \code{rpca} returns a list of class \code{"rpca"} with these components:
#' \describe{
#'   \item{rotation}{A dense \eqn{n} by \eqn{k} matrix of orthonormal
#'     loadings (right singular vectors of the processed data).}
#'   \item{eigvals}{A numeric vector of length \eqn{k} containing the
#'     retained squared singular values divided by \eqn{m-1}. These estimate
#'     the leading eigenvalues of the processed cross-product matrix divided
#'     by \eqn{m-1}, which is the sample covariance matrix for centered data.}
#'   \item{sdev}{A numeric vector of length \eqn{k}: the square roots of
#'     \code{eigvals}.}
#'   \item{var}{The squared Frobenius norm of the processed matrix divided
#'     by \eqn{m-1}, including directions not returned. For centered data
#'     this is the total sample variance.}
#'   \item{x}{A dense \eqn{m} by \eqn{k} matrix of fitted scores, or
#'     \code{NULL} when \code{retx = FALSE}. Rows correspond to complete
#'     observations retained from the input.}
#'   \item{center, scale}{Numeric vectors containing the training column
#'     means and scale factors, respectively, or \code{FALSE} when the
#'     corresponding operation was disabled. Centers can be complex;
#'     scale factors are real.}
#' }
#' The \code{predict} method returns a dense matrix with one row per new
#' observation and one column per retained component. The \code{summary}
#' method returns a numeric matrix with rows for explained variance, standard
#' deviation, proportion of variance, and cumulative proportion. Its entries
#' are rounded to three decimal places. The \code{print} method displays
#' standard deviations, eigenvalues, and loadings, and invisibly returns the
#' displayed loading matrix.
#'
#' @note Singular vectors are defined only up to sign (or a constant of modulus
#' one for complex data). With repeated singular values, any orthonormal
#' basis of the corresponding subspace is valid. Individual loading vectors
#' can therefore differ between runs or implementations without changing
#' the represented subspace.
#'
#' @references \itemize{
#'   \item N. B. Erichson, S. Voronin, S. L. Brunton and J. N. Kutz. 2019.
#'   Randomized Matrix Decompositions Using R.
#'   Journal of Statistical Software, 89(11), 1--48.
#'   \doi{10.18637/jss.v089.i11}.
#'   \item N. Halko, P. G. Martinsson and J. A. Tropp. 2011.
#'   Finding Structure with Randomness: Probabilistic Algorithms for
#'   Constructing Approximate Matrix Decompositions.
#'   SIAM Review, 53(2), 217--288.
#'   \doi{10.1137/090771806}.
#' }
#'
#' @seealso \code{\link{rsvd}}, \code{\link[stats]{prcomp}},
#' \code{\link{plot.rpca}}, \code{\link{ggscreeplot}},
#' \code{\link{ggindplot}}, \code{\link{ggcorplot}},
#' \code{\link{rrpca}} for a low-rank plus sparse decomposition.
#'
#' @author N. Benjamin Erichson
#'
#' @examples # Fit two components and project observations held out from the fit.
#' training <- log(as.matrix(iris[1:145, 1:4]))
#' new_data <- log(as.matrix(iris[146:150, 1:4]))
#' set.seed(42)
#' fit <- rpca(training, k = 2)
#' summary(fit)
#' predict(fit, new_data)
#'
#' # Reconstruct in input coordinates (here, log measurements).
#' reconstructed <- fit$x %*% Conj(t(fit$rotation))
#' reconstructed <- sweep(reconstructed, 2, fit$scale, "*")
#' reconstructed <- sweep(reconstructed, 2, fit$center, "+")
#'
#' # k = NULL uses the deterministic SVD and returns all components.
#' full <- rpca(training)
#' length(full$sdev)
#'
#' # Sparse input: implicit preprocessing, with training scores omitted.
#' A <- Matrix::bdiag(outer(1:12, 1:6), outer(1:12, 6:1))
#' set.seed(42)
#' sparse_fit <- rpca(A, k = 2, retx = FALSE)
#' predict(sparse_fit, A[1:2, , drop = FALSE])
#'
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot(fit, type = "ratio", groups = iris$Species[1:145],
#'        ind_labels = FALSE, ellipse = FALSE)
#' }
#'
#' @export
rpca <- function(A, k=NULL, center=TRUE, scale=TRUE, retx=TRUE, p=10, q=2, rand = TRUE) UseMethod("rpca")

#' @export
rpca.default <- function(A, k=NULL, center=TRUE, scale=TRUE, retx=TRUE, p=10, q=2, rand = TRUE) {
  
  sparse <- inherits(A, "sparseMatrix") && !is.null(k) && rand == TRUE
  if(!sparse) A <- as.matrix(A)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Checks
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (any(is.na(A))) {
    warning("Missing values are omitted: na.omit(A).")
    if(sparse) {
      A <- A[Matrix::rowSums(is.na(A)) == 0, , drop=FALSE]
    } else A <- stats::na.omit(A)
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
  if(sparse) {
    A <- rpca_sparse_numeric(A)
    rpcaObj$center <- if(center == TRUE) rpca_sparse_means(A) else FALSE
    squared_norms <- rpca_sparse_norms(A, rpcaObj$center)
    if(scale == TRUE) {
      rpcaObj$scale <- sqrt(squared_norms / (m-1))
      rpcaObj$scale[rpcaObj$scale < 1e-8] <- 1
      rpcaObj$var <- sum(squared_norms / rpcaObj$scale^2) / (m-1)
    } else {
      rpcaObj$scale <- FALSE
      rpcaObj$var <- sum(squared_norms) / (m-1)
    }
    operators <- rpca_sparse_products(A, rpcaObj$center, rpcaObj$scale)
    if(!is.finite(rpcaObj$var)) {
      rpcaObj$var <- sum(operators$squared_norms()) / (m-1)
    }
  } else {
    if(center == TRUE) {
      rpcaObj$center <- colMeans(A)
      A <- sweep(A, MARGIN = 2, STATS = rpcaObj$center, FUN = "-", check.margin = TRUE)
    } else { rpcaObj$center <- FALSE }

    if(scale == TRUE) {
      squared <- if(is.complex(A)) Mod(A)**2 else A**2
      squared_norms <- colSums(squared)
      rm(squared)
      rpcaObj$scale <- sqrt(squared_norms / (m-1))
      rpcaObj$scale[rpcaObj$scale < 1e-8] <- 1
      A <- sweep(A, MARGIN = 2, STATS = rpcaObj$scale, FUN = "/", check.margin = TRUE)
      rpcaObj$var <- sum(squared_norms / rpcaObj$scale^2) / (m-1)
      if(!is.finite(rpcaObj$var)) {
        rpcaObj$var <- (if(is.complex(A)) {
          norm(Re(A), "F")^2 + norm(Im(A), "F")^2
        } else norm(A, "F")^2) / (m-1)
      }
    } else {
      rpcaObj$scale <- FALSE
      squared_norm <- if(is.complex(A)) {
        norm(Re(A), "F")^2 + norm(Im(A), "F")^2
      } else norm(A, "F")^2
      rpcaObj$var <- squared_norm / (m-1)
    }
  }
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Compute randomized svd / eigen decomposition
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(rand == TRUE) {
    svdalg = 'rsvd'
  }else { 
    svdalg = 'svd' 
  }
  
  out <- switch(svdalg,
                svd = svd(A, nu = if(retx) k else 0, nv = k),
                rsvd = if(sparse) {
                  rsvd_compute(A, k = k, p = p, q = q,
                               nu = if(retx) k else 0, nv = k,
                               operators = operators)
                } else rsvd(A, k = k, p = p, q = q,
                             nu = if(retx) k else 0, nv = k),
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


#' @rdname rpca
#' @method print rpca
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

#' @rdname rpca
#' @method summary rpca
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

#' @rdname rpca
#' @export
predict.rpca <- function( object, newdata, ...)
{
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Predict
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(inherits(newdata, "sparseMatrix") && !is.complex(object$center) &&
     !is.complex(object$rotation)) {
    operators <- rpca_sparse_products(rpca_sparse_numeric(newdata),
                                      object$center, object$scale,
                                      replace_nan = !is.logical(object$scale))
    x <- operators$forward(as.matrix(object$rotation))
    rownames(x) <- rownames(newdata)
    return(x)
  }

  newdata <- as.matrix(newdata)
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

