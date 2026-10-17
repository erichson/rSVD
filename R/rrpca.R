#' @title Randomized robust principal component analysis
#'
#' @description Decompose a real matrix into a low-rank component and an entrywise
#' sparse component using principal component pursuit, with randomized or
#' deterministic singular value decompositions.
#'
#' @details Principal component pursuit separates a matrix into components satisfying
#' \deqn{A = L + S.}{A = L + S.}
#' The objective is to minimize the nuclear norm of \eqn{L} plus
#' \code{lambda} times the entrywise one-norm of \eqn{S}, subject to this
#' equality. The nuclear norm
#' is the sum of singular values; the entrywise one-norm is the sum of
#' absolute entries. These penalties encourage a low-rank component and a
#' component with few nonzero entries. Whether the result recovers an
#' underlying signal depends on the data and the separation assumptions;
#' the decomposition is not a general guarantee of denoising.
#'
#' The inexact augmented Lagrange multiplier method (IALM) alternates
#' entrywise soft thresholding and singular-value thresholding. The
#' low-rank component's rank is selected through thresholding rather than
#' through a fixed user-specified rank. Input and working arrays are dense,
#' including for sparse input. No automatic centering or scaling is applied:
#' any desired preprocessing must be done before calling \code{rrpca}.
#' Use \code{set.seed} before fitting for reproducible random draws when
#' \code{rand = TRUE}.
#'
#' The stopping test requires both the relative reconstruction residual
#' \deqn{\frac{\|A-L-S\|_F}{\|A\|_F}}{norm(A - L - S, "F") / norm(A, "F")}
#' and a stationarity diagnostic to be at most \code{tol}, unless
#' \code{maxiter} is reached first. With the sparse component updated first,
#' the stationarity diagnostic is
#' \deqn{\frac{\mu\|L_{new}-L_{old}\|_F}{\max(1,\|Z\|_F)}}{mu * norm(L_new - L_old, "F") / max(1, norm(Z, "F"))}
#' where \eqn{Z} is the updated Lagrange multiplier and \eqn{\mu} is the
#' penalty used for that iteration. The penalty increases only when this
#' diagnostic is sufficiently small. This avoids stopping on reconstruction
#' feasibility alone and makes the tests invariant to rescaling the input.
#'
#' Randomized singular-value thresholding increases its target rank if all
#' returned singular values exceed the shrinkage threshold. If accuracy
#' stalls, the algorithm makes bounded increases in the number of power
#' iterations before falling back to deterministic SVD. These refinements
#' can increase runtime. The stopping diagnostics and randomized
#' approximations are not a general certificate of optimality. The iteration
#' limit may be reached before both diagnostics satisfy \code{tol}.
#' \code{trace = TRUE} displays both diagnostics; \code{err} records only
#' reconstruction residuals, so a small final \code{err} alone is not a
#' convergence certificate.
#'
#' \code{rrpca} is a batch decomposition of the supplied matrix. Its
#' \code{L} and \code{S} contain fitted entries for those observations, not
#' a reusable PCA rotation matrix, and there is no \code{predict.rrpca}
#' method. Projecting new rows onto a basis extracted from \code{L} gives
#' an ordinary low-rank projection with a generally dense residual; it does
#' not compute a robust sparse component for the new rows. Robust inference
#' with a fixed learned basis requires a separate optimization problem.
#' For ordinary PCA with new-data scores, use \code{\link{rpca}} and its
#' \code{predict} method.
#'
#' @param A A real numeric matrix or data frame to decompose. A sparse matrix from
#' \pkg{Matrix} is accepted and converted to a dense matrix. Rows containing
#' missing values are omitted with a warning. In the dimensions below,
#' \eqn{m} and \eqn{n} are the numbers of rows and columns after omission;
#' both must be positive.
#'
#' @param lambda A positive numeric penalty on the entrywise sparse component. The
#' default, \code{NULL}, uses \code{1 / sqrt(max(m, n))}. Larger values
#' penalize the sparse component more strongly.
#'
#' @param maxiter A positive integer maximum number of iterations (default \code{50}).
#' Reaching this limit does not imply convergence.
#'
#' @param tol A positive numeric tolerance for both the relative reconstruction
#' residual and stationarity diagnostic (default \code{1e-5}).
#'
#' @param p A nonnegative integer oversampling parameter for \code{\link{rsvd}}
#' (default \code{10}); used only in randomized SVD steps.
#'
#' @param q A nonnegative integer initial number of power iterations for
#' \code{\link{rsvd}} (default \code{2}); increased adaptively if randomized
#' accuracy stalls. Ignored when \code{rand = FALSE}.
#'
#' @param trace Logical; print iteration progress, including both convergence
#' diagnostics (default \code{FALSE}).
#'
#' @param rand Logical; start with randomized SVD steps when \code{TRUE} (default),
#' with deterministic SVD used for large target ranks or stalled accuracy.
#' With \code{FALSE}, use the deterministic SVD throughout.
#'
#' @return A list of class \code{"rrpca"} containing:
#' \describe{
#'   \item{L}{A dense \eqn{m} by \eqn{n} matrix: the estimated low-rank
#'     component.}
#'   \item{S}{A dense \eqn{m} by \eqn{n} matrix containing the estimated
#'     sparse component. Sparsity refers to its entries, not its storage
#'     class.}
#'   \item{err}{A numeric vector of relative Frobenius reconstruction
#'     residuals, one per completed iteration. For a zero input, the
#'     function returns zero components and \code{err = 0} immediately.}
#' }
#' Rows of \code{L} and \code{S} correspond to retained input rows after
#' missing rows have been omitted. No fitted loading basis or new-data
#' prediction method is returned.
#'
#' @references \itemize{
#'   \item N. B. Erichson, S. Voronin, S. L. Brunton and J. N. Kutz. 2019.
#'   Randomized Matrix Decompositions Using R.
#'   Journal of Statistical Software, 89(11), 1--48.
#'   \doi{10.18637/jss.v089.i11}.
#'   \item Z. Lin, M. Chen and Y. Ma. 2010.
#'   The Augmented Lagrange Multiplier Method for Exact Recovery of
#'   Corrupted Low-Rank Matrices.
#'   \url{https://arxiv.org/abs/1009.5055}.
#'   \item S. Boyd, N. Parikh, E. Chu, B. Peleato and J. Eckstein. 2011.
#'   Distributed Optimization and Statistical Learning via the Alternating
#'   Direction Method of Multipliers.
#'   Foundations and Trends in Machine Learning, 3(1), 1--122.
#'   \doi{10.1561/2200000016}.
#' }
#'
#' @seealso \code{\link{rpca}}, \code{\link{rsvd}}.
#'
#' @author N. Benjamin Erichson
#'
#' @examples # A small video: each column is one frame, with pixels stored in rows.
#' background <- outer(1:32, 1:32, function(x, y) 0.1 + 0.001 * (x + y))
#' video <- matrix(rep(as.vector(background), 30), nrow = 32 * 32)
#' for (i in 1:30) {
#'   foreground <- matrix(0, 32, 32)
#'   foreground[i:(i + 2), 15:17] <- 0.8
#'   video[, i] <- video[, i] + as.vector(foreground)
#' }
#'
#' set.seed(42)
#' out <- rrpca(video, maxiter = 100, tol = 1e-6)
#' tail(out$err, 1) # Reconstruction residual; see Details for convergence.
#'
#' # Compare one frame and its estimated background and foreground.
#' # Each panel uses its own intensity range.
#' old_par <- par(mfrow = c(1, 4), mar = c(0, 0, 2, 0))
#' image(background, col = gray.colors(64), axes = FALSE,
#'       asp = 1, useRaster = TRUE, main = "Background")
#' image(matrix(video[, 10], 32, 32), col = gray.colors(64),
#'       axes = FALSE, asp = 1, useRaster = TRUE, main = "Frame")
#' image(matrix(out$L[, 10], 32, 32), col = gray.colors(64),
#'       axes = FALSE, asp = 1, useRaster = TRUE, main = "Low-rank part")
#' image(matrix(out$S[, 10], 32, 32), col = gray.colors(64),
#'       axes = FALSE, asp = 1, useRaster = TRUE, main = "Sparse part")
#' par(old_par)
#'
#' @export
rrpca <- function(A, lambda=NULL, maxiter=50, tol=1.0e-5, p=10, q=2, trace=FALSE, rand=TRUE) UseMethod("rrpca")

#' @export
rrpca.default <- function(A, lambda=NULL, maxiter=50, tol=1.0e-5, p=10, q=2, trace=FALSE, rand=TRUE) {
  #*************************************************************************
  #***        Author: N. Benjamin Erichson <nbe@st-andrews.ac.uk>        ***
  #***                              <2016>                               ***
  #***                       License: BSD 3 clause                       ***
  #*************************************************************************
  A <- as.matrix(A)
  if(any(is.na(A))) {
    warning("Missing values are omitted: na.omit(A).")
    A <- stats::na.omit(A)
  }
  m <- nrow(A)
  n <- ncol(A)
  if(m < 1L || n < 1L) {
    stop("Input matrix must have positive dimensions after omitting missing values.")
  }

  rrpcaObj = list(L = NULL,
                  S = NULL,
                  err = NULL)

  froNorm <- norm(A, "F")
  if(froNorm == 0) {
    rrpcaObj$L <- matrix(0, nrow = m, ncol = n)
    rrpcaObj$S <- matrix(0, nrow = m, ncol = n)
    rrpcaObj$err <- 0
    class(rrpcaObj) <- "rrpca"
    return(rrpcaObj)
  }

  # Set target rank
  k <- 1L
  max_rank <- min(m, n)

  # Set lambda, gamma, rho
  if(is.null(lambda)) lambda <- max(m,n)**-0.5
  gamma <- 1.25
  rho <- 1.5

  if(rand == TRUE) {
    svdalg = 'rsvd'
  }else { 
    svdalg = 'svd' 
  }  
  
  # Compute matrix norms
  spectralNorm <- switch(svdalg,
                    svd = norm(A, "2"),
                    rsvd = rsvd(A, k=1, p=10, q=1, nu=0, nv=0)$d,
                    stop("Selected SVD algorithm is not supported!")
  )

  # The dual of the entrywise l1 norm is the largest absolute entry,
  # rather than the matrix infinity norm (maximum absolute row sum).
  infNorm <- norm(A, "M") / lambda
  dualNorm <- max( spectralNorm , infNorm)

  # Initialize Lagrange multiplier
  Z <- A / dualNorm

  # Initialize tuning parameter
  mu <- gamma / spectralNorm
  mubar <- mu * 1e7
  mu <- min( mu * rho , mubar )

  # Init low-rank and sparse matrix
  L = matrix(0, nrow = m, ncol = n)
  S = matrix(0, nrow = m, ncol = n)  
  
  niter <- 1
  err <- 1
  dualErr <- 1
  power_iterations <- q
  refinement <- 0L
  stagnation <- 0L
  best_error <- Inf
  use_exact <- FALSE
  while((err > tol || dualErr > tol) && niter <= maxiter) {

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Update S using soft-threshold
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      epsi = lambda / mu
      temp_S = A - L + Z / mu
      
      S = matrix(0, nrow = m, ncol = n)

      idxL <- which(temp_S < -epsi)
      idxH <- which(temp_S > epsi)
      S[idxL] <- temp_S[idxL] + epsi
      S[idxH] <- temp_S[idxH] - epsi
      remove(temp_S, idxL, idxH)
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #Singular Value Decomposition
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      R <- A - S + Z / mu

      randomized_step <- FALSE
      if(svdalg == 'svd') svd_out <- svd(R)
      if(svdalg == 'rsvd') {
        if(use_exact || k > max_rank / 5) {
          svd_out <- svd(R)
        } else {
          target <- min(k + 10L, max_rank)
          repeat {
            svd_out <- rsvd(R, k = target, p = p, q = power_iterations)
            randomized_step <- TRUE
            # A saturated truncation can omit singular values that survive
            # shrinkage. Increase its rank before accepting this update.
            # This is a coverage safeguard, not an error certificate for
            # the randomized singular values themselves.
            if(target == max_rank || svd_out$d[target] <= 1 / mu) break
            target <- min(2L * target, max_rank)
            if(target > max_rank / 5) {
              svd_out <- svd(R)
              randomized_step <- FALSE
              break
            }
          }
        }
      }
      remove(R)
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Predict optimal rank and update
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      svp = sum(svd_out$d > 1/mu)
      
      if(svp <= k){
        k = min(svp + 1, max_rank)
      } else {
        k = min(svp + round(0.05 * max_rank), max_rank)
      }

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Truncate SVD and update L
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      previous_L <- L
      if(svp > 0L) {
        keep <- seq_len(svp)
        L <- sweep(svd_out$u[, keep, drop=FALSE], 2,
                   svd_out$d[keep] - 1/mu, FUN = "*") %*%
          t(svd_out$v[, keep, drop=FALSE])
      } else {
        L[] <- 0
      }

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Compute error
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Astar = A - L - S
      Z = Z + Astar * mu

      err = norm( Astar , 'F') / froNorm
      # With exact SVD and the S-then-L order, the two proximal subgradients
      # differ by mu * (L - previous_L). Feasibility alone can be exact at
      # a nonoptimal decomposition, so also require this change to be small.
      # Normalize by the dual variable to preserve input-scale invariance.
      dualErr <- mu * norm(L - previous_L, "F") / max(1, norm(Z, "F"))
      remove(previous_L)

      # A fixed sketch accuracy can stall the stationarity test. Increase
      # the power iterations after three steps without a 10% improvement.
      # Try two refinements before using exact SVD for the remaining steps.
      if(randomized_step) {
        current_error <- max(err, dualErr)
        if(current_error < 0.9 * best_error) {
          best_error <- current_error
          stagnation <- 0L
        } else {
          stagnation <- stagnation + 1L
        }
        if(stagnation >= 3L) {
          if(refinement < 2L) {
            power_iterations <- power_iterations + 2L
            refinement <- refinement + 1L
          } else {
            use_exact <- TRUE
          }
          best_error <- Inf
          stagnation <- 0L
        }
      }
      rrpcaObj$err <- c(rrpcaObj$err, err)

      if(trace==TRUE){
        cat('\n', paste0('Iteration: ', niter ), paste0(' predicted rank = ', svp ), paste0(' target rank k = ', k ),  paste0(' Fro. error = ', rrpcaObj$err[niter] ), paste0(' stationarity error = ', dualErr))
        }

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Update mu
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Increasing mu before the proximal steps agree can make a feasible
      # but incorrect split stagnate. Keep the penalty fixed until they do.
      if(dualErr <= tol) mu <- min(mu * rho, mubar)

      niter =  niter + 1

  }# End while loop
  rrpcaObj$L <- L
  rrpcaObj$S <- S
  
  class(rrpcaObj) <- "rrpca"
  return( rrpcaObj )

}



