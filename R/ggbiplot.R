#' @title PCA scores and variable directions in one plot
#'
#' @description Overlay the fitted observation scores and variable
#' directions for two retained principal components from \code{\link{rpca}}.
#' The optional \pkg{ggplot2} package is required.
#'
#' @param rpcaObj An object returned by \code{\link{rpca}}, with
#' real-valued scores and directions. Fit with \code{retx = TRUE} so that
#' \code{rpcaObj$x} is available.
#'
#' @param pcs Integer vector of length two selecting the retained principal
#' components to plot, in horizontal and vertical order. The default is
#' \code{c(1, 2)}; both selected components must be present in the fit.
#'
#' @param loadings Logical; when \code{TRUE} (the default), multiply each
#' selected direction by its component standard deviation, so the coordinates
#' are \code{rotation[, pcs]} scaled columnwise by \code{sdev[pcs]}.
#' When \code{FALSE}, plot the unscaled entries of \code{rotation[, pcs]}.
#'
#' @param groups Optional factor with one entry for each row of \code{rpcaObj$x},
#' in the same order. It controls the point colours and, when requested, the
#' ellipses. The default, \code{NULL}, uses one colour and draws no ellipse.
#' If \code{rpca} omitted incomplete observations, supply groups for the
#' retained rows only.
#'
#' @param alpha Numeric value between zero and one controlling point
#' transparency, with zero fully transparent and one opaque. The default is
#' \code{0.6}.
#'
#' @param ellipse Logical; draw a filled data ellipse for each supplied group
#' when \code{TRUE} (the default). Uses the defaults of
#' \code{\link[ggplot2]{stat_ellipse}}: a multivariate-t ellipse at level
#' \code{0.95}. See Details for its interpretation.
#'
#' @param alpha.ellipse Numeric value between zero and one controlling ellipse
#' transparency. The default is \code{0.2}; ignored when no ellipse is drawn.
#'
#' @param var_labels Logical; label the variable arrows when
#' \code{TRUE} (the default).
#'
#' @param var_labels.names Optional character vector with one label per row of
#' \code{rpcaObj$rotation}, in the original variable order. By default, the
#' rotation row names are used, or sequential variable numbers when names
#' are absent.
#' Labels are displayed only when \code{var_labels = TRUE}.
#'
#' @param ind_labels Logical; show a text label for each observation when
#' \code{TRUE} (the default). Set to \code{FALSE} to reduce clutter.
#'
#' @param ind_labels.names Optional character vector with one label per row of
#' \code{rpcaObj$x}, in the same order. By default, the score row names are
#' used, or sequential row numbers when names are absent. Labels are displayed
#' only when \code{ind_labels = TRUE}.
#'
#' @details Observation coordinates are taken directly from
#' \code{rpcaObj$x}. These are the fitted principal-component scores, with the
#' centering and scaling requested in \code{rpca}; no additional normalization
#' is applied. The \code{loadings} argument changes only the arrow coordinates.
#' Thus, arrows and observations use their respective loading and score
#' scales, rather than an automatically matched visual scale.
#'
#' For centered variables with unit sample variance and nonzero component
#' variance, the default arrow coordinates describe correlations of the
#' variables with the fitted component scores. For other preprocessing choices,
#' they are scaled directions and need not be correlations. Use
#' \code{\link{ggcorplot}} to view the directions separately with a unit circle.
#'
#' Axis labels report the selected components and their shares of the total
#' processed variance stored in \code{rpcaObj$var}, including variance outside
#' the retained components.
#'
#' The ellipses describe the spread of fitted scores in
#' each group under the multivariate-t model used by
#' \code{\link[ggplot2]{stat_ellipse}}. Its multivariate-t computation also
#' requires the recommended package \pkg{MASS}. They are data ellipses, not confidence
#' regions for group means or tests of group separation. Small groups or
#' degenerate score distributions may not support ellipse estimation; use
#' \code{ellipse = FALSE} in that case.
#'
#' @return A \pkg{ggplot2} plot object. Print the object to draw the
#' plot, or add layers, labels, and themes with \pkg{ggplot2}.
#'
#' @seealso \code{\link{rpca}}, \code{\link{ggindplot}},
#' \code{\link{ggcorplot}}, \code{\link{ggscreeplot}},
#' \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{stat_ellipse}}
#'
#' @author N. Benjamin Erichson
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(42)
#'   fit <- rpca(iris[, 1:4], k = 2)
#'   p <- ggbiplot(fit, groups = iris$Species, ellipse = FALSE,
#'                 ind_labels = FALSE)
#'   print(p + ggplot2::labs(title = "Iris scores and variable directions"))
#' }
#'
#' @export
ggbiplot <- function( rpcaObj, pcs = c(1,2), loadings=TRUE, groups = NULL, alpha = 0.6, 
                      ellipse = TRUE, alpha.ellipse=0.2,
                      var_labels=TRUE, var_labels.names=NULL,
                      ind_labels=TRUE, ind_labels.names=NULL
                      )
{
  
  if (!requireNamespace('ggplot2', quietly = TRUE)) {
    stop("The package 'ggplot2' is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Check if retx is provided
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.null(rpcaObj$x)) stop("ggbiplot requires the rotated variables, i.e., set rpca(..., retx = TRUE).")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Check selected pcs
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  stopifnot(length(pcs) == 2)
  if(max(pcs) > ncol(rpcaObj$rotation)) stop("Selected PC is not valid.")
  if(min(pcs) < 1) stop("Selected PC is not valid.")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Dimensions
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  n = nrow(rpcaObj$x)
  p = nrow(rpcaObj$rotation)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Label PCs
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  variance = rpcaObj$sdev**2
  explained_variance_ratio = round(variance / rpcaObj$var,3) * 100
  PC1 = paste("PC", pcs[1], "(", explained_variance_ratio[pcs[1]]  , "% explained var.)", sep="")
  PC2 = paste("PC", pcs[2], "(", explained_variance_ratio[pcs[2]]  , "% explained var.)", sep="")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Scale principal component scores
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Z1 = rpcaObj$x[, pcs[1]] #* (rpcaObj$eigvals[pcs[1]]**0.5)
  Z2 = rpcaObj$x[, pcs[2]] #* (rpcaObj$eigvals[pcs[2]]**0.5)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Create data frame for variables scores
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dfscores <- data.frame(scores=cbind(Z1,Z2), row.names = 1:n)
  colnames(dfscores) <- c( 'a', 'b')
  
  if(is.null(rownames(rpcaObj$x))) {
    dfscores$"indName" <- as.character(1:n)
  } else {
    dfscores$"indName" <- rownames(rpcaObj$x)
  }  
  
  if(!is.null(ind_labels.names)) dfscores$"indName" <- ind_labels.names
  

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  # Create data frame for variables
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  if(loadings==FALSE) rotation = rpcaObj$rotation[,pcs]  
  if(loadings==TRUE)  rotation = t(t(rpcaObj$rotation[,pcs]) * rpcaObj$eigvals[pcs]**0.5)
  
  dfvariables <- data.frame(rotation=rotation, row.names = 1:p)
  colnames(dfvariables) <- c( 'a', 'b')
  
  if(is.null(rownames(rpcaObj$rotation))) {
    dfvariables$"varName" <- as.character(1:p)
  } else {
    dfvariables$"varName" <- rownames(rpcaObj$rotation)
  }
  
  if(!is.null(var_labels.names)) dfvariables$"varName" <- var_labels.names
  
    
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Workaround for CRAN: Nulling
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  x <- NULL # Setting the variables to NULL first
  y <- NULL # Setting the variables to NULL first
  a <- NULL # Setting the variables to NULL first
  b <- NULL # Setting the variables to NULL first
  class <- NULL # Setting the variables to NULL first
  indName <- NULL # Setting the variables to NULL first
  varName <- NULL # Setting the variables to NULL first
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Scores
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.null(groups)) {
    dfscores$"class" <- 'subjects'
  } else{
    dfscores$"class" <- groups
  }
  
  g <- ggplot2::ggplot(data=dfscores, ggplot2::aes(x = a, y = b, colour = class )) +
    ggplot2::geom_point(size = 2, alpha = alpha) +
    ggplot2::theme(legend.position = "none")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Stat ellipse
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.null(groups) == 0 && ellipse==TRUE){
    g <- g + ggplot2::stat_ellipse( geom = "polygon", alpha = alpha.ellipse,
                                    ggplot2::aes(fill = class))
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Label data points
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(ind_labels==TRUE) {
    g <- g + ggplot2::geom_text(data = dfscores,
                                ggplot2::aes(label = indName, x = a, y = b,
                                             angle = 0, hjust = 0.1, vjust = 0.1),
                                color = 'black', size = 4)
  }    
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Variable factor map
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  g <- g + ggplot2::geom_segment(data = dfvariables,
                                 ggplot2::aes(x = 0, y = 0, xend = a, yend = b ),
                                 arrow = grid::arrow(length = grid::unit(0.5, 'picas')),
                                 color = '#de2d26' ,  lwd = 0.5)
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Label the variable axes
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(var_labels==TRUE) {
    g <- g + ggplot2::geom_text(data = dfvariables,
                                ggplot2::aes(label = varName, x = a, y = b,
                                             angle = 0, hjust = 0.1, vjust = 0.1),
                                color = '#de2d26', size = 4)
  }
  
  g <- g + ggplot2::ggtitle('Individuals factor map (PCA)')
  g <- g + ggplot2::xlab(PC1) +  ggplot2::ylab(PC2) 
  g <- g + ggplot2::guides(colour="none")
  g <- g + ggplot2::geom_vline(xintercept=0, linetype="dashed", color = "black")  
  g <- g + ggplot2::geom_hline(yintercept=0, linetype="dashed", color = "black")    
  g <- g + ggplot2::theme_bw() 
  g <- g + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                          panel.grid.minor = ggplot2::element_blank())    
  
  return(g)
}
