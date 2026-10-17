#' @title PCA observation scores
#'
#' @description Display observations in a plane defined by two retained
#' principal components from \code{\link{rpca}}. Each point shows an
#' observation's fitted scores, allowing patterns in the selected components
#' to be compared across observations or supplied groups. Labels and group
#' ellipses can be added. The optional \pkg{ggplot2} package is required.
#'
#' @param rpcaObj An object returned by \code{\link{rpca}}, with
#' real-valued scores. Fit with \code{retx = TRUE} so that \code{rpcaObj$x}
#' is available.
#'
#' @param pcs Integer vector of length two selecting the retained principal
#' components to plot, in horizontal and vertical order. The default is
#' \code{c(1, 2)}; both selected components must be present in the fit.
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
#' when \code{TRUE} (the default). The defaults of
#' \code{\link[ggplot2]{stat_ellipse}} give a multivariate-t ellipse at
#' level \code{0.95}. See Details for its interpretation.
#'
#' @param alpha.ellipse Numeric value between zero and one controlling ellipse
#' transparency. The default is \code{0.2}; ignored when no ellipse is drawn.
#'
#' @param ind_labels Logical; show a text label for each observation when
#' \code{TRUE} (the default). Set to \code{FALSE} to reduce clutter.
#'
#' @param ind_labels.names Optional character vector with one label per row of
#' \code{rpcaObj$x}, in the same order. By default, the score row names are
#' used, or sequential row numbers when names are absent. Labels are displayed
#' only when \code{ind_labels = TRUE}.
#'
#' @details The horizontal and vertical coordinates are the columns of
#' \code{rpcaObj$x} selected by \code{pcs}. These fitted scores retain the
#' centering and scaling used in \code{rpca}; the plotting function does not
#' standardize or recompute them. For randomized fits, they are the stored
#' projected-SVD scores described in \code{\link{rpca}}.
#'
#' Each axis label gives its component's share of \code{rpcaObj$var}, the
#' total variance of the processed input. Variance in omitted components
#' remains in the denominator, so a view of two components need not account
#' for all of the variation.
#'
#' When groups are supplied, point colours identify group membership and
#' data ellipses summarize the spread of fitted scores within each group.
#' The ellipses use the multivariate-t model in
#' \code{\link[ggplot2]{stat_ellipse}}, whose computation also requires
#' the recommended package \pkg{MASS}. They do not estimate confidence
#' regions for group means or test group separation. Small groups or
#' degenerate score distributions may not support ellipse estimation; use
#' \code{ellipse = FALSE} in that case.
#'
#' @return A \pkg{ggplot2} plot object. Print the object to draw the
#' plot, or add layers, labels, and themes with \pkg{ggplot2}.
#'
#' @seealso \code{\link{rpca}}, \code{\link{ggbiplot}},
#' \code{\link{ggcorplot}}, \code{\link{ggscreeplot}},
#' \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{stat_ellipse}}
#'
#' @author N. Benjamin Erichson
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(42)
#'   fit <- rpca(iris[, 1:4], k = 3)
#'   print(ggindplot(fit, pcs = c(1, 3), groups = iris$Species,
#'                   ellipse = FALSE, ind_labels = FALSE))
#' }
#'
#' @export
ggindplot <- function( rpcaObj, pcs = c(1,2), groups = NULL, alpha = 0.6, ellipse = TRUE, alpha.ellipse=0.2,
                       ind_labels=TRUE, ind_labels.names=NULL)
{

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The package 'ggplot2' is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Check if retx is provided
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.null(rpcaObj$x)) stop("ggindplot requires the rotated variables, i.e., set rpca(..., retx = TRUE).")

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
  # Create data frame
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  df <- data.frame(scores=cbind(Z1,Z2), row.names = 1:n)
  colnames(df) <- c( 'a', 'b')
  
  if(is.null(rownames(rpcaObj$x))) {
    df$"indName" <- as.character(1:n)
  } else {
    df$"indName" <- rownames(rpcaObj$x)
  }  
  
  if(!is.null(ind_labels.names)) df$"indName" <- ind_labels.names


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Workaround for CRAN: Nulling
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  x <- NULL # Setting the variables to NULL first
  y <- NULL # Setting the variables to NULL first
  a <- NULL # Setting the variables to NULL first
  b <- NULL # Setting the variables to NULL first
  class <- NULL # Setting the variables to NULL first
  indName <- NULL # Setting the variables to NULL first
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Scores
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.null(groups)) {
        df$"class" <- 'subjects'
  } else{
        df$"class" <- groups
  }

  g <- ggplot2::ggplot(data=df, ggplot2::aes(x = a, y = b, colour = class )) +
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
        g <- g + ggplot2::geom_text(data = df,
                             ggplot2::aes(label = indName, x = a, y = b,
                                          angle = 0, hjust = 0.1, vjust = 0.1),
                             color = 'black', size = 4)
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
