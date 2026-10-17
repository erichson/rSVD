#' @title PCA variable directions and correlation circle
#'
#' @description Plot the variable directions for two retained
#' principal components, optionally scaled by the component standard
#' deviations, together with a unit circle. The optional \pkg{ggplot2} package
#' is required.
#'
#' @param rpcaObj An object returned by \code{\link{rpca}}, with
#' real-valued directions. Scores are not needed, so fits made with
#' \code{retx = FALSE} are supported.
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
#' @param var_labels Logical; show text labels for the variables
#' selected by \code{top.n} when \code{TRUE}. The default is \code{FALSE}.
#'
#' @param var_labels.names Optional character vector with one label per row of
#' \code{rpcaObj$rotation}, in the original variable order. By default, the
#' rotation row names are used, or sequential variable numbers when names
#' are absent.
#' These names identify the endpoint markers and, when requested, the text
#' labels.
#'
#' @param alpha Numeric value between zero and one controlling
#' arrow transparency. The default is \code{1} (opaque).
#'
#' @param top.n Positive integer limiting the number of variables
#' selected for endpoint markers and text labels. Selection is by decreasing
#' arrow length in the plotted two-component plane. The default, \code{NULL},
#' selects all variables; values larger than the number of variables are
#' clamped to that number. See Details for the marker display rule.
#'
#' @details With \code{loadings = TRUE}, the plotted coordinate of
#' variable \eqn{j} on component \eqn{i} is
#' \eqn{V_{ji}\sqrt{\lambda_i}}{V[j, i] * sqrt(lambda[i])}, where \eqn{V} is
#' \code{rotation} and \eqn{\lambda_i}{lambda[i]} is the corresponding
#' eigenvalue. For centered variables
#' with unit sample variance and nonzero component variance, this is the
#' correlation between the variable and the fitted component score. In that
#' setting the unit circle helps show how well the selected component plane
#' represents each variable. With unscaled or uncentered input, the coordinates
#' need not be correlations and the circle is only a visual reference.
#'
#' With \code{loadings = FALSE}, the plotted coordinates are the unscaled
#' principal-component directions. This changes the relative arrow lengths
#' and may change which variables \code{top.n} selects.
#'
#' All variable arrows are drawn. The \code{top.n} selection affects endpoint
#' markers and, if \code{var_labels = TRUE}, text labels. Endpoint markers are
#' drawn only when the effective number selected is less than 50. A variable's
#' ranking measures its arrow length in this selected plane, not its overall
#' importance across all components.
#'
#' The plot uses equal horizontal and vertical scales. Axis labels include
#' the variance shares relative to \code{rpcaObj$var}, the total variance of
#' the processed input, rather than only the retained variance.
#'
#' @return A \pkg{ggplot2} plot object. Print the object to draw the
#' plot, or add layers, labels, and themes with \pkg{ggplot2}.
#'
#' @seealso \code{\link{rpca}}, \code{\link{ggbiplot}},
#' \code{\link{ggindplot}}, \code{\link{ggscreeplot}},
#' \code{\link[ggplot2]{ggplot}}
#'
#' @author N. Benjamin Erichson
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(42)
#'   fit <- rpca(iris[, 1:4], k = 2, retx = FALSE)
#'   # Label the three longest arrows in the selected component plane.
#'   print(ggcorplot(fit, var_labels = TRUE, top.n = 3))
#' }
#'
#' @export
ggcorplot <- function(rpcaObj, pcs=c(1,2),  loadings=TRUE, var_labels=FALSE, var_labels.names=NULL, alpha=1, top.n=NULL) {
  
  if (!requireNamespace('ggplot2', quietly = TRUE)) {
    stop("The package 'ggplot2' is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Number of variables
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  p = nrow(rpcaObj$rotation)  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Check selected pcs
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  stopifnot(length(pcs) == 2)
  if(max(pcs) > ncol(rpcaObj$rotation)) stop("Selected PC is not valid.")
  if(min(pcs) < 1) stop("Selected PC is not valid.")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  # Select PCs
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  PC1 = paste("PC", pcs[1], sep="")
  PC2 = paste("PC", pcs[2], sep="")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  # Generate circle
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  theta <- c(seq(-pi, pi, length = 360))
  circle <- data.frame(x = cos(theta), y = sin(theta))
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  # Create data frame
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  if(loadings==FALSE) rotation = rpcaObj$rotation[,pcs]  
  if(loadings==TRUE)  rotation = t(t(rpcaObj$rotation[,pcs]) * rpcaObj$eigvals[pcs]**0.5)
  
  df <- data.frame(rotation=rotation, row.names = 1:p)
  colnames(df) <- c( 'a', 'b')

    if(is.null(rownames(rpcaObj$rotation))) {
    df$"varName" <- as.character(1:p)
  } else {
    df$"varName" <- rownames(rpcaObj$rotation)
  }
  
  if(!is.null(var_labels.names)) df$"varName" <- var_labels.names
  df$abs <- sqrt(df$a**2 + df$b**2)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Label PCs
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  variance = rpcaObj$sdev**2
  explained_variance_ratio = round(variance / rpcaObj$var, 3) * 100
  PC1 = paste("PC ", pcs[1], "(", explained_variance_ratio[pcs[1]]  , "% explained var.)", sep="")
  PC2 = paste("PC ", pcs[2], "(", explained_variance_ratio[pcs[2]]  , "% explained var.)", sep="")
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  #Workaround for CRAN: Nulling
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  x <- NULL # Setting the variables to NULL first
  y <- NULL # Setting the variables to NULL first
  a <- NULL # Setting the variables to NULL first
  b <- NULL # Setting the variables to NULL first
  varName <- NULL # Setting the variables to NULL first  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  # Create ggplot2::ggplot
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  g <- ggplot2::ggplot( circle , ggplot2::aes( x , y)  ) +
    ggplot2::geom_path( lwd=0.5, colour="black" )

  
  if(is.null(top.n)) top.n <- nrow(df)
  if(top.n>nrow(df)) top.n <- nrow(df)
  if(top.n < 50) {
  g <- g + ggplot2::geom_point(data = df[order(df$abs, decreasing=TRUE)[1:top.n], ], 
                               size = 4, mapping = ggplot2::aes(x = a, y = b, colour = varName ) ) +
    ggplot2::theme(legend.position = "none")
  }
  
  
  g <- g + ggplot2::geom_segment(data = df,
                                 ggplot2::aes(x = 0, y = 0, xend = a, yend = b ),
                                 arrow = grid::arrow(length = grid::unit(0.5, 'picas')),
                                 color = 'black' ,  lwd = 0.5, alpha = alpha)
  
  g <- g + ggplot2::coord_fixed(ratio=1) 
  g <- g + ggplot2::ggtitle('Variables factor map (PCA)')
  g <- g + ggplot2::xlab(PC1) +  ggplot2::ylab(PC2) 
  g <- g + ggplot2::guides(colour=ggplot2::guide_legend(title=NULL))
  
  g <- g + ggplot2::geom_vline(xintercept=0, linetype="dashed", color = "black")  
  g <- g + ggplot2::geom_hline(yintercept=0, linetype="dashed", color = "black")  
  
  # Label the variable axes
  if(var_labels == TRUE) {
    #df$a <- df$a  *1.1
    #df$b <- df$b  * 1.1
    g <- g + ggplot2::geom_text(data = df[order(df$abs, decreasing=TRUE)[1:top.n], ], ggplot2::aes(label = varName, x = a, y = b,
                                                        angle = 0, hjust = 0, vjust = 0),
                                color = 'black', size = 4)
    
  }  
  
  
  
  g <- g + ggplot2::theme_bw() #+ ggplot2::xlim(c(-1.2,1.2)) +  ggplot2::ylim(c(-1.2,1.2)) +
  g <- g + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank())    
   
  return( g )
}
