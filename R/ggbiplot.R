#devtools::use_package("ggplot2", type = "Suggests")
#devtools::use_package("plyr", type = "Suggests")
#devtools::use_package("scales", type = "Suggests")
#devtools::use_package("grid", type = "Suggests")
#library(ggplot2)
#library(plyr)
#library(scales)
#library(grid)

#' Biplot for \code{rPCA} using ggplot2
#'
#' @param  rpcaObj          object containing the \code{sdev} component, such as that returned
#'                          by \code{rpca}
#' @param  pcs              array_like \cr
#'                          an array with two values indicating which two PCs should be plotted,
#'                          by default the first two PCs are used, e.g., \eqn{c(1,2)}. \cr
#' @param  groups           optional, factor \cr
#'                          factor indicating groups
#' @param  ellipse          draw a 1sd data ellipse for each group
#' @param  alpha.ellipse    alpha transparency for ellipse
#' @param  loadings         draw arrows for the variables
#' @param  labels           label variables
#' @param  ...              arguments passed to or from other methods, see \code{\link[ggplot2]{ggplot}}.
#'
#' @seealso \code{\link{rpca}}, \code{\link[ggplot2]{ggplot}}
#'
#' @author N. Benjamin Erichson, \email{nbe@st-andrews.ac.uk}
#' @examples
#' #See ?rpca

#' @export
ggbiplot <- function( rpcaObj, pcs = c(1,2), groups = NULL, ellipse = TRUE, alpha.ellipse=0.2,
                      loadings = TRUE, labels=TRUE, ... )
{

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The package 'ggplot2' is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Check selected pcs
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  stopifnot(length(pcs) == 2)
  if(max(pcs) > ncol(rpcaObj$rotation)) stop("Selected PC is not valid.")
  if(min(pcs) < 1) stop("Selected PC is not valid.")

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # dimensions
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
  # No scaling ?
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Z1 = rpcaObj$x[, pcs[1]] #* (rpcaObj$eigvals**0.5)[pcs[1]]
  Z2 = rpcaObj$x[, pcs[2]] #* (rpcaObj$eigvals**0.5)[pcs[2]]

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Scale prinicipal directions i.e., loadings
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dfrotation = rpcaObj$rotation %*%  diag(rpcaObj$eigvals**0.5)
  dfrotation = rpcaObj$rotation %*%  diag(rpcaObj$eigvals**0.5)

  # Generate circle
  #theta <- c(seq(-pi, pi, length = 360))
  #circle <- data.frame(x = cos(theta), y = sin(theta))

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Create data frame
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dfScores <- data.frame(scores=cbind(Z1,Z2))
  colnames(dfScores) <- c( 'a', 'b')

  dfrotation <- data.frame(rotation=dfrotation)
  rownames(dfrotation) <- rownames(rpcaObj$rotation)
  colnames(dfrotation) <- c( 'a', 'b')
  dfrotation$"varName" <- rownames(rpcaObj$rotation)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Workaround for CRAN: Nulling
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  x <- NULL # Setting the variables to NULL first
  y <- NULL # Setting the variables to NULL first
  a <- NULL # Setting the variables to NULL first
  b <- NULL # Setting the variables to NULL first
  Class <- NULL # Setting the variables to NULL first
  varName <- NULL # Setting the variables to NULL first
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Create ggplot2:: ggplot
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Circle
  #g <- ggplot2::ggplot( circle , ggplot2::aes( x , y)  )  +
  #     ggplot2::geom_path( size=1, colour="#9400d3"  )

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Scores
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(is.null(groups)) {
        dfScores$"Class" <- 'Observations'
      } else{
        dfScores$"Class" <- groups

      }

      g <- ggplot2::ggplot(data=dfScores, ggplot2::aes(x = a, y = b, colour = Class )) +
           ggplot2::geom_point(size = 2) +
           ggplot2::theme(legend.position = "none")

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # stat_ellipse
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(is.null(groups) == 0 && ellipse==TRUE){
      g <- g + ggplot2::stat_ellipse( geom = "polygon", alpha = alpha.ellipse,
                               ggplot2::aes(fill = Class))
      }

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Principal component directions
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(loadings==TRUE){
        g <- g +
          ggplot2::geom_segment(data = dfrotation,
                                ggplot2::aes(x = 0, y = 0, xend = a, yend = b),
                                arrow = grid::arrow(length = grid::unit(1, 'picas')),
                                color = 'black',  size = 1.3)
      }

      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Label the variable axes
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(labels==TRUE) {
        g <- g +
          ggplot2::geom_text(data = dfrotation,
                             ggplot2::aes(label = varName, x = a, y = b,
                                          angle = 1, hjust = 1, vjust = 1),
                             color = 'black', size = 5)
      }
  # g <- g + guides(colour=FALSE)
  g <- g + ggplot2::labs(x = noquote(PC1), y = noquote(PC2)) + ggplot2::guides(colour=FALSE)
  return(g)
}
