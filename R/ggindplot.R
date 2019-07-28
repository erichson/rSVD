#devtools::use_package('ggplot2')
#library('ggplot2')

#' Individual factor map for \code{\link[rsvd]{rpca}} using \code{\link[ggplot2]{ggplot}}.
#'
#' @description Creates a pretty plot which is showing the individual factor map, i.e,
#'    plotting the principal component scores.
#'    
#' @param  rpcaObj        Object returned by the \code{\link[rsvd]{rpca}} function.
#'                          
#' @param  pcs            Array_like. \cr
#'                        An array with two values indicating the two PCs which should be used for plotting. 
#'                        By default the first two PCs are used, e.g., \eqn{c(1,2)}. 
#'                          
#' @param  groups         Factor, optional. \cr
#'                        Factor indicating groups.
#'                          
#' @param  alpha          Scalar, optional. \cr
#'                        Alpha transparency for scatter plot.
#'                                                                                
#' @param  ellipse        Bool (\eqn{TRUE}, \eqn{FALSE}), optional. \cr
#'                        Draw a 1sd data ellipse for each group, if \eqn{TRUE}.
#'  
#' @param  alpha.ellipse  Scalar, optional. \cr
#'                        Alpha transparency for ellipse.
#' 
#' @param  ind_labels         Bool (\eqn{TRUE}, \eqn{FALSE}), optional. \cr
#'                        Plot names for each individual point, if \eqn{TRUE}.
#' 
#' @param ind_labels.names     Array_like, optional. \cr
#'                        User specific labels for the individual points.     
#' 
#' @seealso \code{\link[rsvd]{rpca}}, \code{\link[ggplot2]{ggplot}}
#'
#' @author N. Benjamin Erichson, \email{erichson@berkeley.edu}
#' 
#' @examples
#' #See ?rpca

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
  g <- g + ggplot2::guides(colour=FALSE)
  g <- g + ggplot2::geom_vline(xintercept=0, linetype="dashed", color = "black")  
  g <- g + ggplot2::geom_hline(yintercept=0, linetype="dashed", color = "black")    
  g <- g + ggplot2::theme_bw() 
  g <- g + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank())    
  
  return(g)
}
