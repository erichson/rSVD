#devtools::use_package('ggplot2', type = 'Suggests')
#library('ggplot2')

#' Biplot for \code{\link[rsvd]{rpca}} using \code{\link[ggplot2]{ggplot}}.
#'
#' @description Creates a pretty biplot which is showing the individual factor map overlayed by the
#' variables factor map, i.e, plotting both the principal component scores and directions.
#'    
#' @param  rpcaObj        Object returned by the \code{\link[rsvd]{rpca}} function.
#'                          
#' @param  pcs            Array_like. \cr
#'                        An array with two values indicating the two PCs which should be used for plotting. 
#'                        By default the first two PCs are used, e.g., \eqn{c(1,2)}. 
#'
#' @param loadings   Bool (\eqn{TRUE}, \eqn{FALSE}), optional. \cr
#'                  If \eqn{TRUE}, the eigenvectors
#'                  are unit scaled by the square root of the eigenvalues \eqn{W = W * diag(sqrt(eigvals))}.                       
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
#' @param  var_labels  Bool (\eqn{TRUE}, \eqn{FALSE}), optional. \cr
#'                        Plot variable names, if \eqn{TRUE}.
#' 
#' @param var_labels.names     Array_like, optional. \cr
#'                        User specific labels for the individuals.     
#' 
#' @param  ind_labels  Bool (\eqn{TRUE}, \eqn{FALSE}), optional. \cr
#'                        Plot data point names, if \eqn{TRUE}.
#' 
#' @param ind_labels.names     Array_like, optional. \cr
#'                        User specific labels for data points. 
#' 
#' @seealso \code{\link[rsvd]{rpca}}, \code{\link[ggplot2]{ggplot}}
#'
#' @author N. Benjamin Erichson, \email{erichson@berkeley.edu}
#' 
#' @examples
#' #See ?rpca

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
                                 color = '#de2d26' ,  size = 0.5)  
  
  
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
  g <- g + ggplot2::guides(colour=FALSE)
  g <- g + ggplot2::geom_vline(xintercept=0, linetype="dashed", color = "black")  
  g <- g + ggplot2::geom_hline(yintercept=0, linetype="dashed", color = "black")    
  g <- g + ggplot2::theme_bw() 
  g <- g + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                          panel.grid.minor = ggplot2::element_blank())    
  
  return(g)
}