#devtools::use_package('ggplot2')
#library('ggplot2')

#' @title  Variables factor map for \code{\link[rsvd]{rpca}} using \code{\link[ggplot2]{ggplot}}.
#
#' @description Creates a pretty plot which is showing the correlation of
#'    the original variable with the principal component (PCs).
#'
#' @param rpcaObj  Object returned by the \code{\link[rsvd]{rpca}} function.
#'
#' @param pcs   Array_like. \cr
#'              An array with two values indicating the two PCs which should be used for plotting. 
#'              By default the first two PCs are used, e.g., \eqn{c(1,2)}. 
#'              
#' @param loadings   Bool (\eqn{TRUE}, \eqn{FALSE}), optional. \cr
#'                  If \eqn{TRUE}, the eigenvectors
#'                  are unit scaled by the square root of the eigenvalues \eqn{W = W * diag(sqrt(eigvals))}.
#'              
#' @param var_labels  Bool (\eqn{TRUE}, \eqn{FALSE}), optional. \cr
#'                Plot variable names, if \eqn{TRUE}.
#'                
#' @param var_labels.names  Array_like, optional. \cr
#'                User specific labels for the variables                  
#' 
#' @param  alpha  Scalar, optional. \cr
#'                Alpha transparency of the arrows. 
#' 
#' @param  top.n  Scalar, optional. \cr
#'                Number of (most influencial) variables to label with small circles.  
#' 
#' @seealso \code{\link[rsvd]{rpca}}, \code{\link[ggplot2]{ggplot}}
#'
#' @author N. Benjamin Erichson, \email{erichson@berkeley.edu}
#'
#' @examples #
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
    ggplot2::geom_path( size=0.5, colour="black" )  

  
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
                                 color = 'black' ,  size = 0.5, alpha = alpha)
  
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
