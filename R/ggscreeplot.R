#devtools::use_package('ggplot2')
#library('ggplot2')

#' @title Pretty Screeplot
#
#' @description Creates a pretty screeplpot using \code{\link[ggplot2]{ggplot}}. By default the explained variance is plotted
#' agaings the number of the principal component.
#' Alternatively the explained variance ratio, the cumulative
#' explained variance ratio, or the eigenvalues can be plotted.
#'
#' @param rpcaObj  Object returned by the \code{\link[rsvd]{rpca}} function.
#'
#' @param type   String c('var', 'ratio', 'cum', 'eigenvals'), optional. \cr
#'
#' @seealso \code{\link[rsvd]{rpca}}, \code{\link[ggplot2]{ggplot}}
#'
#' @author N. Benjamin Erichson, \email{erichson@berkeley.edu}
#'
#' @examples #
#'

#'@export
ggscreeplot <- function(rpcaObj, type = c('var', 'ratio', 'cum', 'eigenvals')) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The package 'ggplot2' is needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  
  type <- match.arg(type)
  
  y <- switch(type,
              var = rpcaObj$sdev**2,
              ratio = rpcaObj$sdev**2 / rpcaObj$var,
              cum = cumsum(rpcaObj$sdev**2 / rpcaObj$var ),
              eigenvals = rpcaObj$eigvals,
              stop("Selected plot option is not supported!")
  )
  
  y.label <- switch(type,
                    var = 'Explained variance',
                    ratio = 'Proportion of variance',
                    cum = 'Cummulative proportion',
                    eigenvals = 'Eigenvalues',
                    stop("Selected plot option is not supported!")
  )
  
  df <- data.frame(x = 1:length(rpcaObj$sdev), y = y)
  
  #Workaround for CRAN: Nulling
  x <- NULL # Setting the variables to NULL first
  
  g <- ggplot2::ggplot(data = df, ggplot2::aes(x = x, y = y) ) +
    ggplot2::xlab('Principal components') + ggplot2::ylab( y.label ) +
    ggplot2::geom_line(size=0.5, color = 'black', linetype="dashed") +
    ggplot2::geom_point(size=2, color = '#ef3b2c') + 
    ggplot2::guides(colour=FALSE)
  
  g <- g + ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank())
  
  
  return(g)
}
