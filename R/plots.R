#devtools::use_package("ggplot2")
#library(ggplot2)

#' @title Screeplot
#
#' @description Creates a screeplpot. By default the explained variance is plotted
#' agaings the number of the principal component.
#' Alternatively the eigenvalues, explained variance ratio or the cumulative
#' explained variance ratio can be plotted.
#'
#' @param x   object containing the \code{sdev} component, such as that returned
#'            by \code{rpca}
#'
#' @param type      str c('var', 'ratio', 'cum', 'eigenvals'), optional \cr
#'
#' @param ...     arguments passed to or from other methods, see \code{\link{plot}}.
#'
#' @param ................. .
#'
#' @seealso \code{\link{rpca}}
#'
#' @examples #
#'

#'@export
plot.rpca <- function(x, type = c('var', 'ratio', 'cum', 'eigenvals'), ... ) {

  type <- match.arg(type)

  y <- switch(type,
              var = x$sdev**2,
              ratio = x$sdev**2 / x$var,
              cum = cumsum(x$sdev**2 / x$var ),
              eigenvals = x$eigvals,
              stop("Selected plot option is not supported!")

  )

  y.label <- switch(type,
                    var = 'Explained variance',
                    ratio = 'Proportion of variance',
                    cum = 'Cummulative proportion',
                    eigenvals = 'Eigenvalues',
                    stop("Selected plot option is not supported!")
  )

  df <- data.frame(PC = 1:length(x$sdev), y = y)

  graphics::plot.default(df$PC, df$y, xlab ='PCs', ylab=y.label,
               type = 'b', pch=20, col='red')
}



#' @title Pretty Screeplot
#
#' @description Creates a pretty screeplpot using \code{\link[ggplot2]{ggplot}}. By default the explained variance is plotted
#' agaings the number of the principal component.
#' Alternatively the eigenvalues, explained variance ratio or the cumulative
#' explained variance ratio can be plotted.
#'
#' @param rpcaObj  object containing the \code{sdev} component, such as that returned
#'                 by \code{rpca}
#'
#' @param type   str c('var', 'ratio', 'cum', 'eigenvals'), optional \cr
#'
#' @param ................. .
#'
#' @seealso \code{\link{rpca}}, \code{\link[ggplot2]{ggplot}}
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

  df <- data.frame(PC = 1:length(rpcaObj$sdev), y = y)

  #Workaround for CRAN: Nulling
  PC <- NULL # Setting the variables to NULL first


  ggplot2::ggplot(data = df, ggplot2::aes(x = PC, y = y, color = 'red') ) +
    ggplot2::xlab('PCs') + ggplot2::ylab( y.label ) +
    ggplot2::geom_point(size=5) + ggplot2::geom_line(size=1.2) + ggplot2::guides(colour=FALSE)
}



#' @title Correlation plot
#
#' @description Creates a pretty plot which is showing the correlation of
#'    the original variable with the principal component (PCs).
#'
#' @param rpcaObj    object containing the \code{sdev} component, such as that returned
#'                   by \code{rpca}
#'
#' @param pcs   array_like \cr
#'              an array with two values indicating which two PCs should be plotted,
#'              by default the first two PCs are used, e.g., \eqn{c(1,2)}. \cr
#'
#' @param ................. .
#'
#' @seealso \code{\link{rpca}}, \code{\link[ggplot2]{ggplot}}
#'
#' @examples #
#'


#' @export
ggcorplot <- function( rpcaObj, pcs=c(1,2) ) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The package 'ggplot2' is needed for this function to work. Please install it.",
         call. = FALSE)
  }


  # Check selected pcs
  stopifnot(length(pcs) == 2)
  if(max(pcs) > ncol(rpcaObj$rotation)) stop("Selected PC is not valid.")
  if(min(pcs) < 1) stop("Selected PC is not valid.")

  # Select PCs
  PC1 = paste("PC", pcs[1], sep="")
  PC2 = paste("PC", pcs[2], sep="")

  # Generate circle
  theta <- c(seq(-pi, pi, length = 360))
  circle <- data.frame(x = cos(theta), y = sin(theta))

  # Create data frame
  df <- data.frame(rpcaObj$rotation[ , pcs], labels=row.names(rpcaObj$rotation))
  colnames(df) <- c( 'a', 'b', 'labels' )

  #Workaround for CRAN: Nulling
  x <- NULL # Setting the variables to NULL first
  y <- NULL # Setting the variables to NULL first
  a <- NULL # Setting the variables to NULL first
  b <- NULL # Setting the variables to NULL first

  # Create ggplot2:: ggplot
  g <- ggplot2::ggplot( circle , ggplot2::aes( x , y)  )  + ggplot2::geom_path( size=1, colour="#9400d3"  )
  g <- g + ggplot2::geom_point(data=df, size = 6,  mapping = ggplot2::aes(x = a, y = b, label = labels, colour = labels ) ) +
    ggplot2::coord_fixed(ratio=1) + ggplot2::labs(x = noquote(PC1), y = noquote(PC2)) + ggplot2::guides(colour=ggplot2::guide_legend(title=NULL))

  return( g )
}
