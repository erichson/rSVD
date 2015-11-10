devtools::use_package("ggplot2")


#' @title Screeplot
#
#' @description Creates a screeplpot. By default the explained variance is plotted
#' agaings the number of the principal component.
#' Alternatively the eigenvalues, explained variance ratio or the cumulative
#' explained variance ratio can be plotted.
#'
#' @param x      object containing the \code{sdev} component, such as that returned
#'               by \code{rpca}
#'
#' @param type   str c('var', 'ratio', 'cum', 'eigenvals'), optional \cr
#'
#' @param ................. .
#'
#' @seealso \code{\link{rsvd}}
#'
#' @examples #
#'
#'

plot.rpca <- function(rpcaObj, type = c('var', 'ratio', 'cum', 'eigenvals')) {

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
                    ratio = 'Relative explained variance',
                    cum = 'Cummulative proportion of explained variance',
                    eigenvals = 'Eigenvalues',
                    stop("Selected plot option is not supported!")
  )

  df <- data.frame(PC = 1:length(rpcaObj$sdev), y = y)

  plot.default(df$PC, df$y, xlab ='Principal components number', ylab=y.label,
               type = 'b', pch=20, col='red')
}



#' @title Pretty Screeplot
#
#' @description Creates a pretty screeplpot using \code{ggplot2}. By default the explained variance is plotted
#' agaings the number of the principal component.
#' Alternatively the eigenvalues, explained variance ratio or the cumulative
#' explained variance ratio can be plotted.
#'
#' @param x      object containing the \code{sdev} component, such as that returned
#'               by \code{rpca}
#'
#' @param type   str c('var', 'ratio', 'cum', 'eigenvals'), optional \cr
#'
#' @param ................. .
#'
#' @seealso \code{\link{rsvd}}, \code{ggplot2}
#'
#' @examples #
#'
#'

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
                    ratio = 'Relative explained variance',
                    cum = 'Cummulative proportion of explained variance',
                    eigenvals = 'Eigenvalues',
                    stop("Selected plot option is not supported!")
                    )

  df <- data.frame(PC = 1:length(rpcaObj$sdev), y = y)

  ggplot(data = df, aes(x = PC, y = y, color = 'red') ) +
    xlab('Principal components number') + ylab( y.label ) +
    geom_point(size=5) + geom_line(size=1.2) + guides(colour=FALSE)
}



#' @title Correlation plot
#
#' @description Creates a pretty plot which is showing the correlation of
#'    the original variable with the principal component (PCs).
#'
#' @param x      object containing the \code{sdev} component, such as that returned
#'               by \code{rpca}
#'
#' @param pcs   array_like \cr
#'              an array with two values indicating which two PCs should be plotted,
#'              by default the first two PCs are used, e.g., \eqn{c(1,2)}. \cr
#'
#' @param ................. .
#'
#' @seealso \code{\link{rsvd}}, \code{ggplot2}
#'
#' @examples #
#'
#

ggcorplot <- function( rpcaObj, pcs=c(1,2) ) {

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

  # Create ggplot
  g <- ggplot( circle ,aes( x , y)  )  + geom_path( size=1, colour="red"  )
  g <- g + geom_text(data=df, mapping=aes(x = a, y = b, label = labels, colour = labels ) ) +
      coord_fixed(ratio=1) + labs(x = noquote(PC1), y = noquote(PC2)) + guides(colour=FALSE)

  return( g )
}
