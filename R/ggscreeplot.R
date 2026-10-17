#' @title PCA scree plot and explained variance
#'
#' @description Show how much variation each retained principal component
#' from \code{\link{rpca}} captures in the processed data. Plot individual
#' component variances, their proportions of total variance,
#' or the cumulative proportion to assess how much variation the retained
#' representation captures. The optional \pkg{ggplot2} package is required.
#'
#' @param rpcaObj An object returned by \code{\link{rpca}}.
#' Scores are not needed, so fits made with \code{retx = FALSE} are supported.
#' The variance quantities are real-valued for both real and complex PCA.
#'
#' @param type Character string selecting the vertical-axis
#' quantity. \code{"var"} (the default) plots \code{sdev^2}; \code{"ratio"}
#' plots \code{sdev^2 / var}; \code{"cum"} plots the cumulative sum of those
#' ratios; and \code{"eigenvals"} plots \code{eigvals}. For an \code{rpca}
#' object, \code{"var"} and \code{"eigenvals"} give the same values with
#' different axis labels. Unambiguous partial matches are accepted.
#'
#' @details The horizontal axis follows the retained components in their
#' fitted order. The variance and eigenvalue views show how much variation
#' each component captures, in squared units of the processed input. The
#' ratio view expresses each amount as a proportion of total variance; the
#' cumulative view adds those proportions through each displayed component.
#'
#' Both proportion views use \code{rpcaObj$var} as their denominator. This
#' total includes energy outside the retained components, so the cumulative
#' curve can end below one when the fit retains only a low-rank approximation.
#' The proportions lie between zero and one when total processed variance is
#' positive. They are undefined for an input with zero total processed
#' variance.
#'
#' Interpretation follows the preprocessing used in \code{\link{rpca}}.
#' With \code{center = FALSE}, the quantities use uncentered sums of squares,
#' which include the contribution of any nonzero column means. They need not
#' equal variances computed after subtracting the means of the fitted score
#' columns.
#'
#' @return A \pkg{ggplot2} plot object. Print the object to draw the
#' plot, or add layers, labels, and themes with \pkg{ggplot2}.
#'
#' @seealso \code{\link{rpca}}, \code{\link{ggbiplot}},
#' \code{\link{ggcorplot}}, \code{\link{ggindplot}},
#' \code{\link[ggplot2]{ggplot}}
#'
#' @author N. Benjamin Erichson
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(42)
#'   fit <- rpca(iris[, 1:4], k = 2, retx = FALSE)
#'   # The curve uses all processed variance, including omitted components.
#'   print(ggscreeplot(fit, type = "cum"))
#' }
#'
#' @export
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
                    cum = 'Cumulative proportion',
                    eigenvals = 'Eigenvalues',
                    stop("Selected plot option is not supported!")
  )
  
  df <- data.frame(x = 1:length(rpcaObj$sdev), y = y)
  
  #Workaround for CRAN: Nulling
  x <- NULL # Setting the variables to NULL first
  
  g <- ggplot2::ggplot(data = df, ggplot2::aes(x = x, y = y) ) +
    ggplot2::xlab('Principal components') + ggplot2::ylab( y.label ) +
    ggplot2::geom_line(lwd=0.5, color = 'black', linetype="dashed") +
    ggplot2::geom_point(size=2, color = '#ef3b2c') + 
    ggplot2::guides(colour="none")
  
  g <- g + ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank())
  
  
  return(g)
}
