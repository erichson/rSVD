#devtools::use_package("ggplot2", type = "Suggests")
#library("ggplot2")

#' @title Screeplot
#
#' @description Creates a screeplot, variables and individual factor maps to
#'              summarize the results of the \code{rpca()} function.
#'
#' @param x   Object returned by the \code{rpca()} function.
#'
#' @param ...     Additional arguments passed to the individual plot functions (see below).
#'
#' @param ................. .
#'
#' @seealso \code{\link{ggscreeplot}}, \code{\link{ggcorplot}} , \code{\link{ggindplot}}  
#'
#' @examples #
#'
#'@export
plot.rpca <- function(x, ... ) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The package 'ggplot2' is needed for this function to work. Please install it.",
         call. = FALSE)
  } 
   
  print(ggscreeplot(x, ... ))

  print(ggcorplot(x, ... ))

  print(ggindplot(x, ... ))

}


