#devtools::use_package('ggplot2')
#library('ggplot2')

#' @title PCA summary plots
#
#' @description Creates the available scree, variable, and individual factor
#'              plots for an object returned by \code{\link[rsvd]{rpca}}.
#'
#' @details A scree plot is always drawn. The variable factor map requires at
#' least two retained components. The individual factor map additionally
#' requires fitted scores, obtained with \code{retx = TRUE}.
#'
#' @param x   Object returned by the \code{\link[rsvd]{rpca}} function.
#'
#' @param ... Named arguments for \code{\link{ggscreeplot}},
#' \code{\link{ggcorplot}}, or \code{\link{ggindplot}}. Each argument is passed
#' only to the helpers that accept it; for example, \code{type} selects the
#' scree-plot quantity and \code{pcs} selects components for the factor maps.
#' Argument names must be exact and unique. Unknown names are rejected.
#' An explicit \code{pcs} value must contain two finite numeric component
#' indices within the retained range, even when a factor map is unavailable.
#'
#' @return The last plot drawn, returned invisibly.
#'
#' @seealso \code{\link{ggscreeplot}}, \code{\link{ggcorplot}} , \code{\link{ggindplot}}  
#'
#' @examples #
#'
#'@export
plot.rpca <- function(x, ... ) {
  
  if (!requireNamespace('ggplot2', quietly = TRUE)) {
    stop("The package 'ggplot2' is needed for this function to work. Please install it.",
         call. = FALSE)
  } 
   
  args <- list(...)
  if(length(args)) {
    arg_names <- names(args)
    if(is.null(arg_names) || anyNA(arg_names) || any(arg_names == "")) {
      stop("Additional plot arguments must be named.", call. = FALSE)
    }
    if(anyDuplicated(arg_names)) {
      stop("Each plot argument must be supplied once.", call. = FALSE)
    }
    accepted <- unique(unlist(lapply(list(ggscreeplot, ggcorplot, ggindplot),
      function(helper) setdiff(names(formals(helper)), "rpcaObj"))))
    unknown <- setdiff(arg_names, accepted)
    if(length(unknown)) {
      stop("Unknown plot argument(s): ", paste(unknown, collapse = ", "),
           call. = FALSE)
    }
    if("pcs" %in% arg_names) {
      pcs <- args$pcs
      if(length(pcs) != 2L || !is.numeric(pcs) || is.complex(pcs) ||
         any(!is.finite(pcs)) || any(pcs < 1) ||
         any(pcs > ncol(x$rotation))) {
        stop("Selected PC is not valid.", call. = FALSE)
      }
    }
  }

  draw <- function(helper) {
    selected <- args[names(args) %in% names(formals(helper))]
    print(do.call(helper, c(list(rpcaObj = x), selected)))
  }

  last_plot <- draw(ggscreeplot)
  if(ncol(x$rotation) >= 2L) {
    last_plot <- draw(ggcorplot)
    if(!is.null(x$x)) last_plot <- draw(ggindplot)
  }
  invisible(last_plot)

}


