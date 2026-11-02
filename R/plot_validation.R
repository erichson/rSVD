# Plot indices must agree with the component numbers printed on the axes.
# R otherwise silently truncates fractional numeric subscripts.
check_plot_pcs <- function(rpcaObj, pcs) {
  if(length(pcs) != 2L || !is.numeric(pcs) || is.complex(pcs) ||
     any(!is.finite(pcs)) || any(pcs != floor(pcs)) ||
     any(pcs < 1) || any(pcs > ncol(rpcaObj$rotation))) {
    stop("Selected PC is not valid: use two integer component indices.",
         call. = FALSE)
  }
  invisible(NULL)
}

check_plot_real <- function(values) {
  if(is.complex(values)) {
    stop("PCA factor maps require real-valued scores and directions.",
         call. = FALSE)
  }
  invisible(NULL)
}

# Recycling a shorter vector can silently assign observations to wrong groups
# or label variables incorrectly. Require the documented one-to-one mapping.
check_plot_length <- function(values, expected, name) {
  if(!is.null(values) && length(values) != expected) {
    stop("'", name, "' must have one entry per plotted row (", expected, ").",
         call. = FALSE)
  }
  invisible(NULL)
}
