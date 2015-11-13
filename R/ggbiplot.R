#devtools::use_package("ggplot2", type = "Suggests")
#devtools::use_package("plyr", type = "Suggests")
#devtools::use_package("scales", type = "Suggests")
#devtools::use_package("grid", type = "Suggests")


#library(ggplot2)
#library(plyr)
#library(scales)
#library(grid)


#
# This is a modified version of ggbiplot.r for rPCA
#
#  ggbiplot.r
#
#  Copyright 2011 Vincent Q. Vu.
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

#' Biplot for \code{rPCA} using ggplot2
#'
#' @param rpcObj          object containing the \code{sdev} component,
#' @param pcs             an array with two values indicating which two PCs should be plotted,
#'                        by default the first two PCs are used, e.g., \eqn{c(1,2)}
#' @param scale           covariance biplot (scale = 1), form biplot (scale = 0). When scale = 1, the inner product between the variables approximates the covariance and the distance between the points approximates the Mahalanobis distance.
#' @param obs.scale       scale factor to apply to observations
#' @param var.scale       scale factor to apply to variables
#' @param pc.biplot       for compatibility with biplot.princomp()
#' @param groups          optional factor variable indicating the groups that the observations belong to. If provided the points will be colored according to groups
#' @param ellipse         draw a normal data ellipse for each group
#' @param ellipse.prob    size of the ellipse in Normal probability
#' @param labels          optional vector of labels for the observations
#' @param labels.size     size of the text used for the labels
#' @param alpha           alpha transparency value for the points (0 = transparent, 1 = opaque)
#' @param circle          draw a correlation circle (only applies when prcomp was called with scale = TRUE and when var.scale = 1)
#' @param circle.prob     size of the circe in Normal probability
#' @param var.axes        draw arrows for the variables
#' @param varname.size    size of the text for variable names
#' @param varname.adjust  adjustment factor the placement of the variable names, >= 1 means farther from the arrow
#' @param varname.abbrev  whether or not to abbreviate the variable names
#' @param ...     arguments passed to or from other methods, see \code{\link[ggplot2]{ggplot}}.
#'
#' @seealso \code{\link{rpca}}, \code{\link[ggplot2]{ggplot}}
#'
#' @author The original implementation of \code{ggbiplot} was written by Vincent Q. Vu (2011).
#' @examples
#' #See ?rsvd

#' @export
ggbiplot <- function( rpcObj, pcs = c(1,2), scale = 1, pc.biplot = TRUE,
                     obs.scale = 1 - scale, var.scale = scale,
                     groups = NULL, ellipse = TRUE, ellipse.prob = 0.68,
                     labels = NULL, labels.size = 3, alpha = 1,
                     var.axes = TRUE,
                     circle = TRUE, circle.prob = 0.69,
                     varname.size = 3, varname.adjust = 1.5,
                     varname.abbrev = FALSE, ...)
{

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The package 'ggplot2' is needed for this function to work. Please install it.",
         call. = FALSE)
  }


  choices = pcs
  stopifnot(length(choices) == 2)

  # Get values from the rpca object
    nobs.factor <- sqrt(nrow(rpcObj$x) - 1)
    d <- rpcObj$sdev
    u <- sweep(rpcObj$x, 2, 1 / (d * nobs.factor), FUN = '*')
    v <- rpcObj$rotation


  # Scores
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*'))

  # Directions
  v <- sweep(v, 2, d^var.scale, FUN='*')
  df.v <- as.data.frame(v[, choices])

  names(df.u) <- c('xvar', 'yvar')
  names(df.v) <- names(df.u)

  if(pc.biplot) {
    df.u <- df.u * nobs.factor
  }

  # Scale the radius of the correlation circle so that it corresponds to
  # a data ellipse for the standardized PC scores
  r <- sqrt(stats::qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)

  # Scale directions
  v.scale <- rowSums(v^2)
  df.v <- r * df.v / sqrt(max(v.scale))

  # Change the labels for the axes
  if(obs.scale == 0) {
    u.axis.labs <- paste('standardized PC', choices, sep='')
  } else {
    u.axis.labs <- paste('PC', choices, sep='')
  }

  # Append the proportion of explained variance to the axis labels
  u.axis.labs <- paste(u.axis.labs,
                       sprintf('(%0.1f%% explained var.)',
                               100 * rpcObj$sdev[choices]^2/sum(rpcObj$sdev^2)))

  # Score Labels
  if(!is.null(labels)) {
    df.u$labels <- labels
  }

  # Grouping variable
  if(!is.null(groups)) {
    df.u$groups <- groups
  }

  # Variable Names
  if(varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  } else {
    df.v$varname <- rownames(v)
  }

  # Variables for text label placement
  df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar)) / 2)



  #Workaround for CRAN: Nulling
  xvar <- NULL # Setting the variables to NULL first
  yvar <- NULL # Setting the variables to NULL first
  muted <- NULL # Setting the variables to NULL first
  varname <- NULL # Setting the variables to NULL first
  angle <- NULL # Setting the variables to NULL first
  hjust <- NULL # Setting the variables to NULL first
  x <- NULL # Setting the variables to NULL first
  y <- NULL # Setting the variables to NULL first
  a <- NULL # Setting the variables to NULL first
  b <- NULL # Setting the variables to NULL first




  # Base plot
  g <- ggplot2::ggplot(data = df.u, ggplot2::aes(x = xvar, y = yvar)) +
    ggplot2::xlab(u.axis.labs[1]) + ggplot2::ylab(u.axis.labs[2]) + ggplot2::coord_equal()

  if(var.axes) {
    # Draw circle
    if(circle)
    {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * sin(theta))
      g <- g + ggplot2::geom_path(data = circle, color = muted('white'),
                         size = 1/2, alpha = 1/3)
    }

    # Draw directions
    g <- g +
      ggplot2::geom_segment(data = df.v,
                  ggplot2::aes(x = 0, y = 0, xend = xvar, yend = yvar),
                   arrow = grid::arrow(length = grid::unit(1/2, 'picas')),
                   color = scales::muted('red'))
  }

  # Draw either labels or points
  if(!is.null(df.u$labels)) {
    if(!is.null(df.u$groups)) {
      g <- g + ggplot2::geom_text(ggplot2::aes(label = labels, color = groups),
                         size = labels.size)
    } else {
      g <- g + ggplot2::geom_text(ggplot2::aes(label = labels), size = labels.size)
    }
  } else {
    if(!is.null(df.u$groups)) {
      g <- g + ggplot2::geom_point(ggplot2::aes(color = groups), alpha = alpha)
    } else {
      g <- g + ggplot2::geom_point(alpha = alpha)
    }
  }

  # Overlay a concentration ellipse if there are groups
  if(!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))

    ell <- plyr::ddply(df.u, 'groups', function(x) {
      if(nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- stats::var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(stats::qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'),
                 groups = x$groups[1])
    })
    names(ell)[1:2] <- c('xvar', 'yvar')
    g <- g + ggplot2::geom_path(data = ell, ggplot2::aes(color = groups, group = groups))
  }

  # Label the variable axes
  if(var.axes) {
    g <- g +
      ggplot2::geom_text(data = df.v,
                      ggplot2::aes(label = varname, x = xvar, y = yvar,
                    angle = angle, hjust = hjust),
                color = 'darkred', size = varname.size)
  }


  return(g)
}
