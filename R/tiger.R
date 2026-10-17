#' @title Tiger image
#'
#' @description
#' A grayscale tiger photograph for illustrating low-rank approximation and
#' image reconstruction. Pixel intensities have been divided by 255 to lie
#' between zero and one.
#'
#' @docType data
#' @usage data("tiger")
#' @format A numeric matrix with 1600 rows and 1200 columns. There are no labels
#' or color channels; each entry is a grayscale intensity.
#' @details
#' The stored orientation is suitable for \code{\link[graphics]{image}} as
#' shown below. Matrix rows and columns correspond to the horizontal and
#' vertical axes used by that function, respectively. For a decomposition,
#' use \code{tiger} directly and reconstruct the matrix in the same orientation.
#'
#' The image is a useful visual example rather than an exactly low-rank matrix.
#' Increasing the approximation rank preserves more texture and fine detail,
#' while storing only the factors can require much less space than storing
#' all pixels. See \code{\link{rsvd}}, \code{\link{rid}} and \code{\link{rcur}}
#' for decomposition examples.
#'
#' @keywords datasets
#' @references S. Taheri (2006). Panthera tigris altaica. Online photograph.
#' @source \href{https://en.wikipedia.org/wiki/File:Siberischer_tiger_de_edit02.jpg}{Wikimedia image and attribution}
#' @examples
#' data("tiger")
#' dim(tiger)
#' range(tiger)
#' image(seq_len(nrow(tiger)), seq_len(ncol(tiger)), tiger,
#'       col = gray((0:255) / 255), asp = 1,
#'       axes = FALSE, useRaster = TRUE)

"tiger"
