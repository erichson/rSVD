#'Tiger
#'
#' @description A 1600 by 1200 grayscale image with pixel intensities
#' normalized from 8-bit values to the interval [0, 1].
#'
#' @docType data
#'
#' @usage data('tiger')
#'
#' @format A numeric matrix with 1600 rows and 1200 columns.
#'
#' @keywords image
#'
#' @references S. Taheri (2006). "Panthera tigris altaica", (Online image)
#'
#' @source \href{https://en.wikipedia.org/wiki/File:Siberischer_tiger_de_edit02.jpg}{Wikimedia}
#'
#' @examples
#' \dontrun{
#' library('rsvd')
#' data('tiger')
#'
#' #Display image
#' image(tiger, col = gray((0:255)/255))
#' }
"tiger"
