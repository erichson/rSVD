#'Tiger
#'
#' @description 1600x1200 grayscaled (8 bit [0-255]/255) image.
#'
#' @docType data
#'
#' @usage data('tiger')
#'
#' @format An object of class \code{\link[rsvd]{rsvd}}.
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
