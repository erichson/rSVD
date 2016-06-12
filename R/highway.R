#'Highway
#'
#' @description 176x144 grayscaled (8 bit [0-255]) surveillance video with 200 frames.
#' Each frame is flattened, and stored as a column vector.
#'
#' @docType data
#'
#' @usage data(highway)
#'
#' @format An object of class \code{"rsvd"}.
#'
#' @keywords video
#'
#' @references N. Goyette, P.-M. Jodoin, F. Porikli, J. Konrad, and P. Ishwar,
#'            changedetection.net: A new change detection benchmark dataset,
#'            in Proc. IEEE Workshop on Change Detection (CDW-2012) at CVPR-2012,
#'            Providence, RI, 16-21 Jun., 2012.
#'
#' @source \href{http://changedetection.net/}{changedetection.net}
#'
#' @examples
#' library(rsvd)
#' data(highway)
#'
#' # Reshape and display the 100th frame:
#' frame <- matrix(highway[,100], ncol=144, nrow=176)
#' image(frame, col = gray((0:255)/255))
#'
"highway"
