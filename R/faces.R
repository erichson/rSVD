#'kaggle faces
#'
#' @description Grayscaled faces (8 bit [0-255]), 96x96 size. A few images of several different people
#' and 1500 total images. Each face image is flattened and stored as a column vector.
#'
#' @note This data-set is only a subset of the original kaggle face data-set,
#'  which includes about 7000 faces.
#'
#' @docType data
#'
#' @usage data(faces)
#'
#' @format An object of class \code{"rsvd"}.
#'
#' @keywords faces
#'
#' @references Donated by Dr. Yoshua Bengio for the kaggle facial keypoints detection competion.
#'
#' @source \href{https://www.kaggle.com/c/facial-keypoints-detection}{kaggle}
#'
#' @examples
#' library(rsvd)
#' data(faces)
#'
#' #Display 10th face image
#' img <- matrix(rev(faces[,10]), nrow=96, ncol=96)
#' image(img, col=gray((0:255)/255))
#'
"faces"
