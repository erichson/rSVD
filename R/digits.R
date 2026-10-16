#'Digits
#'
#' @description Subsampled MNIST database of handwritten digits. This smaller dataset has 3000 samples for each
#' of the digits corresponding to the class labels 0,1,2,3. Each 28x28 image patch is stored as a flattened row vector.
#'
#' @docType data
#'
#' @usage data('digits')
#'
#' @format An object of class \code{\link[rsvd]{rsvd}}.
#'
#' @keywords pattern recognition
#'
#' @references Y. LeCun, L. Bottou, Y. Bengio, and P. Haffner. 
#' "Gradient-based learning applied to document recognition." 
#' Proceedings of the IEEE, 86(11):2278-2324, November 1998.
#'
#' @examples
#' \dontrun{
#' library('rsvd')
#' data('digits')
#'
#' #Display first digit
#' digit <- matrix(digits[1,], nrow = 28, ncol = 28)
#' image(digit[,28:1], col = gray(255:0 / 255))
#' }   
#'
"digits"
