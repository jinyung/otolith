#' Select training data dimension
#' 
#' @description 
#' a utility function to trim down the matrix of coefficients (NEF)/ PC scores (GPA)
#' @param train Matrix of EFA coefficients or PC scores result from \code{\link{rNEF}}
#'   or \code{\link{rGPA}} function (called as \code{object$coeff} or 
#'   \code{object$score}) to be trimmed down
#' @param har Numeric. Number of harmonics (for coefficient matrix)
#' @param pc Numeric. Number of pc (for PC score matrix)
#' @return 
#'  coefficients/ PC scores of selected number ready to be used for training
#' @seealso 
#'  input from: \code{\link{rNEF}}, \code{\link{rGPA}}
#'  
#'  optimising \code{pc/har} number: \code{\link{pccv}}, \code{\link{harcv}}
#' @export

selectdim <- function(train, har=NULL, pc=NULL) {
  if (is.null(har) == FALSE & is.null(pc) == FALSE)
    stop ("please give either har or pc value")
  if (is.null(har) == TRUE & is.null(pc) == TRUE)
    stop ("please give either har or pc value")
  if (is.null(har) == FALSE) {
    if (har>(dim(train)[2] / 4 - 1))
      stop("har number too large")
    train <- train[, c(paste0("A", 2:har), paste0("B", 2:har), 
              paste0("C", 2:har), paste0("D", 2:har))]
  }
  if (is.null(pc) == FALSE) {
    if (pc > (dim(train)[2] - 4))
      stop("pc number too large")  
    train <- train[, c(paste0("PC",  1:pc))]
  }
  return(train)  
}