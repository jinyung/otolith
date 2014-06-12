#' Re-arrage semi-landmark configuration
#' 
#' @description To flip or change the order of semi-landmark configuration 
#' @param landdata semi-landmarks data generated from \code{\link{equaldist}} 
#'  functions. A matrix/dataframe consists of two columns of x and y coordinates
#' @param type character of \code{"ori"}, \code{"rev"}, \code{"flip"}, 
#'  or \code{"fliprev"} to set the configuration of the semi-landmarks. See details. 
#' @return matrix of xy coordinates of re-arranged semi-landmarks data. 
#' @details There are four possible configurations of semi-landmark sampling for 
#'  any set of otolith, depending on the side of the otolith and the initial 
#'  orietation of the image. This function can be used to standardize the sampling
#'  of semi-landmarks. 
#'  
#'  \code{type} argument change the semi-landmarks as if the otolith image
#'  was taken from different end on the left (\code{"rev"}), 
#'  or as if the otolith is from the other side (\code{"flip"}), 
#'  or as if the otolith is from the other side and image taken from different 
#'  end on the left (\code{"fliprev"}). \code{"ori"} do nothing on the semi-landmarks.  
#' @importFrom secr flip
#' @export

reland <- function (landdata, type=c("ori", "rev", "flip", "fliprev")) {
  require(secr)
  nd <- dim(landdata)[1] / 2
  landdatanew <- landdataF <- landdataFN <- 
    array(data = NA,  dim=c(dim(landdata)[1], 2))
  landdataF <- flip(landdata, lr=TRUE)
  landdataFN[1, ] <- landdataF[nd + 1, ]
  landdataFN[nd + 1, ] <- landdataF[1, ]
  landdataFN[2:nd, ] <- landdataF[nd:2, ]
  landdataFN[(nd + 2):(nd * 2), ] <- landdataF[(nd * 2):(nd + 2), ]
  if (type == "ori") {
    landdatanew <- landdata
  } else if (type == "flip") {
    landdatanew <- landdataFN
  } else {
    if (type == "fliprev")
      landdata <- landdataFN
    else if (type == "rev")
      landdata <- landdata
    else
      stop("incorrect type, please choose from \"ori\", \"rev\", \"flip\" or \"fliprev\"")
    landdatanew[1, ] <- landdata[nd + 1, ]
    landdatanew[nd + 1, ] <- landdata[1, ]
    landdatanew[2:nd, ] <- landdata [(nd + 2):(nd * 2),  ]
    landdatanew[(nd + 2):(nd * 2), ] <- landdata [2:nd, ]
  }
  return(landdatanew)
}