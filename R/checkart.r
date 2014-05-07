#' Check outline artifact
#' 
#' @description Checking whether the outline containing any pixel artifact 
#'  that will cause error in subsequent sampling of outline
#' @param outdata outline data to be checked. A matrix consists of xy coordinates
#' @return 
#'  \item{artifact}{logical. Whether the outline contains pixel artifact}
#'  \item{index}{vector containg the index of pixel that is problematic}  
#' @seealso \code{\link{orderOutline}}
#' @export

checkart <- function (outdata) {
  neigh <- NULL
  artifact <- FALSE 
  for (i in 1:(dim(outdata)[1] - 1)){
    neigh[i] <- 0
    xtemp <- outdata[outdata[, 1] == outdata[i, 1], ]
    ytemp <- outdata[outdata[, 2] == outdata[i, 2], ]
    if (any(xtemp[, 2] == (outdata[i, 2] + 1)))
      neigh[i] <- neigh[i] + 1
    if (any(xtemp[, 2] == (outdata[i, 2] - 1)))
      neigh[i] <- neigh[i] + 1
    if (any(ytemp[, 1] == (outdata[i, 1] + 1)))
      neigh[i] <- neigh[i] + 1
    if (any(ytemp[, 1] == (outdata[i, 1] - 1)))
      neigh[i] <- neigh[i] + 1
  }
  if (any(neigh > 2,  na.rm=T))
    artifact <- TRUE
  index <- which(neigh > 2)
  return(list("artifact"=artifact,  "index"=index))  
}