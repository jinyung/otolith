#' Order outline pixels
#' 
#' @description Order an outline so that the successive points' indices are 
#'  in successive order. This is important because \code{\link{equaldist}} function 
#'  assumes the coordinates are ordered. 
#' @param outdata outline data that need to be re-ordered. 
#'  A matrix consists of xy coordinates
#' @return re-ordered outline coordinates
#' @seealso \code{\link{checkart}}
#' @export

orderOutline <- function (outdata) {
  new.out <- matrix(data = NA,  dim(outdata)[1] + 1, 2)
  new.out[1, ] <- c(t(outdata[1, ]))
  for (i in 1:dim(outdata)[1]) {
    j <- i + 1
    xeq <- outdata [(outdata[, 1] == new.out[i, 1]), ] 
    yeq <- outdata [(outdata[, 2] == new.out[i, 2]), ] 
    tempout <- xeq[(xeq[, 2] == new.out[i, 2] + 1)|(xeq[, 2] == 
                new.out[i, 2] - 1), ] 
    tempout <- rbind(tempout,  yeq[(yeq[, 1] == new.out[i, 1] + 1)|
                (yeq[, 1] == new.out[i, 1] - 1), ]) 
    if (i == 1) {
      nextneighbor <- tempout[1, ]
    } else {
      nextneighbor <- subset(tempout,  !(tempout[, 1] ==  new.out[i - 1, 1] & 
                        tempout[, 2] ==  new.out[i - 1, 2])) 
    }
    new.out[j, ] <-c(t(nextneighbor))
  }
  new.out <- new.out[!duplicated(new.out), ]
  new.out
}