#' Check outline order
#' 
#' @description Checking whether the outline XY coordinates are in order, i.e. 
#'  whether successive points' indices are in successive order? 
#'  This is important because 'equaldist' function assumes the coordinates 
#'  are ordered.
#' @param test a matrix consists of x and y coordinates, the outline data to be checked.
#' @return logical. Whether the outline is ordered.
#' @export

checkOrder <- function(test) {
  final <- NULL # result for every pixel
  result <- NULL # final overall result 
  for (i in 1:(dim(test)[1] - 1)){ 
    if (test[i, 1] == test[i + 1, 1]) { 
      if (test[i, 2] == test[i + 1, 2] + 1)
        final[i] <- TRUE
      else if (test[i, 2] == test[i + 1, 2] - 1)
        final[i] <- TRUE
      else 
        final[i] <- FALSE
    } else if (test[i, 2] == test[i + 1, 2]) {
      if (test[i, 1] == test[i + 1, 1] + 1)
        final[i] <- TRUE
      else if (test[i, 1] == test[i + 1, 1] - 1)
        final[i] <- TRUE
      else 
        final[i] <- FALSE
    } else {
      final[i] <- FALSE
    }
  }
  if (any(final == FALSE))
    result <- FALSE
  else
    result <- TRUE
  return(result)
}