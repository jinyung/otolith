#' Calculate within species Procrustes distance
#' 
#' @description A wrapper to calculate the within species Procrustes 
#'  (or Riemannian, as used by \code{\link{Morpho}}) distance
#' @param A p x k x n array of landmark configuration
#' @param class a factor giving the species or the grouping of \code{A}, 
#'  must be of same length and order as \code{A}. can be obtained using
#'  \code{\link{getclass}}
#' @return 
#'  \item{stat}{matrix of result of calculated within species Procrustes distance}
#'  \item{meanshape}{array of meanshape configurations of each species/ group}  
#' @importFrom Morpho kendalldist rotonto
#' @seealso 
#'  function that wraps this function: \code{\link{rGPA}}
#'  
#'  which this function wraps: \code{\link{kendalldist}},  \code{\link{rotonto}}, 
#'   \code{\link{mshape}}
#' @export


sprdist <- function(A, class) {
  # require(Morpho)
  if (!is.factor(class))
    class <- factor(class)
  grlevel <- levels(class)
  sdist.mean <- numeric(length=length(grlevel)) # to save dist for each class 
  names(sdist.mean) <- grlevel
  sdist.min <- sdist.mean
  sdist.max <- sdist.mean
  cat("\n--------Calculating Riemannian distance range--------\n")
  meanshape <- array(NA, dim=c(dim(A)[1], 2, length(grlevel)), dimnames=list(NULL, NULL, grlevel))
  for (i in 1:length(grlevel)) {
    cat("calculating ", grlevel[i], "\n")
    flush.console()
    Aindex <- which(class == grlevel[i])
    tempA <- A[, , Aindex]
    rdist <- numeric() # temporary to save all dist of a species
    ms <- mshape(tempA) 
    rdist <- numeric()
    for (j in 1:length(Aindex)) {
      pair <- Morpho::rotonto(x=tempA[, , j], y=ms)
      rdist[j] <- Morpho::kendalldist(pair$X, pair$Y)
    }
    meanshape[, , i] <- ms
    sdist.mean[i] <- mean(rdist)
    sdist.min[i] <- min(rdist)
    sdist.max[i] <- max(rdist)
  }
  cat("-----Calculating Riemannian distance range ended-----\n\n")
  return(list(stat=data.frame(cbind(sdist.mean, sdist.min, sdist.max)),
              meanshape=meanshape))
} 