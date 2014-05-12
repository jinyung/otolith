#' Normalized elliptical Fourier analysis routine
#' 
#' @description a wrapper functions using the codes from J Claude's book,  
#'  for routine use to run normalized elliptical Fourier analysis.
#' @param A p x k x n array of raw semilandmarks configurations
#' @param plotpca Logical. Whether to plot the PCA for quick preliminary view. 
#'  Plot first 3 PCs. 
#' @param class factor giving the species/grouping, used for plotting purpose, 
#'  e.g. using sp value from \code{\link{routine1}}. 
#'  Used only when \code{plotpca=TRUE}
#' @return 
#'  \item{coeff}{matrix of all harmonic coefficients}
#'  \item{expvar}{explained variation of the NEF harmonics 
#'    (first harmonic is removed from the calculation, see reference)}  
#'  \item{size}{size of the semi-major axis of the first fitting ellipse. See reference}
#' @seealso 
#'  \code{\link{rGPA}}
#' @export
#' @references Claude J. (2008). Morphometrics with R. Springer.

rNEF <- function (A, plotpca=FALSE, class) {
  n <- dim(A)[3]
  p <- dim(A)[1]
  coeff <- matrix(NA, n, (p/2) * 4)
  size <- numeric()
  for (i in 1:n) {
    ef <- NEF(A[, , i])
    coeff[i, ] <- c(ef$A, ef$B, ef$C, ef$D)
    size[i] <- ef$size 
  }
  colnames(coeff) <- c(paste0("A", 1:(p/2)), paste0("B", 1:(p/2)), paste0("C", 1:(p/2)), paste0("D", 1:(p/2)))
  vharm <- apply(coeff, 2, var)
  variation <- apply(matrix(vharm, (p/2), 4), 1, sum)
  expvar <- round(cumsum(variation[-1]) / sum(variation[-1]), 3)[1:(p/2 - 1)]
  if (plotpca == TRUE) {
    pca <- prcomp(coeff)
    plotpca(pca, size=size, class=class)
  }
  return(list(coeff=coeff, expvar=expvar, size=size))
}