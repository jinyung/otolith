#' Generalized Procrustes analysis routine
#' 
#' @description a wrapper function for routine use to
#'   run generalized Procrustes analysis.
#' @param A p x k x n array of raw semilandmarks configurations
#' @param fix a numeric vector giving the semi-landmakrs that not to be slided. 
#'  Default is all semilandmarks will be slided (See Note for issue).
#' @param plotpca logical. whether to plot PCA for preliminary assessment of 
#'  the data
#' @param class factor giving the species/grouping, used for plotting purpose
#'   and calculation of within class range of Procrustes distance, e.g. using sp
#'   value from \code{\link{routine1}} or using \code{\link{getclass}}
#' @return 
#'  \item{tanc}{Procrustes aligned configurations projected onto tangent shape
#'    space}
#'  \item{size}{centroid size}  
#'  \item{expvar}{summary of expalined variations by each PC}
#'  \item{score}{matrix contain the PC scores, to be used for training}
#'  \item{rdist}{see \code{\link{sprdist}}}
#'  \item{pca}{the PCA model, to be used by \code{\link{otopred}}}
#'  \item{meanshape}{meanshape of each species/ groups}
#'  \item{mshape}{meanshape of all species/ groups}
#' @importFrom Morpho procSym
#' @seealso 
#'  which this function wraps: \code{\link{procSym}}, \code{\link{sprdist}}
#' @export
#' @note There is an issue with \code{\link{otopred}} function, which cannot
#'  handle sliding semi-landmarks yet and hence it is recommended that all
#'  semi-landmarks should be fixed if the user intend to use the 
#'  \code{\link{otopred}} function (6 May 2014).

rGPA <- function(A, fix=NULL, plotpca=FALSE, class) {
  # check
  # require(Morpho)
  if (missing(class))
    stop("please provide class, you can use getclass()")
  p <- dim(A)[1]
  n <- dim(A)[3]
  slide <- 1:p
  if (!is.null(fix)) {
    if (length(slide[-fix]) == 0)
      slide <- NULL
    else 
      slide <- slide[-fix]
  }
  A <- aligne(A)
  gpa <- Morpho::procSym(A, SMvector=slide, outlines=c(1:p, 1))
  size <- gpa$size
  tanc <- gpa$orpdata
  pca <- prcomp(as.data.frame(t(matrix(tanc, p * 2, n))))
  if (plotpca) {
    plotpca(pca, size=size, class=class)
  }
  rdist <- sprdist(gpa$rotated, class=class)
  invisible(list(tanc=tanc, size=size, expvar= gpa$Variance, 
              score=gpa$PCscores, rdist=rdist$stat, pca=pca, 
              meanshape=rdist$meanshape, mshape=gpa$mshape))
}