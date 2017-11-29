#' Extract outline
#' 
#' @description Extract outline from a fish otolith image
#' @param img an image array read by \code{\link{loadimg}}
#' @param threshold numeric. value of range 0-1 to set the theshold for conversion
#'  into binary image.
#' @param plot \code{no} = no plot; \code{overlay} = plot of outline overlaid 
#'  on image; \code{plain} = plot of plain outline
#' @return a matrix of outline coords. note the outline is in 8-connected format.
#' @importFrom EBImage bwlabel computeFeatures.shape rmObjects ocontour NumberOfFrames
#' @seealso
#' Which this function wraps: \code{\link[EBImage]{ocontour}}
#' 
#' Function that wraps this function: \code{\link{img2landmark}}
#' @export

extractout <- function (img, threshold=0.3,
                            plot=c("no", "overlay", "plain")) {
  # require(EBImage)
  # simple threshoding
  img.mask <- img > threshold
  # segmentation
  img.mask <- EBImage::bwlabel(img.mask)
  # calculate area of each object
  img.par <- EBImage::computeFeatures.shape(img.mask)
  if (dim(img.par)[1] > 1) {
    # remove all object except the one with largest area
    rm.index <- (1:dim(img.par)[1])[- which.max(img.par[, "s.area"])] 
    img.mask <- EBImage::rmObjects(img.mask, rm.index)
  }
  # extract outline
  img.con <- EBImage::ocontour(EBImage::opening(img.mask))
  # show plot
  plot <- match.arg(plot) 
  if (plot != "no") {
    if (plot == "overlay") {
      par(mar=rep(0, 4))
      x <- dim(img)[1]
      y <- dim(img)[2]
      plot(x, y, xlim=c(0, x), ylim=c(0, y), asp=1, type="n")
      rasterImage(img, 0, y, x, 0)
      polygon(img.con[[1]], border="yellow")
    } else if (plot == "plain") {
      par(mar=rep(0, 4)) 
      plot(img.con[[1]], asp=1, type='n', axes=FALSE, xlab="", ylab="")
      polygon(img.con[[1]])
    }
  }
  invisible(img.con[[1]])
}
