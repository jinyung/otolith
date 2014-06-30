#' Plot allometric deformation
#' 
#' @description A wrapper function to draw min, mean, and max shape of a species
#'   out of all configurations of selected group. Min and max shapes are
#'   displayed as TPS deformations from the mean shape
#' @param gpa a \code{\link{rGPA}} object
#' @param group factor giving the groupings, must be in same order as \code{A}
#'   of \code{\link{rGPA}} object
#' @param target character. group that you want to draw, one of the level of
#'   \code{group}
#' @param main character. title of the plot. default is using the \code{target}
#' @param n numeric. no. of grids used in TPS visualization, argument passed to
#'   \code{tps}
#' @param col color of the outline
#' @param saveplot logical. Whether to save the plot. No plot will be displayed
#'   in window if \code{saveplot=TRUE}
#' @param plotsize width of the plot in pixel. Used only when
#'   \code{saveplot=TRUE}
#' @return just the plot
#' @export

plotallo <- function (gpa, group, target, col = 1, main, 
                      n = 40, saveplot = FALSE, plotsize = 1300 
                      ) {
  A <- gpa$tanc
  cens <- gpa$size
  max.spec <- A[, , which(group == target)][, , which(log(cens[which(
              group == target)]) == max(log(cens[which(group == target)])))]
  min.spec <- A[, , which(group == target)][, , which(log(cens[which(
              group == target)]) == min(log(cens[which(group == target)])))]
  ref.spec <- mshape(A[, , which(group == target)])
  if (missing(main))
    main <- target
  if (saveplot == FALSE)    
    dev.new(width=13, height=4)
  if (saveplot == TRUE) {
    filename <- paste0(deparse(substitute(gpa)), "-plotallo(", 
                       Sys.Date(), ").tif")
    tiff(filename, plotsize, (plotsize - plotsize / 13) / 3, res = 144, 
         compression = "lzw")  
  }
  par(mfrow=c(1, 3), mar=c(0, 0, 1, 0) + 0.5, font.main=3, cex.main=2)
    tps2(ref.spec, min.spec, n, A)
      title(paste(main, "(Min)"))
      polygon(min.spec, border=col, lwd=2)
    tps2(ref.spec, ref.spec, n, A)
      title(paste(main, "(Average)"))
      polygon(ref.spec, border=col, lwd=2)
    tps2(ref.spec, max.spec, n, A)
      title(paste(main, "(Max)"))
      polygon(max.spec, border=col, lwd=2)
  if (saveplot == TRUE) {
    dev.off()
    cat("The plot is saved at:", 
        paste(getwd(), filename, sep="/"), "\n\n")
  }
}