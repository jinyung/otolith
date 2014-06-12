#' Plot PCA bubble plot
#' 
#' @description For preliminary visual assessment of data groups in PCA. 
#' @param pca \code{\link{prcomp}} object 
#' @param size optional. numeric. geometric size, e.g. \code{size} value from 
#'  \code{\link{rGPA}} or \code{\link{rNEF}}. size is overlaid as bubble size 
#'  in the plot.
#' @param class factor. giving the species/grouping. NOTE: currently only 
#'  support up to 12 levels. 
#' @param saveplot logical. whether to save the plot
#' @param plotsize numeric. plot size for the plot saved into file 
#'  (used only when \code{saveplot=TRUE}), unit in pixel.
#' @param sizeamp numeric. to adjust the size of the bubble.
#' @return Just the plot, consist of the first 3 PCs. 
#' @export
#' @seealso 
#'  Function that wraps this function: \code{\link{rNEF}}, \code{\link{rGPA}}, 
#'  \code{\link{routine1}}
#' @keywords internal  

plotpca <- function (pca, size=NULL, sizeamp=3, class, saveplot=FALSE, plotsize=1000) {
  if (length(levels(class)) > 12)
    stop("now this function supports up to 12 groups only")
  range01 <- function(x){((x-min(x))/(max(x)-min(x)))}
  if (!is.null(size)) { 
    radius <- sqrt((size) / pi)
    radius <- range01(radius) * sizeamp
  } else {
    radius <- 1
  }
  bg <- rainbow(length(levels(class)))[as.numeric(class)]
  if (saveplot == FALSE) {
    dev.new(width=10,  height=10)
  } else {
    filename <- paste0(bquote(pca), "-pca.tif")
    tiff(filename, plotsize, plotsize, res=172)
  }
  par (mfrow=c(2, 2), mar=c(4, 4, 1, 1))
  plot(pca$x[, c(1, 2)], asp=1,  cex= radius, pch=21, bg=bg, 
       xlab= paste("PC1", "(", round(summary(pca)$
       importance[, 1][2] * 100), "%)"), ylab= paste("PC2", 
       "(", round(summary(pca)$importance[, 2][2] * 100), "%)"))
  abline(v=0, h=0,  col="gray", lty=2)
  plot(pca$x[, c(1, 3)], asp=1, cex= radius,  pch=21, bg=bg, 
       xlab=paste("PC1", "(", round(summary(pca)$importance[, 1][2] * 100),
       "%)"), ylab=paste("PC3", "(", round(summary(pca)$
       importance[, 3][2] * 100), "%)"))
  abline(v=0, h=0, col="gray", lty=2)
  plot(pca$x[, c(2, 3)],  asp=1,  cex= radius,  pch=21, bg=bg, 
       xlab=paste("PC2", "(", round(summary(pca)$importance[, 2][2] * 100), 
       "%)"), ylab=paste("PC3", "(", round(summary(pca)$
       importance[, 3][2] * 100), "%)"))
  abline(v=0, h=0, col="gray", lty=2)
  plot(1, 1, type='n',  axes=FALSE,  frame=FALSE,  ann=FALSE)
  legend("center", fill=1:length(levels(as.factor(class))),  
         title="Species", legend= levels(as.factor(class)))
  if (saveplot == TRUE) {
    dev.off()
    cat("The plot is saved at:",  
        paste(getwd(), filename, sep="/"), "\n\n")
  }
}