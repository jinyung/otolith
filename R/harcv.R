#' Optimize harmonics number
#' 
#' @description modified version of \code{\link{pccv}} function to optimize 
#'  dimension reduction.Different number of harmonics are used in training 
#'  the classification model and the result (misclassification rate) 
#'  is plotted out.
#' @param X matrix of EFA coefficients from \code{\link{rNEF}} 
#'  object (the \code{coeff} value)
#' @param Y a factor giving the grouping, e.g. the \code{sp} value from 
#'  \code{\link{routine1}} object
#' @param har numeric. number of harmonics to be tested.
#' @param saveplot logical. The plot will be saved if \code{TRUE}, 
#'  and no plot is displayed in the window. Note that outliers of
#'   boxplot are not plotted.
#' @param plotsize numeric. Plot size for the plot saved into file 
#'  (used only when \code{saveplot=TRUE}), unit in pixel.
#' @param method argument to be passed to \code{\link{mrkfcv}}
#' @param run argument to be passed to \code{\link{mrkfcv}}
#' @param k argument to be passed to \code{\link{mrkfcv}}
#' @return a matrix and a plot are returned, giving the misclassification rate
#'  of the range of harmonics tested.
#' @seealso 
#'  Similar: \code{\link{pccv}}
#'  
#'  Which this function wraps: \code{\link{mrkfcv}}
#' @export

harcv <- function(X, Y, har, saveplot=FALSE, plotsize=1000, 
                  method="lda", run=30, k=5) {
  if(har > (dim(X)[2] / 4 - 1))
    stop("harmonics number too large", call.=FALSE)
  error <- numeric()
  errorsd <- numeric()
  temp <- NULL
  misclass <- data.frame(matrix(NA,k * run * (har-1),2))
  misclass[, 2] <- factor(rep(paste0("har", 2:har), each=run * k), 
                   levels= paste0("har", 2:har), ordered=TRUE, labels=2:har)
  cat("------ evaluating dimension reduction ------\n")
  for (i in 2:har) {
    cat("harmonics", i, "\n")
    flush.console()
    temp <- mrkfcv(X= selectdim(X, har= i), Y=Y, suppress=TRUE, 
                   method=method, k=k, run=run)
    misclass[((i - 2) * run * k + 1):((i-1) * run * k), 1] <- temp$misclass
    error[i-1] <- 100 - temp$accuracy
    errorsd[i-1] <- temp$accu.sd
  }
  cat("--- dimension reduction evaluation ended---\n\n")
  if (saveplot == TRUE) {
    filename <- "harmonics-optimization.tif"
    tiff(filename, plotsize, plotsize, res=172)
  }  
  boxplot(misclass[, 1] ~ misclass[, 2], outline=FALSE, 
          ylab="Cross-validated misclassification rate",
          xlab="No. of Harmonics used in training")
  if (saveplot == TRUE) {
    dev.off()
    cat("The plot is saved at:", 
        paste(getwd(), filename, sep="/"), "\n\n")
  }
  result <- cbind(error, errorsd)
  rownames(result) <- paste0("har-", 1:har)
  return(result)
}