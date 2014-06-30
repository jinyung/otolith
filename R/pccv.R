#' Optimize PC number
#' 
#' @description A wrapper function to optimize dimension reduction using 
#'  multiple-run k-fold cross validation. Different number of PCs are used 
#'  in training the classification model and the result (misclassification rate) 
#'  is plotted out.
#' @param X matrix of PC scores from \code{\link{rGPA}} object (the \code{score} value)
#' @param Y a factor giving the grouping, e.g. the \code{sp} value from 
#'  \code{\link{routine1}} object
#' @param pc numeric. number of PCs to be tested.
#' @param saveplot logical. The plot will be saved if \code{TRUE}, 
#'  and no plot is displayed in the window. Note that outliers of
#'   boxplot are not plotted.
#' @param plotsize numeric. Plot size for the plot saved into file 
#'  (used only when \code{saveplot=TRUE}), unit in pixel.
#' @param method argument to be passed to \code{\link{mrkfcv}}
#' @param run argument to be passed to \code{\link{mrkfcv}}
#' @param k argument to be passed to \code{\link{mrkfcv}}
#' @return a matrix and a plot are returned, giving the misclassification rate
#'  of the range of PCs tested.
#' @seealso 
#'  Similar: \code{\link{harcv}}, \code{\link{threcv}}
#'  
#'  Which this function wraps: \code{\link{mrkfcv}}
#' @export

pccv <- function(X, Y, pc, saveplot=FALSE, plotsize=1000, 
                 method=c("lda", "tree", "plsda"), run=30, k=5) {
  error <- numeric()
  errorsd <- numeric()
  temp <- NULL
  misclass <- data.frame(matrix(NA, k * run * pc, 2))
  misclass[, 2] <- factor(rep(paste0("PC", 1:pc), each=run * k), 
                          levels= paste0("PC", 1:pc), ordered=TRUE, labels=1:pc)
  for (i in 1:pc) {
    cat ("\r                       (Evaluating PC: 1 - ", i,") | pccv progress: [",  
         round(i / pc * 100), "%]       ", sep="")
    flush.console()
    temp <- mrkfcv(X=data.frame(selectdim(X, pc= i)), Y= Y, suppress="text", 
                   method=method, k=k, run=run) 
    # have to put X into data.frame because kfcv extract info using dim(), 
    # and when pc=1, it will becaome a vector and dim() will give error
    misclass[((i - 1) * run * k + 1):(i * run * k), 1] <- temp$misclass
    error[i] <- 100 - temp$accuracy
    errorsd[i] <- temp$accu.sd
  }
  if (saveplot == TRUE) {
    filename <- "pc-optimization.tif"
    tiff(filename, plotsize, plotsize, res = 172, compression = "lzw")
  }  
  boxplot(misclass[, 1] ~ misclass[, 2], outline=FALSE, 
          ylab="Cross-validated misclassification rate", 
          xlab="No. of PC used in training")
  if (saveplot == TRUE) {
    dev.off()
    cat("\nThe plot is saved at:", 
        paste(getwd(), filename, sep="/"), "\n")
  }
  result <- cbind(error, errorsd)
  rownames(result) <- paste0("PC1-", 1:pc)
  cat("\nEvaluation completed\n")
  return(result)
}