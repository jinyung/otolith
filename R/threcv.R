#' Optimize threshold on posterior probability
#' 
#' @description To estimate the effect of setting a threshold 
#'  on posterior probability to reject "unsure" predictions.
#' @details only \code{\link{lda}} is supported now. \code{lda} decide the class 
#'  (e.g. species) based on the calculated posterior probability. 
#'  New sample is assigned to the class that has the higher posterior probability. 
#'  We can reject the classification result if the posterior probability is low 
#'  (the "unsure" classification). In this way we can increase the overall 
#'  accuracy (correct predictions out of all the reported results), 
#'  but the tradeoff is the decrease in number of reported classification result.
#'  Currently the function test the threshold value of 0.50 0.60 0.65 0.70 0.75 
#'  0.80 0.85 0.90 0.95 0.99. Since thresholding is just an option, 
#'  the choice of threshold value is really arbitrary and depends on one's objectives.
#' @param X matrix of PC scores from \code{\link{rGPA}} object (the \code{score} value)
#' @param Y a factor giving the grouping, e.g. the \code{sp} value from 
#'  \code{\link{routine1}} object or obtain from \code{\link{getclass}}
#' @param saveplot logical. The plot will be saved if \code{TRUE}, 
#'  and no plot is displayed in the window. Note that outliers of
#'   boxplot are not plotted.
#' @param plotsize numeric. Plot size for the plot saved into file 
#'  (used only when \code{saveplot=TRUE}), unit in pixel.
#' @param run argument to be passed to \code{\link{mrkfcv}}
#' @param k argument to be passed to \code{\link{mrkfcv}}
#' @return a matrix and a plot are returned, giving the overall accuracy and 
#'  the total reported prediction percentage over a range of threshold values. 
#' @seealso 
#'  Similar: \code{\link{pccv}}, \code{\link{harcv}}
#'  
#'  Which this function wraps: \code{\link{mrkfcv}}
#' @export
#' @references
#' Beleites, C., & Salzer, R. (2008). Assessing and improving the stability of 
#'  chemometric models in small sample size situations. 
#'  \emph{Analytical and Bioanalytical Chemistry}, 390(5), 1261-1271.

threcv <- function(X, Y, saveplot=FALSE, plotsize=1000, run=30, k=5) {
  accu <- numeric()
  accu.sd <- numeric()
  total <- numeric()
  total.sd <- numeric()
  temp <- NULL
  result <- data.frame(matrix(NA, k * run * 10, 3))
  levellabel <- c(0.5, seq(0.6, 0.95, 0.05), 0.99)
  threlevel <- gl(10, run * k, labels=levellabel, ordered=TRUE) # group level
  result[, 3] <- threlevel
  for (i in 1:10) {
    cat ("\r                       (Evaluating threshold at: ", levellabel[i], 
         ") | threcv progress: [", round(i / 10 * 100), "%]       ", sep="")
    temp <- mrkfcv(X=X, Y= Y, suppress="text", k=k, run=run, 
            threshold=levellabel[i])
    result[((i - 1) * run * k + 1):(i * run * k), 1] <- 100 - temp$misclass
    result[((i - 1) * run * k + 1):(i * run * k), 2] <- temp$total.pred
    accu[i] <- temp$accuracy
    accu.sd[i] <- temp$accu.sd
    total[i] <- temp$total 
    total.sd[i] <- temp$total.sd
    #setTxtProgressBar(probar, i)
  }
  if (saveplot == TRUE) {
    filename <- "threshold-optimization.tif"
    tiff(filename, plotsize, plotsize, res=172, compression = "lzw")
  }  
  boxplot(result[, 1] ~ result[, 3], outline=FALSE, 
          ylab="Percentage", border=3, ylim=c(min(result[, 2]), 100), 
          at=(1:10) - 0.15, boxwex = 0.2, axes=FALSE,
          xlab="Posterior probability threshold value")
  legend("bottomleft", fill=c(3, 4), legend= c("Cross-validated accuracy",
         "Proportion of data predicted"), bty='n')
  axis(2, at=seq(0, 100, 10), las=2)
  axis(1, at=1:10, labels=levellabel)
  boxplot(result[, 2] ~ result[, 3], border=4, outline=F, 
          add=TRUE, at=(1:10) + 0.15, boxwex = 0.2, axes=FALSE)
  box()
  if (saveplot == TRUE) {
    dev.off()
    cat("\n\nThe plot is saved at:", 
        paste(getwd(), filename, sep="/"), "\n")
  }
  result <- data.frame(cbind(accu, accu.sd, total, total.sd))
  result$threshold <- levellabel
  cat("\nEvaluation completed\n")
  return(result)
}