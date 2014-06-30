#' Multiple-run k-fold cross-validation (version-2)
#' 
#' @description modified version of \code{\link{mrkfcv}}, comes with calculation 
#'  of recall, precision, and specificity.
#' @details The calculated by-class statistics (\code{stat.sum}) are average of all values of 
#'  number of \code{k x run} of submodels (\code{NA} values are excluded). 
#' @inheritParams mrkfcv
#' @return 
#'  \item{accuracy}{cross-validated accuracy for the tested classifier, 
#'    resulted from the average of \code{k x run} numbers of accuracy generated 
#'    by the function}
#'  \item{accu.sd}{standard deviation for the accuracy, calculated from 
#'    the \code{k x run} number of results}  
#'  \item{stat.sum}{cross-validated by-class precision, recall and specificity}
#'  \item{conmat}{confusion matrix shown in proportion, average across all 
#'    confusion matrices of \code{k x run} number of submodels. Proportion = 
#'    number correctly or incorrectly predicted divided by the total number 
#'    of that class in training set.}
#'  \item{total}{mean total successful prediction in percent, give \code{NA} if 
#'    \code{threshold} is not given}
#'  \item{total.sd}{sd of total successful prediction in percent, give \code{NA} if 
#'    \code{threshold} is not given}
#'    
#' @seealso 
#'  Similar: \code{\link{mrkfcv}}
#'  
#'  Which this function wraps: \code{\link{kfcv2}}
#' @export

mrkfcv2 <- function(X, Y, method=c("lda", "tree", "plsda"), k=5, run=100, 
                    threshold, ncomp, suppress=FALSE) {
  method <- match.arg(method)
  misclass <- numeric()
  total.pred <- NULL
  class.level <- levels(Y)
  class.length <- length(class.level)
  stat.array <- array (data=NA, dim=c(class.length, 3, run * k), 
                dimnames= list(class.level, c("recall", "precision", 
                "specificity"), NULL))
  stat.sum <- matrix(data=NA, nrow = class.length, ncol= 3)
  stat.sum.sd <- stat.sum
  ConMat <- array (data=NA, dim=c(class.length, class.length, run * k), 
            dimnames= list(class.level, class.level, NULL))
  #progress bar
  if (!suppress) {
    cat("\n **Running ", run, "-runs of ", k, "-fold cross-validation:\n\n", 
        sep="")
    pb <- txtProgressBar(1, run, style=3, char="|")  
  }
  for (p in 1:run) {
    run.p <- kfcv2(X=X, Y=Y, method=method, k=k, threshold=threshold, ncomp=ncomp)
    misclass[(p * k - (k - 1)):(p * k)] <- run.p$misclass
    total.pred[(p * k-(k-1)):(p * k)] <- run.p$total    
    stat.array[, , (p * k - (k - 1)):(p * k)] <- run.p$stat
    ConMat[, , (p * k - (k - 1)):(p * k)] <- run.p$conmat
    if (!suppress)
      setTxtProgressBar(pb, p)
  }
  stat.sum <- round(apply(stat.array, c(1, 2), mean, na.rm=TRUE), 2)
  stat.sum.sd <- round(apply(stat.array, c(1, 2), sd, na.rm=TRUE), 2)
  conmatF <- ConMat[, , 1]
  conmatF <- round(apply(ConMat, c(1, 2), mean, na.rm=TRUE), 2)
  accuracy <- round(100 - mean(misclass), 2)
  accu.sd <- round(sd(misclass), 2)
  if (!missing(threshold)) {
    total <- round(mean(total.pred), 2)
    total.sd <- round(sd(total.pred), 2)
    return.list <- list(accuracy=accuracy, accu.sd=accu.sd, stat.sum=stat.sum, 
                        conmat=conmatF, total=total, total.sd=total.sd)
  } else {
    return.list <- list(accuracy=accuracy, accu.sd=accu.sd, stat.sum=stat.sum, 
                        conmat=conmatF)
  }  
  cat("\n\n")
  return (return.list)
}