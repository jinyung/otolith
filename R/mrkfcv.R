#' Multiple-run k-fold cross-validation
#' 
#' @description run multiple runs of k-fold cross validation, see referece. 
#'  "Use all data" variant is implemented here.
#' @param run number of run to be used in multiple runs of k-fold cross-validation
#' @param suppress suppress the running status in R console
#' @inheritParams kfcv
#' @return 
#'  \item{accuracy}{cross-validated accuracy for the tested classifier, 
#'    resulted from the average of \code{k x ru}n numbers of accuracy generated 
#'    by the function}
#'  \item{accu.sd}{standard deviation for the accuracy, calculated from 
#'    the \code{k x run} number of results}  
#'  \item{total}{mean total successful prediction in percent, give \code{NA} if 
#'    \code{threshold} is not given}
#'  \item{total.sd}{sd of total successful prediction in percent, give \code{NA} if 
#'    \code{threshold} is not given}
#'  \item{misclass}{vector of \code{run x k} number of misclassification rate}
#'  \item{total.pred}{raw prediction success for \code{run x k} number 
#'    of submodels, give \code{NULL} if \code{threshold} is not given}
#' @seealso 
#'  Similar: \code{\link{mrkfcv2}}
#'  
#'  Which this function wraps: \code{\link{kfcv}}
#'  
#'  Function that wraps this function: \code{\link{pccv}}, \code{\link{harcv}}, 
#'    \code{\link{threcv}} 
#' @export
#' @references 
#'  Bouckaert, R.R., 2003. Choosing between two learning algorithms based on 
#'  calibrated tests. In: Fawcett, T., Mishra, N. (Eds.), 
#'  Proceedings of the Twentieth International Conference (ICML 2003) on 
#'  Machine Learning. August 21-24, 2003, AAAI Press, Washington.
#'  
#'  Beleites, C., Baumgartner, R., Bowman, C., Somorjai, R., Steiner, G., 
#'  Salzer, R., & Sowa, M. G. (2005). Variance reduction in estimating 
#'  classification error using sparse datasets. 
#'  Chemometrics and Intelligent Laboratory Systems, 79(1), 91-100.

mrkfcv <- function (X, Y, method="lda", k=5, threshold= NULL, 
                    ncomp=NULL, run=100, suppress=FALSE) {
  misclass <- numeric()
  total.pred <- NULL
  if (suppress == FALSE){
    cat("\n **Running ", run, "-runs of ", k, "-fold cross-validation:\n\n", sep="")
    pb <- txtProgressBar(1, run, style=3, char="|")
  }
  for (p in 1:run){
    run.p <- kfcv(X=X, Y=Y, method=method, k=k, threshold=threshold, ncomp=ncomp)
    misclass[(p * k - (k - 1)):(p * k)] <- run.p$misclass
    total.pred[(p * k - (k - 1)):(p * k)] <- run.p$total
    if (suppress == FALSE)
      setTxtProgressBar(pb, p)      
  }
  accuracy <- round(100 - mean(misclass), 2)
  accu.sd <- round(sd(misclass), 2)
  total <- round(mean(total.pred), 2)
  total.sd <- round(sd(total.pred), 2)
  if (suppress == FALSE)
    cat("\n\n")
  return (list(accuracy=accuracy, accu.sd=accu.sd, total=total, total.sd=total.sd, 
          misclass=misclass, total.pred=total.pred))
}