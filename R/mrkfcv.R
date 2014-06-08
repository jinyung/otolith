#' Multiple-run k-fold cross-validation
#' 
#' @description run multiple runs of k-fold cross validation, see referece. 
#'  "Use all data" variant is implemented here.
#' @param run number of run to be used in multiple runs of k-fold cross-validation
#' @param suppress suppress the running status in R console when \code{TRUE}. the 
#'  option \code{"text"} is used for \code{\link{pccv}/\link{harcv}}wrappers.
#' @inheritParams kfcv
#' @return 
#'  \item{accuracy}{cross-validated accuracy for the tested classifier, 
#'    resulted from the average of \code{k x ru}n numbers of accuracy generated 
#'    by the function}
#'  \item{accu.sd}{standard deviation for the accuracy, calculated from 
#'    the \code{k x run} number of results}  
#'  \item{total}{mean total successful prediction in percent, no returned if 
#'    \code{threshold} value is not given}
#'  \item{total.sd}{sd of total successful prediction in percent, not returned if 
#'    \code{threshold} value is not given}
#'  \item{misclass}{vector of \code{run x k} number of misclassification rate}
#'  \item{total.pred}{raw prediction success for \code{run x k} number 
#'    of submodels, give \code{NULL} if \code{threshold} is not given}
#' @seealso 
#'  Which this function wraps: \code{\link{kfcv}}
#'  
#'  Function that wraps this function: \code{\link{pccv}}, \code{\link{harcv}}, 
#'    \code{\link{threcv}} 
#' @export
#' @references
#'  Bouckaert, R.R., 2003. Choosing between two learning algorithms based on 
#'  calibrated tests. In: Fawcett, T., Mishra, N. (Eds.), 
#'  \emph{Proceedings of the Twentieth International Conference (ICML 2003) on 
#'  Machine Learning}. August 21-24, 2003, AAAI Press, Washington.
#'  
#'  Beleites, C., Baumgartner, R., Bowman, C., Somorjai, R., Steiner, G., 
#'  Salzer, R., & Sowa, M. G. (2005). Variance reduction in estimating 
#'  classification error using sparse datasets. 
#'  \emph{Chemometrics and Intelligent Laboratory Systems}, 79(1), 91-100.

mrkfcv <- function (X, Y, method=c("lda", "tree", "plsda"), k=5, run=100, 
                    threshold, ncomp, suppress=c(FALSE, TRUE, "text")) {
  method <- match.arg(method)
  suppress <- match.arg(suppress)
  misclass <- numeric()
  if (!missing(threshold))
    total.pred <- NULL
  if (suppress == FALSE) {
    cat("\n **Running ", run, "-runs of ", k, "-fold cross-validation:\n\n", sep="")
    pb <- txtProgressBar(1, run, style=3, char="|")
  }
  for (p in 1:run) {
    run.p <- kfcv(X=X, Y=Y, method=method, k=k, threshold=threshold, ncomp=ncomp)
    misclass[(p * k - (k - 1)):(p * k)] <- run.p$misclass
    if (!missing(threshold))
      total.pred[(p * k - (k - 1)):(p * k)] <- run.p$total
    if (suppress == FALSE)
      setTxtProgressBar(pb, p)
    else if (suppress == "text")
      cat("\rmrkfcv progress: [", round(p / run * 100), "%]", sep="")
  }
  accuracy <- round(100 - mean(misclass), 2)
  accu.sd <- round(sd(misclass), 2)
  if (!missing(threshold)) {
    total <- round(mean(total.pred), 2)
    total.sd <- round(sd(total.pred), 2)
    return.list <- list(accuracy=accuracy, accu.sd=accu.sd, total=total, total.sd=total.sd, 
                   misclass=misclass, total.pred=total.pred)
  } else {
    return.list <- list(accuracy=accuracy, accu.sd=accu.sd, misclass=misclass)
  }
  if (suppress == FALSE)
    cat("\n\n")
  return (return.list)
}