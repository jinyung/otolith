#' LDA model aggregation 
#' 
#' @description aggregate LDA models based on iterated k-fold resampling method 
#' @param X matrix/ dataframe of predictors, e.g. EFA coefficients/ PC scores 
#'  selected using \code{\link{selectdim}}
#' @param Y factor/ character giving the class, e.g. value obtained from 
#'   \code{\link{getclass}} or \code{sp} value from \code{\link{routine1}}
#'   object
#' @param newdata matrix/ dataframe of newdata to be predicted. see details.
#' @param type prediction based on the majority vote (\code{"vote"}) from the 
#'  submodels or based on the mean posterior probability from the submodels 
#'  (\code{"post"})
#' @param k the fold number for k-fold resampling
#' @param run the iteration number for iteration of k-fold resampling
#' @param threshold single numeric value of < 1. threshold of the proportion of 
#'   majority vote or the mean posterior probability. predictions with less than
#'   the value will be reported as \code{NA}
#' @param prior the prior used in \code{\link[MASS]{lda}} models, \code{"equal"}
#'   means all classes have same prior, \code{"proportion"} means prior
#'   according to the classes weight in the training data.
#' @param suppress logical. whether to suppress the progress monitoring output
#' @details 
#'   If \code{newdata} is provided, the function is in the \emph{prediction
#'   mode}, the aggregated model will be built from \code{X} and \code{Y} and 
#'   predicition is performed on \code{newdata}. Otherwise, if \code{newdata =
#'   NULL} (default) the function is in \emph{evaluation mode}.
#'   
#'   In \emph{evaluation mode}, overall accuracy of the model and the by-class
#'   statistics are calculated, similar to that of \code{\link{mrkfcv2}}. 
#'   However, the statistics are calculated based on the aggregated prediction.
#'   See reference for explanation on model aggregation and thresholding.  
#' @return
#'  \item{accuracy}{[\emph{evaluation mode}] the overall accuracy in percent}
#'  \item{conmat}{[\emph{evaluation mode}] confusion matrix}
#'  \item{stat}{[\emph{evaluation mode}] matrix containing the statistics of
#'    each class, see details}
#'  \item{total}{[\emph{evaluation mode}] the total percent of reported
#'   prediction after threshold. give \code{NULL} if \code{threshold} is not
#'   given}
#'  \item{ind.prediction}{matrix containing the prediction result on each 
#'    training/ new specimens}  
#' @importFrom MASS lda
#' @seealso \code{\link{kfcv}}, \code{\link{mrkfcv}}
#' @references 
#'  Beleites, C., & Salzer, R. (2008). Assessing and improving the stability of
#'  chemometric models in small sample size situations. \emph{Analytical and
#'  Bioanalytical Chemistry}, 390(5), 1261-1271.
#' @export

agglda <- function (X, Y, newdata = NULL, type = c("vote", "post"), 
                    k = 5, run = 100, threshold, 
                    prior = c("equal", "proportion"), suppress = FALSE) {
  # initialize
  type <- match.arg(type)
  prior <- match.arg(prior)
  Y <- factor(Y)
  class.length <- length(levels(Y))
  if (is.null(newdata)) {
    agg.class.all <- as.data.frame(matrix(data = NA, run, dim(X)[1], 
                                          dimnames = list(paste0("run", 1:run), 
                                                          rownames(X))))
    agg.posterior.all <- array(data = NA, dim = c(dim(X)[1], class.length, 
                         run), dimnames = list(rownames(X), levels(Y), 
                         paste0("run-", 1:run)))
  } else {
    krun.dimname <- paste(rep(paste0("run", 1:run), each = k),
                          paste0("fold", 1:k), sep="-")
    agg.class.all <- as.data.frame(matrix(data = NA, run * k, dim(newdata)[1], 
                                          dimnames = list(krun.dimname, 
                                                          rownames(newdata))))
    agg.posterior.all <- array(data = NA, dim = c(dim(newdata)[1], class.length, 
                         run * k), dimnames = list(rownames(newdata), levels(Y), 
                                                 krun.dimname))
  }
  # progress bar
  if (suppress == FALSE) {
    cat("\n **Running ", run, "iterations ")
    pb <- txtProgressBar(1, run, style=3, char="|")
  }
  # going thru iterations
  for (i in 1:run) {
    run.i <- .kflda(X = X, Y = Y, k = k, newdata = newdata, 
                    prior = prior)
    if (is.null(newdata)) {
      agg.class.all[i, ] <- run.i$agg.class
      agg.posterior.all[, , i] <- run.i$agg.posterior
    } else {
      subindex <- (((i - 1) * k) + 1):(k * i)
      agg.class.all[subindex, ] <- run.i$agg.class
      agg.posterior.all[, , subindex] <- run.i$agg.posterior
    }
    if (suppress == FALSE)
      setTxtProgressBar(pb, i)
  }
  # table vote counts/ mean posterior of each candidates across runs
  vprediction <- switch(type, 
                        vote = vapply(agg.class.all, function (x) 
                          summary(factor(x, levels=levels(Y))), 
                          integer(class.length)), 
                          # this return a matrix of summary (summary on each col)
                        post = t(apply(agg.posterior.all, c(1, 2), 
                                        mean, na.rm = TRUE)))
  # calculate sd for mean posterior too
  if (type == "post")
    posterior.sd <- t(apply(agg.posterior.all, c(1, 2), 
                            sd, na.rm = TRUE))
  # announce the winner
  prediction.class <- levels(Y)[apply(vprediction, 2, which.max)]
  # highlight the correct candidate/ winner's result
  if (is.null(newdata)) { # get the vote %/ mean posterior of correct candidate
    prediction.stat <- mapply(function(x, y) x[which(levels(y) == y)], 
                              data.frame(vprediction), Y)
    if (type == "post")
      posterior.sd <- mapply(function(x, y) x[which(levels(y) == y)], 
                             data.frame(posterior.sd), Y)
  } else { # get whatever the highest (winner) since dont know which is correct
    prediction.stat <- apply(vprediction, 2, max)
    if (type == "post")
      posterior.sd <- mapply(function(x, y) x[y], data.frame(posterior.sd), 
                               apply(vprediction, 2, which.max))
  }
  # then get the proportion (for posterior, value did not change)
  prediction.percent <- prediction.stat / apply(vprediction, 2, sum)
  # thresholding based on the proportion
  if (!missing(threshold)) {
    prediction.class[prediction.percent < threshold] <- NA
    pred.length <- length(prediction.class[!is.na(prediction.class)])
    total.pred <- round(pred.length / dim(X)[1] * 100, 2)
  } else {
    total.pred <- NULL
    pred.length <- dim(X)[1]
  }
  if (is.null(newdata)) {
    ind.prediction <- prediction.class == Y
    accuracy <- sum(ind.prediction, na.rm = TRUE) / pred.length * 100
    result <- data.frame(ID = dimnames(X)[[1]],
                             species = Y, predict = prediction.class, 
                             correct = ind.prediction)
    # confusion matrix and stats
    cm <- as.matrix(as.data.frame.matrix(table(prediction.class, Y)))
    TP <- diag(cm)
    TN <- sapply(c(1:class.length), function(x) sum(cm[ -x, -x], na.rm = TRUE))
    FP <- rowSums(cm, na.rm = TRUE) - TP
    FN <- colSums(cm, na.rm = TRUE) - TP
    recall <- TP / (TP + FN) # recall
    precision <- TP / (TP + FP) # precision
    specificity <- TN / (TN + FP) # specificity
    prevalence <- colSums(cm) / sum(cm, na.rm = TRUE)
    PosPredVal <- (recall * prevalence) / ((recall * prevalence) + 
                    ((1 - specificity) * (1 - prevalence)))
    NegPredVal <- (specificity * (1 - prevalence)) / 
                    (((1 - recall) * prevalence) + 
                    ((specificity) * (1 - prevalence)))
    stat <- round(cbind(recall, precision, specificity, 
                        PosPredVal, NegPredVal), 2)    
    # different % name for different method
    switch(type, vote = {result$vote.percent <- round(prediction.percent, 3)},
           post = {result$mean.posterior <- round(prediction.percent, 3)
                   result$sd.posterior <- round(posterior.sd, 3)})
    result <- list(accuracy = round(accuracy, 2), conmat = cm, stat = stat, 
                   total = total.pred, ind.prediction = result)
  } else {
    result <- data.frame(predict = prediction.class)
    switch(type, vote = {result$vote.percent <- round(prediction.percent, 3)},
           post = {result$mean.posterior <- round(prediction.percent, 3)
                   result$sd.posterior <- round(posterior.sd, 3)})
    rownames(result) <- rownames(newdata)
  }
  if (suppress == FALSE)
    cat("\n\n")
  return(result)
}