#' k-fold cross-validation (version-2)
#' 
#' @description improved \code{\link{kfcv}} function, with the ability to 
#'  calculate the by-class statistics (recall, precision and specificity)  
#' @inheritParams kfcv
#' @importFrom MASS lda
#' @importFrom mixOmics plsda
#' @importFrom tree tree
#' @return
#'  \item{misclass}{vector of k values of misclassification rate in percent 
#'    resulted from each fold of testing}
#'  \item{stat}{\code{k} number of matrix containing the calculated precision, 
#'    sensitivity(recall) and specificity for each class, for each fold. 
#'    May contain \code{NA} values if the class is not present in the fold.}  
#'  \item{conmat}{\code{k} number of confusion matrices, shown as proportion rather 
#'    than counts. Proportion = number correctly or incorrectly predicted divided 
#'    by the total number of that class in training set.}
#'  \item{total}{total number of prediction after excluding the ones below 
#'  threshold, if \code{threhold != NULL}}
#' @details
#'  \itemize{
#'    \item Recall = Sensitivity = tp / (tp + fn)
#'    \item Precision = tp / (tp + fp)
#'    \item Specificity = tn / (tn + fp) 
#'  }
#'  and, tp= true positive, tn=true negative, fp=false positive, fn=false negative
#'  Please refer to reference for detailed explanation.
#' @seealso
#'  Similar: \code{\link{kfcv}}
#'  
#'  Function that wraps this function: \code{\link{mrkfcv2}}
#' @export
#' @references 
#'  Sokolova, M., & Lapalme, G. (2009). A systematic analysis of performance 
#'  measures for classification tasks. \emph{Information Processing & Management}, 
#'  45(4), 427-437.

kfcv2 <- function(X, Y, method=c("lda", "plsda", "tree"),
                  k=5, threshold, ncomp) {
  # create the folds
  alltestingindices <- kfcv.testing(dim(X)[1], k=k)
  # some processing
  method <- match.arg(method)
  Y <- factor(Y)
  dat <- data.frame(X)
  dat$predictor <- as.matrix(X) # matrix can be saved as a variable in dataframe, 
                                # easier to be called later
  dat$class <- Y
  class.level <- levels(Y)
  class.length <- length(class.level)
  # initialize variables 
  misclass <- numeric()
  total.pred <- NULL
  stat <- array (data=NA, dim= c(class.length, 3, k), 
          dimnames=list(class.level, c("recall", "precision", 
          "specificity"), paste0("fold", 1:k)))
  tempstat <- matrix(data=NA, nrow = class.length, ncol= 3, 
              dimnames= list(class.level, c("recall", "precision", 
              "specificity")))
  conmat <- array(data=NA, dim= c(class.length, class.length, k), 
            dimnames=list(class.level, class.level, paste0("fold", 1:k)))
  # iterate through folds
  for (i in 1:k) {
    testingindices <- alltestingindices[[i]]
    train <- dat[-testingindices,]
    test <- dat[testingindices,]
    testlength<- length(levels(test$class))
    trainlength<- length(levels(train$class))
    # build classification sub-models with selected method
    if (method == "lda"){
      require(MASS);
      mod.i <- lda(class ~ predictor, data=train, prior=
                     rep(1/trainlength, trainlength))
      prediction.temp <- predict(mod.i, test)
      prediction.i <- prediction.temp$class
      # add in posterior thresholding 25.4.14
      if (!missing(threshold)) {
        posterior <- prediction.temp$posterior
        posmax <- apply(posterior, 1, max)
        prediction.i[which (posmax < threshold)] <- NA 
          # class predicted with posterior < thereshold are NA-ed.  
        total.pred.i <- (dim(test)[1] - sum(posmax < threshold)) /
                        dim(test)[1] * 100 
                        # total prediction (= percent of total - NA)
      }
    } else if (method == "plsda") {
      require(mixOmics)
      if(missing(ncomp))
        ncomp <- trainlength - 1
      mod.i <- plsda(X=train$predictor, Y=train$class, ncomp=ncomp)
      prediction.temp <- predict(mod.i, test$predictor)
      prediction.i<- levels(train$class)[prediction.temp$class$max.dist[, ncomp]]
      prediction.i <- factor(prediction.i, levels=levels(Y))
    } else if (method == "tree"){
      require(tree)
      mod.i <- tree(class ~ predictor, data=train)
      prediction.i <- predict(mod.i, test, type='class')
    }
    # get the misclassification rate
    t.i <- table(prediction.i, test$class)
    tempcon <- conmat[, , 1]
    t.i <- as.data.frame.matrix(t.i)
    wrongsum.i <- 0
    for (j in 1:testlength) {
      wrongsum.i <- wrongsum.i+ sum(t.i[rownames(t.i) != colnames(t.i)[j], 
                    colnames(t.i)[j]])
    }
    misclass.rate.i <- wrongsum.i / length(alltestingindices[[i]]) * 100
    misclass[i] <- misclass.rate.i 
    if (!missing(threshold))
      total.pred[i] <- total.pred.i    
    tempstat.i <- tempstat
    # calculate the stats
    for (m in 1: class.length) {
      tempstat.i[m, 1] <- t.i[m, m] / (sum(t.i[, m])) # recall
      if (sum(t.i[m, ]) > 0)
        tempstat.i[m, 2] <- t.i[m, m] / (sum(t.i[m,])) # precision
      tempstat.i[m, 3] <- sum(t.i[-m, -m])/ (sum(t.i[-m, -m]) + 
                         sum(t.i[m, ])-t.i[m, m] ) # specificity
      for (r in 1:class.length)
        tempcon[r, m] <- t.i[r, m]/sum(t.i[, m])   
    }
    # remove the stat row /cm column with missing class level in training/
    # test sets
    removeindex1.i <- which(summary(test$class) == 0)
    removeindex2.i <- which(summary(train$class) == 0)
    tempstat.i[removeindex1.i, ] <- NA
    tempstat.i[removeindex2.i, ] <- NA
    tempcon[, removeindex1.i] <- NA
    tempcon[, removeindex2.i] <- NA
    stat[, , i] <- tempstat.i  
    conmat[, , i] <- tempcon
  }
  invisible(list(misclass=misclass, stat=stat, conmat=conmat, total=total.pred))
}