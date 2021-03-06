#' k-fold cross-validation 
#' 
#' @description Run k-fold cross validation to validate a classification model
#' @param X matrix/ dataframe of predictors, e.g. EFA coefficients/ PC scores 
#'  selected using \code{\link{selectdim}}
#' @param Y vector giving the class, e.g. \code{sp} value from 
#'  \code{\link{routine1}} object
#' @param method method \code{"\link{lda}"} for linear discriminant analysis,
#'   \code{"\link{tree}"} for classification tree, 
#' @param k fold number of cross-validation
#' @param threshold optional. A numeric value between 0-1 to set the threshold of 
#'  posterior probility. Any class prediction with posterior probility lower than this 
#'  value will be \code{NA}-ed and not reported.  
#' @param ncomp argument passed to \code{\link[mixOmics]{plsda}}.Used only 
#'  when \code{method="plsda"}.
#' @return 
#'  \item{misclass}{vector of k values of misclassification rate in percent 
#'    resulted from each fold of testing} 
#'  \item{total}{total number of prediction after excluding the ones lower than 
#'  threshold, if \code{threshold} value is given}
#' @seealso 
#'  Similar: \code{\link{kfcv2}}
#'  
#'  Function that wraps this function: \code{\link{mrkfcv}}
#'  
#'  Methods of classifier: \code{\link{lda}}, \code{\link{plsda}}, \code{\link{tree}}
#' @importFrom MASS lda
#' @importFrom tree tree
#' @export

kfcv <- function(X, Y, method=c("lda", "plsda", "tree"), 
                 k=5, threshold, ncomp) {
  alltestingindices <- kfcv.testing(dim(X)[1], k=k)
  method <- match.arg(method)
  Y <- factor(Y)
  dat <- data.frame(X)
  dat$predictor <- as.matrix(X)
  dat$class <- Y
  misclass <- numeric()
  total.pred <- NULL
  ind.prediction <- logical(length(Y)) # new addition of by subject success
  for (i in 1:k) {
    testingindices <- alltestingindices[[i]]
    train <- dat[-testingindices, ]
    test <- dat[testingindices, ]
    testlength <- length(levels(test$class)) # no more important as droplevels() not use here, retained anyway.      
    trainlength <- length(levels(train$class))
    if (method == "lda") {
      # require(MASS)
      mod.i<- MASS::lda(class ~ predictor, data=train,
                        prior= rep(1 / trainlength, trainlength))
      prediction.temp <- predict(mod.i, test)
      prediction.i <- prediction.temp$class
      # add in posterior thresholding 25.4.14
      if (!missing(threshold)) {
        posterior <- prediction.temp$posterior
        posmax <- apply(posterior, 1, max)
        #class predicted with posterior < thereshold are NA-ed.  
        prediction.i[which(posmax < threshold)] <- NA 
        total.pred.i <- (dim(test)[1] - sum(posmax < threshold)) / 
                        dim(test)[1] * 100 
                        # total prediction (= percent of total - NA)
      }
    } else if (method == "plsda") {
      # require(mixOmics)
      # if (missing(ncomp)) 
      #   ncomp<- trainlength-1
      # mod.i <- mixOmics::plsda(X=train$predictor, Y=train$class, ncomp=ncomp)
      # prediction.temp <- predict(mod.i, test$predictor)
      # prediction.i <- levels(train$class)[prediction.temp$class$max.dist[, ncomp]]
      # prediction.i <- factor(prediction.i, levels=levels(Y))
    } else if (method == "tree") {
      # require(tree)
      mod.i <- tree::tree(class ~ predictor, data=train)
      prediction.i <- predict(mod.i, test, type='class')
    }
    t.i <- table(prediction.i, test$class)
    t.i <- as.data.frame.matrix(t.i)
    ind.prediction[testingindices] <- prediction.i == test$class
    wrongsum.i <- 0
    for (j in 1:testlength) {
      wrongsum.i <- wrongsum.i +  
        sum(t.i[rownames(t.i) != colnames(t.i)[j], colnames(t.i)[j]])
    }
    misclass.rate.i <- wrongsum.i / length(alltestingindices[[i]]) * 100
    misclass[i] <- misclass.rate.i
    if (!missing(threshold))
      total.pred[i] <- total.pred.i      
  }
  return(list(misclass=misclass, total=total.pred, ind.prediction = ind.prediction))
}