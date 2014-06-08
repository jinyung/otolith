\name{kfcv / kfcv2}
\alias{kfcv}
\alias{kfcv2}
\title{k-fold cross-validation}
\description{
  Run k-fold cross validation to validate a classification model
}
\usage{
kfcv(X, Y, method = c("lda", "plsda", "tree"), k = 5, threshold, ncomp)
kfcv2(X, Y, method = c("lda", "plsda", "tree"), k = 5, threshold, ncomp)
}
\arguments{
  \item{X}{matrix/ dataframe of predictors, e.g. EFA
  coefficients/ PC scores selected using \code{\link{selectdim}}}

  \item{Y}{vector giving the class, e.g. value obtained from 
  \code{\link{getclass}} or \code{sp} value
  from \code{\link{routine1}} object}

  \item{method}{method \code{"\link[MASS]{lda}"} for linear
  discriminant analysis, \code{"\link[tree]{tree}"} for
  classification tree, \code{"\link[mixOmics]{plsda}"} for partial
  least square-linear discriminant analysis}

  \item{k}{fold number of cross-validation}

  \item{threshold}{optional. A numeric value between 0-1 to
  set the threshold of posterior probility. Any class
  prediction with posterior probility lower than this value
  will be \code{NA}-ed and not reported. See \code{\link{threcv}}}

  \item{ncomp}{argument passed to \code{\link[mixOmics]{plsda}}.Used
  only when \code{method="plsda"}.}
}
\value{
  \item{misclass}{vector of k values of misclassification
  rate in percent resulted from each fold of testing}
    
  \item{total}{total
  number of prediction after excluding the ones lower than
  threshold, if \code{threshold} value is given}
  
  \item{stat}{[\code{kfcv2} only] \code{k} number of matrix containing the
  calculated precision, sensitivity(recall) and specificity
  for each class, for each fold. May contain \code{NA}
  values if the class is not present in the fold.}
  
  \item{conmat}{[\code{kfcv2} only] \code{k} number of confusion matrices, shown
  as proportion rather than counts. Proportion = number
  correctly or incorrectly predicted divided by the total
  number of that class in training set.} 
}
\details{

  \code{kfcv2} is an improved version of \code{kfcv} function, 
  with the ability to calculate the by-class statistics 
  (recall, precision and specificity). It also calculate confusion matrices for 
  each folds. Because of the extra calculations, the 
  execution takes a little bit more time, and hence user might prefer \code{kfcv} 
  in certain circumstances when speed is a concern.

  What the stat values mean:
  \itemize{ 
    \item Recall = Sensitivity = tp / (tp + fn) 
    \item Precision = tp / (tp + fp) 
    \item Specificity = tn / (tn + fp) 
  } and, tp= true positive, tn=true negative, fp=false
  positive, fn=false negative Please refer to reference for
  detailed explanation.
  
}
\references{
  Sokolova, M., & Lapalme, G. (2009). A systematic analysis
  of performance measures for classification tasks.
  \emph{Information Processing & Management}, 45(4), 427-437.  
}
\seealso{
  Function that wraps this function: \code{\link{mrkfcv}}, \code{\link{mrkfcv2}}
}