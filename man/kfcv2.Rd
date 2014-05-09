\name{kfcv2}
\alias{kfcv2}
\title{k-fold cross-validation (version-2)}
\usage{
kfcv2(X, Y, method = c("lda", "plsda", "tree"), k = 5, threshold = NULL,
  ncomp = NULL)
}
\arguments{
  \item{X}{matrix/ dataframe of predictors, e.g. EFA
  coefficients/ PC scores selected using
  \code{\link{selectdim}}}

  \item{Y}{vector giving the class, e.g. \code{sp} value
  from \code{\link{routine1}} object}

  \item{method}{method \code{"\link{lda}"} for linear
  discriminant analysis, \code{"\link{tree}"} for
  classification tree, \code{"\link{plsda}"} for partial
  least square-linear discriminant analysis}

  \item{k}{fold number of cross-validation}

  \item{threshold}{optional. A numeric value between 0-1 to
  set the threshold of posterior probility. Any class
  prediction with posterior probility lower than this value
  will be \code{NA}-ed and not reported.}

  \item{ncomp}{argument passed to \code{\link{plsda}}.Used
  only when \code{method="plsda"}.}
}
\value{
\item{misclass}{vector of k values of misclassification
rate in percent resulted from each fold of testing}
\item{stat}{\code{k} number of matrix containing the
calculated precision, sensitivity(recall) and specificity
for each class, for each fold.  May contain \code{NA}
values if the class is not present in the fold.}
\item{conmat}{\code{k} number of confusion matrices, shown
as proportion rather than counts. Proportion = number
correctly or incorrectly predicted divided by the total
number of that class in training set.} \item{total}{total
number of prediction after excluding the ones below
threshold, if \code{threhold != NULL}}
}
\description{
improved \code{\link{kfcv}} function, with the ability to
calculate the by-class statistics (recall, precision and
specificity)
}
\details{
Recall = Sensitivity = tp / (tp + fn) Precision = tp / (tp
+ fp) Specificity = tn / (tn + fp) and, tp= true positive,
tn=true negative, fp=false positive, fn=false negative
Please refer to reference for detailed explanation.
}
\references{
Sokolova, M., & Lapalme, G. (2009). A systematic analysis
of performance measures for classification tasks.
Information Processing & Management, 45(4), 427-437.
}
\seealso{
Similar: \code{\link{kfcv}}

Function that wraps this function: \code{\link{mrkfcv2}}
}

