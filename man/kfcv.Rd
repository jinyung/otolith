\name{kfcv}
\alias{kfcv}
\title{k-fold cross-validation}
\usage{
kfcv(X, Y, method = c("lda", "plsda", "tree"), k = 5, threshold = NULL,
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
\item{total}{total number of prediction after excluding the
ones below threshold, if \code{threhold != NULL}}
}
\description{
Run k-fold cross validation to validate a classification
model
}
\seealso{
Similar: \code{\link{kfcv2}}

Function that wraps this function: \code{\link{mrkfcv}}

Methods of classifier: \code{\link{lda}},
\code{\link{plsda}}, \code{\link{tree}}
}

