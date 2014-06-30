\name{agglda}
\alias{agglda}
\title{LDA model aggregation}
\usage{
agglda(X, Y, newdata = NULL, type = c("vote", "post"), k = 5, run = 100,
  threshold, prior = c("equal", "proportion"), suppress = FALSE)
}
\arguments{
  \item{X}{matrix/ dataframe of predictors, e.g. EFA
  coefficients/ PC scores selected using
  \code{\link{selectdim}}}

  \item{Y}{factor/ character giving the class, e.g. value
  obtained from \code{\link{getclass}} or \code{sp} value
  from \code{\link{routine1}} object}

  \item{newdata}{matrix/ dataframe of newdata to be
  predicted. see details.}

  \item{type}{prediction based on the majority vote
  (\code{"vote"}) from the submodels or based on the mean
  posterior probability from the submodels (\code{"post"})}

  \item{k}{the fold number for k-fold resampling}

  \item{run}{the iteration number for iteration of k-fold
  resampling}

  \item{threshold}{single numeric value of < 1. threshold
  of the proportion of majority vote or the mean posterior
  probability. predictions with less than the value will be
  reported as \code{NA}}

  \item{prior}{the prior used in \code{\link[MASS]{lda}}
  models, \code{"equal"} means all classes have same prior,
  \code{"proportion"} means prior according to the classes
  weight in the training data.}

  \item{suppress}{logical. whether to suppress the progress
  monitoring output}
}
\value{
\item{accuracy}{[\emph{evaluation mode}] the overall
accuracy in percent} \item{conmat}{[\emph{evaluation mode}]
confusion matrix} \item{stat}{[\emph{evaluation mode}]
matrix containing the statistics of each class, see
details} \item{total}{[\emph{evaluation mode}] the total
percent of reported prediction after threshold. give
\code{NULL} if \code{threshold} is not given}
\item{ind.prediction}{matrix containing the prediction
result on each training/ new specimens}
}
\description{
aggregate LDA models based on iterated k-fold resampling
method
}
\details{
If \code{newdata} is provided, the function is in the
\emph{prediction mode}, the aggregated model will be built
from \code{X} and \code{Y} and predicition is performed on
\code{newdata}. Otherwise, if \code{newdata = NULL}
(default) the function is in \emph{evaluation mode}.

In \emph{evaluation mode}, overall accuracy of the model
and the by-class statistics are calculated, similar to that
of \code{\link{mrkfcv2}}.  However, the statistics are
calculated based on the aggregated prediction.  See
reference for explanation on model aggregation and
thresholding.
}
\references{
Beleites, C., & Salzer, R. (2008). Assessing and improving
the stability of chemometric models in small sample size
situations. \emph{Analytical and Bioanalytical Chemistry},
390(5), 1261-1271.
}
\seealso{
\code{\link{kfcv}}, \code{\link{mrkfcv}}
}

