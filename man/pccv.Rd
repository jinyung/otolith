\name{pccv}
\alias{pccv}
\title{Optimize PC number}
\usage{
pccv(X, Y, pc, saveplot = FALSE, plotsize = 1000, method = "lda",
  run = 30, k = 5)
}
\arguments{
  \item{X}{matrix of PC scores from \code{\link{rGPA}}
  object (the \code{score} value)}

  \item{Y}{a factor giving the grouping, e.g. the \code{sp}
  value from \code{\link{routine1}} object}

  \item{pc}{numeric. number of PCs to be tested.}

  \item{saveplot}{logical. The plot will be saved if
  \code{TRUE}, and no plot is displayed in the window. Note
  that outliers of boxplot are not plotted.}

  \item{plotsize}{numeric. Plot size for the plot saved
  into file (used only when \code{saveplot=TRUE}), unit in
  pixel.}

  \item{method}{argument to be passed to
  \code{\link{mrkfcv}}}

  \item{run}{argument to be passed to \code{\link{mrkfcv}}}

  \item{k}{argument to be passed to \code{\link{mrkfcv}}}
}
\value{
a matrix and a plot are returned, giving the
misclassification rate of the range of PCs tested.
}
\description{
A wrapper function to optimize dimension reduction using
multiple-run k-fold cross validation. Different number of
PCs are used in training the classification model and the
result (misclassification rate) is plotted out.
}
\seealso{
Similar: \code{\link{harcv}}, \code{\link{threcv}}

Which this function wraps: \code{\link{mrkfcv}}
}

