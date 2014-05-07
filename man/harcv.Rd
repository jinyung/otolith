\name{harcv}
\alias{harcv}
\title{Optimize harmonics number}
\usage{
harcv(X, Y, har, saveplot = FALSE, plotsize = 1000, method = "lda",
  run = 30, k = 5)
}
\arguments{
  \item{X}{matrix of EFA coefficients from
  \code{\link{rNEF}} object (the \code{coeff} value)}

  \item{Y}{a factor giving the grouping, e.g. the \code{sp}
  value from \code{\link{routine1}} object}

  \item{har}{numeric. number of harmonics to be tested.}

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
misclassification rate of the range of harmonics tested.
}
\description{
modified version of \code{\link{pccv}} function to optimize
dimension reduction.Different number of harmonics are used
in training the classification model and the result
(misclassification rate) is plotted out.
}
\seealso{
Similar: \code{\link{pccv}}

Which this function wraps: \code{\link{mrkfcv}}
}

