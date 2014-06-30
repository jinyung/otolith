\name{pccv / harcv}
\alias{pccv}
\alias{harcv}
\title{Optimizing PC/ harmonics number}
\usage{
pccv(X, Y, pc, saveplot = FALSE, plotsize = 1000, method = c("lda", "tree", 
  "plsda"), run = 30, k = 5)
  
harcv(X, Y, har, saveplot = FALSE, plotsize = 1000, method = c("lda", "tree", 
  "plsda"), run = 30, k = 5)
}

\arguments{
  \item{X}{[for \code{pccv}] matrix of PC scores from \code{\link{rGPA}}
  object (the \code{score} value) [for \code{harcv}] matrix of EFA 
  coefficients from \code{\link{rNEF}} object (the \code{coeff} value)}

  \item{Y}{a factor giving the groupings, e.g. value from \code{\link{getclass}} 
  or the \code{sp} value from \code{\link{routine1}} object}

  \item{pc}{[for \code{pccv}] numeric. number of PCs to be tested.}
  \item{har}{[for \code{harcv}] numeric. number of harmonics to be tested.}

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
misclassification rate of the range of PCs/ harmonics tested.
}
\description{
A wrapper function to optimize dimension reduction using
multiple-run k-fold cross validation. Different number of
PCs/ harmonics are used in training the classification model and the
result (misclassification rate) is plotted out.
}
\seealso{
Similar: \code{\link{threcv}}

Which this function wraps: \code{\link{mrkfcv}}
}
