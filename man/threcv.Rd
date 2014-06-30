\name{threcv}
\alias{threcv}
\title{Optimize threshold on posterior probability}
\usage{
threcv(X, Y, saveplot = FALSE, plotsize = 1000, run = 30, k = 5, level)
}
\arguments{
  \item{X}{matrix of PC scores from \code{\link{rGPA}}
  object (the \code{score} value)}

  \item{Y}{a factor giving the grouping, e.g. the \code{sp}
  value from \code{\link{routine1}} object or obtain from
  \code{\link{getclass}}}

  \item{saveplot}{logical. The plot will be saved if
  \code{TRUE}, and no plot is displayed in the window. Note
  that outliers of boxplot are not plotted.}

  \item{plotsize}{numeric. Plot size for the plot saved
  into file (used only when \code{saveplot=TRUE}), unit in
  pixel.}

  \item{run}{argument to be passed to \code{\link{mrkfcv}}}

  \item{k}{argument to be passed to \code{\link{mrkfcv}}}

  \item{level}{numeric vector giving the threshold range to
  be tested. see details.}
}
\value{
a matrix and a plot are returned, giving the overall
accuracy and the total reported prediction percentage over
a range of threshold values.
}
\description{
To estimate the effect of setting a threshold on posterior
probability to reject "unsure" predictions.
}
\details{
only \code{\link{lda}} is supported now. \code{lda} decide
the class (e.g. species) based on the calculated posterior
probability. New sample is assigned to the class that has
the higher posterior probability. We can reject the
classification result if the posterior probability is low
(the "unsure" classification). In this way we can increase
the overall accuracy (correct predictions out of all the
reported results), but the tradeoff is the decrease in
number of reported classification result. The function test
the threshold value of 0.50 0.60 0.65 0.70 0.75 0.80 0.85
0.90 0.95 0.99 if \code{level} is not provided. Since
thresholding is just an option, the choice of threshold
value is really arbitrary and depends on one's objectives.
}
\references{
Beleites, C., & Salzer, R. (2008). Assessing and improving
the stability of chemometric models in small sample size
situations.  \emph{Analytical and Bioanalytical Chemistry},
390(5), 1261-1271.
}
\seealso{
Similar: \code{\link{pccv}}, \code{\link{harcv}}

Which this function wraps: \code{\link{mrkfcv}}
}

