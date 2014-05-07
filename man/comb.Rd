\name{comb}
\alias{comb}
\title{Sample 'comb'-style curvilinear semilandmarks}
\usage{
comb(outdata, nd = 50, plot = FALSE, label = FALSE)
}
\arguments{
  \item{outdata}{outline data. A matrix consists of xy
  coordinates}

  \item{nd}{numeric. Number of equally spaced divides along
  the outline for each part of the outline.}

  \item{plot}{logical. Whether to plot}

  \item{label}{logical. Whether to label the plot}
}
\value{
semi-landmarks sampled and the plot of the sampled
semi-landmarks on the outline if \code{plot=TRUE}
}
\description{
Sample semi-landmarks on outline using the 'comb' method.
}
\note{
produce errorneous results for complex outline with high
outline curvature especially at the 'excisura' site.
}
\references{
Ponton, D. (2006). Is geometric morphometrics efficient for
comparing otolith shape of different fish species?. Journal
of Morphology, 267(6), 750-757.

Sheets, H. D., Covino, K. M., Panasiewicz, J. M., & Morris,
S. R. (2006). Comparison of geometric morphometric outline
methods in the discrimination of age-related differences in
feather shape. Frontiers in Zoology, 3(1), 1-12.
}
\seealso{
\code{\link{equaldist}}
}

