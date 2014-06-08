\name{comb / equaldist}
\alias{comb}
\alias{equaldist}
\title{Sample 'comb'/'equally spaced' -style curvilinear semilandmarks}
\usage{
comb(outdata, nd = 50, plot = FALSE, label = FALSE)
equaldist(outdata, nd = 50, plot = FALSE, label = FALSE)
}
\arguments{
  \item{outdata}{outline data. a matrix consists of xy
  coordinates. works on both 4-connected and 8-connected
  outline.}

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
Sample curvillinear semi-landmarks on a closed outline, using either the 'comb' or 
the 'equally-spaced' methods.
}
\details{
  \code{comb} divides the outline into two parts by larget PC axis, a 'comb' is
  then used to guide the sampling of the semi-landmarks. See reference for 'comb' 
  method explanations.
  
  \code{comb} produces errorneous results for complex outline with high
  outline curvature, e.g. at the 'excisura' site of some fish
  otoliths, thus, the \code{equaldist} method is preferable for fish otoliths and 
  is used in the \code{\link{img2landmark}} wrapper.
  
  \code{equaldist} divides the outline in similar way as \code{comb}, then semi
  -landmarks are sampled in equal intervals in each part, along the outline.
}
\references{
  Ponton, D. (2006). Is geometric morphometrics efficient for
  comparing otolith shape of different fish species?.
  \emph{Journal of Morphology}, 267(6), 750-757.
  
  Sheets, H. D., Covino, K. M., Panasiewicz, J. M., & Morris,
  S. R. (2006). Comparison of geometric morphometric outline
  methods in the discrimination of age-related differences in
  feather shape. \emph{Frontiers in Zoology}, 3(1), 1-12.
}