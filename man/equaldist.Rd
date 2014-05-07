\name{equaldist}
\alias{equaldist}
\title{Sample 'equally spaced'-style curvilinear semi-landmarks}
\usage{
equaldist(outdata, nd = 50, plot = FALSE, label = FALSE)
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
Sample points along outline. Outline is divided into two
parts by larget PC axis, semi-landmarks are then sampled in
equal intervals in each part.
}
\seealso{
\code{\link{routine2}}
}

