\name{sprdist}
\alias{sprdist}
\title{Calculate within species Procrustes distance}
\usage{
sprdist(A, group)
}
\arguments{
  \item{A}{A p x k x n array of landmark configuration}

  \item{group}{A factor giving the species or the grouping
  of \code{A}, must be of same length and order as
  \code{A}}
}
\value{
\item{stat}{matrix of result of calculated within species
Procrustes distance} \item{meanshape}{Array of meanshape
configurations of each species/ group}
}
\description{
A wrapper to calculate the within species Procrustes (or
Riemannian, as used by \code{\link{Morpho}}) distance
}
\seealso{
function that wraps this function: \code{\link{rGPA}}

which this function wraps: \code{\link{kendalldist}}
\code{\link{rotonto}} \code{\link{mshape}}
}

