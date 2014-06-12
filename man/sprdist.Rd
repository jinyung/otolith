\name{sprdist}
\alias{sprdist}
\title{Calculate within species Procrustes distance}
\usage{
sprdist(A, class)
}
\arguments{
  \item{A}{p x k x n array of landmark configuration}

  \item{class}{a factor giving the species or the grouping
  of \code{A}, must be of same length and order as
  \code{A}. can be obtained using \code{\link{getclass}}}
}
\value{
\item{stat}{matrix of result of calculated within species
Procrustes distance} \item{meanshape}{array of meanshape
configurations of each species/ group}
}
\description{
A wrapper to calculate the within species Procrustes (or
Riemannian, as used by \code{\link{Morpho}}) distance
}
\seealso{
function that wraps this function: \code{\link{rGPA}}

which this function wraps: \code{\link{kendalldist}},
\code{\link{rotonto}}, \code{\link{mshape}}
}

