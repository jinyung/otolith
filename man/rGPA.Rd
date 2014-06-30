\name{rGPA}
\alias{rGPA}
\title{Generalized Procrustes analysis routine}
\usage{
rGPA(A, fix = NULL, plotpca = FALSE, class)
}
\arguments{
  \item{A}{p x k x n array of raw semilandmarks
  configurations}

  \item{fix}{a numeric vector giving the semi-landmakrs
  that not to be slided.  Default is all semilandmarks will
  be slided (See Note for issue).}

  \item{plotpca}{logical. whether to plot PCA for
  preliminary assessment of the data}

  \item{class}{factor giving the species/grouping, used for
  plotting purpose and calculation of within class range of
  Procrustes distance, e.g. using sp value from
  \code{\link{routine1}} or using \code{\link{getclass}}}
}
\value{
\item{tanc}{Procrustes aligned configurations projected
onto tangent shape space} \item{size}{centroid size}
\item{expvar}{summary of expalined variations by each PC}
\item{score}{matrix contain the PC scores, to be used for
training} \item{rdist}{see \code{\link{sprdist}}}
\item{pca}{the PCA model, to be used by
\code{\link{otopred}}} \item{meanshape}{meanshape of each
species/ groups} \item{mshape}{meanshape of all species/
groups}
}
\description{
a wrapper function for routine use to run generalized
Procrustes analysis.
}
\note{
There is an issue with \code{\link{otopred}} function,
which cannot handle sliding semi-landmarks yet and hence it
is recommended that all semi-landmarks should be fixed if
the user intend to use the \code{\link{otopred}} function
(6 May 2014).
}
\seealso{
which this function wraps: \code{\link{procSym}},
\code{\link{sprdist}}
}

