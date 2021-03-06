\name{rNEF}
\alias{rNEF}
\title{Normalized elliptical Fourier analysis routine}
\usage{
rNEF(A, plotpca = FALSE, class, ...)
}
\arguments{
  \item{A}{p x k x n array of raw semilandmarks
  configurations}

  \item{plotpca}{logical. Whether to plot the PCA for quick
  preliminary view.  Plot first 3 PCs}

  \item{class}{factor giving the species/grouping, used for
  plotting purpose, e.g. using \code{sp} value from
  \code{\link{routine1}} or using the
  \code{\link{getclass}} function. Used only when
  \code{plotpca=TRUE}}

  \item{...}{other arguments passed to
  \code{\link{plotpca}}}
}
\value{
\item{coeff}{matrix of all harmonic coefficients}
\item{expvar}{explained variation of the NEF harmonics
(first harmonic is removed from the calculation, see
reference)} \item{size}{size of the semi-major axis of the
first fitting ellipse. See reference}
}
\description{
A wrapper functions using the codes from J Claude's book,
for routine use to run normalized elliptical Fourier
analysis.
}
\references{
Claude J. (2008). \emph{Morphometrics with R}. Springer.
}
\seealso{
\code{\link{rGPA}}
}

