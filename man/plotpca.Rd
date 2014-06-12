\name{plotpca}
\alias{plotpca}
\title{Plot PCA bubble plot}
\usage{
plotpca(pca, size = NULL, sizeamp = 3, class, saveplot = FALSE,
  plotsize = 1000)
}
\arguments{
  \item{pca}{\code{\link{prcomp}} object}

  \item{size}{optional. numeric. geometric size, e.g.
  \code{size} value from \code{\link{rGPA}} or
  \code{\link{rNEF}}. size is overlaid as bubble size in
  the plot.}

  \item{class}{factor. giving the species/grouping. NOTE:
  currently only support up to 12 levels.}

  \item{saveplot}{logical. whether to save the plot}

  \item{plotsize}{numeric. plot size for the plot saved
  into file (used only when \code{saveplot=TRUE}), unit in
  pixel.}

  \item{sizeamp}{numeric. to adjust the size of the
  bubble.}
}
\value{
Just the plot, consist of the first 3 PCs.
}
\description{
For preliminary visual assessment of data groups in PCA.
}
\seealso{
Function that wraps this function: \code{\link{rNEF}},
\code{\link{rGPA}}, \code{\link{routine1}}
}
\keyword{internal}

