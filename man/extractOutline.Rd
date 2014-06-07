\name{extractOutline}
\alias{extractOutline}
\title{Extract outline}
\usage{
extractOutline(img, threshold = 0.3, plot = c("no", "overlay", "plain"))
}
\arguments{
  \item{img}{an image array read by
  \code{\link{loadImage}}}

  \item{threshold}{numeric. value < 1.}

  \item{plot}{\code{no} = no plot; \code{overlay} = plot of
  outline overlaid on image; \code{plain} = plot of plain
  outline}
}
\value{
a matrix of outline coords. note the outline is in
8-connected format.
}
\description{
Extract outline from a fish otolith image
}
\seealso{
Which this function wraps: \code{\link[EBImage]{ocontour}}
}

