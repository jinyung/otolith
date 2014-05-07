\name{orderOutline}
\alias{orderOutline}
\title{Order outline pixels}
\usage{
orderOutline(outdata)
}
\arguments{
  \item{outdata}{outline data that need to be re-ordered.
  A matrix consists of xy coordinates}
}
\value{
re-ordered outline coordinates
}
\description{
Order an outline so that the successive points' indices are
in successive order. This is important because
\code{\link{equaldist}} function assumes the coordinates
are ordered.
}
\seealso{
\code{\link{checkart}}
}

