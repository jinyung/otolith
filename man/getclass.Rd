\name{getclass}
\alias{getclass}
\title{Extract class name from landmark configurations}
\usage{
getclass(landmark, extract = c(1, 6))
}
\arguments{
  \item{landmark}{p x k x n array of landmark. n dimension
  has to be named.}

  \item{extract}{numeric vector of 2 values, with the first
  giving the starting character for extraction, and second
  for ending character for extraction}
}
\value{
class name
}
\description{
A simple wrapper to extract class name, essentially a
wrapped \code{\link{substr}} function
}

