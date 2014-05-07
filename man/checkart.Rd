\name{checkart}
\alias{checkart}
\title{Check outline artifact}
\usage{
checkart(outdata)
}
\arguments{
  \item{outdata}{outline data to be checked. A matrix
  consists of xy coordinates}
}
\value{
\item{artifact}{logical. Whether the outline contains pixel
artifact} \item{index}{vector containg the index of pixel
that is problematic}
}
\description{
Checking whether the outline containing any pixel artifact
that will cause error in subsequent sampling of outline
}
\seealso{
\code{\link{orderOutline}}
}

