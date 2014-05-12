\name{reland}
\alias{reland}
\title{Re-arrage semi-landmark configuration}
\usage{
reland(landdata, type = c("ori", "rev", "flip", "fliprev"))
}
\arguments{
  \item{landdata}{semi-landmarks data generated from
  \code{\link{equaldist}} functions. A matrix/dataframe
  consists of two columns of x and y coordinates}

  \item{type}{character of \code{"ori"}, \code{"rev"},
  \code{"flip"}, or \code{"fliprev"} to set the
  configuration of the semi-landmarks. See details.}
}
\value{
matrix of xy coordinates of re-arranged semi-landmarks
data.
}
\description{
To flip or change the order of semi-landmark configuration
}
\details{
There are four possible configurations of semi-landmark
sampling for any set of otolith, depending on the side of
the otolith and the initial orietation of the image. This
function can be used to standardize the sampling of
semi-landmarks.

\code{type} argument change the semi-landmarks as if the
otolith image was taken from different end on the left
(\code{"rev"}), or as if the otolith is from the other side
(\code{"flip"}), or as if the otolith is from the other
side and image taken from different end on the left
(\code{"fliprev"}). \code{"ori"} do nothing on the
semi-landmarks.
}

