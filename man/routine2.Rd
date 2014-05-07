\name{routine2}
\alias{routine2}
\title{Process outline and sample semi-landmarks}
\usage{
routine2(imagejfolder, rfolder, write = TRUE, plot = TRUE, nd = 50,
  scale = NULL)
}
\arguments{
  \item{imagejfolder}{Path of the folder containg the
  outline xy coordinate files}

  \item{rfolder}{Path of the output folder (must exist
  already)}

  \item{write}{Whether to save new files of record of
  outline with defects, re-ordered outline coordinates,
  curvilinear semi-landmarks and a combined file
  curvilinear semilandmarks of all configuration in TPS
  format.  TPS format is compatible with many other
  morphometrics software.}

  \item{plot}{Logical. Whether to save the plots of outline
  with semi-landmarks.}

  \item{nd}{Numeric. Number of sampling interval. See
  \code{\link{equaldist}} function.}

  \item{scale}{scale in the unit of um/px, i.e. 1 pixel of
  the image equal to how many um?}
}
\value{
\item{rmindex}{index of the outline with pixel artifacts}
\item{rmfile}{file label of the outline with pixel
artifacts} \item{nlist}{list of file labels after checking}
\item{A}{array of semi-landmarks configurations}
}
\description{
A wrapper function to do the routine of reading the files
from imagej output, then read the outline files, checking
for pixel artifacts, re-order the outline coordinates and
laying curvilinear semi-landmarks
}
\seealso{
similar: \code{\link{routine1}}

which this function wraps: \code{\link{orderOutline}
\link{equaldist}}
}

