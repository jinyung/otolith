\name{routine1}
\alias{routine1}
\title{Process shape descriptor files}
\usage{
routine1(imagejfolder, write = TRUE, plot = FALSE, label = TRUE,
  extract = c(1, 6))
}
\arguments{
  \item{imagejfolder}{path of the folder containg the shape
  descriptor files}

  \item{write}{logical. whether to save a new file of
  combined descriptor.}

  \item{plot}{logical. whether to plot PCA for preliminary
  assessment of descriptor data. Only graphical output of
  first 3 PCs.}

  \item{label}{logical. Whether to extract label from the
  image file name.}

  \item{extract}{a numeric vector of two number telling
  which characters to extract from the file name to become
  the species name, if \code{label=TRUE}. The first is the
  starting position, second is the last position of
  character to extract.}
}
\value{
A dataframe of combined shape descriptor data
}
\description{
A wrapper function to do the routine of reading the files
from imagej output and combine the descriptor files into a
single file.
}
\seealso{
\code{\link{routine2}}
}

