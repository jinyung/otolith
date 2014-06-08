\name{loadimg}
\alias{loadimg}
\title{Load image}
\usage{
loadimg(path, resize = TRUE, width)
}
\arguments{
  \item{path}{path of a single image}

  \item{resize}{logical. whether to resize the picture}

  \item{width}{width to resize to, in pixel. used only if
  \code{resize=TRUE}.  if \code{resize=TRUE} and width is
  not given, image is resized to 1280px in width by
  default}
}
\value{
an image matrix
}
\description{
To read image into R
}
\details{
All images loaded are converted into gray scale. Image
format is automatically detected from the extension.
Supported formats include \code{.jpg}, \code{.tif}, and
\code{.png}. This function is used in place of
\code{\link[EBImage]{readImage}} (from \code{EBImage}) for
this package to speed up the reading of image.
}
\seealso{
Which this function wraps: \code{\link[jpeg]{readJPEG}},
\code{\link[tiff]{readTIFF}}, \code{\link[png]{readPNG}}
}

