\name{img2landmark}
\alias{img2landmark}
\title{Convert image to semi-landmarks}
\usage{
img2landmark(resize = TRUE, opendir, saveoutline = TRUE,
  savelandmark = TRUE, plot = c("no", "overlay", "plain"), savedir,
  threshold = 0.3, nd = 50, suppress = FALSE)
}
\arguments{
  \item{opendir}{directory path. folder containing the
  images to be opened.}

  \item{saveoutline}{logical. whether to save the outline.
  If \code{TRUE}, .csv file(s) of the coordinates of each
  image will be saved}

  \item{savelandmark}{logical. whether to save the sampled
  semi-landmarks. If \code{TRUE}, .csv file(s) of the
  coordinates of each image will be saved and a .tps file
  of all images will be saved.}

  \item{savedir}{directory path. folder for files to be
  saved in if \code{saveoutline=TRUE} or
  \code{savelandmark=TRUE} (default).}

  \item{threshold}{numeric. argument passed to
  \code{\link{extractOutline}}}

  \item{nd}{numeric. argument passed to
  \code{\link{equaldist}}}

  \item{suppress}{logical. whether to supress the system
  messages on running status.}

  \item{plot}{\code{no} = no plot; \code{overlay} = plot of
  outline overlaid on image; \code{plain} = plot of plain
  outline}
}
\value{
if \code{saveoutline=TRUE} or \code{savelandmark=TRUE} or
\code{plot != "no"}, files will be saved into specified
folder. Other values: \item{outline}{list of outline xy
coordinates} \item{landmark}{array of sampled
semi-landmarks xy coordinates}
}
\description{
A wrapper function to convert image files to semi-landmarks
}
\seealso{
Which this function wraps: \code{\link{extractOutline}},
\code{\link{equaldist}}
}

