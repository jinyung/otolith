\name{runmacro}
\alias{runmacro}
\title{Launch ImageJ outline extraction macro}
\usage{
runmacro(imagej = c("built-in", "select"), input = NULL, output = NULL)
}
\arguments{
  \item{imagej}{path to folder that contain the ImageJ's
  "ij.jar" file, see details.}

  \item{input}{path to folder that contain the otolith
  images to be processed, see details.}

  \item{output}{path to folder that the output files will
  be saved, see details.}
}
\value{
For each image, there will be three files generated:
\tabular{ll}{ 1. \code{<image name>_outline.tif}\tab :
image showing the outline for inspection\cr 2. \code{<image
name>_xycoords.txt}\tab : outline XY coordinates\cr 3.
\code{<image name>_shapedescriptors.txt}\tab : containing
calculated shape descriptor indices\cr }
}
\description{
This function launches ImageJ macro from shell, the macro
will extract the outline from otolith images.
}
\details{
When \code{imagej="built-in"} (the default), ImageJ is
launched from the ImageJ bundled within the package (only
works with the installation of package-zip that bundled
with ImageJ), otherwise, when \code{imagej="select"}, GUI
selector will pop out, and user can choose the folder to
locate ImageJ.  Note: when using package-zip without ImageJ
bundled with, user need to put the
\code{batch_outline_extraction_macro.ijm} file into the
"macro" folder of ImageJ or Fiji. When \code{input=NULL} or
\code{output=NULL} (the default), GUI selector will pop
out, and user can choose the folder, otherwise, character
input of path is also accepted.
}
\note{
Because \code{choose.dir()} is used in the function, which
is Windows specific, this function works in Windows only.
The shell command is also specific to Windows, not tested
under other OS.
}
\references{
The Particle 8 ImageJ's plugin at G Landhini's website
(\url{http://www.mecourse.com/landinig/software/software.html})
}

