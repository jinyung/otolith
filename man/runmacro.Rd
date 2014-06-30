\name{runmacro}
\alias{runmacro}
\title{Launch ImageJ shape indices macro}
\usage{
runmacro(imagejdir, opendir, savedir, threshold = 30)
}
\arguments{
  \item{imagejdir}{path to folder that contain the ImageJ's
  "ij.jar" file, see details}

  \item{opendir}{path to folder that contain the otolith
  images to be processed, see details}

  \item{savedir}{path to folder that the output files will
  be saved, see details.}

  \item{threshold}{numeric. threshold value to be passed to
  ImageJ. range within 0-255, which is different from the
  \code{\link{img2landmark}} that use a 0-1 range.}
}
\value{
For each image, a file named \code{<image
name>_shapedescriptors.txt} will be generated, containing
calculated shape descriptor indices
}
\description{
This function launches ImageJ macro from shell, the macro
will calculate dimensionless shape indices from otolith
images.
}
\details{
if \code{imagejdir} is not given, this function will check
whether ImageJ can be launched from the ImageJ bundled
within the package (only works with the installation of
package-zip that bundled with ImageJ), if not, a GUI folder
selector will pop out to prompt user to choose the folder
to locate ImageJ. Note: when using package-zip without
ImageJ bundled with, user need to put the
\code{batch_shape_indices_macro.ijm} file into the "macro"
folder of ImageJ or Fiji, and also put the Particle 8
plugin into the plugin folder, see ref.  When
\code{opendir} or \code{savedir} are not given, a GUI
selector will pop out, and user can choose the folder.
}
\note{
Because \code{choose.dir()} is used in the function, which
is Windows specific, this function works in Windows only.
User is advised to give the path of directories manually if
it doesn't work.  The shell command was tested under
Windows only, may not work in other OS.
}
\references{
The Particle 8 ImageJ's plugin at G Landhini's website
(\url{http://www.mecourse.com/landinig/software/software.html})
}

