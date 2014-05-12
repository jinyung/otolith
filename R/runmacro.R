#' Launch ImageJ outline extraction macro
#' 
#' @description This function launches ImageJ macro from shell, the macro 
#'  will extract the outline from otolith images.
#' @param imagej path to folder that contain the ImageJ's "ij.jar" file, see details.
#' @param input path to folder that contain the otolith images to be processed, see details.
#' @param output path to folder that the output files will be saved, see details.
#' @return 
#' For each image, there will be three files generated:
#' \tabular{ll}{
#' 1. \code{<image name>_outline.tif}\tab : image showing the outline for inspection\cr
#' 2. \code{<image name>_xycoords.txt}\tab : outline XY coordinates\cr
#' 3. \code{<image name>_shapedescriptors.txt}\tab : containing calculated shape descriptor indices\cr
#' }
#' @details When \code{imagej="built-in"} (the default), ImageJ is launched from 
#'  the ImageJ bundled within the package (only works with the installation of 
#'  package-zip that bundled with ImageJ), otherwise, when \code{imagej="select"},
#'  GUI selector will pop out, and user can choose the folder to locate ImageJ. 
#'  Note: when using package-zip without ImageJ bundled with, user need to put the  
#'  \code{batch_outline_extraction_macro.ijm} file into the "macro" folder of ImageJ
#'  or Fiji. When \code{input=NULL} or \code{output=NULL} (the default), 
#'  GUI selector will pop out, and user can choose the folder, otherwise, 
#'  character input of path is also accepted.
#' @export
#' @note Because \code{choose.dir()} is used in the function, 
#'  which is Windows specific, this function works in Windows only. 
#'  The shell command is also specific to Windows, not tested under
#'  other OS.
#' @references The Particle 8 ImageJ's plugin at G Landhini's website 
#'  (\url{http://www.mecourse.com/landinig/software/software.html})

runmacro <- function (imagej=c("built-in","select"), input=NULL, output=NULL) {
  if (imagej == "built-in") {
    # only works with the package-zip that bundled with ImageJ
    imagej <- system.file("ImageJ", package="otolith")
  } else {
    imagej <- choose.dir(caption="Select folder containing ImageJ's \"ij.jar\" file")
    imagej <- gsub("\\", "/", imagej, fix=TRUE)
  }
  if (is.null(input)) { 
    input <- paste0(choose.dir(caption="Select input folder containing the images"),"\\")
    input <- gsub("\\", "/", input, fix=TRUE)
  }
  if (is.null(output)) {
    output <- paste0(choose.dir(caption="Select output folder to save the outlines"),"\\")
    output <- gsub("\\", "/", input, fix=TRUE)
  }
  # the command to be passed to shell
  cmd <- paste("cd", imagej, "&& ij.jar", 
               "-macro batch_outline_extraction_macro", 
               paste(input, output, sep="*")) 
  # "cd" to any directory followed
  # "&&" to pass second command
  # coz in shell everything need to be passed in one line 
  # some messages
  cat("ImageJ is found at:", imagej, 
      "\nThe input folder is:", input, 
      "\nThe output folder is:", output, 
      "\n\nImageJ will be opened now, please wait till the macro finish running\n\n")
  shell(cmd)
}