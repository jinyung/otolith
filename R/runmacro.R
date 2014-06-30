#' Launch ImageJ shape indices macro
#' 
#' @description This function launches ImageJ macro from shell, the macro 
#'  will calculate dimensionless shape indices from otolith images.
#' @param imagejdir path to folder that contain the ImageJ's "ij.jar" file, 
#'  see details
#' @param opendir path to folder that contain the otolith images to be processed, 
#'  see details
#' @param savedir path to folder that the output files will be saved, see details.
#' @param threshold numeric. threshold value to be passed to ImageJ. range within 
#'  0-255, which is different from the \code{\link{img2landmark}} that use a 0-1 range.
#' @return 
#'  For each image, a file named \code{<image name>_shapedescriptors.txt} will be 
#'  generated, containing calculated shape descriptor indices
#' @details if \code{imagejdir} is not given, this function will check whether
#'  ImageJ can be launched from the ImageJ bundled within the package (only works 
#'  with the installation of package-zip that bundled with ImageJ), if not, 
#'  a GUI folder selector will pop out to prompt user to choose the folder to 
#'  locate ImageJ. Note: when using package-zip without ImageJ bundled with, user 
#'  need to put the  
#'  \code{batch_shape_indices_macro.ijm} file into the "macro" folder of ImageJ
#'  or Fiji, and also put the Particle 8 plugin into the plugin folder, see ref.
#'  When \code{opendir} or \code{savedir} are not given, a
#'  GUI selector will pop out, and user can choose the folder.
#' @export
#' @note Because \code{choose.dir()} is used in the function, 
#'  which is Windows specific, this function works in Windows only. 
#'  User is advised to give the path of directories manually if it doesn't work.
#'  The shell command was tested under Windows only, may not work in other OS.
#' @references The Particle 8 ImageJ's plugin at G Landhini's website 
#'  (\url{http://www.mecourse.com/landinig/software/software.html})

runmacro <- function (imagejdir, opendir, savedir, threshold=30) {
  if (missing(imagejdir)) {
    imagej <- system.file("ImageJ", package="otolith")
    # check if the ImageJ folder exist
    if (imagej == "") {
      # if not, let the user to choose the directory
      if (Sys.info()['sysname'] == "Windows") 
        imagej <- choose.dir(caption="Select folder containing ImageJ's \"ij.jar\" file")
      else 
        stop("you need to provide imagejdir (path to folder containing ImageJ's \"ij.jar\")")
    }
  } else {
    if (file.exists(imagejdir))
      imagej <- imagejdir
    else 
      stop("imagejdir path (path to folder to find ImageJ's \"ij.jar\") not exist")
  }
  if (missing(opendir)) { 
    if (Sys.info()['sysname'] == "Windows") 
      opendir <- paste0(choose.dir(caption=
               "Select input folder containing the images"),"\\")
    else
      stop("you need to provide opendir (path to folder containing the images)")
  }
  if (missing(savedir)) {
    if (Sys.info()['sysname'] == "Windows")
      savedir <- paste0(choose.dir(default=opendir, 
               caption="Select output folder to save the shape indices result"),"\\")
    else
      stop("you need to provide savedir (path to folder to save the shape indices result")
  }
  # normalize the path, avoid problem from shell in reading the path
  imagej <- normalizePath(imagej, mustWork = TRUE)
  opendir <- normalizePath(opendir, mustWork = TRUE)
  savedir <- normalizePath(savedir, mustWork = TRUE)
  # the command to be passed to shell
  cmd <- paste("cd", imagej, "&& ij.jar", 
               "-macro batch_shape_indices_macro", 
               paste(opendir, savedir, threshold, sep="*")) 
  # "cd" to any directory followed (in cmd.exe)
  # "&&" to pass second command
  # coz in shell everything need to be passed in one line 
  # some messages
  cat("ImageJ is found at:", imagej, 
      "\nThe input folder is:", opendir, 
      "\nThe output folder is:", savedir, 
      "\n\nImageJ will be opened now, please wait till the macro finish running\n\n")
  shell(cmd)
}