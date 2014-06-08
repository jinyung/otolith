#' Convert image to semi-landmarks
#' 
#' @description A wrapper function to convert image files to semi-landmarks
#' @param opendir directory path. folder containing the images to be opened.
#' @param threshold numeric. argument passed to \code{\link{extractout}}
#' @param resize logical. whether to resize the image when reading the image. 
#'  argument passed to \code{\link{loadimg}}
#' @param saveoutline logical. whether to save the outline. If  
#'  \code{TRUE}, .csv file(s) of the coordinates of each image will be saved 
#' @param savelandmark logical. whether to save the sampled semi-landmarks. If  
#'  \code{TRUE}, .csv file(s) of the coordinates of each image will be saved and 
#'  a .tps file of all images will be saved.
#' @param savedir directory path. folder for files to be saved in if 
#'  \code{saveoutline=TRUE} or \code{savelandmark=TRUE} or \code{plot != "no"}
#'  (default).
#' @param nd numeric. argument passed to \code{\link{equaldist}}
#' @param suppress logical. whether to supress the system messages on running status.
#' @param plot \code{"no"} = no plot; \code{"overlay"} = plot of outline 
#'  and semi-landmarks overlaid on image; \code{"plain"} = plot of plain outline
#' @param start numeric. number of which the 
#' @param extract argument passed to \code{\link{getclass}}
#' @return 
#'  if \code{saveoutline=TRUE} or \code{savelandmark=TRUE} or \code{plot != "no"}, 
#'  files will be saved into specified folder. Other values:
#'  \item{outline}{list of outline xy coordinates}
#'  \item{landmark}{array of sampled semi-landmarks xy coordinates}
#'  \item{class}{factor of class name}
#' @importFrom geomorph writeland.tps
#' @seealso
#' Which this function wraps: \code{\link{extractout}}, \code{\link{equaldist}}
#' @export

img2landmark <- function (opendir, threshold = 0.3, resize = TRUE, nd = 50, 
                          saveoutline = TRUE, savelandmark = TRUE,  
                          plot = c("no", "overlay", "plain"), savedir,                           
                          suppress=FALSE, start, extract = c(1, 6)) {
  # package and arguments check
  require(EBImage)
  if (missing(opendir)) {
    if (Sys.info()['sysname'] == "Windows")
      opendir <- choose.dir(default = getwd(), 
                 caption = "Select folder containing the images")
    else
      stop("you need to provide path to folder containing the images")
  }
  plot <- match.arg(plot)
  if (saveoutline | savelandmark | plot != "no") {
    if (missing(savedir)) {
      if (Sys.info()['sysname'] == "Windows")
        savedir <- choose.dir(default = opendir, 
                   caption = "Select folder to save the outlines")
      else
        stop("you need to provide path to folder to save the file(s)")
    }
  }
  # get the file list and check
  imglist <- list.files(opendir)
  ext <- sapply(imglist, function(x) rev(unlist(strsplit(x, "[.]")))[1])
    # split character at ".", take the last splitted part
    # fixed bug for file name with "." in it and fixed bug for Thumb.db
    # using sapply to avoid loop
  rmindex <- which(ext != "jpg" & ext != "tif" & ext!= "png")
    # remove the files belongs to format not supported
  imglist <- imglist[-rmindex]
  imgname <- sapply(imglist, function(x) substr(x, 1, nchar(x) - 4))
  # initialize 
  outline <- vector("list", length(imglist))
  names(outline) <- imglist
  landmark <- array(data = NA, dim = c(nd * 2, 2, length(imgname)), 
              dimnames = list(NULL, NULL, imgname))
  # running thru the images
  if (missing(start))
    start <- 1
  each <- NULL
  begin <- NULL
  for (i in start:length(imglist)) { 
    if (i == 1)
      begin <- proc.time()
    img <- loadimg(paste(opendir, imglist[i], sep = "/"), resize = resize)
    outline[[i]] <- extractout(img, threshold = threshold, plot = "no")
    landmark[, , i] <- equaldist(outline[[i]], nd = nd)
    # save files
    if (saveoutline) {
      write.csv(outline[[i]], file = paste0(savedir, "/", imgname[[i]], 
                "-outline.csv"), row.names = FALSE)    
    }
    if (savelandmark) {
      write.csv(landmark[, , i], file = paste0(savedir, "/curvillinearLand_", 
                imgname[i], ".csv"), row.names = FALSE)
    }
    # plot
    if (plot != "no") {
      x <- dim(img)[1]
      y <- dim(img)[2]
      tiff(width = x, height = y,   
           filename = paste0(savedir, "/curvillinearLand_", imgname[i], ".tif"),
           compression = "lzw", res = 128)
      par(mar = rep(0, 4))
      if (plot == "overlay") {
        plot(x, y, xlim = c(0, x), ylim = c(0, y), asp = 1, type = "n")
        rasterImage(img, 0, y, x, 0)
        polygon(outline[[i]], border = "green")
      } else if (plot == "plain") {
        plot(outline[[i]], asp = 1, type = 'n')
        polygon(outline[[i]])
      }
      # label landmarks
      points(landmark[- c(1, nd + 1), , i], pch = 21, col = "green", bg = "blue")
      points(landmark[c(1, nd + 1), , i], pch = 21, col = "green", bg = "red")
      dev.off()
      }
      # estimate remaining time
      if (i < 3) {
        remaining <- "estimating..." 
        if (i == 1 & !suppress) {
          cat("The input folder is selected as:", opendir, 
              "\nThe saving options are:\n\t\tsaveoutline\t\t:", saveoutline, 
              "\n\t\tsavelandmark\t:", savelandmark, 
              "\n\t\tplot\t\t\t\t\t\t\t\t\t:", plot, "\n\n")      
        }
      } else {
        each <- (proc.time() - begin)[3] / i
        # using average of all previous runs to estimate remaining time
        remaining <- format(.POSIXct(((length(imglist) - i) * each), tz = "GMT"), 
                     "%H:%M:%S") 
                   # change format to hour:minute:second
                   # ref: http://stackoverflow.com/questions/11017933/format-time-span-to-show-hours-minutes-seconds
      }
      if (!suppress) {
        cat ("\r**Processing: ", imglist[i], " (", i, "/", 
               length(imglist), ") ",  "; Estimated remaining time = ", 
               remaining, "                                 ")
            # \r rewrite current line, blank ensure clean overwrite 
      }
    }
  if (savelandmark) {
    require(geomorph)
    writeland.tps(landmark, paste(savedir, "semi-landmarks.tps", sep = "/"))
  }
  if (!suppress & (saveoutline | savelandmark | plot != "no"))
    cat("\nProcessing completed, the files are saved at: ", savedir)
  invisible(list(outline = outline, landmark = landmark, 
            class = getclass(landmark, extract)))
}