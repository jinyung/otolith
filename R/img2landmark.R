#' Convert image to semi-landmarks
#' 
#' @description A wrapper function to convert image files to semi-landmarks
#' @param opendir directory path. folder containing the images to be opened.
#' @param saveoutline logical. whether to save the outline. If  
#'  \code{TRUE}, .csv file(s) of the coordinates of each image will be saved 
#' @param savelandmark logical. whether to save the sampled semi-landmarks. If  
#'  \code{TRUE}, .csv file(s) of the coordinates of each image will be saved and 
#'  a .tps file of all images will be saved.
#' @param savedir directory path. folder for files to be saved in if 
#'  \code{saveoutline=TRUE} or \code{savelandmark=TRUE} (default).
#' @param threshold numeric. argument passed to \code{\link{extractOutline}}
#' @param nd numeric. argument passed to \code{\link{equaldist}}
#' @param suppress logical. whether to supress the system messages on running status.
#' @param plot \code{no} = no plot; \code{overlay} = plot of outline overlaid 
#'  on image; \code{plain} = plot of plain outline
#' @return 
#'  if \code{saveoutline=TRUE} or \code{savelandmark=TRUE} or \code{plot != "no"}, 
#'  files will be saved into specified folder. Other values:
#'  \item{outline}{list of outline xy coordinates}
#'  \item{landmark}{array of sampled semi-landmarks xy coordinates}
#' @importFrom geomorph writeland.tps
#' @seealso
#' Which this function wraps: \code{\link{extractOutline}}, \code{\link{equaldist}}
#' @export

img2landmark <- function (resize = TRUE, 
                          opendir, 
                          saveoutline = TRUE, savelandmark = TRUE,  
                          plot = c("no", "overlay", "plain"), savedir,                           
                          threshold = 0.3, nd = 50, suppress=FALSE) {
  # package and arguments check
  require(EBImage)
  if (missing(opendir)) {
    if (Sys.info()['sysname'] == "Windows")
      opendir <- choose.dir(default=getwd(), caption="Select folder containing the images")
    else
      stop("you need to provide path to folder containing the images")
  }
  if (saveoutline | savelandmark | plot != "no") {
    if (missing(savedir)) {
      if (Sys.info()['sysname'] == "Windows")
        savedir <- choose.dir(default=opendir, caption="Select folder to save the outlines")
      else
        stop("you need to provide path to folder to save the file(s)")
    }
  }
  #switch(match.arg(type), 
  #       jpg = require(jpeg),
  #       tif = require(tiff),
  #       png = require(png))
  # file names contain in opendir
  imglist <- list.files(opendir)
  imgname <- unlist(strsplit(imglist, "[.]")) # split character at "."
  imgname <- imgname[c(TRUE, FALSE)] # get rid of the even characters
  # initialize 
  outline <- vector("list", length(imglist))
  names(outline) <- imglist
  landmark <- array(data=NA, dim=c(nd * 2, 2, length(imglist)), dimnames=imglist)
  # running thru the images
  for (i in 1:length(imglist)) {
    if (!suppress) 
      cat("\nprocessing...", imglist[i], " (", i, "/", length(imglist), ") ")
    flush.console()
    img <- loadImage(paste(opendir, imglist[i], sep="/"), resize=resize)
    #switch(match.arg(type), 
    #       jpg = assign("img", readJPEG(paste(opendir, imglist[i], sep="/"))), 
    #       tif = assign("img", readTIFF(paste(opendir, imglist[i], sep="/"))),
    #       png = assign("img", readPNG(paste(opendir, imglist[i], sep="/"))))
    # resize
    #if (dim(img)[2] > 1280)
    #  img <- resize(img, h = 1280)
    # greyscale
    #if (getNumberOfFrames(img) == 3) 
    # img <- img[, , 1] * 0.3 + img[, , 2] * 0.59 + img[, , 3] * 0.11
    # correct read format
    #img <- t(img)
    #img <- EBImage::flip(img)
    # extract outline and sample semi-landmark
    outline[[i]] <- extractOutline(img, threshold=threshold, plot="no")
    landmark[, , i] <- equaldist(outline[[i]], nd=nd)
    # save files
    if (saveoutline) {
      write.csv(outline[[i]], file=paste0(savedir, "/", imgname[[i]], 
                "-outline.csv"), row.names=FALSE)    
    }
    if (savelandmark) {
      write.csv(landmark[, , i], file=paste0(savedir, "/curvillinearLand_", 
                imgname[i], ".csv"), row.names=FALSE)
    }
    # plot
    if (plot != "no") {
      x <- dim(img)[1]
      y <- dim(img)[2]
      tiff(width=x, height=y,   
           filename=paste0(savedir, "/curvillinearLand_", imgname[i], ".tif"),
           compression="lzw", res=128)
      par(mar=rep(0, 4))
      if (plot == "overlay") {
        plot(x, y, xlim=c(0, x), ylim=c(0, y), asp=1, type="n")
        rasterImage(img, 0, y, x, 0)
        polygon(outline[[i]], border="green")           
      } else if (plot == "plain") {
        plot(outline[[i]], asp=1, type='n')
        polygon(outline[[i]])
      }
      # label landmarks
      points(landmark[- c(1, nd + 1), , i], pch=21, col="green", bg="blue")
      points(landmark[c(1, nd + 1), , i], pch=21, col="green", bg="red")
      dev.off()
    }
  }
  if (savelandmark){
    require(geomorph)
    writeland.tps(landmark, paste(savedir, "semi-landmarks.tps", sep="/"))
  }
  if (!suppress)
    cat("\nConversion completed, the files are saved at: ", savedir)
  invisible(list(outline=outline, landmark=landmark))
}