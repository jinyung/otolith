#' Load image
#' 
#' @description To read image into R 
#' @param path path of a single image
#' @param resize logical. whether to resize the picture
#' @param width width to resize to, in pixel. used only if \code{resize=TRUE}.
#'  if \code{resize=TRUE} and width not given, image is resize to 1280px in width
#'  by default
#' @return an image matrix
#' @importFrom jpeg readJPEG
#' @importFrom tiff readTIFF
#' @importFrom png readPNG
#' @details All images loaded are converted into gray scale. Image format is 
#'  automatically detected from the extension. Supported formats include \code{.jpg}, 
#'  \code{.tif}, and \code{.png}. This function is used in place of 
#'  \code{\link[EBImage]{writeImage}} for this package to speed up the reading of
#'  image.
#' @seealso 
#'  Which this function wraps: \code{\link[jpeg]{readJPEG}}, 
#'  \code{\link[tiff]{readTIFF}}, \code{\link[png]{readPNG}} 
#' @export

load_image <- function (path, resize=TRUE, width) {
  imgname <- unlist(strsplit(path, "[.]")) # split at "."
  imgtype <- rev(imgname)[1] # get the last str segment after "."
  if (imgtype == "jpg") {
    require(jpeg)
    img <- readJPEG(path)
  } else if (imgtype == "tif") {
    require(tiff)
    img <- readTIFF(path)
  } else if (imgtype == "png") {
    require(png)
    img <- readPNG(path)
  } else {
    stop("file type not supported, only .jpg, .tif or .png are supported")
  }
  if (getNumberOfFrames(img) == 3) 
      img <- img[, , 1] * 0.3 + img[, , 2] * 0.59 + img[, , 3] * 0.11
  img <- t(img)
  img <- EBImage::flip(img)
  x <- dim(img)[1]
  y <- dim(img)[2]
  if (resize == TRUE) {
    if (missing(width)) {
      if ( x > 1280)
        img <- resize(img, w = 1280)
      else 
        img <- resize(img, w = x)
    } else {
      img <- resize(img, w = width)
    }
  } else {
    img <- resize(img, w=x)
  }
  invisible(img)
} 