#' Load image
#' 
#' @description To read image into R 
#' @param path path of a single image
#' @param resize logical. whether to resize the picture
#' @param width width to resize to, in pixel. used only if \code{resize=TRUE}.
#'  if \code{resize=TRUE} and width is not given, image is resized to 
#'  1280px in width by default
#' @return an image matrix
#' @importFrom jpeg readJPEG
#' @importFrom tiff readTIFF
#' @importFrom png readPNG
#' @importFrom EBImage resize
#' @details All images loaded are converted into gray scale. Image format is 
#'  automatically detected from the extension. Supported formats include \code{.jpg}, 
#'  \code{.tif}, and \code{.png}. This function is used in place of 
#'  \code{\link[EBImage]{readImage}} (from \code{EBImage}) for this package to 
#'  speed up the reading of image.
#' @seealso 
#'  Which this function wraps: \code{\link[jpeg]{readJPEG}}, 
#'  \code{\link[tiff]{readTIFF}}, \code{\link[png]{readPNG}} 
#' @export

loadimg <- function (path, resize=TRUE, width) {
#   imgtype <- .ext(path) # get the last str segment after "."
  if (grepl("jpe?g", path, ignore.case = TRUE)) {
    # require(jpeg)
    img <- jpeg::readJPEG(path)
  } else if (grepl("tif?f", path, ignore.case = TRUE)) {
    # require(tiff)
    img <- tiff::readTIFF(path)
  } else if (grepl("png", path, ignore.case = TRUE)) {
    # require(png)
    img <- png::readPNG(path)
  } else {
    stop("file type not supported, only .jpg, .tif or .png are supported")
  }
  # require(EBImage)
  if (EBImage::numberOfFrames(img) == 3) 
      img <- img[, , 1] * 0.3 + img[, , 2] * 0.59 + img[, , 3] * 0.11
  img <- t(img)
  img <- EBImage::flip(img)
  x <- dim(img)[1]
  y <- dim(img)[2]
  if (resize == TRUE) {
    # require(EBImage)
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