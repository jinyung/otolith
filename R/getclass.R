#' Extract class name from landmark configurations
#' @description A simple wrapper to extract class name, essentially a 
#'  wrapped \code{\link{substr}} function
#' @param landmark p x k x n array of landmark. n dimension has to be named.
#' @param extract numeric vector of 2 values, with the first giving the starting
#'  character for extraction, and second for ending character for extraction
#' @return class name
#' @export

getclass <- function(landmark, extract) {
  if (is.null(dimnames(landmark)[[3]]))
    stop("landmark array do not have names")
  if (!is.numeric(extract) & length(extract == 2))
    stop("extract should be a numeric vector consist of 2 values, first the start and second the end character to extract")
  factor(substr(dimnames(landmark)[[3]], extract[1], extract[2]))
}