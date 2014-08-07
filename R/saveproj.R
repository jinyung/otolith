#' Save a project
#' 
#' @description this function save the processed semi-landmarks, shape descriptor 
#'  data, GPA/ EFA transformed data, optimized parameters into a \code{.rds} file.
#' @param name character or path. file name/ path of the project, should end with \code{".rds"}
#' @param landmark p x k x n array of semi-landmarks, e.g. \code{landmark} value of 
#'  \code{\link{img2landmark}}.
#' @param des a \code{\link{routine1}} object 
#' @param nef a \code{\link{rNEF}} object
#' @param gpa a \code{\link{rGPA}} object
#' @param pc numeric. number of optimized PC range determined by \code{\link{pccv}}
#' @param har numeric. number of optimized harmonics range determined by 
#'  \code{\link{harcv}}
#' @return a \code{.rds} file saved
#' @Note The saved project could be read into R using \code{\link{readRDS}} 
#'  function.
#' @seealso 
#'  Similar: \code{\link{readRDS}}, \code{\link{updateproj}}
#'  
#'  Which this function wraps: \code{\link{saveRDS}}
#' @export

saveproj <- function(name, landmark = NULL, des = NULL, gpa = NULL, nef = NULL, 
                     pc = NULL, har = NULL, class = NULL) {
  if (is.null(class) & !is.null(des$sp))
    class <- des$sp
  temp <- list(landmark = landmark, des = des, nef = nef, gpa = gpa, pc = pc, 
               har = har, class = class, name=gsub(".rds", "", name))
  saveRDS(temp, file=name)
  cat("The project is saved at:", 
      paste(getwd(), name, sep="/"), "\n\n")
}