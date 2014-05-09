#' Update the project
#' 
#' @description update the project previously saved using \code{\link{saveproj}}
#' @param name character. file name of the project to be updated.
#' @inheritParams saveproj 
#' @return the selected values are updated (overwrited) with the new addtions.
#' @details the objects assigned with new values will be overwrited with new ones,
#'  the rest will remain the same.
#' @note the current version only support overwrite(!), not append.  
#' @seealso 
#'  \code{\link{saveproj}}, \code{\link{loadproj}}
#' @export

updateproj <- function(name, landmark=NULL, des=NULL, gpa=NULL, nef=NULL, pc=NULL, 
                      har=NULL, des.threshold=NULL, nef.threshold=NULL, 
                      gpa.threshold=NULL) {
  temp <- readRDS(name)
  if (is.null(landmark)) landmark <- temp$landmark
  if (is.null(des)) des <- temp$des
  if (is.null(gpa)) gpa <- temp$gpa
  if (is.null(nef)) nef <- temp$nef
  if (is.null(pc)) pc <- temp$pc
  if (is.null(har)) har <- temp$har
  if (is.null(des.threshold)) des.threshold <- temp$des.threshold
  if (is.null(nef.threshold)) nef.threshold <- temp$nef.threshold
  if (is.null(gpa.threshold)) gpa.threshold <- temp$gpa.threshold
  temp2 <- list (landmark=landmark, des=des, nef=nef, gpa=gpa, pc=pc, 
           har=har, class=des$sp)
  saveRDS(temp2, file=paste0(name))
  cat("The updated project is saved at:", 
      paste(getwd(), name, sep="/"),"\n\n")
}