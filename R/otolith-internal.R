#' @keywords internal 

# function to get extension from a list of full/ half path
# how it works: split character at ".", take the last splitted part as ext
# fixed bug for file name with "." in it 
# using vapply to avoid loop and specify output type
# @param pathlist list of path of length >=1
.ext <- function (pathlist) {
  vapply(pathlist, function(x) rev(unlist(strsplit(x, "[.]")))[1], 
         character(1), USE.NAMES = FALSE)
}

# function to get file names from a list of full path
# how it works: split character at "//" OR "/" (regular expression accepts logical
#  conditions), take the last part, and if remove.ext=TRUE, substitute the 
#  ".<ext>" with empty string
# @param pathlist list of path of length >=1
# @param remove.ext specify whether to remove the extension
# @reference http://stackoverflow.com/questions/19424709/r-gsub-pattern-vector-and-replacement-vector
.getfilename <- function (pathlist, remove.ext = TRUE) {
  output <- vapply(pathlist, function(x) rev(unlist(strsplit(x, 
            split = "[\\]|/")))[1], character(1), USE.NAMES = FALSE)
  if(remove.ext)
    output <- mapply(gsub, paste0(".", .ext(output)), "", output, 
                     USE.NAMES = FALSE)
  return(output)
}

# function to get directory of where the file is from a full path name/ list
# @param path a single path name or a path list
# @note if a path list is provided, extraction is based on the first file on 
#   the list only
# @return path
.getdir <- function(path) {  
  if (length(path) > 1) path <- path[1]
  paste(rev(rev(unlist(strsplit(path, "[\\]|/")))[-1]), collapse="/") 
}

# function to determine if an object is a otolith project
# @param x object to be determined
# @return logical
.isproj <- function (x) {
  temp <- c("landmark","des", "nef", "gpa", "pc", "har", "class") %in% unlist(attributes(x))
  !any(temp == FALSE)
}