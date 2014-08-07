#' Process shape descriptor files
#' 
#' @description A wrapper function to do the routine of reading the files from 
#'  imagej output (from \code{\link{runmacro}}) and combine the descriptor files 
#'  into a single file. 
#' @param imagejdir path of the folder containing the shape descriptor files. 
#'  if not given, a GUI folder selector will prompt user to select (if run in Windows)
#' @param write logical. whether to save a new file of combined descriptor.
#' @param plot logical. whether to plot PCA for preliminary assessment of 
#'  descriptor data. Only graphical output of first 3 PCs.
#' @param label logical. Whether to extract label from the image file name. 
#' @param extract a numeric vector of two number telling which characters to extract
#'  from the file name to become the species name, if \code{label=TRUE}. The first
#'  is the starting position, second is the last position of character to extract. 
#' @return A dataframe of combined shape descriptor data
#' @seealso \code{\link{runmacro}}
#' @export

routine1 <- function(imagejdir, write=TRUE, plot=FALSE, label=TRUE, 
                     extract=c(1, 6)) {
  if (missing(imagejdir)) {
    if (Sys.info()['sysname'] == "Windows")
      imagejdir <- choose.dir(default = getwd(), 
                   caption = "Select folder containing the shapedescriptors.txt files")
    else
      stop("you need to provide path to folder containing the shape descriptor files")
  }
  deslist <- list.files(imagejdir, pattern= "_shapedescriptors.txt")  
  nlist <- gsub("_shapedescriptors.txt", "", deslist) 
  # create shapedescriptor dataframe
  desdata <- NULL
  cat("------Reading shape descriptor files------\n")
  for (i in 1:length(deslist)) {
    temp <- read.table(paste(imagejdir, 
                             deslist[i], sep="/"), comment.char="")
    desdata <- rbind(desdata, temp)
    cat(paste(i, "  ", nlist[i], "\n"))
  }
  cat("---Reading shape descriptor files ended---\n")  
  # write the combined files into a new csv file
  desdatan <- desdata[, -c(2:22, 23:33, 37:39, 42, 48)]
  rownames(desdatan) <- desdata[, "Label"]
  if (label) {
    if (length(extract) > 2 || !is.numeric(extract))
      stop("please provide two number for file name extraction only")
    sp <- factor(substr(desdatan$Label, extract[1], extract[2]))
    desdatan$sp <- sp
  }
  if (write == TRUE) {
    filename1 <- paste0("combined_descriptors_complete(", Sys.Date(), ").csv")
    filename2 <- paste0("combined_descriptors(", Sys.Date(), ").csv")
    write.csv(desdata, file=filename1, row.names=FALSE)
    write.csv(desdatan, file=filename2, row.names=FALSE)
    cat("\nThe combined file is saved at:\n1) ", 
        paste(getwd(), filename1, sep="/"),
        "\n2) ", paste(getwd(), filename2, sep="/"), "\n\n")
  }
  if (plot == TRUE) {
    pca <- prcomp(desdatan[, 2:13], scale.=TRUE)
    plotpca(pca, class=desdatan$sp)
  } 
  invisible(desdatan)
}