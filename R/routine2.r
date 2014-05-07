#' Process outline and sample semi-landmarks
#' 
#' @description A wrapper function to do the routine of reading the files from 
#'  imagej output, then read the outline files, checking for pixel artifacts, 
#'  re-order the outline coordinates and laying curvilinear semi-landmarks
#' @param imagejfolder Path of the folder containg the outline xy coordinate files
#' @param rfolder Path of the output folder (must exist already)
#' @param write Whether to save new files of record of outline with defects, 
#'  re-ordered outline coordinates, curvilinear semi-landmarks and a combined 
#'  file curvilinear semilandmarks of all configuration in TPS format. 
#'  TPS format is compatible with many other morphometrics software.
#' @param plot Logical. Whether to save the plots of outline with semi-landmarks.
#' @param nd Numeric. Number of sampling interval. See \code{\link{equaldist}} function. 
#' @param scale scale in the unit of um/px, 
#'  i.e. 1 pixel of the image equal to how many um?
#' @return 
#'  \item{rmindex}{index of the outline with pixel artifacts}
#'  \item{rmfile}{file label of the outline with pixel artifacts}  
#'  \item{nlist}{list of file labels after checking}
#'  \item{A}{array of semi-landmarks configurations}
#' @importFrom geomorph writeland.tps
#' @seealso 
#'  similar: \code{\link{routine1}}
#'  
#'  which this function wraps: \code{\link{orderOutline} \link{equaldist}}
#' @export

routine2 <- function(imagejfolder, rfolder, write=TRUE, 
                     plot=TRUE, nd=50, scale=NULL) {
  # check directory validity first
  setwd(rfolder) 
  setwd(imagejfolder)
  coordlist <- list.files(pattern="_xycoords.txt")
  nlist <- gsub("_xycoords.txt", "", coordlist) 
  coorddatalist <- vector("list", length(coordlist)) 
  names(coorddatalist) <- nlist 
  #reading the outline data from imagej
  cat("------Reading outline coordinate files------\n")
  for (i in 1:length(nlist)) {
    coorddatalist[[i]] <- read.table(paste(imagejfolder,
                          coordlist[i], sep="/"))
    # remove third column as the output in this column from ImageJ 
    coorddatalist[[i]] <- coorddatalist[[i]][, -3] 
    cat(paste(i, nlist[i], "\n"))
    flush.console()
  }
  cat("---Reading outline coordinate files ended---\n\n")
  #check for artefacts
  cat("------Checking for pixel artifacts------\n")
  checkresult <- NULL
  for(i in 1:length(nlist)) {
    cat(paste("Image (", i, "/", length(nlist), ")\nFile Name: ",
              nlist[i], "\nChecking for artifact..."))
    flush.console()
    temp <- checkart(coorddatalist[[i]])
    cat("done\n\n")
    checkresult[i] <- temp$artifact
  }
  cat("---Checking for pixel artifacts ended---\n\n")
  #remove the problematic ones
  if (length(which(checkresult == TRUE)) > 0) {
    setwd(rfolder)
    rmindex <- which(checkresult== TRUE)
    rmfile <- nlist[rmindex]
    if (write == TRUE)
      write.table(cbind(rmindex, rmfile), file="defect.txt", 
                  quote=FALSE, row.names=FALSE, sep="\t")
    nlistN <- nlist[-rmindex]
    rmlistR <- rev(rmindex)
    for (j in 1:length(rmindex)) {
      coorddatalist[[(rmindex[j])]] <- NULL
    }
  }
  if (length(which(checkresult == T)) == 0) {
    rmindex <- NULL
    rmfile <- NULL
    nlistN <- nlist
  }
  # set the file name pattern for 3 types of output
  savenames <- paste0("curvillinearLand_",nlistN,".tif")
  savename2 <- paste0("orderedLand_",nlistN,".csv")
  savename3 <- paste0("curvillinearLand_",nlistN,".csv")
  # order outline batch processing + semilandmarks laying
  require(secr)
  A <- array(data=NA, dim=c(nd * 2, 2, length(nlistN)), 
             dimnames=list(paste0("LM", 1:(nd * 2)), c("x", "y"), nlistN))
  for (i in 1:length(nlistN)) {
    cat(paste("Image (", i, "/", length(nlistN), ")\nFile Name: ", 
        nlistN[i], "\nRe-order outline coords sequences..."))
    flush.console()
    # order outline
    coorddatalist[[i]] <- orderOutline(coorddatalist[[i]])
    cat("done\nLaying semilandmarks...")
    flush.console()
    # semilandmarks
    if (plot == TRUE) {
      setwd(rfolder)
      tiff(savenames[i],width=800, height=600,res = 144)
      par(mar=c(rep(0,4)))
      landtemp <- equaldist(coorddatalist[[i]], 
                  plot=plot, label=T, nd=nd)
      dev.off()
    } else {
      landtemp <- equaldist(coorddatalist[[i]], plot=plot, label=T, nd=nd)
    }
    A[, , i] <- landtemp
    if (write == TRUE) {
      setwd(rfolder)
      write.csv(coorddatalist[[i]], file=savename2[i],  row.names=FALSE)
      write.csv(landtemp, file=savename3[i], row.names=FALSE)
    }
    cat("done\n\n")
    flush.console()
  }
  if (!is.null(scale))
    A <- A * scale
  if (write == TRUE) {
    writeland.tps(A, file= "semilandmarks.tps")
    cat("\nThe sampled semilandmark file is saved at:",  
        paste(rfolder, "semilandmarks.tps", sep="/"), "\n")
  }
  invisible(list(rmindex=rmindex, rmfile=rmfile, nlist=nlistN, A=A))
}