#' Sample 'equally spaced'-style curvilinear semi-landmarks
#' 
#' @description Sample points along outline. Outline is divided into two parts 
#'  by larget PC axis, semi-landmarks are then sampled in equal intervals in 
#'  each part. 
#' @param outdata outline data. A matrix consists of xy coordinates
#' @param nd numeric. Number of equally spaced divides along the outline for each part of the outline. 
#' @param plot logical. Whether to plot
#' @param label logical. Whether to label the plot
#' @return semi-landmarks sampled and the plot of the sampled 
#'  semi-landmarks on the outline if \code{plot=TRUE} 
#' @seealso \code{\link{routine2}}
#' @export

equaldist <- function (outdata,  nd=50,  plot=FALSE,  label=FALSE){
  # aligne the configuration to largest axis 
  alout <- aligne2(outdata)
  index <- c(which(alout[,1] == min(alout[, 1])), 
             which(alout[, 1] == max(alout[, 1])))
  # get the end points
  endpoints <- outdata[c(index), ]
  if (index[1] < index[2]){
    whofirst <-1;whonext <-2 
  } else if (index[2] < index[1]){
    whofirst <-2;whonext <-1} 
  toplandindex <- round(1:(nd - 1) * abs(
                  index[1]-index[2])/nd + index[whofirst])
  botlandindex <- round(1:(nd - 1) * (dim(outdata)[1] 
                  - index[whonext] + index[whofirst]) / nd + index[whonext])
  botlandindex <- c(botlandindex[botlandindex <= dim(
                  outdata)[1]], botlandindex[botlandindex > dim(
                  outdata)[1]] - dim(outdata)[1])
  curviland <- rbind(endpoints[whofirst, ], 
               outdata[toplandindex, ], endpoints[whonext, ], 
               outdata[botlandindex, ])
  # plot the result
  if (plot) {
    par(mar=c(0, 0, 0, 0))  
    plot(outdata, type='n',asp=1)
    polygon(outdata)
    points(endpoints, pch=21, bg="red")
    lines(endpoints, col="red",lty=3)
    points(outdata[toplandindex, ], pch=21, bg='blue') 
    points(outdata[botlandindex, ], pch=21, bg='blue')    
    if (label) {
      text(outdata[toplandindex, ], pos=3,  label=2:nd,  cex=0.5)
      text(outdata[botlandindex, ], pos=1,
           label=(1:(nd - 1)) + nd + 1, cex=0.5)
      text(outdata[index[c(whofirst, whonext)], ], 
           labels=c(1, nd + 1), pos=c(2, 4), cex=0.5)
    }
  }
  return(curviland)
}