#' Sample 'comb'-style curvilinear semilandmarks 
#' 
#' @description Sample semi-landmarks on outline using the 'comb' method.
#' @inheritParams equaldist
#' @return semi-landmarks sampled and the plot of the sampled 
#'  semi-landmarks on the outline if \code{plot=TRUE} 
#' @seealso \code{\link{equaldist}}
#' @export
#' @references
#' Ponton, D. (2006). Is geometric morphometrics efficient for comparing 
#' otolith shape of different fish species?. Journal of Morphology, 
#' 267(6), 750-757.
#' 
#' Sheets, H. D., Covino, K. M., Panasiewicz, J. M., & Morris, S. R. (2006). 
#' Comparison of geometric morphometric outline methods in the discrimination 
#' of age-related differences in feather shape. Frontiers in Zoology, 3(1), 1-12.
#' 
#' @note produce errorneous results for complex outline with  high outline curvature 
#'  especially at the 'excisura' site.

comb <- function (outdata, nd=50, plot=FALSE, label=FALSE) {
  # get the end points
  alout <- aligne2(outdata)
  index <- c(which(alout[,1] == min(alout[, 1])), 
             which(alout[, 1] == max(alout[, 1])))
  endpoints <- outdata[c(index), ]  
  x1 <- endpoints[1, 1]
  x2 <- endpoints[2, 1]
  y1 <- endpoints[1, 2]
  y2 <- endpoints[2, 2]
  # get the slope and intercept
  slope <- (y2 - y1) / (x2 - x1) #slope 
  ceptf <- y1 - slope * x1 #intercept 
  # get the dividing points
  divide.x <- NULL
  divide.y <- NULL
  for (j in 1:(nd - 1)) {
    if (slope < 0) {
      divide.x[j] <- max(x2, x1) - j * abs((x2 - x1) / nd)
    }
    else if (slope > 0) {
      divide.x[j] <- min(x2, x1) + j * abs((x2 - x1) / nd)
    }
    divide.y[j] <- min(y2, y1) + j * abs((y2 - y1) / nd)
  }
  # dividing the outline into two parts
  yf <- slope * outdata[, 1] + ceptf
  topindex <- which(outdata[, 2] > yf)
  bottomindex <- which (outdata[, 2] < yf)
  outdataTop <- outdata[topindex, ]
  outdataBottom <- outdata[bottomindex, ]
  # searching the points on outline 
  combPointIndexTop <- NULL
  combPointIndexBottom <- NULL  
  for (k in 1:(nd - 1)) {
    slopePTop <- (outdataTop[, 2] - divide.y[k]) / (outdataTop[, 1] 
                 - divide.x[k])
    slopedifftop <- abs(slopePTop - (-1 / slope))
    slopePBottom <- (outdataBottom[, 2] - divide.y[k]) / (outdataBottom[, 1] 
                     - divide.x[k])
    slopediffbottom <- abs(slopePBottom - (-1/slope))
    combPointIndexTop[k] <- topindex[which(slopedifftop == min(slopedifftop))]
    combPointIndexBottom[k] <- bottomindex[which(slopediffbottom == 
                               min(slopediffbottom))]
  }
  # organize the semi-landmarks correctly
  topsemiland <- outdata[combPointIndexTop, ]
  botsemiland <- outdata[combPointIndexBottom, ]
  botsemilandR <- botsemiland[nrow(botsemiland):1, ] # reverse the row of the matrix
  landmark <- rbind(endpoints[1, ], topsemiland, endpoints[2, ], botsemilandR)
  # plotting output
  if (plot) {
    par(mar=c(0, 0, 0, 0))  
    plot(outdata, type="n", asp=1)
    polygon(outdata)
    points(endpoints, pch=21, bg="red")
    lines(endpoints, col="red", lty=3)
    points(rbind(topsemiland, botsemiland), pch=21, bg="blue")
    for (i in 1:(nd-1))
      lines(rbind(topsemiland[i, ], botsemiland[i, ]),  lty=3, col="blue")
    if (label) {
      text(topsemiland, labels=2:nd,  pos=3,  cex=0.5)
      text(botsemilandR, labels=(1:(nd-1)) + nd + 1, pos=1, cex=0.5)
      text(endpoints, labels=c(1, nd+1), pos=c(2, 4), cex=0.5)
    }
  }
  return(landmark)
}