#' Otosearch - internal function
#' 
#' @description search one specimen against whole database, with exhaustive 
#'  search using all possible semi-landmarks configuration arrangements
#' @param specimen matrix of configuration of query 
#' @param project a project object, see \code{\link{saveproj}}
#' @param show how many search results to be shown (in order of ranking)
#' @return search result
#' @importFrom Morpho rotonto kendalldist bindArr
#' @seealso 
#'  Functions that wraps this function: \code{\link{otosearch}}
#'  
#'  Which this function wraps: \code{\link{otosearch2}}, \code{\link{reland}}
#' @keywords internal
#' @export

otosearch3 <- function(specimen, project, show = 5) {
  # require(Morpho)
  meanshape <- project$gpa$meanshape
  # first get specimen to (by-species) meanshape dist
  ms.rdist <- numeric() 
  mslevel <- dim(meanshape)[3]  
  for (i in 1:mslevel) {  
    temp <- Morpho::rotonto(x = specimen, y = meanshape[, , i], 
                            reflection = FALSE, scale = TRUE)
    ms.rdist[i] <- Morpho::kendalldist(temp$X, temp$Y)
    temp2 <- Morpho::rotonto(x = reland(specimen, "rev"), y = meanshape[, , i], 
                             reflection = FALSE, scale = TRUE)
    ms.rdist[mslevel + i] <- Morpho::kendalldist(temp2$X, temp2$Y)
    temp3 <- Morpho::rotonto(x = reland(specimen, "flip"), y = meanshape[, , i], 
                             reflection = FALSE, scale = TRUE)
    ms.rdist[mslevel * 2 + i] <- Morpho::kendalldist(temp3$X, temp3$Y)
    temp4 <- Morpho::rotonto(x = reland(specimen, "fliprev"), 
                             y = meanshape[, , i], 
                             reflection = FALSE, scale = TRUE)
    ms.rdist[mslevel * 3 + i] <- Morpho::kendalldist(temp4$X, temp4$Y)
  }
  distmax <- project$gpa$rdist[, 3]
  # initialize
  result <- NULL
  # see if any specimen to meanshape dist less than the max dist to 
  # (contd) meanshape present in the databse
  if (any(ms.rdist < distmax)) {
    # if yes search only these 
    index <- which(ms.rdist < distmax)
     for (i in 1:length(index)) {
      if (index[i] <= mslevel) {
        # first the original configuration
        clevel <- levels(project$class)[index[i]]
        Aindex <- which(project$class == clevel)
        tresult <- otosearch2(specimen, project$landmark[, , Aindex], 
                   show=length(Aindex))
        tresult$species <- rep(clevel, length(Aindex))        
        tresult$inside <- rep(TRUE, length(Aindex))
        tresult$orient <- rep("ori", length(Aindex))
        result <- rbind(result, tresult)
      } else if (index[i] <= mslevel * 2) {
        # then the reversed configuration
        clevel <- levels(project$class)[index[i] - mslevel]
        Aindex <- which(project$class == clevel)
        tresult <- otosearch2(reland(specimen, "rev"), 
                   project$landmark[, , Aindex], show = length(Aindex))
        tresult$species <- rep(clevel, length(Aindex))
        tresult$inside <- rep(TRUE, length(Aindex))
        tresult$orient <- rep("rev", length(Aindex))
        result <- rbind(result, tresult)
      } else if (index[i] <= mslevel * 3) {
        # then the flipped configuration
        clevel <- levels(project$class)[index[i] - mslevel * 2]
        Aindex <- which(project$class == clevel)
        tresult <- otosearch2(reland(specimen, "flip"),
                   project$landmark[, , Aindex], show = length(Aindex))
        tresult$species <- rep(clevel, length(Aindex))
        tresult$inside <- rep(TRUE, length(Aindex))
        tresult$orient <- rep("flip", length(Aindex))
        result <- rbind(result, tresult)
      } else {
        # last the flipped and reversed configuration
        clevel <- levels(project$class)[index[i] - mslevel * 3]
        Aindex <- which(project$class == clevel)
        tresult <- otosearch2(reland(specimen, "fliprev"), 
                   project$landmark[, , Aindex], show = length(Aindex))
        tresult$species <- rep(clevel, length(Aindex))
        tresult$inside <- rep(TRUE, length(Aindex))
        tresult$orient <- rep("fliprev", length(Aindex))
        result <- rbind(result,tresult)
      }    
    }
  } else {
    # if not search the closest five species only
    # below is the same as above 
    index <- order(ms.rdist - distmax)[1:5]
    for (i in 1:5) {
      if (index[i] <= mslevel) {
        clevel <- levels(project$class)[index[i]]
        Aindex <- which(project$class == clevel)
        tresult <- otosearch2(specimen, project$landmark[, , Aindex], 
                   show=length(Aindex))
        tresult$species <- rep(clevel, length(Aindex))
        tresult$inside <- rep(FALSE, length(Aindex))
        tresult$orient <- rep("ori", length(Aindex))
        result <- rbind(result, tresult)
      } else if (index[i] <= mslevel * 2) {
        clevel <- levels(project$class)[index[i] - mslevel]
        Aindex <- which(project$class == clevel)
        tresult <- otosearch2(reland(specimen, "rev"),
                   project$landmark[, , Aindex], show = length(Aindex))
        tresult$species <- rep(clevel, length(Aindex))
        tresult$inside <- rep(FALSE, length(Aindex))
        tresult$orient <- rep("rev", length(Aindex))
        result <- rbind(result, tresult)
      } else if (index[i] <= mslevel * 3) {
        clevel <- levels(project$class)[index[i] - mslevel * 2]
        Aindex <-which(project$class == clevel)
        tresult <- otosearch2(reland(specimen, "flip"),
                   project$landmark[, , Aindex], show = length(Aindex))
        tresult$species <- rep(clevel, length(Aindex))
        tresult$inside <- rep(FALSE, length(Aindex))
        tresult$orient <- rep("flip", length(Aindex))
        result <- rbind(result, tresult)
      } else {
        clevel <- levels(project$class)[index[i] - mslevel * 3]
        Aindex <- which(project$class == clevel)
        tresult <- otosearch2(reland(specimen, "fliprev"), 
                   project$landmark[, , Aindex], show = length(Aindex))
        tresult$species <- rep(clevel, length(Aindex))
        tresult$inside <- rep(FALSE, length(Aindex))
        tresult$orient <- rep("fliprev", length(Aindex))
        result <- rbind(result, tresult)
      }
    }
  }
  if (show > dim(result)[1])
    show <- dim(result)[1]  # so that show will always < results
  result.index <- order(result[, 1])[1:show]
  result <- result[result.index, ]
  rownames(result) <- paste0("rank", 1:dim(result)[1])
  return(result)
}