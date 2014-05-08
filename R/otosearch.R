#' Search new specimens' configuration against project database 
#' 
#' @description search the semi-landmarks configuration(s) of new, unknown specimens
#'  against the configurations saved in the project.
#'  
#' @details the search is based on the Procrustes (Riemannian) distance between the 
#'  query and the database. the lower the distance (\code{rdist}) between the query and the project's 
#'  configuration, the higher their ranking in the search result. Perfect match will
#'  have 0 distance.
#'  
#'  To shorten the time of searching, the query is first searched against the meanshape
#'  configurations of each species present in the project. For each species in the project,
#'  maximum distance between individuals and the species's meanshape was already
#'  calculated (by \code{\link{sprdist}}, wrapped within \code{\link{rGPA}}). 
#'  First the search will determine if the query to each species's meanshape 
#'  distance is within this range, if say, the query is within the range of 
#'  2 species' range, then searching will be continued with the individuals 
#'  of these 2 species only. Otherwise, searching will be continued with 
#'  the individuals of the top 5 closest species only. 
#'  
#'  The search has taken the semi-landmarks configuration arrangements into consideration
#'  already, hence the user doesn't need to know the side/ direction of the query.
#'   
#' @param query p x k matrix or p x k x n array. the query configuration(s) 
#' @param project a project saved using \code{\link{saveproj}}. query is searched 
#'  against the whole array of configurations saved within the project.
#' @param show numeric. how many search results to be shown (in order of ranking).
#' @param save logical. whether to save the result.
#' @param name character. optional. the file name, if \code{save=TRUE}. 
#' @return The search result. if \code{save=TRUE}, the result is written into a 
#'  \code{.txt} file
#' @seealso 
#'  Similar: \code{\link{otopred}}
#'  
#'  Which this function wraps: \code{\link{otosearch3}}
#' @export

otosearch <- function(query, project, show=5, save=FALSE, name=NULL) {
  if (is.matrix(query)) # matrix, single specimen
    result <- otosearch3(specimen=query, project=project, show=show)
  else if(is.array(query)){ # array, multiple specimens
    n <- dim(query)[3]
    result <- vector("list", n) # result shown in list for multiple queries
    if (!is.null(dimnames(query)[[3]]))
      names(result) <- dimnames(query)[[3]]
    else 
      names(result) <- paste("Specimen", 1:n)  
    # progress bar
    cat("\n**Searching in progress (", dim(query)[3], "queries):\n\n")  
    pb <- txtProgressBar(1, n, style=3, char="|")  
    for (i in 1:n) {
      result[[i]] <- otosearch3(specimen=query[, , i], project=project, show=show)
      #label, for pretty presentation in txt file
      result[[i]][dim(result[[i]])[1] + 1, 1] <- dimnames(query)[[3]][i]
      rownames(result[[i]])[dim(result[[i]])[1]] <- "Specimen:"
      result[[i]][dim(result[[i]])[1] +  1, ] <- NA
      result[[i]][is.na(result[[i]])] <- ""
      rownames(result[[i]])[dim(result[[i]])[1]] <- ""
      setTxtProgressBar(pb, i)
    }
  }
  if (save) {
    if (is.null(name))
      name <- paste0("otosearchresult(", Sys.Date(), ")")
    lapply(result, capture.output, file=paste0(name, ".txt"), append=TRUE)
    cat("\n\n**The search result is saved at:", 
        paste(getwd(), paste0(name, ".txt"), sep="/"))
  }
  cat("\n\n")
  return(result)
}