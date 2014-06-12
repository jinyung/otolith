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
#' @param project path to a project (\code{.rds} file) to be read, saved using 
#'  \code{\link{saveproj}}, or a project object already read into R.
#'  if not given, interactive file selector will pop out to prompt user to select 
#'  a \code{.rds} file (Windows only)
#' @param query path(s) to otolith images/ path(s) of folder containing the otolith 
#'  images/ code{.tps} file containing the semi-landmark configurations/ p x k matrix 
#'  or p x k x n array of semi-landmark configuration(s) to be searched. If none is 
#'  given, interactive file selector will pop out to prompt user to select images 
#'  to be searched (Windows only)
#' @param show numeric. how many search results to be shown (in order of ranking).
#' @param saveresult logical. whether to save the result.
#' @param name character. optional. the file name, if \code{saveresult=TRUE}. 
#' @return The search result. if \code{saveresult=TRUE}, the result is written into a 
#'  \code{.txt} file. The results include:
#'    \itemize{ 
#'    \item \code{rdist}: distance between the query and the match 
#'    \item \code{label}: label of the match from the database
#'    \item \code{species}: species of the match
#'    \item \code{inside}: whether the \code{rdist} is within the range of that 
#'      species in the database
#'    \item \code{orient}: type of configuration that match 
#'    }
#' @importFrom geomorph readland.tps
#' @seealso 
#'  Similar: \code{\link{otopred}}
#' @export

otosearch <- function(project, query, show = 5, saveresult = FALSE, name) {
  # get the project databse
  if (missing(project)) {
    if (Sys.info()['sysname'] == "Windows") 
    project <- readRDS(choose.files(getwd(), caption = 
               "Select project containing the database", multi = FALSE, 
               filters = c("Project (*.rds)", "*.rds; *.RDS")))
    else 
      stop("Please provide the project")
  } else {
    if (.isproj(project))
      project <- project
    else if (is.character(project)) {
      if (file.exists(project) & .ext(project) == "rds" | .ext(project) == "RDS")
        project <- readRDS(project)
      else
        stop("Please provide right format of project")
    } else 
      stop("Please provide right format of project")
  }
  # get the query/ determine the query type
  if (missing(query)) { # interactive selection if path not specified
    query <- img2landmark(type = "file", saveoutline = FALSE, 
                          savelandmark = FALSE)$landmark  
  } else if (is.character(query)) { # is it even a character (path)?
    if (!any(file.exists(query))) # ok, so is it a legit path?
      stop("path(s) given do not exist, please check again")
    if (any(.ext(query) == "tps")) {
      if (length(query) > 1)
        stop("multiple .tps files are not allowed")
      require(geomorph)
      query <- readland.tps(query, specID="ID")
    }
    if (file.info(query)$isdir) # ok, is the path given dir?
      pathtype <- "dir"
    else 
      pathtype <- "file"
    query <- img2landmark(query, type = pathtype, saveoutline = FALSE, 
                          savelandmark = FALSE)$landmark    
  } else if (is.matrix(query) | is.array(query)) { # if not path, is it config?
    query <- query  
  } else { # not path nor config, error
    stop(paste("query should be path(in character) OR",
         "matrix/ array of semi-landmark configuration(s) OR",
         "a .tps file containing the configurations"))
  }
  # determine one specimen or multi specimens
  if (is.matrix(query)) # matrix, single specimen
    result <- otosearch3(specimen = query, project = project, show = show)
  else if (is.array(query)) { # array, multiple specimens
    n <- dim(query)[3]
    if (n == 1) {
      result <- list(otosearch3(specimen = query[, , 1], project = project, 
                           show = show))
      names(result) <- "query"
    } else {
      result <- vector("list", n) # result shown in list for multiple queries
      if (!is.null(dimnames(query)[[3]]))
        names(result) <- dimnames(query)[[3]]
      else 
        names(result) <- paste("query", 1:n)  
      # progress bar
      cat("\n**Search in progress (", dim(query)[3], "queries):\n\n")  
      pb <- txtProgressBar(1, n, style = 3, char = "|")  
      for (i in 1:n) {
        result[[i]] <- otosearch3(specimen = query[, , i], 
                                  project = project, show = show)
        # label, for pretty presentation in txt file only
        if (saveresult) {
          result[[i]][dim(result[[i]])[1] + 1, 1] <- dimnames(query)[[3]][i]
          rownames(result[[i]])[dim(result[[i]])[1]] <- "Specimen:"
          result[[i]][dim(result[[i]])[1] +  1, ] <- NA
          result[[i]][is.na(result[[i]])] <- ""
          rownames(result[[i]])[dim(result[[i]])[1]] <- ""
        }
        setTxtProgressBar(pb, i)
      }
    }
  }
  if (saveresult) {
    if (missing(name))
      name <- paste0("otosearchresult(", Sys.Date(), ")")
    lapply(result, capture.output, file = paste0(name, ".txt"), append = TRUE)
    cat("\n\n**The search result is saved at:", 
        paste(getwd(), paste0(name, ".txt"), sep = "/"))
  }
  cat("\n\n")
  return(result)
}