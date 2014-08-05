#' Search new specimens' configuration against project database 
#' 
#' @description search the semi-landmarks configuration(s) of new, unknown
#'   specimens against the configurations saved in the project.
#'  
#' @details 
#'   The search is based on the Procrustes (Riemannian) distance between
#'   the query and the database. the lower the distance (\code{rdist}) between
#'   the query and the project's configuration, the higher their ranking in the
#'   search result. Perfect match will have 0 distance (for non-sliding landmark
#'   method).
#'   
#'   To shorten the time of searching, the query is first searched against the
#'   meanshape configurations of each species present in the project. For each
#'   species in the project, maximum distance between individuals and the
#'   species's meanshape was already calculated (by \code{\link{sprdist}},
#'   wrapped within \code{\link{rGPA}}). First the search will determine if the
#'   query to each species's meanshape distance is within this range, if say,
#'   the query is within the range of 2 species' range, then searching will be
#'   continued with the individuals of these 2 species only. Otherwise,
#'   searching will be continued with the individuals of the top 5 closest
#'   species only.
#'  
#'   The search has taken the semi-landmarks configuration arrangements into
#'   consideration already, hence the user doesn't need to know the side/
#'   direction of the query.
#'   
#'   \bold{Note:} sliding semi-landmark is currently not supported due to the 
#'   difficulty to achieve short searching time and correct configuration 
#'   guessing at the same time. Therefore the \code{gpa} object saved in the 
#'   project should be built with fixed semi-landmarks. However, this does not
#'   affect \code{\link{otopred}}, thus if sliding is preferred, \code{otopred}
#'   will still give good predictions.
#' 
#' @param project path to a project (\code{.rds} file) to be read, saved using 
#'   \code{\link{saveproj}}, or a project object already read into R. if not
#'   given, interactive file selector will pop out to prompt user to select a
#'   \code{.rds} file (Windows only)
#' @param query path(s) to otolith images/ path(s) of folder containing the
#'   otolith images/ \code{.tps} file containing the semi-landmark
#'   configurations/ p x k matrix or p x k x n array of semi-landmark
#'   configuration(s) to be searched. If none is given, interactive file
#'   selector will pop out to prompt user to select images to be searched
#'   (Windows only)
#' @param show integer. how many search results to be shown (in order of 
#'   ranking). Note that if number of specimens in database searched against is
#'   less than show, the result shown will be less than the given number
#' @param saveresult logical. whether to save the result
#' @param showplot logical. whether to show plot of query and their matches
#' @param savename character. optional. the file name, if \code{saveresult=TRUE}
#' @return The search result. if \code{saveresult=TRUE}, the result is written 
#'   into a \code{.txt} file. The results include:
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
# hierarchy for otosearch: otosearch > otosearch3 > otosearch2  

otosearch <- function(project, query, show = 5, saveresult = FALSE, 
                      showplot = TRUE, savename) {
  # Note: sliding is currently unsupported
  # get the project database
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
  # error message
  if (is.null(dimnames(project$landmark)[[3]]))
    stop("landmarks of project should be named")
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
    } else {
      if (file.info(query)$isdir) # ok, is the path given dir?
        pathtype <- "dir"
      else 
        pathtype <- "file"
      query <- img2landmark(query, type = pathtype, saveoutline = FALSE, 
                          savelandmark = FALSE)$landmark    
    }
  } else if (is.matrix(query) | is.array(query)) { # if not path, is it config?
    query <- query  
  } else { # not path nor config, error
    stop(paste("query should be path(in character) OR",
         "matrix/ array of semi-landmark configuration(s) OR",
         "a .tps file containing the configurations"))
  }
  # determine one specimen or multi specimens
  p <- dim(query)[1]
  if (is.matrix(query)) { # matrix, single specimen
    query <- array(query, dim = c(p, 2, 1), dimnames = list(NULL, NULL, "Query"))
  }   
  if (is.array(query)) { # array, multiple specimens
    n <- dim(query)[3]    
    result <- vector("list", n) # result shown in list for multiple queries
    if (!is.null(dimnames(query)[[3]]))
      names(result) <- dimnames(query)[[3]]
    else 
      names(result) <- paste("query", 1:n)  
    # progress bar
    cat("\n**Search in progress (", n, "queries):\n\n")
    if (n > 2)
      pb <- txtProgressBar(1, n, style = 3, char = "|")  
    # looping thru the queries
    for (i in 1:n) {
        result[[i]] <- otosearch3(specimen = query[, , i],
                                  project = project, show = show)
      # some result may be shorter than show
      if (show > dim(result[[i]])[1])
        shown <- dim(result[[i]])[1] # shown the variable to loop later
      else 
        shown <- show
      # plot
      if (showplot) {
        # determine plot number and plot area
        nplot <- shown + 1
        if (nplot <= 5) {
          rn <- 1; cn <- nplot
        } else if (nplot <= 10) {
          rn <- 2; cn <- ceiling(nplot / 2)
        } else if (nplot <= 14) {
          rn <- 3; cn <- ceiling((nplot - 4) / 2)
        } else {
          shown <- 14; rn <- 3; cn <- 5 
        }  
        dev.new(width = cn * 3, height = rn * 3) # new dev for each query
        par(mfrow = c(rn, cn), mar = c(0, 0, 0, 0))
        # plot the query first
        if (!is.null(dimnames(query)[[3]]))
          query.label <- dimnames(query)[[3]][i]
        else
          query.label <- ""
        # set the xlim and ylim of the plots
        r1 <- range(query[, 1, i])
        r2 <- range(query[, 2, i])
        r1.1<- abs(r1[2] - r1[1])
        r2.1<- abs(r2[2] - r2[1])
        xlim<- c(r1[1] - r1.1/20, r1[2] + r1.1/20)
        ylim<- c(r2[1] - r2.1/20, r2[2] + r2.1/20)
        plot(query[, , i], asp = 1, type = "n", axes = FALSE, xlim = xlim, 
             ylim = ylim)
        polygon(query[, , i], border = "blue")
        text(mean(query[,1 , i]), mean(query[,2 , i]), 
               c(paste0("Query #", i, "\n", query.label)), bty = "n")
        # insert the index into result too
        result[[i]]$index.in.db <- NA 
        # then plot the matches 1 by 1
        for (j in 1:shown) { 
          index <- which(dimnames(project$landmark)[[3]] == 
                             result[[i]]$label[j])
          result[[i]]$index.in.db[j] <- index
            plotorient <- result[[i]]$orient[j]
          match <- project$landmark[, , index]
          matchlab <- switch(plotorient, ori = "original", rev = "diff direction", 
                 flip = "diff side", fliprev = "diff direction + side")
          match <- rotonto(x = query[, , i], y = reland(match, plotorient), 
                                scale = TRUE, reflection = FALSE)$yrot
          plot(match, asp = 1, type = "n", axes = FALSE, xlim = xlim, ylim = ylim)
          polygon(match)
          text(mean(match[,1]), mean(match[,2]), 
               paste0(paste0("Match #", j), "\n", result[[i]]$label[j], 
                        paste0("\nSp:", result[[i]]$species[j]), 
                        paste0("\n[", matchlab, "]")), 
                 bty = "n")
        }
      }
      # label, for pretty presentation in txt file only
      if (saveresult) {
        result[[i]][dim(result[[i]])[1] + 1, 1] <- dimnames(query)[[3]][i]
        rownames(result[[i]])[dim(result[[i]])[1]] <- "Specimen:"
        result[[i]][dim(result[[i]])[1] +  1, ] <- NA
        result[[i]][is.na(result[[i]])] <- ""
        rownames(result[[i]])[dim(result[[i]])[1]] <- ""
      }
      if (n > 2)
        setTxtProgressBar(pb, i)
    }
  }
  if (saveresult) {
    if (missing(savename))
      savename <- paste0("otosearchresult(", Sys.Date(), ")")
    lapply(result, capture.output, file = paste0(savename, ".txt"), append = TRUE)
    cat("\n\n**The search result is saved at:", 
        paste(getwd(), paste0(savename, ".txt"), sep = "/"))
  }
  cat("\n\n")
  return(result)
}