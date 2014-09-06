#' Predict new samples based on classification models
#' 
#' @description Predict new, unknown samples using either shape decriptors or 
#'  semi-landmark configurations.
#' @details The functions include the \code{\link{otosearch}} algorithm as a
#'   means to guess the semi-landmarks arragement (the four types of
#'   arrangements, see \code{\link{reland}}), so that the user can predict the
#'   new, unknown samples even if the side and the direction of the otoliths are
#'   unknown. This could be turned off using \code{reland = FALSE} to speed up the
#'   prediction if the side and direction of the query is known, and is already
#'   in the same arrangements as the dataset in project.
#'   
#'   Combination of multiple views of otoliths for prediction can be used by 
#'   setting \code{multiview = TRUE}. Under \emph{multiview} mode, input of
#'   multiple \code{project} and \code{query} can be done by using list. For 
#'   example, by setting \code{project = list(medial = project_medial, 
#'   anterior = project_anterior)} and \code{query = list(medial = query_medial, 
#'   anterior = query_anterior)}. The names of list between project and query 
#'   should match. Objects in the list for \code{project} should be projects 
#'   created using \code{\link{saveproj}}. Objects in the list for query should
#'   be 3-dimensional array consist of the semi-landmark configurations(for
#'   \code{gpa}/ \code{nef} methods) or matrix containing the shape indices.
#' @param project path to a project (\code{.rds} file) to be read, saved using 
#'   \code{\link{saveproj}}, or a project object already read into R. if not
#'   given, interactive file selector will pop out to prompt user to select a
#'   \code{.rds} file (Windows only)
#' @param query path(s) to otolith images/ path(s) of folder containing the
#'   otolith images/ code{.tps} file containing the semi-landmark
#'   configurations/ p x k matrix or p x k x n array of semi-landmark
#'   configuration(s) to be predicted. If none is given, interactive file
#'   selector will pop out to prompt user to select images to be searched
#'   (Windows only)
#' @param multiview logical. Turn on mode of combination of different views for
#'   prediction. see details
#' @param type type of data to predict
#' @param method classification method. see Note
#' @param har numeric. optional. By default \code{har} range saved in the
#'   project is used. a different value could be set using this argument
#' @param pc numeric. optional. By default \code{pc} range saved in the project
#'   is used. a different value could be set using this argument
#' @param threshold numeric. optional. threshold value on posterior probalility 
#'  to reject the prediction. see \code{\link{threcv}}. Currently works with 
#'  \code{lda} only
#' @param reland logical. whether to do automatic re-arrangement of landmark- 
#'   configuration
#' @param tol numeric. max limit of distance (see \code{\link{sprdist}}) to
#'   which the automatic re-arrangement of landmark configuration should be
#'   carried out
#' @param fix numeric. for \code{gpa} method. index of landmarks that should
#'   stay fixed. default is \code{NULL}, i.e. all semi-landmarks are slid. see
#'   Note
#' @param mode  when \code{="search+pred"}, searching is included in addition to
#'   prediction using \code{\link{otosearch}}. automatically changed to 
#'   \code{"search+pred"} when \code{reland=TRUE}
#' @param saveresult logical. whether to save the result
#' @param search.plot logical. whether to plot the search results. used only when
#'  \code{reland = TRUE} or \code{mode = "search+pred"}
#' @param ... other arguments passed to \code{\link{agglda}}
#' @note 
#'  Currently the \code{method} supported are limited to \code{\link[MASS]{lda}}
#'  and \code{\link{agglda}} only.
#'  
#'  Because sliding semi-landmark method is not supported by 
#'  \code{\link{otosearch}}, thus the \code{project} used should contain 
#'  \code{gpa} object from non-sliding GPA transformation. Sliding will be 
#'  performed by \code{otopred} instead if \code{gpa} is the preferred method
#'  and \code{fix = NULL}.
#' @return matrix of prediction class and posterior probablity.
#' @importFrom MASS lda
#' @importFrom tree tree
#' @importFrom mixOmics plsda
#' @importFrom geomorph readland.tps
#' @importFrom Morpho relaxLM
#' @seealso 
#'  Which this function wraps: \code{\link{otosearch}}
#'  
#'  Methods of classifier: \code{\link{lda}}, \code{\link{plsda}}, 
#'    \code{\link{tree}}, \code{\link{agglda}}
#' @export

otopred <- function(project, query, multiview = FALSE, 
                    type = c("nef", "gpa", "des"),  
                    method = c("lda", "agglda", "tree", "plsda"), har, pc, 
                    threshold, reland = TRUE, tol = NULL, fix = NULL, 
                    mode = c("pred", "search+pred"), saveresult = FALSE, 
                    search.plot = FALSE, ...) {
  type <- match.arg(type)
  method <- match.arg(method)
  mode <- match.arg(mode)
  # support for combined views, but with less features
  if (multiview) { 
    if (any((names(project) == names(query)) == FALSE))
      stop("project list should match the query list")
    if (reland | mode == "search+pred")
      stop("search mode is not supported with multiview now, ", 
           "please set reland = FALSE")
    if (type == "nef" | type == "gpa") {
      nq <- sapply(query, function(x) dim(x)[3])
      if (diff(range(nq)) != 0) # do all query views have same n?
        stop("query from different views do not have same n")
      np <- sapply(project, function(x) dim(x$landmark)[3])
      if (diff(range(np)) != 0) # do all project views have same n?
        stop("landmark from different projects do not have same n") 
      if (type == "nef") { # check if har is null
        harp <- sapply(project, function(x) x$har)
        if (any(is.null(harp)))
          stop("all projects must contain har number")
      } else if (type == "gpa") { # check if pc is null
        pcp <- sapply(project, function(x) x$pc)
        if (any(is.null(pcp)))
          stop("all projects must contain pc number")
      }
      pq <- sapply(query, function(x) dim(x)[1])
      pp <- sapply(project, function(x) dim(x$landmark)[1])
      if (any((pq == pp) == FALSE))
        stop("project and query do not have the same number of landmarks")
    }
    view.names <- names(project)
    view.length <- length(view.names)
    class <- project[[1]]$class    
  } else {
    class <- project$class
    # get the project / query, copied from otosearch
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
          if (file.exists(project) & .ext(project) == 
                "rds" | .ext(project) == "RDS")
            project <- readRDS(project)
          else
            stop("Please provide right format of project")
        } else 
            stop("Please provide right format of project")
      }
  }
  # get the query/ determine the query type
  if (!multiview) {
    if (type == "nef" || type == "gpa") {
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
    } else if (type == "des") {
        if (!is.data.frame(query))
          stop("Please provide a dataframe containing shape descriptors")
    }
  }
  # taking care of orientation issue with otosearch()
  if (reland)
    mode <- "search+pred"
  else
    mode <- match.arg(mode)
  if (mode == "search+pred") {
    if (type == "nef" || type == "gpa") {
      if (is.matrix(query))
        query <- array(data=query, dim=c(dim(query)[1], 2, 1))
      # putting the search result in  
      searchresult <- otosearch(query = query, project = project, 
                                showplot = search.plot) 
      print(searchresult)
      if (reland) {
        for (i in 1:dim(query)[3]) {
          orient <- searchresult[[i]]$orient[1]
          if (!is.null(tol))
            if (searchresult[[i]]$rdist[1] > tol)
              orient <- "ori"      
          query[, , i] <- reland(query[, , i], orient)
        }
      }
    }
  }
  # initialize variables for multiview option
  if (multiview) {
    newdata <- vector("list", view.length)
    moddata <- vector("list", view.length)
  }
  # going thru type of choice
  if (type == "nef") {
    if(!multiview) {
      if (missing(har)) {
        if (!is.null(project$har)) 
          har <- project$har
        else 
          stop("please provide har value")
      }
      newdata <- selectdim(rNEF(query)$coeff, har=har)
      moddata <- selectdim(project$nef$coeff, har=har)
    } else {
      for (i in 1:view.length) {
        newdata[[i]] <- selectdim(rNEF(query[[i]])$coeff, har=project[[i]]$har)
      }
      for (i in 1:view.length) {
        moddata[[i]] <- selectdim(project[[i]]$nef$coeff, har = project[[i]]$har)
      }
    } 
  } else if (type == "gpa") {
      if (!multiview) {
        temp <- .GPApred(project=project, query = query, fix = fix, pc = pc)
        newdata <- temp$newdata
        moddata <- temp$moddata    
      } else {
        for (i in 1:view.length) {
          temp <- .GPApred(project[[i]], query[[i]], fix = fix, 
                                  pc = project[[i]]$pc)
          newdata[[i]] <- temp$newdata
          moddata[[i]] <- temp$moddata
        }
      }
  } else if (type == "des") {
      des.predictors <- c("AspRatio", "Circ", "Roundness", "Compactness", 
                    "Solidity", "Convexity", "Shape", "RFactor", "ModRatio")
      if (!multiview) {
        newdata <- query[, des.predictors]
        moddata <- project$des[, des.predictors]
      } else {
        for (i in 1:view.length) {
          newdata[[i]] <- query[[i]][, des.predictors]
          moddata[[i]] <- project[[i]]$des[, des.predictors]
        }
      }
  }
  # cbind for newdata(query) and moddata(training) for multiview
  if (multiview) {
    newdata <- do.call(cbind, newdata)
    moddata <- do.call(cbind, moddata)    
  }
  # going through method of choice
  if (method == "lda") {
    mod <- lda(x=moddata, class)
    prediction <- predict(mod, newdata=newdata)
    result <- data.frame(prediction$class)
    result$posterior <- round(apply(prediction$posterior, 1, max), 2)
    if (!missing(threshold))
      result$prediction.class[which(result$posterior < threshold)] <- NA    
  } else if (method == "agglda") {
    require(MASS)
    result <- agglda(X = moddata, Y = class, 
                     newdata = as.matrix(newdata), 
                     threshold = threshold, suppress = TRUE, ...)
    colnames(result)[1] <- c("prediction.class")
  } else if (method == "tree") {
    require(tree)
    # to be added 
  } else if (method == "plsda") {
    require(mixOmics)
    # to be added
  }
  # change the row names of result
  if (multiview)
    query <- query[[1]]
  if (type == "nef" || type == "gpa") {
    if (!is.null(dimnames(query)[[3]]))
      rownames(result) <- dimnames(query)[[3]]
    else
      rownames(result) <- paste0("Query", 1:dim(result)[1])
  } else if (type == "des") {
    rownames(result) <- query[, "Label"]
  }
  #save result
  if (saveresult) {
    name <- paste0("otopred-result(", Sys.Date(), ").txt")
    write.table(result, file = name, sep = "\t", quote = FALSE)
    cat("\n\n**The prediction result is saved at:", 
    paste(getwd(), name, sep = "/"), "\n\n")
  }
  return(result)
}