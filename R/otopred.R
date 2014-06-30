#' Predict new samples based on classification models
#' 
#' @description Predict new, unknown samples using either shape decriptors or 
#'  semi-landmark configurations.
#' @details The functions include the \code{\link{otosearch}} algorithm as a
#'   means to guess the semi-landmarks arragement (the four types of
#'   arrangements, see \code{\link{reland}}), so that the user can predict the
#'   new, unknown samples even if the side and the direction of the otoliths are
#'   unknown. This could be turned off using \code{reland=FALSE} to speed up the
#'   prediction if the side and direction of the query is known, and is already
#'   in the same arrangements as the dataset in project.
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
#' @param write logical. whether to save the result
#' @param search.plot logical. whether to plot the search results. used only when
#'  \code{reland = TRUE} or \code{mode = "search+pred"}
#' @param ... other arguments passed to \code{\link{agglda}}
#' @note 
#'  Currently the \code{method} supported are limited to \code{"lda"} and
#'  \code{\link{agglda}} only. 
#'  
#'  Because sliding semi-landmark method is not supported by
#'  \code{\link{otosearch}}, thus the \code{project} used should contain \code{gpa}
#'  object from non-sliding GPA transformation. Sliding will be performed by
#'  \code{otopred} instead if \code{gpa} is the preferred method and 
#'  \code{fix = NULL}. 
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

otopred <- function(project, query, type = c("nef", "gpa", "des"),  
                    method = c("lda", "agglda", "tree", "plsda"), har, pc, 
                    threshold, reland = TRUE, tol = 0.2, fix = NULL, 
                    mode = c("pred", "search+pred"), write = FALSE, 
                    search.plot = FALSE, ...) {
  type <- match.arg(type)
  method <- match.arg(method)
  # get the project / query, copied from otosearch
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
          if (searchresult[[i]]$rdist[1] > tol)
            orient <- "ori"      
          query[, , i] <- reland(query[, , i], orient)
        }
      }
    }
  }  
  if (type == "nef") {
    if (missing(har)) {
      if (!is.null(project$har)) 
        har <- project$har
      else 
        stop("please provide har value")
    }
    newdata <- selectdim(rNEF(query)$coeff, har=har)
    if (method == "lda")
      mod <- lda(x=selectdim(project$nef$coeff, har=har), project$class)      
  } else if (type == "gpa") {
    # set sliding param
    p <- dim(query)[1]
    n <- dim(query)[3]
    slide <- 1:p
    dorelax <- TRUE
    if (!is.null(fix)) {
      if (length(slide[-fix]) == 0)
        dorelax <- FALSE
      else {
        slide <- slide[-fix]
        outline <- slide
      }
    } else {
      outline <- c(slide, 1)
    }
    # superimpose on the database meanshape
    rot <- array (data=NA, dim= c(p, 2, n))
    if (dorelax) {
      sink("NUL") # suppress procSym cat()
      gpanew <- rGPA(project$landmark, fix = fix, class = project$class)
      sink() # suppress procSym cat()
    } else
      gpanew <- project$gpa
    for (i in 1:n) {
      if (dorelax) {
        sink("NUL") # suppress relaxLM cat(), only works in Windows
        rot[, , i] <- relaxLM(query[, , i], gpanew$mshape, SMvector = slide, 
                              outlines = outline)
        sink() # suppress relaxLM cat()
      } else {
        rot[, , i] <- query[, , i]
      }
      rot[, , i] <- rotonto(gpanew$mshape, rot[, , i], scale = TRUE)$yrot
    }
    query.score <- predict(gpanew$pca, 
                           as.data.frame(t(matrix(rot, p * 2, n))))
    if (missing(pc)) {
      if (!is.null(project$pc))
        pc <- project$pc
      else
        stop("please provide pc value")
    }
    newdata <- selectdim(query.score, pc=pc)    
    if (method == "lda")
      mod <- lda(x=selectdim(gpanew$score, pc=pc), project$class)  
  } else if (type == "des") {
    if (is.array(query))
      stop ("query data type wrong")
    newdata <- query[, c("AspRatio", "Circ", "Roundness", "Compactness", 
               "Solidity", "Convexity", "Shape", "RFactor", "ModRatio")]
    if (method == "lda") {
      require(MASS)
      mod <- lda(x=project$des[, c("AspRatio", "Circ", "Roundness","Compactness",
             "Solidity", "Convexity", "Shape", "RFactor", "ModRatio")], 
             project$class)
    }
  }  
  if (method == "lda") {
    prediction <- predict(mod, newdata=newdata)
    result <- data.frame(prediction$class)
    result$posterior <- round(apply(prediction$posterior, 1, max), 2)
    if (!missing(threshold))
      result$prediction.class[which(result$posterior < threshold)] <- NA
    if (type == "nef" || type == "gpa") {
      if(!is.null(dimnames(query)[[3]]))
        rownames(result) <- dimnames(query)[[3]]
      else
        rownames(result) <- paste0("Query", 1:dim(result)[1])
    } else if (type == "des") {
      rownames(result) <- query[,"Label"]
    }
  } else if (method == "agglda") {
    require(MASS)
    traindat <- switch(type, gpa = selectdim(gpanew$score, pc=pc),
                       nef = selectdim(project$nef$coeff, har=har), 
                       des = project$des[, c("AspRatio", "Circ", 
                                             "Roundness","Compactness",
                                             "Solidity", "Convexity", 
                                             "Shape", "RFactor", "ModRatio")]
                       )
    result <- agglda(X = traindat, Y = project$class, newdata = newdata, 
                     threshold = threshold, suppress = TRUE, ...)
  } else if (method == "tree") {
    require(tree)
    # to be added 
  } else if (method == "plsda") {
    require(mixOmics)
    # to be added
  }
  if (write == TRUE) {
    name <- paste0("otopred-result(", Sys.Date(), ").txt")
    write.csv(result, file=name, append=TRUE)
    cat("\n\n**The prediction result is saved at:", 
    paste(getwd(), name, sep="/"))
  }
  return(result)
}