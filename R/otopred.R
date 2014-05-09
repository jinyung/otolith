#' Predict new samples based on classification models
#' 
#' @description Predict new, unknown samples using either shape decriptors or 
#'  semi-landmark configurations.
#' @details The functions include the \code{\link{otosearch}} algorithm as a mean 
#'  to guess the semi-landmarks arragement (the four types of arrangements, see
#'  \code{\link{reland}}), so that the user can predict the new, unknown samples
#'  even if the side and the direction of the otoliths are unknown. This could be 
#'  turned off using \code{reland=FALSE} to speed up the prediction if the side 
#'  and direction of the query is known, and is already in the same arrangements 
#'  as the dataset in project.
#' @param query dataframe of shape descriptor data or p x k x n array of 
#'  semi-landmarks
#' @param project the project data, see \code{\link{saveproj}}
#' @param type type of data to predict.
#' @param method classification method. see note.
#' @param har numeric. optional. By default \code{har} range saved in the project is 
#'  used. a different value could be set using this argument.  
#' @param pc numeric. optional. By default \code{pc} range saved in the project is 
#'  used. a different value could be set using this argument.  
#' @param threshold numeric. optional. threshold value on posterior probalility 
#'  to reject the prediction. see \code{\link{threcv}}. Currently works with 
#'  \code{lda} only.
#' @param reland logical. whether to do automatic re-arrangement of landmark-
#'  configuration.
#' @param tol numeric. max limit of distance (see \code{\link{sprdist}}) to which 
#'  the automatic re-arrangement of landmark configuration should be carried out.    
#' @param write logical. whether to save the result.
#' @note Currently the \code{method} supported is limited to \code{"lda"} only. 
#' @return matrix of prediction class and posterior probablity.
#' @importFrom MASS lda
#' @importFrom tree tree
#' @importFrom mixOmics plsda
#' @seealso 
#'  Which this function wraps: \code{\link{otosearch}}
#'  
#'  Methods of classifier: \code{\link{lda}}, \code{\link{plsda}}, 
#'    \code{\link{tree}}
#' @export

otopred <- function(query, project, type=c("nef", "gpa", "des"),  
                    method=c("lda", "tree", "plsda"), har=NULL, pc=NULL, 
                    threshold=NULL, reland=TRUE, tol=0.1, write=FALSE) {
  # solve the single specimen (matrix) problem
  if (reland == TRUE) {
    if (type == "nef" || type == "gpa") {
      if (is.matrix(query))
        query <- array(data=query, dim=c(dim(query)[1], 2, 1))
      # putting the search result in  
      searchresult <- otosearch(query=query, project=project) 
      for (i in 1:dim(query)[3]) {
        orient <- searchresult[[i]]$orient[1]
        if (searchresult[[i]]$rdist[1] > tol)
          orient <- "ori"      
        query[, , i] <- reland(query[, , i], orient)
      }
    }
  }  
  if (type == "des") {
    if (is.array(query))
      stop ("query data type wrong")
    newdata <- query[, c("AspRatio", "Circ", "Roundness", "Compactness", 
               "Solidity", "Convexity", "Shape", "RFactor", "ModRatio")]
    if (method == "lda")
      mod <- lda(x=project$des[, c("AspRatio", "Circ", "Roundness","Compactness",
             "Solidity", "Convexity", "Shape", "RFactor", "ModRatio")], 
             project$class)
  } else if (type == "gpa") {
    p <- dim(query)[1]
    n <- dim(query)[3]
    # superimpose on the database meanshape
    rot <- array (data=NA, dim= c(p, 2, n))
    for (i in 1:n)
      rot[, , i] <- fPsup(query[, , i], project$gpa$mshape)$Mp1
    query.score <- predict(project$gpa$pca, as.data.frame(t(matrix(rot, p * 2, n))))
    if (is.null(har)) 
      pc <- project$pc
    newdata <- selectdim(query.score, pc=pc)    
    if (method == "lda")
      mod <- lda(x=selectdim(project$gpa$score, pc=pc), project$class)  
  } else if (type == "nef") {
    if (is.null(har)) 
      har <- project$har
    newdata <- selectdim(rNEF(query)$coeff, har=har)
    if (method == "lda")
      mod <- lda(x=selectdim(project$nef$coeff, har=har), project$class)      
  }  
  if (method == "lda") {
    prediction <- predict(mod, newdata=newdata)
    result <- data.frame(prediction$class)
    result$posterior <- round(apply(prediction$posterior, 1, max), 2)
    if (!is.null(threshold))
      result$prediction.class[which(result$posterior < threshold)] <- NA
    if (type == "nef" || type == "gpa") {
      if(!is.null(dimnames(query)[[3]]))
        rownames(result) <- dimnames(query)[[3]]
      else
        rownames(result) <- paste0("Query", 1:dim(result)[1])
    } else if (type == "des") {
      rownames(result) <- query[,"Label"]
    }
  } else if (method == "tree") {
    # to be added 
  } else if (method == "plsda") {
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