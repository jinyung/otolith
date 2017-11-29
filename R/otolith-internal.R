#' @keywords internal 

#---------------------------------------------------------

# function to get extension from a list of full/ half path
# how it works: split character at ".", take the last splitted part as ext
# fixed bug for file name with "." in it 
# using vapply to avoid loop and specify output type
# @param pathlist list of path of length >=1
.ext <- function (pathlist) {
  vapply(pathlist, function(x) rev(unlist(strsplit(x, "[.]")))[1], 
         character(1), USE.NAMES = FALSE)
}

#---------------------------------------------------------

# function to get file names from a list of full path how it works: split
# character at "//" OR "/" (regular expression accepts logical conditions), take
# the last part, and if remove.ext=TRUE, substitute the ".<ext>" with empty
# string
# @param pathlist list of path of length >=1
# @param remove.ext specify whether to remove the extension
# @reference
#   http://stackoverflow.com/questions/19424709
# actually same as basename
.getfilename <- function (pathlist, remove.ext = TRUE) {
  output <- vapply(pathlist, function(x) rev(unlist(strsplit(x, 
            split = "[\\]|/")))[1], character(1), USE.NAMES = FALSE)
  if(remove.ext)
    output <- mapply(gsub, paste0(".", .ext(output)), "", output, 
                     USE.NAMES = FALSE)
  return(output)
}

#---------------------------------------------------------

# function to get directory of where the file is from a full path name/ list
# @param path a single path name or a path list
# @note if a path list is provided, extraction is based on the first file on 
#   the list only
# @return path
# actually same as dirname
.getdir <- function(path) {  
  if (length(path) > 1) path <- path[1]
  paste(rev(rev(unlist(strsplit(path, "[\\]|/")))[-1]), collapse="/") 
}

#---------------------------------------------------------

# function to determine if an object is a otolith project
# @param x object to be determined
# @return logical
.isproj <- function (x) {
  temp <- c("landmark","des", "nef", 
            "gpa", "pc", "har", "class") %in% unlist(attributes(x))
  !any(temp == FALSE)
}

#---------------------------------------------------------

# internal function wrapped under agglda(), do k-fold sampling aggregation

.kflda <- function(X, Y, newdata = NULL, k = 5, 
                   prior = c("equal", "proportion")) {
  # require(MASS)
  Y <- factor(Y)
  dat <- data.frame(class=Y)
  dat$predictor <- as.matrix(X)   
  prior <- match.arg(prior)
  alltestingindices <- kfcv.testing(dim(X)[1], k = k)  
  # initiate for aggregating class and posterior 
  if (is.null(newdata)) {
    agg.class <- factor(levels = levels(Y))
    agg.posterior <- matrix(data = NA, dim(X)[1], length(levels(Y)), 
                            dimnames = list(rownames(X), levels(Y)))
  } else {
    agg.class <- as.data.frame(matrix(data = NA, k, dim(newdata)[1], dimnames = 
                               list(paste0("fold", 1:k), rownames(newdata))))
    agg.posterior <- array(data = NA, dim = c(dim(newdata)[1], 
                                              length(levels(Y)), k), 
                           dimnames = list(rownames(newdata), levels(Y), 
                                           paste0("fold", 1:k)))
  }
  # going thru the folds
  for (i in 1:k) {
    testingindices <- alltestingindices[[i]]
    train <- dat[-testingindices, ]
    test <- dat[testingindices, ]
    trainlength <- length(levels(train$class))
    switch(prior, equal = assign("prior.", rep(1 / trainlength, trainlength)), 
           proportion = assign("prior.", sapply(summary(train$class), 
                                                function (x) x / dim(train)[1])))
    mod.i <- MASS::lda(class ~ predictor, data = train, prior = prior.)
    if (!is.null(newdata)) {
      temp <- data.frame(newdata)
      temp$predictor <- as.matrix(newdata)
      test <- temp 
    }
    prediction.temp <- predict(mod.i, test)
    prediction.i <- prediction.temp$class
    posterior <- prediction.temp$posterior
    if (is.null(newdata))
      agg.class[testingindices] <- prediction.i
    else 
      agg.class[i, ] <- prediction.i
    # if any class did not appear in the training set in this sub-model
    if (dim(posterior)[2] < dim(agg.posterior)[2]) {
      posmatindex <- sapply(colnames(posterior), function(x) 
                     which(x == levels(train$class)))
      if (is.null(newdata)) # those without the predictions remain NA
        agg.posterior[testingindices, posmatindex] <- posterior 
      else
        agg.posterior[, posmatindex, i] <- posterior
    } else {
      if (is.null(newdata))
        agg.posterior[testingindices, ] <- posterior 
      else
        agg.posterior[, , i] <- posterior
    }     
  }
  return(list(agg.class = agg.class, agg.posterior = agg.posterior))
}

#---------------------------------------------------------

# predicting based on sliding GPA


.GPApred <- function(project, query, fix = NULL, pc) {
  # require(Morpho)  
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
    rot <- array (data = NA, dim = c(p, 2, n))
    if (dorelax) {
      sink("NUL") # suppress procSym cat()
      gpanew <- rGPA(project$landmark, fix = fix, class = project$class)
      sink() # suppress procSym cat()
    } else
      gpanew <- project$gpa
    for (i in 1:n) {
      if (dorelax) {
        sink("NUL") # suppress relaxLM cat(), only works in Windows
        rot[, , i] <- Morpho::relaxLM(query[, , i], gpanew$mshape, 
                                      SMvector = slide, outlines = outline)
        sink() # suppress relaxLM cat()
      } else {
        rot[, , i] <- query[, , i]
      }
      rot[, , i] <- Morpho::rotonto(gpanew$mshape, rot[, , i], scale = TRUE)$yrot
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
    moddata <- selectdim(gpanew$score, pc=pc) 
    sink()
    return(list(newdata = newdata, moddata = moddata))
}