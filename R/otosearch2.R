#' [Otosearch - internal function] 
#' 
#' @description search one specimen against whole database 
#' @param specimen matrix of configuration of query 
#' @param database array of configuration to search against
#' @param species species/ grouping of \code{database} in the same order
#' @param show how many search results to be shown (in order of ranking)
#' @return search result
#' @importFrom Morpho rotonto kendalldist
#' @seealso 
#'  Functions that wraps this function: \code{\link{otosearch3}}
#'  
#'  Which this function wraps: \code{\link{rotonto}}, \code{\link{kendalldist}}

otosearch2 <- function(specimen, database, species=NULL, show=5) {
  require(Morpho)
  if (dim(specimen)[1] != dim(database)[1])
    stop("the number of semi-landmarks differs between 
         the specimen and the database")
  rdist <- numeric()
  n <- dim(database)[3]
  for (i in 1: n) {
    temp <- rotonto(x=specimen, y=database[, , i])
    rdist[i] <- kendalldist(temp$X, temp$Y)
  }
  result <- round(data.frame(rdist), 5)
  result$label <- dimnames(database)[[3]]
  result$species <- species
  result <- result[order(result$rdist), ]
  rownames(result) <- paste0("rank", 1:n)
  return(result = result[1:show, ])
}