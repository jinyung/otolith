#' Otosearch - internal 
#' 
#' @description search one specimen against whole database 
#' @param specimen matrix of configuration of query 
#' @param database array of configuration to search against
#' @param species species/ grouping of \code{database} in the same order
#' @param show how many search results to be shown (in order of ranking)
#' @return search result
#' @importFrom Morpho rotonto kendalldist
#' @export
#' @seealso 
#'  Functions that wraps this function: \code{\link{otosearch3}}
#'  
#'  Which this function wraps: \code{\link{rotonto}}, \code{\link{kendalldist}}
#' @keywords internal

otosearch2 <- function(specimen, database, 
                       species = NULL, show = 5) {
  require(Morpho)
  if (dim(specimen)[1] != dim(database)[1])
    stop("the number of semi-landmarks differs between 
         the specimen and the reference specimens in database")
  n <- dim(database)[3]
  rdist <- numeric()
  for (i in 1:n) {
    temp <- rotonto(x = specimen, y = database[, , i], scale = TRUE, 
                    reflection = FALSE)
    rdist[i] <- kendalldist(temp$yrot, specimen)
  } # in rotonto, x is fixed, y is rotated
  result <- round(data.frame(rdist), 5)
  result$label <- dimnames(database)[[3]]
  result$species <- species
  result <- result[order(result$rdist), ]
  rownames(result) <- paste0("rank-", 1:n)
  return(result[1:show, ])
}